"""Genome coordinate model and shared plotting helpers.

Every genome-wide plot in cnplot maps bins or sites onto a single horizontal axis.
This module owns that mapping. Callers build a :class:`GenomeAxis` once and pass it
to the drawing functions in the other modules, so a figure's panels always share
one coordinate system.

Two layouts are supported:

- "collapsed": concatenates whitelist segments only, so centromeres and blacklisted
  regions occupy no width. Segment joins are recorded in ``seg_coords`` and are
  normally drawn as dashed lines.
- "linear": lays out full chromosomes end to end using a chromosome-sizes table.
  Gaps keep their real width and ``seg_coords`` is empty.
"""

from collections import OrderedDict
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .cnplot_colormap import BLACK

__all__ = [
    "GenomeAxis",
    "FigureSaver",
    "adaptive_dot_size",
    "build_genome_axis",
    "decorate_genome_axis",
    "draw_chr_boundaries",
    "draw_segment_boundaries",
    "get_chr_sizes",
    "get_transparency",
    "read_bed_by_chr",
    "shade_regions",
]


# =============================================================================
# Genome coordinate model
# =============================================================================


@dataclass
class GenomeAxis:
    """Mapping from genomic bins or sites to positions on a single plot axis.

    Built by :func:`build_genome_axis`; do not construct directly. Per-bin arrays
    are aligned to the row order of the ``coords_df`` passed to the builder, and
    hold NaN for rows that fall outside the layout (a bin in no whitelist segment,
    or on a chromosome missing from the sizes table).

    Attributes:
        chrs: Chromosome names in plotting order.
        positions: (n,) bin midpoints in axis coordinates.
        abs_starts: (n,) bin left edges in axis coordinates.
        abs_ends: (n,) bin right edges in axis coordinates.
        ch_coords: len(chrs)+1 chromosome start offsets; the last entry is the
            genome end, so chromosome i spans ch_coords[i] to ch_coords[i + 1].
        seg_coords: Within-chromosome segment joins (collapsed layout only),
            where a centromere or blacklisted region was removed.
        axis_start: Left axis limit, including any leading ``chr_shift`` pad.
        axis_end: Right axis limit, including a trailing pad of the same width.
        chr_offsets: Per-chromosome start offset, for translating raw genomic
            coordinates onto the axis. Only meaningful in the linear layout,
            where the mapping is a single shift; None in the collapsed layout.
        x_edges: (m+1,) cell edges for ``pcolormesh``; None in the linear layout.
        col_bin_ids: (m,) row index of ``coords_df`` for each cell in ``x_edges``,
            with -1 for filler cells covering gaps. None in the linear layout.
    """

    chrs: list
    positions: np.ndarray
    abs_starts: np.ndarray
    abs_ends: np.ndarray
    ch_coords: list
    seg_coords: list
    axis_start: float
    axis_end: float
    chr_offsets: dict | None = None
    x_edges: np.ndarray | None = None
    col_bin_ids: list | None = field(default=None, repr=False)

    @property
    def chr_end(self) -> float:
        """Right edge of the last chromosome, excluding the trailing pad."""
        return self.ch_coords[-1]

    @property
    def chr_boundaries(self) -> list:
        """Internal chromosome boundaries, excluding both genome ends."""
        return self.ch_coords[1:-1]

    @property
    def xtick_chrs(self) -> list:
        """Chromosome midpoints, for placing one tick per chromosome."""
        return [
            (self.ch_coords[i] + self.ch_coords[i + 1]) / 2
            for i in range(len(self.chrs))
        ]

    @property
    def xlab_chrs(self) -> list:
        """Chromosome tick labels, aligned to :attr:`xtick_chrs`."""
        return list(self.chrs)


def build_genome_axis(
    coords_df: pd.DataFrame,
    wl_segments: pd.DataFrame | None = None,
    chrom_sizes: dict | None = None,
    contain: bool = False,
    chr_shift: float = 0.0,
) -> GenomeAxis:
    """Map bins or sites onto a single genome-wide axis.

    Pass ``wl_segments`` for the collapsed layout or ``chrom_sizes`` for the
    linear layout; exactly one is required. Chromosome order is taken from the
    order of first appearance in ``coords_df``, which is assumed sorted.

    Args:
        coords_df: Bins or sites with a "#CHR" column plus either "START" and
            "END" (intervals) or "POS" (points). Any aggregation level works,
            since no value columns are read. Row order defines the output arrays.
        wl_segments: Whitelist segments with "#CHR", "START", "END", sorted
            within each chromosome. Selects the collapsed layout.
        chrom_sizes: Chromosome name to length, e.g. from :func:`get_chr_sizes`.
            Selects the linear layout. Chromosomes absent from this table are
            dropped from the axis and their bins map to NaN.
        contain: Collapsed layout only. If True, a bin joins a segment only when
            fully inside it, matching integer-CN profiles where a partly covered
            bin has no well-defined state. If False, overlapping bins join and
            are clipped to the segment, which is what scatter and heatmap panels
            use. Ignored in the linear layout.
        chr_shift: Blank pad added before the first chromosome and after the
            last, in axis units. HATCHet's 1D plots use 10e6; profiles use 0.

    Returns:
        A :class:`GenomeAxis` for the requested layout.

    Raises:
        ValueError: If neither or both of ``wl_segments`` and ``chrom_sizes`` are
            given, or if ``coords_df`` lacks the required columns.
    """
    if (wl_segments is None) == (chrom_sizes is None):
        raise ValueError("pass exactly one of wl_segments (collapsed) or chrom_sizes")
    if "#CHR" not in coords_df.columns:
        raise ValueError("coords_df must have a '#CHR' column")
    has_interval = {"START", "END"}.issubset(coords_df.columns)
    if not has_interval and "POS" not in coords_df.columns:
        raise ValueError("coords_df must have START/END or POS columns")

    if wl_segments is not None:
        if not has_interval:
            raise ValueError("collapsed layout needs START/END, not POS")
        return _build_collapsed_axis(coords_df, wl_segments, contain, float(chr_shift))
    return _build_linear_axis(coords_df, chrom_sizes, float(chr_shift))


def _build_collapsed_axis(
    coords_df: pd.DataFrame,
    wl_segments: pd.DataFrame,
    contain: bool,
    chr_shift: float,
) -> GenomeAxis:
    """Concatenate whitelist segments, giving removed regions zero width.

    Walks each chromosome's segments in order, accumulating an offset so that a
    segment of length L occupies exactly L axis units regardless of where it sits
    on the chromosome. Also accumulates the ``pcolormesh`` grid: one cell per bin
    plus filler cells (id -1) wherever bins do not tile the segment.

    A segment join is recorded in ``seg_coords`` when it is not the chromosome's
    last segment, or when the first segment starts past position 0 (a removed
    leading telomere).

    Args:
        coords_df: Bins with "#CHR", "START", "END".
        wl_segments: Whitelist segments with "#CHR", "START", "END".
        contain: Require bins to be fully inside a segment rather than merely
            overlapping it.
        chr_shift: Leading and trailing pad in axis units.

    Returns:
        A :class:`GenomeAxis` with ``x_edges`` and ``col_bin_ids`` populated and
        ``chr_offsets`` None, the mapping being piecewise rather than a shift.
    """
    coords_df = coords_df.reset_index(drop=True)
    chrs = list(coords_df["#CHR"].unique())
    wl_chs = wl_segments.groupby("#CHR", sort=False)
    bins_chs = coords_df.groupby("#CHR", sort=False, observed=True)

    n = len(coords_df)
    positions = np.full(n, np.nan)
    abs_starts = np.full(n, np.nan)
    abs_ends = np.full(n, np.nan)

    x_edges = [chr_shift]
    col_bin_ids = []
    ch_coords = []
    seg_coords = []
    ch_offset = chr_shift

    for ch in chrs:
        ch_coords.append(ch_offset)
        wl_ch = wl_chs.get_group(ch)
        bins_ch = bins_chs.get_group(ch)
        n_seg = len(wl_ch)

        for si in range(n_seg):
            wl_row = wl_ch.iloc[si]
            wl_s, wl_e = wl_row["START"], wl_row["END"]
            seg_start = ch_offset
            seg_end = ch_offset + (wl_e - wl_s)
            is_join = si < n_seg - 1 or (si == 0 and wl_s > 0)

            if contain:
                in_seg = bins_ch[(bins_ch["START"] >= wl_s) & (bins_ch["END"] <= wl_e)]
            else:
                in_seg = bins_ch[(bins_ch["START"] < wl_e) & (bins_ch["END"] > wl_s)]

            if in_seg.empty:
                if seg_end > x_edges[-1]:
                    col_bin_ids.append(-1)
                    x_edges.append(seg_end)
                ch_offset = seg_end
                if is_join:
                    seg_coords.append(ch_offset)
                continue

            # clipping is a no-op when contain=True, so one path serves both modes
            bin_starts = (
                np.maximum(in_seg["START"], wl_s) - wl_s + ch_offset
            ).to_numpy(float)
            bin_ends = (np.minimum(in_seg["END"], wl_e) - wl_s + ch_offset).to_numpy(
                float
            )
            bin_ids = in_seg.index.to_numpy()

            abs_starts[bin_ids] = bin_starts
            abs_ends[bin_ids] = bin_ends
            positions[bin_ids] = (bin_starts + bin_ends) / 2

            ch_offset = seg_end
            if is_join:
                seg_coords.append(ch_offset)

            cur = seg_start
            if seg_start > x_edges[-1]:
                col_bin_ids.append(-1)
                x_edges.append(seg_start)
            for s, e, bid in zip(bin_starts, bin_ends, bin_ids, strict=True):
                if s > cur:
                    col_bin_ids.append(-1)
                    x_edges.append(s)
                    cur = s
                if e > cur:
                    col_bin_ids.append(int(bid))
                    x_edges.append(e)
                    cur = e
            if cur < seg_end:
                col_bin_ids.append(-1)
                x_edges.append(seg_end)

    ch_coords.append(ch_offset)
    return GenomeAxis(
        chrs=chrs,
        positions=positions,
        abs_starts=abs_starts,
        abs_ends=abs_ends,
        ch_coords=ch_coords,
        seg_coords=seg_coords,
        axis_start=chr_shift,
        axis_end=ch_offset + chr_shift,
        chr_offsets=None,
        x_edges=np.asarray(x_edges, dtype=float),
        col_bin_ids=col_bin_ids,
    )


def _build_linear_axis(
    coords_df: pd.DataFrame,
    chrom_sizes: dict,
    chr_shift: float,
) -> GenomeAxis:
    """Lay out full chromosomes end to end, preserving gap widths.

    Each chromosome contributes its full length from ``chrom_sizes``, so the
    mapping is a per-chromosome shift and centromeres keep their real width.

    Args:
        coords_df: Bins with "#CHR", "START", "END", or sites with "#CHR", "POS".
        chrom_sizes: Chromosome name to length. Chromosomes absent here are
            dropped and their rows map to NaN.
        chr_shift: Leading and trailing pad in axis units.

    Returns:
        A :class:`GenomeAxis` with ``chr_offsets`` populated, ``seg_coords``
        empty, and no ``pcolormesh`` grid.
    """
    coords_df = coords_df.reset_index(drop=True)
    seen = list(dict.fromkeys(coords_df["#CHR"].tolist()))

    chrs = []
    chr_offsets = {}
    ch_coords = []
    cum = chr_shift
    for ch in seen:
        size = chrom_sizes.get(ch)
        if size is None:
            continue
        chrs.append(ch)
        chr_offsets[ch] = cum
        ch_coords.append(cum)
        cum += size
    ch_coords.append(cum)

    ch = coords_df["#CHR"].to_numpy()
    offsets = np.array([chr_offsets.get(c, np.nan) for c in ch], dtype=float)
    if "POS" in coords_df.columns:
        pos = coords_df["POS"].to_numpy(dtype=float)
        starts = ends = pos
    else:
        starts = coords_df["START"].to_numpy(dtype=float)
        ends = coords_df["END"].to_numpy(dtype=float)

    abs_starts = offsets + starts
    abs_ends = offsets + ends
    return GenomeAxis(
        chrs=chrs,
        positions=(abs_starts + abs_ends) / 2,
        abs_starts=abs_starts,
        abs_ends=abs_ends,
        ch_coords=ch_coords,
        seg_coords=[],
        axis_start=chr_shift,
        axis_end=cum + chr_shift,
        chr_offsets=chr_offsets,
    )


# =============================================================================
# Axis decoration
# =============================================================================


def draw_chr_boundaries(
    ax: plt.Axes,
    axis: GenomeAxis,
    include_ends: bool = False,
    ymax: float = 1.0,
    color=BLACK,
    linewidth: float = 1.0,
    alpha: float = 1.0,
) -> None:
    """Draw vertical lines at chromosome boundaries.

    Args:
        ax: Axes to draw on.
        axis: Genome axis giving the boundary positions.
        include_ends: Also draw at the genome start and end.
        ymax: Line top in axes fraction. Values above 1 extend past the axes and
            disable clipping, which is how CNP profile rows are separated.
        color: Line color.
        linewidth: Line width.
        alpha: Line opacity.
    """
    xs = axis.ch_coords if include_ends else axis.chr_boundaries
    for x in xs:
        line = ax.vlines(
            x,
            ymin=0,
            ymax=ymax,
            transform=ax.get_xaxis_transform(),
            colors=color,
            linewidth=linewidth,
            alpha=alpha,
        )
        if ymax > 1:
            line.set_clip_on(False)


def draw_segment_boundaries(
    ax: plt.Axes,
    axis: GenomeAxis,
    color=BLACK,
    linewidth: float = 0.8,
    linestyle: str = "--",
    alpha: float = 0.6,
) -> None:
    """Draw dashed lines where the collapsed layout removed a region.

    No-op in the linear layout, where nothing was removed.

    Args:
        ax: Axes to draw on.
        axis: Genome axis giving the segment joins.
        color: Line color.
        linewidth: Line width.
        linestyle: Line style.
        alpha: Line opacity.
    """
    for x in axis.seg_coords:
        ax.axvline(
            x, color=color, linewidth=linewidth, linestyle=linestyle, alpha=alpha
        )


def decorate_genome_axis(
    ax: plt.Axes,
    axis: GenomeAxis,
    plot_chrname: bool = True,
    label_pos: str = "bottom",
    rotation: float = 60,
    fontsize: int = 8,
    hide_spines: bool = False,
    pad_axis: bool = False,
) -> None:
    """Apply genome-wide x-axis limits, ticks, and chromosome labels.

    Args:
        ax: Axes to decorate.
        axis: Genome axis supplying limits and tick positions.
        plot_chrname: Draw one chromosome label per tick. If False, ticks and
            labels are cleared, which is what stacked panels use above the last
            row.
        label_pos: "bottom" or "top". CNP profiles label the top so the rows
            below read as a stack.
        rotation: Label rotation in degrees.
        fontsize: Label font size.
        hide_spines: Hide all four spines, as CNP profiles do.
        pad_axis: Use ``axis_start``/``axis_end`` including the ``chr_shift``
            pad, instead of tight limits at the genome ends.

    Raises:
        ValueError: If ``label_pos`` is not "bottom" or "top".
    """
    if label_pos not in ("bottom", "top"):
        raise ValueError("label_pos must be 'bottom' or 'top'")

    ax.grid(False)
    if pad_axis:
        ax.set_xlim(axis.axis_start, axis.axis_end)
    else:
        ax.set_xlim(axis.ch_coords[0], axis.chr_end)
    ax.set_xlabel("")
    if hide_spines:
        for spine in ax.spines.values():
            spine.set_visible(False)

    if not plot_chrname:
        ax.set_xticks([])
        ax.set_xticklabels([])
        return

    ax.set_xticks(axis.xtick_chrs)
    ax.set_xticklabels(axis.xlab_chrs, rotation=rotation, fontsize=fontsize)
    on_top = label_pos == "top"
    ax.tick_params(
        axis="x",
        labeltop=on_top,
        labelbottom=not on_top,
        top=False,
        bottom=False,
    )


def shade_regions(
    ax: plt.Axes,
    axis: GenomeAxis,
    regions_by_chr: dict,
    color: str = "lightblue",
    alpha: float = 0.15,
    zorder: int = 0,
) -> None:
    """Shade genomic intervals, e.g. callable regions or a blacklist.

    Args:
        ax: Axes to draw on.
        axis: Genome axis in the linear layout.
        regions_by_chr: Chromosome to list of (start, end) in raw genomic
            coordinates, as returned by :func:`read_bed_by_chr`.
        color: Fill color.
        alpha: Fill opacity.
        zorder: Draw order; keep below the data.

    Raises:
        ValueError: If ``axis`` uses the collapsed layout, where a raw interval
            has no single axis span.
    """
    if axis.chr_offsets is None:
        raise ValueError("shade_regions needs a linear-layout GenomeAxis")
    for ch, offset in axis.chr_offsets.items():
        for start, end in regions_by_chr.get(ch, []):
            ax.axvspan(
                offset + start, offset + end, color=color, alpha=alpha, zorder=zorder
            )


# =============================================================================
# Marker styling
# =============================================================================


def get_transparency(
    df: pd.DataFrame,
    by: str,
    cols: tuple = ("RD", "BAF"),
    one_tail: float = 0.25,
    tail_alpha: float = 0.2,
    nontail_alpha: float = 1.0,
) -> pd.Series:
    """Fade points lying in the tails of their group's distribution.

    Thresholds are per group and per column, taken as order statistics at
    ``one_tail`` and ``1 - one_tail``. A row is faded when it is at or beyond a
    tail on any column, which dims cluster fringes and leaves cores opaque.

    Args:
        df: Rows to score.
        by: Grouping column, e.g. a cluster id.
        cols: Value columns tested for tail membership.
        one_tail: Fraction of each group in one tail, in [0, 0.5].
        tail_alpha: Alpha for tail rows.
        nontail_alpha: Alpha for the rest.

    Returns:
        Alpha per row, named "transparency" and indexed like ``df``.
    """
    thresholds = {}
    for key, grp in df.groupby(by, sort=False):
        per_col = {}
        for col in cols:
            arr = np.sort(grp[col].to_numpy())
            lo = arr[min(max(0, int(len(arr) * one_tail)), len(arr) - 1)]
            hi = arr[min(max(0, int(len(arr) * (1 - one_tail))), len(arr) - 1)]
            per_col[col] = (lo, hi)
        thresholds[key] = per_col

    def alpha_for(row):
        """Score one row against its group's thresholds.

        Args:
            row: Row of ``df``.

        Returns:
            ``tail_alpha`` if the row is beyond a tail on any column, else
            ``nontail_alpha``.
        """
        per_col = thresholds[row[by]]
        for col in cols:
            lo, hi = per_col[col]
            if row[col] <= lo or row[col] >= hi:
                return tail_alpha
        return nontail_alpha

    alphas = df.apply(alpha_for, axis=1)
    alphas.name = "transparency"
    return alphas


def adaptive_dot_size(
    n_points: int,
    s_base: float = 4,
    s_min: float = 0.5,
    s_max: float = 10,
    n_ref: int = 5000,
) -> float:
    """Scale scatter marker size inversely with point count.

    At ``n_ref`` points the size is ``s_base``; denser plots get smaller markers
    so genome-wide scatters stay readable across bin resolutions.

    Args:
        n_points: Number of points to be drawn.
        s_base: Size at ``n_ref`` points.
        s_min: Lower clamp.
        s_max: Upper clamp.
        n_ref: Reference point count.

    Returns:
        Marker size in points squared, clamped to [``s_min``, ``s_max``].
    """
    if n_points <= 0:
        return float(s_base)
    return float(np.clip(s_base * n_ref / n_points, s_min, s_max))


# =============================================================================
# IO
# =============================================================================


def get_chr_sizes(sz_file: str) -> OrderedDict:
    """Read a two-column chromosome-sizes file.

    Args:
        sz_file: Path to a whitespace-separated file of (chromosome, length),
            such as a ``.fai``-derived ``chrom.sizes``.

    Returns:
        Chromosome name to length, in file order.
    """
    chr_sizes = OrderedDict()
    with open(sz_file) as rfd:
        for line in rfd:
            if not line.strip():
                continue
            ch, size = line.split()[:2]
            chr_sizes[ch] = int(size)
    return chr_sizes


def read_bed_by_chr(bed_path: str | None) -> dict:
    """Read a BED file into per-chromosome interval lists.

    Args:
        bed_path: Path to a 3+ column BED file, or None.

    Returns:
        Chromosome to list of (start, end). Empty when ``bed_path`` is None, so
        callers can pass an optional path straight through.
    """
    if bed_path is None:
        return {}
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["#CHR", "START", "END"],
        dtype={"#CHR": str},
    )
    return {
        ch: list(zip(grp["START"], grp["END"], strict=True))
        for ch, grp in df.groupby("#CHR", sort=False)
    }


class FigureSaver:
    """Multi-figure writer, a drop-in for ``PdfPages``.

    Use as a context manager. "pdf" accumulates one multi-page
    ``{out_base}.pdf``; other formats are written per page as
    ``{out_base}.p{i}.{ext}``. Several formats can be written at once.

    Attributes:
        out_base: Output path without extension.
        img_types: Formats being written.
        dpi: Raster resolution.
        transparent: Whether figures are saved with a transparent background.
    """

    def __init__(
        self,
        out_base: str,
        img_type: str | list = "pdf",
        dpi: int = 300,
        transparent: bool = False,
    ):
        """Open the writer, creating the PDF immediately if one is requested.

        Args:
            out_base: Output path without extension.
            img_type: One format or a list, from "pdf", "png", "svg".
            dpi: Raster resolution.
            transparent: Save with a transparent background.
        """
        self.out_base = out_base
        self.img_types = [img_type] if isinstance(img_type, str) else list(img_type)
        self.dpi = dpi
        self.transparent = transparent
        self._page = 0
        self._pdf = PdfPages(f"{out_base}.pdf") if "pdf" in self.img_types else None

    def savefig(self, fig: plt.Figure, *args, **kwargs) -> None:
        """Write one figure to every configured format.

        Args:
            fig: Figure to write.
            *args: Ignored; accepted for ``PdfPages.savefig`` compatibility.
            **kwargs: Ignored; dpi and transparency come from the writer.
        """
        if self._pdf is not None:
            self._pdf.savefig(
                fig, dpi=self.dpi, bbox_inches="tight", transparent=self.transparent
            )
        for img_type in self.img_types:
            if img_type == "pdf":
                continue
            fig.savefig(
                f"{self.out_base}.p{self._page}.{img_type}",
                dpi=self.dpi,
                bbox_inches="tight",
                transparent=self.transparent,
            )
        self._page += 1

    def close(self) -> None:
        """Close the multi-page PDF, if one was opened."""
        if self._pdf is not None:
            self._pdf.close()

    def __enter__(self) -> "FigureSaver":
        """Enter the context manager.

        Returns:
            This writer.
        """
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        """Close the writer on context exit.

        Args:
            exc_type: Exception type, if one was raised.
            exc: Exception instance, if one was raised.
            tb: Traceback, if an exception was raised.
        """
        self.close()
