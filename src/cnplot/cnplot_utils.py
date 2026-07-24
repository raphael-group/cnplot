"""Shared plotting helpers: column names, axis decoration, marker styling.

These sit on top of the coordinate model in :mod:`cnplot.cnplot_genome_axis`; the
drawing functions take a :class:`~cnplot.cnplot_genome_axis.GenomeAxis` and render it,
while the rest are small utilities the plotting modules share.

The column-name constants describe the seg.ucn layout every plotting module
reads, so they live here rather than in any one of them. Reference readers live
in :mod:`cnplot.cnplot_io_utils`.
"""

import logging

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .cnplot_genome_axis import GenomeAxis

logger = logging.getLogger(__name__)

__all__ = [
    "CN_PREFIX",
    "EXP_PREFIX",
    "FigureSaver",
    "MARKER_SIZE_LARGE",
    "MARKER_SIZE_SMALL",
    "MAX_NDOTS",
    "NORMAL_CLONE",
    "PI_VIOL_COL",
    "SAMPLE_COL",
    "U_PREFIX",
    "adaptive_dot_size",
    "decorate_genome_axis",
    "draw_chr_boundaries",
    "draw_segment_boundaries",
    "get_clone_ylabels",
    "get_transparency",
    "resolve_colors",
    "resolve_marker_size",
    "resolve_ylim",
    "resolve_ylim_scaled",
    "shade_regions",
]


# =============================================================================
# Column names
# =============================================================================

CN_PREFIX = "cn_"
U_PREFIX = "u_"
EXP_PREFIX = "exp_"
NORMAL_CLONE = "normal"
PI_VIOL_COL = "PI_VIOL"
SAMPLE_COL = "SAMPLE"

# Marker sizes for sparse and dense scatters, switched at MAX_NDOTS. The two
# values are what HATCHet (bulk, 16-22k bins) and Copytyping (pseudobulk, far
# fewer) each settled on independently; MAX_NDOTS sits between them and matches
# the n_ref already used by adaptive_dot_size.
MARKER_SIZE_SMALL = 2.0
MARKER_SIZE_LARGE = 20.0
MAX_NDOTS = 5000


# =============================================================================
# Axis decoration
# =============================================================================


def draw_chr_boundaries(
    ax: plt.Axes,
    genome_axis: GenomeAxis,
    include_ends: bool = False,
    ymax: float = 1.0,
    color="black",
    linewidth: float = 1.0,
    alpha: float = 1.0,
) -> None:
    """Draw vertical lines at chromosome boundaries.

    Args:
        ax: Axes to draw on.
        genome_axis: Genome axis giving the boundary positions.
        include_ends: Also draw at the genome start and end.
        ymax: Line top in axes fraction. Above 1 extends past the axes and
            disables clipping.
        color: Line color.
        linewidth: Line width.
        alpha: Line opacity.
    """
    xs = genome_axis.ch_coords if include_ends else genome_axis.chr_boundaries
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
    genome_axis: GenomeAxis,
    include_edges: bool = False,
    color="black",
    linewidth: float = 0.8,
    linestyle: str = "--",
    alpha: float = 0.6,
    shade: str | None = None,
    shade_alpha: float = 0.15,
) -> None:
    """Mark the stretches no whitelist segment covers.

    A dashed line stands in for width that was removed, so only a collapsed gap
    gets one. When gaps are kept the missing stretch is already visible as empty
    axis, and a line there would mark something the reader can see - optionally
    shade it instead.

    Args:
        ax: Axes to draw on.
        genome_axis: Genome axis giving the gaps.
        include_edges: Also mark leading and trailing gaps. Off by default since
            collapsed ones sit under the solid chromosome boundary.
        color: Line color.
        linewidth: Line width.
        linestyle: Line style.
        alpha: Line opacity.
        shade: Fill color for gaps that kept their width, or None to leave them
            blank. Ignored for collapsed gaps, which have nothing to fill.
        shade_alpha: Opacity of the fill.
    """
    gaps = genome_axis.gaps if include_edges else genome_axis.interior_gaps
    for gap in gaps:
        if gap.is_collapsed:
            ax.axvline(
                gap.start,
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                alpha=alpha,
            )
        elif shade is not None:
            ax.axvspan(gap.start, gap.end, color=shade, alpha=shade_alpha, zorder=0)


def decorate_genome_axis(
    ax: plt.Axes,
    genome_axis: GenomeAxis,
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
        genome_axis: Genome axis supplying limits and tick positions.
        plot_chrname: Draw one chromosome label per tick; False clears them,
            as stacked panels do above the last row.
        label_pos: "bottom" or "top". CNP profiles label the top.
        rotation: Label rotation in degrees.
        fontsize: Label font size.
        hide_spines: Hide all four spines.
        pad_axis: Use the padded limits instead of tight genome ends.

    Raises:
        ValueError: If ``label_pos`` is not "bottom" or "top".
    """
    if label_pos not in ("bottom", "top"):
        raise ValueError("label_pos must be 'bottom' or 'top'")

    ax.grid(False)
    if pad_axis:
        ax.set_xlim(genome_axis.axis_start, genome_axis.axis_end)
    else:
        ax.set_xlim(genome_axis.ch_coords[0], genome_axis.chr_end)
    ax.set_xlabel("")
    if hide_spines:
        for spine in ax.spines.values():
            spine.set_visible(False)

    if not plot_chrname:
        ax.set_xticks([])
        ax.set_xticklabels([])
        return

    ax.set_xticks(genome_axis.xtick_chrs)
    ax.set_xticklabels(genome_axis.xlab_chrs, rotation=rotation, fontsize=fontsize)
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
    genome_axis: GenomeAxis,
    regions: pd.DataFrame,
    color: str = "lightblue",
    alpha: float = 0.15,
    zorder: int = 0,
) -> None:
    """Shade genomic intervals, e.g. callable regions or a blacklist.

    Intervals are intersected with each axis segment, so a region spanning a
    shrunk gap shades as the two pieces actually drawn rather than one block
    across the join. Works on either layout.

    Args:
        ax: Axes to draw on.
        genome_axis: Genome axis to map the intervals through.
        regions: Intervals with "#CHR", "START", "END", as returned by
            :func:`cnplot.cnplot_io_utils.read_bed`.
        color: Fill color.
        alpha: Fill opacity.
        zorder: Draw order; keep below the data.
    """
    by_chrom = {ch: grp for ch, grp in regions.groupby("#CHR", sort=False)}
    for seg in genome_axis.segments:
        grp = by_chrom.get(seg.chrom)
        if grp is None:
            continue
        shift = seg.axis_start - seg.raw_start
        for start, end in zip(grp["START"], grp["END"], strict=True):
            lo = max(float(start), seg.raw_start)
            hi = min(float(end), seg.raw_end)
            if hi > lo:
                ax.axvspan(
                    lo + shift, hi + shift, color=color, alpha=alpha, zorder=zorder
                )


def get_clone_ylabels(
    clones: list,
    clone_ploidies: dict | None = None,
    clone_props: dict | None = None,
) -> list:
    """Build y-tick labels for a clone-stacked panel.

    Returned bottom-to-top, matching ascending y, so ``clones[0]`` ends up as the
    top row. Clone names are used verbatim. Each optional map adds a line only for
    the clones it contains.

    Args:
        clones: Clone names in stacking order, top row first.
        clone_ploidies: Optional {clone: ploidy}, adding a "ploidy X" line.
        clone_props: Optional {clone: proportion} in [0, 1], adding "prop X%".

    Returns:
        One label string per clone, top row first.
    """
    ylabels = []
    for name in reversed(clones):
        lines = [str(name)]
        if clone_ploidies is not None and name in clone_ploidies:
            lines.append(f"ploidy {round(clone_ploidies[name], 2)}")
        if clone_props is not None and name in clone_props:
            lines.append(f"prop {round(clone_props[name] * 100, 2)}%")
        ylabels.append("\n".join(lines))
    return ylabels


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

    Thresholds are per group and column, at ``one_tail`` and ``1 - one_tail``. A
    row fades if it is beyond a tail on any column, dimming cluster fringes.

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

    At ``n_ref`` points the size is ``s_base``, so denser plots get smaller
    markers.

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


def resolve_marker_size(
    n_points: int,
    small: float = MARKER_SIZE_SMALL,
    large: float = MARKER_SIZE_LARGE,
    max_ndots: int = MAX_NDOTS,
) -> float:
    """Choose a marker size from how many points will be drawn.

    A two-step switch rather than a ramp, so a figure has one of two looks
    instead of a size that drifts with the data. Resolve it once per figure from
    the total point count: applied per axes, two rows either side of
    ``max_ndots`` would differ tenfold. Use :func:`adaptive_dot_size` when a
    continuous scale is wanted instead.

    Args:
        n_points: Number of points to be drawn.
        small: Size at or above ``max_ndots``.
        large: Size below ``max_ndots``.
        max_ndots: Switch point.

    Returns:
        Marker size in points squared.
    """
    return float(large) if n_points < max_ndots else float(small)


def _finite(*arrays) -> np.ndarray:
    """Concatenate the finite values of several optional arrays.

    Args:
        *arrays: Arrays or None, each flattened and filtered.

    Returns:
        1-D array of every finite value, empty if there are none.
    """
    parts = []
    for arr in arrays:
        if arr is None:
            continue
        flat = np.asarray(arr, dtype=float).ravel()
        parts.append(flat[np.isfinite(flat)])
    return np.concatenate(parts) if parts else np.empty(0)


def resolve_ylim(
    values,
    expected=None,
    windows=((-2, 2), (-5, 5)),
) -> tuple:
    """Pick the tightest preset window holding the data.

    Windows are tried in order, so list them narrowest first; the last is the
    fallback when nothing fits. Expected values are taken into account so an
    overlay line cannot fall off the axes.

    Args:
        values: Observed values; non-finite entries ignored.
        expected: Expected values to keep visible, or None.
        windows: Candidate (lo, hi) limits, narrowest first.

    Returns:
        The first window containing every finite value, else the last one.
    """
    vals = _finite(values, expected)
    if vals.size == 0:
        return tuple(windows[-1])
    lo, hi = float(vals.min()), float(vals.max())
    for window in windows:
        if lo >= window[0] and hi <= window[1]:
            return tuple(window)
    return tuple(windows[-1])


def resolve_ylim_scaled(
    values,
    expected=None,
    scale: float = 1.1,
    min_top: float = 2.0,
    max_top: float | None = 6.0,
    floor_frac: float = -0.05,
) -> tuple:
    """Derive limits from the data for a ratio-like quantity.

    For values with a meaningful zero and no natural ceiling, e.g. RDR or
    fractional copy number. The top is the largest value with headroom, held to
    at least ``min_top`` so a flat diploid profile is not blown up, and capped at
    ``max_top`` so a few outliers cannot flatten everything else.

    Args:
        values: Observed values; non-finite entries ignored.
        expected: Expected values to keep visible, or None.
        scale: Headroom factor applied to the maximum.
        min_top: Smallest acceptable top.
        max_top: Largest acceptable top, or None for uncapped.
        floor_frac: Bottom as a fraction of the top, giving a little room below
            zero.

    Returns:
        (lo, hi) limits.
    """
    vals = _finite(values, expected)
    top = max(float(vals.max()) * scale, min_top) if vals.size else min_top
    if max_top is not None:
        top = min(top, max_top)
    return (floor_frac * top, top)


def resolve_colors(
    n: int,
    hue=None,
    palette=None,
    colors=None,
    alphas=None,
    default="steelblue",
) -> np.ndarray:
    """Resolve per-point colors to an RGBA array.

    Precedence is ``colors``, then ``hue`` through ``palette``, then ``default``.
    Doing this up front is what lets the 1D and 2D plots share one coloring: the
    alternative, reading resolved colors back off a drawn collection, forces one
    plot to run before the other.

    Args:
        n: Number of points.
        hue: (n,) categorical label per point, or None.
        palette: {label: color} map, or a color sequence taken in order of first
            appearance in ``hue``. Unmapped labels fall back to ``default``.
        colors: Explicit color, or one per point. Wins over ``hue``.
        alphas: (n,) alpha per point, written into the resolved colors, or None.
        default: Color for points with nothing else to go on.

    Returns:
        (n, 4) RGBA array.

    Raises:
        ValueError: If a resolved array does not have one color per point, or if
            ``hue`` is given without a ``palette``.
    """
    if colors is not None:
        rgba = mcolors.to_rgba_array(colors)
        if len(rgba) == 1:
            rgba = np.tile(rgba, (n, 1))
    elif hue is not None:
        if palette is None:
            raise ValueError("hue needs a palette")
        hue_arr = np.asarray(hue)
        lut = palette
        if not isinstance(palette, dict):
            keys = list(dict.fromkeys(hue_arr.tolist()))
            lut = dict(zip(keys, palette, strict=False))
        rgba = mcolors.to_rgba_array([lut.get(h, default) for h in hue_arr.tolist()])
    else:
        rgba = np.tile(mcolors.to_rgba_array(default), (n, 1))

    if len(rgba) != n:
        raise ValueError(f"resolved {len(rgba)} colors for {n} points")
    if alphas is not None:
        rgba = rgba.copy()
        rgba[:, 3] = np.asarray(alphas, dtype=float)
    return rgba


# =============================================================================
# Figure output
# =============================================================================


class FigureSaver:
    """Multi-figure writer, a drop-in for ``PdfPages``.

    Use as a context manager. "pdf" accumulates one multi-page
    ``{out_base}.pdf``; other formats write per page as
    ``{out_base}.p{i}.{ext}``. Several formats can be written at once.

    Attributes:
        out_base: Output path without extension.
        img_types: Formats being written.
        dpi: Raster resolution.
        transparent: Whether figures use a transparent background.
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
