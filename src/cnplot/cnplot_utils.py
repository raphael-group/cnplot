"""Shared plotting helpers: column names, axis decoration, marker styling.

These sit on top of the coordinate model in :mod:`cnplot.cnplot_genome_axis`; the
drawing functions take a :class:`~cnplot.cnplot_genome_axis.GenomeAxis` and render it,
while the rest are small utilities the plotting modules share.

The column-name constants describe the seg.ucn layout every plotting module
reads, so they live here rather than in any one of them. Reference readers live
in :mod:`cnplot.cnplot_io_utils`.
"""

import logging
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .cnplot_genome_axis import GenomeAxis

logger = logging.getLogger(__name__)

__all__ = [
    "CN_PREFIX",
    "FigureSaver",
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
    "shade_regions",
]


# =============================================================================
# Column names
# =============================================================================

CN_PREFIX = "cn_"
U_PREFIX = "u_"
NORMAL_CLONE = "normal"
PI_VIOL_COL = "PI_VIOL"
SAMPLE_COL = "SAMPLE"


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

    A collapsed gap gets one dashed line; a gap with width gets two bounding
    lines and optional shading.

    Args:
        ax: Axes to draw on.
        genome_axis: Genome axis giving the gaps.
        include_edges: Also mark leading and trailing gaps. Off by default since
            collapsed ones sit under the solid chromosome boundary.
        color: Line color.
        linewidth: Line width.
        linestyle: Line style.
        alpha: Line opacity.
        shade: Fill color for gaps with width, or None for lines only.
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
            continue
        if shade is not None:
            ax.axvspan(gap.start, gap.end, color=shade, alpha=shade_alpha, zorder=0)
        for x in (gap.start, gap.end):
            ax.axvline(
                x,
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                alpha=alpha,
            )


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
    plot_clone_name: bool = True,
    clone_ploidies: dict | None = None,
    clone_props: dict | None = None,
) -> list:
    """Build y-tick labels for a clone-stacked panel.

    Returned top-to-bottom, reversing ``clones``, since the first clone sits at
    the bottom. A name like "clone2" renders as "Clone 2"; anything else is used
    verbatim. Each optional map adds a line only for the clones it contains.

    Args:
        clones: Clone names in stacking order, bottom row first.
        plot_clone_name: Write "Clone N" rather than a bare "N".
        clone_ploidies: Optional {clone: ploidy}, adding a "ploidy X" line.
        clone_props: Optional {clone: proportion} in [0, 1], adding "prop X%".

    Returns:
        One label string per clone, top row first.
    """
    ylabels = []
    for name in reversed(clones):
        m = re.fullmatch(r"clone(\d+)", str(name))
        if m:
            head = f"Clone {m.group(1)}" if plot_clone_name else m.group(1)
        else:
            head = str(name)
        lines = [head]
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
