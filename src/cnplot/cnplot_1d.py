"""1D genome scatter primitive: draw one value column on a given axes.

The multi-sample panel that stacks these rows over a shared profile lives in
:mod:`cnplot.cnplot_figures`; this module is just the per-axes renderer.

The renderer takes values that are already computed. A column named by ``ycol``
holds whatever is being plotted - RDR, log2RDR, BAF, mhBAF, fractional copy
number - and no scaling, log, or allele split happens here. Expected values
arrive the same way, in their own table.

That contract is what lets one function serve both bulk and single-cell callers:
a bulk expectation is a proportion-weighted mixture over all clones, a
single-cell one is a single pure clone, and neither belongs in a plotting
library.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection

from .cnplot_genome_axis import GenomeAxis
from .cnplot_utils import (
    EXP_PREFIX,
    decorate_genome_axis,
    draw_chr_boundaries,
    draw_segment_boundaries,
    resolve_colors,
    resolve_marker_size,
)

logger = logging.getLogger(__name__)

__all__ = [
    "plot_scatter_1d",
]


# =============================================================================
# Helpers
# =============================================================================


def _merge_exp_lines(
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    chroms: np.ndarray,
    tol: float = 1.0,
) -> list:
    """Join runs of equal expected value into single line segments.

    Consecutive bins sharing a value on one chromosome become one horizontal
    line, so a segment-wide expectation draws as one stroke rather than one per
    bin. Bins with a non-finite position or value are dropped first, which also
    breaks a run.

    A run also breaks where the next bin does not start where the last one ended.
    Merging across such a break would assert an expected copy number over
    sequence no bin covers - visible whenever gaps keep their width, where the
    line would otherwise run straight through a centromere.

    Args:
        starts: (n,) left edges in axis coordinates.
        ends: (n,) right edges in axis coordinates.
        values: (n,) expected value per bin.
        chroms: (n,) chromosome per bin; runs never cross one.
        tol: Largest gap in axis units still counted as abutting, absorbing
            rounding when a shrunk axis makes neighbours meet exactly.

    Returns:
        List of [(x0, y), (x1, y)] pairs for a ``LineCollection``.
    """
    values = np.asarray(values, dtype=float)
    starts = np.asarray(starts, dtype=float)
    ends = np.asarray(ends, dtype=float)
    keep = np.isfinite(starts) & np.isfinite(ends) & np.isfinite(values)
    starts, ends, values = starts[keep], ends[keep], values[keep]
    chrom_arr = np.asarray(chroms)[keep]

    lines = []
    i = 0
    n = len(values)
    while i < n:
        j = i + 1
        while (
            j < n
            and values[j] == values[i]
            and chrom_arr[j] == chrom_arr[i]
            and starts[j] - ends[j - 1] <= tol
        ):
            j += 1
        lines.append([(starts[i], values[i]), (ends[j - 1], values[i])])
        i = j
    return lines


def _expected_column(ycol: str, group=None) -> str:
    """Name the expected-value column for one group.

    Args:
        ycol: Observed value column.
        group: Group id, or None for a table holding a single unnamed group.

    Returns:
        "exp_<ycol>_<group>", or "exp_<ycol>" when ``group`` is None.
    """
    return f"{EXP_PREFIX}{ycol}" + (f"_{group}" if group is not None else "")


def _expected_lines(
    expected_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    ycol: str,
    group=None,
) -> list:
    """Map an expected-value table onto the axis and merge it into lines.

    The table may be coarser than the observations: expected values are constant
    across a segment, so passing segments rather than bins is both cheaper and
    exactly right, and mapping through the same axis keeps the two aligned.

    Args:
        expected_df: Table with "#CHR", "START", "END" and the expected column.
        genome_axis: Axis to map through.
        ycol: Observed value column the expectation belongs to.
        group: Group id selecting the column.

    Returns:
        Line segments for a ``LineCollection``, empty when the column is absent,
        which is how a group with no expectation drops out.
    """
    col = _expected_column(ycol, group)
    if col not in expected_df.columns:
        return []
    coords = genome_axis.build_coordinates(expected_df)
    return _merge_exp_lines(
        coords.starts,
        coords.ends,
        expected_df[col].to_numpy(dtype=float),
        expected_df["#CHR"].to_numpy(),
    )


# =============================================================================
# 1D scatter
# =============================================================================


def plot_scatter_1d(
    ax: plt.Axes,
    obs_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    ycol: str,
    expected_df: pd.DataFrame | None = None,
    group=None,
    hue=None,
    palette=None,
    colors=None,
    alphas=None,
    keep_col: str | None = None,
    filtered_color="red",
    show_legend: bool = True,
    markersize: float | None = None,
    edgecolor: str = "none",
    edgewidth: float = 0.0,
    ylim: tuple | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    href: float | None = None,
    reverse_y: bool = False,
    plot_chrname: bool = True,
    mb_ticks: bool = False,
    mb_tick_step: float = 50_000_000,
    show_gaps: bool = True,
    exp_color="black",
    exp_linewidth: float = 1.5,
    bd_linewidth: float = 1.0,
    rasterized: bool = True,
) -> None:
    """Scatter one value column along the genome on a single axes.

    Bins are placed by coordinate through ``genome_axis``, so ``obs_df`` needs no
    relationship to whatever the axis was built from, and rows falling outside
    the plotted regions are dropped with a warning rather than misplaced.

    Args:
        ax: Axes to draw on.
        obs_df: Bins with "#CHR", "START", "END" and the ``ycol`` column.
        genome_axis: Axis placing the bins and supplying the decoration.
        ycol: Column to plot; its values are used as given.
        expected_df: Expected values with "#CHR", "START", "END" and an
            ``exp_<ycol>_<group>`` column, drawn as a step overlay. None skips
            it.
        group: Group id selecting the expected column.
        hue: (n,) categorical label per bin, or a column name in ``obs_df``.
        palette: {label: color} for ``hue``.
        colors: Explicit per-bin colors, taking precedence over ``hue``.
        alphas: (n,) alpha per bin, or a column name in ``obs_df``.
        keep_col: Bool column marking which bins passed upstream filtering. False
            rows are still drawn, in ``filtered_color``, so a filter can be judged
            against the data it removed. Kept rows keep their ``hue`` colors.
            Bins outside the plotted regions are a separate matter and are dropped
            by the axis whatever this says.
        filtered_color: Color for the rows ``keep_col`` marks False.
        show_legend: Label the kept and filtered counts. Only with ``keep_col``.
        markersize: Marker size; None picks one from the point count.
        edgecolor: Marker edge color.
        edgewidth: Marker edge width.
        ylim: Fixed (lo, hi), or None to leave the autoscale alone.
        ylabel: Y-axis label.
        title: Axes title, drawn left-aligned.
        href: Y value for a dotted grey reference line, e.g. 0.5 on a BAF row.
        reverse_y: Invert the y-axis.
        plot_chrname: Draw chromosome labels under the axes.
        mb_ticks: Label x in Mb every ``mb_tick_step`` bp within each chromosome
            instead of one name-tick per chromosome; see
            :func:`~cnplot.cnplot_utils.decorate_genome_axis`.
        mb_tick_step: Mb-tick spacing in base pairs. Only with ``mb_ticks``.
        show_gaps: Mark uncovered stretches, dashed where they are collapsed.
        exp_color: Color of the expected overlay.
        exp_linewidth: Width of the expected overlay.
        bd_linewidth: Width of the chromosome boundary lines.
        rasterized: Rasterize the point collections, keeping vector output small.

    Raises:
        ValueError: If ``ycol`` is missing from ``obs_df``.
    """
    if ycol not in obs_df.columns:
        raise ValueError(f"obs_df has no column {ycol!r}")
    if isinstance(hue, str):
        hue = obs_df[hue].to_numpy()
    if isinstance(alphas, str):
        alphas = obs_df[alphas].to_numpy()

    coords = genome_axis.build_coordinates(obs_df)
    values = obs_df[ycol].to_numpy(dtype=float)
    keep = coords.mapped & np.isfinite(values)

    rgba = resolve_colors(
        len(obs_df), hue=hue, palette=palette, colors=colors, alphas=alphas
    )
    if markersize is None:
        markersize = resolve_marker_size(int(keep.sum()))

    if keep_col is None:
        drawn = [(keep, rgba[keep], None)]
    else:
        if keep_col not in obs_df.columns:
            raise ValueError(f"obs_df has no column {keep_col!r}")
        flag = obs_df[keep_col].to_numpy(dtype=bool)
        kept, dropped = keep & flag, keep & ~flag
        n_kept, n_dropped = int(kept.sum()), int(dropped.sum())
        # filtered first so kept points sit on top
        drawn = [
            (dropped, filtered_color, f"filtered ({n_dropped})"),
            (kept, rgba[kept], f"kept ({n_kept})"),
        ]
        logger.info("%s: %d kept, %d filtered", keep_col, n_kept, n_dropped)

    for mask, color, label in drawn:
        if not mask.any():
            continue
        ax.scatter(
            coords.positions[mask],
            values[mask],
            s=markersize,
            c=color,
            edgecolors=edgecolor,
            linewidths=edgewidth,
            label=label,
        )
    if keep_col is not None and show_legend:
        ax.legend(loc="upper right", fontsize=8, markerscale=2, framealpha=0.7)
    if rasterized:
        for coll in ax.collections:
            coll.set_rasterized(True)

    if href is not None:
        ax.axhline(href, color="grey", linestyle=":", linewidth=bd_linewidth, zorder=0)

    if expected_df is not None:
        lines = _expected_lines(expected_df, genome_axis, ycol, group)
        if lines:
            ax.add_collection(
                LineCollection(lines, linewidth=exp_linewidth, colors=exp_color)
            )

    draw_chr_boundaries(ax, genome_axis, linewidth=bd_linewidth)
    if show_gaps:
        draw_segment_boundaries(ax, genome_axis)

    if ylim is not None:
        ax.set_ylim(*ylim)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    if title is not None:
        ax.set_title(title, fontsize=12, fontweight="bold", loc="left")
    decorate_genome_axis(
        ax,
        genome_axis,
        plot_chrname=plot_chrname,
        mb_ticks=mb_ticks,
        mb_tick_step=mb_tick_step,
    )
    if reverse_y:
        ax.invert_yaxis()
