"""Shared plotting helpers: column names, axis decoration, marker styling.

These sit on top of the coordinate model in :mod:`cnplot.cnplot_genome_axis`; the
drawing functions take a :class:`~cnplot.cnplot_genome_axis.GenomeAxis` and render it,
while the rest are small utilities the plotting modules share.

The column-name constants describe the seg.ucn layout every plotting module
reads, so they live here rather than in any one of them. Reference readers live
in :mod:`cnplot.cnplot_io_utils`.
"""

import contextlib
import functools
import logging
import math
import threading

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends import backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import ScaledTranslation

from .cnplot_genome_axis import GenomeAxis
from .cnplot_style import FS_TICK

logger = logging.getLogger(__name__)

__all__ = [
    "BAF_LIM",
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
    "chrname_tiers",
    "decorate_genome_axis",
    "genome_label_room",
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

# Default limits for a BAF-like axis on [0, 1]. The pad keeps points at 0 or 1
# off the spine; ticks still run 0.0 to 1.0, the pad being too small to label.
BAF_LIM = (-0.01, 1.01)

# Display names for plotted quantities. A log2 ratio still reads as "RDR" on an
# axis, which is how the heatmap already labels its colorbar; the 1D rows follow
# so the two figures name the same quantity the same way.
AXIS_LABELS = {"log2RDR": "RDR"}


def display_label(name: str | None) -> str | None:
    """Canonical axis text for a quantity name.

    Args:
        name: Column or label name, or None.

    Returns:
        The display name from :data:`AXIS_LABELS`, or ``name`` unchanged when it
        has no canonical form. None passes through.
    """
    return name if name is None else AXIS_LABELS.get(name, name)


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


def _mb_tick_positions(genome_axis: GenomeAxis, step: float) -> tuple:
    """Intra-chromosome Mb tick positions in axis coordinates, and their labels.

    Ticks fall at every ``step`` base pairs within each chromosome, mapped through
    the axis segments so they land correctly under either gap layout; the tick at a
    chromosome's own start (genomic 0) is dropped.

    Args:
        genome_axis: Axis supplying the segment transform.
        step: Tick spacing in base pairs.

    Returns:
        (positions, labels): axis-coordinate x positions and their Mb integer
        labels, in axis order.
    """
    positions, labels = [], []
    for seg in genome_axis.segments:
        shift = seg.axis_start - seg.raw_start
        first = int(np.ceil(seg.raw_start / step)) * int(step)
        for p in range(first, int(seg.raw_end), int(step)):
            if p == 0:
                continue
            positions.append(p + shift)
            labels.append(p // 1_000_000)
    return positions, labels


# Mean glyph advance of the figure face, as a fraction of the font size.
_GLYPH_ADVANCE = 0.55
# Characters in the longest label of each kind: "chr22", "200".
_CHRNAME_CHARS = 5
_MB_CHARS = 3
# Mb labels sit this far under the chromosome-name size.
_MB_SHRINK = 3
# Space one line of text needs across its height, in font sizes.
_LINE = 1.1


def _label_extent_pt(n_chars: int, fontsize: float, rotation: float) -> float:
    """How far a rotated tick label reaches away from its axis, in points.

    Args:
        n_chars: Characters in the label.
        fontsize: Font size in points.
        rotation: Rotation in degrees.

    Returns:
        The extent perpendicular to the axis.
    """
    length = _GLYPH_ADVANCE * fontsize * n_chars
    angle = math.radians(rotation)
    return length * abs(math.sin(angle)) + fontsize * abs(math.cos(angle))


def _chrname_pitch(fontsize: float, rotation: float) -> float:
    """Spacing along the axis two neighbouring chromosome names need, in points.

    Parallel labels clear each other either side by side, a line apart across
    their height, or end to end, a label length apart along their baseline;
    the smaller of the two spacings is enough.

    Args:
        fontsize: Name font size in points.
        rotation: Label rotation in degrees.

    Returns:
        The spacing in points.
    """
    angle = math.radians(rotation)
    length = _GLYPH_ADVANCE * fontsize * (_CHRNAME_CHARS + 1)
    side_by_side = _LINE * fontsize / max(abs(math.sin(angle)), 1e-6)
    end_to_end = length / max(abs(math.cos(angle)), 1e-6)
    return min(side_by_side, end_to_end)


def genome_label_room(
    mb_ticks: bool = False,
    chrnames: bool = True,
    fontsize: float = FS_TICK,
    rotation: float = 60,
    pad: float = 10.0,
    tiers: int = 1,
) -> float:
    """Height the genome-axis labels between two stacked axes need, in inches.

    The gap between a row and the profile under it holds the profile's
    chromosome names and, in Mb mode, the row's Mb tick labels. A layout sizes
    that gap from this instead of from a fixed fraction, so the two sets of
    labels clear each other at any font size.

    Args:
        mb_ticks: The upper axes labels its x axis in Mb along its bottom.
        chrnames: The lower axes writes chromosome names along its top.
        fontsize: Chromosome-name font size; Mb labels are derived from it as
            :func:`decorate_genome_axis` does.
        rotation: Label rotation in degrees, of the names and the Mb labels.
        pad: Clear space added around the labels, in points.
        tiers: Rows of chromosome names, from :func:`chrname_tiers`.

    Returns:
        The room in inches.
    """
    room = pad
    if chrnames:
        room += tiers * _label_extent_pt(_CHRNAME_CHARS, fontsize, rotation)
        room += pad / 2
    if mb_ticks:
        room += _label_extent_pt(_MB_CHARS, max(fontsize - _MB_SHRINK, 4), rotation)
        room += pad / 2
    return room / 72.0


def _positions_pt(ax: plt.Axes, positions) -> np.ndarray:
    """Axis x positions as points along the axes, for spacing labels.

    Args:
        ax: Axes with its x limits already set.
        positions: X positions in data coordinates.

    Returns:
        The positions in points from the axes' left edge.
    """
    x0, x1 = ax.get_xlim()
    width_pt = ax.get_position().width * ax.figure.get_figwidth() * 72.0
    return (np.asarray(positions, dtype=float) - x0) / (x1 - x0) * width_pt


def chrname_tiers(
    genome_axis: GenomeAxis,
    axes_width_in: float,
    fontsize: float = FS_TICK,
    rotation: float = 60,
) -> int:
    """How many rows the chromosome names need on an axes this wide: 1 or 2.

    Rotated names need a line of clearance each. Names closer than that are
    nudged apart along the axis (:func:`decorate_genome_axis`), which works as
    long as the axes is wide enough to hold them all in one row; past that,
    every second crowded name goes on a second row and a layout has to leave
    room for it.

    Args:
        genome_axis: Axis supplying the name positions.
        axes_width_in: Width of the axes in inches.
        fontsize: Name font size in points.
        rotation: Name rotation in degrees.

    Returns:
        2 where one row cannot hold every name, else 1.
    """
    n = len(genome_axis.xtick_chrs)
    return 2 if n * _chrname_pitch(fontsize, rotation) > axes_width_in * 72.0 else 1


def _spread(at: np.ndarray, need: float) -> np.ndarray:
    """Positions moved the least that leaves neighbours ``need`` apart.

    The mean of a left-to-right and a right-to-left sweep: each sweep is the
    closest feasible placement that only pushes one way, and their mean stays
    feasible while sharing the displacement between both neighbours.

    Args:
        at: Positions, in any order.
        need: Smallest spacing wanted between neighbours.

    Returns:
        The adjusted positions, in the input's order.
    """
    order = np.argsort(at)
    fwd = at[order].astype(float)
    back = fwd.copy()
    for i in range(1, len(fwd)):
        fwd[i] = max(fwd[i], fwd[i - 1] + need)
    for i in range(len(back) - 2, -1, -1):
        back[i] = min(back[i], back[i + 1] - need)
    out = np.empty(len(at), dtype=float)
    out[order] = (fwd + back) / 2.0
    return out


def decorate_genome_axis(
    ax: plt.Axes,
    genome_axis: GenomeAxis,
    plot_chrname: bool = True,
    label_pos: str = "bottom",
    rotation: float = 60,
    fontsize: int = FS_TICK,
    hide_spines: bool = False,
    pad_axis: bool = False,
    mb_ticks: bool = False,
    mb_tick_step: float = 50_000_000,
) -> None:
    """Apply genome-wide x-axis limits, ticks, and chromosome labels.

    Two labeling modes. By default one chromosome-name tick sits at each
    chromosome midpoint. With ``mb_ticks`` the ticks instead mark intra-chromosome
    position every ``mb_tick_step`` base pairs (labeled in Mb, reset per
    chromosome), and chromosome names are drawn as compact text off the axis - the
    denser mode for genome-wide scatter QC.

    Args:
        ax: Axes to decorate.
        genome_axis: Genome axis supplying limits and tick positions.
        plot_chrname: Draw the chromosome labels; False clears them, as stacked
            panels do above the last row.
        label_pos: "bottom" or "top". CNP profiles label the top.
        rotation: Tick label rotation in degrees, applied to the chromosome names
            or, with ``mb_ticks``, to the Mb labels; rotating keeps the dense Mb
            labels from overlapping. Chromosome names stay upright in that mode.
            Chromosome names closer together than one line at this angle keep
            the angle and are nudged apart along the axis; on an axes too
            narrow to hold them in one row, every second crowded name goes on a
            second row. Mb labels too close to the previous one are left blank.
        fontsize: Label font size; defaults to the module-level FS_TICK.
        hide_spines: Hide all four spines.
        pad_axis: Use the padded limits instead of tight genome ends.
        mb_ticks: Mark position in Mb every ``mb_tick_step`` bp within each
            chromosome instead of one name-tick per chromosome. Default off keeps
            the chromosome-name behavior.
        mb_tick_step: Mb-tick spacing in base pairs. Only with ``mb_ticks``.

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

    on_top = label_pos == "top"

    if mb_ticks:
        positions, labels = _mb_tick_positions(genome_axis, mb_tick_step)
        ax.set_xticks(positions)
        # NB: a label nearer than one rotated line to the last one kept is left blank
        need = max(fontsize - _MB_SHRINK, 4) / max(
            math.sin(math.radians(rotation)), 0.5
        )
        kept = -math.inf
        labels = list(labels)
        for i, x in enumerate(_positions_pt(ax, positions)):
            if x - kept < need:
                labels[i] = ""
            else:
                kept = x
        ax.set_xticklabels(
            labels,
            fontsize=max(fontsize - _MB_SHRINK, 4),
            rotation=rotation,
            ha="left" if on_top else "right",
            rotation_mode="anchor",
        )
        ax.tick_params(axis="x", length=2, labeltop=on_top, labelbottom=not on_top)
        if plot_chrname:
            # NB: offset in points, so the names clear the Mb labels at any axes height
            mb_fs = max(fontsize - _MB_SHRINK, 4)
            clear = _label_extent_pt(_MB_CHARS, mb_fs, rotation) + 8.0
            for x, name in zip(
                genome_axis.xtick_chrs, genome_axis.xlab_chrs, strict=True
            ):
                ax.annotate(
                    str(name).upper().removeprefix("CHR"),
                    xy=(x, 1.0 if on_top else 0.0),
                    xycoords=ax.get_xaxis_transform(),
                    xytext=(0, clear if on_top else -clear),
                    textcoords="offset points",
                    ha="center",
                    va="bottom" if on_top else "top",
                    fontsize=fontsize,
                    annotation_clip=False,
                )
        return

    if not plot_chrname:
        ax.set_xticks([])
        ax.set_xticklabels([])
        return

    ax.set_xticks(genome_axis.xtick_chrs)
    ax.set_xticklabels(genome_axis.xlab_chrs, rotation=rotation, fontsize=fontsize)
    at = _positions_pt(ax, genome_axis.xtick_chrs)
    need = _chrname_pitch(fontsize, rotation)
    if len(at) > 1 and float(np.diff(np.sort(at)).min()) < need:
        ticks = ax.xaxis.get_major_ticks()
        width_pt = ax.get_position().width * ax.figure.get_figwidth() * 72.0
        if len(at) * need <= width_pt:
            # NB: nudged apart along the axis; no tick mark is drawn to betray it
            for tick, dx in zip(ticks, _spread(at, need) - at, strict=True):
                shift = ScaledTranslation(dx / 72.0, 0, ax.figure.dpi_scale_trans)
                for label in (tick.label1, tick.label2):
                    label.set_transform(label.get_transform() + shift)
        else:
            # NB: one row cannot hold them; alternate crowded names go up a row
            lift = _label_extent_pt(_CHRNAME_CHARS, fontsize, rotation) + 4.0
            last = [-math.inf, -math.inf]
            base = plt.rcParams["xtick.major.pad"]
            for i in np.argsort(at):
                tier = 0 if at[i] - last[0] >= need else 1
                last[tier] = at[i]
                ticks[i].set_pad(base + tier * lift)
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


@contextlib.contextmanager
def _compact_rasterized_images():
    """Make the PDF backend keep a compact copy of each rasterized image.

    Works around a matplotlib memory bug. For every rasterized artist,
    ``backend_mixed.MixedModeRenderer`` renders into a buffer the size of the
    whole figure (``width_in * dpi`` by ``height_in * dpi``) and then hands the
    PDF backend ``img[slice_y, slice_x]``, a *view* that keeps the whole buffer
    alive. ``PdfFile`` holds every image it is given until the file is
    finalized, so a multi-page PDF pins one full-figure buffer per rasterized
    artist: peak memory grows as ``n_artists * fig_w_in * fig_h_in * dpi ** 2``
    and is only released on close.

    Copying the view on the way in drops each buffer as soon as the next artist
    is drawn, leaving only the cropped pixels. The bytes written are unchanged.

    Copy on insert, not per page: one page can hold many rasterized artists, so
    deferring the copy still lets a page's worth of buffers accumulate.
    """
    original = backend_pdf.PdfFile.imageObject

    def imageObject(self, image):
        if getattr(image, "base", None) is not None:
            image = np.ascontiguousarray(image)
        return original(self, image)

    backend_pdf.PdfFile.imageObject = imageObject
    try:
        yield
    finally:
        backend_pdf.PdfFile.imageObject = original


class _TrueColorImage:
    """Stand-in for ``PIL.Image`` while ``PdfFile._writeImg`` runs.

    Every image it builds reports more colours than a palette can hold, so the
    backend takes its continuous-tone path: 8 bits per component, DeviceRGB,
    Flate with a PNG predictor. Everything else is the real module's.

    Attributes:
        calls: Images built so far, so the hook can tell it is still consulted.
    """

    def __init__(self, module):
        """Wrap the real ``PIL.Image`` module.

        Args:
            module: The ``PIL.Image`` module the backend imported.
        """
        self._module = module
        self.calls = 0

    def __getattr__(self, name):
        """Defer every other attribute to the real module.

        Args:
            name: Attribute requested by the backend.

        Returns:
            The real module's attribute.
        """
        return getattr(self._module, name)

    def fromarray(self, *args, **kwargs):
        """Build an image that declines to enumerate its colours.

        Args:
            *args: Passed to ``PIL.Image.fromarray``.
            **kwargs: Passed to ``PIL.Image.fromarray``.

        Returns:
            The image, with ``getcolors`` answering None.
        """
        self.calls += 1
        img = self._module.fromarray(*args, **kwargs)
        img.getcolors = lambda *a, **k: None
        return img


def _write_truecolor_images() -> bool:
    """Stop the PDF backend storing rasterized layers as palette images.

    With ``pdf.compression`` on, ``PdfFile._writeImg`` stores an RGB raster of
    256 or fewer distinct colours as an ``/Indexed`` image, and the index stream
    is packed at 1, 2 or 4 bits per pixel when the palette is small. Vector
    editors read sub-byte packed indices unreliably and draw the layer in the
    wrong colours. A copy-number raster has a few state colours, so it always
    takes that branch. matplotlib tracks the symptom in issue 27893 and has no
    setting for it.

    The backend decides by asking Pillow for the image's colours. For the
    duration of each ``_writeImg`` call the backend's ``Image`` module is swapped
    for :class:`_TrueColorImage`, which sends the raster down the backend's own
    continuous-tone path. Nothing of matplotlib's is reimplemented, compression
    stays on for every stream, and the rendered pixels are unchanged.

    Installed on import rather than inside :class:`FigureSaver`: the backend
    writes images when the file is finalized, not when a figure is saved, and a
    caller that opens its own ``PdfPages`` or calls ``Figure.savefig`` must get
    the same output.

    ``_writeImg`` and the backend's ``Image`` name are private to matplotlib. If
    either is missing the hook is skipped with a warning and output is whatever
    that matplotlib writes. If both exist but a compressed image is written
    without the backend building it through ``Image``, the hook no longer
    decides the encoding and says so once. The swap is held under a lock, so
    figures saved from several threads get the same encoding.

    References:
        https://github.com/matplotlib/matplotlib/issues/27893

    Returns:
        True if the hook was installed; False if it already was, or if this
        matplotlib lacks the internals it relies on.
    """
    original = getattr(backend_pdf.PdfFile, "_writeImg", None)
    pil_image = getattr(backend_pdf, "Image", None)
    if original is None or pil_image is None:
        logger.warning(
            "matplotlib PDF backend changed; rasterized layers may be written "
            "as palette images"
        )
        return False
    if getattr(original, "_cnplot_truecolor", False):
        return False
    proxy = _TrueColorImage(pil_image)
    lock = threading.RLock()
    warned = []

    @functools.wraps(original)
    def _writeImg(self, *args, **kwargs):
        with lock:
            seen = proxy.calls
            backend_pdf.Image = proxy
            try:
                return original(self, *args, **kwargs)
            finally:
                backend_pdf.Image = pil_image
                bypassed = plt.rcParams["pdf.compression"] and proxy.calls == seen
                if bypassed and not warned:
                    warned.append(True)
                    logger.warning(
                        "matplotlib PDF backend no longer builds images through "
                        "PIL.Image; rasterized layers may be written as palette "
                        "images"
                    )

    _writeImg._cnplot_truecolor = True
    backend_pdf.PdfFile._writeImg = _writeImg
    return True


_write_truecolor_images()


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
            with _compact_rasterized_images():
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
        """Close the multi-page PDF, if one is open. Safe to call again."""
        if self._pdf is not None:
            self._pdf.close()
            self._pdf = None

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
