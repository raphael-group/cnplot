"""Integer copy-number profiles and legends over a seg.ucn table.

The plotting functions take the profile table directly, in the layout HATCHet
and Copytyping already write: "#CHR", "START", "END", a ``cn_<clone>`` column of
``"a|b"`` strings per clone, a matching ``u_<clone>`` proportion, and an optional
``PI_VIOL`` flag. Clone names, ordering, proportions, and the normal clone are all
derived from the columns, so callers pass one DataFrame rather than assembling
parallel arrays.

The supported surface is :func:`plot_cnv_profile`, which draws its own legend when
given a second axes. Everything else is support code: the seg.ucn readers are
importable for callers who want the parsed pieces, but they are not part of the
interface and may change.
"""

import logging
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from .cnplot_colormap import get_ascn_colors, get_cn_colors
from .cnplot_genome_axis import GenomeAxis
from .cnplot_utils import (
    CN_PREFIX,
    NORMAL_CLONE,
    PI_VIOL_COL,
    SAMPLE_COL,
    U_PREFIX,
    decorate_genome_axis,
    draw_chr_boundaries,
    draw_segment_boundaries,
    get_clone_ylabels,
)

__all__ = [
    "plot_cnv_profile",
]

logger = logging.getLogger(__name__)


# =============================================================================
# Helpers
# =============================================================================


def _clone_sort_key(name: str) -> tuple:
    """Order clone names: normal first, then clone1, clone2, ..., then the rest.

    Args:
        name: Clone name without the ``cn_`` prefix.

    Returns:
        Sort key tuple.
    """
    if name == NORMAL_CLONE:
        return (0, 0, "")
    m = re.fullmatch(r"clone(\d+)", name)
    return (1, int(m.group(1)), "") if m else (2, 0, name)


def _iter_bins(coords):
    """Walk the bins that landed on the genome_axis.

    Args:
        coords: :class:`~cnplot.cnplot_genome_axis.BinCoords` from
            :meth:`~cnplot.cnplot_genome_axis.GenomeAxis.build_coordinates`.

    Yields:
        (row, x0, width) for each mapped bin.
    """
    for row in np.flatnonzero(coords.mapped):
        x0 = float(coords.starts[row])
        w = float(coords.ends[row]) - x0
        if w > 0:
            yield int(row), x0, w


def _finish_clone_axis(
    ax: plt.Axes,
    genome_axis: GenomeAxis,
    clones: list,
    height: float,
    plot_chrname: bool,
    clone_ploidies: dict | None,
    clone_props: dict | None,
    show_clone_name: bool,
    ylabel: str | None,
    title: str | None,
    label_fontweight: str | None,
    tick_pad: float,
) -> None:
    """Apply the shared x/y axis styling for a clone-stacked profile.

    Args:
        ax: Axes to style.
        genome_axis: Genome axis supplying limits and chromosome ticks.
        clones: Clone names in stacking order, top row first.
        height: Total height of the stack in axes fraction.
        plot_chrname: Draw chromosome labels along the top.
        clone_ploidies: Optional {clone: ploidy} map.
        clone_props: Optional {clone: proportion} map.
        show_clone_name: Write "Clone N" rather than a bare "N".
        ylabel: Optional y-axis label.
        title: Optional axes title.
        label_fontweight: Font weight for clone labels, or None for default.
        tick_pad: Padding for the major y ticks, raised when A/B minor ticks
            occupy the space next to the genome_axis.
    """
    num_clones = len(clones)
    h = height / num_clones
    decorate_genome_axis(
        ax, genome_axis, plot_chrname=plot_chrname, label_pos="top", hide_spines=True
    )
    ax.set_yticks([h * (i + 0.5) for i in range(num_clones)])
    labels = get_clone_ylabels(
        clones,
        plot_clone_name=show_clone_name,
        clone_ploidies=clone_ploidies,
        clone_props=clone_props,
    )
    if label_fontweight is None:
        ax.set_yticklabels(labels, fontsize=8, va="center")
    else:
        ax.set_yticklabels(labels, fontsize=8, fontweight=label_fontweight, va="center")
    ax.set_ylim(0, height)
    ax.tick_params(
        axis="y", which="major", left=True, right=False, length=4, pad=tick_pad
    )
    if ylabel is not None:
        ax.set_ylabel(ylabel, rotation=0, ha="right", va="center")
    if title is not None:
        ax.set_title(title)


def _mirrored_bins(states: np.ndarray) -> np.ndarray:
    """Flag bins whose clones disagree on which allele is amplified.

    A bin is mirrored when some clone has a > b and another has a < b. A bin
    where every clone leans the same way is not mirrored: the integer-CN palette
    gives (a, b) and (b, a) one color, so a uniform lean records nothing but
    which allele was labelled A.

    Args:
        states: (n_bins, n_clones, 2) array from :func:`get_clone_states`.

    Returns:
        (n_bins,) bool array.
    """
    a, b = states[..., 0], states[..., 1]
    return (a > b).any(axis=1) & (a < b).any(axis=1)


def _draw_chevrons(
    ax: plt.Axes,
    x_start: float,
    unit: float,
    gap: float,
    y_high: float,
    y_low: float,
    direction: int,
    n_chev: int = 2,
    transform=None,
) -> None:
    """Draw the stacked chevrons that mark a clone in a mirrored bin.

    Args:
        ax: Axes to draw on.
        x_start: Left edge of the chevron block.
        unit: Width of one chevron.
        gap: Horizontal gap between chevrons.
        y_high: Top of the chevrons.
        y_low: Bottom of the chevrons.
        direction: 1 points the apex right (a > b), -1 points it left (a < b).
        n_chev: Number of chevrons.
        transform: Optional transform, e.g. ``ax.get_xaxis_transform()``.
    """
    y_mid = (y_high + y_low) / 2.0
    kwargs = {"transform": transform} if transform is not None else {}
    for i in range(n_chev):
        cx_left = x_start + i * (unit + gap)
        cx_right = cx_left + unit
        xs = (
            [cx_left, cx_right, cx_left]
            if direction > 0
            else [
                cx_right,
                cx_left,
                cx_right,
            ]
        )
        ax.plot(
            xs,
            [y_high, y_mid, y_low],
            color="black",
            linewidth=1.2,
            alpha=0.5,
            solid_capstyle="round",
            **kwargs,
        )


def _draw_mirror_swatch(
    ax: plt.Axes,
    leg_x: float,
    pair_w: float = 2.0,
    pair_h: float = 0.6,
    fontsize: int = 10,
) -> float:
    """Draw the mirrored-CNA swatch: a "/"-hatched white box labeled "mirrored".

    Args:
        ax: Axes to draw on.
        leg_x: Left edge of the swatch.
        pair_w: Swatch width.
        pair_h: Swatch height.
        fontsize: Label font size.

    Returns:
        The x position just past the swatch.
    """
    ax.add_patch(
        Rectangle(
            (leg_x, 0.0),
            pair_w,
            pair_h,
            facecolor="white",
            edgecolor="black",
            hatch="/",
        )
    )
    ax.text(
        leg_x + pair_w / 2.0,
        -0.2,
        "mirrored",
        ha="center",
        va="top",
        fontsize=fontsize,
    )
    return leg_x + pair_w


# =============================================================================
# seg.ucn parsing
# =============================================================================


def select_sample(
    seg_df: pd.DataFrame,
    sample_id: str | None = None,
    col: str = SAMPLE_COL,
) -> pd.DataFrame:
    """Take one sample's rows from a profile table.

    A seg.ucn file may hold several samples stacked in one table. Profiles are
    per sample, so one has to be chosen; without ``sample_id`` the first in file
    order wins.

    Args:
        seg_df: Profile table, with or without a ``SAMPLE`` column.
        sample_id: Sample to keep. None takes the first one present.
        col: Column holding the sample id.

    Returns:
        The matching rows, re-indexed. The input is returned unchanged when it
        carries no ``SAMPLE`` column.

    Raises:
        ValueError: If ``sample_id`` is given but absent, either from the table
            or because there is no ``SAMPLE`` column to match on.
    """
    if col not in seg_df.columns:
        if sample_id is not None:
            raise ValueError(f"seg_df has no {col!r} column to select {sample_id!r}")
        return seg_df
    samples = list(dict.fromkeys(seg_df[col]))
    wanted = samples[0] if sample_id is None else sample_id
    if wanted not in samples:
        raise ValueError(f"sample {sample_id!r} not in {col}; have {samples}")
    if sample_id is None and len(samples) > 1:
        logger.info("seg_df holds %d samples; plotting %r", len(samples), wanted)
    return seg_df[seg_df[col] == wanted].reset_index(drop=True)


def get_clone_names(seg_df: pd.DataFrame, normal: str | None = NORMAL_CLONE) -> tuple:
    """Split a profile's clones into the normal one and the tumor ones.

    Clone names come from the ``cn_*`` columns, ordered normal first then
    clone1, clone2, ... numerically, with any other name last. The normal clone
    is identified by name, so a profile that calls it something else needs
    ``normal`` set accordingly.

    Args:
        seg_df: Profile table with ``cn_<clone>`` columns.
        normal: Name of the normal clone, or None to treat every clone as tumor.

    Returns:
        (normal_clone, tumor_clones); ``normal_clone`` is None when the profile
        has no normal column.

    Raises:
        ValueError: If no ``cn_*`` column is present.
    """
    names = [c[len(CN_PREFIX) :] for c in seg_df.columns if c.startswith(CN_PREFIX)]
    if not names:
        raise ValueError(f"seg_df has no {CN_PREFIX}* columns")
    names.sort(key=_clone_sort_key)
    normal_clone = normal if normal in names else None
    return normal_clone, [n for n in names if n != normal_clone]


def get_clone_states(seg_df: pd.DataFrame, clones: list) -> np.ndarray:
    """Parse the ``cn_<clone>`` columns into integer copy-number states.

    Args:
        seg_df: Profile table.
        clones: Clone names to read, in the order they should be stacked.

    Returns:
        (n_bins, n_clones, 2) int array of (a, b) pairs.

    Raises:
        ValueError: If a column is missing or a value is not an "a|b" pair.
    """
    out = np.empty((len(seg_df), len(clones), 2), dtype=int)
    for k, clone in enumerate(clones):
        col = CN_PREFIX + clone
        if col not in seg_df.columns:
            raise ValueError(f"seg_df has no column {col!r}")
        parts = seg_df[col].astype(str).str.split("|", expand=True)
        if parts.shape[1] != 2:
            raise ValueError(f"{col!r} values are not 'a|b' pairs")
        try:
            out[:, k, :] = parts.astype(int).to_numpy()
        except ValueError as exc:
            raise ValueError(f"{col!r} has non-integer copy numbers") from exc
    return out


def get_clone_proportions(seg_df: pd.DataFrame, clones: list) -> dict:
    """Read per-clone proportions from the ``u_<clone>`` columns.

    Proportions are constant down a profile, so the first row is used.

    Args:
        seg_df: Profile table.
        clones: Clone names to read.

    Returns:
        {clone: proportion}, empty when the table carries no ``u_*`` columns.
    """
    cols = {c: U_PREFIX + c for c in clones if U_PREFIX + c in seg_df.columns}
    return {c: float(seg_df[u].iloc[0]) for c, u in cols.items()}


def get_pi_viol(seg_df: pd.DataFrame, col: str = PI_VIOL_COL) -> np.ndarray | None:
    """Read the optional per-bin pure-integer violation flags.

    Args:
        seg_df: Profile table.
        col: Column to read, matched case-insensitively.

    Returns:
        (n_bins,) bool array, or None when the column is absent.
    """
    for candidate in seg_df.columns:
        if candidate.lower() == col.lower():
            return seg_df[candidate].to_numpy(dtype=bool)
    return None


def has_mirror(states: np.ndarray) -> bool:
    """Test whether any bin is mirrored, for deciding on a legend swatch.

    Args:
        states: (n_bins, n_clones, 2) array from :func:`get_clone_states`.

    Returns:
        True if any bin has clones leaning opposite ways; see
        :func:`_mirrored_bins`.
    """
    return bool(_mirrored_bins(states).any())


# =============================================================================
# Integer CN profiles
# =============================================================================


def plot_cnv_profile(
    ax: plt.Axes,
    seg_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    ax_leg: plt.Axes | None = None,
    sample_id: str | None = None,
    normal: str | None = NORMAL_CLONE,
    clones: list | None = None,
    height: float = 1.0,
    title: str | None = None,
    ylabel: str | None = None,
    plot_chrname: bool = True,
    show_clone_name: bool = True,
    show_prop: bool = True,
    show_pi_viol: bool = True,
    clone_ploidies: dict | None = None,
    show_mirror: bool = True,
    clone_separators: bool = True,
    contain: bool = True,
    legend_kwargs: dict | None = None,
    rasterized: bool = True,
) -> None:
    """Draw a stacked integer copy-number profile, one row per clone.

    Each bin becomes a rectangle per clone, colored by its joint (a, b) state,
    with clone 1 as the top row. The selected sample's bins are matched onto
    ``genome_axis`` by coordinate, taking only those that lie entirely inside a segment.

    Args:
        ax: Axes to draw on.
        seg_df: Profile table with "#CHR", "START", "END", ``cn_<clone>``
            columns, optional ``u_<clone>`` proportions, and an optional
            ``PI_VIOL`` flag.
        genome_axis: Genome axis to draw on, built once from the reference and reusable
            across samples and panels. Bins are matched to it by coordinate.
        ax_leg: Axes for the legend, drawn in place. None omits it. The mirrored
            swatch appears only when a mirrored bin was actually hatched, so the
            legend never advertises a state the figure does not contain.
        sample_id: Sample to draw when the table holds several. None takes the
            first one present.
        normal: Name of the normal clone, excluded from the stack. None draws
            every clone, for a profile with no normal component.
        clones: Tumor clones to draw, in stacking order. None uses every
            non-normal clone.
        height: Total height of the clone stack, in axes fraction.
        title: Optional axes title.
        ylabel: Optional y-axis label.
        plot_chrname: Draw chromosome labels along the top.
        show_clone_name: Write "Clone N" rather than a bare "N".
        show_prop: Add each clone's proportion to its row label when the table
            carries ``u_<clone>`` columns.
        show_pi_viol: Draw the ``PI_VIOL`` overlay when the column is present:
            a line above each bin, red where True and green where False.
        clone_ploidies: Optional {clone: ploidy} map for the row labels.
        show_mirror: Hatch the flipped clones in mirrored bins, i.e. bins whose
            clones disagree on which allele is amplified. A bin where every clone
            leans the same way is not hatched, since the palette gives (a, b) and
            (b, a) one color and a uniform lean carries no information.
        clone_separators: Draw horizontal rules between adjacent clone rows.
        contain: Require a row to lie entirely inside a region to be drawn,
            which is exact when the profile was called against the same regions
            being plotted. False clips overlapping rows instead. Leaving it True
            means a reference mismatch shows up as a loud unmapped-row warning
            rather than as a quietly clipped, plausible-looking figure.
        legend_kwargs: Extra styling forwarded to the legend, ignored without
            ``ax_leg``.
        rasterized: Rasterize the rectangles, keeping vector output small.

    Raises:
        ValueError: If ``seg_df`` is empty or has no clone left to draw.
    """
    if seg_df.empty:
        raise ValueError("seg_df is empty")
    if clones is None:
        _, clones = get_clone_names(seg_df, normal=normal)
    if not clones:
        raise ValueError("no clone to draw")

    seg_df = select_sample(seg_df, sample_id)
    coords = genome_axis.build_coordinates(seg_df, contain=contain)
    states = get_clone_states(seg_df, clones)
    props = get_clone_proportions(seg_df, clones) if show_prop else None
    pi_viol = get_pi_viol(seg_df) if show_pi_viol else None

    state_style, _ = get_cn_colors()
    mirrored = _mirrored_bins(states) if show_mirror else None
    num_clones = len(clones)
    h = height / num_clones

    for row, x0, w in _iter_bins(coords):
        for k in range(num_clones):
            cna, cnb = states[row, num_clones - k - 1]
            y0 = k * h
            ax.add_patch(
                Rectangle(
                    (x0, y0),
                    w,
                    h,
                    facecolor=state_style.get((cna, cnb), state_style["default"]),
                    edgecolor="none",
                    linewidth=0,
                    antialiased=False,
                    rasterized=rasterized,
                    transform=ax.get_xaxis_transform(),
                )
            )
            if mirrored is not None and mirrored[row] and cnb > cna:
                ax.add_patch(
                    Rectangle(
                        (x0, y0),
                        w,
                        h,
                        facecolor="none",
                        edgecolor="black",
                        hatch="/",
                        linewidth=0,
                        rasterized=rasterized,
                        transform=ax.get_xaxis_transform(),
                    )
                )
        if pi_viol is not None:
            ax.plot(
                [x0, x0 + w],
                [height, height],
                color="#d62728" if bool(pi_viol[row]) else "#2ca02c",
                linewidth=3,
                solid_capstyle="butt",
                transform=ax.get_xaxis_transform(),
            )

    draw_segment_boundaries(ax, genome_axis, linewidth=0.5, alpha=1.0)
    draw_chr_boundaries(ax, genome_axis, ymax=1.15)
    if clone_separators:
        for k in range(1, num_clones):
            ax.hlines(
                k * h,
                xmin=genome_axis.ch_coords[0],
                xmax=genome_axis.chr_end,
                transform=ax.get_xaxis_transform(),
                linewidth=0.8,
                colors="black",
            )

    _finish_clone_axis(
        ax,
        genome_axis,
        clones,
        height,
        plot_chrname,
        clone_ploidies,
        props,
        show_clone_name,
        ylabel,
        title,
        label_fontweight="bold",
        tick_pad=4,
    )
    if ax_leg is not None:
        _plot_cnv_legend(
            ax_leg,
            has_mirror=mirrored is not None and bool(mirrored.any()),
            **(legend_kwargs or {}),
        )


def _plot_cnv_legend(
    ax: plt.Axes,
    has_mirror: bool = True,
    pair_w: float = 2.0,
    pair_h: float = 0.6,
    fontsize: int = 10,
) -> None:
    """Draw the integer copy-number legend, grouped by total copy number.

    Mirrored pairs share a swatch since (1, 2) and (2, 1) share a color; the "/"
    hatch distinguishes them. A trailing box covers total CN above 7.

    Args:
        ax: Axes to draw on. Its own frame is turned off.
        has_mirror: Append the mirrored-CNA swatch. :func:`plot_cnv_profile`
            passes True only when it actually hatched a mirrored bin.
        pair_w: Width of each swatch.
        pair_h: Height of each swatch.
        fontsize: Swatch label font size.

    """
    state_style, tcn_states = get_cn_colors()
    ax.axis("off")

    gap_pairs = 0.0
    gap_groups = 0.75 * pair_w
    mirror_gap = 1.5 * pair_w
    leg_x = 0.0

    for _total, states in sorted(tcn_states.items()):
        uniq_pairs = sorted({tuple(sorted(s, reverse=True)) for s in states})
        group_x0 = leg_x
        for i, pair in enumerate(uniq_pairs):
            x0 = group_x0 + i * (pair_w + gap_pairs)
            ax.add_patch(
                Rectangle(
                    (x0, 0.0),
                    pair_w,
                    pair_h,
                    facecolor=state_style[pair],
                    edgecolor="black",
                )
            )
            ax.text(
                x0 + pair_w / 2.0,
                -0.2,
                f"{pair}",
                ha="center",
                va="top",
                fontsize=fontsize,
            )
        group_w = len(uniq_pairs) * pair_w + (len(uniq_pairs) - 1) * gap_pairs
        leg_x = group_x0 + group_w + gap_groups

    ax.add_patch(
        Rectangle(
            (leg_x, 0.0),
            pair_w,
            pair_h,
            facecolor=state_style["default"],
            edgecolor="black",
        )
    )
    ax.text(leg_x + pair_w / 2.0, -0.2, ">7", ha="center", va="top", fontsize=fontsize)
    leg_x += pair_w + gap_groups

    if has_mirror:
        leg_x = _draw_mirror_swatch(
            ax, leg_x + mirror_gap, pair_w=pair_w, pair_h=pair_h, fontsize=fontsize
        )
        leg_x += mirror_gap

    ax.text(
        -0.5,
        pair_h / 2.0,
        "CNA",
        fontsize=fontsize + 2,
        fontweight="bold",
        ha="right",
        va="center",
    )
    ax.set_xlim(-2.0, leg_x)
    ax.set_ylim(-0.8, pair_h + 0.8)
    ax.set_aspect("auto")


# =============================================================================
# Legacy
# =============================================================================
# Allele-specific panels, kept for the existing HATCHet figures that use them.
# Not exported, and not covered by the supported interface.


def plot_ascn_profile(
    ax: plt.Axes,
    seg_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    ax_leg: plt.Axes | None = None,
    sample_id: str | None = None,
    normal: str | None = NORMAL_CLONE,
    clones: list | None = None,
    height: float = 1.0,
    title: str | None = None,
    ylabel: str | None = None,
    plot_chrname: bool = True,
    show_clone_name: bool = True,
    show_prop: bool = True,
    clone_ploidies: dict | None = None,
    show_mirror: bool = True,
    outline_alleles: bool = True,
    contain: bool = True,
    legend_kwargs: dict | None = None,
    rasterized: bool = True,
) -> None:
    """Draw an allele-specific copy-number profile, two sub-bars per clone.

    Legacy. Retained for the HATCHet figures that already use it; not part of the
    supported interface. It shows per-allele copy number, which
    :func:`plot_cnv_profile` does not - its joint palette gives (3, 0) and (2, 1)
    distinct colors but never the A and B values separately.

    Each clone slot splits into a lower B bar and an upper A bar, coloured by
    its own integer copy number. Non-zero values render at alpha 0.5, so CN 1
    reads as mid gray against the white CN 0. Bins are matched onto ``genome_axis`` by
    coordinate, as in :func:`plot_cnv_profile`.

    Args:
        ax: Axes to draw on.
        seg_df: Profile table with "#CHR", "START", "END", ``cn_<clone>``
            columns, and optional ``u_<clone>`` proportions.
        genome_axis: Genome axis to draw on, built once from the reference and reusable
            across samples and panels. Bins are matched to it by coordinate.
        ax_leg: Axes for the legend, drawn in place. None omits it. The mirrored
            swatch appears only when a mirrored bin was actually hatched, so the
            legend never advertises a state the figure does not contain.
        sample_id: Sample to draw when the table holds several. None takes the
            first one present.
        normal: Name of the normal clone, excluded from the stack. None draws
            every clone.
        clones: Tumor clones to draw, in stacking order. None uses every
            non-normal clone.
        height: Total height of the clone stack, in axes fraction.
        title: Optional axes title.
        ylabel: Optional y-axis label.
        plot_chrname: Draw chromosome labels along the top.
        show_clone_name: Write "Clone N" rather than a bare "N".
        show_prop: Add each clone's proportion to its row label.
        clone_ploidies: Optional {clone: ploidy} map for the row labels.
        show_mirror: Mark mirrored bins - those whose clones disagree on which
            allele is amplified - with chevrons pointing toward the larger
            allele. Same rule as :func:`plot_cnv_profile`.
        outline_alleles: Outline each allele row across the genome.
        contain: Require a row to lie entirely inside a region to be drawn.
            Same rule as :func:`plot_cnv_profile`.
        legend_kwargs: Extra styling forwarded to the legend, ignored without
            ``ax_leg``.
        rasterized: Rasterize the rectangles.

    Raises:
        ValueError: If ``seg_df`` is empty or has no clone left to draw.
    """
    if seg_df.empty:
        raise ValueError("seg_df is empty")
    if clones is None:
        _, clones = get_clone_names(seg_df, normal=normal)
    if not clones:
        raise ValueError("no clone to draw")

    seg_df = select_sample(seg_df, sample_id)
    coords = genome_axis.build_coordinates(seg_df, contain=contain)
    states = get_clone_states(seg_df, clones)
    props = get_clone_proportions(seg_df, clones) if show_prop else None

    state_style, _ = get_ascn_colors()
    mirrored = _mirrored_bins(states) if show_mirror else None
    num_clones = len(clones)
    h = height / num_clones
    clone_gap = 0.10 * h
    h_pair = h - clone_gap
    h_sub = h_pair / 2
    y_gap = clone_gap / 2

    for row, x0, w in _iter_bins(coords):
        bin_states = states[row]
        for k in range(num_clones):
            cna, cnb = bin_states[num_clones - k - 1]
            y_b = k * h + y_gap
            y_a = y_b + h_sub
            for cn, y in ((cnb, y_b), (cna, y_a)):
                ax.add_patch(
                    Rectangle(
                        (x0, y),
                        w,
                        h_sub,
                        facecolor=state_style.get(cn, state_style["default"]),
                        edgecolor="none",
                        linewidth=0,
                        alpha=1.0 if cn == 0 else 0.5,
                        rasterized=rasterized,
                        transform=ax.get_xaxis_transform(),
                    )
                )
            if mirrored is not None and mirrored[row] and cna != cnb:
                direction = 1 if cna > cnb else -1
                unit = w * 0.12
                gap = w * 0.04
                total_w = 2 * unit + gap
                _draw_chevrons(
                    ax,
                    x0 + (w - total_w) / 2.0,
                    unit,
                    gap,
                    y_b + 2 * h_sub,
                    y_b,
                    direction,
                    transform=ax.get_xaxis_transform(),
                )

    # per-clone dashed lines so a removed region does not cross the clone gaps
    for x in genome_axis.seg_coords:
        for k in range(num_clones):
            ax.vlines(
                x,
                ymin=k * h + y_gap,
                ymax=k * h + y_gap + h_pair,
                transform=ax.get_xaxis_transform(),
                linewidth=0.5,
                colors="black",
                linestyles="dashed",
            )
    draw_chr_boundaries(ax, genome_axis, ymax=1.15)

    if outline_alleles:
        span = genome_axis.chr_end - genome_axis.ch_coords[0]
        for k in range(num_clones):
            y_b_k = k * h + y_gap
            for y0 in (y_b_k, y_b_k + h_sub):
                ax.add_patch(
                    Rectangle(
                        (genome_axis.ch_coords[0], y0),
                        span,
                        h_sub,
                        facecolor="none",
                        edgecolor="black",
                        linewidth=0.5,
                        transform=ax.get_xaxis_transform(),
                    )
                )

    minor_positions = []
    minor_labels = []
    for k in range(num_clones):
        minor_positions += [k * h + y_gap + h_sub * 0.5, k * h + y_gap + h_sub * 1.5]
        minor_labels += ["B", "A"]
    ax.set_yticks(minor_positions, minor=True)
    ax.set_yticklabels(minor_labels, minor=True, fontsize=6)
    ax.tick_params(axis="y", which="minor", left=False, right=False, pad=2)

    _finish_clone_axis(
        ax,
        genome_axis,
        clones,
        height,
        plot_chrname,
        clone_ploidies,
        props,
        show_clone_name,
        ylabel,
        title,
        label_fontweight=None,
        tick_pad=20,
    )
    if ax_leg is not None:
        _plot_ascn_legend(
            ax_leg,
            show_mirror=mirrored is not None and bool(mirrored.any()),
            **(legend_kwargs or {}),
        )


def _plot_ascn_legend(
    ax: plt.Axes,
    box_w: float = 1.2,
    box_h: float = 0.4,
    tick_len: float = 0.08,
    label_fontsize: int = 12,
    show_mirror: bool = True,
) -> None:
    """Draw the allele copy-number legend as a horizontal ramp.

    Legacy, paired with :func:`plot_ascn_profile`.

    One box per copy number 0 to 6 plus "7+", at the same alpha the profiles use.

    Args:
        ax: Axes to draw on. Its own frame is turned off.
        box_w: Width of each box.
        box_h: Height of each box.
        tick_len: Length of the tick below each box.
        label_fontsize: Font size for labels and the title.
        show_mirror: Append the mirrored chevron swatch.
            :func:`plot_ascn_profile` passes True only when it actually drew
            chevrons.

    """
    state_style, tcn_states = get_ascn_colors()
    boxes = list(tcn_states) + ["7+"]
    ax.axis("off")

    for i, label in enumerate(boxes):
        color = state_style["default"] if label == "7+" else state_style[label]
        ax.add_patch(
            Rectangle(
                (i * box_w, 0.0),
                box_w,
                box_h,
                facecolor=color,
                edgecolor="black",
                alpha=1.0 if label == 0 else 0.5,
            )
        )
        xc = i * box_w + box_w / 2.0
        ax.plot([xc, xc], [-tick_len, 0.0], color="black", linewidth=0.8)
        ax.text(
            xc,
            -tick_len - 0.04,
            str(label),
            ha="center",
            va="top",
            fontsize=label_fontsize,
            fontweight="bold",
        )

    total_w = len(boxes) * box_w
    ax.text(
        -0.3,
        box_h / 2.0,
        "Allele copy number",
        fontsize=label_fontsize,
        fontweight="bold",
        ha="right",
        va="center",
    )

    right = total_w
    if show_mirror:
        swatch_w = box_w * 0.7
        chev_box_x = total_w + 1.0
        ax.add_patch(
            Rectangle(
                (chev_box_x, 0.0),
                swatch_w,
                box_h,
                facecolor="white",
                edgecolor="black",
            )
        )
        unit = swatch_w * 0.30
        gap = swatch_w * 0.10
        chev_total_w = 2 * unit + gap
        _draw_chevrons(
            ax,
            chev_box_x + (swatch_w - chev_total_w) / 2.0,
            unit,
            gap,
            box_h,
            0.0,
            direction=1,
        )
        ax.text(
            chev_box_x + swatch_w / 2.0,
            -tick_len - 0.04,
            "mirrored",
            ha="center",
            va="top",
            fontsize=label_fontsize,
            fontweight="bold",
        )
        right = chev_box_x + swatch_w + 0.5

    ax.set_xlim(-2.0, right)
    ax.set_ylim(-0.5, box_h + 0.2)
    ax.set_aspect("auto")
