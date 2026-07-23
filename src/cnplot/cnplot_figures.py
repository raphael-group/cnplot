"""High-level figure builders that compose the per-axes primitives.

The rest of the package draws on an ``ax`` you own and returns None. The builders
here own the layout instead: they create a figure, place a ``GridSpec``, call the
primitives (``plot_scatter_1d``, ``plot_cnv_profile``, ...), and hand back a
finished figure. They sit at the top of the dependency graph and import across the
domain modules, which is why they live apart from any one of them.

Built here: ``plot_scatter_1d_multisample`` (stacked 1D panels + shared profile)
and ``plot_heatmap_cnp`` (heatmap + CN profile + legend page). Planned: the
multi-solution CNP panels (``plot_cnp_panel`` / ``pool_cnp``, M9).
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D

from .cnplot_1d import plot_scatter_1d
from .cnplot_genome_axis import GenomeAxis
from .cnplot_heatmap import plot_heatmap
from .cnplot_intcnp import (
    get_clone_names,
    get_clone_proportions,
    plot_cnv_profile,
    select_sample,
)
from .cnplot_utils import (
    NORMAL_CLONE,
    SAMPLE_COL,
    format_clone_name,
    resolve_marker_size,
)

logger = logging.getLogger(__name__)

__all__ = [
    "make_row_spec",
    "plot_heatmap_cnp",
    "plot_scatter_1d_multisample",
]


# =============================================================================
# Helpers
# =============================================================================


def _draw_prop_legend(
    ax: plt.Axes,
    seg_df: pd.DataFrame,
    group,
    group_col: str,
    normal: str | None,
    display_min_clone_prop: float | None,
) -> None:
    """Label one group's clone proportions beside its row.
    Args:
        ax: Axes to attach the legend to; it is placed just outside the right edge.
        seg_df: Profile table carrying ``u_<clone>`` columns per group.
        group: Group to read. None takes the first.
        group_col: Column identifying a group.
        normal: Normal clone name, always listed however small.
        display_min_clone_prop: Hide tumor clones below this proportion.
    """
    try:
        sub = select_sample(seg_df, group, col=group_col)
    except ValueError:
        logger.info("no %r row in seg_df; skipping its proportion legend", group)
        return
    normal_clone, tumor = get_clone_names(sub, normal=normal)
    names = ([normal_clone] if normal_clone else []) + tumor
    props = get_clone_proportions(sub, names)
    handles = [
        Line2D([0], [0], alpha=0, label=f"{format_clone_name(n)}: {p:.3f}")
        for n, p in props.items()
        if n == normal_clone
        or display_min_clone_prop is None
        or p >= display_min_clone_prop
    ]
    if not handles:
        return
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        fontsize="small",
        fancybox=True,
        framealpha=0.7,
        handlelength=0,
        handletextpad=0,
    )


def make_row_spec(
    ycol: str,
    ylabel: str | None = None,
    ylim: tuple | None = None,
    href: float | None = None,
    reverse_y: bool = False,
    hue_col: str | None = None,
) -> dict:
    """Describe one row of a stacked 1D panel.

    A panel row differs from its neighbours only in which column it draws and
    how that column is scaled, so rows are data rather than code; the panel
    supplies the shared styling and decides which row carries the chromosome
    labels.

    Args:
        ycol: Observed value column for this row.
        ylabel: Y-axis label; defaults to ``ycol``.
        ylim: Fixed (lo, hi), or None to leave matplotlib's autoscale alone.
        href: Y value for a dotted grey reference line, or None.
        reverse_y: Invert the y-axis, for a mirrored allele row.
        hue_col: Column to color this row by, overriding the panel default.

    Returns:
        A row spec for :func:`plot_scatter_1d_multisample`.
    """
    return {
        "ycol": ycol,
        "ylabel": ycol if ylabel is None else ylabel,
        "ylim": ylim,
        "href": href,
        "reverse_y": reverse_y,
        "hue_col": hue_col,
    }


# =============================================================================
# 1D multi-sample panel
# =============================================================================


def plot_scatter_1d_multisample(
    obs_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    row_specs: list,
    groups: list | None = None,
    group_col: str = SAMPLE_COL,
    expected_df: pd.DataFrame | None = None,
    hue=None,
    palette=None,
    alphas=None,
    keep_col: str | None = None,
    filtered_color="red",
    seg_df: pd.DataFrame | None = None,
    show_props: bool = True,
    normal: str | None = NORMAL_CLONE,
    display_min_clone_prop: float | None = None,
    titles: dict | None = None,
    row_width: float = 20.0,
    row_height: float = 2.0,
    markersize: float | None = None,
    intra_group_hspace: float = 0.25,
    inter_group_hspace: float = 0.15,
    profile_hspace: float = 0.2,
    rasterized: bool = True,
    **profile_kwargs,
) -> plt.Figure:
    """Stack per-group rows over a shared copy-number profile.

    One block of ``row_specs`` rows per group, then the clonal CN profile and its
    legend once at the bottom. A group is a bulk sample or a cell group - the
    layout does not distinguish them, and neither host repo's version of this
    figure ever did.

    Only the bottom row carries chromosome labels, and the marker size is
    resolved once for the whole figure so no two rows disagree.

    Args:
        obs_df: Bins for every group, with "#CHR", "START", "END", the value
            columns named by ``row_specs``, and ``group_col``.
        genome_axis: Axis shared by every row.
        row_specs: Rows per group, from :func:`make_row_spec`.
        groups: Groups to draw, top to bottom. None takes them in order of
            appearance.
        group_col: Column identifying a group. Absent means one unnamed group.
        expected_df: Expected values for every group, with one
            ``exp_<ycol>_<group>`` column each.
        hue: Column name in ``obs_df`` to color points by, or an (n,) array
            aligned to it.
        palette: {label: color} for ``hue``.
        alphas: Column name or (n,) array of per-point alpha.
        keep_col: Bool column in ``obs_df`` splitting kept from filtered rows;
            see :func:`~cnplot.cnplot_1d.plot_scatter_1d`.
        filtered_color: Color for the rows ``keep_col`` marks False.
        seg_df: Copy-number profile for the bottom strip. None omits the strip
            and its legend. Also supplies the per-group proportions.
        show_props: Label each group's clone proportions beside its rows rather
            than on the profile. The profile is shared across groups while the
            proportions are not, so printing them there would attribute one
            group's mixture to all of them. Pass ``show_prop=True`` through
            ``profile_kwargs`` to override.
        normal: Normal clone name, always listed in the proportion legend.
        display_min_clone_prop: Hide tumor clones below this proportion from the
            proportion legend.
        titles: {group: title}; a group without one gets its own name.
        row_width: Figure width in inches.
        row_height: Height of one row in inches.
        markersize: Marker size; None resolves one from the busiest group.
        intra_group_hspace: Space between rows of one group.
        inter_group_hspace: Space between groups, as a fraction of a group's
            height; a group is several rows tall, so keep this small.
        profile_hspace: Space between the last row and the profile, which needs
            room for the chromosome labels it draws along its top.
        rasterized: Rasterize the point collections.
        **profile_kwargs: Passed to
            :func:`~cnplot.cnplot_intcnp.plot_cnv_profile`.

    Returns:
        The figure. The caller saves and closes it.

    Raises:
        ValueError: If ``row_specs`` is empty or a named group is absent.
    """
    if not row_specs:
        raise ValueError("row_specs is empty")
    profile_kwargs.setdefault("show_prop", False)
    if group_col in obs_df.columns:
        present = list(dict.fromkeys(obs_df[group_col]))
        groups = present if groups is None else list(groups)
        missing = [g for g in groups if g not in present]
        if missing:
            raise ValueError(f"groups absent from {group_col}: {missing}")
    else:
        groups = [None]

    if isinstance(hue, str):
        hue = obs_df[hue].to_numpy()
    if isinstance(alphas, str):
        alphas = obs_df[alphas].to_numpy()
    hue = None if hue is None else np.asarray(hue)
    alphas = None if alphas is None else np.asarray(alphas)

    masks = {
        g: (
            np.ones(len(obs_df), dtype=bool)
            if g is None
            else (obs_df[group_col] == g).to_numpy()
        )
        for g in groups
    }
    if markersize is None:
        markersize = resolve_marker_size(max(int(m.sum()) for m in masks.values()))

    n_rows = len(row_specs)
    n_groups = len(groups)
    has_profile = seg_df is not None
    rows_h = row_height * n_rows * n_groups
    profile_h = 2.0 if has_profile else 0.0
    fig_h = rows_h + profile_h + 0.5
    fig = plt.figure(figsize=(row_width, fig_h))

    # Nested so the three gaps stay independent: rows within a group, group to
    # group, and rows to profile. One flat grid would tie the last two together.
    if has_profile:
        outer = GridSpec(
            2,
            1,
            figure=fig,
            height_ratios=[rows_h, profile_h],
            hspace=profile_hspace,
            top=0.97,
        )
        groups_gs = GridSpecFromSubplotSpec(
            n_groups, 1, subplot_spec=outer[0], hspace=inter_group_hspace
        )
    else:
        groups_gs = GridSpec(
            n_groups, 1, figure=fig, hspace=inter_group_hspace, top=0.97
        )

    for gi, group in enumerate(groups):
        inner = GridSpecFromSubplotSpec(
            n_rows,
            1,
            subplot_spec=groups_gs[gi],
            hspace=intra_group_hspace,
        )
        mask = masks[group]
        sub = obs_df[mask]
        last_group = gi == n_groups - 1
        for ri, spec in enumerate(row_specs):
            ax = fig.add_subplot(inner[ri])
            last_row = ri == n_rows - 1
            if spec["hue_col"] is not None:
                row_hue = sub[spec["hue_col"]].to_numpy()
            elif hue is not None:
                row_hue = hue[mask]
            else:
                row_hue = None
            plot_scatter_1d(
                ax,
                sub,
                genome_axis,
                spec["ycol"],
                expected_df=expected_df,
                group=group,
                hue=row_hue,
                palette=palette,
                alphas=None if alphas is None else alphas[mask],
                keep_col=keep_col,
                filtered_color=filtered_color,
                show_legend=ri == 0,
                markersize=markersize,
                ylim=spec["ylim"],
                ylabel=spec["ylabel"],
                title=(
                    (titles or {}).get(group, group if group is not None else None)
                    if ri == 0
                    else None
                ),
                href=spec["href"],
                reverse_y=spec["reverse_y"],
                plot_chrname=last_row and last_group and not has_profile,
                rasterized=rasterized,
            )
            if last_row and show_props and seg_df is not None:
                _draw_prop_legend(
                    ax, seg_df, group, group_col, normal, display_min_clone_prop
                )

    if has_profile:
        # outer is [all groups, profile]; the groups are nested inside outer[0]
        inner_bot = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[1], height_ratios=[3, 1], hspace=0.15
        )
        ax_profile = fig.add_subplot(inner_bot[0])
        ax_legend = fig.add_subplot(inner_bot[1])
        plot_cnv_profile(
            ax_profile, seg_df, genome_axis, ax_leg=ax_legend, **profile_kwargs
        )
    return fig


# =============================================================================
# Single-cell heatmap page
# =============================================================================


def plot_heatmap_cnp(
    matrix: np.ndarray,
    coords_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    seg_df: pd.DataFrame,
    sample_id: str | None = None,
    title: str | None = None,
    figsize: tuple = (15, 8),
    height_ratios: tuple = (10, 2, 1.5),
    profile_hspace: float = 0.35,
    top: float = 0.9,
    profile_kwargs: dict | None = None,
    **heatmap_kwargs,
) -> plt.Figure:
    """Compose the single-cell heatmap page: mesh, CN profile, and legend.

    Three stacked rows - the value heatmap (with its colorbar, side strips, and
    legends), the integer-CN profile, and the profile legend. The chromosome
    names sit on top of the profile, between it and the heatmap above, and label
    both; the heatmap itself is drawn without them.

    Everything that shapes the mesh - the row-reduced ``matrix``, its ``cmap`` /
    ``norm``, and the ``strips`` - is a ``plot_heatmap`` argument passed straight
    through; building those inputs stays with the caller.

    Args:
        matrix: (n_rows, n_bins) values, columns aligned to ``coords_df`` rows.
        coords_df: The bins, with "#CHR", "START", "END".
        genome_axis: Axis shared by the heatmap and the profile.
        seg_df: seg.ucn profile drawn under the heatmap.
        sample_id: Sample to draw from ``seg_df`` when it holds several.
        title: Figure title, bold, above the heatmap. None omits it.
        figsize: Figure size in inches.
        height_ratios: Row heights for (heatmap, profile, legend).
        profile_hspace: Space between the heatmap and the profile below it; holds
            the chromosome labels drawn on top of the profile.
        top: Top of the axes block in figure fraction; the space above holds the
            title.
        profile_kwargs: Extra arguments for
            :func:`~cnplot.cnplot_intcnp.plot_cnv_profile`. ``plot_chrname``
            defaults on so the profile carries the chromosome labels.
        **heatmap_kwargs: Passed to :func:`~cnplot.cnplot_heatmap.plot_heatmap`;
            ``show_colorbar`` and hidden block labels default on for the page, and
            ``plot_chrname`` is forced off so only the profile labels chromosomes.

    Returns:
        The figure. The caller saves and closes it.
    """
    heatmap_kwargs.setdefault("show_colorbar", True)
    heatmap_kwargs.setdefault("show_block_labels", False)
    heatmap_kwargs["plot_chrname"] = False
    profile_kwargs = dict(profile_kwargs or {})
    profile_kwargs.setdefault("plot_chrname", True)

    fig, axes = plt.subplots(
        3, 1, figsize=figsize, gridspec_kw={"height_ratios": list(height_ratios)}
    )
    fig.subplots_adjust(top=top, hspace=profile_hspace)

    plot_heatmap(axes[0], matrix, coords_df, genome_axis, **heatmap_kwargs)
    plot_cnv_profile(
        axes[1],
        seg_df,
        genome_axis,
        ax_leg=axes[2],
        sample_id=sample_id,
        **profile_kwargs,
    )
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)
    return fig
