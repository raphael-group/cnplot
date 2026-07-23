"""2D scatter: one value against another, with expected copy-number landmarks.

Like the 1D module, both plotted quantities arrive precomputed in columns named
by ``xcol`` and ``ycol``, and the expected coordinates arrive in their own table.

The landmark table uses the seg.ucn layout the profile plots already read, so the
``cn_<clone>`` columns are present and one clone subset switches between the two
conventions the host repos use: a bulk caller passes every clone and gets one
landmark per distinct joint state, a single-cell caller passes the one clone
behind the group and gets one per distinct (a, b). The ``(a, b)`` label needs
those columns anyway, so a bare table of coordinates would not do.
"""

import contextlib
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.lines import Line2D

from .cnplot_intcnp import (
    get_clone_names,
    get_clone_proportions,
    get_clone_states,
    select_sample,
)
from .cnplot_utils import (
    EXP_PREFIX,
    NORMAL_CLONE,
    SAMPLE_COL,
    resolve_colors,
    resolve_marker_size,
)

logging.getLogger("adjustText").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

__all__ = [
    "get_landmarks",
    "plot_scatter_2d",
]


# =============================================================================
# Helpers
# =============================================================================


def _visible_clones(
    clones: list,
    props: dict,
    display_min_clone_prop: float | None,
) -> list:
    """Keep the clones big enough to speak for the sample.

    Whether a state reads as clonal depends on which clones are actually present
    in *this* sample, so a clone below the threshold is not allowed to make a
    state look subclonal. Falls back to every clone when the threshold would hide
    them all, rather than producing an empty label.

    Args:
        clones: Clone names being labelled.
        props: {clone: proportion} for this sample; a clone missing from it
            counts as visible.
        display_min_clone_prop: Smallest proportion that counts, or None to keep
            every clone.

    Returns:
        The clones that count toward the label, in the given order.
    """
    if display_min_clone_prop is None:
        return list(clones)
    visible = [c for c in clones if props.get(c, 1.0) >= display_min_clone_prop]
    return visible or list(clones)


def _state_label(states: np.ndarray) -> tuple:
    """Label one landmark from the copy-number states behind it.

    Args:
        states: (n_clones, 2) array of (a, b) pairs, already restricted to the
            visible clones.

    Returns:
        (text, clonal): the distinct states joined by commas in clone order, and
        whether every clone agreed on one state.
    """
    seen = list(dict.fromkeys(f"({int(a)},{int(b)})" for a, b in states))
    return ",".join(seen), len(seen) == 1


def _annotate_landmarks(
    ax: plt.Axes,
    landmarks: list,
    markersize: float = 30,
    linewidth: float = 1.0,
    fontsize: int = 9,
    color="black",
) -> None:
    """Draw hollow markers with repelled ``(a, b)`` labels.

    Labels are pushed off each other by ``adjust_text`` and tied back with a thin
    connector, since landmarks of neighbouring states often overlap.

    Args:
        ax: Axes to draw on.
        landmarks: Dicts with "x", "y", "label" and optional "clonal" and
            "balanced". A clonal label is drawn bold and a balanced state gets a
            square rather than a circle. Every landmark gets a marker; one whose
            "label" is empty gets no text.
        markersize: Marker size in points squared.
        linewidth: Marker edge width.
        fontsize: Label font size.
        color: Marker edge and connector color.
    """
    if not landmarks:
        return
    xs = np.array([float(p["x"]) for p in landmarks])
    ys = np.array([float(p["y"]) for p in landmarks])
    balanced = np.array([bool(p.get("balanced", False)) for p in landmarks])

    for mask, marker in ((~balanced, "o"), (balanced, "s")):
        if mask.any():
            ax.scatter(
                xs[mask],
                ys[mask],
                facecolors="none",
                edgecolors=color,
                s=markersize,
                linewidth=linewidth,
                marker=marker,
                zorder=10,
            )
    labelled = [p for p in landmarks if p.get("label")]
    texts = [
        ax.text(
            p["x"],
            p["y"],
            p["label"],
            fontsize=fontsize,
            fontweight="bold" if p.get("clonal") else "normal",
            zorder=11,
        )
        for p in labelled
    ]
    if not texts:
        return
    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
        adjust_text(
            texts,
            x=np.array([float(p["x"]) for p in labelled]),
            y=np.array([float(p["y"]) for p in labelled]),
            ax=ax,
            arrowprops=dict(arrowstyle="-", color=color, lw=0.5),
        )


def get_landmarks(
    expected_df: pd.DataFrame,
    xcol: str,
    ycol: str,
    group=None,
    clones: list | None = None,
    normal: str | None = NORMAL_CLONE,
    group_col: str = SAMPLE_COL,
    display_min_clone_prop: float | None = None,
) -> list:
    """Collect the expected (x, y) landmarks for one group.

    Rows sharing a copy-number state collapse to one landmark, which is what
    both host repos do: many segments carry the same state and would otherwise
    stack identical markers on one point.

    Args:
        expected_df: seg.ucn-layout table with ``cn_<clone>`` columns, the
            precomputed ``exp_<xcol>`` and ``exp_<ycol>`` coordinates, and
            optional ``is_filtered`` / ``is_balanced`` flags.
        xcol: Observed x column, naming ``exp_<xcol>``.
        ycol: Observed y column, naming ``exp_<ycol>``.
        group: Group to take when the table holds several. None takes the first.
        clones: Clones whose states define a landmark. None uses every clone in
            the table, giving one landmark per distinct joint state.
        normal: Normal clone to exclude when ``clones`` is None, since a
            landmark labels the tumor states. None keeps every clone.
        group_col: Column identifying a group.
        display_min_clone_prop: Clones below this proportion do not count toward
            the label or the clonal test. Proportions are per group, so the same
            state can read clonal in one sample and subclonal in another. None
            counts every clone.

    Returns:
        Dicts with "x", "y", "label", "clonal", "states", "visible_clones" and
        "balanced", in first-seen order. Rows flagged ``is_filtered`` are left
        out, and a repeated label is blanked so only the first carries text.

    Raises:
        ValueError: If an expected coordinate column is missing.
    """
    df = select_sample(expected_df, group, col=group_col)
    x_exp, y_exp = f"{EXP_PREFIX}{xcol}", f"{EXP_PREFIX}{ycol}"
    for col in (x_exp, y_exp):
        if col not in df.columns:
            raise ValueError(f"expected_df has no column {col!r}")

    if clones is None:
        _, clones = get_clone_names(df, normal=normal)
    states = get_clone_states(df, clones)
    visible = _visible_clones(
        clones, get_clone_proportions(df, clones), display_min_clone_prop
    )
    vis_idx = [clones.index(c) for c in visible]

    if "is_filtered" in df.columns:
        keep = ~df["is_filtered"].to_numpy(dtype=bool)
    else:
        keep = np.ones(len(df), dtype=bool)
    balanced = (
        df["is_balanced"].to_numpy(dtype=bool)
        if "is_balanced" in df.columns
        else np.zeros(len(df), dtype=bool)
    )
    xs = df[x_exp].to_numpy(dtype=float)
    ys = df[y_exp].to_numpy(dtype=float)

    landmarks = []
    seen_states = set()
    seen_labels = set()
    for i in range(len(df)):
        if not keep[i] or not (np.isfinite(xs[i]) and np.isfinite(ys[i])):
            continue
        key = states[i].tobytes()
        if key in seen_states:
            continue
        seen_states.add(key)
        label, clonal = _state_label(states[i][vis_idx])
        # Distinct states can share a label once hidden clones drop out. Both
        # keep their marker; only the first carries the text.
        if label in seen_labels:
            label = ""
        else:
            seen_labels.add(label)
        landmarks.append(
            {
                "x": xs[i],
                "y": ys[i],
                "label": label,
                "clonal": clonal,
                "states": states[i],
                "visible_clones": visible,
                "balanced": bool(balanced[i]),
            }
        )
    return landmarks


# =============================================================================
# 2D scatter
# =============================================================================


def plot_scatter_2d(
    obs_df: pd.DataFrame,
    xcol: str,
    ycol: str,
    expected_df: pd.DataFrame | None = None,
    group=None,
    group_col: str = SAMPLE_COL,
    clones: list | None = None,
    normal: str | None = NORMAL_CLONE,
    hue=None,
    palette=None,
    colors=None,
    alphas=None,
    markersize: float | None = None,
    landmark_size: float = 30,
    xlim: tuple | None = None,
    ylim: tuple | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    refline_x: float | None = 0.5,
    refline_y: float | None = None,
    show_marginals: bool = True,
    show_props: bool = True,
    display_min_clone_prop: float | None = None,
    rasterized: bool = True,
) -> sns.JointGrid:
    """Scatter two value columns against each other with marginal densities.

    A ``JointGrid`` carries the scatter plus one unfilled KDE curve per hue level
    on each margin. Expected copy-number landmarks are drawn as hollow markers
    with ``(a, b)`` labels, squares where a state is flagged balanced.

    Args:
        obs_df: Rows with the ``xcol`` and ``ycol`` columns.
        xcol: Column on the x-axis; values used as given.
        ycol: Column on the y-axis; values used as given.
        expected_df: seg.ucn-layout landmark table; see :func:`get_landmarks`.
            None draws no landmarks.
        group: Group to plot. Selects rows of ``obs_df`` when ``group_col`` is
            present, and the matching landmarks. None takes the first.
        group_col: Column identifying a group.
        clones: Clones defining a landmark; see :func:`get_landmarks`.
        normal: Normal clone to exclude when ``clones`` is None, since a
            landmark labels the tumor states. None keeps every clone.
        hue: (n,) categorical label per point, or a column name in ``obs_df``.
            Drives both the point colors and the marginal curves.
        palette: {label: color} for ``hue``.
        colors: Explicit per-point colors, taking precedence over ``hue``.
        alphas: (n,) alpha per point, or a column name in ``obs_df``.
        markersize: Marker size; None picks one from the point count.
        landmark_size: Landmark marker size.
        xlim: Fixed x limits, or None to autoscale.
        ylim: Fixed y limits, or None to autoscale.
        xlabel: X-axis label; defaults to ``xcol``.
        ylabel: Y-axis label; defaults to ``ycol``.
        title: Figure title.
        refline_x: X value for a reference line, or None. 0.5 suits a BAF axis.
        refline_y: Y value for a reference line, or None.
        show_marginals: Draw the marginal KDE curves.
        show_props: Add a clone-proportion legend when ``expected_df`` carries
            ``u_<clone>`` columns.
        display_min_clone_prop: Hide clones below this proportion from the
            landmark labels and the legend; see :func:`get_landmarks`.
        rasterized: Rasterize the point collection.

    Returns:
        The ``JointGrid``. The caller saves and closes ``grid.figure``.

    Raises:
        ValueError: If ``xcol`` or ``ycol`` is missing from ``obs_df``.
    """
    for col in (xcol, ycol):
        if col not in obs_df.columns:
            raise ValueError(f"obs_df has no column {col!r}")
    if group_col in obs_df.columns:
        obs_df = select_sample(obs_df, group, col=group_col)
    if isinstance(hue, str):
        hue = obs_df[hue].to_numpy()
    if isinstance(alphas, str):
        alphas = obs_df[alphas].to_numpy()

    x = obs_df[xcol].to_numpy(dtype=float)
    y = obs_df[ycol].to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    rgba = resolve_colors(
        len(obs_df), hue=hue, palette=palette, colors=colors, alphas=alphas
    )
    if markersize is None:
        markersize = resolve_marker_size(int(keep.sum()))

    hue_keep = None if hue is None else np.asarray(hue)[keep]
    grid = sns.JointGrid(
        x=x[keep], y=y[keep], hue=hue_keep, palette=palette, xlim=xlim, ylim=ylim
    )
    grid.plot_joint(sns.scatterplot, s=markersize, legend=False, edgecolor="none")
    if show_marginals:
        grid.plot_marginals(sns.kdeplot, common_norm=False, linewidth=0.8, fill=False)
    if refline_x is not None or refline_y is not None:
        grid.refline(x=refline_x, y=refline_y)

    scatter = grid.ax_joint.collections[0]
    scatter.set_rasterized(rasterized)
    # Override seaborn's own resolution so hue, explicit colors, and alphas all
    # route through one place and match the 1D plot point for point.
    scatter.set_alpha(None)
    scatter.set_facecolors(rgba[keep])

    if expected_df is not None:
        landmarks = get_landmarks(
            expected_df,
            xcol,
            ycol,
            group=group,
            clones=clones,
            normal=normal,
            group_col=group_col,
            display_min_clone_prop=display_min_clone_prop,
        )
        _annotate_landmarks(grid.ax_joint, landmarks, markersize=landmark_size)

        if show_props:
            sub = select_sample(expected_df, group, col=group_col)
            names = clones
            if names is None:
                normal_clone, tumor = get_clone_names(sub, normal=normal)
                names = ([normal_clone] if normal_clone else []) + tumor
            props = get_clone_proportions(sub, names)
            # The normal clone always stays listed, however small: it is the
            # baseline the tumor proportions are read against.
            shown = _visible_clones(
                [n for n in names if n != normal], props, display_min_clone_prop
            )
            if props:
                handles = [
                    Line2D([0], [0], alpha=0, label=f"{name}: {prop:.3f}")
                    for name, prop in props.items()
                    if name == normal or name in shown
                ]
                grid.ax_joint.legend(
                    handles=handles,
                    loc="best",
                    fontsize="small",
                    fancybox=True,
                    framealpha=0.7,
                    handlelength=0,
                    handletextpad=0,
                )

    grid.set_axis_labels(
        xlabel=xcol if xlabel is None else xlabel,
        ylabel=ycol if ylabel is None else ylabel,
        fontsize=12,
    )
    if title is not None:
        grid.figure.suptitle(title, fontsize=12, fontweight="bold")
    grid.figure.tight_layout()
    return grid
