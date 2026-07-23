"""Single-cell copy-number heatmap: a per-row matrix along the genome.

Copytyping-only; there is no second consumer, so this ports the drawing rather
than deduplicating across repos. The split follows M5: cnplot draws a matrix that
is already reduced to rows, and the reducers stay upstream.

``plot_heatmap`` takes an ``(n_rows, n_bins)`` matrix and a bin table, places the
columns through a :class:`~cnplot.cnplot_genome_axis.GenomeAxis`, and draws a
``pcolormesh`` with the shared chromosome decoration. Building that matrix -
aggregating cells into rows, pooling counts, turning them into RDR or BAF - needs
``copytyping.inference.model_utils`` and stays in Copytyping (``_row_layout``,
``prepare_rdr`` / ``prepare_baf`` / ``prepare_pi_gk``, ``plot_cnv_heatmap``).

The categorical side strips and their legends are general, so they live here.
"""

import logging

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle

from .cnplot_genome_axis import GenomeAxis
from .cnplot_utils import (
    decorate_genome_axis,
    draw_chr_boundaries,
    draw_segment_boundaries,
)

logger = logging.getLogger(__name__)

__all__ = [
    "plot_column_strips",
    "plot_heatmap",
    "plot_strip_legend",
]


# =============================================================================
# Heatmap
# =============================================================================


def plot_heatmap(
    ax: plt.Axes,
    matrix: np.ndarray,
    coords_df: pd.DataFrame,
    genome_axis: GenomeAxis,
    row_labels: np.ndarray | None = None,
    height: float = 1.0,
    plot_chrname: bool = True,
    title: str | None = None,
    ylabel: str | None = None,
    cmap=None,
    norm=None,
    show_block_labels: bool = True,
    show_gaps: bool = True,
    show_colorbar: bool = False,
    cbar_label: str | None = None,
    cbar_ticks=None,
    cbar_frac: float = 0.2,
    strips: list | None = None,
    show_strip_legend: bool = True,
    legend_x0: float | None = None,
    rasterized: bool = True,
) -> tuple:
    """Draw a per-row value matrix as a genome-wide heatmap.

    Columns are placed by coordinate through ``genome_axis``: each bin becomes one
    mesh cell, and stretches no bin covers become masked filler so nothing is
    stretched across a gap. Rows are drawn bottom to top, matching ``pcolormesh``.

    When ``strips`` is given the side strips and their legends are drawn too, so
    a single call produces the full annotated figure.

    Args:
        ax: Axes to draw on.
        matrix: (n_rows, n_bins) values, columns aligned to ``coords_df`` rows.
        coords_df: The bins, with "#CHR", "START", "END"; one column per bin.
        genome_axis: Axis placing the columns and supplying the decoration.
        row_labels: (n_rows,) label per row, bottom to top. Contiguous runs get a
            bounding box and a y-tick; None draws no row labels.
        height: Mesh height in axes fraction.
        plot_chrname: Draw chromosome labels along the top.
        title: Axes title.
        ylabel: Y-axis label.
        cmap: Colormap for the mesh.
        norm: Normalization for the mesh.
        show_block_labels: Write "<label> (<pct>%)" at each row block; else the
            labels live in the side strips.
        show_gaps: Dash the collapsed gaps, as the profile does.
        show_colorbar: Draw a colorbar for the mesh values, right of the axes.
        cbar_label: Colorbar label.
        cbar_ticks: Explicit colorbar ticks, or None to autoscale.
        cbar_frac: Colorbar height as a fraction of the mesh height, matching
            Copytyping's short bar.
        strips: Ordered list of side-strip specs drawn left of the heatmap, the
            first closest to it; see :func:`plot_column_strips` for the spec
            format. None draws no strips.
        show_strip_legend: Draw one legend per strip to the right.
        legend_x0: Left edge of the legends in figure fraction; defaults to just
            right of the axes.
        rasterized: Rasterize the mesh, keeping vector output small.

    Returns:
        (x_edges, y_edges, masked_matrix): the mesh geometry and the masked array
        drawn, aligned to the strips.

    Raises:
        ValueError: If ``matrix`` columns do not match ``coords_df`` rows.
    """
    n_rows, n_bins = matrix.shape
    if n_bins != len(coords_df):
        raise ValueError(f"matrix has {n_bins} columns for {len(coords_df)} bins")

    x_edges, col_bin_ids = genome_axis.grid(coords_df)
    x_edges = np.asarray(x_edges, dtype=float)
    ext = np.full((n_rows, len(col_bin_ids)), np.nan, dtype=float)
    for j, bid in enumerate(col_bin_ids):
        if bid >= 0:
            ext[:, j] = matrix[:, bid]

    y_edges = np.linspace(0.0, height, n_rows + 1)
    masked = np.ma.masked_invalid(ext)
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        masked,
        cmap=cmap,
        norm=norm,
        shading="flat",
        rasterized=rasterized,
    )

    draw_chr_boundaries(ax, genome_axis)
    if show_gaps:
        draw_segment_boundaries(ax, genome_axis)

    if row_labels is not None:
        _draw_row_blocks(
            ax, np.asarray(row_labels), y_edges, genome_axis, show_block_labels
        )
    else:
        ax.set_yticks([])
    ax.tick_params(axis="y", length=0)

    ax.set_ylim(0.0, height)
    if ylabel is not None:
        ax.set_ylabel(ylabel, rotation=90, va="center")
    if title is not None:
        ax.set_title(title)
    decorate_genome_axis(ax, genome_axis, plot_chrname=plot_chrname, label_pos="top")
    if plot_chrname:
        plt.setp(ax.get_xticklabels(), fontweight="bold")

    # right-side chrome, left to right: colorbar, then one legend per strip
    right = None
    if show_colorbar:
        ax.figure.canvas.draw()
        box = ax.get_position()
        # short bar (Copytyping uses height / 5), bottom-aligned to the mesh
        cax = ax.figure.add_axes([box.x1 + 0.012, box.y0, 0.01, box.height * cbar_frac])
        cbar = ax.figure.colorbar(mesh, cax=cax)
        if cbar_ticks is not None:
            cbar.set_ticks(cbar_ticks)
        if cbar_label is not None:
            # label above the bar rather than alongside it
            cax.set_title(cbar_label, fontsize=10, fontweight="bold", pad=4)
        right = box.x1 + 0.06

    if strips:
        info = plot_column_strips(ax.figure, ax, y_edges, strips)
        if show_strip_legend:
            if legend_x0 is None:
                ax.figure.canvas.draw()
                right = right if right is not None else ax.get_position().x1 + 0.02
                legend_x0 = right
            plot_strip_legend(ax.figure, ax, info, x0=legend_x0)
    return x_edges, y_edges, masked


def _draw_row_blocks(
    ax: plt.Axes,
    row_labels: np.ndarray,
    y_edges: np.ndarray,
    genome_axis: GenomeAxis,
    show_block_labels: bool,
) -> None:
    """Box each contiguous run of equal row labels and tick it at its center.

    Args:
        ax: Axes to draw on.
        row_labels: (n_rows,) label per row, bottom to top.
        y_edges: (n_rows + 1,) row edges in data coordinates.
        genome_axis: Axis giving the horizontal span for the boxes.
        show_block_labels: Write "<label> (<pct>%)"; else only draw the boxes.
    """
    x0 = genome_axis.ch_coords[0]
    span = genome_axis.chr_end - x0
    n = len(row_labels)
    starts = np.flatnonzero(np.r_[True, row_labels[1:] != row_labels[:-1]])
    ends = np.r_[starts[1:], n]

    yticks, yticklabels = [], []
    for start, end in zip(starts, ends, strict=True):
        y0, y1 = y_edges[start], y_edges[end]
        ax.add_patch(
            Rectangle(
                (x0, y0), span, y1 - y0, fill=False, edgecolor="black", linewidth=1.0
            )
        )
        yticks.append(0.5 * (y0 + y1))
        pct = round(100 * (end - start) / n, 1)
        yticklabels.append(f"{row_labels[start]} ({pct}%)")

    if show_block_labels:
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=11, fontweight="bold")
    else:
        ax.set_yticks([])


# =============================================================================
# Side strips
# =============================================================================


def _resolve_strips(strips: list) -> list:
    """Fill spec defaults: categorical cmaps and display names.

    Categorical strips missing a ``cmap`` are colored together through
    :func:`~cnplot.cnplot_colormap.get_multiclass_cmap` so shared values agree
    across strips; the first such strip is the primary.

    Args:
        strips: The strip specs, each a categorical (``values``), distribution
            (``matrix``), or continuous (``scalar``) dict.

    Returns:
        Shallow copies of ``strips`` with ``cmap`` and ``display_name`` set.
    """
    strips = [dict(s) for s in strips]
    missing = {
        s["name"]: s["values"]
        for s in strips
        if "values" in s and s.get("cmap") is None
    }
    if missing:
        from .cnplot_colormap import get_multiclass_cmap

        built = get_multiclass_cmap(missing, primary_label=next(iter(missing)))
        for s in strips:
            if s["name"] in missing:
                s["cmap"] = built[s["name"]]
    for s in strips:
        s.setdefault("display_name", s["name"])
    return strips


def _draw_dist_strip(ax: plt.Axes, spec: dict, y_edges: np.ndarray) -> tuple:
    """Stack each row's distribution into a horizontal bar.

    Instead of a solid color per row, a distribution strip (e.g. copy-typing
    posteriors) shows assignment confidence rather than just the hard call.

    Args:
        ax: Strip axes to draw on.
        spec: Distribution spec with ``matrix`` (n_rows, len(order)), ``order``,
            ``cmap`` {category: color}, and optional ``props`` {category:
            fraction} for the legend.
        y_edges: Row edges from :func:`plot_heatmap`.

    Returns:
        (display_name, {category: color}, {category: fraction}) for the legend.
    """
    matrix = np.asarray(spec["matrix"])
    order = spec["order"]
    colors = spec["cmap"]
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    heights = np.diff(y_edges)
    left = np.zeros(len(y_centers))
    for k, val in enumerate(order):
        w = matrix[:, k]
        ax.barh(
            y_centers,
            w,
            left=left,
            height=heights,
            color=colors[val],
            edgecolor="none",
            align="center",
            rasterized=True,
        )
        left = left + w
    ax.set_xlim(0.0, 1.0)
    props = spec.get("props") or dict.fromkeys(order, 0.0)
    return spec["display_name"], {v: colors[v] for v in order}, props


def _draw_cont_strip(ax: plt.Axes, spec: dict, y_edges: np.ndarray) -> None:
    """Draw a per-row continuous scalar as a colored column, e.g. purity.

    Args:
        ax: Strip axes to draw on.
        spec: Continuous spec with ``scalar`` (n_rows,), ``cmap`` (a colormap),
            and optional ``norm``.
        y_edges: Row edges from :func:`plot_heatmap`.

    Returns:
        None - a continuous strip carries no categorical legend.
    """
    scalar = np.asarray(spec["scalar"], dtype=float)[:, None]
    ax.pcolormesh(
        np.array([0.0, 1.0]),
        y_edges,
        scalar,
        cmap=spec["cmap"],
        norm=spec.get("norm"),
        shading="flat",
        rasterized=True,
    )
    return None


def _draw_cat_strip(ax: plt.Axes, spec: dict, y_edges: np.ndarray) -> tuple:
    """Draw a solid color per row from the row's categorical value.

    Args:
        ax: Strip axes to draw on.
        spec: Categorical spec with ``values`` (n_rows,) and ``cmap`` {value:
            color}.
        y_edges: Row edges from :func:`plot_heatmap`.

    Returns:
        (display_name, {value: color}, {value: fraction}) for the legend, where
        fraction is each value's share of rows.
    """
    values = np.array([str(v) for v in spec["values"]])
    n_rows = max(len(values), 1)
    color_dict = spec["cmap"]
    order = list(color_dict)
    prop_dict = {v: int((values == v).sum()) / n_rows for v in color_dict}
    codes = np.array([order.index(v) for v in values], dtype=float)[:, None]
    strip_cmap = mcolors.ListedColormap([color_dict[v] for v in order])
    strip_norm = mcolors.BoundaryNorm(np.arange(len(order) + 1) - 0.5, strip_cmap.N)
    ax.pcolormesh(
        np.array([0.0, 1.0]),
        y_edges,
        codes,
        cmap=strip_cmap,
        norm=strip_norm,
        shading="flat",
        rasterized=True,
    )
    return spec["display_name"], color_dict, prop_dict


def _finish_strip_axis(ax: plt.Axes, label: str, base_ax: plt.Axes) -> None:
    """Title a strip along its top and align it to the heatmap rows.

    Args:
        ax: Strip axes to finish.
        label: Text written vertically above the strip.
        base_ax: Heatmap axes whose y-limits the strip matches.
    """
    ax.set_xticks([0.5])
    ax.set_xticklabels([label], rotation=90, fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", labeltop=True, labelbottom=False, top=False, bottom=False)
    ax.set_yticks([])
    ax.set_ylim(base_ax.get_ylim())


def plot_column_strips(
    fig: plt.Figure,
    base_ax: plt.Axes,
    y_edges: np.ndarray,
    strips: list,
    strip_width: float = 0.012,
    gap: float = 0.004,
) -> list:
    """Draw vertical annotation strips left of ``base_ax``, in list order.

    Each entry of ``strips`` is one strip, drawn outward from the heatmap in the
    given order (the first is closest). A spec is one of three kinds, told apart by
    its data key:

    - categorical (``values``): ``{"name", "values": (n_rows,), "cmap"?,
      "display_name"?}`` - a solid color per row from the row's value. ``cmap`` may
      be omitted to auto-color; see :func:`_resolve_strips`.
    - distribution (``matrix``): ``{"name", "matrix": (n_rows, k), "order":
      [k cats], "cmap", "props"?, "display_name"?}`` - each row's distribution
      stacked into a bar, e.g. copy-typing posteriors.
    - continuous (``scalar``): ``{"name", "scalar": (n_rows,), "cmap", "norm"?,
      "display_name"?}`` - a per-row scalar colored through ``cmap`` / ``norm``,
      e.g. purity. Carries no legend.

    Args:
        fig: Figure holding ``base_ax``; strips are added as new axes.
        base_ax: Heatmap axes the strips align to.
        y_edges: Row edges from :func:`plot_heatmap`.
        strips: Ordered list of strip specs, the first closest to the heatmap.
        strip_width: Strip width in figure fraction.
        gap: Space between strips in figure fraction.

    Returns:
        [(display_name, {value: color}, {value: fraction})], one per strip that
        carries a legend (categorical and distribution), in order, for
        :func:`plot_strip_legend`; continuous strips contribute none.
    """
    strips = _resolve_strips(strips)
    fig.canvas.draw()
    bbox = base_ax.get_position()
    x_cursor = bbox.x0 - gap

    legends_info = []
    for spec in strips:
        x_cursor -= strip_width
        ax = fig.add_axes([x_cursor, bbox.y0, strip_width, bbox.height])
        if "scalar" in spec:
            legend = _draw_cont_strip(ax, spec, y_edges)
        elif "matrix" in spec:
            legend = _draw_dist_strip(ax, spec, y_edges)
        else:
            legend = _draw_cat_strip(ax, spec, y_edges)
        _finish_strip_axis(ax, spec["display_name"], base_ax)
        x_cursor -= gap
        if legend is not None:
            legends_info.append(legend)
    return legends_info


def plot_strip_legend(
    fig: plt.Figure,
    base_ax: plt.Axes,
    legends_info: list,
    x0: float,
    entry_h: float = 0.038,
    gap: float = 0.06,
) -> float:
    """Stack one categorical legend per strip, top-aligned, at figure-x ``x0``.

    Args:
        fig: Figure to add the legends to.
        base_ax: Heatmap axes the legends align to at the top.
        legends_info: The list :func:`plot_column_strips` returns.
        x0: Left edge of the legends, in figure fraction.
        entry_h: Height of one legend entry, in figure fraction.
        gap: Space between legends, in figure fraction.

    Returns:
        The y of the last legend's bottom, for anything stacked below.
    """
    fig.canvas.draw()
    bbox = base_ax.get_position()
    y_top = bbox.y1
    for name, color_dict, prop_dict in legends_info:
        handles = [
            Patch(facecolor=col, label=f"{v}: {prop_dict.get(v, 0.0) * 100:.2f}%")
            for v, col in color_dict.items()
        ]
        leg = fig.legend(
            handles=handles,
            title=name,
            loc="upper left",
            bbox_to_anchor=(x0, y_top),
            frameon=False,
            fontsize=13,
            title_fontsize=15,
        )
        leg.get_title().set_fontweight("bold")
        fig.add_artist(leg)
        y_top -= entry_h * (len(handles) + 1) + gap
    return y_top
