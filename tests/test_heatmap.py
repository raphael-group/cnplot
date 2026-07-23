"""Single-cell heatmap: mesh, masked gaps, side strips, legends."""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import TwoSlopeNorm

from cnplot import (
    build_label_cmaps,
    get_baf_cmap,
    plot_column_strips,
    plot_heatmap,
    plot_heatmap_cnp,
    plot_strip_legend,
)


def test_heatmap_renders(axis, sim, saved):
    cmap, norm = get_baf_cmap()
    fig, ax = plt.subplots(figsize=(12, 4))
    x_edges, y_edges, C = plot_heatmap(
        ax,
        sim.heatmap_baf,
        sim.bins,
        axis,
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        ylabel="cells",
    )
    assert len(y_edges) == sim.heatmap_baf.shape[0] + 1
    assert x_edges[0] == axis.axis_start and x_edges[-1] == axis.chr_end
    saved(fig, "heatmap")


def test_filler_columns_are_masked(axis, sim):
    cmap, norm = get_baf_cmap()
    fig, ax = plt.subplots()
    _, _, C = plot_heatmap(ax, sim.heatmap_baf, sim.bins, axis, cmap=cmap, norm=norm)
    _, col_bin_ids = axis.grid(sim.bins)
    fill = np.array([b < 0 for b in col_bin_ids])
    masked = np.array([bool(C[:, j].mask.all()) for j in range(C.shape[1])])
    assert np.array_equal(masked, fill)
    plt.close(fig)


def test_column_mismatch_raises(axis, sim):
    with pytest.raises(ValueError):
        fig, ax = plt.subplots()
        plot_heatmap(ax, sim.heatmap_baf[:, :-1], sim.bins, axis)
    plt.close("all")


def test_row_blocks_labelled(axis, sim):
    fig, ax = plt.subplots()
    plot_heatmap(
        ax,
        sim.heatmap_rdr,
        sim.bins,
        axis,
        row_labels=sim.heatmap_labels,
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
    )
    labels = [t.get_text() for t in ax.get_yticklabels()]
    n_blocks = len(np.unique(sim.heatmap_labels))
    assert len(labels) == n_blocks
    boxes = [p for p in ax.patches if p.get_fill() is False]
    assert len(boxes) == n_blocks
    plt.close(fig)


def test_side_strips_and_legends(axis, sim):
    cmap, norm = get_baf_cmap()
    fig, ax = plt.subplots(figsize=(12, 4))
    _, y_edges, _ = plot_heatmap(
        ax,
        sim.heatmap_baf,
        sim.bins,
        axis,
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        show_block_labels=False,
    )
    celltype = np.where(sim.heatmap_labels == "normal", "normal", "tumor")
    row_label_map = {"clone": sim.heatmap_labels, "cell_type": celltype}
    cmaps = build_label_cmaps(row_label_map, primary_label="clone")

    n_before = len(fig.axes)
    info = plot_column_strips(fig, ax, y_edges, row_label_map, cmaps)
    assert len(fig.axes) == n_before + len(row_label_map)
    for _name, color_dict, prop_dict in info:
        assert abs(sum(prop_dict.get(v, 0.0) for v in color_dict) - 1.0) < 1e-9

    n_leg = len(fig.legends)
    plot_strip_legend(fig, ax, info, x0=0.92)
    assert len(fig.legends) == n_leg + len(info)
    plt.close(fig)


def test_plot_heatmap_draws_strips_and_legends(axis, sim, saved):
    cmap, norm = get_baf_cmap()
    celltype = np.where(sim.heatmap_labels == "normal", "normal", "tumor")
    strip_label_map = {"clone": sim.heatmap_labels, "cell_type": celltype}
    fig, ax = plt.subplots(figsize=(12, 4))
    # one call: mesh + strips + legends, colors auto-built
    plot_heatmap(
        ax,
        sim.heatmap_baf,
        sim.bins,
        axis,
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        show_block_labels=False,
        strip_label_map=strip_label_map,
        display_names={"clone": "Copy-typing", "cell_type": "Cell-Type"},
    )
    # base axes + one per strip
    assert len(fig.axes) == 1 + len(strip_label_map)
    assert len(fig.legends) == len(strip_label_map)
    saved(fig, "heatmap_full")


def test_dist_strip_draws_and_legends(axis, sim):
    cmap, norm = get_baf_cmap()
    order = list(np.unique(sim.heatmap_labels))
    color_dict = {v: c for v, c in zip(order, ["C0", "C1", "C2"], strict=False)}
    # per-row posterior peaked at the true clone
    post = np.array(
        [[1.0 if v == lab else 0.0 for v in order] for lab in sim.heatmap_labels]
    )
    prop = {v: float((sim.heatmap_labels == v).mean()) for v in order}
    dist = ("Copy-typing", post, order, color_dict, prop)

    fig, ax = plt.subplots(figsize=(12, 4))
    n_before = len(fig.axes)
    plot_heatmap(
        ax,
        sim.heatmap_baf,
        sim.bins,
        axis,
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        show_block_labels=False,
        dist_strip=dist,
    )
    # one strip axes for the distribution, and one legend for it
    assert len(fig.axes) == n_before + 1
    assert len(fig.legends) == 1
    # the distribution legend uses the "{label}: xx.xx%" format
    texts = [t.get_text() for t in fig.legends[0].get_texts()]
    assert all(": " in t and t.endswith("%") for t in texts), texts
    plt.close(fig)


def test_colorbar_adds_axes(axis, sim):
    cmap, norm = get_baf_cmap()
    fig, ax = plt.subplots(figsize=(12, 4))
    n_before = len(fig.axes)
    plot_heatmap(
        ax,
        sim.heatmap_baf,
        sim.bins,
        axis,
        cmap=cmap,
        norm=norm,
        show_colorbar=True,
        cbar_label="BAF",
        cbar_ticks=[0, 0.5, 1],
    )
    assert len(fig.axes) == n_before + 1  # the colorbar axes
    plt.close(fig)


def test_plot_heatmap_cnp_page(axis, sim, saved):
    import matplotlib.pyplot as plt
    from matplotlib.collections import QuadMesh
    from matplotlib.patches import Rectangle

    cmap, norm = get_baf_cmap()
    celltype = np.where(sim.heatmap_labels == "normal", "normal", "tumor")
    fig = plot_heatmap_cnp(
        sim.heatmap_baf,
        sim.bins,
        axis,
        sim.seg_ucn,
        sample_id="S1",
        title="page",
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        cbar_label="BAF",
        cbar_ticks=[0, 0.5, 1],
        strip_label_map={"cell_type": celltype},
    )
    assert isinstance(fig, plt.Figure)
    heatmap_ax, profile_ax = fig.axes[0], fig.axes[1]
    # heatmap mesh on the top axes
    assert any(isinstance(c, QuadMesh) for c in heatmap_ax.collections)
    # CN profile rectangles on the middle axes (two clones x segments)
    n_seg = len(sim.seg_ucn[sim.seg_ucn["SAMPLE"] == "S1"])
    assert len([p for p in profile_ax.patches if isinstance(p, Rectangle)]) >= 2 * n_seg
    # colorbar + one strip axes were added beyond the three base axes
    assert len(fig.axes) >= 3 + 2
    # only the heatmap carries chromosome labels; the profile does not
    prof_labels = [t.get_text() for t in profile_ax.get_xticklabels() if t.get_text()]
    assert prof_labels == []
    saved(fig, "heatmap_cnp_page")
