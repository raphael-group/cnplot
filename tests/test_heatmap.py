"""Single-cell heatmap: mesh, masked gaps, side strips, legends."""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import TwoSlopeNorm

from cnplot import (
    get_baf_cmap,
    plot_column_strips,
    plot_heatmap,
    plot_heatmap_cnp,
    plot_strip_legend,
)


def test_heatmap_renders(axis, sim, saved):
    cmap, norm, _ = get_baf_cmap()
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
    cmap, norm, _ = get_baf_cmap()
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
    cmap, norm, _ = get_baf_cmap()
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
    strips = [
        {"name": "clone", "values": sim.heatmap_labels},
        {"name": "cell_type", "values": celltype},
    ]

    n_before = len(fig.axes)
    info = plot_column_strips(fig, ax, y_edges, strips)
    assert len(fig.axes) == n_before + len(strips)
    for _name, color_dict, prop_dict in info:
        assert abs(sum(prop_dict.get(v, 0.0) for v in color_dict) - 1.0) < 1e-9

    n_leg = len(fig.legends)
    plot_strip_legend(fig, ax, info, x0=0.92)
    assert len(fig.legends) == n_leg + len(info)
    plt.close(fig)


def test_plot_heatmap_draws_strips_and_legends(axis, sim, saved):
    cmap, norm, _ = get_baf_cmap()
    celltype = np.where(sim.heatmap_labels == "normal", "normal", "tumor")
    strips = [
        {"name": "clone", "values": sim.heatmap_labels, "display_name": "Copy-typing"},
        {"name": "cell_type", "values": celltype, "display_name": "Cell-Type"},
    ]
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
        strips=strips,
    )
    # base axes + one per strip
    assert len(fig.axes) == 1 + len(strips)
    assert len(fig.legends) == len(strips)
    saved(fig, "heatmap_full")


def test_dist_strip_draws_and_legends(axis, sim):
    cmap, norm, _ = get_baf_cmap()
    order = list(np.unique(sim.heatmap_labels))
    color_dict = {v: c for v, c in zip(order, ["C0", "C1", "C2"], strict=False)}
    # per-row posterior peaked at the true clone
    post = np.array(
        [[1.0 if v == lab else 0.0 for v in order] for lab in sim.heatmap_labels]
    )
    prop = {v: float((sim.heatmap_labels == v).mean()) for v in order}
    dist = {
        "name": "Copy-typing",
        "matrix": post,
        "order": order,
        "cmap": color_dict,
        "props": prop,
    }

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
        strips=[dist],
    )
    # one strip axes for the distribution, and one legend for it
    assert len(fig.axes) == n_before + 1
    assert len(fig.legends) == 1
    # the distribution legend uses the "{label}: xx.xx%" format
    texts = [t.get_text() for t in fig.legends[0].get_texts()]
    assert all(": " in t and t.endswith("%") for t in texts), texts
    plt.close(fig)


def test_continuous_strip_draws_without_legend(axis, sim):
    from matplotlib.colors import Normalize

    cmap, norm, _ = get_baf_cmap()
    rng = np.random.default_rng(0)
    purity = rng.random(len(sim.heatmap_labels))
    # continuous strip closest to the heatmap, then a categorical strip
    strips = [
        {
            "name": "purity",
            "scalar": purity,
            "cmap": "magma_r",
            "norm": Normalize(0, 1),
        },
        {"name": "clone", "values": sim.heatmap_labels},
    ]
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
        strips=strips,
    )
    # two strip axes, but only the categorical strip carries a legend
    assert len(fig.axes) == n_before + 2
    assert len(fig.legends) == 1
    plt.close(fig)


def test_continuous_strip_colorbar(axis, sim):
    from matplotlib.colors import Normalize

    cmap, norm, _ = get_baf_cmap()
    rng = np.random.default_rng(0)
    purity = rng.random(len(sim.heatmap_labels))
    strips = [
        {
            "name": "purity",
            "scalar": purity,
            "cmap": "magma_r",
            "norm": Normalize(0, 1),
            "show_cbar": True,
            "cbar_ticks": [0, 0.5, 1],
        },
    ]
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
        strips=strips,
    )
    # the strip axes + its colorbar axes; no swatch legend
    assert len(fig.axes) == n_before + 2
    assert len(fig.legends) == 0
    plt.close(fig)


def test_mb_ticks_on_heatmap(axis, sim, saved):
    cmap, norm, _ = get_baf_cmap()
    fig, ax = plt.subplots(figsize=(12, 4))
    plot_heatmap(
        ax,
        sim.heatmap_baf,
        sim.bins,
        axis,
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        mb_ticks=True,
        mb_tick_step=20_000_000,
    )
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels and all(lbl.isdigit() for lbl in labels)  # Mb integers
    assert "0" not in labels  # chromosome-start tick dropped
    # chromosome names drawn as bold off-axis text, one per chromosome
    names = [t for t in ax.texts if t.get_text()]
    assert len(names) == len(axis.chrs)
    assert all(t.get_fontweight() == "bold" for t in names)
    saved(fig, "heatmap_mb_ticks")


def test_colorbar_adds_axes(axis, sim):
    cmap, norm, _ = get_baf_cmap()
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

    cmap, norm, _ = get_baf_cmap()
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
        strips=[{"name": "cell_type", "values": celltype}],
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
    # chromosome labels sit on the profile (top), not the heatmap
    prof_labels = [t.get_text() for t in profile_ax.get_xticklabels() if t.get_text()]
    heat_labels = [t.get_text() for t in heatmap_ax.get_xticklabels() if t.get_text()]
    assert prof_labels and heat_labels == []
    saved(fig, "heatmap_cnp_page")


def test_plot_heatmap_cnp_ascn_profile(axis, sim):
    from matplotlib.patches import Rectangle

    cmap, norm, _ = get_baf_cmap()
    fig = plot_heatmap_cnp(
        sim.heatmap_baf,
        sim.bins,
        axis,
        sim.seg_ucn,
        sample_id="S1",
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        profile="ascn",
    )
    # allele-specific profile splits each clone into A and B bars: >= 4x segments
    n_seg = len(sim.seg_ucn[sim.seg_ucn["SAMPLE"] == "S1"])
    profile_ax = fig.axes[1]
    rects = [p for p in profile_ax.patches if isinstance(p, Rectangle)]
    assert len(rects) >= 4 * n_seg
    plt.close(fig)

    with pytest.raises(ValueError):
        plot_heatmap_cnp(sim.heatmap_baf, sim.bins, axis, sim.seg_ucn, profile="bad")
    plt.close("all")
