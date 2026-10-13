"""Measured layout works on a canvas that has no Agg renderer of its own."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cnplot import get_baf_cmap, plot_cnv_profile, plot_heatmap_cnp


@pytest.fixture(params=["pdf", "svg"])
def backend(request):
    """Switch to a vector backend for one test, then back."""
    previous = plt.get_backend()
    plt.switch_backend(request.param)
    yield request.param
    plt.switch_backend(previous)


def test_profile_key_draws(backend, axis, sim, tmp_path):
    fig, (ax, key) = plt.subplots(2, 1, figsize=(12, 3), height_ratios=[3, 1])
    plot_cnv_profile(ax, sim.seg_ucn, axis, sample_id="S1", ax_leg=key)
    assert key.texts
    fig.savefig(tmp_path / f"profile.{backend}", bbox_inches="tight")
    plt.close(fig)


def test_heatmap_page_draws(backend, axis, sim, tmp_path):
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
    assert fig.legends
    fig.savefig(tmp_path / f"page.{backend}", bbox_inches="tight")
    plt.close(fig)
