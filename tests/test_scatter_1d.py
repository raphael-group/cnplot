"""1D genome scatter: single axes, multi-sample panel, filtering, gap handling."""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import PathCollection

from cnplot.cnplot_1d import _expected_lines, plot_scatter_1d
from cnplot.cnplot_figures import make_row_spec, plot_scatter_1d_multisample
from cnplot.cnplot_intcnp import plot_cnv_profile


def _scatters(ax):
    return [c for c in ax.collections if isinstance(c, PathCollection)]


def test_scatter_draws_mapped_points(axis, sim, state_palette, saved):
    obs = sim.obs[sim.obs["SAMPLE"] == "S1"]
    fig, ax = plt.subplots(figsize=(12, 2))
    assert (
        plot_scatter_1d(
            ax, obs, axis, "RD", hue="state", palette=state_palette, ylim=(0, 3)
        )
        is None
    )
    n = int((axis.build_coordinates(obs).mapped & np.isfinite(obs["RD"])).sum())
    assert len(_scatters(ax)[0].get_offsets()) == n
    saved(fig, "scatter_1d")


def test_missing_column_raises(axis, sim):
    with pytest.raises(ValueError):
        fig, ax = plt.subplots()
        plot_scatter_1d(ax, sim.obs, axis, "NOPE")
    plt.close("all")


def test_expected_column_hue_needs_palette(axis, sim):
    obs = sim.obs[sim.obs["SAMPLE"] == "S1"]
    with pytest.raises(ValueError):
        fig, ax = plt.subplots()
        plot_scatter_1d(ax, obs, axis, "RD", hue="state", palette=None)
        # hue as array with no palette
        plot_scatter_1d(ax, obs, axis, "RD", hue=np.array(["x"] * len(obs)))
    plt.close("all")


def test_keep_col_splits_without_hiding(axis, sim, state_palette):
    obs = sim.obs[sim.obs["SAMPLE"] == "S1"].reset_index(drop=True)
    mapped = axis.build_coordinates(obs).mapped & np.isfinite(obs["RD"])
    flag = obs["pass_qc"].to_numpy()
    fig, ax = plt.subplots()
    plot_scatter_1d(
        ax, obs, axis, "RD", keep_col="pass_qc", hue="state", palette=state_palette
    )
    counts = [len(c.get_offsets()) for c in _scatters(ax)]
    assert sum(counts) == int(mapped.sum())
    assert counts[0] == int((mapped & ~flag).sum())  # filtered first
    assert counts[1] == int((mapped & flag).sum())
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == [f"filtered ({counts[0]})", f"kept ({counts[1]})"]
    plt.close(fig)


def test_multisample_shared_profile(axis, sim2, saved):
    rows = [make_row_spec("RD", ylim=(0, 3), href=1.0), make_row_spec("BAF", href=0.5)]
    fig = plot_scatter_1d_multisample(
        sim2.obs, axis, rows, expected_df=sim2.expected_1d, seg_df=sim2.seg_ucn
    )
    axs = fig.get_axes()
    # 2 samples x 2 rows + profile + legend
    assert len(axs) == 2 * 2 + 2
    # profile carries no per-clone proportion, each sample's last row does
    prof_labels = [t.get_text() for t in axs[4].get_yticklabels()]
    assert all("prop" not in t for t in prof_labels)
    legs = [axs[i].get_legend() for i in (1, 3)]
    assert all(lg is not None for lg in legs)
    texts = [[t.get_text() for t in lg.get_texts()] for lg in legs]
    assert texts[0] != texts[1]
    saved(fig, "multisample_1d")


def _n_dashes(ax):
    return len([ln for ln in ax.lines if ln.get_linestyle() != "-"])


def test_dashes_only_when_gaps_collapsed(axis, axis_keep, sim):
    obs = sim.obs[sim.obs["SAMPLE"] == "S1"]
    n_interior = len(axis.interior_gaps)
    assert n_interior == len(axis_keep.interior_gaps) >= 1

    fig, ax = plt.subplots()
    plot_scatter_1d(ax, obs, axis, "RD")
    assert _n_dashes(ax) == n_interior
    plt.close(fig)

    fig, ax = plt.subplots()
    plot_scatter_1d(ax, obs, axis_keep, "RD")
    assert _n_dashes(ax) == 0
    plt.close(fig)


def test_expected_lines_break_at_gaps(axis, axis_keep, sim):
    exp = sim.expected_1d
    gaps_keep = [(g.start, g.end) for g in axis_keep.interior_gaps if g.width > 0]

    def crossings(ga):
        lines = _expected_lines(exp, ga, "RD", "S1")
        gaps = [(g.start, g.end) for g in ga.interior_gaps if g.width > 0]
        return sum(
            any(ln[0][0] < ge and ln[1][0] > gs for gs, ge in gaps) for ln in lines
        )

    assert gaps_keep, "kept layout must have gaps with width"
    assert crossings(axis) == 0
    assert crossings(axis_keep) == 0


def test_profile_and_scatter_share_gap_rule(axis_keep, sim):
    obs = sim.obs[sim.obs["SAMPLE"] == "S1"]
    fig, ax = plt.subplots()
    plot_cnv_profile(ax, sim.seg_ucn, axis_keep, sample_id="S1")
    assert _n_dashes(ax) == 0
    plt.close(fig)
    fig, ax = plt.subplots()
    plot_scatter_1d(ax, obs, axis_keep, "RD")
    assert _n_dashes(ax) == 0
    plt.close(fig)
