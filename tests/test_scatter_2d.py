"""2D scatter and the copy-number landmarks it draws."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cnplot import get_landmarks, plot_scatter_2d


def test_landmarks_dedupe_on_state(sim):
    lm = get_landmarks(sim.expected_2d, "BAF", "RD", group="S1")
    states = {s["states"].tobytes() for s in lm}
    assert len(lm) == len(states)
    # normal excluded by default, so labels are the tumor states
    assert any(s["clonal"] for s in lm)


def test_clones_subset_switches_convention(sim):
    joint = get_landmarks(sim.expected_2d, "BAF", "RD", group="S1")
    single = get_landmarks(sim.expected_2d, "BAF", "RD", group="S1", clones=["clone1"])
    # a single-clone label is one "(a,b)"; a joint label separates states with "),("
    assert all("),(" not in s["label"] for s in single)
    # chr3 disagrees (2|1 vs 1|2), so the joint set has a multi-state label
    assert any("),(" in s["label"] for s in joint)
    assert len(joint) >= len(single)


def test_clonal_flag_matches_agreement(sim):
    lm = get_landmarks(sim.expected_2d, "BAF", "RD", group="S1")
    for s in lm:
        distinct = {tuple(ab) for ab in s["states"]}
        assert s["clonal"] == (len(distinct) == 1)


def test_scatter_renders(sim, state_palette, saved):
    grid = plot_scatter_2d(
        sim.obs,
        "BAF",
        "RD",
        expected_df=sim.expected_2d,
        group="S1",
        hue="state",
        palette=state_palette,
        xlim=(0, 1),
        ylim=(0, 3),
    )
    pts = grid.ax_joint.collections[0].get_offsets()
    obs = sim.obs[sim.obs["SAMPLE"] == "S1"]
    n = int((np.isfinite(obs["BAF"]) & np.isfinite(obs["RD"])).sum())
    assert len(pts) == n
    # every landmark gets a marker
    lm = get_landmarks(sim.expected_2d, "BAF", "RD", group="S1")
    n_marks = sum(len(c.get_offsets()) for c in grid.ax_joint.collections[1:])
    assert n_marks == len(lm)
    saved(grid.figure, "scatter_2d")


def test_prop_legend_from_u_columns(sim):
    grid = plot_scatter_2d(
        sim.obs, "BAF", "RD", expected_df=sim.expected_2d, group="S1", show_props=True
    )
    leg = grid.ax_joint.get_legend()
    assert leg is not None and len(leg.get_texts()) == 3  # normal + 2 clones
    plt.close(grid.figure)


def test_missing_value_column_raises(sim):
    with pytest.raises(ValueError):
        plot_scatter_2d(sim.obs, "BAF", "NOPE")
    plt.close("all")
