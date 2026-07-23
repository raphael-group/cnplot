"""Integer copy-number profile and its self-drawn legend."""

import matplotlib.pyplot as plt
import pandas as pd

from cnplot.cnplot_intcnp import (
    get_clone_names,
    get_clone_states,
    has_mirror,
    plot_cnv_profile,
)


def _cnv_swatch(ax_leg):
    return any(p.get_hatch() == "/" for p in ax_leg.patches)


def _mirrored_profile(mirror=True):
    """A tiny two-segment chr1 profile; the first segment mirrors when asked."""
    c1, c2 = ("2|1", "1|2") if mirror else ("2|1", "2|1")
    return pd.DataFrame(
        {
            "SAMPLE": ["M", "M"],
            "#CHR": ["chr1", "chr1"],
            "START": [0, 10_000_000],
            "END": [10_000_000, 20_000_000],
            "cn_normal": ["1|1", "1|1"],
            "u_normal": [0.3, 0.3],
            "cn_clone1": [c1, "1|1"],
            "u_clone1": [0.3, 0.3],
            "cn_clone2": [c2, "1|1"],
            "u_clone2": [0.4, 0.4],
        }
    )


def test_profile_renders(axis, sim, saved):
    fig, ax = plt.subplots(figsize=(12, 2))
    assert plot_cnv_profile(ax, sim.seg_ucn, axis, sample_id="S1") is None
    n_seg = len(sim.seg_ucn[sim.seg_ucn["SAMPLE"] == "S1"])
    # two tumor clones (normal excluded) -> two stacked rows of rectangles
    assert len(ax.patches) >= 2 * n_seg
    saved(fig, "profile")


def test_clone_names_split_normal(sim):
    normal, tumor = get_clone_names(sim.seg_ucn)
    assert normal == "normal"
    assert tumor == ["clone1", "clone2"]


def test_mirror_rule(sim):
    # the real profile has no mirrored segment
    seg = sim.seg_ucn[sim.seg_ucn["SAMPLE"] == "S1"]
    assert not has_mirror(get_clone_states(seg, ["clone1", "clone2"]))
    # a clone1=2|1 vs clone2=1|2 segment is mirrored
    mirrored = _mirrored_profile(mirror=True)
    assert has_mirror(get_clone_states(mirrored, ["clone1", "clone2"]))


def test_legend_swatch_tracks_data(axis):
    fig, (ax, lg) = plt.subplots(2, 1)
    plot_cnv_profile(ax, _mirrored_profile(mirror=True), axis, ax_leg=lg)
    assert _cnv_swatch(lg), "mirrored data must show the swatch"
    plt.close(fig)

    fig, (ax, lg) = plt.subplots(2, 1)
    plot_cnv_profile(ax, _mirrored_profile(mirror=False), axis, ax_leg=lg)
    assert not _cnv_swatch(lg)
    plt.close(fig)


def test_show_mirror_false_suppresses_swatch(axis):
    fig, (ax, lg) = plt.subplots(2, 1)
    plot_cnv_profile(
        ax, _mirrored_profile(mirror=True), axis, ax_leg=lg, show_mirror=False
    )
    assert not _cnv_swatch(lg)
    plt.close(fig)


def test_pi_viol_overlay_present(axis, sim):
    fig, ax = plt.subplots()
    plot_cnv_profile(ax, sim.seg_ucn, axis, sample_id="S1", show_pi_viol=True)
    # PI_VIOL draws colored horizontal segments above the stack
    colors = {
        tuple(ln.get_color()) if not isinstance(ln.get_color(), str) else ln.get_color()
        for ln in ax.lines
    }
    assert colors, "expected PI_VIOL overlay lines"
    plt.close(fig)


def test_unknown_sample_raises(axis, sim):
    import pytest

    with pytest.raises(ValueError):
        fig, ax = plt.subplots()
        plot_cnv_profile(ax, sim.seg_ucn, axis, sample_id="NOPE")
    plt.close("all")
