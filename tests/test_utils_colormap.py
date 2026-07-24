"""Palettes and the shared resolvers - pure functions, no rendering."""

import matplotlib.colors as mcolors
import numpy as np

from cnplot import (
    MARKER_SIZE_LARGE,
    MARKER_SIZE_SMALL,
    get_baf_cmap,
    get_categorical_cmap,
    get_cn_cmap,
    get_log2rdr_cmap,
    get_multiclass_cmap,
    resolve_colors,
    resolve_marker_size,
    resolve_ylim,
    resolve_ylim_scaled,
)


def test_cn_palette_is_allele_symmetric():
    style, _ = get_cn_cmap()
    assert style[(2, 1)] == style[(1, 2)]
    assert "default" in style


def test_marker_size_switch():
    assert resolve_marker_size(100) == MARKER_SIZE_LARGE
    assert resolve_marker_size(4999) == MARKER_SIZE_LARGE
    assert resolve_marker_size(5000) == MARKER_SIZE_SMALL


def test_ylim_picks_tightest_window():
    assert resolve_ylim([0.1, -0.4]) == (-2, 2)
    assert resolve_ylim([0.1, -4.0]) == (-5, 5)
    assert resolve_ylim([0.1], expected=[3.0]) == (-5, 5)
    assert resolve_ylim([np.nan]) == (-5, 5)


def test_ylim_scaled_caps_outliers():
    lo, hi = resolve_ylim_scaled([0.5, 1.0])
    assert (lo, hi) == (-0.1, 2.0)
    assert resolve_ylim_scaled([100.0])[1] == 6.0


def test_resolve_colors_hue_and_alpha():
    rgba = resolve_colors(3, hue=["a", "b", "a"], palette={"a": "red", "b": "blue"})
    assert rgba.shape == (3, 4)
    assert np.allclose(rgba[0], rgba[2]) and not np.allclose(rgba[0], rgba[1])
    faded = resolve_colors(2, colors="red", alphas=[0.2, 0.9])
    assert np.allclose(faded[:, 3], [0.2, 0.9])


def test_resolve_colors_unmapped_falls_back():
    out = resolve_colors(1, hue=["z"], palette={}, default="steelblue")
    assert np.allclose(out[0], mcolors.to_rgba("steelblue"))


def test_baf_cmap_builds():
    cmap, norm, ticks = get_baf_cmap()
    assert cmap is not None and norm is not None
    assert ticks == [0.0, 0.25, 0.5, 0.75, 1.0]


def test_log2rdr_cmap_builds():
    cmap, norm, ticks = get_log2rdr_cmap()
    assert cmap is not None and norm is not None
    assert norm.vcenter == 0 and ticks == [-1.0, -0.5, 0.0, 0.5, 1.0]


def test_label_cmaps_cover_every_value(sim):
    row_labels = sim.heatmap_labels
    maps = get_multiclass_cmap({"clone": row_labels}, primary_label="clone")
    assert set(maps["clone"]) >= set(np.unique(row_labels).astype(str))


def test_get_categorical_cmap():
    cmap = get_categorical_cmap(["d1", "d2", "NA", "d1"], "Set1")
    assert cmap["NA"] == "darkgray"
    assert cmap["d1"] != cmap["d2"] and cmap["d1"] != "darkgray"
    assert set(cmap) == {"d1", "d2", "NA"}
    # single-label maps are self-contained; normal-like -> normal gray
    clones = get_categorical_cmap(["normal", "clone1", "clone2"], "tab10")
    assert clones["normal"] == "lightgray"
    assert clones["clone1"] != clones["clone2"]
