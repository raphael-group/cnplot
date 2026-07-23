"""Palettes and the shared resolvers - pure functions, no rendering."""

import matplotlib.colors as mcolors
import numpy as np

from cnplot.cnplot_colormap import (
    build_label_cmaps,
    get_baf_cmap,
    get_cn_cmap,
)
from cnplot.cnplot_utils import (
    MARKER_SIZE_LARGE,
    MARKER_SIZE_SMALL,
    format_clone_name,
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


def test_format_clone_name():
    assert format_clone_name("normal") == "Normal"
    assert format_clone_name("clone2") == "Clone 2"
    assert format_clone_name("clone2", plot_clone_name=False) == "2"
    assert format_clone_name("T_cell") == "T_cell"


def test_baf_cmap_builds():
    cmap, norm = get_baf_cmap()
    assert cmap is not None and norm is not None


def test_label_cmaps_cover_every_value(sim):
    row_labels = sim.heatmap_labels
    maps = build_label_cmaps({"clone": row_labels}, primary_label="clone")
    assert set(maps["clone"]) >= set(np.unique(row_labels).astype(str))
