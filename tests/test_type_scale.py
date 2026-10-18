"""Text is drawn on the shared type scale, and the CN key wraps instead of shrinking."""

import math

import matplotlib.pyplot as plt
import pytest

from cnplot import FS_ANNOT, FS_LEGEND, FS_TICK, chrname_tiers, decorate_genome_axis
from cnplot.cnplot_intcnp import _plot_cnv_legend, _wrap_groups


def _key_rows(width):
    fig, ax = plt.subplots(figsize=(width, 1.0))
    _plot_cnv_legend(ax, has_mirror=True)
    labels = [t for t in ax.texts if t.get_text() != "CNA"]
    rows = {round(t.get_position()[1], 6) for t in labels}
    plt.close(fig)
    return labels, rows


def test_wrap_groups_keeps_order_and_balances():
    rows = _wrap_groups([2, 4, 4, 6, 6, 8, 8, 2], 2)
    assert [i for row in rows for i in row] == list(range(8))
    assert len(rows) == 2
    assert _wrap_groups([2, 4], 5) == [[0], [1]]


@pytest.mark.parametrize("width", [60.0, 16.0, 8.0])
def test_cn_key_keeps_its_type_size_and_notation(width):
    labels, _ = _key_rows(width)
    assert {t.get_fontsize() for t in labels} == {FS_LEGEND}
    states = [t.get_text() for t in labels if t.get_text().startswith("(")]
    assert states and all(" " not in s and "|" not in s for s in states)


def test_cn_key_wraps_when_the_axes_is_narrow():
    _, wide = _key_rows(60.0)
    _, narrow = _key_rows(8.0)
    assert len(wide) == 1
    assert len(narrow) > 1


def test_rc_defaults_are_on_the_scale():
    for key in ("font.size", "xtick.labelsize", "ytick.labelsize", "legend.fontsize"):
        assert plt.rcParams[key] >= FS_ANNOT


def _line_gaps_pt(fig, labels, rotation):
    """Clearance across the text height between neighbouring parallel labels."""
    fig.canvas.draw()
    boxes = [t.get_window_extent() for t in labels if t.get_text()]
    x_pt = sorted((b.x0 + b.x1) / 2.0 * 72.0 / fig.dpi for b in boxes)
    sine = math.sin(math.radians(rotation))
    return [(b - a) * sine for a, b in zip(x_pt, x_pt[1:], strict=False)]


@pytest.mark.parametrize("width", [20.0, 12.0, 10.0])
def test_chromosome_names_keep_their_angle_and_clear_each_other(axis, width):
    fig, ax = plt.subplots(figsize=(width, 3))
    decorate_genome_axis(ax, axis, label_pos="top", rotation=60)
    labels = [t.label2 for t in ax.xaxis.get_major_ticks()]
    assert {t.get_fontsize() for t in labels} == {FS_TICK}
    assert {t.get_rotation() for t in labels} == {60.0}
    assert min(_line_gaps_pt(fig, labels, 60)) >= FS_TICK * 0.99
    plt.close(fig)


def test_narrow_axes_needs_a_second_row_of_names(axis):
    assert chrname_tiers(axis, 20.0) == 1
    assert chrname_tiers(axis, 4.0) == 2


def test_kept_mb_labels_are_a_line_apart(axis):
    fig, ax = plt.subplots(figsize=(12, 3))
    decorate_genome_axis(ax, axis, mb_ticks=True, plot_chrname=False, rotation=60)
    fig.canvas.draw()
    kept = [t for t in ax.get_xticklabels() if t.get_text()]
    assert len(kept) > len(axis.chrs) / 2
    # rotated labels clear each other when their baselines are a line apart
    x_pt = sorted(ax.transData.transform((t.get_position()[0], 0))[0] for t in kept)
    gap_pt = min(b - a for a, b in zip(x_pt, x_pt[1:], strict=False)) * 72 / fig.dpi
    assert gap_pt * math.sin(math.radians(60)) >= kept[0].get_fontsize() * 0.99
    plt.close(fig)
