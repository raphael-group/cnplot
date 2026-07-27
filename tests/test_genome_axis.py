"""Coordinate model: gaps, mapping, containment, and the mesh - no matplotlib."""

import numpy as np

from cnplot import GapKind


def test_chrom_order_follows_sizes(axis):
    assert axis.chrs[:3] == ["chr1", "chr2", "chr3"]
    assert axis.chrs[-1] == "chr22"
    assert len(axis.chrs) == 22


def test_gap_kinds_present(axis):
    kinds = {g.kind for g in axis.gaps}
    # T2T arm regions give interior (centromere) and leading gaps
    assert GapKind.INTERIOR in kinds
    assert GapKind.LEADING in kinds


def test_shrunk_gaps_are_zero_width(axis):
    assert all(g.width == 0 for g in axis.gaps)
    assert all(g.raw_width > 0 for g in axis.interior_gaps)


def test_kept_gaps_have_width(axis_keep):
    assert all(g.width > 0 for g in axis_keep.interior_gaps)
    assert not any(g.is_collapsed for g in axis_keep.interior_gaps)


def test_shrink_axis_is_contiguous(axis):
    segs = axis.segments
    assert all(
        abs(b.axis_start - a.axis_end) < 1e-6
        for a, b in zip(segs[:-1], segs[1:], strict=True)
    )


def test_all_bins_map(axis, sim):
    coords = axis.build_coordinates(sim.bins)
    assert coords.n_unmapped == 0
    assert coords.n_mapped == len(sim.bins)


def test_bins_land_inside_axis(axis, sim):
    coords = axis.build_coordinates(sim.bins)
    assert coords.positions.min() >= axis.ch_coords[0]
    assert coords.positions.max() <= axis.chr_end


def test_contain_matches_overlap_when_aligned(axis, sim):
    # seg.ucn rows are the regions themselves, so both rules map them all
    seg = sim.seg_ucn[sim.seg_ucn["SAMPLE"] == "S1"]
    assert axis.build_coordinates(seg, contain=True).n_mapped == len(seg)
    assert axis.build_coordinates(seg, contain=False).n_mapped == len(seg)


def test_grid_covers_axis_end_to_end(axis, sim):
    x_edges, col_bin_ids = axis.grid(sim.bins)
    assert x_edges[0] == axis.axis_start
    assert x_edges[-1] == axis.chr_end
    assert sum(b >= 0 for b in col_bin_ids) == len(sim.bins)
    assert np.all(np.diff(x_edges) > 0)


def test_excluded_chroms_drop(reference):
    from cnplot import GenomeAxis

    kept = GenomeAxis(*reference, excluded_chroms=["chr21", "chr22"])
    assert "chr21" not in kept.chrs and "chr22" not in kept.chrs


def test_region_bed_none_is_whole_genome(reference, sim):
    from cnplot import GenomeAxis

    _, chrom_sizes = reference
    whole = GenomeAxis(None, chrom_sizes, collapse_gaps=False)
    # every chromosome kept and no uncovered stretch: a plain full-genome axis
    assert len(whole.chrs) == 22
    assert whole.gaps == []
    assert whole.build_coordinates(sim.bins).n_unmapped == 0


def test_mb_ticks_labels_are_mb_and_drop_zero(axis_keep):
    import matplotlib.pyplot as plt

    from cnplot import decorate_genome_axis

    fig, ax = plt.subplots()
    decorate_genome_axis(ax, axis_keep, mb_ticks=True, mb_tick_step=20_000_000)
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels and all(lbl.isdigit() for lbl in labels)  # Mb integers
    assert "0" not in labels  # chromosome-start tick dropped
    # chromosome names drawn as off-axis text, one per chromosome
    names = [t.get_text() for t in ax.texts if t.get_text()]
    assert len(names) == len(axis_keep.chrs)
    plt.close(fig)
