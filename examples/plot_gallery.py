#!/usr/bin/env python3
# Render the README gallery figures into examples/, from the test simulator.
#
# Runpeng Luo (2026-07-23)
#
# Dependencies:
#   cnplot, and the tests/ simulator (tests/simulate.py + tests/data/)
#
# Usage:
#   python examples/plot_gallery.py
#
# Outputs:
#   examples/{profile,scatter_1d,scatter_2d,fcn_ab,heatmap_rdr,heatmap_baf}.png

import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))

from cnplot import (  # noqa: E402
    GenomeAxis,
    get_baf_cmap,
    get_log2rdr_cmap,
    get_mixcn_cmap,
    get_multiclass_cmap,
    make_row_spec,
    plot_cnv_profile,
    plot_heatmap_cnp,
    plot_scatter_1d_multisample,
    plot_scatter_2d,
)
from simulate import CELL_POSTERIORS, load_dataset, reference  # noqa: E402


def save(fig, name):
    path = os.path.join(HERE, f"{name}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(path)


def suptitle(fig, text, y=1.0):
    fig.suptitle(text, fontsize=14, fontweight="bold", y=y)


genome_axis = GenomeAxis(*reference())
sim = load_dataset(["S1", "S2"])

# color scatter points by the joint CNP, so a subclonal segment stays distinct
# from a clonal one at the same total copy number
palette = get_mixcn_cmap(
    list(sim.obs["cnp"].unique()),
    [float(sim.seg_ucn[f"u_{c}"].iloc[0]) for c in sim.clones],
)

# 1. integer copy-number profile
fig, (ax, lg) = plt.subplots(2, 1, figsize=(12, 3), height_ratios=[3, 1])
plot_cnv_profile(
    ax, sim.seg_ucn, genome_axis, sample_id="S1", ax_leg=lg, show_pi_viol=False
)
suptitle(fig, "Integer copy-number profile", y=1.18)
save(fig, "profile")

# 2. genome-wide RDR + BAF, multi-sample, over the shared profile
row_specs = [
    make_row_spec("RD", ylabel="RDR", ylim=(0, 3), href=1.0),
    make_row_spec("BAF", ylabel="mhBAF", ylim=(-0.05, 1.05), href=0.5),
]
fig = plot_scatter_1d_multisample(
    sim.obs,
    genome_axis,
    row_specs,
    expected_df=sim.expected_1d,
    hue="cnp",
    palette=palette,
    seg_df=sim.seg_ucn,
    mb_ticks=True,
    show_pi_viol=False,
)
suptitle(fig, "Genome-wide RDR and BAF (multi-sample)", y=1.0)
save(fig, "scatter_1d")

# 3. RDR vs BAF joint scatter, with copy-number landmarks
grid = plot_scatter_2d(
    sim.obs,
    "BAF",
    "RD",
    expected_df=sim.expected_2d,
    group="S1",
    hue="cnp",
    palette=palette,
    xlim=(0, 1),
    ylim=(0, 3),
    xlabel="mhBAF",
    ylabel="RDR",
)
suptitle(grid.figure, "RDR vs BAF joint scatter", y=1.0)
save(grid.figure, "scatter_2d")

# 4. allele-specific fractional copy number: FCN-A major, FCN-B minor mirrored
# (FCN-A / FCN-B and their exp_ overlays are carried by the bbc.ucn / expected_1d)
rows_ab = [
    make_row_spec("FCN-A", ylabel="FCN-A", ylim=(0, 3.5), href=1.0),
    make_row_spec("FCN-B", ylabel="FCN-B", ylim=(0, 3.5), href=1.0, reverse_y=True),
]
fig = plot_scatter_1d_multisample(
    sim.obs,
    genome_axis,
    rows_ab,
    expected_df=sim.expected_1d,
    hue="cnp",
    palette=palette,
    seg_df=sim.seg_ucn,
    mb_ticks=True,
    show_pi_viol=False,
)
suptitle(fig, "Allele-specific fractional copy number (FCN-A / FCN-B)", y=1.0)
save(fig, "fcn_ab")

# 5 & 6. single-cell heatmap page: mesh + colorbar + strips + CN profile + legend
cells = pd.read_table(CELL_POSTERIORS)
cell_labels = cells["clone"].to_numpy()
celltype = np.where(cell_labels == "normal", "normal", "tumor")
clone_cmap = get_multiclass_cmap({"clone": cell_labels}, "clone")["clone"]
order = list(clone_cmap)
post = cells[order].to_numpy()
post_prop = {v: float((cell_labels == v).mean()) for v in order}
# posterior distribution strip closest to the heatmap, then a categorical strip
strips = [
    {
        "name": "Copy-typing",
        "matrix": post,
        "order": order,
        "cmap": clone_cmap,
        "props": post_prop,
    },
    {"name": "cell_type", "values": celltype, "display_name": "Cell-type"},
]
rdr_cmap, rdr_norm, rdr_ticks = get_log2rdr_cmap()
baf_cmap, baf_norm, baf_ticks = get_baf_cmap()

for matrix, cmap, norm, label, ticks, out in [
    (sim.heatmap_rdr, rdr_cmap, rdr_norm, "RDR", rdr_ticks, "heatmap_rdr"),
    (sim.heatmap_baf, baf_cmap, baf_norm, "BAF", baf_ticks, "heatmap_baf"),
]:
    fig = plot_heatmap_cnp(
        matrix,
        sim.bins,
        genome_axis,
        sim.seg_ucn,
        sample_id="S1",
        title=f"Single-cell copy-number heatmap ({label})",
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        cbar_label=label,
        cbar_ticks=ticks,
        strips=strips,
        mb_ticks=True,
        profile_kwargs={"show_pi_viol": False},
    )
    save(fig, out)
