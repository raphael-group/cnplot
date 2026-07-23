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
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))

from cnplot.cnplot_2d import plot_scatter_2d  # noqa: E402
from cnplot.cnplot_colormap import (  # noqa: E402
    build_label_cmaps,
    build_mixture_cn_cmap,
    get_baf_cmap,
)
from cnplot.cnplot_figures import (  # noqa: E402
    make_row_spec,
    plot_heatmap_cnp,
    plot_scatter_1d_multisample,
)
from cnplot.cnplot_genome_axis import GenomeAxis  # noqa: E402
from cnplot.cnplot_intcnp import plot_cnv_profile  # noqa: E402
from simulate import reference, simulate  # noqa: E402


def save(fig, name):
    path = os.path.join(HERE, f"{name}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(path)


def suptitle(fig, text, y=1.0):
    fig.suptitle(text, fontsize=14, fontweight="bold", y=y)


axis = GenomeAxis(*reference())
sim = simulate(seed=0, samples=["S1", "S2"])

# color scatter points by the joint CNP, so a subclonal segment stays distinct
# from a clonal one at the same total copy number
palette = build_mixture_cn_cmap(
    list(sim.obs["cnp"].unique()),
    [float(sim.seg_ucn[f"u_{c}"].iloc[0]) for c in sim.clones],
)

# 1. integer copy-number profile
fig, (ax, lg) = plt.subplots(2, 1, figsize=(12, 3), height_ratios=[3, 1])
plot_cnv_profile(ax, sim.seg_ucn, axis, sample_id="S1", ax_leg=lg)
suptitle(fig, "Integer copy-number profile", y=1.18)
save(fig, "profile")

# 2. genome-wide RDR + BAF, multi-sample, over the shared profile
rows = [
    make_row_spec("RD", ylabel="RDR", ylim=(0, 3), href=1.0),
    make_row_spec("BAF", ylabel="mhBAF", ylim=(-0.05, 1.05), href=0.5),
]
fig = plot_scatter_1d_multisample(
    sim.obs,
    axis,
    rows,
    expected_df=sim.expected_1d,
    hue="cnp",
    palette=palette,
    seg_df=sim.seg_ucn,
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
obs_ab = sim.obs.copy()
obs_ab["FCN-A"] = 2 * obs_ab["RD"] * (1 - obs_ab["BAF"])
obs_ab["FCN-B"] = 2 * obs_ab["RD"] * obs_ab["BAF"]
exp_ab = sim.expected_1d[["#CHR", "START", "END"]].copy()
for s in sim.samples:
    rd, baf = sim.expected_1d[f"exp_RD_{s}"], sim.expected_1d[f"exp_BAF_{s}"]
    exp_ab[f"exp_FCN-A_{s}"] = 2 * rd * (1 - baf)
    exp_ab[f"exp_FCN-B_{s}"] = 2 * rd * baf
rows_ab = [
    make_row_spec("FCN-A", ylabel="FCN-A", ylim=(0, 3.5), href=1.0),
    make_row_spec("FCN-B", ylabel="FCN-B", ylim=(0, 3.5), href=1.0, reverse_y=True),
]
fig = plot_scatter_1d_multisample(
    obs_ab,
    axis,
    rows_ab,
    expected_df=exp_ab,
    hue="cnp",
    palette=palette,
    seg_df=sim.seg_ucn,
)
suptitle(fig, "Allele-specific fractional copy number (FCN-A / FCN-B)", y=1.0)
save(fig, "fcn_ab")

# 5 & 6. single-cell heatmap page: mesh + colorbar + strips + CN profile + legend
celltype = np.where(sim.heatmap_labels == "normal", "normal", "tumor")
clone_cmap = build_label_cmaps({"clone": sim.heatmap_labels}, "clone")["clone"]
order = list(clone_cmap)
rng = np.random.default_rng(1)
post = np.zeros((len(sim.heatmap_labels), len(order)))
for i, lab in enumerate(sim.heatmap_labels):
    p = rng.dirichlet([1.0] * len(order))
    p[order.index(lab)] += 2.0
    post[i] = p / p.sum()
post_prop = {v: float((sim.heatmap_labels == v).mean()) for v in order}
dist_strip = ("Copy-typing", post, order, clone_cmap, post_prop)
strips = {"cell_type": celltype}
names = {"cell_type": "Cell-type", "Copy-typing": "Copy-typing"}
baf_cmap, baf_norm = get_baf_cmap()

for matrix, cmap, norm, label, ticks, out in [
    (
        sim.heatmap_rdr,
        "coolwarm",
        TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
        "RDR",
        [-1, -0.5, 0, 0.5, 1],
        "heatmap_rdr",
    ),
    (
        sim.heatmap_baf,
        baf_cmap,
        baf_norm,
        "BAF",
        [0, 0.25, 0.5, 0.75, 1],
        "heatmap_baf",
    ),
]:
    fig = plot_heatmap_cnp(
        matrix,
        sim.bins,
        axis,
        sim.seg_ucn,
        sample_id="S1",
        title=f"Single-cell copy-number heatmap ({label})",
        row_labels=sim.heatmap_labels,
        cmap=cmap,
        norm=norm,
        cbar_label=label,
        cbar_ticks=ticks,
        strip_label_map=strips,
        dist_strip=dist_strip,
        display_names=names,
    )
    save(fig, out)
