#!/usr/bin/env python3
# Deterministic toy copy-number data in every cnplot input format, built from a
# real 22-chromosome CNP profile.
#
# Runpeng Luo (2026-07-23)
#
# Dependencies:
#   numpy, pandas (no cnplot import at module load; produces the plain tables
#   cnplot consumes)
#
# Usage:
#   from tests.simulate import simulate, reference
#   sim = simulate(seed=0)          # bundle of tables + arrays
#   bed, sizes = reference()        # vendored T2T region BED + chrom sizes
#
# Inputs (vendored under tests/data/)
#   sample.seg.ucn.tsv        HATCHet seg.ucn: 69 segments over 22 autosomes,
#                             normal + clone1 + clone2, states 1|0..2|2.
#   T2T-CHM13v2.0.regions.bed arm-level whitelist the profile was called against.
#   T2T-CHM13v2.0.sizes       chromosome lengths.
# Parameters
#   The seg.ucn is the ground-truth CNP. Fine bins are cut inside each segment;
#   observed 1D/2D values are the proportion-weighted clone mixture (a bulk
#   sample) plus Gaussian noise, and heatmap rows are per-cell single-clone
#   truths plus noise. Both derive from the same profile so expected overlays
#   line up.
# Notes/References:
#   Real HATCHet seg.ucn profile (SAMPLE anonymised to "sample") over the
#   T2T-CHM13v2.0 reference, recast into cnplot's DataFrame contracts (seg.ucn
#   columns, exp_<col>_<group>).

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
REGION_BED = os.path.join(DATA_DIR, "T2T-CHM13v2.0.regions.bed")
CHROM_SIZES = os.path.join(DATA_DIR, "T2T-CHM13v2.0.sizes")
SEG_UCN = os.path.join(DATA_DIR, "sample.seg.ucn.tsv")

BIN_SIZE = 2_000_000
SIGMA_RDR = 0.12
SIGMA_BAF = 0.05
NUM_CELLS = 120


@dataclass
class Sim:
    """A bundle of one dataset in every cnplot input format.

    Attributes:
        bins: Bin table with "#CHR", "START", "END", one row per bin.
        seg_ucn: seg.ucn profile: "SAMPLE", "#CHR/START/END", ``cn_<clone>``,
            ``u_<clone>``, and a "PI_VIOL" flag. One block per sample.
        obs: Per-bin observations for the scatter plots: the bin columns plus
            "SAMPLE", "RD", "BAF", "log2RDR", "state" (dominant clone's ``a|b``),
            "cnp" (joint ``;``-joined CNP over all clones, for
            :func:`~cnplot.cnplot_colormap.get_mixcn_cmap`), "CLUSTER"
            (segment id), and "pass_qc".
        expected_1d: Segment-resolution expected values, "#CHR/START/END" plus
            one ``exp_RD_<sample>`` / ``exp_BAF_<sample>`` column per sample.
        expected_2d: seg.ucn layout plus per-row "exp_BAF" / "exp_RD" landmark
            coordinates.
        clones: Clone names, normal first.
        samples: Sample ids present in the multi-sample tables.
        heatmap_rdr: (n_cells, n_bins) per-cell log2RDR.
        heatmap_baf: (n_cells, n_bins) per-cell BAF.
        heatmap_labels: (n_cells,) clone label per row, grouped bottom to top.
    """

    bins: pd.DataFrame
    seg_ucn: pd.DataFrame
    obs: pd.DataFrame
    expected_1d: pd.DataFrame
    expected_2d: pd.DataFrame
    clones: list
    samples: list
    heatmap_rdr: np.ndarray
    heatmap_baf: np.ndarray
    heatmap_labels: np.ndarray


def reference() -> tuple:
    """Paths to the vendored T2T region BED and chromosome-sizes file.

    Returns:
        (region_bed_path, chrom_sizes_path).
    """
    return REGION_BED, CHROM_SIZES


def _clone_names(profile: pd.DataFrame) -> list:
    """Clone names from the ``cn_`` columns, normal first then clone1, clone2.

    Args:
        profile: seg.ucn table.

    Returns:
        Clone names in stacking order.
    """
    names = [c[3:] for c in profile.columns if c.startswith("cn_")]

    def key(n):
        if n == "normal":
            return (0, 0)
        return (1, int(n[5:])) if n.startswith("clone") else (2, 0)

    return sorted(names, key=key)


def _states(profile: pd.DataFrame, clones: list) -> np.ndarray:
    """Parse ``cn_<clone>`` strings into an (n_seg, n_clones, 2) int array.

    Args:
        profile: seg.ucn table.
        clones: Clone names to read.

    Returns:
        Per-segment (a, b) pairs per clone.
    """
    out = np.empty((len(profile), len(clones), 2), dtype=int)
    for k, c in enumerate(clones):
        parts = profile[f"cn_{c}"].astype(str).str.split("|", expand=True)
        out[:, k, :] = parts.astype(int).to_numpy()
    return out


def _mixture(states: np.ndarray, props: np.ndarray) -> tuple:
    """Proportion-weighted (RD, BAF) per segment.

    Args:
        states: (n_seg, n_clones, 2) copy numbers.
        props: (n_clones,) proportions.

    Returns:
        (rd, baf): expected RD (mixture total CN over 2) and BAF (mixture B over
        total), one value per segment.
    """
    a_mix = (states[:, :, 0] * props).sum(axis=1)
    b_mix = (states[:, :, 1] * props).sum(axis=1)
    c_mix = a_mix + b_mix
    return c_mix / 2.0, np.where(c_mix > 0, b_mix / c_mix, 0.5)


def _build_bins(profile: pd.DataFrame) -> tuple:
    """Cut each segment into ``BIN_SIZE`` bins.

    Args:
        profile: Single-sample seg.ucn table.

    Returns:
        (bins, seg_idx): the "#CHR/START/END" bin table and the segment id per
        bin, indexing ``profile``.
    """
    rows, seg_idx = [], []
    for si, seg in profile.reset_index(drop=True).iterrows():
        start, end = int(seg["START"]), int(seg["END"])
        n = max(1, round((end - start) / BIN_SIZE))
        edges = np.linspace(start, end, n + 1, dtype=np.int64)
        for i in range(n):
            rows.append((seg["#CHR"], int(edges[i]), int(edges[i + 1])))
            seg_idx.append(si)
    bins = pd.DataFrame(rows, columns=["#CHR", "START", "END"])
    return bins, np.asarray(seg_idx)


def simulate(seed: int = 0, samples: list | None = None) -> Sim:
    """Generate one dataset in every cnplot input format from the real profile.

    Args:
        seed: Seed for the observation and cell noise.
        samples: Sample ids for the multi-sample tables. The first keeps the
            profile's own proportions; others perturb them so a shared profile
            meets differing mixtures. None keeps the single real sample.

    Returns:
        A :class:`Sim` bundle.
    """
    rng = np.random.default_rng(seed)
    profile = pd.read_table(SEG_UCN)
    clones = _clone_names(profile)
    base_sample = str(profile["SAMPLE"].iloc[0])
    if samples is None:
        samples = [base_sample]

    base_props = {c: float(profile[f"u_{c}"].iloc[0]) for c in clones}
    tumor = [c for c in clones if c != "normal"]
    props_by_sample = {samples[0]: dict(base_props)}
    for i, s in enumerate(samples[1:], start=1):
        shift = min(0.1 * i, base_props.get("normal", 0.2) - 0.02)
        p = dict(base_props)
        if "normal" in p and tumor:
            p["normal"] -= shift
            p[tumor[0]] += shift
        props_by_sample[s] = p

    states = _states(profile, clones)  # (n_seg, n_clones, 2)
    dominant = max(tumor, key=lambda c: base_props[c]) if tumor else clones[0]
    dom_k = clones.index(dominant)

    bins, seg_idx = _build_bins(profile)
    coords = bins.copy()

    # per-bin observations: bulk mixture (base props) plus noise
    base_props_arr = np.array([base_props[c] for c in clones])
    rd_seg, baf_seg = _mixture(states, base_props_arr)
    rd_bin, baf_bin = rd_seg[seg_idx], baf_seg[seg_idx]
    obs_blocks = []
    for sample in samples:
        df = coords.copy()
        df["SAMPLE"] = sample
        df["RD"] = np.clip(rd_bin + rng.normal(0, SIGMA_RDR, len(df)), 0, None)
        df["BAF"] = np.clip(baf_bin + rng.normal(0, SIGMA_BAF, len(df)), 0, 1)
        df["log2RDR"] = np.log2(np.maximum(df["RD"], 1e-6))
        df["state"] = [f"{states[s, dom_k, 0]}|{states[s, dom_k, 1]}" for s in seg_idx]
        # joint CNP string over all clones, for get_mixcn_cmap (keeps a
        # subclonal segment distinct from a clonal one at the same total CN)
        cnp_by_seg = [
            ";".join(
                f"{states[si, k, 0]}|{states[si, k, 1]}" for k in range(len(clones))
            )
            for si in range(len(profile))
        ]
        df["cnp"] = [cnp_by_seg[s] for s in seg_idx]
        df["CLUSTER"] = seg_idx
        df["pass_qc"] = rng.random(len(df)) > 0.15
        obs_blocks.append(df)
    obs = pd.concat(obs_blocks, ignore_index=True)

    # multi-sample seg.ucn and expected tables (shared CN, per-sample proportions)
    seg_blocks = []
    exp1d = profile[["#CHR", "START", "END"]].copy()
    a = states[..., 0]
    b = states[..., 1]
    for sample in samples:
        block = profile.copy()
        block["SAMPLE"] = sample
        for c in clones:
            block[f"u_{c}"] = props_by_sample[sample][c]
        p = np.array([props_by_sample[sample][c] for c in clones])
        c_mix = (a * p).sum(1) + (b * p).sum(1)
        # PI_VIOL: flag segments whose mixture total CN is far from an integer
        block["PI_VIOL"] = np.abs(c_mix - np.round(c_mix)) > 0.15
        seg_blocks.append(block)
        rd_s, baf_s = _mixture(states, p)
        exp1d[f"exp_RD_{sample}"] = rd_s
        exp1d[f"exp_BAF_{sample}"] = baf_s
    seg_ucn = pd.concat(seg_blocks, ignore_index=True)

    expected_2d = seg_ucn.copy()
    e_rd, e_baf = [], []
    n_seg = len(profile)
    for i in range(len(expected_2d)):
        sample = expected_2d["SAMPLE"].iloc[i]
        p = np.array([props_by_sample[sample][c] for c in clones])
        rd_s, baf_s = _mixture(states[i % n_seg][None], p)
        e_rd.append(float(rd_s[0]))
        e_baf.append(float(baf_s[0]))
    expected_2d["exp_RD"] = e_rd
    expected_2d["exp_BAF"] = e_baf

    # per-cell heatmap: single-clone truth plus noise. Rows go bottom to top in
    # pcolormesh, so reverse the clone order to read normal, clone1, ... top down.
    counts = {c: max(1, int(round(NUM_CELLS * base_props[c]))) for c in clones}
    cell_labels = np.concatenate([[c] * counts[c] for c in reversed(clones)])
    rdr_rows, baf_rows = [], []
    for c in cell_labels:
        k = clones.index(c)
        aa = states[seg_idx, k, 0].astype(float)
        bb = states[seg_idx, k, 1].astype(float)
        tot = aa + bb
        rdr_rows.append(
            np.log2(np.maximum(tot / 2.0, 1e-6)) + rng.normal(0, SIGMA_RDR, len(tot))
        )
        baf_rows.append(
            np.clip(
                np.where(tot > 0, bb / tot, 0.5) + rng.normal(0, SIGMA_BAF, len(tot)),
                0,
                1,
            )
        )

    return Sim(
        bins=coords.reset_index(drop=True),
        seg_ucn=seg_ucn,
        obs=obs,
        expected_1d=exp1d,
        expected_2d=expected_2d,
        clones=clones,
        samples=list(samples),
        heatmap_rdr=np.vstack(rdr_rows),
        heatmap_baf=np.vstack(baf_rows),
        heatmap_labels=cell_labels,
    )


def _demo() -> None:
    """Render a couple of demo figures into the project .tmp, as a sanity check."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from cnplot import (
        GenomeAxis,
        make_row_spec,
        plot_cnv_profile,
        plot_scatter_1d_multisample,
    )

    out = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".tmp"
    )
    os.makedirs(out, exist_ok=True)
    axis = GenomeAxis(*reference())
    sim = simulate(seed=0, samples=["S1", "S2"])
    print(f"{len(sim.bins)} bins, {len(sim.seg_ucn) // len(sim.samples)} segments")

    fig, ax = plt.subplots(figsize=(16, 2))
    plot_cnv_profile(ax, sim.seg_ucn, axis, sample_id="S1")
    fig.savefig(os.path.join(out, "sim_profile.png"), dpi=90, bbox_inches="tight")
    plt.close(fig)

    rows = [make_row_spec("RD", ylabel="RDR", href=1.0), make_row_spec("BAF", href=0.5)]
    fig = plot_scatter_1d_multisample(
        sim.obs, axis, rows, expected_df=sim.expected_1d, seg_df=sim.seg_ucn
    )
    fig.savefig(os.path.join(out, "sim_multisample.png"), dpi=90, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote demo figures to {out}")


if __name__ == "__main__":
    _demo()
