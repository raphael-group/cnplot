# Plotting inventory across three CNV codebases

Cross-reference of plotting code in three sibling repos under
`proj-copy_number_variation/`:

- **HATCHet3** - `hatchet-long-read/src/hatchet/plot/`
- **Copytyping** - `Copytyping-dev/src/copytyping/plot/`
- **UGP** (Universal-Genotyping-Pipeline) - `universal-genotyping/workflow/scripts/plot_utils.py`
  (+ `hmm_genotype_utils.plot_qc`)

Roles differ by pipeline stage:

- **UGP** = preprocessing QC (read-depth correction, segmentation, SNP/allele QC). Stops
  before CN calling, so it has no integer-CN profile/legend/colormap.
- **HATCHet3** = bulk allele/clone-specific CN calling + plotting (CNP-string driven).
- **Copytyping** = single-cell / spatial CN typing (cell-matrix driven), plus heatmaps and
  Visium.

---

## 1. Per-codebase inventory

### HATCHet3 (`src/hatchet/plot/`)

| File | Public functions | Role |
|---|---|---|
| plot_cn.py | `run` | `plot-cn` CLI; builds FCN/BAF + FCN-A/B 1D and 2D per solution |
| plot_cnp_panel.py | `run` | `plot-panel` CLI; multi-solution CNP grid |
| plot_1d2d.py | `plot_1d`, `plot_2d`, `plot_clusters`, `plot_rdr_baf`, `get_transparency`, `get_abs_positions_ignore_gap/_keep_gap`, `_is_multimodal` | 1D/2D scatter primitives + cluster diagnostics |
| plot_cn_utils.py | `plot_cnv_profile`, `plot_cnv_legend`, `plot_ascn_profile`, `plot_ascn_legend`, `get_cn_colors`, `get_ascn_colors`, `build_cnp_palette`, `_distinct_subclonal_colors` | Integer-CN profiles, legends, colormaps |
| plot_utils.py | `plot_combined_1d`, `_row_def`, `get_expected_baf_fcn`, `override_solution`, `load_gammas`, `set_palette` | Combined 1D builder + compute-cn glue |
| plot_cn_tree.py | `render_cnt_tree` (+ helpers) | Clone tree (cnt_cd) |
| plot_pool.py | `plot_pool_cnp` | Solution-pool CNP panels |
| plot_scaling_2d.py | `plot_scaling_2d` | Scaling-factor inference viz |
| plot_common.py | `plot_summary_pdf`, `plot_bars` | Purity/ploidy bar summaries |

### Copytyping (`src/copytyping/plot/`)

| File | Public functions | Role |
|---|---|---|
| plot_copynumber.py | `plot_cnv_profile`, `plot_cnv_legend`, `plot_ascn_profile`, `plot_ascn_legend`, `get_cn_colors`, `get_ascn_colors`, `cnp_has_mirror` (+ layout helpers) | Integer-CN profiles/legends/colors |
| plot_scatter_1d.py | `plot_scatter_1d_pseudobulk`, `plot_rdr_baf_1d_pseudobulk`, `_build_ch_boundary`, `_merge_exp_lines` | 1D genome scatter |
| plot_scatter_2d.py | `plot_rdr_baf_2d_pseudobulk`, `plot_scatter_2d_per_cell`, `_annotate_cn_landmarks` | 2D scatter |
| plot_heatmap.py | `plot_heatmap`, `plot_cnv_heatmap`, `plot_label_strips`, `draw_label_legends`, `prepare_rdr/baf/pi_gk`, `_aggregate_columns`, `_row_layout`, `_mode` | Single-cell CN heatmaps |
| plot_spatial.py | `plot_visium_all/panel/loh_baf/iters`, `build_visium_slices`, `compute_loh_baf`, `blend_purity_rgba`, `set_label_colors`, `build_legend` | Visium/spatial |
| plot_common.py | `build_label_colors`, `build_label_color_maps`, `build_categorical_color_map`, `make_baf_cmap`, `build_wl_coords`, `plot_loss`, `_clone_order_key`, `_is_normal_like`, ... | Label color machinery + training loss |

### UGP (`workflow/scripts/plot_utils.py`, `hmm_genotype_utils.py`)

| Function | Role |
|---|---|
| `plot_rd_2d_kde` | RD-correction QC: before/after 2D KDE vs GC/mappability/repliseq |
| `plot_rd_1d_scatter` | RD-correction QC: before/after 1D genome scatter |
| `plot_1d_sample`, `plot_1d_multi_sample` | Generic 1D genome-wide scatter |
| `plot_rdr_baf` | Genome-wide RDR+BAF, one page per sample, two rows |
| `plot_segmentation_qc` | Segmentation QC histograms (combine_counts output) |
| `plot_snp_depth_histogram`, `plot_allele_freqs` | SNP total-depth / allele-frequency QC |
| `hmm_genotype_utils.plot_qc` | Per-chrom genotyping QC (alpha/beta/total, het calls) |
| `_genome_coords`, `_add_chrom_decorations`, `_parse_bed_by_chr`, `_hist_with_stats`, `_seg_gene_counts`, `_extract_col`, `_plot_cov_panel` | Helpers (coords, chrom decoration, histograms) |

---

## 2. Union by capability

Grouped by what is drawn; cells mark which repo implements it.

| Capability | HATCHet3 | Copytyping | UGP |
|---|---|---|---|
| Integer-CN profile (per-clone bars) | `plot_cnv_profile` | `plot_cnv_profile` | - |
| Allele-specific CN profile (A/B) | `plot_ascn_profile` | `plot_ascn_profile` | - |
| CN legend (total / allele) | `plot_cnv_legend` / `plot_ascn_legend` | `plot_cnv_legend` / `plot_ascn_legend` | - |
| Integer-CN colormap | `get_cn_colors` / `get_ascn_colors` | `get_cn_colors` / `get_ascn_colors` | - |
| Clonal/subclonal scatter palette | `build_cnp_palette`, `_distinct_subclonal_colors` | - | - |
| 1D genome scatter primitive | `plot_1d` | `plot_scatter_1d_pseudobulk` | `plot_1d_sample` / `plot_1d_multi_sample` |
| Genome RDR+BAF, per-sample 2-row | `plot_combined_1d` | `plot_rdr_baf_1d_pseudobulk` | `plot_rdr_baf` |
| 2D BAF-vs-RDR joint scatter | `plot_2d` | `plot_rdr_baf_2d_pseudobulk` | - |
| Per-cluster 2D diagnostics | `plot_clusters` | `plot_scatter_2d_per_cell` | - |
| Bin -> genome x-coord mapping | `get_abs_positions_ignore_gap/_keep_gap` | `_build_ch_boundary`, `build_wl_coords` | `_genome_coords` |
| Chromosome-boundary decoration | inline | `_draw_chr_boundaries`, `_decorate_cnp_xaxis` | `_add_chrom_decorations` |
| Expected-value line merging | inline in `plot_1d` | `_merge_exp_lines` | inline |
| CN-state centroid labels | inline in `plot_2d` | `_annotate_cn_landmarks` | - |
| Categorical / cluster palette | `set_palette` | `build_label_colors`, `build_categorical_color_map` | (label colors inline) |
| Clone tree | `render_cnt_tree` | - | - |
| Solution-pool CNP panels | `plot_pool_cnp` | - | - |
| Scaling-factor viz | `plot_scaling_2d` | - | - |
| Purity/ploidy summary bars | `plot_summary_pdf`, `plot_bars` | - | - |
| Single-cell CN heatmap | - | `plot_cnv_heatmap` (+ suite) | - |
| Spatial / Visium | - | `plot_visium_*` (suite) | - |
| BAF diverging colormap | - | `make_baf_cmap` | - |
| Training-loss curve | - | `plot_loss` | - |
| RD-correction QC (KDE / 1D) | - | - | `plot_rd_2d_kde`, `plot_rd_1d_scatter` |
| Segmentation QC histograms | - | - | `plot_segmentation_qc` |
| SNP depth / allele-freq QC | - | - | `plot_snp_depth_histogram`, `plot_allele_freqs` |
| Genotyping QC | - | - | `hmm_genotype_utils.plot_qc` |

---

## 3. Redundancy

### Hard duplication (identical or lightly drifted; dedup candidates)

HATCHet3 <-> Copytyping share the integer-CN plotting stack:

| Function | Status |
|---|---|
| `get_cn_colors` | **byte-identical** (both use `khaki` for (2,1)) |
| `get_ascn_colors` | identical / near-identical |
| `plot_cnv_profile` | same purpose, drifted signature (`bin_info,regions` vs `seg_cnprofile,wl_segments`) |
| `plot_cnv_legend` | drifted (Copytyping adds `has_mirror` + `cnp_has_mirror`/`_draw_mirror_swatch`) |
| `plot_ascn_profile` | same, drifted signature |
| `plot_ascn_legend` | same |

These are copy-pasted and drift independently -> maintenance hazard. UGP does not
participate (no CN calling).

### Conceptual duplication (same idea, diverged implementations)

Present in all three or two, but shaped to each data model:

| Concept | HATCHet3 | Copytyping | UGP |
|---|---|---|---|
| Genome RDR+BAF per-sample 2-row plot | `plot_combined_1d` | `plot_rdr_baf_1d_pseudobulk` | `plot_rdr_baf` |
| 1D genome scatter primitive | `plot_1d` | `plot_scatter_1d_pseudobulk` | `plot_1d_sample` |
| Bin -> genome-coord mapping | `get_abs_positions_*` | `_build_ch_boundary`/`build_wl_coords` | `_genome_coords` |
| 2D BAF-vs-RDR joint scatter | `plot_2d` | `plot_rdr_baf_2d_pseudobulk` | - |
| Chrom-boundary decoration | inline | `_draw_chr_boundaries` | `_add_chrom_decorations` |

The genome RDR+BAF two-row plot and the bin->coord mapping are the strongest cross-repo
concept duplication (all three). HATCHet is CNP-string driven, Copytyping is cell-matrix
driven, UGP is QC/preprocessing driven, so unifying needs a refactor, not a copy.

### Not redundant (single-repo)

- **HATCHet3 only**: `build_cnp_palette`, `_distinct_subclonal_colors`, `plot_combined_1d`/`_row_def`,
  `get_expected_baf_fcn`, `override_solution`, `load_gammas`, `plot_scaling_2d`, `plot_pool_cnp`,
  `plot_summary_pdf`, `render_cnt_tree`.
- **Copytyping only**: heatmap suite, Visium/spatial suite, label-color machinery, `make_baf_cmap`,
  `plot_loss`, mirror helpers.
- **UGP only**: RD-correction QC, segmentation QC, SNP/allele-freq QC, genotyping QC.

---

## 4. Recommendations

1. **Extract a shared `cn_colors` + CNP-profile module** for the 6 HATCHet3<->Copytyping
   duplicated functions (start with the identical `get_cn_colors`/`get_ascn_colors`). A single
   source prevents silent drift; the `khaki` change happens to match today but the next edit
   to one won't propagate.
2. **Leave the 1D/2D scatter concept duplication as-is** unless a shared plotting library is
   worth building - the three implementations are bound to different data models (CNP string /
   cell matrix / QC arrays).
3. **UGP stays separate** - it is preprocessing QC and shares only the genome-coordinate and
   1D-scatter concepts, not the CN-calling visuals.
