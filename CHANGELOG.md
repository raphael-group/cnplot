# Changelog

## 0.1.0

First release of `cnplot`, the shared copy-number plotting library.

### cnplot package
- `GenomeAxis`: a genome coordinate transform built from a region BED and chromosome
  sizes, mapping any bin table by `(#CHR, START, END)`.
- Colormaps: `get_cn_cmap`, `get_ascn_cmap`, `get_baf_cmap`, `get_log2rdr_cmap`,
  `get_mixcn_cmap`, `get_multiclass_cmap`, `set_palette`.
- Integer copy-number profiles: `plot_cnv_profile` over the seg.ucn table.
- 1D genome scatter: `plot_scatter_1d` and the `plot_scatter_1d_multisample` builder,
  with `hue`/palette and an optional `expected_df` overlay.
- 2D RDR/BAF scatter: `plot_scatter_2d` (seaborn `JointGrid`) with copy-number landmarks.
- Single-cell heatmap: `plot_heatmap` / `plot_heatmap_cnp` with ordered categorical and
  posterior side `strips` and an optional value colorbar.
