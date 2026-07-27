# Changelog

## 0.1.1

Backward-compatible additions; existing calls and output are unchanged.

### cnplot package
- `decorate_genome_axis` gains `mb_ticks` / `mb_tick_step`: label x in Mb every
  `mb_tick_step` bp within each chromosome (reset per chromosome), with chromosome names
  drawn as compact off-axis text - the denser mode for genome-wide scatter QC. Plumbed
  through `plot_scatter_1d`, `plot_scatter_1d_multisample`, `plot_heatmap`, and
  `plot_heatmap_cnp`. Default off, 50Mb step; the Mb labels are rotated by `rotation`
  (default 60 degrees) so the dense ticks do not overlap.
- Version lives in the repo-root `VERSION` file: `pyproject.toml` builds from it
  (`dynamic version = {file = "VERSION"}`) and `cnplot.__version__` reads it back, falling
  back to the installed metadata outside a source checkout. One edit per bump.
- `GenomeAxis(region_bed=None)`: build one full-length region per chromosome from
  `chrom_sizes`, i.e. a plain whole-genome axis with no gaps, for callers that only have a
  sizes file.

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
