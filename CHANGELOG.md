# Changelog

All notable changes to `cnplot` are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/); versions follow semantic versioning.

## 0.1.0 - Initial release

First standalone release. `cnplot` factors the shared copy-number plotting code out of
HATCHet3, Copytyping, and UGP into one PyPI-installable library (numpy, pandas, matplotlib,
seaborn, adjustText; no scipy). Two tiers: primitives take an `ax` and return None; figure
builders own the layout and return a `Figure`.

### Added

- Genome coordinate model - `GenomeAxis`, a coordinate transform built from a region BED and
  chromosome sizes that maps any bin table by `(#CHR, START, END)`, so one axis serves every
  sample, bin resolution, and panel. `build_coordinates(contain=)`, `grid()` for meshes,
  `Segment` / `Gap` / `BinCoords` records, and `read_bed` / `read_chr_sizes` readers.
- Colormaps (`cnplot_colormap`) - `get_cn_cmap`, `get_ascn_cmap`, `get_baf_cmap`,
  `get_log2rdr_cmap` (all cmap builders return `(cmap, norm, ticks)`), `get_mixcn_cmap` for
  joint-CNP scatter coloring, `get_multiclass_cmap` for shared-palette annotation strips, and
  `set_palette`. Naming convention: built-in `get_*`, data-derived `build_*`, mutators `set_*`.
- Integer copy-number profiles (`cnplot_intcnp`) - `plot_cnv_profile` driven by the seg.ucn
  table (`cn_<clone>`, `u_<clone>`, optional `PI_VIOL` / `SAMPLE`); draws its own legend via
  `ax_leg=`, with a legacy allele-specific pair kept for HATCHet.
- 1D genome scatter (`cnplot_1d`) - `plot_scatter_1d` primitive: caller-supplied `ycol`,
  `hue` + palette, optional `expected_df` overlay (segment or bin resolution), and a
  kept/filtered two-color mode.
- 2D RDR/BAF scatter (`cnplot_2d`) - `plot_scatter_2d` (seaborn `JointGrid`) with seg.ucn-style
  landmarks, per-sample clonal/subclonal labelling, proportion legend, and `get_landmarks`.
- Single-cell heatmap (`cnplot_heatmap`) - `plot_heatmap` places a reduced `(n_rows, n_bins)`
  matrix through `GenomeAxis.grid`, masking filler columns so nothing bleeds across a gap. One
  ordered `strips` list draws categorical and posterior-distribution side strips
  (`plot_column_strips`, `plot_strip_legend`), plus an optional value colorbar.
- Figure builders (`cnplot_figures`) - `plot_scatter_1d_multisample` (one row per group over a
  shared profile), `plot_heatmap_cnp` (heatmap + colorbar + strips + CN profile page), and
  `make_row_spec`.
- Shared helpers (`cnplot_utils`) - axis decoration, marker-size / color / ylim resolvers,
  clone-name formatting, region shading, and column-name constants.
- Public API - the full surface re-exported from `cnplot` with an explicit `__all__`; no
  imports from the sibling repos.
- Tests and examples - a deterministic simulator over a real 22-chromosome T2T-CHM13v2.0
  profile emits every input format; the suite renders through `Agg` and asserts structural
  invariants. `examples/plot_gallery.py` renders the six README figures.

### Notes (deliberate behavior changes vs the origin repos)

Downstream migration (see `TODO.md`) must absorb these; they are intentional, not bugs.

- Mirror rule corrected - a bin is mirrored only when clones disagree on which allele is
  amplified, not when any single clone has `b > a`. Removes false hatching on uniform `1|2`
  bins and `1|1` + `1|2` mixes. The mirror swatch is drawn only when a bin was actually hatched.
- Gaps are first-class and found by comparing regions to the chromosome, not counting segments.
  A collapsed-gap dash marks removed width; under `collapse_gaps=False` the width is visible, so
  no line is drawn (both host repos always drew it). The same rule breaks expected lines and
  masks heatmap columns across uncovered sequence.
- `contain=True` is the profile default - a reference mismatch becomes a loud unmapped-row
  warning instead of a plausible-looking wrong figure.
- Coordinates are matched by `(#CHR, START, END)`, not indexed by row position, removing the
  length check against the table being drawn.
- Chromosome order always follows the sizes file, never sorted, so `chr10` stays after `chr9`.
- Drawing functions that take an `ax` return None; only figure builders return a `Figure` /
  `JointGrid`. Profile/scatter legends are private, reached via `ax_leg=`.
