# cnplot TODO

Goal: a PyPI-installable library (`cnplot`) for copy-number plotting, factored out of the
sibling repos inventoried in [plan.md](plan.md):

- HATCHet3 - `hatchet-long-read/src/hatchet/plot/` (bulk clone CN, CNP-string driven)
- Copytyping - `Copytyping-dev/src/copytyping/plot/` (single-cell CN, cell-matrix driven)
- UGP - `universal-genotyping/workflow/scripts/plot_utils.py` (reference only, not a consumer;
  still a useful third data point for the M1 coordinate model)

## Scope

In: genome-coordinate CN visualization - colormaps, integer-CN profiles and legends,
1D/2D scatter, single-cell heatmap. Consumers are HATCHet3 and Copytyping.

Out: preprocessing QC (stays in UGP), spatial/Visium (stays in Copytyping), and any function
whose core computation needs scipy. This keeps the dependency set at numpy, pandas,
matplotlib, seaborn, adjustText - no optional scientific dependencies.

Priority: M1 -> M2 -> M3 was the critical path; it retires the hard duplication in plan.md
section 3. M0-M6 are done. M7+ are additive.

Ground rule: cnplot never modifies or deletes files in the sibling repos, and downstream
migration is not interleaved with the build. Finish the package first (M0-M8), then do the
migration in one pass - see M9. Until then the host repos keep their own copies and nothing
outside `cnplot/` is touched.

---

## Module layout

| Module | Contents |
|---|---|
| `cnplot_io_utils` | `read_chr_sizes`, `read_bed` - reference readers, pandas only |
| `cnplot_colormap` | integer-CN, allele-CN, and categorical label palettes |
| `cnplot_genome_axis` | `GenomeAxis`, `Gap`, `Segment`, `BinCoords` - pure geometry |
| `cnplot_utils` | column-name constants, axis decoration, marker/color/limit resolution, `FigureSaver` |
| `cnplot_intcnp` | `plot_cnv_profile` - draws its own legend (+ legacy ascn pair) |
| `cnplot_1d` | `plot_scatter_1d` - the per-axes 1D primitive |
| `cnplot_2d` | `plot_scatter_2d`, `get_landmarks` |
| `cnplot_heatmap` | `plot_heatmap`, `plot_column_strips`, `plot_strip_legend` |
| `cnplot_figures` | figure builders that compose the primitives: `plot_scatter_1d_multisample`, `make_row_spec` (M7 panels + M6/M9 `plot_heatmap_cnp` land here too) |

Two tiers, kept in separate modules: **primitives** take an `ax` you own and return
None (`plot_scatter_1d`, `plot_cnv_profile`, `plot_heatmap`); **figure builders** own the
layout, create a figure, call the primitives, and return it (`cnplot_figures`). This mirrors
HATCHet's own split (primitives in `plot_1d2d.py`, `plot_combined_1d` builder in
`plot_utils.py`). `plot_scatter_2d` returns a `JointGrid` but stays a primitive - it is one
plot that owns a figure only because seaborn forces it, not a composition.

Dependencies run one way: `io_utils` and `colormap` are leaves, `genome_axis` -> `io_utils`,
`utils` -> `genome_axis`, `intcnp` -> all three, `1d` / `2d` / `heatmap` sit on top of the
primitives, and `figures` sits above everything (imports `1d`, `intcnp`, `utils`).
`io_utils` and `genome_axis` are matplotlib-free, so the coordinate logic is testable
without a plotting backend. Verified acyclic by AST scan of the relative imports.

## M0. Package scaffolding - DONE

Verified 2026-07-21: `python -m build` produces `cnplot-0.1.0.dev0-py3-none-any.whl`;
installing it into an isolated target resolves `cnplot.__version__`. `ruff format --check`
and `ruff check` clean.

- [x] `pyproject.toml`: dist name `cnplot` (matches the import name; `py-cnplot` was also
      free but would split the two), `requires-python = ">=3.10"`, setuptools>=77 backend,
      `package-dir = {"" = "src"}`, version dynamic via `attr = "cnplot.__version__"`.
- [x] Deps: `numpy`, `pandas`, `matplotlib`, `seaborn`, `adjustText`. Only extra is `dev`.
      `seaborn` is needed by M2 (`sns.color_palette`) and M4/M5 (`scatterplot`, `kdeplot`,
      `move_legend`); `adjustText` by M5/M7 label placement.
- [x] `LICENSE`: MIT, declared via PEP 639 (`license = "MIT"` + `license-files`), which is
      why the build needs setuptools>=77. No sibling repo had a LICENSE to match.
- [x] `[project]` metadata: authors, keywords, urls (raphael-group/cnplot), classifiers.
- [x] `.gitignore`. Verified with `git status --untracked-files=all` that `build/`, `dist/`,
      `*.egg-info/`, and `.DS_Store` stay hidden after a real build.
- [x] README: fixed "Pythong" typo, install name, added install-from-source.
- [x] `[tool.ruff]` (line-length 88, py310, double quotes, `E,F,I,UP,B`) and `[tool.pytest]`.
- [x] `src/cnplot/__init__.py`: docstring + `__version__`. Re-exports pending (M8).
- [x] No CLI entry point - API-only library; plotting CLIs stay in the host repos.

`twine check dist/*` moved to M8 (publish), where the upload actually happens.

## M1. `cnplot_genome_axis.py` - genome coordinate model - DONE

The strongest cross-repo duplication: all three repos mapped bins to genome x-coordinates
with incompatible implementations, each returning a different tuple.

`GenomeAxis` is a **coordinate transform, not a parallel array**. It is constructed from
the reference alone and then maps any bin table by `(#CHR, START, END)`:

    genome_axis = GenomeAxis(region_bed, chrom_sizes)      # once
    coords = genome_axis.build_coordinates(seg_df, contain=True)

One axis therefore serves several samples, several bin resolutions, and every panel of a
figure. An earlier version indexed per-bin arrays by row position, which forced a length
check against the table being drawn; matching by coordinate removed that coupling
entirely.

- [x] Subsumes `get_abs_positions_ignore_gap` / `get_abs_positions_keep_gap` (HATCHet),
      `build_wl_coords` (Copytyping), and `_genome_coords` (UGP).
- [x] A plain class: `__init__(region_bed, chrom_sizes, excluded_chroms, chr_shift,
      collapse_gaps)` reads both files itself via `read_bed` / `read_chr_sizes`. There is
      no separate builder function and no derive-a-whitelist fallback - the reference is
      always supplied.
- [x] `Segment` records carry the transform; `chrs`, `chr_end`, `chr_boundaries`,
      `xtick_chrs`, `xlab_chrs`, `chr_offsets` are properties, so nothing can drift out of
      sync with `ch_coords`.
- [x] `build_coordinates(df, contain=)` returns `BinCoords`; `contain=True` takes bins
      fully inside a segment (integer-CN profiles), False clips overlapping bins (scatter,
      heatmap). `grid(df)` returns the `pcolormesh` mesh for M6.
- [x] Both input shapes: `#CHR/START/END` bins and `#CHR/POS` sites.
- [x] `excluded_chroms`, default `("chrX", "chrY", "chrM")`, filters both the sizes and the
      regions before building. Matched ignoring case and any `chr` prefix, so an
      unprefixed reference works too. Pass `[]` to keep everything.
- [x] Unmapped bins are reported: `build_coordinates` warns with a split between bins on
      chromosomes the axis does not carry and bins falling inside gaps. Only mapped bins
      are drawn.
- [x] `decorate_genome_axis`, `draw_chr_boundaries`, `draw_segment_boundaries`, and
      `shade_regions` in `cnplot_utils` replace the three inline decorators. They take
      `genome_axis`, not `axis`, to keep it distinct from matplotlib's `ax`.

Behavior changes vs the originals, all documented in the docstrings:

- **Gaps are first-class, and shrinking them is optional.** `collapse_gaps=True` (default)
  removes uncovered stretches, so the axis is the regions concatenated and every `Gap` is
  zero-width - what both host repos did. `collapse_gaps=False` lays chromosomes out at
  full length, so bins keep true spacing and every `Gap` keeps its real width. Neither host
  repo could do the latter: HATCHet's `get_abs_positions_keep_gap` preserves spacing but
  discards the whitelist, so it cannot say where the gaps are.

  The same bin gets **different coordinates** under the two, so one axis must not be shared
  across them. `Gap` keeps `raw_start`/`raw_end` either way, so a shrunk dash can still
  report the size of the region it stands for.

- **Gaps are found by comparing segments, not counting them.** Both host repos placed
  dashed lines by *counting* (`si < n_seg - 1`, plus a broken special case in Copytyping
  that could only fire on single-segment chromosomes and recorded the segment END rather
  than the leading boundary). cnplot compares consecutive regions against the chromosome,
  so a gap is reported only where material is actually missing:

  | Gap | Detected from | Lands on |
  |---|---|---|
  | leading | first region `START > 0` | chromosome start offset |
  | interior | `next.START > this.END` | strictly inside the chromosome |
  | trailing | last region `END < chrom_size` | chromosome end offset |

  A dashed line stands in for **width that was removed**, so it is drawn only for a
  collapsed gap. Under `collapse_gaps=False` the missing stretch is already visible as empty
  axis and gets no line at all - previously it got two, bounding the gap, which marked
  something the reader could already see. `shade=` fills those instead. Fixed 2026-07-22
  after the M4 renders made the doubled lines obvious.

  Consequences the count-based rule got wrong: adjacent regions leave no gap and are no
  longer marked; trailing gaps are detectable at all, which is why `chrom_sizes` is
  required; a chromosome with no regions is dropped when gaps are shrunk and kept as one
  full-length gap when they are not; a region overrunning its chromosome raises.

  Verified across the full matrix: 0/1/2/3 regions per chromosome, each with and without
  leading, interior, and trailing gaps, adjacent vs separated regions, and both
  `collapse_gaps` settings side by side.

- **Chromosome order always follows the sizes file**, never sorted, so `chr10` does not
  land between `chr1` and `chr2`. Order is the user's knob: reorder that file and the axis
  follows. Earlier the two layouts took order from *different* files, so flipping
  `collapse_gaps` could silently reorder the genome.

- `shade_regions` intersects each interval with each segment, so it works on either layout
  and a region crossing a shrunk gap shades as the two pieces actually drawn. It takes the
  region table, the same shape the axis is built from; the second `{chrom: [(s, e)]}`
  representation and its reader are gone.
- `get_transparency` takes `cols=("RD", "BAF")` instead of hardcoding those columns.
- UGP drew chromosome labels via `ax.text`; cnplot uses real ticks.
- `get_clone_ylabels` (was Copytyping's `_clone_ylabels`) lives in `cnplot_utils`, takes
  clone *names* rather than a count, and types `clone_ploidies` as `dict` - which is what
  the body always required, though Copytyping annotated it `list | None`.

## M2. `cnplot_colormap.py` - palettes - DONE

- [x] `get_cn_cmap()` - was byte-identical in both repos; now single-source in cnplot
      (the host-repo copies stay put until their owners retire them). Palette and
      state list lifted to module constants so the table is visible without reading the
      function. Verified against the original 20 entries, both orderings.
- [x] `get_ascn_cmap()` - same treatment; values verified.
- [x] `build_mixture_cn_cmap` + `_distinct_subclonal_colors` (HATCHet).
- [x] `get_baf_cmap` (Copytyping). Note it uses `ListedColormap` + `BoundaryNorm`, not
      `TwoSlopeNorm` as the earlier inventory claimed.
- [x] Label-color machinery: `build_label_cmaps` (the live one, used by the heatmap
      strips), `_clone_order_key`, `_is_normal_like`, `_is_colored_label`. Removed
      2026-07-23 as dead code: `build_label_colors` (a clone-indexed color *list* for
      Copytyping's `plot_scatter_2d_per_cell`, which stays upstream) and
      `build_categorical_cmap` (per-annotation named palette, e.g. datasets in Set1), plus
      the now-orphaned `_label_color_index`. Neither is subsumed by `build_label_cmaps`,
      which only draws from the shared clone palette; re-add them in M9 if the dataset /
      cell-type strips or the per-cell diagnostic migrate onto cnplot.
- [x] Naming convention (2026-07-23): built-in palettes `get_*`, data-derived `build_*`,
      global mutators `set_*`; suffix by return shape - `_cmap` for a colormap-like
      lookup or `(Colormap, Norm)`, `_cmaps` for several, `_colors` for a list.
- [x] `set_palette` ported from HATCHet; docstring flags that it mutates seaborn globals.
- [x] Resolved: `INVALID_LABELS` / `NA_CELLTYPE` ship as module defaults but every function
      that consults them takes `invalid_labels=` / `na_labels=` overrides. Callers keep
      control of what "missing" means without having to pass it every call.

## M3. `cnplot_intcnp.py` - integer CN profiles and legends - DONE

The 6-function hard-duplication block (plan.md section 3), now single-sourced.

Not needed here - covered elsewhere: `_cnp_segment_geometry` (subsumed by
`GenomeAxis.build_coordinates`), `_draw_chr_boundaries` / `_decorate_cnp_xaxis` (M1),
`_clone_ylabels` (now `get_clone_ylabels` in `cnplot_utils`), and the palettes (M2).

Public surface is one function; everything else is support code:

    genome_axis = GenomeAxis(region_bed, chrom_sizes)
    plot_cnv_profile(ax, seg_df, genome_axis, ax_leg=ax2, sample_id="HT941")

- [x] The interface is the seg.ucn table itself, not parallel arrays: "#CHR"/"START"/"END",
      `cn_<clone>` columns of `"a|b"` strings, matching `u_<clone>` proportions, an
      optional `PI_VIOL` flag, and an optional `SAMPLE` column. Clone names, stacking
      order, proportions, the normal clone, and the violation flags all come from the
      columns.
- [x] Table readers: `select_sample`, `get_clone_names`, `get_clone_states` (n, k, 2 int
      array), `get_clone_proportions`, `get_pi_viol`, `has_mirror`. The CNP-string helpers
      are gone - a `";"`-joined string was only ever an intermediate the loaders built on
      top of these same columns.
- [x] `sample_id=` selects one sample from a multi-sample table; None takes the first and
      logs which.
- [x] `normal="normal"` names the clone to exclude, by column name. None draws every clone,
      and `clones=` overrides the set and stacking order outright. The earlier structural
      detector is gone: with named columns the name is the identifier, so a genome-wide
      diploid tumor clone is no longer ambiguous.
- [x] `plot_cnv_profile` owns the legend. It takes an optional `ax_leg` and calls the
      legend itself, so `__all__` is one name and a caller cannot pair a profile with a
      legend that disagrees with it. Both legends are private (`_plot_cnv_legend`,
      `_plot_ascn_legend`); styling still reaches them through `legend_kwargs=`.
- [x] The mirrored swatch is drawn only when a mirrored bin was actually hatched. The flag
      passed to the legend is the same `mirrored` array that drives the hatching, not a
      recomputation, so `show_mirror=False` suppresses both and the legend can never
      advertise a state the figure does not contain. Same wiring on the ascn path, whose
      indicator is chevrons rather than a hatched box.
- [x] Functions handed an `ax` draw on it and return None; only a function that creates the
      figure returns anything (`plot_scatter_1d_multisample` -> Figure, `plot_scatter_2d` ->
      JointGrid). Returning the caller's own axes "for chaining" was never used by anything
      in the package or either host repo, and left `cnplot_utils`'s decorators as the only
      helpers following the rule. All 18 drawing functions audited.
- [x] `PI_VIOL` overlay, picked up automatically when the column is present (matched
      case-insensitively), suppressible with `show_pi_viol=False`.
- [x] Mirror rule corrected. A bin is mirrored when its clones **disagree** on which allele
      is amplified - some clone with a > b and another with a < b - not when any single
      clone has b > a. Both host repos use the per-clone test, which fires on a uniformly
      `1|2` bin (a pure allele-labelling artefact, since the palette gives (a, b) and
      (b, a) one color) and on `1|1` + `1|2` (nothing disagrees). One `_mirrored_bins`
      helper now drives the hatch, the chevrons, and the legend swatch.
- [x] `plot_ascn_profile` / `plot_ascn_legend` moved to a Legacy section: unexported,
      docstring-flagged, kept for the HATCHet figures that use them. They are not
      superseded - the joint palette never shows per-allele copy number - so the note says
      so rather than pointing at a replacement.

Verified by reading all 7 benchmark `seg.ucn` files off disk and plotting with no manual
assembly, across 3- and 4-clone solutions: clone names and order, state array shape,
proportions, y-tick count/order/labels, every rectangle's facecolor against
`get_cn_cmap`, `normal=None`, a `clones=` subset, custom and missing normal names,
`show_prop=False`, `PI_VIOL` pickup in both cases, multi-sample selection by default and
by id, and the empty-table / missing-column / unknown-sample errors.

Reconciled drift, all now explicit parameters rather than per-repo behavior:

- HATCHet's `show_prop` survives as a flag; values come from the `u_<clone>` columns.
  Copytyping had no proportion line.
- Copytyping's clone separators become `clone_separators` (default on).
- Bin membership (containment vs overlap) lives in the axis, so neither repo's choice is
  hardcoded.
- Copytyping bolded CNV row labels, HATCHet did not; cnplot bolds them.
- The mirror-rule fix above changes output for both repos - flag it in M9.

`contain` is now an explicit parameter rather than a hardcoded True, which the M1 notes
always said it should be. **The default stays True.**

Retracted 2026-07-22: I first flipped the default to False after measuring that
`contain=True` drew only 14-38 of 48-72 segments across the benchmark datasets, and recorded
that as a bug. It was not. I had built the axis from `hg38.regions.bed` while the benchmark
data is **T2T-CHM13v2.0** - the seg.ucn rows match the T2T BED line for line
(`chr1 0 116796047`, `chr1 147241659 248387328`). Against the right reference the two
settings are identical, all 48-72 segments mapping either way, and bin-level unmapped counts
go from 213-383 to **0** in all 7 datasets. The lesson is worth keeping: `contain=True`
turns a reference mismatch into a loud unmapped-row warning, where `contain=False` would
quietly clip and produce a plausible-looking wrong figure. That is a reason to keep it.

Left for M3: nothing. The one open thread is the differential test, which still targets the
pre-class `build_genome_axis` signature and needs rewriting against `GenomeAxis` - tracked
under M8.

## M4/M5 shared design - DONE

Settled 2026-07-22 after reading all five source files. Both milestones build one
single-axes core; the page builders sit on top. The three repos already converged on the
same three tiers, which is why one core per dimension is enough:

| Tier | HATCHet | Copytyping | UGP |
|---|---|---|---|
| 3 page builder | `plot_combined_1d` + `_row_def`; `plot_rdr_baf` | `plot_rdr_baf_1d_pseudobulk`, `plot_rdr_baf_2d_pseudobulk`, `plot_scatter_2d_per_cell` | `plot_rdr_baf`, `plot_1d_multi_sample`, `plot_1d_sample`, `plot_rd_1d_scatter` |
| 2 single-axes | `plot_1d`, `plot_2d` | `plot_scatter_1d_pseudobulk` (2D not factored out) | inline `ax.scatter` |
| 1 coordinates | `get_abs_positions_ignore_gap` / `_keep_gap` | `build_wl_coords` / `_build_ch_boundary` | `_genome_coords` |

Tier 1 is already retired by `GenomeAxis`. Copytyping's `_merge_exp_lines` and
`_annotate_cn_landmarks` are the factored-out forms of HATCHet's inline blocks - the latter's
docstring says so outright ("Matches the hatchet `plot_2d` style").

**Callers supply plotted values; cnplot never computes them.** Observations arrive in a
DataFrame with `#CHR/START/END`, a group column, and one column per plotted axis named by
`xcol` / `ycol`. FCN, mhBAF, log2RDR, RDR are all just a column name: cnplot applies no
gamma scaling, no log2, no allele split. This is forced by the deepest divergence between
the repos - HATCHet's expectation is a proportion-weighted mixture over all clones
(`get_expected_baf_fcn`, plot_utils.py:324), Copytyping's is a single pure clone
(`clone_C / sum(base_props * clone_C)`, plot_scatter_1d.py:346). Bulk mixes clones,
single-cell does not. That is a data-model fact and stays upstream.

**Sample and cell-group are the same axis.** One row/page per group, whatever a group
means: a bulk sample in HATCHet, a cell group in Copytyping. The layouts are already
identical - `plot_combined_1d`'s own docstring says "Layout mirrors copytyping
`plot_rdr_baf_1d_pseudobulk`". One `groups` parameter and one group column (default
`SAMPLE`, reusing the M3 constant), no per-repo branch.

Expected values therefore key on the **group**, not the clone. In Copytyping the group is a
clone so the two coincide; in HATCHet the group is a sample and the clone mixture is already
collapsed into it. A missing column means no expected overlay, which is exactly how
Copytyping's non-clone groups (`"NA"`, unmatched labels) drop out today.

- [x] `MARKER_SIZE_SMALL = 2.0`, `MARKER_SIZE_LARGE = 20.0`, `MAX_NDOTS = 5000` in
      `cnplot_utils`, all overridable per call. Below `MAX_NDOTS` dots use the large marker,
      at or above it the small one. Measured evidence for the threshold, from
      `copytyping-benchmark/datasets/bulk_cnprofiles` (7 datasets):

      | Level | Dots | Repo default |
      |---|---|---|
      | HATCHet bulk `bbc.ucn` bins/sample | 15,701 - 21,760 | `s=2` |
      | HATCHet `bbc_phases` bins/sample | 7,977 - 11,452 | - |
      | HATCHet `seg.ucn` segments | 48 - 72 | - |
      | Copytyping `cnp_bin` | ~500 - 3,000 (estimated) | `s=20` |

      Copytyping's `cnp_bin` grid merges bbc rows until pseudobulk SNP count reaches
      `min_snp_count=300` or the span exceeds `max_bin_length=5_000_000`, never crossing a
      segment (count_data.py:265-278). The 5 Mb cap alone floors the genome at ~500-600
      bins; the SNP cap binds far earlier in bulk (<=6.7k-9.1k if single-cell allele
      coverage matched bulk `#SNPS`, which it does not) so real G sits near the floor.
      **Not directly measured** - no Copytyping run output is on disk. Re-check against one
      real run before locking the constant. 5000 also matches the `n_ref` already used by
      the ported `adaptive_dot_size`, so the two agree by construction.
- [x] Resolve the size **once per figure** from the total dot count, not per axes. A binary
      switch applied per row would give two rows of the same figure a 10x marker difference
      at 4,999 vs 5,001 dots. `adaptive_dot_size` (continuous, already in `cnplot_utils`)
      stays available for callers who prefer no discontinuity at all.
- [x] Default `get_transparency` (already ported and generalized to arbitrary `cols`) on for
      both repos, not just HATCHet, so tail fading is a package style rather than a
      per-repo accident. Copytyping has no equivalent today.
- [x] One adaptive-ylim helper, applied for both. Both repos already do this, in different
      places and by different rules: Copytyping inside the plot function (log2 tightens to
      (-2, 2) when observed and expected both fit, else (-5, 5); linear takes the data max
      clamped to 6.0 - plot_scatter_1d.py:355-371), HATCHet in the caller (`lim_fcn` from
      the data max clamped by `maxlim_fcn`, `lim_baf` as (0, 1) or (0, 0.55) -
      plot_cn.py:160-168). Fold into `cnplot_utils` and drive both off the observed and
      expected values actually passed in.

## M4. `cnplot_1d.py` + `cnplot_figures.py` - 1D genome scatter - DONE

Split 2026-07-23: `plot_scatter_1d` (per-axes primitive) stays in `cnplot_1d`;
`plot_scatter_1d_multisample` + `make_row_spec` (the figure builder) moved to the new
`cnplot_figures` module. See the Module layout note on the primitive / figure-builder tiers.


- [x] `plot_scatter_1d(ax, obs_df, genome_axis, ycol, ...)`. Once tier 1 is removed,
      `plot_1d` and `plot_scatter_1d_pseudobulk` are the same renderer step for step:
      scatter -> rasterize -> chromosome vlines -> expected `LineCollection` -> ylim/ylabel/
      title -> xticks at chromosome midpoints rotated 60. Only four real differences, all
      parameters: `hue`+`palette` vs a per-bin `colors` array; per-point `alphas`; marker
      style; expected lines derived from a column vs passed in.
- [x] `expected_df` (default None): `#CHR/START/END` plus one `exp_{ycol}_{group}` column
      per group. May be at **segment** resolution while `obs_df` is at bin resolution, since
      both map through the same `GenomeAxis` - expected values are piecewise-constant per
      segment, so that is the honest representation.
- [x] Port `_merge_exp_lines` (Copytyping): join equal-valued neighbours into one stroke.
      Extended past both host repos' version - a run also breaks where the next row does not
      start where the last ended (1 bp tolerance). Without it a kept-gap layout draws the
      expected line straight through every centromere, asserting a copy number over sequence
      no bin covers. Measured on HT973N1: 17/17 gaps bridged before, 0 after; the same
      break also caught 2 within-region holes the shrunk layout hid.
- [x] Gap dashes are drawn only when the gap was collapsed (`draw_segment_boundaries`). A
      dash stands in for width that was removed; under `collapse_gaps=False` the width is
      already visible as empty axis, so a line there marks nothing - `shade=` fills it
      instead. Both host repos always drew the line; cnplot ties it to `Gap.is_collapsed`.
      This is the same principle as the expected-line break: never draw across sequence that
      is not on the axis.
- [x] Hue is a column name plus a palette dict. HATCHet colors by joint CNP cluster string
      (`build_mixture_cn_cmap`, plot_cn_utils.py:383), Copytyping by that group's (A, B) state
      (`get_cn_cmap`). Different meaning, identical rendering call - so the core takes a
      categorical label per row and a palette, and neither semantics is baked in.
- [x] Break the `g0_colors` coupling. HATCHet runs `plot_2d` first and feeds its returned
      face colors into `plot_1d` as `colors=` (plot_cn.py:230 -> plot_utils.py:156), making
      1D depend on 2D having been drawn. Replace with a shared
      `resolve_colors(hue, palette, alphas) -> RGBA` used by both cores, so ordering is free.
- [x] Two-row RDR+BAF page builder - implemented three times (`plot_combined_1d` + `_row_def`,
      `plot_rdr_baf_1d_pseudobulk`, UGP `plot_rdr_baf`). One row-spec driven function;
      HATCHet `_row_def` is the closest model.
- [x] Multi-sample variant (UGP `plot_1d_multi_sample`): one row per group, shared axis,
      as `plot_scatter_1d_multisample`. Its BED shading is covered by `shade_regions` and
      `adaptive_dot_size` is ported.
- [x] Clone proportions live beside each group's rows, not on the shared profile
      (`show_props`, default on; the profile gets `show_prop=False` unless overridden).
      Several samples routinely share one clonal CN solution while mixing those clones in
      different amounts, so printing proportions against the shared profile attributes one
      sample's mixture to all of them. This is HATCHet's arrangement - `plot_combined_1d`
      passes `show_prop=False` (plot_utils.py:230) and builds a per-sample legend on each
      sample's bottom row (:198-218) - and cnplot follows it, including
      `display_min_clone_prop` and always listing normal.
- [x] `format_clone_name` in `cnplot_utils` renders clone names for display, shared by the
      row labels and both proportion legends so a profile cannot say "Clone 1" while its
      legend says "clone1".
- [x] Kept/filtered two-color via an optional bool `keep_col` in `obs_df`: False rows are
      still drawn, in `filtered_color`, with kept/filtered counts in the legend, so a
      filter can be judged against what it removed. Two deliberate differences from UGP's
      `mask`:

      - Kept rows keep their `hue`; UGP flattens both classes to blue/red and throws the
        CN-state coloring away for no reason.
      - The legend appears on the top row only, not once per row.

      Worth knowing what the flag can mean. UGP's `snp_mask` is a conjunction of region
      inclusion, blacklist exclusion, RNA coverage, and `exon_only`
      (aggregation_utils.py:699-716, phase_and_concat_nonbulk.py:108-136). Only the first
      overlaps `GenomeAxis`, which already drops out-of-region bins and warns - red would
      draw nothing for those. The useful reading is "inside the regions yet excluded for a
      data-quality reason", which is the other three components. Also note the UGP path is
      dead code: `snp_mask` defaults to None (plot_utils.py:844), the sole forward is
      :908, and no caller in the repo ever sets it.

## M5. `cnplot_2d.py` - 2D RDR/BAF scatter - DONE

- [x] `plot_scatter_2d(obs_df, xcol, ycol, ...)`. Same JointGrid skeleton in both repos:
      `sns.JointGrid(hue, palette)` -> `plot_joint(scatterplot)` -> `plot_marginals(kdeplot,
      common_norm=False, fill=False)` -> refline -> hollow black circle and `adjust_text`
      label per landmark. HATCHet-only extras to parameterize: the `clone_props` legend,
      `alphas`, and `filtered_ids`/`balanced_ids` (skip a landmark; square vs circle marker).
- [x] **2D expected values are indexed differently from 1D**, and the right table is the
      seg.ucn layout M3 already reads. 1D wants one value per (bin, group), a line along the
      genome. 2D wants one (x, y) landmark per distinct CN identity within a group: HATCHet
      keys landmarks by CNP cluster (`exp_labels` = `unique_labels`, one per distinct joint
      state across all clones), Copytyping by distinct (a, b) of the single group's clone.
      Both are "dedup on the `cn_<clone>` columns you care about", so `expected_df` is a
      seg.ucn-layout table - `SAMPLE`, `cn_<clone>`, `u_<clone>` - plus precomputed
      `exp_{xcol}` / `exp_{ycol}` columns. Consequences:

      - The M3 `clones=` subset parameter is the same knob that switches HATCHet's
        joint-state dedup to Copytyping's single-clone dedup. No per-repo branch.
      - `select_sample`, `get_clone_names`, `get_clone_states`, `get_clone_proportions` in
        `cnplot_intcnp` are reused as-is; M5 adds no parsing path.
      - The `(a, b)` landmark label has to come from the `cn_` columns - precomputed
        `exp_` coordinates alone cannot produce it. This is why the CN columns must be
        present and a bare (x, y, label) table would not do.
      - `u_<clone>` gives HATCHet's `clone_props` legend directly.
      - HATCHet's `filtered_ids` / `balanced_ids` become optional `is_filtered` /
        `is_balanced` boolean columns (both already exist in HATCHet's segs table -
        compute_cn.py:65, scaling.py:51), driving landmark skip and circle-vs-square.
        Absent for Copytyping.
      - Genome coordinates stay **optional**. A landmark is determined by CN identity, not
        position; `#CHR/START/END` (seg.ucn) and `#ID` (HATCHet's cluster-level `.seg`,
        keyed `(#ID, SAMPLE)`) are both accepted for provenance, neither is required.
- [x] Unify `_annotate_cn_landmarks` (Copytyping) with HATCHet's inline CN-centroid labels.
- [x] A clonal state - every visible clone agreeing on one (a, b) - is labelled **bold**,
      subclonal states plain (HATCHet plot_1d2d.py:374). The normal clone defaults out of
      the label, as HATCHet does by starting at `range(1, len(cid))`; leaving it in makes
      its constant `1|1` pollute every label so nothing is ever clonal.
- [x] `display_min_clone_prop` drops low-proportion clones from the label and the clonal
      test (`_tumor_visible`), falling back to all clones if it would hide every one.
      Proportions are per group, so the same joint state can read clonal in one sample and
      subclonal in another - clonality is a property of the sample, not the state alone.
      The proportion legend honours the same threshold; normal always stays listed, being
      the baseline the tumor fractions are read against.
- [x] `seen_labels`: distinct states can collapse to the same label once hidden clones drop
      out. Both keep their marker, only the first carries text - matching HATCHet, where
      the `continue` skips the text while the circles come from a separate `is_vis` mask.
      Needs 3+ tumor clones to arise, so the n3 benchmark never exercised it; covered by a
      synthetic fixture instead.
- [x] Per-cell diagnostics (`plot_scatter_2d_per_cell`) deliberately **not** ported. It
      stays in Copytyping and is refactored there to call cnplot instead - tracked as the
      last M9 item. It is a single-consumer diagnostic with no cross-repo duplication to
      retire, so moving it would grow the library without serving its purpose.
- [x] Excluded: HATCHet `plot_clusters` + `_is_multimodal`. Their scipy use is load-bearing,
      not incidental - `beta.pdf`/`norm.logpdf` draw the theoretical densities,
      `beta.ppf`/`norm.ppf` the QQ lines, and `gaussian_kde` + `find_peaks` are the
      multimodality test itself. Reopen only if HATCHet needs them served from cnplot.
- [x] Excluded: UGP `plot_rd_2d_kde`. Despite the name it is GC/mappability/repliseq
      bias-correction QC, not an RDR/BAF scatter. UGP has no 2D RDR/BAF plot at all.

Verified 2026-07-22 end to end against real HATCHet data: hg38 arm-level `regions.bed` +
`chrom.sizes`, and the HT973N1 `bbc.ucn` (18,253 bins) / `seg.ucn` (48 segments) pair, with
the bulk mixture expectation computed upstream in the test the way a caller would. Checked:
point count drawn equals mapped-and-finite bins (18,017, with 236 correctly dropped into
centromere gaps), every x inside the axis span, auto marker size resolving to SMALL for a
16k-bin sample and staying identical across panel rows, 31 merged expected-line segments
from 48 input segments, an unknown group silently drawing no overlay, chromosome labels on
the profile row only, 4 distinct joint landmark states vs 3 for `clones=['clone1']` (the two
host conventions off one table), `is_filtered` / `is_balanced` honoured, and the
proportion legend picking up `u_<clone>`. Plus unit checks on the four new helpers.

Left for M4/M5: nothing. Both cores ship, the page builder ships, and the dedup that
justified these milestones is done. `plot_scatter_2d_per_cell` stays in Copytyping and
becomes an M9 consumer of the library rather than a port.

## M4/M5 resolved question

- [x] **`region_bed` stays required in `GenomeAxis`.** HATCHet's cluster-diagnostic path
      calls `plot_1d(..., regions=None, ignore_gap=False)` (plot_1d2d.py:822) and has no
      region BED, so the option was to weaken the cnplot API or fix the caller. Resolved
      2026-07-22 in favour of the caller: that path should read a region BED like every
      other one. One coordinate contract, no conditional-required argument, and the
      diagnostic plots gain the shrink-gap layout they currently cannot use. The plumbing
      is a HATCHet change, so it belongs to M9 - see the `cluster_bins` item there.

## M6. `cnplot_heatmap.py` - single-cell CN heatmap - DONE

Copytyping-only, no dedup pressure. The line between what ported and what stayed follows M5:
cnplot draws a matrix already reduced to rows; the reducers - which need
`copytyping.inference.model_utils` and the `anns` table - stay upstream. Copytyping's
imports made the reuse map explicit: everything `plot_heatmap` pulled from
`plot_common` / `plot_copynumber` (`build_wl_coords`, `get_baf_cmap`, `BAF_COLORS`,
`build_label_cmaps`, `PURITY_CMAP`, `plot_cnv_profile`) already exists in cnplot.

- [x] `plot_heatmap(ax, matrix, coords_df, genome_axis, row_labels=, ...)`. Takes a
      precomputed `(n_rows, n_bins)` matrix; columns are placed by `GenomeAxis.grid`, so the
      whole coordinate block, the chromosome vlines, and the centromere dashes are the M1
      helpers rather than the inline versions. Returns `(x_edges, y_edges, masked)` so a
      caller can align side strips to the same rows.
- [x] Filler columns are masked, so no value bleeds across a gap - the mesh analogue of the
      expected-line and dash fixes. Verified: masked columns equal exactly the `-1` filler
      columns from `grid`.
- [x] Categorical side strips and their legends are `plot_column_strips` and
      `plot_strip_legend` (renamed from Copytyping's `plot_label_strips` /
      `draw_label_legends`). General chrome, no data-model dependency; `_display_name`'s
      hardcoded label map became a `display_names=` argument. `plot_heatmap` calls both
      itself when given `strip_label_map=`, mirroring the profile's `ax_leg=`, so one call
      draws mesh + strips + legends and the user need not orchestrate the pieces.
- [x] `plot_heatmap` also draws the value **colorbar** (`show_colorbar=`, `cbar_label=`,
      `cbar_ticks=`) right of the axes, before the legends. Copytyping had this as a separate
      `add_colorbar` helper in the page builder.
- [x] Reuses M1 `grid()` for the `pcolormesh` mesh; the collapsed-gap dashing reuses
      `draw_segment_boundaries`, so it inherits the M4 fix (dash only where width was
      removed).
- [x] Dense-only; Copytyping's heatmap has no sparse handling to port.
- [x] Verified against the M6 tests (real profile / T2T sim): mesh geometry, filler==masked,
      one boxed y-tick per label run, chr labels on top, the column-mismatch guard, aligned
      side strips whose per-value fractions sum to 1, both legends, and the colorbar axes.

Coloring note (not an M6 bug): a scatter colored by a **single** clone's state collapses a
subclonal segment onto a clonal one at the same total CN. The joint-CNP coloring is
`build_mixture_cn_cmap` (M2), which gives subclonal states their own colors while clonal states
keep the integer-CN color - that is what the gallery uses (`hue="cnp"`), and the palette has
no color collisions (each distinct `(a, b)` is unique).

Stays in Copytyping (upstream, needs `model_utils` / `anns`), tracked for the M9 refactor:
`_row_layout`, `_aggregate_columns`, `_mode`, `prepare_rdr` / `prepare_baf` / `prepare_pi_gk`,
and the `plot_cnv_heatmap` page builder that pools counts and calls all of the above.

## M7. Additional modules (not yet stubbed)

- [ ] `cnplot_tree.py`: HATCHet `render_cnt_tree` (clone tree). Check its drawing backend
      first - may pull a non-matplotlib dependency, which would put it out of scope.
- [ ] Multi-solution CNP grids into **`cnplot_figures.py`** (not a separate panel module -
      it is the one home for figure builders): HATCHet `plot_cnp_panel.run`, `plot_pool_cnp`,
      `plot_summary_pdf` / `plot_bars`, `plot_scaling_2d`. Strip the solution-loading glue
      (`override_solution`, `load_gammas`, `get_expected_baf_fcn`) and keep only the drawing.
- [ ] `plot_heatmap_cnp` into `cnplot_figures.py` (M6/M9): the heatmap page - heatmap + CN
      profile + legend + colorbar + strips - currently composed by hand in the gallery.
      Should expose its own `profile_hspace`-style spacing knob so callers do not set
      `subplots_adjust` themselves.

## M8. Public API and quality

- [ ] `__init__.py`: explicit re-exports + `__all__`. No lazy-import machinery needed.
- [ ] Remove all `from hatchet.* import` / `from copytyping.* import` / `import *` usages.
- [ ] Type hints on all public functions; Google-style docstrings with Args/Returns.
- [x] Tests: `tests/` built ahead of M7 (2026-07-23). One `simulate.py` generates a
      deterministic dataset in **every** input format - bins, seg.ucn, per-bin observations,
      `exp_<col>_<group>` overlays, and per-cell heatmap matrices - from a **real 22-chromosome
      CNP profile** (69 segments, states `1|0`..`2|2`) over the **real T2T-CHM13v2.0**
      arm reference. The three files are vendored under `tests/data/` (seg.ucn ~7 KB, BED, and
      sizes), so tests are self-contained and deterministic while exercising real gap structure
      (17 interior + 5 leading gaps) and real CN diversity. Fine bins are cut inside each
      segment (~1,270 bins); observed values are the bulk clone mixture plus noise, heatmap
      rows are per-cell single-clone truths plus noise. The profile has no mirrored segment, so the
      mirror-swatch and mirror-rule tests build a small synthetic `2|1` vs `1|2` profile inline.
      46 tests over `pytest` + `Agg`: geometry and containment (no backend), palettes and the
      resolvers, the profile with its conditional mirror swatch and PI_VIOL overlay, 1D scatter
      / multi-sample / `keep_col` / the gap-dash and expected-line-break rules, 2D landmarks
      with the joint-vs-single-clone convention switch, and the heatmap (mesh + masked filler +
      the one-call side strips + legends). Plotting tests render to `tmp_path` and assert the
      figure is non-empty plus structural invariants. `pythonpath = ["src", "tests"]` and
      `known-first-party` wire it up. Pending: `pytest-cov` is in the `dev` extra, so run
      coverage there; the original M1/M2 differential could still be folded in as a golden
      check.
- [ ] Docs: README API section + an `examples/` notebook per module.
- [ ] `CHANGELOG.md` (empty): start at 0.1.0.
- [ ] Publish: `twine check dist/*` (twine is in the `dev` extra, not `base`), then
      TestPyPI first, then PyPI trusted publishing via GitHub Actions on tag.

## M9. Downstream migration - LAST, only after M0-M8 are finished

Deliberately deferred to the end. Nothing here starts while the package is still being
built, and none of it edits or deletes files in the sibling repos on cnplot's behalf.

- [ ] Precondition: M0-M8 complete, tests in place, and a version published.
- [ ] Per-function equivalence report: for every function cnplot replaces, confirm the
      cnplot version matches the host-repo original on real data, and list the deliberate
      divergences - see the M1 list. Flag one to Copytyping specifically: cnplot fixes the
      `seg_coords` leading-gap bug, so adopting it removes a stray dashed line from any
      figure with a single-segment chromosome whose whitelist starts past 0.
- [ ] Signature changes downstream must absorb, none of them mechanical:
      - `plot_cnv_legend` / `plot_ascn_legend` are private. Callers pass `ax_leg=` to the
        profile instead of calling the legend themselves. HATCHet does this in two places
        (plot_utils.py:244-252) and Copytyping in one (plot_scatter_1d.py:437-440).
      - Profile and scatter functions return None, not the axes.
      - `contain` is now a parameter, default unchanged at True.
      - The mirrored swatch disappears from any figure whose data has no mirrored bin.
- [ ] Propose re-export shims for HATCHet3 and Copytyping (`from cnplot... import ...` in
      place of the local copy) so their call sites keep working unchanged.
- [ ] **Propose `--region_bed` for HATCHet `cluster_bins`.** Its two `plot_rdr_baf` calls
      (cluster_bins.py:163 and :456) are the only ones reaching `plot_1d` with
      `regions=None`, because `cluster_bins` never receives a region BED - `add_arguments_
      cluster_bins` has no such argument, only `genome_size`. `--region_bed` already exists
      on `compute_cn`, `plot_cn`, and `plot_panel`, `required=True` on all three
      (hatchet_parser.py:622, 705, 827), so any real pipeline run already supplies the file;
      `cluster_bins` simply is not wired to it. Adding the argument and passing it through
      lets those diagnostics use `GenomeAxis` unchanged. Note it is a CLI change for anyone
      invoking `hatchet cluster_bins` standalone - recommend `required=True` to match the
      three neighbours, but that is the repo owner's call, not cnplot's.
- [ ] Hand the shims to each repo's owner. Removing the superseded local copies is their
      call, not cnplot's.
- [ ] Last: refactor Copytyping's `plot_scatter_2d_per_cell` to build on cnplot -
      `get_landmarks` / `_annotate_landmarks` for the per-clone CN landmarks, the shared
      color and marker-size helpers for the grid cells - rather than reimplementing them.
      The cross-tab layout and the BAF-histogram fallback stay Copytyping's.
- [ ] Only after a repo has adopted cnplot: revisit anything descoped for lack of a second
      consumer, e.g. HATCHet's `plot_clusters` (M5) and the clone tree (M7).

---

## Open decisions

- [ ] **Column-name contract.** Everything assumes `#CHR/START/END` and the `CNP` string
      format. Bake this in as the documented input schema, or add an adapter layer? Baking
      it in is simpler and matches both consumers today.
- [ ] **Axes-in vs figure-out.** HATCHet/Copytyping primitives take `ax`; UGP functions
      create figures and write files. Standardize on `ax`-taking primitives plus thin page
      builders.
- [x] **Migration sequencing.** Resolved: build the whole package first, migrate once at
      the end. See M9.
