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
section 3. M0-M3 are done. M4+ are additive.

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
| `cnplot_utils` | column-name constants, axis decoration, marker styling, `FigureSaver` |
| `cnplot_intcnp` | `plot_cnv_profile`, `plot_cnv_legend` (+ legacy ascn pair) |

Dependencies run one way: `io_utils` and `colormap` are leaves, `genome_axis` -> `io_utils`,
`utils` -> `genome_axis`, `intcnp` -> all three. `io_utils` and `genome_axis` are
matplotlib-free, so the coordinate logic is testable without a plotting backend.

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
- [ ] `twine check dist/*` before first upload (twine is in the `dev` extra, not `base`).

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

- [x] `get_cn_colors()` - was byte-identical in both repos; now single-source in cnplot
      (the host-repo copies stay put until their owners retire them). Palette and
      state list lifted to module constants so the table is visible without reading the
      function. Verified against the original 20 entries, both orderings.
- [x] `get_ascn_colors()` - same treatment; values verified.
- [x] `build_cnp_palette` + `_distinct_subclonal_colors` (HATCHet).
- [x] `make_baf_cmap` (Copytyping). Note it uses `ListedColormap` + `BoundaryNorm`, not
      `TwoSlopeNorm` as the earlier inventory claimed.
- [x] Label-color machinery: `build_label_colors`, `build_label_color_maps`,
      `build_categorical_color_map`, `_clone_order_key`, `_is_normal_like`,
      `_is_colored_label`, `_label_color_index`.
- [x] `set_palette` ported from HATCHet; docstring flags that it mutates seaborn globals.
- [x] Resolved: `INVALID_LABELS` / `NA_CELLTYPE` ship as module defaults but every function
      that consults them takes `invalid_labels=` / `na_labels=` overrides. Callers keep
      control of what "missing" means without having to pass it every call.

## M3. `cnplot_intcnp.py` - integer CN profiles and legends - DONE

The 6-function hard-duplication block (plan.md section 3), now single-sourced.

Not needed here - covered elsewhere: `_cnp_segment_geometry` (subsumed by
`GenomeAxis.build_coordinates`), `_draw_chr_boundaries` / `_decorate_cnp_xaxis` (M1),
`_clone_ylabels` (now `get_clone_ylabels` in `cnplot_utils`), and the palettes (M2).

Public surface is two functions; everything else is support code:

    axis = GenomeAxis(region_bed, chrom_sizes)
    plot_cnv_profile(ax, seg_df, axis, sample_id="HT941")
    plot_cnv_legend(ax2, has_mirror=...)

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
- [x] `plot_cnv_profile`, `plot_cnv_legend`.
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
`get_cn_colors`, `normal=None`, a `clones=` subset, custom and missing normal names,
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

Left for M3: nothing. The one open thread is the differential test, which still targets the
pre-class `build_genome_axis` signature and needs rewriting against `GenomeAxis` - tracked
under M8.

## M4. `cnplot_1d.py` - 1D genome scatter

- [ ] Core primitive from HATCHet `plot_1d` (richest: exp lines, hue/palette, alphas, gap
      handling, legend), reshaped to take precomputed arrays like Copytyping
      `plot_scatter_1d_pseudobulk` so it is data-model agnostic.
- [ ] Port `_merge_exp_lines` (Copytyping) - HATCHet and UGP do this inline.
- [ ] Two-row genome RDR+BAF page builder - implemented three times (`plot_combined_1d` +
      `_row_def`, `plot_rdr_baf_1d_pseudobulk`, UGP `plot_rdr_baf`). One row-spec driven
      function; HATCHet `_row_def` is the closest model.
- [ ] Multi-sample variant (UGP `plot_1d_multi_sample`): one row per sample, shared axis.
- [ ] Decide array-first vs DataFrame-first entry points; likely both.

## M5. `cnplot_2d.py` - 2D RDR/BAF scatter

- [ ] Core from HATCHet `plot_2d` and Copytyping `plot_rdr_baf_2d_pseudobulk`.
- [ ] Unify `_annotate_cn_landmarks` (Copytyping) with HATCHet's inline CN-centroid labels.
- [ ] Per-cell diagnostics: Copytyping `plot_scatter_2d_per_cell`.
- [ ] Excluded: HATCHet `plot_clusters` + `_is_multimodal`. Their scipy use is load-bearing,
      not incidental - `beta.pdf`/`norm.logpdf` draw the theoretical densities,
      `beta.ppf`/`norm.ppf` the QQ lines, and `gaussian_kde` + `find_peaks` are the
      multimodality test itself. Reopen only if HATCHet needs them served from cnplot.

## M6. `cnplot_heatmap.py` - single-cell CN heatmap

Copytyping-only, no dedup pressure.

- [ ] Port `plot_heatmap`, `plot_cnv_heatmap`, `plot_label_strips`, `draw_label_legends`.
- [ ] Port `_row_layout`, `_aggregate_columns`, `_mode`, `_display_name`.
- [ ] `prepare_rdr` / `prepare_baf` / `prepare_pi_gk` depend on
      `copytyping.inference.model_utils` (`cell_rdr_matrix`, `cell_baf_matrix`) - accept
      precomputed matrices rather than vendoring the reducers.
- [ ] Reuse M1 `x_edges` / `col_bin_ids` for the `pcolormesh` grid.
- [ ] Dense-only; Copytyping's heatmap has no sparse handling to port.

## M7. Additional modules (not yet stubbed)

- [ ] `cnplot_tree.py`: HATCHet `render_cnt_tree` (clone tree). Check its drawing backend
      first - may pull a non-matplotlib dependency, which would put it out of scope.
- [ ] `cnplot_panel.py`: multi-solution CNP grids - HATCHet `plot_cnp_panel.run`,
      `plot_pool_cnp`, `plot_summary_pdf` / `plot_bars`, `plot_scaling_2d`. Strip the
      solution-loading glue (`override_solution`, `load_gammas`, `get_expected_baf_fcn`)
      and keep only the drawing.
- [ ] `plot_loss` (training-loss curve) - decide if it belongs in a CN plotting library.

## M8. Public API and quality

- [ ] `__init__.py`: explicit re-exports + `__all__`. No lazy-import machinery needed.
- [ ] Remove all `from hatchet.* import` / `from copytyping.* import` / `import *` usages.
- [ ] Type hints on all public functions; Google-style docstrings with Args/Returns.
- [ ] Tests: golden-image or numeric-invariant per module (matplotlib `Agg`), with synthetic
      fixtures for region BEDs, sizes files, and seg.ucn tables. The differential script
      that checked M1/M2 against the inlined originals still targets the pre-class
      `build_genome_axis` signature and needs rewriting against `GenomeAxis`; fold it into
      `tests/` so the equivalence keeps being checked rather than being a one-off result.
- [ ] Docs: README API section + an `examples/` notebook per module.
- [ ] `CHANGELOG.md` (empty): start at 0.1.0.
- [ ] Publish: TestPyPI first, then PyPI trusted publishing via GitHub Actions on tag.

## M9. Downstream migration - LAST, only after M0-M8 are finished

Deliberately deferred to the end. Nothing here starts while the package is still being
built, and none of it edits or deletes files in the sibling repos on cnplot's behalf.

- [ ] Precondition: M0-M8 complete, tests in place, and a version published.
- [ ] Per-function equivalence report: for every function cnplot replaces, confirm the
      cnplot version matches the host-repo original on real data, and list the deliberate
      divergences - see the M1 list. Flag one to Copytyping specifically: cnplot fixes the
      `seg_coords` leading-gap bug, so adopting it removes a stray dashed line from any
      figure with a single-segment chromosome whose whitelist starts past 0.
- [ ] Propose re-export shims for HATCHet3 and Copytyping (`from cnplot... import ...` in
      place of the local copy) so their call sites keep working unchanged.
- [ ] Hand the shims to each repo's owner. Removing the superseded local copies is their
      call, not cnplot's.
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
