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

Priority: M1 -> M2 -> M3 is the critical path; it retires the hard duplication in plan.md
section 3. M1 and M2 are done, so M3 is next and is where the duplication actually gets
deleted from the host repos. M4+ are additive.

---

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

## M1. `cnplot_utils.py` - genome coordinate model - DONE

The strongest cross-repo duplication: all three repos mapped bins to genome x-coordinates
with incompatible implementations. Verified by a differential test against the four
originals inlined verbatim; all fields match exactly (see "Verification" below).

- [x] `build_genome_axis()` subsumes `get_abs_positions_ignore_gap` /
      `get_abs_positions_keep_gap` (HATCHet), `build_wl_coords` (Copytyping), and
      `_genome_coords` (UGP).
- [x] Returns a `GenomeAxis` dataclass instead of 3-to-9 tuple returns. Copytyping's dict
      was the superset; `chr_end`, `chr_boundaries`, `xtick_chrs`, `xlab_chrs` are
      properties, so nothing can drift out of sync with `ch_coords`.
- [x] Both layouts in one call: `wl_segments=` (collapsed) or `chrom_sizes=` (linear),
      with `chr_shift` padding either. Passing both or neither raises.
- [x] `contain=` selects bin membership: True = fully inside (HATCHet CNV profiles),
      False = overlap and clip (Copytyping scatter/heatmap). One code path serves both,
      since clipping is a no-op under containment.
- [x] Both input shapes: `#CHR/START/END` bins and `#CHR/POS` sites.
- [x] `decorate_genome_axis()`, `draw_chr_boundaries()`, `draw_segment_boundaries()`,
      and `shade_regions()` replace the three inline decorators.
- [x] Ported `get_chr_sizes`, `get_transparency`, `adaptive_dot_size`, `FigureSaver`,
      and `read_bed_by_chr`. No host-repo imports remain.

Two deliberate behavior changes, both documented in the docstrings:

- `seg_coords` follows Copytyping's rule, which also marks a removed leading telomere
  (first whitelist segment starting past 0). HATCHet marked only internal joins, so a
  leading-gap whitelist now yields one extra dashed line. Confirmed by the test.
- `get_transparency` takes `cols=("RD", "BAF")` instead of hardcoding those columns.
- `shade_regions` raises on a collapsed axis, where a raw interval has no single span.
  UGP drew chromosome labels via `ax.text`; cnplot uses real ticks.

## M2. `cnplot_colormap.py` - palettes - DONE

- [x] `get_cn_colors()` - was byte-identical in both repos; now single-source. Palette and
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

## M3. `cnplot_cnp.py` - integer CN profiles and legends

The 6-function hard-duplication block (plan.md section 3).

- [ ] Unify `plot_cnv_profile`: HATCHet `(bin_info, regions)` vs Copytyping
      `(seg_cnprofile, wl_segments)` - same semantics, drifted names. Build on `GenomeAxis`.
- [ ] Unify `plot_ascn_profile` (same drift).
- [ ] Unify `plot_cnv_legend` / `plot_ascn_legend`; keep Copytyping's `has_mirror` extension
      (`cnp_has_mirror`, `_draw_mirror_swatch`) opt-in so HATCHet callers are unaffected.
- [ ] Keep Copytyping's optional `PI_VIOL` per-bin violation overlay behind a flag.
- [ ] Port `_cnp_segment_geometry`, `_clone_ylabels`.
- [ ] Document the CNP string contract: `";"`-joined per clone, `"a|b"` per clone, first
      field is normal, `PROPS` gives clone proportions.
- [ ] Migration shims in HATCHet3 and Copytyping, then delete the local copies. This is the
      payoff - do not stop at M3 without the deletion.

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
      fixtures for CNP strings, whitelist segments, and cell matrices. M1/M2 were verified
      by a throwaway differential script against the inlined originals; fold that approach
      into `tests/` so the equivalence keeps being checked.
- [ ] Docs: README API section + an `examples/` notebook per module.
- [ ] `CHANGELOG.md` (empty): start at 0.1.0.
- [ ] Publish: TestPyPI first, then PyPI trusted publishing via GitHub Actions on tag.

---

## Open decisions

- [ ] **Column-name contract.** Everything assumes `#CHR/START/END` and the `CNP` string
      format. Bake this in as the documented input schema, or add an adapter layer? Baking
      it in is simpler and matches both consumers today.
- [ ] **Axes-in vs figure-out.** HATCHet/Copytyping primitives take `ax`; UGP functions
      create figures and write files. Standardize on `ax`-taking primitives plus thin page
      builders.
- [ ] **Migration sequencing.** Land 0.1.0 with M1-M3, cut over both host repos, then
      continue - rather than porting everything before any repo depends on it.
