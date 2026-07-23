# cnplot TODO

`cnplot` is a PyPI-installable library for copy-number plotting, factored out of three
sibling repos:

- HATCHet3 - `hatchet-long-read/src/hatchet/plot/` (bulk clone CN, CNP-string driven)
- Copytyping - `Copytyping-dev/src/copytyping/plot/` (single-cell CN, cell-matrix driven)
- UGP - `universal-genotyping/workflow/scripts/plot_utils.py` (reference only, not a consumer)

The build (M0-M7) is complete and staged as 0.1.0 - see [CHANGELOG.md](CHANGELOG.md) for the
delivered surface and the deliberate behavior changes vs the origin repos. What remains, in
order: the downstream migration (M9), then publishing (M8). **Publish last** - only after a
real HATCHet/Copytyping run drives the API do we freeze it on PyPI, since a signature or
behavior gap is cheap to fix before the first upload and expensive after. The migration is the
current focus.

Ground rule (unchanged): cnplot never modifies or deletes files in the sibling repos. The
migration produces shims and proposals; removing superseded local copies is each repo owner's
call.

---

## Module layout (for the migration)

| Module | Public surface |
|---|---|
| `cnplot_io_utils` | `read_chr_sizes`, `read_bed` |
| `cnplot_colormap` | `get_cn_cmap`, `get_ascn_cmap`, `get_baf_cmap`, `get_log2rdr_cmap`, `get_mixcn_cmap`, `get_multiclass_cmap`, `set_palette` |
| `cnplot_genome_axis` | `GenomeAxis`, `Gap`, `Segment`, `BinCoords` |
| `cnplot_utils` | column constants, axis decoration, marker/color/ylim resolvers, `format_clone_name`, `shade_regions`, `FigureSaver` |
| `cnplot_intcnp` | `plot_cnv_profile` (draws its own legend) |
| `cnplot_1d` | `plot_scatter_1d` |
| `cnplot_2d` | `plot_scatter_2d`, `get_landmarks` |
| `cnplot_heatmap` | `plot_heatmap`, `plot_column_strips`, `plot_strip_legend` |
| `cnplot_figures` | `plot_scatter_1d_multisample`, `plot_heatmap_cnp`, `make_row_spec` |

Everything is re-exported from the top-level `cnplot` package, so a call site migrates to
`from cnplot import ...` regardless of which module a name lives in.

## M9. Migration into HATCHet3 and Copytyping - do first

Precondition: the build (M0-M7) is done; runs against the local `0.1.0.dev0` (editable
install or path), before publishing. The question this milestone answers is **how a host repo
adopts cnplot without breaking its call sites**. Three steps per repo, in order.

### Step 1 - shim the local copies (call sites unchanged)

For each function cnplot replaces, replace the local body with a re-export:

    # Copytyping-dev/src/copytyping/plot/plot_common.py
    from cnplot import get_baf_cmap, get_multiclass_cmap, plot_cnv_profile  # was local defs

Keep the local module name and the imported symbols, so nothing downstream of the shim
changes yet. This is the reversible step and where the equivalence check runs.

- [ ] **Equivalence report, per replaced function.** Render the host original and the cnplot
      version on the same real data and diff the figure. Confirm match or record the change as
      one of the deliberate divergences below. Flag to Copytyping specifically: cnplot fixes
      the `seg_coords` leading-gap bug, so a single-segment chromosome whose whitelist starts
      past 0 loses a stray dashed line. This subsumes the old M1/M2 differential idea (a
      unit-level parity check against `build_genome_axis` / the original palettes) - keep any
      still-useful cases as golden tests in `tests/`, rewritten against `GenomeAxis`; drop the
      rest.

Reuse map (what each repo's local code becomes):

| Host local | cnplot replacement |
|---|---|
| HATCHet `get_abs_positions_*`, Copytyping `build_wl_coords`, UGP `_genome_coords` | `GenomeAxis` |
| HATCHet/Copytyping `get_cn_cmap`, `get_ascn_cmap`, `get_baf_cmap` | same names in `cnplot_colormap` |
| HATCHet `build_mixture_cn_cmap` | `get_mixcn_cmap` |
| Copytyping `build_label_cmaps` | `get_multiclass_cmap` |
| the 6-function integer-CN profile block (plan.md sec 3) | `plot_cnv_profile` |
| HATCHet `plot_1d` / `plot_scatter_1d_pseudobulk` | `plot_scatter_1d` |
| HATCHet `plot_2d`, Copytyping 2D pseudobulk | `plot_scatter_2d` |
| Copytyping `plot_label_strips` / `draw_label_legends` | `plot_column_strips` / `plot_strip_legend` |
| the RDR+BAF page builders (all three repos) | `plot_scatter_1d_multisample` |

### Step 2 - absorb the intentional signature / behavior changes

None of these are mechanical; each is a real edit at the call site or a visible figure change.

- [ ] Legends are private. Callers pass `ax_leg=` to `plot_cnv_profile` instead of calling a
      legend function themselves - HATCHet in two places (`plot_utils.py:244-252`), Copytyping
      in one (`plot_scatter_1d.py:437-440`).
- [ ] Profile and scatter functions return None, not the axes. Remove any `ax = plot_...(...)`
      chaining (none found in the host repos, but confirm).
- [ ] `contain` is a parameter, default unchanged at True.
- [ ] The mirrored swatch disappears from any figure whose data has no mirrored bin, and the
      mirror rule itself changed (clones-disagree, not single-clone `b > a`). Expect diffs on
      figures with uniform `1|2` or `1|1`+`1|2` bins.
- [ ] Heatmap side strips take one ordered `strips` list, not `strip_label_map` + `dist_*`
      kwargs. Copytyping's `plot_cnv_heatmap` page builder must build the list.
- [ ] cmap builders now return `(cmap, norm, ticks)`; a caller unpacking `(cmap, norm)` breaks.

### Step 3 - hand off removal + wiring proposals

- [ ] Propose the re-export shims to each repo owner; removing the superseded local modules is
      their decision.
- [ ] **Propose `--region_bed` for HATCHet `cluster_bins`.** Its two `plot_rdr_baf` calls
      (`cluster_bins.py:163`, `:456`) are the only ones reaching `plot_1d` with `regions=None`,
      because `cluster_bins` has no region-BED argument (only `genome_size`). `--region_bed` is
      already `required=True` on `compute_cn` / `plot_cn` / `plot_panel`
      (`hatchet_parser.py:622, 705, 827`), so real pipelines already supply the file. Adding the
      argument lets those diagnostics use `GenomeAxis` unchanged. It is a CLI change for anyone
      running `hatchet cluster_bins` standalone - recommend `required=True` to match the
      neighbours, but that is the owner's call.

### Deferred ports (single-consumer; migrate with the repo, not speculatively)

- [ ] Copytyping `plot_scatter_2d_per_cell` - refactor to build on `get_landmarks` /
      `_annotate_landmarks` and the shared color/marker helpers; the cross-tab layout and the
      BAF-histogram fallback stay Copytyping's. It has no cross-repo duplication, so it moves
      only when Copytyping adopts cnplot.
- [ ] Copytyping heatmap reducers (`_row_layout`, `_aggregate_columns`, `prepare_rdr/baf/pi_gk`,
      `plot_cnv_heatmap`) - stay upstream (need `model_utils` / `anns`); the page builder rewires
      to call `plot_heatmap_cnp`.
- [ ] HATCHet multi-solution CNP grids into `cnplot_figures.py` (`plot_cnp_panel.run`,
      `plot_pool_cnp`, `plot_summary_pdf` / `plot_bars`, `plot_scaling_2d`) - strip the
      solution-loading glue, keep only the drawing.
- [ ] HATCHet `render_cnt_tree` (clone tree) as `cnplot_tree.py` - check the drawing backend
      first; a non-matplotlib dependency would put it out of scope.
- [ ] Re-scoped-out items to revisit only after a repo adopts cnplot: HATCHet `plot_clusters`
      (`_is_multimodal`) - needs scipy, reopen only if served from cnplot.

## M8. Publish - last, after M9 verifies downstream

Done during the build (M0-M7), kept here as the release checklist:

- [x] `__init__.py`: explicit re-exports + `__all__`.
- [x] No `from hatchet.*` / `from copytyping.*` / `import *` in `src/` (verified).
- [x] Type hints + Google-style docstrings on the public surface.
- [x] Tests under `tests/` (real T2T profile simulator; renders through `Agg`).
- [x] `CHANGELOG.md` started at 0.1.0.
- [x] Docs: README gallery (6 figures) + `docs/reference.md` input-format tables +
      `examples/plot_gallery.py`. Per-module notebooks were dropped in favour of the one
      gallery script; revisit only if a notebook is actually requested.

Remaining, gated on M9:

- [ ] Precondition: at least one host repo migrated onto the local build and its figures
      verified (M9 equivalence report), so any signature or behavior gap is fixed before the
      first upload freezes the API.
- [ ] Bump `__version__` off `.dev0`.
- [ ] `twine check dist/*` (twine is in the `dev` extra, not `base`).
- [ ] Publish: TestPyPI first, then PyPI trusted publishing via GitHub Actions on tag.

---

## Open decisions

- [ ] **Column-name contract.** Everything assumes `#CHR/START/END` and the CNP-string format.
      Bake it in as the documented schema (now in `docs/reference.md`) or add an adapter layer?
      Baking it in is simpler and matches both consumers.
- [x] **Axes-in vs figure-out.** Resolved: `ax`-taking primitives plus thin figure builders.
- [x] **Migration sequencing.** Resolved: build the whole package first, migrate once at the end.
