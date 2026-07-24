# cnplot TODO

`cnplot` is the copy-number plotting library factored out of three sibling repos:

- HATCHet3 - `hatchet-long-read/src/hatchet/plot/` (bulk clone CN)
- Copytyping - `Copytyping-dev/src/copytyping/plot/` (single-cell CN)
- UGP - `universal-genotyping/workflow/scripts/plot_utils.py` (reference only, not a consumer)

The build (0.1.0) is done - see [CHANGELOG.md](CHANGELOG.md) for the delivered surface and the
deliberate behavior changes vs the origin repos. What remains: finish the migration (M9), then
publish (M8). **Publish last** - only after a real host run drives the API do we freeze it.

Ground rule: cnplot never modifies sibling-repo files on its own. Migration produces edits in
those repos only when that repo's migration is the active task; removing superseded local copies
is each owner's call.

---

## M9. Migration into the host repos

### Copytyping - done (uncommitted, commits separately)

Migrated onto cnplot: `plot_heatmap`, `plot_scatter_1d`, `plot_scatter_2d`, `plot_common`,
`inference`, `utils`, parser; `plot_copynumber.py` deleted (fully redundant). ascn-profile
support dropped throughout. Reducers that need `model_utils` / `anns` stay upstream.

- [ ] `plot_spatial.py` - last un-migrated module; already imports several cnplot helpers.
      Verify it needs nothing further, then it is done.
- [ ] Equivalence report on real data, then commit the Copytyping side.

### HATCHet - not started

Reuse map (host local -> cnplot):

| Host local | cnplot replacement |
|---|---|
| `get_abs_positions_*` | `GenomeAxis` |
| `get_cn_cmap`, `get_ascn_cmap`, `get_baf_cmap` | same names in `cnplot_colormap` |
| `build_mixture_cn_cmap` | `get_mixcn_cmap` |
| the 6-function integer-CN profile block | `plot_cnv_profile` |
| `plot_1d` / `plot_scatter_1d_pseudobulk` | `plot_scatter_1d` |
| `plot_2d` | `plot_scatter_2d` |
| the RDR+BAF page builders | `plot_scatter_1d_multisample` |

Signature/behavior changes to absorb at the call sites (the *why* is in CHANGELOG Notes):

- [ ] Legends are private: pass `ax_leg=` to `plot_cnv_profile` (HATCHet `plot_utils.py:244-252`).
- [ ] Profile/scatter functions return None, not the axes - drop any `ax = plot_...()` chaining.
- [ ] cmap builders return `(cmap, norm, ticks)`; a caller unpacking `(cmap, norm)` breaks.
- [ ] Heatmap side strips take one ordered `strips` list, not `strip_label_map` + `dist_*` kwargs.
- [ ] Equivalence report per replaced function on real data; keep still-useful parity cases as
      golden tests under `tests/`, rewritten against `GenomeAxis`.
- [ ] Propose `--region_bed` for `cluster_bins` - its two `plot_rdr_baf` calls
      (`cluster_bins.py:163`, `:456`) are the only ones reaching `plot_1d` with `regions=None`.
      Already `required=True` on `compute_cn` / `plot_cn` / `plot_panel`, so pipelines supply it.

### Deferred ports (single-consumer; move only when that repo adopts cnplot)

- [ ] HATCHet multi-solution CNP grids -> `cnplot_figures.py` (strip the solution-loading glue).
- [ ] HATCHet `render_cnt_tree` -> `cnplot_tree.py` (check the drawing backend first; a
      non-matplotlib dep puts it out of scope).
- [ ] HATCHet `plot_clusters` (`_is_multimodal`) - needs scipy; reopen only if served from cnplot.

## M8. Publish - after M9 verifies downstream

PyPI is the source of truth; conda-forge is generated from the PyPI sdist. Strict order:
**M8a must be live before M8b**. M8c wires the host repos onto the release.

### M8a. PyPI (setup ~half a day; each release ~minutes)

- [ ] Precondition: at least one host repo migrated onto the local build and figures verified,
      so any signature/behavior gap is fixed before the first upload freezes the API.
- [ ] Confirm every runtime dep (numpy, pandas, matplotlib, seaborn, adjustText) has a
      conda-forge feedstock covering `python >=3.10`. Gates M8b; catch a PyPI-only dep now.
- [ ] Bump `__version__` off `.dev0` to `0.1.0`.
- [ ] `python -m build`, then `twine check dist/*` (twine is in the `dev` extra).
- [ ] Register on TestPyPI + PyPI as GitHub trusted publishers (OIDC, no stored token).
- [ ] Release workflow: build -> `twine check` -> TestPyPI on pre-release tag, PyPI on `v*` tag.
      Smoke-test by installing from TestPyPI and running `examples/plot_gallery.py`.
- [ ] Tag `v0.1.0`; installable via `pip` within minutes of the workflow finishing.

### M8b. conda-forge (first PR days to ~2 weeks; updates same-day)

- [ ] `grayskull pypi cnplot` from the published sdist; hand-check the `run:` deps + `python`
      floor. Set `noarch: python` (pure-Python).
- [ ] PR to `conda-forge/staged-recipes`. Bot lint + CI is minutes; **human review is the long
      pole - a few days, up to ~2 weeks on first submission**.
- [ ] On merge the feedstock is auto-created and builds; verify
      `conda install -c conda-forge cnplot`. Thereafter `regro-cf-autotick-bot` opens an update
      PR within hours of each PyPI release - same-day merge, no manual recipe edits.

### M8c. Wire into HATCHet / Copytyping (~1 day per repo + their bioconda review)

- [ ] pip metadata: add `cnplot>=0.1,<1.0` to each repo's `pyproject.toml`.
- [ ] bioconda recipe: add `- cnplot >=0.1,<1.0` under `requirements: run:` in each `meta.yaml`.
- [ ] Never let cnplot arrive via `pip` alongside conda-managed numpy/pandas/matplotlib - the
      conda solver does not see pip installs and can clobber them. Once M8b is live, cnplot
      resolves through conda in every conda/container build.

---

## Open decisions

- [ ] **Column-name contract.** Everything assumes `#CHR/START/END` and the CNP-string format,
      now documented in `docs/reference.md`. Bake it in as the schema, or add an adapter layer?
      Baking in is simpler and matches both consumers.
