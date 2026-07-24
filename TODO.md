# cnplot TODO

`cnplot` is the copy-number plotting library factored out of three sibling repos:

- HATCHet3 - `hatchet-long-read/src/hatchet/plot/` (bulk clone CN)
- Copytyping - `Copytyping-dev/src/copytyping/plot/` (single-cell CN)
- UGP - `universal-genotyping/workflow/scripts/plot_utils.py` (reference only, not a consumer)

The build (0.1.0) is done - see [CHANGELOG.md](CHANGELOG.md) for the delivered surface and the
deliberate behavior changes vs the origin repos. The migration (M9) is done and committed in both
consumer repos. What remains: publish (M8).

Ground rule: cnplot never modifies sibling-repo files on its own.

---

## M9. Migration into the host repos - done

Both repos delegate to cnplot on their own branches; verify their figures on a real run before
M8 freezes the API (the M8a precondition).

- **Copytyping** (committed, `730220e` / `3cdb6b2`): every plot module delegates to cnplot;
  `plot_copynumber.py` deleted; ascn-profile support dropped; `plot_spatial` on cnplot.
- **HATCHet** (committed, `7c85ea8`): `plot_utils` builds a `GenomeAxis`; `plot_cn`,
  `plot_cluster_bins`, `plot_cnp_panel`, `plot_compute_cn` delegate; `--region_bed` added to
  `cluster_bins`; legends via `ax_leg=`; `(cmap, norm, ticks)` unpacking; heatmap `strips` list.

Kept upstream by design (single-consumer or non-portable deps): HATCHet `_is_multimodal` cluster
gate (scipy), `render_cnt_tree` clone-tree drawing, and the `plot_cnp_panel` multi-solution grid
glue (which already reuses `plot_cnv_profile`).

- [ ] Optional future port, only if a second consumer appears: `render_cnt_tree` ->
      `cnplot_tree.py` (matplotlib-based, so portable).

## M8. Publish

PyPI is the source of truth; conda-forge is generated from the PyPI sdist. Strict order:
**M8a must be live before M8b**. M8c wires the host repos onto the release.

### M8a. PyPI (setup ~half a day; each release ~minutes)

- [ ] Precondition: host figures verified on a real run, so any signature/behavior gap is fixed
      before the first upload freezes the API.
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

- [ ] **Column-name contract.** Everything assumes `#CHR/START/END` and the per-clone
      `cn_<clone>` columns, now documented in `docs/reference.md`. Bake it in as the schema, or
      add an adapter layer? Baking in is simpler and matches both consumers.
