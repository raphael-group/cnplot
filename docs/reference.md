# Input / Output Reference

---

## `region_bed`

Headerless, tab-separated file. Read by `GenomeAxis` (via `read_bed`); the whitelist regions
the genome axis is drawn from.

| Column | Required | Type | Description |
|---|---|---|---|
| `#CHR` | Yes | str | Chromosome. |
| `START` | Yes | int | Region start (0-based). |
| `END` | Yes | int | Region end. |
| name | No | str | 4th column (e.g. arm id); read but not used. |

## `chrom_sizes`

Headerless, whitespace-separated file. Read by `GenomeAxis`.

| Column | Required | Type | Description |
|---|---|---|---|
| chromosome | Yes | str | Chromosome name. Row order sets the axis chromosome order. |
| length | Yes | int | Chromosome length. |

## `seg_df` (HATCHet-style seg.ucn profile)

Integer copy-number profile, one row per segment. Consumed by `plot_cnv_profile` and, in the
2D form below, by `plot_scatter_2d`.

| Column | Required | Type | Description |
|---|---|---|---|
| `SAMPLE` | For multi-sample | str | Sample id; needed only when the table stacks several samples (select with `sample_id=`). |
| `#CHR` | Yes | str | Chromosome. |
| `START` | Yes | int | Start coordinate (0-indexed, inclusive). |
| `END` | Yes | int | End coordinate (exclusive). |
| `cn_(\w+)` | Yes (>= 1) | str | Per-clone allele copy number as `"a\|b"`, one column per clone (`cn_normal`, `cn_clone1`, ...). |
| `u_(\w+)` | No | float | Per-clone proportion in `[0, 1]`; drives the proportion legend. Constant down the column. |
| `PI_VIOL` | No | bool | Per-segment pure-integer-violation flag; drawn as a red/green overlay. Matched case-insensitively. |

## `obs_df` (Bin-level observations)

Bin-level observations. Consumed by `plot_scatter_1d` /
`plot_scatter_1d_multisample` (`ycol`) and `plot_scatter_2d` (`xcol`, `ycol`).

| Column | Required | Type | Description |
|---|---|---|---|
| `#CHR` | Yes | str | Chromosome. |
| `START` | Yes | int | Start coordinate (0-indexed, inclusive). |
| `END` | Yes | int | End coordinate (exclusive). |
| group | For multi-sample | str | Group / sample id; column name set by `group_col=` (default `SAMPLE`). |
| value | Yes (>= 1) | float | The plotted quantity, named freely and passed as `ycol` / `xcol` (e.g. `RD`, `BAF`, `FCN-A`). Used as given - no scaling or log. |
| hue | No | any | Categorical column to color points by, named via `hue=` (e.g. `cnp` holding the joint `"a\|b;a\|b"` state). |
| `keep_col` | No | bool | Marks bins that passed upstream QC; False bins are drawn in `filtered_color`. Column name set by `keep_col=`. |

## `expected_df` (1D overlay)

Expected value per segment (or bin) for the step overlay on `plot_scatter_1d`. May be coarser
than `obs_df` - it maps through the same axis.

| Column | Required | Type | Description |
|---|---|---|---|
| `#CHR` | Yes | str | Chromosome. |
| `START` | Yes | int | Start coordinate (0-indexed, inclusive). |
| `END` | Yes | int | End coordinate (exclusive). |
| `exp_(\w+)_(\w+)` | Yes | float | Expected value, one `exp_<value>_<group>` column per value per group (e.g. `exp_RD_S1`). A missing column simply draws no overlay for that group. |

## `expected_df` (2D landmarks)

Expected copy-number landmarks for `plot_scatter_2d`; a seg.ucn-style table with precomputed
landmark coordinates.

| Column | Required | Type | Description |
|---|---|---|---|
| `SAMPLE` | For multi-sample | str | Group id; selected with `group=`. |
| `#CHR`, `START`, `END` | No | str/int | Provenance only; |
| `cn_(\w+)` | Yes | str | Per-clone `"a\|b"` CN, one column per clone; sets each landmark's `(a, b)` label and dedup key. |
| `u_(\w+)` | No | float | Per-clone proportion; drives the proportion legend. |
| `exp_<xcol>` | Yes | float | Landmark x coordinate (matches the observed `xcol`). |
| `exp_<ycol>` | Yes | float | Landmark y coordinate (matches the observed `ycol`). |
| `is_filtered` | No | bool | Drop this landmark from the plot. |
| `is_balanced` | No | bool | Draw a square marker instead of a circle. |

## Heatmap inputs

`plot_heatmap` / `plot_heatmap_cnp` take a reduced matrix, not a DataFrame, plus a bin table.

| Input | Required | Type | Description |
|---|---|---|---|
| `matrix` | Yes | ndarray | `(n_rows, n_bins)` values; columns aligned to `coords_df` rows, rows bottom-to-top. |
| `coords_df` | Yes | DataFrame | The bins: `#CHR`, `START`, `END`, one row per matrix column. |
| `row_labels` | No | ndarray | `(n_rows,)` label per row; contiguous runs get a boxed y-tick. |
| `strips` | No | list[dict] | Ordered annotation strips left of the heatmap, the first closest to it. Each dict is one of the two strip specs below. |

### Strip specs (entries of `strips`)

A categorical strip colors each row by its value; a distribution strip stacks each
row's distribution (e.g. copy-typing posteriors) into a bar; a continuous strip colors
each row by a scalar (e.g. purity). The kind is told apart by the data key: `values`
marks categorical, `matrix` distribution, `scalar` continuous.

| Key | Kind | Required | Type | Description |
|---|---|---|---|---|
| `name` | all | Yes | str | Strip key; also the default legend title and shared-palette identity. |
| `display_name` | all | No | str | Title shown on the strip and legend; defaults to `name`. |
| `values` | categorical | Yes | ndarray | `(n_rows,)` value per row. |
| `cmap` | categorical | No | dict | `{value: color}`; auto-built (shared across strips) when omitted. |
| `matrix` | distribution | Yes | ndarray | `(n_rows, len(order))` per-row distribution, rows summing to 1. |
| `order` | distribution | Yes | list | Column order of `matrix`, i.e. the stacked categories. |
| `cmap` | distribution | Yes | dict | `{category: color}` for the stacked bars. |
| `props` | distribution | No | dict | `{category: fraction}` for the strip's legend. |
| `scalar` | continuous | Yes | ndarray | `(n_rows,)` scalar per row. |
| `cmap` | continuous | Yes | Colormap | Colormap (name or object) for the scalar. |
| `norm` | continuous | No | Normalize | Scalar normalization; matplotlib autoscales when omitted. |
| `show_cbar` | continuous | No | bool | Draw a colorbar for the scalar in the legend column (default off). |
| `cbar_ticks` | continuous | No | list | Explicit colorbar ticks when `show_cbar` is set. |

Continuous strips carry no legend unless `show_cbar` is set (then a colorbar).
