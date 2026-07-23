# CN-Plot

## Installation
`Python >=3.10` is required for installation.
### Install PyPI package:
```sh
pip install cnplot
```

### Install from source:
```sh
git clone https://github.com/raphael-group/cnplot.git
cd cnplot
pip install -e ".[dev]"
```

## API Usage

Build a `GenomeAxis` once from a region BED and a chromosome-sizes file, then reuse it
across every plot:

```python
from cnplot.cnplot_genome_axis import GenomeAxis

axis = GenomeAxis("regions.bed", "chrom.sizes")
```

Each plot takes a DataFrame in the `seg.ucn` layout (`#CHR`, `START`, `END`, `cn_<clone>`,
`u_<clone>`) or per-bin observations, and maps them onto the shared axis by coordinate. See
the Galleries below for the essential call behind each figure.

## Galleries

Every figure is produced by [`examples/plot_gallery.py`](examples/plot_gallery.py) from a
small bundled dataset; the snippets below show the essential call.

### 1. Integer copy-number profile

```python
import matplotlib.pyplot as plt
from cnplot.cnplot_intcnp import plot_cnv_profile

fig, (ax, ax_leg) = plt.subplots(2, 1, height_ratios=[3, 1])
plot_cnv_profile(ax, seg_df, axis, sample_id="S1", ax_leg=ax_leg)
```

![profile](examples/profile.png)

### 2. Genome-wide RDR and BAF (multi-sample)

```python
from cnplot.cnplot_figures import make_row_spec, plot_scatter_1d_multisample

rows = [
    make_row_spec("RD", ylabel="RDR", href=1.0),
    make_row_spec("BAF", ylabel="mhBAF", href=0.5),
]
fig = plot_scatter_1d_multisample(
    obs_df, axis, rows, expected_df=expected_df,
    hue="cnp", palette=palette, seg_df=seg_df,
)
```

![scatter_1d](examples/scatter_1d.png)

### 3. RDR vs BAF joint scatter

```python
from cnplot.cnplot_2d import plot_scatter_2d

grid = plot_scatter_2d(
    obs_df, "BAF", "RD", expected_df=seg_df, group="S1",
    hue="cnp", palette=palette,
)
```

![scatter_2d](examples/scatter_2d.png)

### 4. Allele-specific fractional copy number (FCN-A / FCN-B)

The minor allele (FCN-B) is mirrored below the major (FCN-A) with `reverse_y=True`:

```python
rows = [
    make_row_spec("FCN-A", href=1.0),
    make_row_spec("FCN-B", href=1.0, reverse_y=True),
]
fig = plot_scatter_1d_multisample(
    obs_df, axis, rows, expected_df=expected_df,
    hue="cnp", palette=palette, seg_df=seg_df,
)
```

![fcn_ab](examples/fcn_ab.png)

### 5. Single-cell copy-number heatmap

One call draws the mesh, its colorbar, the categorical / posterior side strips, and the
integer-CN profile with its legend. Shown below for RDR and BAF:

```python
from cnplot.cnplot_figures import plot_heatmap_cnp

fig = plot_heatmap_cnp(
    matrix, bins, axis, seg_df, sample_id="S1",
    row_labels=cell_labels, cmap="coolwarm", norm=norm,
    cbar_label="RDR", strip_label_map={"cell_type": celltype},
    dist_strip=("Copy-typing", posteriors, clones, clone_cmap, props),
)
```

![heatmap_rdr](examples/heatmap_rdr.png)

![heatmap_baf](examples/heatmap_baf.png)