# CN-Plot
`cnplot` is an allele-specific copy-number profile visualization python package implemented based on matplotlib. This package is developed for the purpose of easily and beautifully plotting copy-number profiles as well as read-depth ratio, B-allele frequency (and others) observations across multiple samples with a common reference coordinates with both single-cell and (pseudo-)bulk option.


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
See **[docs/reference.md](docs/reference.md)** for input format descriptions.

## Galleries

Every figure is produced by [`examples/plot_gallery.py`](examples/plot_gallery.py) from a
small bundled dataset; the snippets below show the essential call.

### 0. Load the data and build the axis

```python
import pandas as pd
from cnplot import GenomeAxis, get_mixcn_cmap

# build a common coordinate axis from the reference genome
genome_axis = GenomeAxis("regions.bed", "chrom.sizes")

seg_df = pd.read_table("sample.seg.ucn.tsv")   # the copy-number profile
obs_df = pd.read_table("bins.tsv")             # bin-level RDR / BAF observations

# color points by the joint CNP, so a subclonal segment stays distinct from a
# clonal one at the same total copy number
palette = get_mixcn_cmap(obs_df["cnp"].unique())
```

### 1. Integer copy-number profile

```python
import matplotlib.pyplot as plt
from cnplot import plot_cnv_profile

fig, (ax, ax_leg) = plt.subplots(2, 1, figsize=(12, 3), height_ratios=[3, 1])
plot_cnv_profile(ax, seg_df, genome_axis, sample_id="S1", ax_leg=ax_leg)
```

![profile](examples/profile.png)

### 2. Genome-wide RDR and BAF (multi-sample)

```python
from cnplot import make_row_spec, plot_scatter_1d_multisample

rows = [
    make_row_spec("RD", ylabel="RDR", ylim=(0, 3), href=1.0),
    make_row_spec("BAF", ylabel="mhBAF", ylim=(-0.05, 1.05), href=0.5),
]
fig = plot_scatter_1d_multisample(
    obs_df, genome_axis, rows, expected_df=expected_1d,
    hue="cnp", palette=palette, seg_df=seg_df,
)
```

![scatter_1d](examples/scatter_1d.png)

### 3. RDR vs BAF joint scatter

```python
from cnplot import plot_scatter_2d

grid = plot_scatter_2d(
    obs_df, "BAF", "RD", expected_df=expected_2d, group="S1",
    hue="cnp", palette=palette,
    xlim=(0, 1), ylim=(0, 3), xlabel="mhBAF", ylabel="RDR",
)
```

![scatter_2d](examples/scatter_2d.png)

### 4. Allele-specific fractional copy number (FCN-A / FCN-B)

Derive the major/minor fractional copy numbers, then mirror the minor allele (FCN-B)
below the major (FCN-A) with `reverse_y=True`:

```python
obs_df["FCN-A"] = 2 * obs_df["RD"] * (1 - obs_df["BAF"])
obs_df["FCN-B"] = 2 * obs_df["RD"] * obs_df["BAF"]
rows = [
    make_row_spec("FCN-A", ylabel="FCN-A", ylim=(0, 3.5), href=1.0),
    make_row_spec("FCN-B", ylabel="FCN-B", ylim=(0, 3.5), href=1.0, reverse_y=True),
]
fig = plot_scatter_1d_multisample(
    obs_df, genome_axis, rows, expected_df=expected_ab,
    hue="cnp", palette=palette, seg_df=seg_df,
)
```

![fcn_ab](examples/fcn_ab.png)

### 5. Single-cell copy-number heatmap

One call draws the mesh, its colorbar, the categorical / posterior side strips, and the
integer-CN profile with its legend. Shown below for RDR and BAF:

```python
from cnplot import get_log2rdr_cmap, plot_heatmap_cnp

cmap, norm, ticks = get_log2rdr_cmap()   # or get_baf_cmap() for BAF
fig = plot_heatmap_cnp(
    matrix, bins, genome_axis, seg_df, sample_id="S1",
    row_labels=cell_labels, cmap=cmap, norm=norm,
    cbar_label="RDR", cbar_ticks=ticks,
    strips=[
        # posterior distribution strip, then a categorical strip
        {"name": "Copy-typing", "matrix": posteriors, "order": clones,
         "cmap": clone_cmap, "props": props},
        {"name": "cell_type", "values": celltype, "display_name": "Cell-type"},
    ],
)
```

![heatmap_rdr](examples/heatmap_rdr.png)

![heatmap_baf](examples/heatmap_baf.png)