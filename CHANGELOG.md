# Changelog

## Unreleased

### cnplot package
- Clone names are set in the regular weight: the clone labels beside the integer-CN profile
  and the `<label> (<pct>%)` row-block labels on the heatmap. The CN key's `CNA` title is
  regular as well.
- `BAF_LIM = (-0.01, 1.01)`: default limits for a BAF-like axis on [0, 1], exported from
  `cnplot`. The pad keeps points at exactly 0 or 1 off the spine while the ticks still run
  0.0 to 1.0. Used by the gallery/README BAF row (`ylim`) and the RDR-vs-BAF scatter
  (`xlim`), which previously used `(-0.05, 1.05)` and `(0, 1)`.
- Fonts are configured in one module, `cnplot.cnplot_style`: the family list
  (`FONT_FAMILY`), the embedding types (`FONT_TYPE`), the `FS_*` sizes, and `apply_style()`,
  which writes them into rcParams on import. No other module names a family, a size or an
  embedding type. `FONT_FAMILY` and `apply_style` are exported from `cnplot`.
- One type scale for every figure. The heatmap page, the integer-CN profile, the CN key,
  the proportion legends and the 2D scatter no longer carry their own literal font sizes
  (6 to 15 pt); they take `FS_TICK`, `FS_LEGEND`, `FS_AXLABEL`, `FS_AXTITLE` and `FS_TITLE`,
  and the same sizes are set as rc defaults on import so text drawn without an explicit
  size (a colorbar's ticks) follows. `FS_ANNOT` is new, for labels written on the data.
- Layout measured instead of assumed, so it holds at any font size:
  - the CN key labels each state `(a,b)` and wraps its groups onto further rows until the
    widest label fits under one swatch; it never shrinks the type or the notation;
  - `plot_strip_legend` stacks legends by their rendered height. Its `entry_h` argument is
    removed and `gap` defaults to half a legend line;
  - the heatmap colorbar grows until neighbouring tick labels clear each other, and the
    strip legends start to the right of its rendered tick labels;
  - the `plot_heatmap_cnp` title is anchored above the mesh rather than at a fixed
    figure height.
- Genome-axis labels are placed so they clear each other at the type scale:
  - chromosome names keep their `rotation` (60 degrees by default); names closer together
    than a line at that angle are nudged apart along the axis, and on an axes too narrow
    to hold them in one row, alternate crowded names go on a second row (`chrname_tiers`
    tells a layout which case applies);
  - in Mb mode an Mb label within a line of the previous one is left blank (its tick mark
    stays), the chromosome names under the Mb labels are offset in points rather than as a
    fraction of the axes, and on a stacked page only the last row labels its Mb ticks;
  - the gap between a row and the profile under it is sized in inches from the labels it
    holds (`genome_label_room`) by `plot_scatter_1d_multisample` and `plot_heatmap_cnp`,
    and `plot_heatmap_cnp` no longer opens the same gap between the profile and its key.
- The README gallery is drawn 20 in wide, the width the type scale assumes.
- `FigureSaver.close()` can be called twice on every supported matplotlib.
- Importing cnplot stops matplotlib's PDF backend storing rasterized layers as `/Indexed`
  palette images. With `pdf.compression` on, the backend converts any RGB raster of 256 or
  fewer colours to a palette and packs the indices at 1, 2 or 4 bits per pixel when the
  palette is small, which vector editors misread, so a copy-number raster opened in the
  wrong colours. Rasters are now written at 8 bits per component in DeviceRGB, still
  Flate-compressed. The hook is on `PdfFile._writeImg`, so it applies to every PDF the
  process writes, through `FigureSaver`, `PdfPages` or `Figure.savefig`. Rendered pixels
  are unchanged; a PDF holding such rasters is larger.
  This works around a known matplotlib bug, wrong colours when such a PDF is opened in
  Adobe Illustrator ([matplotlib#27893](https://github.com/matplotlib/matplotlib/issues/27893)),
  for which matplotlib has no setting; remove the hook once that is fixed upstream.
- `FigureSaver` no longer holds a whole-figure raster buffer per rasterized artist.
  matplotlib renders each rasterized artist into a buffer the size of the entire figure
  and passes the PDF backend a view of it, and the backend keeps every image until the
  file is finalized, so a multi-page PDF pinned `n_artists * fig_w_in * fig_h_in * dpi^2`
  bytes until close. `FigureSaver.savefig` now copies each image on the way in, so only
  the cropped pixels are retained. Written output is byte-identical.

## 0.1.1

Backward-compatible additions; existing calls and output are unchanged.

### cnplot package
- `decorate_genome_axis` gains `mb_ticks` / `mb_tick_step`: label x in Mb every
  `mb_tick_step` bp within each chromosome (reset per chromosome), with chromosome names
  drawn as compact off-axis text - the denser mode for genome-wide scatter QC. Plumbed
  through `plot_scatter_1d`, `plot_scatter_1d_multisample`, `plot_heatmap`, and
  `plot_heatmap_cnp`. Default off, 50Mb step; the Mb labels are rotated by `rotation`
  (default 60 degrees) so the dense ticks do not overlap.
- Version lives in the repo-root `VERSION` file: `pyproject.toml` builds from it
  (`dynamic version = {file = "VERSION"}`) and `cnplot.__version__` reads it back, falling
  back to the installed metadata outside a source checkout. One edit per bump.
- `GenomeAxis(region_bed=None)`: build one full-length region per chromosome from
  `chrom_sizes`, i.e. a plain whole-genome axis with no gaps, for callers that only have a
  sizes file.

## 0.1.0

First release of `cnplot`, the shared copy-number plotting library.

### cnplot package
- `GenomeAxis`: a genome coordinate transform built from a region BED and chromosome
  sizes, mapping any bin table by `(#CHR, START, END)`.
- Colormaps: `get_cn_cmap`, `get_ascn_cmap`, `get_baf_cmap`, `get_log2rdr_cmap`,
  `get_mixcn_cmap`, `get_multiclass_cmap`, `set_palette`.
- Integer copy-number profiles: `plot_cnv_profile` over the seg.ucn table.
- 1D genome scatter: `plot_scatter_1d` and the `plot_scatter_1d_multisample` builder,
  with `hue`/palette and an optional `expected_df` overlay.
- 2D RDR/BAF scatter: `plot_scatter_2d` (seaborn `JointGrid`) with copy-number landmarks.
- Single-cell heatmap: `plot_heatmap` / `plot_heatmap_cnp` with ordered categorical and
  posterior side `strips` and an optional value colorbar.
