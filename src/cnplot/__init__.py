"""cnplot: plotting library for copy-number data.

The public API is re-exported here, so ``from cnplot import GenomeAxis,
plot_cnv_profile`` works without reaching into submodules. Each name comes from
one module's ``__all__``; the seg.ucn readers and other private helpers stay in
their modules.
"""

from .cnplot_1d import plot_scatter_1d
from .cnplot_2d import annotate_landmarks, get_landmarks, plot_scatter_2d
from .cnplot_colormap import (
    BAF_COLORS,
    CELLTYPE_CMAP,
    DATASET_CMAP,
    NA_COLOR,
    NA_LABELS,
    NORMAL_COLOR,
    POSTERIOR_CMAP,
    PURITY_CMAP,
    get_ascn_cmap,
    get_baf_cmap,
    get_categorical_cmap,
    get_cn_cmap,
    get_log2rdr_cmap,
    get_mixcn_cmap,
    get_multiclass_cmap,
    set_palette,
)
from .cnplot_figures import (
    make_row_spec,
    plot_heatmap_cnp,
    plot_scatter_1d_multisample,
)
from .cnplot_genome_axis import (
    DEFAULT_EXCLUDED_CHROMS,
    BinCoords,
    Gap,
    GapKind,
    GenomeAxis,
    Segment,
)
from .cnplot_heatmap import plot_column_strips, plot_heatmap, plot_strip_legend
from .cnplot_intcnp import plot_cnv_profile
from .cnplot_io_utils import read_bed, read_chr_sizes
from .cnplot_utils import (
    CN_PREFIX,
    EXP_PREFIX,
    MARKER_SIZE_LARGE,
    MARKER_SIZE_SMALL,
    MAX_NDOTS,
    NORMAL_CLONE,
    PI_VIOL_COL,
    SAMPLE_COL,
    U_PREFIX,
    FigureSaver,
    adaptive_dot_size,
    decorate_genome_axis,
    draw_chr_boundaries,
    draw_segment_boundaries,
    get_clone_ylabels,
    get_transparency,
    resolve_colors,
    resolve_marker_size,
    resolve_ylim,
    resolve_ylim_scaled,
    shade_regions,
)

__version__ = "0.1.0"

__all__ = [
    "BAF_COLORS",
    "CELLTYPE_CMAP",
    "CN_PREFIX",
    "DATASET_CMAP",
    "DEFAULT_EXCLUDED_CHROMS",
    "EXP_PREFIX",
    "MARKER_SIZE_LARGE",
    "MARKER_SIZE_SMALL",
    "MAX_NDOTS",
    "NA_COLOR",
    "NA_LABELS",
    "NORMAL_CLONE",
    "NORMAL_COLOR",
    "PI_VIOL_COL",
    "POSTERIOR_CMAP",
    "PURITY_CMAP",
    "SAMPLE_COL",
    "U_PREFIX",
    "BinCoords",
    "FigureSaver",
    "Gap",
    "GapKind",
    "GenomeAxis",
    "Segment",
    "__version__",
    "adaptive_dot_size",
    "annotate_landmarks",
    "decorate_genome_axis",
    "draw_chr_boundaries",
    "draw_segment_boundaries",
    "get_ascn_cmap",
    "get_baf_cmap",
    "get_categorical_cmap",
    "get_clone_ylabels",
    "get_cn_cmap",
    "get_landmarks",
    "get_log2rdr_cmap",
    "get_mixcn_cmap",
    "get_multiclass_cmap",
    "get_transparency",
    "make_row_spec",
    "plot_cnv_profile",
    "plot_column_strips",
    "plot_heatmap",
    "plot_heatmap_cnp",
    "plot_scatter_1d",
    "plot_scatter_1d_multisample",
    "plot_scatter_2d",
    "plot_strip_legend",
    "read_bed",
    "read_chr_sizes",
    "resolve_colors",
    "resolve_marker_size",
    "resolve_ylim",
    "resolve_ylim_scaled",
    "set_palette",
    "shade_regions",
]
