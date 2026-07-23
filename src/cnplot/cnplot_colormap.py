"""Color palettes for copy-number and categorical label plotting.

Three independent palettes:

- joint CN states (a, b) -> :func:`get_cn_cmap`, shared by profiles, legends,
  and scatter points so all three agree.
- per-allele CN -> :func:`get_ascn_cmap`, a sequential ramp for A/B panels.
- categorical labels -> :func:`build_label_cmaps`, which colors several
  annotation strips from one shared palette and keeps "normal" and missing
  labels neutral.

Functions taking ``invalid_labels`` or ``na_labels`` default to the module
constants but accept caller-supplied sets, since what counts as "missing" belongs
to the caller's data.
"""

import re

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

__all__ = [
    "BAF_COLORS",
    "CELLTYPE_CMAP",
    "DATASET_CMAP",
    "INVALID_LABELS",
    "NA_CELLTYPE",
    "NA_COLOR",
    "NORMAL_COLOR",
    "POSTERIOR_CMAP",
    "PURITY_CMAP",
    "build_mixture_cn_cmap",
    "build_label_cmaps",
    "get_ascn_cmap",
    "get_cn_cmap",
    "get_baf_cmap",
    "set_palette",
]


# =============================================================================
# Palette constants
# =============================================================================

NORMAL_COLOR = "lightgray"
NA_COLOR = "darkgray"

INVALID_LABELS = frozenset({"Doublet", "doublet", "Unknown", "NA"})
NA_CELLTYPE = frozenset({"Unknown", "NA"})

PURITY_CMAP = "magma_r"
POSTERIOR_CMAP = "viridis"
CELLTYPE_CMAP = "Set2"
DATASET_CMAP = "Set1"

# Diverging BAF palette, 10 discrete bins: blue = A-skewed, gray = balanced,
# red = B-skewed.
BAF_COLORS = [
    "#1f77b4",
    "#3b8bc6",
    "#67a9cf",
    "#90c4d6",
    "#b8d6da",
    "#d9d9d9",
    "#fddbc7",
    "#f4a582",
    "#d6604d",
    "#b2182b",
]

_TUMOR_COLORS = [mcolors.to_hex(c) for c in plt.get_cmap("tab10").colors]
_TUMOR_COLORS_20 = [mcolors.to_hex(c) for c in plt.get_cmap("tab20").colors]

# Integer CN palette, keyed by (major, minor). Hue tracks total CN and shade
# tracks allelic imbalance, so LOH reads darker than a balanced state of the
# same total.
_CN_PALETTE = {
    (0, 0): "darkblue",
    (1, 0): "lightblue",
    (1, 1): "lightgray",
    (2, 0): "dimgray",
    (2, 1): "khaki",
    (3, 0): "gold",
    (2, 2): "navajowhite",
    (3, 1): "orange",
    (4, 0): "darkorange",
    (3, 2): "salmon",
    (4, 1): "red",
    (5, 0): "darkred",
    (3, 3): "plum",
    (4, 2): "orchid",
    (5, 1): "purple",
    (6, 0): "indigo",
    (4, 3): "#c0b7f0",
    (5, 2): "#a485f4",
    (6, 1): "#6f42c1",
    (7, 0): "#4b0082",
}

# teal-green, distinct from every CN <= 7 color above
_CN_DEFAULT_COLOR = "#00cc99"

_CN_COPY_STATES = [
    (1, 0),
    (0, 1),
    (0, 2),
    (1, 1),
    (2, 0),
    (0, 3),
    (1, 2),
    (2, 1),
    (3, 0),
    (0, 4),
    (1, 3),
    (2, 2),
    (3, 1),
    (4, 0),
    (0, 5),
    (1, 4),
    (2, 3),
    (3, 2),
    (4, 1),
    (5, 0),
    (0, 6),
    (1, 5),
    (2, 4),
    (3, 3),
    (4, 2),
    (5, 1),
    (6, 0),
    (0, 7),
    (1, 6),
    (2, 5),
    (3, 4),
    (4, 3),
    (5, 2),
    (6, 1),
    (7, 0),
]


# =============================================================================
# Helpers
# =============================================================================


def _distinct_subclonal_colors(n: int, avoid_rgb: np.ndarray, min_dist: float = 0.22):
    """Pick colors that stay clear of the integer-CN palette.

    Greedily walks a husl pool, rejecting candidates within ``min_dist`` (RGB
    Euclidean) of an avoided or already-picked color, then relaxes to the
    furthest remaining if the pool runs short.

    Args:
        n: Number of colors to return.
        avoid_rgb: (c, 3) RGB colors to stay away from.
        min_dist: Minimum RGB distance from avoided and already-picked colors.

    Returns:
        List of ``n`` RGB tuples, or an empty list when ``n`` <= 0.
    """
    if n <= 0:
        return []
    pool = sns.color_palette("husl", n_colors=max(3 * n, 24))
    cand = np.array([mcolors.to_rgb(c) for c in pool])

    def min_dist_to(rgb, ref):
        """Measure the closest approach from one color to a reference set.

        Args:
            rgb: (3,) RGB color to measure from.
            ref: (k, 3) reference colors.

        Returns:
            Smallest Euclidean RGB distance, or infinity if ``ref`` is empty.
        """
        if len(ref) == 0:
            return np.inf
        return float(np.sqrt(((ref - rgb) ** 2).sum(axis=1)).min())

    picked, picked_rgb = [], np.empty((0, 3))
    for i in range(len(cand)):
        if min_dist_to(cand[i], avoid_rgb) >= min_dist and (
            min_dist_to(cand[i], picked_rgb) >= min_dist
        ):
            picked.append(pool[i])
            picked_rgb = np.vstack([picked_rgb, cand[i]])
            if len(picked) == n:
                return picked

    remaining = sorted(
        (i for i in range(len(cand)) if pool[i] not in picked),
        key=lambda i: -min_dist_to(cand[i], avoid_rgb),
    )
    for i in remaining:
        picked.append(pool[i])
        if len(picked) == n:
            break
    return picked


def _is_normal_like(label: str) -> bool:
    """Test whether a label names the normal reference population.

    Args:
        label: Label to test, e.g. "normal" or "Normal_cell".

    Returns:
        True if the label starts with "normal", case-insensitively.
    """
    return str(label).lower().startswith("normal")


def _is_colored_label(label: str, invalid_labels, na_labels) -> bool:
    """Test whether a label consumes a palette slot.

    Args:
        label: Label to test.
        invalid_labels: Labels treated as invalid.
        na_labels: Labels treated as missing.

    Returns:
        True unless the label is normal-like, invalid, or missing, all of which
        take fixed grays instead.
    """
    return (
        label not in invalid_labels
        and label not in na_labels
        and not _is_normal_like(label)
    )


def _clone_order_key(label: str, invalid_labels) -> tuple:
    """Build a sort key ordering clone labels for legends.

    Args:
        label: Label to order.
        invalid_labels: Labels treated as invalid.

    Returns:
        Sort key placing normal first, then clone1, clone2, ... numerically,
        then other labels alphabetically, then invalid labels last.
    """
    if label == "normal":
        return (0, 0, "")
    m = re.match(r"clone(\d+)$", label)
    if m:
        return (1, int(m.group(1)), "")
    if label in invalid_labels:
        return (3, 0, label)
    return (2, 0, label)


# =============================================================================
# Copy-number palettes
# =============================================================================


def get_cn_cmap() -> tuple:
    """Return the integer allele-specific copy-number color map and legend order.

    Keyed by (a, b) under both orderings, so (2, 1) and (1, 2) render alike.
    States above total CN 7 fall through to the default color.

    Returns:
        (state_style, tcn_states): the {(a, b): color} map with a "default"
        key, and {total_cn: [(a, b), ...]} in the order legends draw them.
    """
    tcn_states = {}
    for a, b in _CN_COPY_STATES:
        tcn_states.setdefault(int(a + b), []).append((a, b))

    state_style = {}
    for (major, minor), color in _CN_PALETTE.items():
        state_style[(major, minor)] = color
        state_style[(minor, major)] = color
    state_style["default"] = _CN_DEFAULT_COLOR
    return state_style, tcn_states


def get_ascn_cmap() -> tuple:
    """Return the per-allele copy-number color map and legend order.

    A sequential ramp: 0 white, 1 black, 2 upward sampling inferno from 0.65 to
    0.97, which is perceptually uniform and colorblind-safe. Callers draw
    non-zero colors at alpha 0.5, so CN 1 reads as mid gray.

    Returns:
        (state_style, tcn_states): the {cn: color} map for 0-6 with a "default"
        key for cn >= 7, and its sorted integer keys.
    """
    state_style = {
        0: "#FFFFFF",  # white
        1: "#000000",  # black, mid gray at alpha 0.5
        2: "#EA632A",  # inferno 0.65
        3: "#F57D15",  # inferno 0.71
        4: "#FB9B06",  # inferno 0.78
        5: "#FBBA1F",  # inferno 0.84
        6: "#F5D949",  # inferno 0.91
    }
    state_style["default"] = "#F3F68A"  # inferno 0.97, cn >= 7
    tcn_states = sorted(k for k in state_style if isinstance(k, int))
    return state_style, tcn_states


def build_mixture_cn_cmap(
    clone_states,
    clone_props=None,
    display_min_clone_prop: float | None = None,
) -> dict:
    """Map each joint-clone CNP string to a scatter color.

    States where every visible tumor clone shares one (a, b) take that state's
    integer-CN color, so points match the profile and legend. Subclonal states
    get distinct colors kept away from that palette. If the threshold hides every
    clone, the check falls back to all tumor clones.

    Args:
        clone_states: CNP strings "n|n;a|b;..." with the normal clone first.
        clone_props: Per-clone proportions aligned to the CNP fields, index 0
            normal. None treats every clone as visible.
        display_min_clone_prop: Minimum proportion for a tumor clone to count as
            visible. None treats every clone as visible.

    Returns:
        {cnp_string: color} covering every entry of ``clone_states``.
    """
    state_style, _ = get_cn_cmap()

    def visible(j):
        """Test whether clone j is above the display threshold.

        Args:
            j: Clone index into the CNP fields, where 0 is normal.

        Returns:
            True if the clone's proportion is at least
            ``display_min_clone_prop``, or if no threshold applies.
        """
        return (
            display_min_clone_prop is None
            or clone_props is None
            or clone_props[j] >= display_min_clone_prop
        )

    palette = {}
    subclonal = []
    for cs in clone_states:
        tumor = [
            (int(x.split("|")[0]), int(x.split("|")[1])) for x in str(cs).split(";")[1:]
        ]
        vis = [tumor[j - 1] for j in range(1, len(tumor) + 1) if visible(j)]
        if not vis:
            vis = tumor
        if vis and all(p == vis[0] for p in vis):
            palette[cs] = state_style.get(vis[0], state_style["default"])
        else:
            subclonal.append(cs)

    cn_rgb = np.array([mcolors.to_rgb(c) for c in set(state_style.values())])
    subclonal_colors = _distinct_subclonal_colors(len(subclonal), cn_rgb)
    for cs, color in zip(subclonal, subclonal_colors, strict=False):
        palette[cs] = color
    return palette


def get_baf_cmap() -> tuple:
    """Build the discrete diverging BAF colormap.

    Ten equal bins over [0, 1]: blue A-skewed, gray balanced, red B-skewed.
    Binning keeps small BAF differences from reading as real structure. NaN is
    white.

    Returns:
        Tuple of (cmap, norm) ready to pass to ``imshow`` or ``pcolormesh``.
    """
    cmap = mcolors.ListedColormap(BAF_COLORS, name="baf_disc")
    cmap.set_bad("white")
    norm = mcolors.BoundaryNorm(np.linspace(0, 1, 11), cmap.N, clip=True)
    return cmap, norm


# =============================================================================
# Categorical label palettes
# =============================================================================


def build_label_cmaps(
    row_label_map: dict,
    primary_label: str | None,
    invalid_labels=INVALID_LABELS,
    na_labels=NA_CELLTYPE,
) -> dict:
    """Build color maps for several label sets that share one palette.

    A value shared across sets gets the same color everywhere, so stacked
    annotation strips read against each other. The primary label is visited
    first in clone order; the rest follow alphabetically and reuse its colors.

    Args:
        row_label_map: {label_name: values} for each annotation to color.
        primary_label: Name of the label set that gets first pick, usually the
            clone assignment. None or an unknown name simply skips the priority.
        invalid_labels: Labels drawn in the missing-data gray.
        na_labels: Labels treated as missing.

    Returns:
        {label_name: {value: color}}, one entry per key of ``row_label_map``.
    """
    names = ([primary_label] if primary_label in row_label_map else []) + [
        n for n in row_label_map if n != primary_label
    ]

    def cats_for(name):
        """List one label set's distinct values in assignment order.

        Args:
            name: Key into ``row_label_map``.

        Returns:
            Distinct values, in clone order for the primary label and
            alphabetical otherwise.
        """
        uniq = {str(v) for v in row_label_map[name]}
        if name == primary_label:
            return sorted(uniq, key=lambda c: _clone_order_key(c, invalid_labels))
        return sorted(uniq)

    # first encounter in visit order fixes a value's color
    ordered_values = []
    seen = set()
    for name in names:
        for c in cats_for(name):
            if _is_colored_label(c, invalid_labels, na_labels) and c not in seen:
                seen.add(c)
                ordered_values.append(c)

    palette = (
        _TUMOR_COLORS if len(ordered_values) <= len(_TUMOR_COLORS) else _TUMOR_COLORS_20
    )
    value_color = {c: palette[i % len(palette)] for i, c in enumerate(ordered_values)}

    cmaps = {}
    for name in names:
        cmap = {}
        for c in cats_for(name):
            if c in invalid_labels or c in na_labels:
                cmap[c] = NA_COLOR
            elif _is_normal_like(c):
                cmap[c] = NORMAL_COLOR
            else:
                cmap[c] = value_color[c]
        cmaps[name] = cmap
    return cmaps


# =============================================================================
# Global style
# =============================================================================


def set_palette(num_colors: int = 8, style: str = "whitegrid") -> list:
    """Set the global seaborn style and categorical palette.

    Mutates seaborn global state, so call once per figure, not inside a drawing
    primitive.

    Args:
        num_colors: Distinct colors needed; above 8 switches Set2 to husl.
        style: Seaborn style name.

    Returns:
        The palette that was set, as a list of RGB tuples.
    """
    sns.set_style(style)
    if num_colors > 8:
        palette = sns.color_palette("husl", n_colors=num_colors)
    else:
        palette = sns.color_palette("Set2", n_colors=num_colors)
    sns.set_palette(palette)
    return palette
