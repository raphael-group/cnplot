"""Typography for every cnplot figure: the one place fonts are configured.

Nothing else in the package names a font family, a font size or a font
embedding type. The drawing modules import the ``FS_*`` sizes from here, and
:func:`apply_style` writes the same values into matplotlib's rcParams when the
package is imported, so text drawn without an explicit size - a colorbar's
ticks, a caller's own axes on the same figure - follows too. To restyle the
package, edit the constants below.

The values follow manuscript figure rules: text in Arial or Helvetica, nothing
under 18 pt because a panel is drawn wide and reproduced a few inches wide,
and fonts embedded so text stays text in a vector editor.

``FONT_FAMILY`` is a preference list and the first installed family is used;
DejaVu Sans closes it so a machine without the others falls back instead of
raising. ``FONT_TYPE`` replaces matplotlib's Type 3 default, which a vector
editor imports as outlines. ``FS_ANNOT`` is the one size under the floor, for
labels written on the data, which share the plotting area with the marks.
"""

import contextlib

import matplotlib.font_manager
import matplotlib.pyplot as plt

__all__ = [
    "FONT_FAMILY",
    "FONT_TYPE",
    "FS_ANNOT",
    "FS_AXLABEL",
    "FS_AXTITLE",
    "FS_LEGEND",
    "FS_TICK",
    "FS_TITLE",
    "apply_style",
]

# NB: preference order, first installed wins; DejaVu Sans ships with matplotlib
FONT_FAMILY = ("Arial", "Helvetica", "DejaVu Sans")

# NB: Type 42 embeds TrueType; "none" keeps SVG text as text
FONT_TYPE = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}

# NB: points
FS_TICK = 18
FS_AXLABEL = 22
FS_AXTITLE = 18
FS_TITLE = 24
FS_LEGEND = 18
# NB: labels written on the data
FS_ANNOT = 15

_RC_SIZES = {
    "font.size": FS_TICK,
    "axes.labelsize": FS_AXLABEL,
    "axes.titlesize": FS_AXTITLE,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "legend.fontsize": FS_LEGEND,
    "legend.title_fontsize": FS_LEGEND,
    "figure.titlesize": FS_TITLE,
}


def apply_style() -> None:
    """Write the family, sizes and embedding types above into rcParams.

    Mathtext is set in the same face as the text, so an exponent in an axis
    label is not in a different face from the words beside it; that needs the
    name of a family that is actually installed, which is resolved here.

    Called once when cnplot is imported, because a figure drawn by cnplot may be
    saved by a caller that never touches rcParams. Call it again after anything
    that resets rcParams.
    """
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = list(FONT_FAMILY)
    plt.rcParams.update(_RC_SIZES)
    plt.rcParams.update(FONT_TYPE)
    with contextlib.suppress(Exception):
        installed = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
        family = next((f for f in FONT_FAMILY if f in installed), FONT_FAMILY[-1])
        plt.rcParams["mathtext.fontset"] = "custom"
        plt.rcParams["mathtext.rm"] = family
        plt.rcParams["mathtext.it"] = f"{family}:italic"
        plt.rcParams["mathtext.bf"] = f"{family}:bold"
        plt.rcParams["mathtext.default"] = "regular"


apply_style()
