"""Readers for the reference files the plots are built against.

Kept separate from the plotting modules so both the coordinate model and the
drawing helpers can read input without importing each other. Pandas only, with
no matplotlib dependency.
"""

import logging
from collections import OrderedDict

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "read_chr_sizes",
    "read_bed",
]


def read_chr_sizes(sz_file: str) -> OrderedDict:
    """Read a two-column chromosome-sizes file.

    Args:
        sz_file: Path to a whitespace-separated (chromosome, length) file.

    Returns:
        Chromosome name to length, in file order.
    """
    chr_sizes = OrderedDict()
    with open(sz_file) as rfd:
        for line in rfd:
            if not line.strip():
                continue
            ch, size = line.split()[:2]
            chr_sizes[ch] = int(size)
    return chr_sizes


def read_bed(bed_file: str, names: list | None = None, **kwargs) -> pd.DataFrame:
    """Read a headerless region BED into a table.

    The one BED representation the package uses: the genome axis is built from
    it, and :func:`cnplot.cnplot_utils.shade_regions` shades from it.

    Args:
        bed_file: Path to a tab-separated BED file.
        names: Column names, defaulting to "#CHR", "START", "END", "NAME".
        **kwargs: Passed through to ``pandas.read_table``.

    Returns:
        The regions, in file order.
    """
    if names is None:
        names = ["#CHR", "START", "END", "NAME"]
    return pd.read_table(bed_file, sep="\t", header=None, names=names, **kwargs)
