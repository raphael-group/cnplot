"""Genome coordinate model: mapping bins and sites onto a plot axis.

A :class:`GenomeAxis` is a coordinate transform, not a parallel array. It is
constructed from the reference - a region BED and a chromosome-sizes file - and
then maps any table of bins onto the axis by ``(#CHR, START, END)``. One axis can
therefore serve several samples, several bin resolutions, and every panel of a
figure.

Pure geometry over numpy and pandas, with no matplotlib dependency; the drawing
helpers live in :mod:`cnplot.cnplot_utils`.

Stretches no region covers can be handled two ways:

- shrunk (``collapse_gaps=True``): gaps are removed, so the axis is the regions
  concatenated and each gap is zero-width.
- kept (``collapse_gaps=False``): chromosomes keep full length, so bins keep true
  spacing and each gap keeps its real width.

Both report gaps as :class:`Gap` intervals. The same bin lands at different
coordinates under each, so an axis must not be shared across the two.
"""

import logging
from collections import OrderedDict
from dataclasses import dataclass
from enum import Enum

import numpy as np
import pandas as pd

from .cnplot_io_utils import read_bed, read_chr_sizes

logger = logging.getLogger(__name__)

# Sex and mitochondrial chromosomes, excluded by default: their copy number is
# not comparable to the autosomes without ploidy assumptions the plots do not make.
DEFAULT_EXCLUDED_CHROMS = ("chrX", "chrY", "chrM")

__all__ = [
    "DEFAULT_EXCLUDED_CHROMS",
    "BinCoords",
    "Gap",
    "GapKind",
    "GenomeAxis",
    "Segment",
]


# =============================================================================
# Helpers
# =============================================================================


def _norm_chrom(name: str) -> str:
    """Normalize a chromosome name for comparison.

    Args:
        name: Chromosome name, with or without a "chr" prefix.

    Returns:
        Upper-cased name without the prefix, so "chrX", "chrx", and "X" match.
    """
    return str(name).upper().removeprefix("CHR")


def _bin_edges(coords_df: pd.DataFrame) -> tuple:
    """Extract bin edges, treating a point table as zero-length bins.

    Args:
        coords_df: Bins with "START"/"END", or sites with "POS".

    Returns:
        (starts, ends) float arrays.

    Raises:
        ValueError: If the required columns are missing.
    """
    if "#CHR" not in coords_df.columns:
        raise ValueError("coords_df must have a '#CHR' column")
    if {"START", "END"}.issubset(coords_df.columns):
        return (
            coords_df["START"].to_numpy(dtype=float),
            coords_df["END"].to_numpy(dtype=float),
        )
    if "POS" in coords_df.columns:
        pos = coords_df["POS"].to_numpy(dtype=float)
        return pos, pos
    raise ValueError("coords_df must have START/END or POS columns")


def _chrom_gaps(regions_ch: pd.DataFrame, size: float) -> list:
    """Derive a chromosome's uncovered stretches in raw genomic coordinates.

    Works off the regions themselves rather than their count, so a gap is
    reported only where consecutive regions fail to meet.

    Args:
        regions_ch: This chromosome's regions, sorted and non-overlapping.
        size: Chromosome length.

    Returns:
        List of (start, end, :class:`GapKind`) in raw genomic coordinates,
        ordered along the chromosome. Zero-length gaps are omitted.
    """
    out = []
    if regions_ch.iloc[0]["START"] > 0:
        out.append((0.0, float(regions_ch.iloc[0]["START"]), GapKind.LEADING))
    for si in range(len(regions_ch) - 1):
        this_end = regions_ch.iloc[si]["END"]
        next_start = regions_ch.iloc[si + 1]["START"]
        if next_start > this_end:
            out.append((float(this_end), float(next_start), GapKind.INTERIOR))
    last_end = regions_ch.iloc[-1]["END"]
    if last_end < size:
        out.append((float(last_end), float(size), GapKind.TRAILING))
    return out


# =============================================================================
# Coordinate model
# =============================================================================


class GapKind(str, Enum):
    """Where an uncovered stretch sits relative to a chromosome's regions.

    Subclasses ``str``, so ``kind == "interior"`` works alongside
    ``kind is GapKind.INTERIOR``.
    """

    LEADING = "leading"
    INTERIOR = "interior"
    TRAILING = "trailing"


@dataclass(frozen=True)
class Gap:
    """A stretch of a chromosome that no region covers.

    Shrinking zeroes the axis extent but not the genomic one, so a shrunk gap can
    still report the size of the region it stands for.

    Attributes:
        chrom: Chromosome the gap belongs to.
        start: Left edge in axis coordinates.
        end: Right edge in axis coordinates; equals ``start`` when shrunk.
        kind: Where the gap sits, as :class:`GapKind`.
        raw_start: Left edge in base pairs.
        raw_end: Right edge in base pairs.
    """

    chrom: str
    start: float
    end: float
    kind: GapKind
    raw_start: float
    raw_end: float

    @property
    def width(self) -> float:
        """Width in axis units, zero once shrunk."""
        return self.end - self.start

    @property
    def raw_width(self) -> float:
        """Length in base pairs, which survives shrinking."""
        return self.raw_end - self.raw_start

    @property
    def is_collapsed(self) -> bool:
        """Whether the gap was shrunk to a single position."""
        return self.end == self.start


@dataclass(frozen=True)
class Segment:
    """One drawn stretch of a chromosome, and where it starts on the axis.

    Segments tile the axis in order and carry the whole transform: a position
    ``p`` inside one maps to ``axis_start + (p - raw_start)``.

    Attributes:
        chrom: Chromosome the segment belongs to.
        raw_start: Left edge in base pairs.
        raw_end: Right edge in base pairs.
        axis_start: Left edge in axis coordinates.
    """

    chrom: str
    raw_start: float
    raw_end: float
    axis_start: float

    @property
    def axis_end(self) -> float:
        """Right edge in axis coordinates."""
        return self.axis_start + (self.raw_end - self.raw_start)


@dataclass(frozen=True)
class BinCoords:
    """Where a table of bins landed on a :class:`GenomeAxis`.

    Arrays follow the row order of the mapped table and hold NaN for bins that
    no segment covers.

    Attributes:
        positions: (n,) bin midpoints in axis coordinates.
        starts: (n,) bin left edges in axis coordinates.
        ends: (n,) bin right edges in axis coordinates.
    """

    positions: np.ndarray
    starts: np.ndarray
    ends: np.ndarray

    @property
    def mapped(self) -> np.ndarray:
        """(n,) bool mask of the bins that landed on the axis."""
        return np.isfinite(self.starts)

    @property
    def n_mapped(self) -> int:
        """How many bins landed on the axis."""
        return int(self.mapped.sum())

    @property
    def n_unmapped(self) -> int:
        """How many bins no segment covered, and so are not drawn."""
        return int((~self.mapped).sum())


class GenomeAxis:
    """A genome-wide coordinate transform built from the reference.

    The axis depends only on the reference, not on any bin table, so construct it
    once and map every sample and panel through it with
    :meth:`build_coordinates`.

    Attributes:
        segments: :class:`Segment` stretches in axis order, carrying the
            transform.
        ch_coords: len(chrs)+1 chromosome start offsets, last entry the genome
            end, so chromosome i spans ch_coords[i] to ch_coords[i + 1].
        gaps: :class:`Gap` intervals, ordered along the genome.
        axis_start: Left axis limit, including the ``chr_shift`` pad.
        axis_end: Right axis limit, including the ``chr_shift`` pad.
        collapsed: Whether gaps were removed from the axis.
        chrom_sizes: Chromosome lengths read from ``chrom_sizes``, after
            exclusions.
        regions: Regions read from ``region_bed``, after exclusions.
        excluded_chroms: Chromosomes that were dropped.
    """

    def __init__(
        self,
        region_bed: str,
        chrom_sizes: str,
        excluded_chroms: list | tuple | None = DEFAULT_EXCLUDED_CHROMS,
        chr_shift: float = 0.0,
        collapse_gaps: bool = True,
    ):
        """Build the transform from a region BED and a chromosome-sizes file.

        Chromosome order always follows the sizes file, so it does not change
        with ``collapse_gaps`` and is controlled by ordering that file. A
        chromosome with no regions has no extent when gaps are shrunk and is
        dropped; when they are kept it stays, entirely a gap.

        Args:
            region_bed: Path to the regions to draw, sorted and non-overlapping
                within each chromosome.
            chrom_sizes: Path to the chromosome-sizes file. Required either way:
                it is what makes a trailing gap detectable, since regions alone
                cannot say where a chromosome ends.
            excluded_chroms: Chromosomes to drop from both inputs before
                building, matched ignoring case and any "chr" prefix. Defaults to
                :data:`DEFAULT_EXCLUDED_CHROMS`; pass an empty list to keep
                everything.
            chr_shift: Blank pad before the first and after the last chromosome.
            collapse_gaps: Remove uncovered stretches from the axis. False keeps
                true genomic spacing and each gap's real width.

        Raises:
            ValueError: If a region ends past its chromosome's length, or names
                a chromosome absent from the sizes file.
        """
        sizes = read_chr_sizes(chrom_sizes)
        regions = read_bed(region_bed)
        drop = {_norm_chrom(c) for c in (excluded_chroms or ())}
        if drop:
            removed = [ch for ch in sizes if _norm_chrom(ch) in drop]
            if removed:
                logger.info("excluding %s from the axis", removed)
            sizes = OrderedDict(
                (ch, sz) for ch, sz in sizes.items() if _norm_chrom(ch) not in drop
            )
            regions = regions[
                ~regions["#CHR"].map(lambda c: _norm_chrom(c) in drop)
            ].reset_index(drop=True)
        self.excluded_chroms = tuple(excluded_chroms or ())
        self.chrom_sizes = sizes
        self.regions = regions
        self.collapsed = bool(collapse_gaps)
        self.axis_start = float(chr_shift)

        segments: list = []
        ch_coords: list = []
        gaps: list = []
        offset = float(chr_shift)

        regions_chs = self.regions.groupby("#CHR", sort=False)
        unknown = [ch for ch in regions_chs.groups if ch not in self.chrom_sizes]
        if unknown:
            raise ValueError(f"{unknown} in the region BED but not the sizes file")

        for ch, size in self.chrom_sizes.items():
            has_regions = ch in regions_chs.groups
            if self.collapsed and not has_regions:
                logger.debug("no region for %s; dropping it from the axis", ch)
                continue

            ch_coords.append(offset)
            if not has_regions:
                logger.warning("no region for %s; all of it is a gap", ch)
                segments.append(Segment(ch, 0.0, float(size), offset))
                gaps.append(
                    Gap(ch, offset, offset + size, GapKind.LEADING, 0.0, float(size))
                )
                offset += size
                continue

            regions_ch = regions_chs.get_group(ch)
            if self.collapsed:
                # position 0 is the chromosome start; a region's END is where it
                # ended, so a gap starting there lands on the running offset
                anchor = {0.0: offset}
                for _, row in regions_ch.iterrows():
                    raw_s, raw_e = float(row["START"]), float(row["END"])
                    if raw_e > size:
                        raise ValueError(
                            f"region {ch}:{raw_s:.0f}-{raw_e:.0f} runs past the "
                            f"chromosome length {size}"
                        )
                    segments.append(Segment(ch, raw_s, raw_e, offset))
                    offset += raw_e - raw_s
                    anchor[raw_e] = offset
                for gs, ge, kind in _chrom_gaps(regions_ch, size):
                    at = anchor[gs]
                    gaps.append(Gap(ch, at, at, kind, gs, ge))
            else:
                segments.append(Segment(ch, 0.0, float(size), offset))
                for gs, ge, kind in _chrom_gaps(regions_ch, size):
                    gaps.append(Gap(ch, offset + gs, offset + ge, kind, gs, ge))
                offset += size

        ch_coords.append(offset)
        self.segments = segments
        self.ch_coords = ch_coords
        self.gaps = gaps
        self.axis_end = offset + float(chr_shift)

    def __repr__(self) -> str:
        """Summarize the transform.

        Returns:
            A string with the chromosome, segment, and gap counts.
        """
        return (
            f"GenomeAxis({len(self.chrs)} chrs, {len(self.segments)} segments, "
            f"{len(self.gaps)} gaps, "
            f"{'shrunk' if self.collapsed else 'kept'} gaps, "
            f"span={self.chr_end - self.ch_coords[0]:.0f})"
        )

    @property
    def chrs(self) -> list:
        """Chromosome names in plotting order."""
        return list(dict.fromkeys(s.chrom for s in self.segments))

    @property
    def chr_end(self) -> float:
        """Right edge of the last chromosome, excluding the trailing pad."""
        return self.ch_coords[-1]

    @property
    def chr_boundaries(self) -> list:
        """Internal chromosome boundaries, excluding both genome ends."""
        return self.ch_coords[1:-1]

    @property
    def xtick_chrs(self) -> list:
        """Chromosome midpoints, for placing one tick per chromosome."""
        return [
            (self.ch_coords[i] + self.ch_coords[i + 1]) / 2
            for i in range(len(self.chrs))
        ]

    @property
    def xlab_chrs(self) -> list:
        """Chromosome tick labels, aligned to :attr:`xtick_chrs`."""
        return self.chrs

    @property
    def chr_offsets(self) -> dict:
        """Each chromosome's start offset on the axis."""
        return {ch: self.ch_coords[i] for i, ch in enumerate(self.chrs)}

    @property
    def seg_coords(self) -> list:
        """Left edge of every gap, which is its position once shrunk."""
        return [g.start for g in self.gaps]

    @property
    def interior_gaps(self) -> list:
        """Gaps falling between two regions of one chromosome."""
        return [g for g in self.gaps if g.kind is GapKind.INTERIOR]

    @property
    def edge_gaps(self) -> list:
        """Gaps before a chromosome's first region or after its last."""
        return [g for g in self.gaps if g.kind is not GapKind.INTERIOR]

    def build_coordinates(
        self, coords_df: pd.DataFrame, contain: bool = False
    ) -> BinCoords:
        """Place a table of bins or sites on this axis.

        Bins are matched by coordinate, not by row position, so the table needs
        no relationship to whatever the axis was built from. Bins covered by no
        segment map to NaN and are not drawn; a warning reports how many, split
        into those on chromosomes the axis does not carry and those falling in a
        gap.

        Args:
            coords_df: Bins with "#CHR", "START", "END", or sites with "#CHR",
                "POS". Row order defines the output arrays.
            contain: If True, a bin joins a segment only when fully inside it,
                as integer-CN profiles need. If False, overlapping bins join and
                are clipped to the segment.

        Returns:
            A :class:`BinCoords` aligned to ``coords_df`` row for row.

        Raises:
            ValueError: If ``coords_df`` lacks the required columns.
        """
        starts, ends = _bin_edges(coords_df)
        n = len(coords_df)
        out_start = np.full(n, np.nan)
        out_end = np.full(n, np.nan)
        chrom = coords_df["#CHR"].to_numpy()

        for seg in self.segments:
            on_chrom = chrom == seg.chrom
            if not on_chrom.any():
                continue
            if contain:
                hit = on_chrom & (starts >= seg.raw_start) & (ends <= seg.raw_end)
            else:
                hit = on_chrom & (starts < seg.raw_end) & (ends > seg.raw_start)
            if not hit.any():
                continue
            shift = seg.axis_start - seg.raw_start
            out_start[hit] = np.maximum(starts[hit], seg.raw_start) + shift
            out_end[hit] = np.minimum(ends[hit], seg.raw_end) + shift

        unmapped = ~np.isfinite(out_start)
        if unmapped.any():
            off_axis = int((unmapped & ~np.isin(chrom, self.chrs)).sum())
            in_gap = int(unmapped.sum()) - off_axis
            logger.warning(
                "%d/%d bins are not drawn: %d on chromosomes outside the axis, "
                "%d inside gaps",
                int(unmapped.sum()),
                n,
                off_axis,
                in_gap,
            )
        return BinCoords(
            positions=(out_start + out_end) / 2, starts=out_start, ends=out_end
        )

    def grid(self, coords_df: pd.DataFrame, contain: bool = False) -> tuple:
        """Build a ``pcolormesh`` grid covering the axis.

        Cells run left to right without overlap: one per mapped bin, plus filler
        cells wherever bins do not tile the axis.

        Args:
            coords_df: Bins to place, as for :meth:`build_coordinates`.
            contain: Bin membership rule, as for :meth:`build_coordinates`.

        Returns:
            (x_edges, col_bin_ids): the (m+1,) cell edges, and the ``coords_df``
            row index for each of the m cells, -1 for filler.
        """
        coords = self.build_coordinates(coords_df, contain=contain)
        order = np.argsort(coords.starts, kind="stable")
        order = order[np.isfinite(coords.starts[order])]

        x_edges = [self.axis_start]
        col_bin_ids = []
        cur = self.axis_start
        for row in order:
            s, e = float(coords.starts[row]), float(coords.ends[row])
            if s > cur:
                col_bin_ids.append(-1)
                x_edges.append(s)
                cur = s
            if e > cur:
                col_bin_ids.append(int(row))
                x_edges.append(e)
                cur = e
        if cur < self.chr_end:
            col_bin_ids.append(-1)
            x_edges.append(self.chr_end)
        return np.asarray(x_edges, dtype=float), col_bin_ids
