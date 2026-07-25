"""Shared fixtures: a headless backend, the reference on disk, and the toy data.

The datasets are loaded from the materialized fixture files in ``tests/data``
(written by ``simulate.write_dataset``), so the tests read the same TSVs the
README and gallery do, and stay deterministic.
"""

import matplotlib
import pytest

matplotlib.use("Agg")

from cnplot import GenomeAxis
from simulate import load_dataset
from simulate import reference as _reference


@pytest.fixture(scope="session")
def reference():
    """(region_bed, chrom_sizes): the vendored T2T reference.

    Returns:
        Paths to the region BED and chromosome-sizes files.
    """
    return _reference()


@pytest.fixture
def axis(reference):
    """A shrunk-gap :class:`GenomeAxis` over the toy reference."""
    return GenomeAxis(*reference)


@pytest.fixture
def axis_keep(reference):
    """A kept-gap :class:`GenomeAxis`, for the gap-layout tests."""
    return GenomeAxis(*reference, collapse_gaps=False)


@pytest.fixture
def sim():
    """The single-sample dataset, sample id "S1"."""
    return load_dataset(["S1"])


@pytest.fixture
def sim2():
    """The two-sample toy dataset (shared profile, differing proportions)."""
    return load_dataset(["S1", "S2"])


@pytest.fixture
def state_palette(sim):
    """{state string: color} covering the "state" column, for hue coloring."""
    from cnplot import get_cn_cmap

    style, _ = get_cn_cmap()
    out = {}
    for s in sim.obs["state"].unique():
        a, b = (int(x) for x in s.split("|"))
        out[s] = style.get((a, b), style["default"])
    return out


@pytest.fixture
def saved(tmp_path):
    """Save a figure to the test's tmp dir and assert it rendered.

    Returns:
        A ``save(fig, name) -> Path`` callable asserting the file is non-empty.
    """
    import matplotlib.pyplot as plt

    def _save(fig, name):
        path = tmp_path / f"{name}.png"
        fig.savefig(path, dpi=60)
        plt.close(fig)
        assert path.exists() and path.stat().st_size > 0
        return path

    return _save
