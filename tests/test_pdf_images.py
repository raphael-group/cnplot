"""Rasterized layers are written as 8-bit DeviceRGB, whichever writer is used."""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.backends import backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap

from cnplot import FigureSaver
from cnplot.cnplot_utils import _write_truecolor_images

# NB: 2, 3 and 12 colours pack at 1, 2 and 4 bits in an unpatched backend
N_COLORS = [2, 3, 12]


def _raster(n_colors):
    fig, ax = plt.subplots(figsize=(4, 2))
    values = np.arange(40 * 200).reshape(40, 200) % n_colors
    cmap = ListedColormap(plt.cm.tab20.colors[:n_colors])
    ax.imshow(values, cmap=cmap, aspect="auto", interpolation="nearest")
    ax.images[0].set_rasterized(True)
    return fig


def _write(writer, out, n_colors):
    """Write one raster page to ``{out}.pdf`` through the named writer."""
    fig = _raster(n_colors)
    path = f"{out}.pdf"
    if writer == "savefig":
        fig.savefig(path, dpi=100)
    elif writer == "pdfpages":
        with PdfPages(path) as pdf:
            pdf.savefig(fig, dpi=100)
    else:
        with FigureSaver(str(out), "pdf", dpi=100) as saver:
            saver.savefig(fig)
    plt.close(fig)
    return path


def _image_encodings(path):
    """(has a palette, bits per component of every image) read from the PDF bytes."""
    raw = Path(path).read_bytes()
    assert b"/Subtype /Image" in raw
    depths = [int(d) for d in re.findall(rb"/BitsPerComponent (\d+)", raw)]
    return b"/Indexed" in raw, depths


@pytest.mark.parametrize("n_colors", N_COLORS)
@pytest.mark.parametrize("writer", ["savefig", "pdfpages", "figuresaver"])
def test_few_colour_raster_is_8bit_rgb(tmp_path, writer, n_colors):
    indexed, depths = _image_encodings(_write(writer, tmp_path / "page", n_colors))
    assert not indexed
    assert depths and set(depths) == {8}


def test_hook_is_installed_once_and_leaves_pillow_in_place(tmp_path):
    assert _write_truecolor_images() is False
    _write("savefig", tmp_path / "page", 3)
    assert backend_pdf.Image.__name__ == "PIL.Image"


def test_hook_warns_when_the_backend_bypasses_it(monkeypatch, caplog):
    monkeypatch.setattr(backend_pdf.PdfFile, "_writeImg", lambda self: None)
    assert _write_truecolor_images() is True
    with caplog.at_level("WARNING", logger="cnplot.cnplot_utils"):
        backend_pdf.PdfFile._writeImg(None)
        backend_pdf.PdfFile._writeImg(None)
    assert [r.message for r in caplog.records].count(caplog.records[0].message) == 1
    assert "palette" in caplog.records[0].message


def test_stock_backend_writes_a_palette(tmp_path, monkeypatch):
    monkeypatch.setattr(
        backend_pdf.PdfFile, "_writeImg", backend_pdf.PdfFile._writeImg.__wrapped__
    )
    indexed, _ = _image_encodings(_write("savefig", tmp_path / "stock", 12))
    assert indexed


def test_pixels_match_the_uncompressed_raster(tmp_path):
    fitz = pytest.importorskip("fitz")

    def render(path):
        return fitz.open(path)[0].get_pixmap(dpi=100).samples

    patched = render(_write("savefig", tmp_path / "patched", 12))
    with plt.rc_context({"pdf.compression": 0}):
        exact = render(_write("savefig", tmp_path / "exact", 12))
    assert patched == exact
