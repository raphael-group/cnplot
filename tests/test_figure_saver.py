"""FigureSaver output paths and the rasterized-image memory workaround."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cnplot import FigureSaver


def _page(rasterized, figsize=(4, 8)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(np.linspace(0, 1, 50), np.linspace(0, 1, 50), rasterized=rasterized)
    return fig


def _retained_images(saver):
    """The (image, name, obj) entries the open PDF is still holding."""
    return list(saver._pdf._file._images.values())


def _pinned_bytes(image):
    """Bytes kept alive by ``image``, following it back to the buffer it views."""
    while getattr(image, "base", None) is not None:
        image = image.base
    return image.nbytes


def test_rasterized_pages_do_not_pin_full_figure_buffers(tmp_path):
    out = tmp_path / "pages"
    with FigureSaver(str(out), "pdf") as saver:
        for _ in range(3):
            fig = _page(rasterized=True)
            saver.savefig(fig)
            plt.close(fig)
        retained = _retained_images(saver)
        assert retained
        # matplotlib hands the backend a view into a whole-figure raster buffer
        # and holds every image until close; unpatched, each page pins a buffer
        # far larger than the pixels it contributes
        pinned = sum(_pinned_bytes(img) for img, _, _ in retained)
        cropped = sum(img.nbytes for img, _, _ in retained)
        assert pinned == cropped
    assert out.with_suffix(".pdf").exists()


def test_rasterized_pixels_survive_the_copy(tmp_path):
    fitz = pytest.importorskip("fitz")
    rendered = {}
    for name, rasterized in (("raster", True), ("vector", False)):
        out = tmp_path / name
        with FigureSaver(str(out), "pdf", dpi=100) as saver:
            fig = _page(rasterized=rasterized, figsize=(4, 4))
            saver.savefig(fig)
            plt.close(fig)
        page = fitz.open(out.with_suffix(".pdf"))[0].get_pixmap(dpi=100)
        rendered[name] = np.frombuffer(page.samples, np.uint8).reshape(
            page.height, page.width, page.n
        )
    assert rendered["raster"].shape == rendered["vector"].shape
    diff = np.abs(rendered["raster"].astype(int) - rendered["vector"].astype(int))
    assert diff.mean() < 1.0


def test_close_is_idempotent_and_writes_side_formats(tmp_path):
    out = tmp_path / "pages"
    saver = FigureSaver(str(out), ["pdf", "png"], dpi=50)
    fig = _page(rasterized=True, figsize=(2, 2))
    saver.savefig(fig)
    plt.close(fig)
    saver.close()
    saver.close()
    assert out.with_suffix(".pdf").exists()
    assert (tmp_path / "pages.p0.png").exists()
