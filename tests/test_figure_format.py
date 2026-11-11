"""Figures are exported as PDF (vector) unless another format is requested.

A PNG is a pixel grid: text and hairlines degrade as soon as a reader zooms, and a
journal that wants vector art needs a re-export. PDF keeps the drawing commands, so
one file serves on-screen reading and print. These tests pin that default, and pin
that asking for a raster format still gets exactly that format -- the change must not
make PNG unreachable.
"""
from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from easyepi._plots import (DEFAULT_FIGURE_FORMAT, FIGURE_FORMATS,  # noqa: E402
                            save_figure)
from easyepi.pubfig import finalize_figure  # noqa: E402


def _figure():
    fig, ax = plt.subplots(figsize=(4, 2))
    ax.plot([0, 1, 2], [0, 1, 4])
    ax.set_title("t")
    return fig


def test_default_format_is_pdf():
    assert DEFAULT_FIGURE_FORMAT == "pdf"
    # PDF is offered first, so `--help` leads with the vector option
    assert FIGURE_FORMATS[0] == "pdf"
    assert "png" in FIGURE_FORMATS


def test_save_figure_writes_pdf_by_default(tmp_path):
    written = save_figure(_figure(), tmp_path, "default")
    assert [p.rsplit("/", 1)[-1].rsplit("\\", 1)[-1] for p in written] == ["default.pdf"]
    from pathlib import Path
    assert Path(written[0]).read_bytes().startswith(b"%PDF")


def test_the_pdf_carries_no_raster_image(tmp_path):
    """A vector PDF must not embed a bitmap; matplotlib emits no image XObject here."""
    from pathlib import Path

    written = save_figure(_figure(), tmp_path, "vector")
    data = Path(written[0]).read_bytes()
    assert b"/Subtype /Image" not in data
    assert b"/Subtype/Image" not in data


def test_an_explicit_raster_format_is_still_honoured(tmp_path):
    from pathlib import Path

    written = save_figure(_figure(), tmp_path, "raster", formats=["png"])
    assert written[0].endswith(".png")
    assert Path(written[0]).read_bytes().startswith(b"\x89PNG")


def test_several_formats_in_one_call(tmp_path):
    written = save_figure(_figure(), tmp_path, "both", formats=["pdf", "png"])
    assert [w[-4:] for w in written] == [".pdf", ".png"]


def test_finalize_figure_without_extension_defaults_to_pdf(tmp_path):
    written = finalize_figure(_figure(), str(tmp_path / "out" / "figure"))
    assert written[0].endswith(".pdf")
