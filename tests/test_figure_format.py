"""Figures are exported as PDF (vector) unless another format is requested.

A PNG is a pixel grid: text and hairlines degrade as soon as a reader zooms, and a
journal that wants vector art needs a re-export. PDF keeps the drawing commands, so
one file serves on-screen reading and print. These tests pin that default, and pin
that asking for a raster format still gets exactly that format -- the change must not
make PNG unreachable.
"""
from __future__ import annotations

import pytest

from pathlib import Path

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


# --------------------------------------------------- command-line figure defaults
def test_no_flag_writes_a_pdf_and_no_figures_skips_it(tmp_path):
    """The CLI default exports; `--no-figures` is the only way to skip writing.

    Driven through `render_module_figures` with a module that needs no optional
    dependency (`motif_enrichment` needs biopython, which the CI test job does not
    install). What is under test is the figure plumbing, not any one analysis: the
    bug this pins was that the flag plumbing was fine while the render path raised
    inside a best-effort wrapper, so a run looked successful and wrote nothing.
    """
    import json
    from easyepi._plots import render_module_figures

    def run(formats, name):
        work = tmp_path / name
        work.mkdir()
        # A qc_metrics table is enough: the QC branch renders a bar chart from it.
        (work / "qc_metrics.tsv").write_text(
            "sample_id\tmetric\nS1\t0.5\nS2\t0.7\n", encoding="utf-8")
        made = render_module_figures("qc_metrics", work, formats)
        return sorted(Path(m).name for m in made)

    assert run(("pdf",), "default") == ["qc_metrics.pdf"]
    assert run((), "off") == []
    assert run(("png",), "png") == ["qc_metrics.png"]
    # and the file really is what the extension says
    assert (tmp_path / "default" / "figures" / "qc_metrics.pdf").read_bytes()[:4] == b"%PDF"


def test_matrix_input_does_not_break_a_bar_figure(tmp_path):
    """`render_module_figures` passes some tables as arrays, not DataFrames.

    `_bar`/`_distribution` used to reach `DataFrame.to_numpy` on an ndarray and raise
    AttributeError, which the best-effort wrapper logged and swallowed -- so the
    figure silently never appeared.
    """
    import numpy as np
    from easyepi._plots import _bar, _distribution

    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    assert _bar(arr, tmp_path, "bar", formats=("pdf",)) == [str(tmp_path / "bar.pdf")]
    assert _distribution(arr, tmp_path, "dist", formats=("pdf",)) \
        == [str(tmp_path / "dist.pdf")]


def test_motif_enrichment_json_renders_with_numpy_input(tmp_path):
    """`np_log10` returns an ndarray; the later `.iloc` calls need a Series."""
    import json
    from easyepi._plots import _motif_json

    payload = {"motifs": [
        {"motif": "A", "padj": 1e-5, "odds_ratio": 2.0},
        {"motif": "B", "padj": 0.2, "odds_ratio": 1.1},
    ]}
    path = tmp_path / "motif_enrichment.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    made = _motif_json(path, tmp_path / "fig", formats=("pdf",))
    assert made and made[0].endswith("motif_enrichment.pdf")
