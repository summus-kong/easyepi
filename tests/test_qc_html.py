"""The HTML QC report: self-contained, and honest about absent values.

Two properties are load-bearing and easy to break silently:

* **No external references.** The report is opened from an email attachment or a
  network share, so one stray ``src="figures/x.png"`` makes it render blank where
  a reviewer will assume the analysis failed. Every test asserts the file's only
  resources are inline ``data:`` URIs.
* **"Not measured" never becomes zero.** A reader skimming a table cannot tell a
  real ``0`` from an absent metric once both look like ``0``, and the whole point
  of the threshold layer is that those mean different things.
"""

from __future__ import annotations

import re
from html.parser import HTMLParser

from easyepi.qc_html import fmt_value, render_merged_html, render_qc_html

#: 1x1 PNG, so a figure can be inlined without needing matplotlib in this test.
_TINY_PNG = bytes.fromhex(
    "89504e470d0a1a0a0000000d49484452000000010000000108060000001f15c4"
    "890000000a49444154789c6300010000050001" "0d0a2db40000000049454e44ae426082")


class _Check(HTMLParser):
    """Collect external references and unclosed-tag hints from a document."""

    def __init__(self):
        super().__init__()
        self.external = []
        self.tags = 0

    def handle_starttag(self, tag, attrs):
        self.tags += 1
        for name, value in attrs:
            if not value:
                continue
            if name in ("src", "href", "data", "action", "srcset"):
                if value.startswith(("http://", "https://", "//")):
                    self.external.append((tag, name, value))
                elif name != "href" and not value.startswith(
                        ("data:", "#", "mailto:")):
                    self.external.append((tag, name, value))


def _audit(path):
    parser = _Check()
    text = path.read_text(encoding="utf-8")
    parser.feed(text)
    return parser, text


def _payload(**over):
    """A payload shaped like the CLI's JSON, including the awkward cases."""
    base = {
        "assay": "ATAC-seq",
        "verdict_summary": {"samples": {"pass": 1, "warn": 1,
                                        "not_evaluated": 0},
                            "warning_metrics": {"frip": 1}},
        "samples": {
            "good": {
                "sample_id": "good", "assay": "ATAC-seq",
                "bam_qc": {"mapping_rate": 0.9871, "n_reads": 12_000_000},
                "assay_qc": {"frip": 0.4231, "tss_enrichment": 12.5,
                             "median_fragment_length": 165,
                             "tss_enrichment_status": "measured"},
                "qc_verdict": {"overall": "pass", "warnings": [],
                               "core_coverage": 1.0,
                               "checks": [{"metric": "frip", "status": "pass"}]},
            },
            "sparse": {
                "sample_id": "sparse", "assay": "ATAC-seq",
                "bam_qc": {"mapping_rate": 0.41},
                "assay_qc": {"frip": 0.02,
                             # deliberately absent, not zero
                             "tss_enrichment": None,
                             "tss_enrichment_status": "background_zero"},
                "qc_verdict": {"overall": "warn", "warnings": ["frip"],
                               "core_coverage": 0.75,
                               "checks": [{"metric": "frip",
                                           "status": "breached"}]},
            },
        },
    }
    base.update(over)
    return base


# ---------------------------------------------------------------- basic shape
def test_report_is_written_and_parses(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    assert out.exists()
    parser, text = _audit(out)
    assert parser.tags > 10
    assert text.startswith("<!DOCTYPE html>")
    assert "good" in text and "sparse" in text


def test_report_has_no_external_references(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    parser, _ = _audit(out)
    assert parser.external == [], parser.external


def test_figures_are_inlined_as_data_uris(tmp_path):
    fig_dir = tmp_path / "figures"
    fig_dir.mkdir()
    (fig_dir / "frip.png").write_bytes(_TINY_PNG)
    out = render_qc_html(_payload(), tmp_path / "qc.html", figure_dir=tmp_path)
    parser, text = _audit(out)
    assert parser.external == []
    assert "data:image/png;base64," in text
    assert 'alt="frip"' in text


def test_pdf_only_figures_do_not_break_the_page(tmp_path):
    """A PDF cannot be shown in an <img>; it must be skipped, not emitted broken."""
    fig_dir = tmp_path / "figures"
    fig_dir.mkdir()
    (fig_dir / "frip.pdf").write_bytes(b"%PDF-1.4 fake")
    out = render_qc_html(_payload(), tmp_path / "qc.html", figure_dir=tmp_path)
    _, text = _audit(out)
    assert "<img" not in text
    assert "data:application/pdf" not in text


def test_missing_figure_dir_is_not_an_error(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html",
                         figure_dir=tmp_path / "nope")
    assert out.exists()


# --------------------------------------------------- absent values stay absent
def test_absent_metric_is_a_dash_not_a_zero(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    # tss_enrichment is None for "sparse"; it must not be printed as 0. Checked in
    # the metric table specifically, because the verdict tiles legitimately print
    # a count of 0 samples.
    # Every metric table in the document, not just the first one: the sample that
    # has the absent value is not necessarily the first card.
    tables = re.findall(r'<table class="sortable">.*?</table>', text, re.S)
    assert tables
    assert any("\u2014" in t for t in tables), "no dash rendered anywhere"
    # And the row exists rather than being silently dropped.
    assert any("TSS enrichment" in t for t in tables)


def test_fmt_value_reports_none_and_nan_as_not_measured():
    for value in (None, float("nan")):
        text, css = fmt_value(value)
        assert text == "\u2014"
        assert css == "na"
    # A real zero is a number, and must stay distinguishable from the above.
    text, css = fmt_value(0.0)
    assert text == "0"
    assert css == "num"


def test_status_strings_are_shown_verbatim(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    # background_zero is the whole point of the four-state TSS reporting: it must
    # reach the report, not be reduced to a missing value.
    assert "background_zero" in text


def test_two_tss_definitions_are_reported_as_disagreeing(tmp_path):
    payload = _payload()
    payload["samples"]["good"]["assay_qc"]["tss_enrichment_bigwig"] = 42.4
    payload["samples"]["good"]["assay_qc"]["tss_enrichment_disagreement"] = {
        "bam": 1220.08, "bigwig": 42.4, "ratio": 28.8,
        "note": "different windows",
    }
    out = render_qc_html(payload, tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    assert "disagree" in text.lower()
    assert "28.8" in text


def test_synthesised_zero_is_shown_as_absent(tmp_path):
    """`tss_enrichment = 0.0` with status background_zero is not a score of zero.

    `assay_metrics` writes 0.0 for a TSS ratio it could not form, so that a floor
    threshold still fires. Rendering that 0 in a report says the library scored
    zero, which is a different claim from "the ratio was undefined" -- and the
    status is right there in the payload saying which it is.
    """
    payload = _payload()
    payload["samples"]["sparse"]["assay_qc"].update(
        {"tss_enrichment": 0.0, "tss_enrichment_status": "background_zero"})
    out = render_qc_html(payload, tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    # One card per sample; the asserting one is the sample whose payload carries
    # the sentinel, which is not necessarily the first card in the document.
    cards = re.findall(r"<section class=\"card\">(.*?)</section>", text, re.S)
    sparse_card = next(c for c in cards if "background_zero" in c)
    row = re.search(r'<td class="metric">TSS enrichment</td><td class="(\w+)"[^>]*>([^<]*)<',
                    sparse_card)
    assert row, sparse_card
    assert row.group(1) == "na", row.group(0)
    assert row.group(2) == "\u2014"


def test_a_real_zero_is_still_shown_as_zero(tmp_path):
    """The guard above must not swallow genuine zeros, e.g. mitochondrial fraction."""
    payload = _payload()
    payload["samples"]["good"]["bam_qc"]["mitochondrial_fraction"] = 0.0
    out = render_qc_html(payload, tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    assert re.search(r'mitochondrial fraction</td><td class="num"[^>]*>0<', text)


def test_empty_metric_table_says_why(tmp_path):
    payload = _payload()
    payload["samples"] = {"noassay": {"sample_id": "noassay", "assay_qc": {},
                                      "qc_verdict": {"overall": "not_evaluated"}}}
    out = render_qc_html(payload, tmp_path / "qc.html")
    # Case-insensitive: the sentence is a sentence ("No assay metrics were
    # computed...") and the assertion is about the message existing, not its
    # capitalisation.
    assert "no assay metrics" in out.read_text(encoding="utf-8").lower()


# ------------------------------------------------------------- verdict fidelity
def test_verdict_badges_match_the_json(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    assert '<span class="badge ok">pass</span>' in text
    assert '<span class="badge warn">warn</span>' in text


def test_no_invented_pass_or_fail(tmp_path):
    """The report must not add a vocabulary the threshold layer does not use."""
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    lowered = out.read_text(encoding="utf-8").lower()
    for invented in (">fail<", ">failed<", ">error<"):
        assert invented not in lowered


def test_breached_rows_are_marked(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    assert 'class="breach"' in text


def test_warning_names_are_listed(tmp_path):
    out = render_qc_html(_payload(), tmp_path / "qc.html")
    assert "frip" in out.read_text(encoding="utf-8")


# --------------------------------------------------------------- escaping
def test_html_in_sample_ids_is_escaped(tmp_path):
    payload = _payload()
    payload["samples"] = {
        '<script>alert("x")</script>': payload["samples"]["good"]}
    out = render_qc_html(payload, tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    assert "<script>alert" not in text
    assert "&lt;script&gt;" in text


def test_html_in_notes_is_escaped(tmp_path):
    payload = _payload()
    payload["samples"]["sparse"]["qc_verdict"]["warnings"] = ['<b>frip</b>']
    out = render_qc_html(payload, tmp_path / "qc.html")
    text = out.read_text(encoding="utf-8")
    assert "<b>frip</b>" not in text
    assert "&lt;b&gt;frip" in text


# ------------------------------------------------------------- merged report
def _merged():
    return {
        "module": "merged-qc-report",
        "n_reports": 2, "n_samples": 3,
        "batches": {"batchA": 2, "batchB": 1},
        "verdict_summary": {"samples": {"pass": 2, "warn": 1},
                            "warning_metrics": {"frip": 1}},
        "samples": {
            "batchA::s1": {"batch": "batchA", "assay": "ATAC-seq",
                           "qc_verdict": {"overall": "pass", "warnings": []}},
            "batchA::s2": {"batch": "batchA", "assay": "ATAC-seq",
                           "qc_verdict": {"overall": "warn",
                                          "warnings": ["frip"]}},
            "batchB::s3": {"batch": "batchB", "assay": "ChIP-seq",
                           "qc_verdict": {"overall": "pass", "warnings": []}},
        },
        "problems": ["batchC/qc_report.json: unreadable"],
    }


def test_merged_report_lists_every_sample_and_batch(tmp_path):
    out = render_merged_html(_merged(), tmp_path / "batch.html")
    parser, text = _audit(out)
    assert parser.external == []
    for name in ("batchA::s1", "batchA::s2", "batchB::s3"):
        assert name in text
    assert "batchA" in text and "batchB" in text


def test_merged_report_surfaces_skipped_reports(tmp_path):
    out = render_merged_html(_merged(), tmp_path / "batch.html")
    text = out.read_text(encoding="utf-8")
    # A silently dropped report is the failure mode this section exists for.
    assert "skipped" in text.lower()
    assert "batchC/qc_report.json" in text


def test_merged_report_with_no_samples_is_still_valid(tmp_path):
    empty = {"n_reports": 0, "n_samples": 0, "samples": {},
             "verdict_summary": {}, "problems": []}
    out = render_merged_html(empty, tmp_path / "batch.html")
    text = out.read_text(encoding="utf-8")
    assert "no samples" in text
    assert text.startswith("<!DOCTYPE html>")
