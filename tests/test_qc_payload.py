"""Direct tests for the shared QC payload builder and the entry-point helper.

``build_sample_payload`` is the single source of truth behind three QC command
paths, so it deserves direct coverage rather than only indirect exercise through
the CLI. ``_entry`` is covered because every ``easyepi.universal.*`` script now
depends on it.
"""

from __future__ import annotations

import json
import sys

import pytest

pytest.importorskip("pysam")
pytest.importorskip("pyBigWig")

import pysam  # noqa: E402
import pyBigWig  # noqa: E402

from easyepi import qc_payload  # noqa: E402
from easyepi.qc_thresholds import load_rules  # noqa: E402

CHROM = "chr1"
LENGTH = 100_000


def _bam(path, records):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": LENGTH}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, isize) in enumerate(sorted(records)):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 60
            a.query_qualities = pysam.qualitystring_to_array("I" * 60)
            a.flag = 99 if i % 2 == 0 else 147
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = "60M"
            a.template_length = isize if i % 2 == 0 else -isize
            out.write(a)
    pysam.index(str(path))


def _bed(path, intervals):
    with open(path, "w", encoding="utf-8") as fh:
        for s, e in intervals:
            fh.write(f"{CHROM}\t{s}\t{e}\tpeak\n")


def _bigwig(path, values):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(CHROM, len(values))])
    starts = list(range(len(values)))
    bw.addEntries([CHROM] * len(values), starts,
                  ends=[s + 1 for s in starts],
                  values=[float(v) for v in values])
    bw.close()


# ---------------------------------------------------------------- basics
def test_no_assay_returns_empty_payload_and_no_verdict():
    out = qc_payload.build_sample_payload(None)
    assert out["assay_qc"] == {}
    assert out["qc_verdict"] is None


def test_unknown_assay_returns_empty_but_verdict_present():
    out = qc_payload.build_sample_payload("RNA-seq")
    assert out["assay_qc"] == {}
    # A verdict is still produced so callers can see "everything not_evaluated".
    assert out["qc_verdict"] is not None
    assert out["qc_verdict"]["overall"] == "not_evaluated"


# ---------------------------------------------------------------- paths
def test_missing_files_are_tolerated(tmp_path):
    """Nonexistent inputs must not raise: metrics just come back empty."""
    out = qc_payload.build_sample_payload(
        "ChIP-seq",
        bam=tmp_path / "nope.bam",
        peak=tmp_path / "nope.bed",
        bigwig=tmp_path / "nope.bw",
    )
    assert isinstance(out["assay_qc"], dict)
    assert out["qc_verdict"] is not None


def test_existing_path_is_resolved(tmp_path):
    bam = tmp_path / "r.bam"
    _bam(bam, [(1000, 150), (2000, 160)])
    assert qc_payload._existing(bam) == str(bam)
    assert qc_payload._existing(tmp_path / "missing") is None
    assert qc_payload._existing(None) is None


# ---------------------------------------------------------------- metrics
def test_chip_payload_includes_frip_and_verdict(tmp_path):
    bam = tmp_path / "r.bam"
    _bam(bam, [(1000, 150), (1010, 150), (50_000, 150)])
    peak = tmp_path / "p.bed"
    _bed(peak, [(950, 1200)])
    out = qc_payload.build_sample_payload("ChIP", bam=bam, peak=peak)
    assert "frip" in out["assay_qc"]
    assert out["qc_verdict"]["overall"] in {"pass", "warn"}


def test_atac_alias_resolves_to_atac_profile(tmp_path):
    bam = tmp_path / "r.bam"
    _bam(bam, [(1000, 150), (2000, 200)])
    out = qc_payload.build_sample_payload("ATAC", bam=bam)
    # ATAC profile carries fragment-length structure that ChIP does not.
    assert ("median_fragment_length" in out["assay_qc"]
            or "fraction_nfr" in out["assay_qc"])


def test_bam_qc_folded_into_verdict(tmp_path):
    bam = tmp_path / "r.bam"
    _bam(bam, [(1000, 150)])
    bam_qc = {"total_reads": 100, "mapping_rate": 0.10}
    out = qc_payload.build_sample_payload("ChIP-seq", bam=bam, bam_qc=bam_qc)
    breached = {c["metric"] for c in out["qc_verdict"]["checks"] if c["status"] == "breached"}
    assert "bam_mapping_rate" in breached


def test_custom_rules_are_honoured(tmp_path):
    bam = tmp_path / "r.bam"
    _bam(bam, [(1000, 150), (1010, 150)])
    peak = tmp_path / "p.bed"
    _bed(peak, [(950, 1200)])
    relaxed = load_rules({"frip": {"op": "lt", "value": 0.0}})
    out = qc_payload.build_sample_payload("ChIP-seq", bam=bam, peak=peak, rules=relaxed)
    assert "frip" not in out["qc_verdict"]["warnings"]


# ---------------------------------------------------------------- methylation helper
def test_methylation_beta_values_handles_missing_column():
    class FakeTable:
        columns = ["pos"]
        def __getitem__(self, key):
            raise KeyError(key)
    assert qc_payload.methylation_beta_values(FakeTable()) is None
    assert qc_payload.methylation_beta_values(None) is None


def test_methylation_beta_values_extracts_floats():
    import pandas as pd
    frame = pd.DataFrame({"beta": [0.1, 0.9, None]})
    assert qc_payload.methylation_beta_values(frame) == pytest.approx([0.1, 0.9])


# ---------------------------------------------------------------- merge_into_row
def test_merge_into_row_merges_payload():
    row = {"sample_id": "s1", "assay_qc": {"frip": 0.5}}
    result = {"assay_qc": {"tss_enrichment": 8.0}, "qc_verdict": {"overall": "pass"}}
    merged = qc_payload.merge_into_row(row, result)
    assert merged["assay_qc"] == {"frip": 0.5, "tss_enrichment": 8.0}
    assert merged["qc_verdict"]["overall"] == "pass"


def test_merge_into_row_keeps_existing_verdict_when_none():
    row = {"sample_id": "s1", "qc_verdict": {"overall": "warn"}}
    qc_payload.merge_into_row(row, {"assay_qc": {}, "qc_verdict": None})
    assert row["qc_verdict"]["overall"] == "warn"


# ---------------------------------------------------------------- _entry helper
def test_output_directory_treats_file_suffix_as_file(tmp_path):
    from easyepi._entry import output_directory
    assert output_directory(tmp_path / "x.tsv") == tmp_path
    assert output_directory(tmp_path / "x.json") == tmp_path
    assert output_directory(tmp_path / "dir") == tmp_path / "dir"


def test_redirect_out_to_file_only_when_directory(tmp_path):
    from easyepi._entry import redirect_out_to_file

    class A:
        out = str(tmp_path / "outdir")
    args = A()
    redirect_out_to_file(args, "result.tsv")
    assert args.out.endswith("result.tsv")

    class B:
        out = str(tmp_path / "explicit.tsv")
    args_b = B()
    redirect_out_to_file(args_b, "result.tsv")
    assert args_b.out.endswith("explicit.tsv")


def test_run_entry_writes_provenance(tmp_path):
    """run_entry must produce the same command/metrics/status trio as before."""
    from easyepi._entry import run_entry

    class Args:
        out = str(tmp_path / "run")
        verbose = False
        quiet = True
        log_file = None

    with run_entry("TEST-MODULE", Args(), samples=3) as run:
        (tmp_path / "run" / "artifact.txt").write_text("x", encoding="utf-8")

    outdir = tmp_path / "run"
    assert (outdir / "status.json").exists()
    metrics = json.loads((outdir / "metrics.json").read_text(encoding="utf-8"))
    assert metrics.get("samples") == 3
