"""BAM -> BigWig: bin accumulation, Tn5 shift, normalisation, CLI.

These tests need pysam to build a tiny BAM, so they skip cleanly where it is
unavailable. The arithmetic is checked against hand-computed expectations rather
than against the implementation's own output.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from easyepi.bam2bigwig import _sum_to_bins
from easyepi.bam2bw_cli import (ASSAY_DEFAULTS, _sample_names, resolve_defaults)

pysam = pytest.importorskip("pysam", reason="pysam is required to write test BAMs")


# --------------------------------------------------------------------------
# Bin arithmetic (no BAM needed)
# --------------------------------------------------------------------------
def test_sum_to_bins_aligns_chunk_to_the_global_grid():
    """A chunk starting mid-bin must credit the partial bin, not shift the grid."""
    depth = np.ones(10, dtype=float)      # 10 bases of depth 1
    # bin_size 4, chunk starts at 2 -> bins [2..6), [6..10), [10..14)
    out = _sum_to_bins(depth, offset=2, bin_size=4, n_bins=4)
    assert out[0] == pytest.approx(2.0)   # bases 2,3
    assert out[1] == pytest.approx(4.0)   # bases 4-7
    assert out[2] == pytest.approx(4.0)   # bases 8-11 -> only 8,9,10,11 present
    assert out.sum() == pytest.approx(10.0)


def test_sum_to_bins_preserves_total_signal_across_chunks():
    depth = np.arange(1, 21, dtype=float)
    total = float(depth.sum())
    for offset in (0, 1, 3, 7):
        n_bins = 20
        out = _sum_to_bins(depth, offset, bin_size=5, n_bins=n_bins)
        assert out.sum() == pytest.approx(total), f"offset {offset} lost signal"


def test_sum_to_bins_handles_empty_input():
    out = _sum_to_bins(np.zeros(0), 0, 10, 5)
    assert out.shape == (5,)
    assert out.sum() == 0.0


# --------------------------------------------------------------------------
# Assay defaults
# --------------------------------------------------------------------------
def test_atac_defaults_use_the_tn5_shift():
    d = resolve_defaults("atac")
    assert d["shift"] == "tn5"
    assert d["normalize"] == "cpm"


def test_chip_defaults_use_the_read_body():
    assert resolve_defaults("chip")["shift"] == "none"
    assert resolve_defaults("histone")["shift"] == "none"


def test_fallback_defaults_do_not_alter_the_signal():
    d = resolve_defaults(None)
    assert d["shift"] == "none" and d["normalize"] == "raw"
    assert d["source"] == "fallback"


def test_unknown_assay_lists_the_choices():
    with pytest.raises(ValueError, match="unknown assay"):
        resolve_defaults("rna-seq")


def test_every_assay_default_is_a_valid_combination():
    for assay, d in ASSAY_DEFAULTS.items():
        assert d["shift"] in ("none", "tn5"), assay
        assert d["normalize"] in ("raw", "cpm", "rpkm"), assay
        assert d["bin_size"] >= 1, assay


# --------------------------------------------------------------------------
# Sample naming
# --------------------------------------------------------------------------
def test_sample_names_default_to_the_file_stem():
    names = _sample_names(["a/rep1.bam", "b/rep2.bam"], None)
    assert names == {"a/rep1.bam": "rep1", "b/rep2.bam": "rep2"}


def test_name_override_changes_only_the_named_file():
    names = _sample_names(["a/rep1.bam", "b/rep2.bam"], ["a/rep1.bam=ctrla"])
    assert names["a/rep1.bam"] == "ctrla"
    assert names["b/rep2.bam"] == "rep2"


def test_name_override_must_refer_to_a_given_bam():
    with pytest.raises(ValueError, match="not among the BAM arguments"):
        _sample_names(["a.bam"], ["zzz.bam=x"])


def test_duplicate_names_are_rejected():
    with pytest.raises(ValueError, match="duplicate sample name"):
        _sample_names(["x/s1.bam", "y/s1.bam"], None)


# --------------------------------------------------------------------------
# Real BAM round trip
# --------------------------------------------------------------------------
def _write_bam(path: Path, reads, chrom: str = "chr1", length: int = 1000,
               paired: bool = False) -> Path:
    """Create a coordinate-sorted, indexed BAM.

    ``reads`` is a list of ``(start, is_reverse)``; read length is 10 bp.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": chrom, "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, is_reverse) in enumerate(sorted(reads)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 10
            a.flag = 16 if is_reverse else 0
            if paired:
                a.flag |= 1 | 64          # paired, read1
            a.reference_id = 0
            a.reference_start = int(start)
            a.mapping_quality = 60
            a.cigartuples = [(0, 10)]
            out.write(a)
    pysam.index(str(path))
    return path


def test_bam_to_bigwig_writes_a_readable_track(tmp_path):
    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(100, False), (300, True)])
    bw_path = tmp_path / "s.bw"
    meta = bam_to_bigwig(bam, bw_path, bin_size=10, shift="none", normalize="raw")

    assert bw_path.exists() and bw_path.stat().st_size > 0
    assert meta["bin_size"] == 10 and meta["shift"] == "none"
    assert meta["mapped_reads"]["mapped"] == 2

    import pyBigWig
    bw = pyBigWig.open(str(bw_path))
    try:
        assert bw.chroms()["chr1"] == 1000
        # 10 reads x ... one read of 10 bp at 100 covers 100..109 -> 10 bases
        vals = bw.stats("chr1", 0, 1000, type="sum", nBins=100)
        assert sum(v for v in vals if v) > 0
    finally:
        bw.close()


def test_shift_none_and_tn5_place_signal_differently(tmp_path):
    """The Tn5 path moves the forward read +4 and the reverse read -5."""
    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(500, False), (600, True)])
    plain = tmp_path / "plain.bw"
    shifted = tmp_path / "shifted.bw"
    bam_to_bigwig(bam, plain, bin_size=1, shift="none", normalize="raw")
    bam_to_bigwig(bam, shifted, bin_size=1, shift="tn5", normalize="raw")

    import pyBigWig

    def first_nonzero(path):
        bw = pyBigWig.open(str(path))
        try:
            vals = bw.values("chr1", 0, 1000)
            for i, v in enumerate(vals):
                if v:
                    return i
        finally:
            bw.close()
        return None

    # read body starts at 500; the Tn5 forward cut site is 504
    assert first_nonzero(plain) == 500
    assert first_nonzero(shifted) == 504


def test_tn5_reverse_read_uses_the_minus_five_offset(tmp_path):
    from easyepi.bam2bigwig import bam_to_bigwig

    # reverse read 600..609, so reference_end is 610 -> 610-5 = 605
    bam = _write_bam(tmp_path / "s.bam", [(600, True)])
    bw_path = tmp_path / "s.bw"
    bam_to_bigwig(bam, bw_path, bin_size=1, shift="tn5", normalize="raw")

    import pyBigWig
    bw = pyBigWig.open(str(bw_path))
    try:
        vals = bw.values("chr1", 0, 1000)
        nonzero = [i for i, v in enumerate(vals) if v]
    finally:
        bw.close()
    assert nonzero == [605]


def test_cpm_scaling_is_recorded_and_applied(tmp_path):
    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(100, False)] * 4)
    bw_path = tmp_path / "s.bw"
    meta = bam_to_bigwig(bam, bw_path, bin_size=10, shift="none", normalize="cpm")
    # 4 mapped reads -> scale = 1e6/4
    assert meta["scale_factor"] == pytest.approx(1e6 / 4)
    assert "per million" in meta["signal_units"]


def test_unknown_chromosome_is_reported(tmp_path):
    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(100, False)])
    with pytest.raises(ValueError, match="not in the BAM header"):
        bam_to_bigwig(bam, tmp_path / "s.bw", chrom="chrZ")


def test_invalid_parameters_are_rejected(tmp_path):
    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(100, False)])
    with pytest.raises(ValueError, match="normalize must be"):
        bam_to_bigwig(bam, tmp_path / "a.bw", normalize="zscore")
    with pytest.raises(ValueError, match="shift must be"):
        bam_to_bigwig(bam, tmp_path / "b.bw", shift="weird")
    with pytest.raises(ValueError, match="bin_size"):
        bam_to_bigwig(bam, tmp_path / "c.bw", bin_size=0)


def test_metadata_sidecar_records_every_value_affecting_the_signal(tmp_path):
    import json

    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(100, False)])
    bw_path = tmp_path / "s.bw"
    bam_to_bigwig(bam, bw_path, bin_size=25, shift="tn5", normalize="cpm",
                  min_mapq=20)
    meta = json.loads(Path(str(bw_path) + ".json").read_text())
    for key in ("bam", "bigwig", "bin_size", "normalize", "shift", "min_mapq",
                "scale_factor", "mapped_reads", "signal_units"):
        assert key in meta
    assert meta["bin_size"] == 25 and meta["shift"] == "tn5"


def test_blacklist_regions_are_zeroed(tmp_path):
    from easyepi.bam2bigwig import bam_to_bigwig

    bam = _write_bam(tmp_path / "s.bam", [(100, False)])
    bw_path = tmp_path / "s.bw"
    bam_to_bigwig(bam, bw_path, bin_size=10, shift="none",
                  blacklist=[("chr1", 0, 200)])
    import pyBigWig
    bw = pyBigWig.open(str(bw_path))
    try:
        vals = bw.stats("chr1", 0, 200, type="sum", nBins=20)
    finally:
        bw.close()
    assert all(v in (0, 0.0, None) for v in vals)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_is_registered_with_expected_defaults():
    from easyepi import cli
    parser = cli.build_parser()
    a = parser.parse_args(["bam2bw", "x.bam", "-o", "tracks/"])
    assert a.bam == ["x.bam"] and a.out == "tracks/"
    assert a.shift is None and a.bin_size is None   # filled from --assay


def _cli_args(tmp_path, **overrides):
    """Namespace for run_bam2bw; built from defaults so each test states only
    what it changes. A class body cannot close over the test's locals, which is
    what made the earlier version raise NameError."""
    from types import SimpleNamespace

    base = dict(bam=[], out=str(tmp_path), output_file=None, assay="atac",
                shift=None, bin_size=None, normalize=None, chrom=None,
                min_mapq=30, max_depth=None, name=None, manifest=None)
    base.update(overrides)
    return SimpleNamespace(**base)


def test_cli_runs_end_to_end(tmp_path):
    from easyepi.bam2bw_cli import run_bam2bw

    bam = _write_bam(tmp_path / "rep1.bam", [(100, False), (300, True)])
    out = tmp_path / "tracks"
    args = _cli_args(tmp_path, bam=[str(bam)], out=str(out),
                     manifest=str(tmp_path / "tracks.json"))
    assert run_bam2bw(args) == 0
    assert (out / "rep1.bw").exists()
    assert (tmp_path / "tracks.json").exists()


def test_cli_writes_one_track_per_bam(tmp_path):
    from easyepi.bam2bw_cli import run_bam2bw

    a = _write_bam(tmp_path / "a.bam", [(100, False)])
    b = _write_bam(tmp_path / "b.bam", [(200, False)])
    out = tmp_path / "tracks"
    assert run_bam2bw(_cli_args(tmp_path, bam=[str(a), str(b)],
                                out=str(out))) == 0
    assert (out / "a.bw").exists() and (out / "b.bw").exists()


def test_cli_rejects_output_file_with_several_bams(tmp_path):
    from easyepi.bam2bw_cli import run_bam2bw

    bam = _write_bam(tmp_path / "a.bam", [(100, False)])
    args = _cli_args(tmp_path, bam=[str(bam), str(bam)],
                     output_file=str(tmp_path / "one.bw"))
    with pytest.raises(ValueError, match="single BigWig"):
        run_bam2bw(args)


def test_cli_reports_a_missing_bam(tmp_path):
    from easyepi.bam2bw_cli import run_bam2bw

    args = _cli_args(tmp_path, bam=[str(tmp_path / "absent.bam")])
    with pytest.raises(FileNotFoundError, match="input not found"):
        run_bam2bw(args)

