"""GC-matched background: FASTA parsing, window GC, matching, exclusions."""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.gc_background import (gc_fraction, genome_windows, iter_fasta,
                                   match_background, write_bed)


def _fa(tmp_path, name, records):
    p = tmp_path / name
    p.write_text("".join(f">{h}\n{s}\n" for h, s in records), encoding="utf-8")
    return p


# --------------------------------------------------------------------------
# GC fraction
# --------------------------------------------------------------------------
def test_gc_fraction_basic():
    assert gc_fraction("GCGC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert gc_fraction("ACGT") == 0.5


def test_gc_fraction_ignores_n_and_is_case_insensitive():
    assert gc_fraction("ggNNcc") == 1.0
    assert gc_fraction("Nnnn") == 0.0
    assert gc_fraction("") == 0.0


# --------------------------------------------------------------------------
# FASTA reader
# --------------------------------------------------------------------------
def test_iter_fasta_yields_uppercased_sequences(tmp_path):
    p = _fa(tmp_path, "g.fa", [("chr1", "acgt"), ("chr2 desc", "GGCC")])
    got = list(iter_fasta(p))
    assert got == [("chr1", "ACGT"), ("chr2", "GGCC")]


def test_iter_fasta_handles_wrapped_lines(tmp_path):
    p = tmp_path / "w.fa"
    p.write_text(">chr1\nACGT\nACGT\nACGT\n", encoding="utf-8")
    assert list(iter_fasta(p)) == [("chr1", "ACGTACGTACGT")]


# --------------------------------------------------------------------------
# Windows
# --------------------------------------------------------------------------
def test_genome_windows_computes_gc_per_window(tmp_path):
    p = _fa(tmp_path, "g.fa", [("chr1", "G" * 100 + "A" * 100)])
    w = genome_windows(p, window=100)
    assert len(w["chr1"]) == 2
    assert w["chr1"][0] == (0, 100, 1.0)
    assert w["chr1"][1] == (100, 200, 0.0)


def test_genome_windows_drops_n_dominated_windows(tmp_path):
    p = _fa(tmp_path, "g.fa", [("chr1", "G" * 100 + "N" * 100)])
    w = genome_windows(p, window=100)
    assert len(w["chr1"]) == 1


def test_genome_windows_filters_chromosomes(tmp_path):
    p = _fa(tmp_path, "g.fa", [("chr1", "G" * 50), ("chr2", "A" * 50)])
    w = genome_windows(p, window=50, chromosomes=["chr2"])
    assert set(w) == {"chr2"}


def test_genome_windows_validates_window_and_empties(tmp_path):
    p = _fa(tmp_path, "g.fa", [("chr1", "G" * 50)])
    with pytest.raises(ValueError, match="window"):
        genome_windows(p, window=0)
    empty = _fa(tmp_path, "n.fa", [("chr1", "N" * 50)])
    with pytest.raises(ValueError, match="no usable windows"):
        genome_windows(empty, window=50)


# --------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------
def _two_band_genome(tmp_path):
    """chr1: a GC-rich half then an AT-rich half, 100 bp windows."""
    gc_half = ("GC" * 50)    # 100 bp, GC = 1.0
    at_half = ("AT" * 50)    # 100 bp, GC = 0.0
    return {"chr1": [(0, 100, 1.0), (100, 200, 0.0)]}


def test_match_prefers_the_same_gc_band():
    """With two GC-rich windows and one AT-rich window, a GC-rich target must
    receive a GC-rich control and never itself."""
    windows = {"chr1": [(0, 100, 1.0), (5000, 5100, 1.0), (9000, 9100, 0.0)]}
    targets = [("chr1", 0, 100)]   # GC-rich target
    controls, summary = match_background(targets, windows, tolerance=0.1,
                                         min_distance=0, seed=1)
    assert summary["n_matched"] == 1
    _chrom, cs, ce = controls[0]
    assert (cs, ce) != (0, 100)               # never the target itself
    assert abs(_gc_of(windows, cs, ce) - 1.0) <= 0.1


def test_match_reports_no_candidate_when_only_the_target_matches():
    """A genome with a single GC-rich window equal to the target offers no match."""
    windows = _two_band_genome(None)
    controls, summary = match_background([("chr1", 0, 100)], windows,
                                         tolerance=0.1, min_distance=0, seed=1)
    assert controls == []
    assert summary["n_unmatched"] == 1


def test_match_respects_min_distance():
    windows = _two_band_genome(None)
    targets = [("chr1", 0, 100)]
    # every window is within 500 bp of the target, so nothing can be matched
    controls, summary = match_background(targets, windows, tolerance=0.5,
                                         min_distance=500, seed=1)
    assert controls == []
    assert summary["n_unmatched"] == 1


def test_match_never_returns_the_same_window_twice():
    windows = {"chr1": [(i * 100, i * 100 + 100, 1.0) for i in range(6)]}
    targets = [("chr1", 0, 100)] * 8   # more targets than candidates
    controls, summary = match_background(targets, windows, tolerance=0.1,
                                         min_distance=0, seed=2,
                                         max_tries_per_target=400)
    assert len(controls) == len(set(controls))
    assert summary["n_matched"] <= 6
    assert summary["n_matched"] + summary["n_unmatched"] == len(targets)


def test_match_reports_unmatched_targets_instead_of_dropping_them():
    windows = {"chr1": [(0, 100, 0.0)]}
    targets = [("chr1", 0, 100), ("chrX", 0, 100)]
    controls, summary = match_background(targets, windows, tolerance=0.1,
                                         min_distance=0, seed=0)
    assert summary["n_unmatched"] == 2
    reasons = {u["reason"] for u in summary["unmatched"]}
    assert any("not in genome windows" in r for r in reasons)


def test_match_excludes_regions_from_candidates():
    windows = {"chr1": [(i * 100, i * 100 + 100, 1.0) for i in range(5)]}
    targets = [("chr1", 0, 100)]
    exclude = [("chr1", i * 100, i * 100 + 100) for i in range(1, 5)]
    controls, summary = match_background(targets, windows, tolerance=0.1,
                                         min_distance=0, seed=3, exclude=exclude,
                                         max_tries_per_target=400)
    assert controls == []


def test_match_validates_parameters():
    windows = {"chr1": [(0, 100, 0.5)]}
    with pytest.raises(ValueError, match="tolerance"):
        match_background([("chr1", 0, 100)], windows, tolerance=0.0)
    with pytest.raises(ValueError, match="min_distance"):
        match_background([("chr1", 0, 100)], windows, min_distance=-5)


def test_match_summary_fraction_is_consistent():
    windows = {"chr1": [(i * 1000, i * 1000 + 1000, 0.5) for i in range(50)]}
    targets = [("chr1", i * 1000, i * 1000 + 1000) for i in range(10)]
    _controls, summary = match_background(targets, windows, tolerance=0.2,
                                          min_distance=0, seed=7)
    assert summary["frac_matched"] == pytest.approx(
        summary["n_matched"] / summary["n_targets"])


def test_write_bed_roundtrip(tmp_path):
    out = write_bed([("chr1", 10, 20), ("chr1", 30, 40)], tmp_path / "bg.bed")
    rows = out.read_text().strip().splitlines()
    assert rows[0].split("\t")[:3] == ["chr1", "10", "20"]
    assert rows[0].split("\t")[3] == "bg_1"


def _gc_of(windows, s, e):
    for ws, we, gc in windows["chr1"]:
        if ws == s and we == e:
            return gc
    return None
