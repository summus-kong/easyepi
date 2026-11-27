"""Background construction and its report: the numbers a reader must be able to check."""
from __future__ import annotations

import json

import numpy as np
import pytest

from easyepi.background_report import (background_report, compare_distributions)
from easyepi.background_cli import _read_bed


# --------------------------------------------------------------------------
# Distribution summary
# --------------------------------------------------------------------------
def test_summary_reports_quartiles_and_not_a_fabricated_zero():
    from easyepi.background_report import _summarize

    s = _summarize([0.4, 0.5, 0.6])
    assert s["n"] == 3
    assert s["median"] == pytest.approx(0.5)
    assert s["min"] == 0.4 and s["max"] == 0.6

    empty = _summarize([])
    assert empty["n"] == 0 and empty["mean"] is None


def test_ks_flags_an_unmatched_background():
    """A GC-composition difference must be visible, not averaged away."""
    rng = np.random.default_rng(0)
    target = rng.normal(0.60, 0.03, 200)
    bg = rng.normal(0.40, 0.03, 200)
    d = compare_distributions(target, bg, [200] * 200, [200] * 200)
    assert d["gc"]["ks"]["pvalue"] < 1e-6
    assert d["gc"]["target"]["mean"] > d["gc"]["background"]["mean"]


def test_ks_is_quiet_when_the_background_really_matches():
    rng = np.random.default_rng(1)
    target = rng.normal(0.50, 0.03, 200)
    bg = rng.normal(0.50, 0.03, 200)
    d = compare_distributions(target, bg, [200] * 200, [200] * 200)
    assert d["gc"]["ks"]["pvalue"] > 0.01


def test_ks_reports_why_it_did_not_run():
    d = compare_distributions([], [0.5], [200], [200])
    assert d["gc"]["ks"]["statistic"] is None
    assert "note" in d["gc"]["ks"]


# --------------------------------------------------------------------------
# Report assembly
# --------------------------------------------------------------------------
def _targets(n=5):
    return [("chr1", i * 1000, i * 1000 + 200) for i in range(n)]


def test_report_without_sequences_says_gc_was_not_measured():
    """Reporting a GC comparison that never happened is the failure mode."""
    t = _targets()
    c = [("chr1", i * 1000 + 5000, i * 1000 + 5200) for i in range(5)]
    rep = background_report(t, c)
    assert rep["gc_source"] == "not measured"
    assert rep["distributions"]["gc"]["status"] == "not measured"
    # length is still compared, since it needs no sequences
    assert rep["distributions"]["length"]["target"]["median"] == 200


def test_report_measures_gc_when_sequences_are_supplied():
    t = _targets()
    c = [("chr1", i * 1000 + 5000, i * 1000 + 5200) for i in range(5)]
    reps = ["GCGC" * 50] * 5
    rep = background_report(t, c, target_seqs=reps, control_seqs=["ATAT" * 50] * 5)
    assert rep["gc_source"] == "measured on sequences"
    assert rep["distributions"]["gc"]["target"]["mean"] == 1.0
    assert rep["distributions"]["gc"]["background"]["mean"] == 0.0


def test_report_states_what_was_not_matched():
    rep = background_report(_targets(), _targets())
    for item in ("mappability", "repeat content", "chromatin state",
                 "accessibility"):
        assert item in rep["not_matched"]
    assert rep["composition_matched_only"] is True


def test_low_match_rate_produces_a_warning():
    """The bias that matters is the targets that silently got no control."""
    summary = {"n_matched": 40, "n_unmatched": 60, "frac_matched": 0.4,
               "tolerance": 0.05, "min_distance": 2000}
    rep = background_report(_targets(100), [], match_summary=summary)
    assert "warning" in rep
    assert "40%" in rep["warning"]


def test_high_match_rate_does_not_warn():
    summary = {"n_matched": 95, "n_unmatched": 5, "frac_matched": 0.95,
               "tolerance": 0.05, "min_distance": 2000}
    rep = background_report(_targets(100), [], match_summary=summary)
    assert "warning" not in rep


# --------------------------------------------------------------------------
# BED reading
# --------------------------------------------------------------------------
def test_read_bed_skips_track_lines_and_requires_regions(tmp_path):
    p = tmp_path / "p.bed"
    p.write_text("track name=x\n# comment\nchr1\t10\t20\nchr1\t30\t40\n",
                 encoding="utf-8")
    assert _read_bed(p) == [("chr1", 10, 20), ("chr1", 30, 40)]

    empty = tmp_path / "e.bed"
    empty.write_text("track name=y\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no regions"):
        _read_bed(empty)


# --------------------------------------------------------------------------
# End to end: does a real run produce a matched set on a synthetic genome?
# --------------------------------------------------------------------------
def test_background_end_to_end_matches_gc(tmp_path):
    """Known-truth check: targets drawn from a GC-rich half must get GC-rich controls."""
    from easyepi.background_cli import run_background
    from easyepi.gc_background import gc_fraction, iter_fasta

    # 40 kb genome: first half GC-rich, second half AT-rich, 1 kb blocks
    seq = ("GC" * 500) * 20 + ("AT" * 500) * 20
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + seq + "\n", encoding="utf-8")

    # targets inside the GC-rich half
    peaks = tmp_path / "p.bed"
    rows = ["chr1\t%d\t%d" % (i * 1000, i * 1000 + 200) for i in range(5)]
    peaks.write_text("\n".join(rows) + "\n", encoding="utf-8")

    out = tmp_path / "bg.bed"

    class A:
        pass

    a = A()
    a.peaks = str(peaks)
    a.genome_fasta = str(fa)
    a.window = 1000
    a.tolerance = 0.05
    a.min_distance = 0
    a.exclude = None
    a.any_chrom = False
    a.chromosomes = None
    a.seed = 0
    a.out = str(out)
    a.verbose = False
    a.quiet = True
    a.log_file = None

    assert run_background(a) == 0
    assert out.exists()

    controls = _read_bed(out)
    assert controls, "no controls were written"

    genome = dict(iter_fasta(fa))["chr1"]
    ctrl_gc = [gc_fraction(genome[s:e]) for _c, s, e in controls]
    # every control must come from the GC-rich half
    assert all(g > 0.9 for g in ctrl_gc), ctrl_gc

    report = json.loads(out.with_suffix(".bed.json").read_text())
    assert report["distributions"]["gc"]["target"]["mean"] > 0.9
    assert report["distributions"]["gc"]["background"]["mean"] > 0.9
    assert report["params"]["seed"] == 0
