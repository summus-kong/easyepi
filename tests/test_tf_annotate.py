"""TF annotation: exact PWM thresholds, coordinate mapping, enrichment, CLI."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from easyepi.motifs.formats import Motif
from easyepi.tf_annotate import (annotate_sequences, enrichment_test,
                                 load_motif_database, pwm_score_distribution,
                                 pwm_threshold, read_peak_intervals,
                                 site_pvalue)


def _motif(name: str, columns: list[dict[str, float]], mid: str = "MA0000.1") -> Motif:
    """Build a Motif from per-column base counts given as {base: count}."""
    mat = np.zeros((4, len(columns)))
    for j, col in enumerate(columns):
        for i, base in enumerate("ACGT"):
            mat[i, j] = col.get(base, 0.0)
    return Motif(name=name, matrix=mat, id=mid, source="test")


# A fully degenerate "AAAA" motif: one possible score, probability 0.25**4 each.
_AAAA = _motif("AAAA", [{"A": 1.0}] * 4)
# A perfectly informative "ACGT" motif.
_ACGT = _motif("ACGT", [{"A": 1.0}, {"C": 1.0}, {"G": 1.0}, {"T": 1.0}])
# A two-base variety so the score distribution has repeated values.
_GATA = _motif("GATA", [{"G": 1.0}, {"A": 1.0}, {"T": 1.0}, {"A": 1.0}])


# --------------------------------------------------------------------------
# Exact score distribution
# --------------------------------------------------------------------------
def test_score_distribution_is_a_normalised_pmf():
    probs_sum = pwm_score_distribution(_GATA.to_log_odds())[1].sum()
    assert probs_sum == pytest.approx(1.0, abs=1e-12)


def test_degenerate_motif_outcomes_are_ranked_by_column_penalty():
    """A single-base motif has L+1 outcomes: 0..L bases matching the consensus.

    Non-consensus bases are floored at 1e-9, so each mismatch contributes the same
    large negative score. Binomial counts mean 1-of-4 matching is the modal outcome
    while the exact consensus site (all L matching) is the rarest and highest.
    """
    scores, probs = pwm_score_distribution(_AAAA.to_log_odds())
    assert scores.size == _AAAA.length + 1
    assert scores[-1] == scores.max()
    # the exact consensus site: probability 0.25**L, and it is the top score
    assert probs[-1] == pytest.approx(0.25 ** _AAAA.length, rel=1e-6)
    assert probs.sum() == pytest.approx(1.0, abs=1e-12)
    # scores increase with the number of matching columns
    assert np.all(np.diff(scores) > 0)


def test_informative_motif_has_the_same_outcome_structure():
    scores, probs = pwm_score_distribution(_ACGT.to_log_odds())
    assert scores.size == _ACGT.length + 1
    # the exact "ACGT" site is the best and has probability 0.25**4
    assert probs[-1] == pytest.approx(0.25 ** _ACGT.length, rel=1e-6)
    assert scores[-1] > scores[0]


def test_distribution_matches_brute_force_convolution():
    """Independent check: enumerate all 4**L sequences for a short PWM."""
    motif = _motif("mix", [{"A": 3.0, "C": 1.0}, {"G": 2.0, "T": 2.0},
                           {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0}])
    pwm = motif.to_log_odds()
    from itertools import product
    observed = np.zeros(pwm.shape[1] * 4)
    # every base at every position, unweighted -> uniform 0.25 per base
    scores = []
    for combo in product(range(4), repeat=pwm.shape[1]):
        score = 0.0
        for j, b in enumerate(combo):
            score += pwm[b, j]
        scores.append(round(score, 6))
    uniq, counts = np.unique(scores, return_counts=True)
    expected_probs = counts / counts.sum()
    got_scores, got_probs = pwm_score_distribution(pwm)
    # The score axis is binned (as FIMO does), so scores agree to within one bin
    # width rather than exactly; the probabilities are what must match.
    bin_width = float(np.ptp(pwm, axis=0).sum()) / 20_000
    assert np.allclose(uniq, got_scores, atol=2 * bin_width)
    assert np.allclose(expected_probs, got_probs, atol=1e-9)
    assert got_probs.sum() == pytest.approx(1.0, abs=1e-12)


# --------------------------------------------------------------------------
# Thresholds
# --------------------------------------------------------------------------
def test_threshold_tail_is_at_most_requested_fpr():
    pwm = _GATA.to_log_odds()
    thr, achieved = pwm_threshold(pwm, 0.05)
    assert 0 < achieved <= 0.05
    scores, probs = pwm_score_distribution(pwm)
    tail = np.cumsum(probs[::-1])[::-1]
    assert tail[np.searchsorted(scores, thr, side="left")] == pytest.approx(achieved)


def test_threshold_rejects_unreachable_fpr():
    """A 4 bp motif's best site has probability 0.25**4 = 0.0039.

    So fpr=1e-4 is unreachable and must raise, not silently pass everything.
    """
    with pytest.raises(ValueError, match="no PWM score reaches"):
        pwm_threshold(_AAAA.to_log_odds(), 1e-4)
    # a reachable rate succeeds and reports the achieved tail
    thr, achieved = pwm_threshold(_AAAA.to_log_odds(), 0.01)
    assert achieved == pytest.approx(0.25 ** 4, rel=1e-6)


def test_threshold_validates_fpr_range():
    for bad in (0.0, 1.0, -0.1):
        with pytest.raises(ValueError, match="fpr must be"):
            pwm_threshold(_GATA.to_log_odds(), bad)


def test_site_pvalue_is_monotone_in_score():
    pwm = _GATA.to_log_odds()
    scores, _ = pwm_score_distribution(pwm)
    ps = [site_pvalue(pwm, s) for s in scores]
    assert all(a >= b for a, b in zip(ps, ps[1:]))


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------
def test_annotation_finds_the_planted_site_and_reports_forward_coordinates():
    seq = "TTTT" + "GATA" + "CCCC"
    ann = annotate_sequences([seq], [_GATA], fpr=0.25, top=5)
    rows = [r for r in ann["per_peak"] if r["seq_index"] == 0]
    assert rows, "planted GATA site was not detected"
    # the motif is at forward offset 4 regardless of which strand scored best
    assert rows[0]["offset_start"] == 4
    assert rows[0]["offset_end"] == 8
    assert rows[0]["strand"] in {"+", "-"}


def test_annotation_reports_nothing_when_no_site_passes():
    ann = annotate_sequences(["CCCCCCCCCC"], [_GATA], fpr=0.25, top=5)
    assert ann["per_peak"] == []
    assert ann["summary"][0]["n_peaks_with_site"] == 0


def test_annotation_respects_top():
    seqs = ["GATA" * 6]
    ann = annotate_sequences(seqs, [_GATA, _ACGT, _AAAA], fpr=0.25, top=2)
    assert len(ann["per_peak"]) <= 2
    assert all(r["rank"] <= 2 for r in ann["per_peak"])


def test_annotation_rejects_empty_input():
    with pytest.raises(ValueError, match="no sequences"):
        annotate_sequences([], [_GATA])
    with pytest.raises(ValueError, match="top must be"):
        annotate_sequences(["ACGT"], [_GATA], top=0)


def test_lowercase_and_n_bases_do_not_crash():
    ann = annotate_sequences(["gatannnnGATA".lower()], [_GATA], fpr=0.25, top=3)
    assert ann["n_sequences"] == 1


# --------------------------------------------------------------------------
# Peaks, coordinates, enrichment
# --------------------------------------------------------------------------
def test_read_peak_intervals_uses_bed_name_and_skips_bad_rows(tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text("track name=x\n"
                   "chr1\t100\t200\tpeakA\t10\t.\n"
                   "# comment\n"
                   "chr1\tbad\t300\n"
                   "chr2\t5\t5\n"
                   "chr2\t50\t70\n", encoding="utf-8")
    recs = read_peak_intervals(bed)
    assert [r["name"] for r in recs] == ["peakA", "peak_6"]
    assert recs[1]["chrom"] == "chr2" and recs[1]["start"] == 50


def test_read_peak_intervals_errors_on_empty_file(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("# nothing\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no intervals"):
        read_peak_intervals(bed)


def test_enrichment_marks_a_clearly_enriched_motif_and_corrects_pvalues():
    ann = {"summary": [
        {"tf": "GATA", "motif_id": "M1", "n_peaks_with_site": 40,
         "pct_peaks_with_site": 80.0},
        {"tf": "NOPE", "motif_id": "M2", "n_peaks_with_site": 5,
         "pct_peaks_with_site": 10.0},
    ]}
    rows = enrichment_test(ann, n_peaks=50, bg_present={"GATA": 5, "NOPE": 5}, n_bg=50)
    by_tf = {r["tf"]: r for r in rows}
    assert by_tf["GATA"]["pvalue"] < by_tf["NOPE"]["pvalue"]
    assert by_tf["GATA"]["odds_ratio"] > 1
    assert all(r["padj"] is not None for r in rows)
    # BH never makes a p-value smaller than the raw one
    assert all(r["padj"] >= r["pvalue"] - 1e-12 for r in rows)


def test_enrichment_skips_motifs_absent_from_both_sets():
    ann = {"summary": [{"tf": "Z", "motif_id": "MZ", "n_peaks_with_site": 0,
                        "pct_peaks_with_site": 0.0}]}
    assert enrichment_test(ann, n_peaks=10, bg_present={}, n_bg=10) == []


# --------------------------------------------------------------------------
# Motif database loading
# --------------------------------------------------------------------------
def test_load_motif_database_merges_comma_separated_files(tmp_path):
    a = tmp_path / "a.pfm"
    a.write_text(">MA0001.1 Foo\n1 4 0 0\n0 0 0 4\n0 0 4 0\n4 0 0 0\n", encoding="utf-8")
    b = tmp_path / "b.pfm"
    b.write_text(">MA0002.1 Bar\n4 0 0 0\n0 4 0 0\n0 0 4 0\n0 0 0 4\n", encoding="utf-8")
    motifs = load_motif_database(f"{a},{b}")
    assert [m.name for m in motifs] == ["Foo", "Bar"]


# --------------------------------------------------------------------------
# CLI integration (no pysam needed: the FASTA path is faked)
# --------------------------------------------------------------------------
def test_cli_writes_per_peak_summary_and_json(tmp_path, monkeypatch):
    from easyepi import tf_annotate

    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\n" + "N" * 100 + "GATA" + "N" * 200 + "\n")
    motif_file = tmp_path / "m.pfm"
    motif_file.write_text(">MA0001.1 GATA\n0 4 0 4\n0 0 0 0\n4 0 0 0\n0 0 4 0\n",
                          encoding="utf-8")
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t90\t120\tpeakA\n", encoding="utf-8")

    monkeypatch.setattr(tf_annotate, "fetch_sequences",
                        lambda records, fasta: ["N" * 10 + "GATA" + "N" * 16])

    class Args:
        peaks = str(bed)
        genome_fasta = str(genome)
        motifs = str(motif_file)
        fpr = 0.3
        background = None
        top = 5
        out = str(tmp_path / "out")
        figure_formats = None

    assert tf_annotate.run_tf_annotate(Args) == 0
    out = tmp_path / "out"
    assert (out / "tf_per_peak.tsv").exists()
    assert (out / "tf_summary.tsv").exists()
    payload = json.loads((out / "tf_annotate.json").read_text())
    assert payload["n_peaks_used"] == 1
    assert payload["score_model"].startswith("exact PWM convolution")
    rows = (out / "tf_per_peak.tsv").read_text().strip().splitlines()
    assert rows[0].split("\t")[:5] == ["peak_id", "chrom", "start", "end", "rank"]
    # the written site coordinates are genomic, i.e. peak start + offset
    if len(rows) > 1:
        fields = rows[1].split("\t")
        assert int(fields[11]) == 90 + 10


def test_short_motif_is_skipped_and_recorded_not_silently_empty():
    """A 4 bp motif cannot reach fpr=1e-4; it must be reported as skipped."""
    ann = annotate_sequences(["GATAGATAGATA"], [_GATA], fpr=1e-4)
    assert ann["summary"] == []
    assert "GATA" in ann["skipped_motifs"]
    rec = ann["skipped_motifs"]["GATA"]
    assert rec["best_achievable_fpr"] == pytest.approx(0.25 ** 4, rel=1e-6)
    assert rec["length"] == 4


def test_cli_reports_skipped_motifs_in_the_json(tmp_path, monkeypatch):
    from easyepi import tf_annotate

    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\n" + "N" * 400 + "\n")
    motif_file = tmp_path / "m.pfm"
    motif_file.write_text(">MA0001.1 GATA\n0 4 0 4\n0 0 0 0\n4 0 0 0\n0 0 4 0\n",
                          encoding="utf-8")
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t0\t100\tpeakA\n", encoding="utf-8")
    monkeypatch.setattr(tf_annotate, "fetch_sequences",
                        lambda records, fasta: ["N" * 100])

    class Args:
        peaks = str(bed)
        genome_fasta = str(genome)
        motifs = str(motif_file)
        fpr = 1e-4
        background = None
        top = 5
        out = str(tmp_path / "out")
        figure_formats = None

    assert tf_annotate.run_tf_annotate(Args) == 0
    payload = json.loads((tmp_path / "out" / "tf_annotate.json").read_text())
    assert "GATA" in payload["skipped_motifs"]
    assert payload["summary"] == []
    # the summary file still exists with only its header
    assert (tmp_path / "out" / "tf_summary.tsv").read_text().strip().count("\n") == 0


def test_cli_errors_when_every_sequence_is_empty(tmp_path, monkeypatch):
    from easyepi import tf_annotate

    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    motif_file = tmp_path / "m.pfm"
    motif_file.write_text(">MA0001.1 GATA\n1 0 0 0\n0 0 0 1\n0 0 0 1\n0 1 0 0\n",
                          encoding="utf-8")
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t0\t10\tpeakA\n", encoding="utf-8")
    monkeypatch.setattr(tf_annotate, "fetch_sequences",
                        lambda records, fasta: [""])

    class Args:
        peaks = str(bed)
        genome_fasta = str(genome)
        motifs = str(motif_file)
        fpr = 0.3
        background = None
        top = 5
        out = str(tmp_path / "out")
        figure_formats = None

    with pytest.raises(ValueError, match="every peak sequence is empty"):
        tf_annotate.run_tf_annotate(Args)
