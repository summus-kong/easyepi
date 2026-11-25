"""`easyepi.motifs.scan_sequences` finds the same best site as FIMO.

The harness that prints the full table is `tests/oracle/compare_to_fimo.py`; this
module asserts the properties that must hold every run. Skips cleanly when
`fimo` is absent, matching how `tests/oracle/` handles R.

What is asserted is the best window and the ranking, because those are what a
motif call depends on. The absolute scores are deliberately **not** compared:
measured on four fixed 7-mers they differ by an amount whose sign changes with
base composition (GATATCA +0.40, GATATCG +0.11, AAAACGT +0.26, GGGGGGG -1.14),
and a constant offset, a pseudocount and background smoothing were each ruled
out. FIMO's `score` column is its own statistic, not `read_pwm`'s log2-odds sum.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests" / "oracle"))

from easyepi.motifs import read_pwm, scan_sequences  # noqa: E402


def _harness():
    """Import the oracle harness for its fixture and writers."""
    import importlib.util

    path = REPO / "tests" / "oracle" / "compare_to_fimo.py"
    spec = importlib.util.spec_from_file_location("cmp_fimo", str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _fimo():
    cand = Path(sys.executable).parent / "fimo"
    if cand.exists():
        return str(cand)
    found = shutil.which("fimo")
    if found:
        return found
    pytest.skip("fimo not installed; the FIMO oracle needs MEME Suite")


@pytest.fixture(scope="module")
def compared(tmp_path_factory):
    _fimo()                            # skip before doing any work
    h = _harness()
    tmp = tmp_path_factory.mktemp("fimo")
    labels, seqs = h.fixture_sequences()
    meme, jaspar, fasta = h.write_inputs(tmp, labels, seqs)

    pwm = read_pwm(jaspar)
    our = scan_sequences(pwm, seqs)
    rows = h.run_fimo(meme, fasta, thresh=1.0)
    best = h.best_per_sequence(rows)
    records = h.compare(labels, seqs, our, best)
    return h, pwm, [r for r in records if r["had_hit"]]


def test_the_matrix_is_not_truncated(compared):
    """A 4-row JASPAR file for a 7-column motif must stay 7 wide.

    The first version of the harness transposed COUNTS the wrong way and every
    comparison silently ran against a 4 bp motif.
    """
    _h, pwm, _records = compared
    assert pwm.shape == (4, 7), pwm.shape


def test_the_best_window_matches_for_planted_sites(compared):
    """Where a consensus was planted, both tools must pick the same bp and strand."""
    _h, _pwm, records = compared
    planted = [r for r in records if r["label"] in ("seq0", "seq1", "seq3",
                                                    "seq4", "seq6", "seq7",
                                                    "seq9", "seq10")]
    assert planted, "no planted-site records"
    for r in planted:
        assert r["our_pos"] == r["fimo_pos"], (
            "%s: easyepi %d vs FIMO %d" % (r["label"], r["our_pos"], r["fimo_pos"]))
        assert r["our_strand"] == r["fimo_strand"], (
            "%s strand: easyepi %s vs FIMO %s"
            % (r["label"], r["our_strand"], r["fimo_strand"]))


def test_the_two_rank_the_sequences_identically(compared):
    """Ranking is what a threshold-free comparison actually uses."""
    _h, _pwm, records = compared
    if len(records) < 3:
        pytest.skip("not enough records to rank")
    a = np.array([r["our_score"] for r in records])
    b = np.array([r["fimo_score"] for r in records])
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    rho = float(np.corrcoef(ra, rb)[0, 1])
    assert rho == pytest.approx(1.0), "ranking disagrees, rho=%.3f" % rho


def test_at_most_one_sequence_disagrees(compared):
    """One noise-only sequence may tie differently; more would be a real defect."""
    _h, _pwm, records = compared
    mismatched = [r["label"] for r in records
                  if r["our_pos"] != r["fimo_pos"]
                  or r["our_strand"] != r["fimo_strand"]]
    assert len(mismatched) <= 1, "disagreed on: %s" % mismatched
