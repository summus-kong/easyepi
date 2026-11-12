"""The PyDESeq2 adapter must get the sign, the refusals and the reporting right.

The sign check is the important one. PyDESeq2's `contrast=[factor, first, second]`
reports log2(second / first) -- the opposite of R DESeq2's `resultsNames` convention --
so an adapter that passes the operands in the natural reading order reports every
fold-change upside down. P-values stay identical, which means the mistake survives any
check that only looks at significance. It was found here by building data where the
direction is known, not by reading the API.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pydeseq2")

from easyepi.differential.deseq2 import deseq2_available, run_deseq2  # noqa: E402


def _neutral_counts(seed=4, n=40, per=4):
    """Counts where a quarter go up, a quarter go down, so size factors stay ~1.

    Without that balance a one-sided change shifts the median-of-ratios size factors
    and the expected log2FC is no longer +/-1 -- which makes the sign check unable to
    distinguish a wrong sign from an absorbed normalisation shift.
    """
    rng = np.random.default_rng(seed)
    samples = [f"A{i}" for i in range(per)] + [f"B{i}" for i in range(per)]
    design = {s: {"batch": "b1" if int(s[1]) % 2 == 0 else "b2"} for s in samples}
    mask = np.array([True] * per + [False] * per)
    base = rng.integers(200, 400, n)
    counts = np.zeros((n, len(samples)), dtype=int)
    for j, s in enumerate(samples):
        mu = base.astype(float).copy()
        if s.startswith("B"):
            mu[0:10] *= 2.0          # up
            mu[10:20] *= 0.5         # down
        counts[:, j] = rng.poisson(mu)
    return counts, samples, design, mask


def test_log2fc_points_from_denominator_to_numerator():
    """A feature doubled in group B must report a positive log2FC for 'A vs B'."""
    counts, samples, design, mask = _neutral_counts()
    res = run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                     design, mask, ("A", "B"), quiet=True)
    lfc = res["log2FoldChange"]

    assert np.nanmedian(lfc[0:10]) == pytest.approx(1.0, abs=0.25), "up set sign/shape"
    assert np.nanmedian(lfc[10:20]) == pytest.approx(-1.0, abs=0.25), "down set sign/shape"
    assert abs(float(np.nanmedian(lfc[20:40]))) < 0.25, "null set not centred on zero"
    assert res["info"]["contrast_reported_as"] == "log2(A / B)"


def test_signal_is_recovered_and_null_is_not():
    counts, samples, design, mask = _neutral_counts()
    res = run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                     design, mask, ("A", "B"), quiet=True)
    sig = np.nan_to_num(res["pvalue"], nan=1.0) < 0.05
    assert sig[0:10].sum() >= 8, f"only {sig[0:10].sum()}/10 up features called"
    assert sig[10:20].sum() >= 8, f"only {sig[10:20].sum()}/10 down features called"
    assert sig[20:40].sum() <= 2, f"{sig[20:40].sum()}/20 null features called"


def test_model_provenance_is_reported():
    """A model is only reproducible with the version and design that produced it."""
    counts, samples, design, mask = _neutral_counts()
    res = run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                     design, mask, ("A", "B"), quiet=True)
    info = res["info"]
    assert info["pydeseq2_version"]
    assert info["design"].startswith("~")
    assert info["contrast"] == ["easyepi_group", "A", "B"]
    assert info["replicates_a"] == info["replicates_b"] == 4
    assert "treated_as_na_before_filtering" in info


def test_a_batch_term_is_actually_modelled():
    counts, samples, design, mask = _neutral_counts()
    res = run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                     design, mask, ("A", "B"), design_factors=("batch",), quiet=True)
    assert "batch" in res["info"]["design"]


def test_a_constant_batch_is_refused():
    """A factor with one level cannot be estimated; silently dropping it hides a mistake."""
    counts, samples, design, mask = _neutral_counts()
    flat = {s: {"batch": "same"} for s in samples}
    with pytest.raises(ValueError, match="single level"):
        run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                   flat, mask, ("A", "B"), design_factors=("batch",), quiet=True)


def test_a_missing_factor_value_is_refused():
    counts, samples, design, mask = _neutral_counts()
    broken = dict(design)
    broken[samples[0]] = {}
    with pytest.raises(ValueError, match="missing for sample"):
        run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                   broken, mask, ("A", "B"), design_factors=("batch",), quiet=True)


@pytest.mark.parametrize("bad, message", [
    (lambda c: c + 0.5, "integer counts"),
    (lambda c: c - 10_000, "negative values"),
    (lambda c: np.where(c > 500, np.nan, c), "NaN or inf"),
])
def test_inputs_deseq2_cannot_model_are_refused(bad, message):
    counts, samples, design, mask = _neutral_counts()
    with pytest.raises(ValueError, match=message):
        run_deseq2([f"f{i}" for i in range(counts.shape[0])], bad(counts), samples,
                   design, mask, ("A", "B"), quiet=True)


def test_too_few_samples_is_refused():
    counts, samples, design, mask = _neutral_counts()
    with pytest.raises(ValueError, match="at least 4 samples"):
        run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts[:, :3],
                   samples[:3], design, np.array([True, False, True]), ("A", "B"),
                   quiet=True)


def test_one_replicate_per_group_is_refused():
    counts, samples, design, mask = _neutral_counts()
    unbalanced = np.array([True] + [False] * (counts.shape[1] - 1))
    with pytest.raises(ValueError, match="at least 2 replicates"):
        run_deseq2([f"f{i}" for i in range(counts.shape[0])], counts, samples,
                   design, unbalanced, ("A", "B"), quiet=True)


def test_availability_probe_does_not_import():
    assert deseq2_available() is True
