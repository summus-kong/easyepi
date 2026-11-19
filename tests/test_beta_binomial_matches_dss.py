"""Check the pure-Python estimator against DSS, when R happens to be available.

The product no longer needs R: ``easyepi.dmr_beta_binomial`` reproduces
``DSS::DMLtest`` in numpy/scipy. That claim is only worth making if it is
checked, and checked at every step: a matching p-value can hide a wrong
dispersion that happens to cancel, so ``tests/oracle/dss_oracle.R`` emits the
intermediates (priors, both smoothed levels, phi, both variances, se, diff,
pval) and each column is compared separately.

The oracle also asserts internally that its own hand-rolled chain agrees with
``DMLtest`` to 1e-12, so a green Python comparison really is a comparison with
DSS and not just with the test's idea of it.

These tests skip cleanly when ``Rscript`` or ``DSS`` is missing, which is the
normal case for an installed package -- the point of the rewrite is that the
package works without them.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from easyepi.dmr_beta_binomial import beta_binomial_test
from easyepi.dmr_r import find_rscript

ORACLE = Path(__file__).parent / "oracle" / "dss_oracle.R"

#: Sites per chromosome. Enough for a prior (>50 usable) and for windows that
#: hold several sites, small enough that DSS stays a few seconds.
PER_CHR = 150
N_A = 3
N_B = 3
#: First this many sites per chromosome carry a real difference.
N_DIFF = 40


def _fixture(n_a: int = N_A, n_b: int = N_B):
    """Deterministic, irregularly spaced, multi-chromosome counts.

    Three properties are deliberate:

    * irregular spacing, so an index-based window cannot pass;
    * two chromosomes, so the smoothing must not cross the boundary;
    * scatter *within* each group, so the dispersion is identifiable instead
      of being pushed to the upper bound.
    """
    rng = np.random.default_rng(11)
    n_sites = 2 * PER_CHR
    chrom = np.array(["chr1"] * PER_CHR + ["chr2"] * PER_CHR)
    gaps = rng.integers(1, 120, n_sites)
    pos = np.concatenate([
        np.cumsum(gaps[:PER_CHR]) + 1000,
        np.cumsum(gaps[PER_CHR:]) + 5000,
    ]).astype(np.int64)

    n_samples = n_a + n_b
    total = rng.integers(12, 45, (n_samples, n_sites)).astype(np.int64)
    base = np.empty((n_samples, n_sites), dtype=float)
    base[:] = 0.30
    base[n_a:, :N_DIFF] = 0.65
    base[n_a:, PER_CHR:PER_CHR + N_DIFF] = 0.65
    # Between-replicate scatter: without it every phi is pinned to the bound.
    level = np.clip(rng.normal(base, 0.06), 0.01, 0.99)
    meth = rng.binomial(total, level).astype(np.int64)
    design = np.concatenate([np.zeros(n_a, dtype=int), np.ones(n_b, dtype=int)])
    return meth, total, design, chrom, pos


def _write_payload(meth, total, design, chrom, pos, path: Path) -> None:
    n_sites = meth.shape[1]
    n_samples = meth.shape[0]
    payload = "\n".join([
        str(n_sites),
        str(n_samples),
        " ".join(str(c) for c in chrom),
        " ".join(str(int(v)) for v in pos),
        " ".join(str(int(v)) for v in meth.ravel()),
        " ".join(str(int(v)) for v in total.ravel()),
        " ".join(str(int(v)) for v in design),
    ]) + "\n"
    path.write_text(payload, encoding="utf-8")


def _run_oracle(tmp: Path, smoothing: bool = True, equal_disp: bool = False,
                direct: bool = False, **fixture_args) -> pd.DataFrame:
    meth, total, design, chrom, pos = _fixture(**fixture_args)
    payload = tmp / "counts.txt"
    out = tmp / "dss.tsv"
    _write_payload(meth, total, design, chrom, pos, payload)
    proc = subprocess.run(
        [find_rscript(), "--vanilla", str(ORACLE), str(payload), str(out),
         "TRUE" if smoothing else "FALSE",
         "TRUE" if equal_disp else "FALSE",
         "TRUE" if direct else "FALSE"],
        capture_output=True, text=True, timeout=3600,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    return pd.read_csv(out, sep="\t")


requires_r = pytest.mark.skipif(
    find_rscript() is None or not ORACLE.exists(),
    reason="Rscript with DSS not available")


def _compare(smoothing: bool, equal_disp: bool, direct: bool = False,
             **fixture_args) -> None:
    """Compare every intermediate against DSS for one (smoothing, equal_disp) mode."""
    tmp_ctx = tempfile.TemporaryDirectory(prefix="easyepi-dss-oracle-")
    with tmp_ctx as tmp:
        ref = _run_oracle(Path(tmp), smoothing=smoothing,
                          equal_disp=equal_disp, direct=direct, **fixture_args)

    meth, total, design, chrom, pos = _fixture(**fixture_args)
    mine = beta_binomial_test(meth, total, design, smoothing_span=500,
                              chrom=chrom, pos=pos, equal_disp=equal_disp,
                              smoothing=smoothing)
    tested = np.asarray(mine["tested"], dtype=bool)
    assert tested.sum() == len(ref), (
        f"site filter disagrees with DSS: {tested.sum()} kept, DSS kept {len(ref)}")
    assert np.array_equal(np.flatnonzero(tested), ref["rn"].to_numpy(int) - 1)

    # DSS returns its sites in (chr, pos) order, which is the fixture order
    # here only if the two chromosomes are contiguous and sorted; sort both
    # sides on the returned coordinates to be independent of that.
    ref = ref.sort_values(["chr", "pos"]).reset_index(drop=True)
    order = np.lexsort((pos[tested], chrom[tested].astype(str)))
    mine_idx = np.flatnonzero(tested)[order]
    got = {
        "e1": mine["level1"][mine_idx],
        "e2": mine["level2"][mine_idx],
        "e1w": mine["level1_weighted"][mine_idx],
        "e2w": mine["level2_weighted"][mine_idx],
        "phi1": mine["phi1"][mine_idx],
        "phi2": mine["phi2"][mine_idx],
        "var1": mine["var1"][mine_idx],
        "var2": mine["var2"][mine_idx],
        "diff": mine["diff"][mine_idx],
        "se": mine["se"][mine_idx],
        "stat": mine["stat"][mine_idx],
        "pval": mine["pvalue"][mine_idx],
    }

    # Smoothed levels are exact transcriptions of the same C kernel.
    for key in ("e1", "e2", "e1w", "e2w"):
        np.testing.assert_allclose(got[key], ref[key].to_numpy(float),
                                   rtol=1e-10, atol=1e-12, err_msg=key)

    # The dispersion is optimised, not solved; DSS's own tolerance is 1e-3.
    for key in ("phi1", "phi2"):
        np.testing.assert_allclose(got[key], ref[key].to_numpy(float),
                                   rtol=2e-3, atol=1e-6, err_msg=key)

    for key in ("var1", "var2", "diff", "se", "stat"):
        np.testing.assert_allclose(got[key], ref[key].to_numpy(float),
                                   rtol=1e-6, atol=1e-9, err_msg=key)

    np.testing.assert_allclose(got["pval"], ref["pval"].to_numpy(float),
                               rtol=1e-6, atol=1e-30, err_msg="pval")

    assert mine["dispersion_method"] == "dss_logN_shrinkage"
    assert "DSS 2.58.0" in mine["reference"]
    assert mine["smoothing"] is smoothing
    # `equal_disp` is what was *used*, not what was asked for: with
    # smoothing=FALSE and a group of one replicate, DSS shares a dispersion on
    # its own, so the reported flag is True even though the caller did not ask.
    forced_share = (not smoothing) and (fixture_args.get("n_a", N_A) == 1
                                        or fixture_args.get("n_b", N_B) == 1)
    assert mine["equal_disp"] is (equal_disp or forced_share)
    # And the same sites are called.
    dss_p = ref["pval"].to_numpy(float)
    assert ((got["pval"] < 0.05) == (dss_p < 0.05)).all()


@requires_r
def test_default_mode_reproduces_dss():
    """``DMLtest(smoothing = TRUE, equal.disp = FALSE)``, every intermediate.

    The tolerances are the measured ones, not a guess: the dispersion is a
    bounded 1-D optimisation whose tolerance DSS fixes at ``tol = 0.001`` on
    ``log(phi)``, so two correct implementations can differ there by that much
    and no further.
    """
    _compare(smoothing=True, equal_disp=False)


@requires_r
def test_equal_dispersion_matches_dss():
    """``DMLtest(equal.disp = TRUE)``: one dispersion, compared not just checked."""
    _compare(smoothing=True, equal_disp=True)


@requires_r
def test_no_smoothing_matches_dss():
    """``DMLtest(smoothing = FALSE)``: a separate path with its own variance floor."""
    _compare(smoothing=False, equal_disp=False)


@requires_r
def test_no_smoothing_equal_disp_matches_dss():
    """``smoothing = FALSE`` with ``equal.disp = TRUE``."""
    _compare(smoothing=False, equal_disp=True)


@requires_r
def test_no_smoothing_one_replicate_matches_dss():
    """One replicate per group with ``smoothing = FALSE``.

    ``DMLtest`` refuses this combination outright ("There is no biological
    replicates in at least one condition"), but ``DSS:::DMLtest.noSmooth``
    handles it by sharing one dispersion. The oracle therefore calls the
    internal function, and this test pins that Python makes the same choice.
    """
    _compare(smoothing=False, equal_disp=False, direct=True, n_a=1, n_b=3)


def test_transcription_constants_are_the_documented_ones():
    """The numbers DSS hard-codes, pinned so a refactor cannot drift silently."""
    from easyepi import dmr_beta_binomial as m

    assert m._RHO == 0.8                      # compute.var.Smooth
    assert m._COVERAGE_SPAN == 200            # hasCoverage: ws <- 200
    assert (m._VAR_LOW_TRIGGER, m._VAR_LOW_SET) == (1e-5, 1e-4)
    assert (m._VAR_HIGH_TRIGGER, m._VAR_HIGH_SET) == (1 - 1e-5, 1 - 1e-4)
    # The no-smoothing path floors to 1e-05, one-sided. Same trigger, different
    # replacement, and no cap: two guards, deliberately not shared.
    assert m._NO_SMOOTH_VAR_FLOOR == 1e-5
    assert m._NO_SMOOTH_VAR_FLOOR != m._VAR_LOW_SET
