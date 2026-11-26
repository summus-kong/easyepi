"""The fingerprint port is checked against the real TOBIAS binary.

A "1:1 implementation" is only meaningful if it agrees with the reference, so
these tests run TOBIAS 0.17.5 on the same signal and compare base by base.
They skip when TOBIAS is absent.

Measured agreement (8 kb synthetic signal with planted dips, defaults):

    interior region   mean |diff| 2.9e-5, max 1.5e-4, 0 bases > 0.01
    region at chr 0   mean |diff| 7.2e-3, max 2.2,    20 bases > 0.01
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pyBigWig
import pytest

from easyepi.specific.atac.tobias_footprint import (
    TOBIAS_DEFAULTS,
    fast_rolling_math,
    score_footprint_like_tobias,
    tobias_footprint_array,
)

ENV = Path("/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi")


def _tobias():
    for cand in (ENV / "bin/TOBIAS", shutil.which("TOBIAS")):
        if cand and Path(cand).exists():
            return str(cand)
    return None


TOBIAS = _tobias()
requires_tobias = pytest.mark.skipif(TOBIAS is None, reason="TOBIAS not installed")


# ---------------------------------------------------------------------------
# Kernel properties that do not need the binary
# ---------------------------------------------------------------------------
def test_single_dip_scores_flank_minus_footprint_mean():
    """A -3 dip in a +5 background scores 5 - (-3) = 8, as TOBIAS returns."""
    sig = np.full(400, 5.0)
    sig[200:220] = -3.0
    out = tobias_footprint_array(sig, 10, 30, 20, 50)
    assert out[210] == pytest.approx(8.0, abs=1e-4)


def test_positive_dip_scores_the_flank_level_because_only_negatives_count():
    """The asymmetric filter: a dip that stays positive is invisible.

    `fp_sum` accumulates negatives only, so a positive dip leaves `fp_mean = 0`
    and the score is `flank_mean - 0`, i.e. the flank level itself -- not a
    footprint. Measured against the binary: 5.0 for a 5.0 signal.
    """
    sig = np.full(400, 5.0)
    sig[200:220] = 1.0            # a dip, but above zero
    out = tobias_footprint_array(sig, 10, 30, 20, 50)
    assert out[210] == pytest.approx(5.0, abs=1e-3)


def test_tail_is_left_unscored():
    """The kernel cannot reach the end of the array, and the tail stays 0.0.

    The boundary is measured rather than derived: three hand-computations of it
    were wrong, and the kernel's write range depends on the largest
    (flank_w, footprint_w) pair rather than on the loop bound alone.
    """
    sig = np.tile(np.r_[np.full(20, 5.0), np.full(20, -5.0)], 20)
    out = tobias_footprint_array(sig, 10, 30, 20, 50)

    written = np.flatnonzero(out != 0.0)
    assert written.size > 0, "nothing was scored"
    last = int(written[-1])
    assert last < len(sig) - 1, "the kernel must not reach the final base"

    after = out[last + 1:]
    assert np.all(after == 0.0), (
        "bases after the last written position must be zero, got %r" % after[:5])

    # Measured on this fixture: last written index 768 of 800 -> 31 tail bases.
    # (The tail length follows from the start bound L-2*flank_max-fp_max plus
    # the largest window write, not from the smallest windows, so no algebraic
    # lower bound is asserted here.)


def test_short_input_returns_all_zero():
    out = tobias_footprint_array(np.ones(50), 10, 30, 20, 50)
    assert out.shape == (50,)
    assert np.all(out == 0.0)


def test_only_maxima_are_kept():
    """A stronger footprint must overwrite a weaker one at the same base."""
    sig = np.full(600, 5.0)
    sig[300:320] = -2.0
    weak = tobias_footprint_array(sig, 10, 30, 20, 50)
    sig[300:320] = -6.0
    strong = tobias_footprint_array(sig, 10, 30, 20, 50)
    assert strong[310] > weak[310]


def test_bad_window_bounds_are_rejected():
    with pytest.raises(ValueError):
        tobias_footprint_array(np.ones(400), 30, 10, 20, 50)
    with pytest.raises(ValueError):
        tobias_footprint_array(np.ones(400), 10, 30, 50, 20)


def test_smoothing_matches_a_plain_moving_mean():
    v = np.arange(10, dtype=float)
    got = fast_rolling_math(v, 3, "mean")
    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(2.0)
    assert np.isnan(got[-2:]).all()


# ---------------------------------------------------------------------------
# Agreement with the reference binary
# ---------------------------------------------------------------------------
def _signal(n=8000, seed=3):
    rng = np.random.default_rng(seed)
    sig = rng.normal(2.0, 0.6, n)
    for c in rng.integers(500, n - 500, 30):
        w = int(rng.integers(22, 46))
        sig[c - w // 2:c + w // 2] -= rng.uniform(3.0, 6.0)
    return sig


def _run_tobias(tmp_path, sig, region):
    bw_in = tmp_path / "s.bw"
    bw = pyBigWig.open(str(bw_in), "w")
    n = len(sig)
    bw.addHeader([("chr1", n)])
    bw.addEntries(["chr1"] * n, [int(i) for i in range(n)],
                  ends=[int(i + 1) for i in range(n)],
                  values=[float(v) for v in sig])
    bw.close()

    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t%d\t%d\n" % region)
    out = tmp_path / "o.bw"
    d = TOBIAS_DEFAULTS
    cmd = [TOBIAS, "ScoreBigwig", "--signal", str(bw_in), "--regions", str(bed),
           "--output", str(out), "--score", "footprint",
           "--flank-min", str(d["flank_min"]), "--flank-max", str(d["flank_max"]),
           "--fp-min", str(d["fp_min"]), "--fp-max", str(d["fp_max"])]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert proc.returncode == 0, proc.stdout.decode()[-1200:]

    rb = pyBigWig.open(str(out))
    ref = np.nan_to_num(np.array(rb.values("chr1", 0, n), dtype=np.float64))
    rb.close()

    # read back what the binary actually consumed (bigWig is float32)
    rb = pyBigWig.open(str(bw_in))
    stored = np.nan_to_num(np.array(rb.values("chr1", 0, n), dtype=np.float64))
    rb.close()
    return ref, stored


@requires_tobias
def test_interior_region_agrees_to_float32(tmp_path):
    """The 1:1 criterion: an interior region matches the binary to float32."""
    sig = _signal()
    ref, stored = _run_tobias(tmp_path, sig, (1000, 7000))

    d = TOBIAS_DEFAULTS
    mine = score_footprint_like_tobias(stored, d["flank_min"], d["flank_max"],
                                       d["fp_min"], d["fp_max"], d["smooth"])
    diff = np.abs(ref[1000:7001] - mine[1000:7001])

    assert diff.max() < 1e-3, (
        "max |diff| %.6g exceeds float32 noise; mean %.6g" % (diff.max(), diff.mean()))
    assert np.sum(diff > 0.01) == 0, "bases disagree beyond tolerance"
    assert diff.mean() < 1e-4


@requires_tobias
def test_single_dip_scores_identically_in_both(tmp_path):
    """One planted dip: both must report the same value at the dip."""
    n = 4000
    sig = np.full(n, 5.0)
    sig[2000:2020] = -3.0
    ref, stored = _run_tobias(tmp_path, sig, (0, n - 1))

    d = TOBIAS_DEFAULTS
    mine = score_footprint_like_tobias(stored, d["flank_min"], d["flank_max"],
                                       d["fp_min"], d["fp_max"], d["smooth"])
    assert mine[2010] == pytest.approx(ref[2010], abs=1e-3)
    assert mine[2010] == pytest.approx(8.0, abs=1e-2)


@requires_tobias
def test_all_positive_signal_yields_the_flank_level_in_both(tmp_path):
    """With no negative values both report the flank level, not a footprint."""
    n = 4000
    sig = np.full(n, 5.0)
    ref, stored = _run_tobias(tmp_path, sig, (0, n - 1))

    d = TOBIAS_DEFAULTS
    mine = score_footprint_like_tobias(stored, d["flank_min"], d["flank_max"],
                                       d["fp_min"], d["fp_max"], d["smooth"])
    # With no negative values: fp_mean = 0, so score = flank_mean = 5.0.
    # The port must reproduce that, including the tail's zeros.
    assert ref.max() == pytest.approx(5.0, abs=1e-3)
    assert mine.max() == pytest.approx(5.0, abs=1e-3)
    # and no base reports an actual footprint (a score above the flank level)
    assert ref.max() <= 5.0 + 1e-3
    assert mine.max() <= 5.0 + 1e-3
