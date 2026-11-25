"""TOBIAS's footprint score is easyepi's score on a different scale.

The two tools answer the same question with the same two quantities, but publish
them in different units. Read from TOBIAS 0.17.5's own Cython source
(`tobias/utils/signals.pyx`, `tobias_footprint_array`):

    fp_sum     = sum of the NEGATIVE values inside the footprint window
    fp_mean    = fp_sum / footprint_w                 (<= 0)
    flank_mean = (left_sum + right_sum) / (2*flank_w)
    footprint  = flank_mean - fp_mean                 additive, raw-signal units

`easyepi.footprint.footprint_profile` reports

    score = (shoulder_baseline - centre) / shoulder_baseline      fractional

so with the same windows the two are related by exactly

    TOBIAS == flank_level * easyepi

This test pins that relation on a synthetic ATACorrect-style track (positive
flanks, negative footprint), which is the case TOBIAS's additive score assumes:
its `fp_sum` accumulates only negative values, so a purely positive dip makes it
report the flank level unchanged.

Measured ratio on seven sites: 1.000000 (min = max = mean).
"""
from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pyBigWig
import pytest

from easyepi.footprint import footprint_profile

FLANK_LEVEL = 10.0
DEPTH = -5.0
FLANK_MIN, FLANK_MAX = 15, 25
FP_MIN, FP_MAX = 8, 12


def _tobias():
    from pathlib import Path
    import sys

    for cand in (Path(sys.executable).parent / "TOBIAS", shutil.which("TOBIAS")):
        if cand and Path(cand).exists():
            return str(cand)
    return None


TOBIAS = _tobias()
requires_tobias = pytest.mark.skipif(TOBIAS is None, reason="TOBIAS not installed")


def _track(tmp_path, centres):
    n = 40_000
    values = np.full(n, FLANK_LEVEL)
    for c in centres:
        values[c - 5:c + 6] = DEPTH
    bw_path = tmp_path / "cuts.bw"
    bw = pyBigWig.open(str(bw_path), "w")
    bw.addHeader([("chr1", n)])
    bw.addEntries(["chr1"] * n, [int(i) for i in range(n)],
                  ends=[int(i + 1) for i in range(n)],
                  values=[float(v) for v in values])
    bw.close()
    return bw_path, np.array(values, dtype=float)


def _score_with_tobias(tmp_path, bw_path, centres):
    bed = tmp_path / "sites.bed"
    bed.write_text("".join("chr1\t%d\t%d\ts%d\n" % (c - 200, c + 200, i)
                           for i, c in enumerate(centres)), encoding="utf-8")
    out = tmp_path / "fp.bw"
    cmd = [TOBIAS, "ScoreBigwig", "--signal", str(bw_path),
           "--regions", str(bed), "--output", str(out), "--score", "footprint",
           "--flank-min", str(FLANK_MIN), "--flank-max", str(FLANK_MAX),
           "--fp-min", str(FP_MIN), "--fp-max", str(FP_MAX)]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert proc.returncode == 0, proc.stdout.decode()[-1500:]
    assert out.exists(), "TOBIAS produced no bigWig"

    bw = pyBigWig.open(str(out))
    vals = [bw.values("chr1", c, c + 1)[0] for c in centres]
    bw.close()
    return vals


@requires_tobias
def test_tobias_score_equals_flank_level_times_easyepi(tmp_path):
    centres = [4000, 9000, 14000, 19000, 24000, 29000, 34000]
    bw_path, arr = _track(tmp_path, centres)
    tobias_vals = _score_with_tobias(tmp_path, bw_path, centres)

    for c, tob in zip(centres, tobias_vals):
        prof = footprint_profile(arr, [c], center=0, flank=100)
        easy = prof["score"]
        bg = prof["background_profile"]
        base = float(bg[len(bg) // 2])
        assert tob == pytest.approx(base * easy, rel=1e-6, abs=1e-9), (
            "site %d: TOBIAS %.6f vs flank(%.4f) * easyepi(%.6f) = %.6f"
            % (c, tob, base, easy, base * easy))


@requires_tobias
def test_purely_positive_dip_leaves_tobias_at_the_flank_level(tmp_path):
    """TOBIAS accumulates only negative values, so this must not change.

    This is a property of TOBIAS's definition, not a bug in either tool: it is
    why a comparison has to use an ATACorrect-style (negative) footprint.
    """
    centres = [5000, 15000]
    n = 20_000
    values = np.full(n, 10.0)
    for c in centres:
        values[c - 5:c + 6] = 3.0          # positive dip, above zero
    bw_path = tmp_path / "pos.bw"
    bw = pyBigWig.open(str(bw_path), "w")
    bw.addHeader([("chr1", n)])
    bw.addEntries(["chr1"] * n, [int(i) for i in range(n)],
                  ends=[int(i + 1) for i in range(n)],
                  values=[float(v) for v in values])
    bw.close()

    vals = _score_with_tobias(tmp_path, bw_path, centres)
    for c, v in zip(centres, vals):
        assert v == pytest.approx(10.0, abs=1e-6), \
            "site %d: expected the flank level, got %.4f" % (c, v)


def test_easyepi_score_is_zero_when_there_is_no_dip():
    """The denominator guard: a flat track must give 0.0, not a division error."""
    arr = np.full(5_000, 7.0)
    prof = footprint_profile(arr, [2500], center=0, flank=50)
    assert prof["score"] == pytest.approx(0.0, abs=1e-12)
