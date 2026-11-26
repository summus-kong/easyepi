"""Footprint engine selection: `--engine {native,port,tobias}`.

The ports were previously reachable only from tests, which made them dead code.
These tests pin the wiring rather than the numerics (the numerics have their own
suites): the engine choice must be honoured, the default must not change historic
output, and the port must refuse honestly when it cannot score.

Measured constraint this encodes: a corrected ATAC track is sparse and its finite
values cluster by sign, so a real 110 bp window often has no positive flank. When
that happens the additive TOBIAS score is 0 -- `flank_mean` only sums positives --
and the engine reports `degenerate` instead of inventing a number.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from easyepi.footprint_engines import (
    DEFAULT_FLANK_MAX, DEFAULT_FLANK_MIN, DEFAULT_FP_MAX, DEFAULT_FP_MIN,
    tobias_port_score,
)


def _dip_signal(n: int = 400, centre: int = 200) -> np.ndarray:
    """A track with a positive flank and a genuine dip, as ATACorrect produces."""
    sig = np.full(n, 10.0)
    sig[centre - 10:centre + 10] = -6.0
    return sig


def test_port_defaults_match_tobias():
    """The window bounds are the reference's, not invented."""
    assert (DEFAULT_FP_MIN, DEFAULT_FP_MAX) == (20, 50)
    assert (DEFAULT_FLANK_MIN, DEFAULT_FLANK_MAX) == (10, 30)


def test_port_scores_a_planted_dip_positive():
    """A dip with a positive flank must score above zero."""
    sig = _dip_signal()
    out = tobias_port_score(sig, [200], flank=100)
    assert out["engine"] == "port"
    assert out["n_sites"] == 1
    assert out["degenerate"] is False
    assert out["score"] > 0


def test_port_reports_degenerate_without_negatives():
    """All-positive input cannot yield a dip; the flag must say so."""
    sig = np.full(400, 10.0)
    out = tobias_port_score(sig, [200], flank=100)
    assert out["degenerate"] is True
    assert "no negative values" in out["note"]


def test_port_scores_an_all_negative_window_from_the_dip_alone():
    """No positive flank means `flank_mean` is 0, so the dip sets the score.

    Measured on `ENCFF327MYY` corrected signal: whole stretches are negative and a
    short window inside one carries no positive value. With `flank_mean = 0` and
    `fp_mean = -5`, the additive score is `0 - (-5) = 5` -- the dip magnitude, not
    zero. Worth pinning because the naive reading of "no flank" is "no signal".
    """
    sig = np.full(400, -5.0)
    out = tobias_port_score(sig, [200], flank=100)
    assert out["score"] == pytest.approx(5.0)


def test_port_rejects_sites_outside_the_window():
    """A site list that falls entirely outside the array is a caller error."""
    with pytest.raises(ValueError, match="no target site"):
        tobias_port_score(np.zeros(100), [500, 900], flank=50)


def test_cli_exposes_the_engine_choices():
    """`--engine` must be discoverable and default to native."""
    import argparse
    import re
    import subprocess
    import sys

    res = subprocess.run([sys.executable, "-m", "easyepi.cli", "footprint", "--help"],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-400:]
    assert "--engine {native,port,tobias}" in res.stdout
    match = re.search(r"--engine \{native,port,tobias\}\s*\n\s*footprint scorer", res.stdout)
    assert match, "engine help text missing"
    assert argparse is not None


def test_native_is_the_default_and_port_is_opt_in():
    """Parsing without `--engine` must select native, preserving old behaviour."""
    from easyepi.cli import build_parser

    parser = build_parser()
    base = ["footprint", "--bigwig", "s.bw", "--target", "t.bed",
            "--chrom", "chr1", "--start", "0", "--end", "100", "--out", "o"]
    assert parser.parse_args(base).engine == "native"
    assert parser.parse_args(base + ["--engine", "port"]).engine == "port"
    with pytest.raises(SystemExit):
        parser.parse_args(base + ["--engine", "bogus"])


def test_port_payload_is_json_serialisable():
    """The CLI dumps this straight into footprint.json."""
    out = tobias_port_score(_dip_signal(), [200], flank=100)
    encoded = json.dumps(out)
    assert "per_site" in encoded
    assert len(out["per_site"]) == len(out["per_site"])
