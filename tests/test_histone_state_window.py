"""Histone-state windows must respect a non-zero --start.

Both halves of the command assumed the window began at genomic coordinate 0:

  * `_read_signal_per_bin` binned `start // resolution` from coordinate 0, so a run
    with `--start 1000000` folded the genome from 0 into a window declared to begin at
    1,000,000. Signal inside the requested window was looked up at the wrong offset and
    mostly fell outside `n_bins`, so those bins came back zero.
  * `states_to_bed` wrote `start * resolution` from 0, so the emitted BED intervals
    were offset by the window origin relative to the genome.

Verified against a synthetic bedGraph: with `origin=1_000_000` the covered bins read
`[2, 2, 2, 2]`; with the old `origin=0` behaviour they read `[0, 0, 0, 0]`.
"""

from __future__ import annotations

import numpy as np
import pytest

from easyepi.histone_states import states_to_bed
from easyepi.histone_states_cli import _read_signal_per_bin

RES = 10
ORIGIN = 1_000_000
N_BINS = 4


@pytest.fixture()
def bedgraph(tmp_path):
    p = tmp_path / "s.bedGraph"
    p.write_text(
        f"chr1\t{ORIGIN}\t{ORIGIN + RES * N_BINS}\t2.0\n"
        f"chr1\t{ORIGIN - 100}\t{ORIGIN}\t7.0\n"      # entirely before the window
        f"chr1\t{ORIGIN + 10_000}\t{ORIGIN + 10_040}\t9.0\n",  # after the window
        encoding="utf-8",
    )
    return p


def test_signal_inside_a_shifted_window_is_found(bedgraph):
    """The whole point: signal at --start must land in bin 0 of the window."""
    got = _read_signal_per_bin(str(bedgraph), "chr1", N_BINS, RES, origin=ORIGIN)
    assert list(got) == [2.0, 2.0, 2.0, 2.0]


def test_origin_zero_would_miss_it(bedgraph):
    """Guards the regression: without the origin the signal is silently dropped."""
    got = _read_signal_per_bin(str(bedgraph), "chr1", N_BINS, RES, origin=0)
    assert list(got) == [0.0, 0.0, 0.0, 0.0]


def test_signal_outside_the_window_is_excluded(bedgraph):
    """Flanks must not bleed in: totals equal only the in-window contribution."""
    got = _read_signal_per_bin(str(bedgraph), "chr1", N_BINS, RES, origin=ORIGIN)
    assert pytest.approx(float(got.sum())) == 2.0 * N_BINS


def test_states_to_bed_offsets_by_origin():
    """A state path over a shifted window must be emitted at genomic coordinates."""
    paths = np.array([0, 0, 1], dtype=np.int64)
    rows = states_to_bed(paths, "chr1", RES, ["E1", "E2"], origin=ORIGIN)
    assert rows == [
        ("chr1", ORIGIN + 0, ORIGIN + 2 * RES, "E1"),
        ("chr1", ORIGIN + 2 * RES, ORIGIN + 3 * RES, "E2"),
    ]


def test_states_to_bed_defaults_to_unshifted():
    """Callers that do not pass an origin keep the previous coordinate basis."""
    paths = np.array([0, 1], dtype=np.int64)
    rows = states_to_bed(paths, "chr1", RES, ["E1", "E2"])
    assert rows[0][1] == 0
    assert rows[0][2] == RES
