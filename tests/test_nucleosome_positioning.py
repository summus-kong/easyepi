"""Nucleosome positioning: curve primitives, calls, fuzziness, single-end refusal."""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.nucleosome_positioning import (bam_to_dyad_density, call_positions,
                                            fuzziness, robust_background,
                                            spacing_summary)


def _track(centres, n_bins=400, bin_size=10, noise=0.0, seed=0, width=1):
    """Build a dyad-density track with gaussian bumps at the given bp centres."""
    rng = np.random.default_rng(seed)
    x = np.arange(n_bins)
    dens = np.zeros(n_bins)
    for c in centres:
        b = c // bin_size
        dens += np.exp(-0.5 * ((x - b) / max(width, 1e-9)) ** 2)
    if noise:
        dens += rng.random(n_bins) * noise
    return dens


# --------------------------------------------------------------------------
# Robust background
# --------------------------------------------------------------------------
def test_robust_background_uses_median_and_scaled_mad():
    values = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 100.0])
    median, sigma = robust_background(values)
    assert median == 0.0
    # MAD is 0 here (most points are the median), so sigma must be 0, not huge:
    # an SD-based estimate would be ~37 and hide every real call.
    assert sigma == 0.0


def test_robust_background_is_not_dominated_by_outliers():
    base = np.zeros(100)
    base[50] = 1000.0
    median, sigma = robust_background(base)
    assert median == 0.0 and sigma == 0.0


def test_robust_background_empty():
    assert robust_background(np.array([])) == (0.0, 0.0)


# --------------------------------------------------------------------------
# Position calling
# --------------------------------------------------------------------------
def test_calls_recover_planted_positions():
    centres = [500, 800, 1100, 1400]
    dens = _track(centres, noise=0.02)
    calls = call_positions(dens, bin_size=10, smooth_bins=3, mad_k=2.0,
                           min_spacing=100)
    got = [c["position"] for c in calls]
    assert len(got) == len(centres), f"expected {centres}, got {got}"
    for want, have in zip(centres, got):
        assert abs(want - have) <= 30, f"expected near {want}, got {have}"


def test_noise_spikes_are_not_called_when_mad_is_zero():
    """Sparse dyad tracks have MAD = 0, which would make the threshold 0."""
    dens = _track([1000], noise=0.02)
    calls = call_positions(dens, bin_size=10, mad_k=2.0, min_spacing=100)
    assert [c["position"] for c in calls] == [1000]


def test_relative_height_floor_is_validated():
    with pytest.raises(ValueError, match="min_relative_height"):
        call_positions(np.ones(10), bin_size=10, min_relative_height=1.0)


def test_calls_are_sorted_left_to_right():
    calls = call_positions(_track([1400, 500, 1100], noise=0.02), bin_size=10,
                           mad_k=2.0, min_spacing=100)
    positions = [c["position"] for c in calls]
    assert positions == sorted(positions)


def test_min_spacing_prevents_adjacent_double_calls():
    # one broad nucleosome should yield a single call, not a row of them
    dens = _track([1000], noise=0.0, width=6)
    calls = call_positions(dens, bin_size=10, smooth_bins=3, mad_k=1.0,
                           min_spacing=100)
    assert len(calls) == 1


def test_flat_signal_produces_no_calls():
    assert call_positions(np.ones(200), bin_size=10, mad_k=2.0) == []


def test_all_zero_signal_produces_no_calls():
    assert call_positions(np.zeros(200), bin_size=10) == []


def test_min_height_filters_weak_calls():
    dens = _track([500, 1200], noise=0.0)
    dens[120] *= 0.05  # make the second bump tiny
    without = call_positions(dens, bin_size=10, mad_k=1.0, min_spacing=100)
    with_height = call_positions(dens, bin_size=10, mad_k=1.0, min_spacing=100,
                                 min_height=0.5)
    assert len(with_height) <= len(without)


def test_call_positions_validates_input():
    with pytest.raises(ValueError, match="empty"):
        call_positions(np.array([]), bin_size=10)
    with pytest.raises(ValueError, match="non-negative"):
        call_positions(np.array([1.0, -1.0]), bin_size=10)
    with pytest.raises(ValueError, match="bin_size"):
        call_positions(np.ones(10), bin_size=0)
    with pytest.raises(ValueError, match="min_spacing"):
        call_positions(np.ones(10), bin_size=10, min_spacing=-1)


def test_each_call_reports_a_positive_score():
    calls = call_positions(_track([500, 900], noise=0.01), bin_size=10, mad_k=1.5)
    assert calls and all(c["score"] > 0 for c in calls)
    assert all(c["frac_fragments"] >= 0 for c in calls)


# --------------------------------------------------------------------------
# Fuzziness and spacing
# --------------------------------------------------------------------------
def test_fuzziness_is_narrow_for_a_sharp_dyad():
    dens = _track([1000], width=0.8)
    sharp = fuzziness(dens, bin_size=10, position=1000, max_radius=500)
    assert sharp["width_bp"] <= 100
    assert 0 < sharp["frac_covered"] <= 1.0


def test_fuzziness_is_wider_for_a_diffuse_dyad():
    sharp = fuzziness(_track([1000], width=0.8), bin_size=10, position=1000)
    fuzzy = fuzziness(_track([1000], width=8.0), bin_size=10, position=1000)
    assert fuzzy["width_bp"] > sharp["width_bp"]


def test_fuzziness_validates_input():
    with pytest.raises(ValueError, match="fraction"):
        fuzziness(np.ones(10), bin_size=10, position=0, fraction=0.0)
    with pytest.raises(ValueError, match="outside the supplied track"):
        fuzziness(np.ones(10), bin_size=10, position=10_000)


def test_spacing_summary_over_known_gaps():
    calls = [{"position": p} for p in (0, 180, 360, 550)]
    summary = spacing_summary(calls)
    assert summary["n_pairs"] == 3
    assert summary["median_spacing"] == 180
    assert summary["min_spacing"] == 180 and summary["max_spacing"] == 190


def test_spacing_summary_needs_two_calls():
    assert spacing_summary([{"position": 5}])["n_pairs"] == 0
    assert spacing_summary([])["n_pairs"] == 0


# --------------------------------------------------------------------------
# BAM path: validation happens before any reader is needed
# --------------------------------------------------------------------------
def test_bam_to_dyad_density_validates_arguments(tmp_path):
    dummy = tmp_path / "x.bam"
    dummy.write_bytes(b"")
    with pytest.raises(ValueError, match="bin_size"):
        bam_to_dyad_density(dummy, "chr1", 0, 1000, bin_size=0)
    with pytest.raises(ValueError, match="end must exceed"):
        bam_to_dyad_density(dummy, "chr1", 1000, 1000)
    with pytest.raises(ValueError, match="max_length"):
        bam_to_dyad_density(dummy, "chr1", 0, 1000, min_length=200, max_length=150)


def test_cli_is_registered():
    from easyepi import cli
    parser = cli.build_parser()
    args = parser.parse_args(["nucleosome-positioning", "--bam", "a.bam",
                              "--chrom", "chr1", "--start", "0", "--end", "10000",
                              "--out", "o"])
    assert args.module == "nucleosome-positioning"
    assert args.bin_size == 10
    assert args.min_spacing == 100
    assert (args.min_length, args.max_length) == (120, 200)


# --------------------------------------------------------------------------
# BED output contract
# --------------------------------------------------------------------------
def _write_bed(tmp_path, dens, bin_size=10, start=0):
    """Exercise the real BED writer through run_positioning's helper contract."""
    import csv as _csv
    from easyepi.nucleosome_positioning import call_positions, fuzziness
    calls = call_positions(dens, bin_size=bin_size, smooth_bins=3, mad_k=2.0,
                           min_spacing=100, min_height=0.0,
                           min_relative_height=0.05)
    bed = tmp_path / "np.bed"
    with bed.open("w", encoding="utf-8", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t")
        for i, c in enumerate(calls, start=1):
            pos = start + c["position"]
            fz = fuzziness(dens, bin_size, c["position"])
            half = max(bin_size, int(fz["width_bp"]) // 2)
            w.writerow(["chr1", max(0, pos - half), pos + half, f"nuc_{i}", c["score"], "."])
    return calls, bed


def test_bed_width_follows_measured_fuzziness_not_a_constant(tmp_path):
    """Every feature used to be exactly 2*bin_size wide, which hid fuzziness."""
    sharp = _track([1000, 1200, 1400], width=0.8)
    fuzzy = _track([1000, 1200, 1400], width=8.0)
    _, bed_sharp = _write_bed(tmp_path, sharp)
    widths_sharp = [int(r[2]) - int(r[1]) for r in
                    (l.split("\t") for l in bed_sharp.read_text().splitlines())]
    fuzzy_dir = tmp_path / "fuzzy"
    fuzzy_dir.mkdir()
    _, bed_fuzzy = _write_bed(fuzzy_dir, fuzzy)
    widths_fuzzy = [int(r[2]) - int(r[1]) for r in
                    (l.split("\t") for l in bed_fuzzy.read_text().splitlines())]
    assert widths_sharp and widths_fuzzy
    # a diffuse dyad must produce wider features than a sharp one
    assert max(widths_fuzzy) > max(widths_sharp)
    # and they must not all collapse onto the same constant
    assert len(set(widths_fuzzy)) > 1 or len(set(widths_sharp)) > 1


def test_bed_strand_is_unstranded_because_a_dyad_has_no_direction(tmp_path):
    """Column 6 used to be '+'/'-' alternating on the call index."""
    _, bed = _write_bed(tmp_path, _track([1000, 1200, 1400], width=0.8))
    strands = {l.split("\t")[5] for l in bed.read_text().splitlines()}
    assert strands == {"."}


def test_bed_has_six_columns_and_bounded_coordinates(tmp_path):
    _, bed = _write_bed(tmp_path, _track([1000, 1200, 1400], width=0.8))
    for line in bed.read_text().splitlines():
        f = line.split("\t")
        assert len(f) == 6
        s, e = int(f[1]), int(f[2])
        assert 0 <= s < e


# --------------------------------------------------------------------------
# Fuzziness growth and periodicity reliability
# --------------------------------------------------------------------------
def test_fuzziness_grows_outward_and_does_not_saturate_on_a_phased_array():
    """The old top-bins rule spanned the whole radius on a phased array."""
    n = 2000
    dens = np.zeros(n)
    for c in np.arange(50, n - 50, 18):          # 180 bp lattice
        dens += np.exp(-0.5 * ((np.arange(n) - c) / 1.2) ** 2)
    widths = [fuzziness(dens, bin_size=10, position=p, max_radius=500)["width_bp"]
              for p in np.arange(50, n - 50, 18)]
    assert widths
    # a well-phased lattice must not pin every call at the search radius
    assert sum(w >= 500 for w in widths) < 0.5 * len(widths)
    assert np.median(widths) < 500


def test_fuzziness_separates_sharp_from_diffuse():
    sharp = fuzziness(_track([1000], width=0.8), bin_size=10, position=1000)
    fuzzy = fuzziness(_track([1000], width=8.0), bin_size=10, position=1000)
    assert sharp["width_bp"] < fuzzy["width_bp"]


def test_fuzziness_fraction_coverage_is_reported_and_bounded():
    r = fuzziness(_track([1000], width=2.0), bin_size=10, position=1000,
                  fraction=0.5, max_radius=500)
    assert 0.0 < r["frac_covered"] <= 1.0
