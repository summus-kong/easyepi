"""DMR calling: site counting, Wald test, region merging, CLI contract."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from easyepi.dmr import (_split_bams, _wald_log_odds, differential_sites,
                         group_counts, merge_dmr_regions, site_counts)

pysam = pytest.importorskip("pysam", reason="pysam is required to build test BAMs")


def _meth_bam(path, sites, chrom="chr1", length=2000):
    """Write a BAM whose XM tags encode the requested calls.

    ``sites`` is a list of ``(pos, methylated)``; each becomes one read of length 1
    at that position, which keeps the fixture's arithmetic obvious.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": chrom, "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (pos, meth) in enumerate(sorted(sites)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "C"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = int(pos)
            a.mapping_quality = 60
            a.cigartuples = [(0, 1)]
            a.query_qualities = pysam.qualitystring_to_array("I")
            a.set_tag("XM", "Z" if meth else "z")
            out.write(a)
    pysam.index(str(path))
    return path


# --------------------------------------------------------------------------
# Wald test arithmetic
# --------------------------------------------------------------------------
def test_wald_is_zero_when_levels_are_equal():
    diff, stat, p = _wald_log_odds([10], [20], [10], [20])
    assert diff[0] == pytest.approx(0.0)
    assert stat[0] == pytest.approx(0.0)
    assert p[0] == pytest.approx(1.0)


def test_wald_difference_is_group_b_minus_group_a():
    """The sign and magnitude must follow B - A after the continuity correction.

    The correction is (m + 0.5) / (n + 1), so 2/10 vs 8/10 gives 0.590909..., not
    the raw 0.6. Asserting the exact corrected value keeps the test tied to the
    documented estimator rather than to the uncorrected proportion.
    """
    diff, _stat, p = _wald_log_odds([2], [10], [8], [10])
    assert diff[0] == pytest.approx(8.5 / 11 - 2.5 / 11)
    assert diff[0] > 0
    assert p[0] < 0.05


def test_wald_sign_symmetry():
    a = _wald_log_odds([2], [10], [8], [10])
    b = _wald_log_odds([8], [10], [2], [10])
    assert a[0][0] == pytest.approx(-b[0][0])
    assert a[2][0] == pytest.approx(b[2][0])


def test_zero_cells_stay_finite():
    """Without the continuity correction this would be an infinite statistic."""
    _diff, stat, p = _wald_log_odds([0], [10], [10], [10])
    assert np.isfinite(stat[0])
    assert 0.0 <= p[0] <= 1.0


def test_more_evidence_gives_a_smaller_p_value():
    _d1, _s1, p_small_n = _wald_log_odds([1], [5], [4], [5])
    _d2, _s2, p_large_n = _wald_log_odds([100], [500], [400], [500])
    assert p_large_n[0] < p_small_n[0]


# --------------------------------------------------------------------------
# Site counts from a BAM
# --------------------------------------------------------------------------
def test_site_counts_separates_strands(tmp_path):
    bam = _meth_bam(tmp_path / "a.bam", [(100, True), (100, False)])
    frame, stats = site_counts(bam)
    # both reads are forward-strand, so they share one (pos, strand) cell
    assert stats["encoding"] == "XM"
    assert len(frame) == 1
    assert frame.iloc[0]["n_meth"] == 1
    assert frame.iloc[0]["n_total"] == 2


def test_site_counts_refuses_a_bam_without_tags(tmp_path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 1000}]}
    bam = tmp_path / "plain.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = "r"; a.query_sequence = "C"; a.flag = 0
        a.reference_id = 0; a.reference_start = 10
        a.mapping_quality = 60; a.cigartuples = [(0, 1)]
        out.write(a)
    with pytest.raises(ValueError, match="neither an MM/ML tag"):
        site_counts(bam)


def test_site_counts_rejects_unknown_chromosome(tmp_path):
    bam = _meth_bam(tmp_path / "a.bam", [(100, True)])
    with pytest.raises(ValueError, match="not in"):
        site_counts(bam, chroms=["chrZ"])


def test_site_counts_rejects_unknown_context(tmp_path):
    bam = _meth_bam(tmp_path / "a.bam", [(100, True)])
    with pytest.raises(ValueError, match="unknown methylation context"):
        site_counts(bam, contexts=("CpX",))


# --------------------------------------------------------------------------
# Differentially tested sites
# --------------------------------------------------------------------------
def _counts(rows_a, rows_b):
    idx = pd.MultiIndex.from_tuples([(r[0], r[1], r[2]) for r in rows_a],
                                    names=["chrom", "pos", "strand"])
    a = pd.DataFrame({"n_meth": [r[3] for r in rows_a],
                      "n_total": [r[4] for r in rows_a]}, index=idx)
    idx_b = pd.MultiIndex.from_tuples([(r[0], r[1], r[2]) for r in rows_b],
                                      names=["chrom", "pos", "strand"])
    b = pd.DataFrame({"n_meth": [r[3] for r in rows_b],
                      "n_total": [r[4] for r in rows_b]}, index=idx_b)
    return a, b


def test_differential_sites_flags_a_clear_difference():
    a, b = _counts(
        [("chr1", p, "+", 1, 20) for p in range(100, 110)],
        [("chr1", p, "+", 19, 20) for p in range(100, 110)])
    tested, stats = differential_sites(a, b, min_coverage=5)
    assert stats["sites_tested"] == 10
    assert (tested["significant"]).all()
    # group B is more methylated, so the diff must be positive
    assert (tested["diff"] > 0).all()


def test_low_coverage_sites_are_dropped_per_group():
    a, b = _counts([("chr1", 100, "+", 1, 2)], [("chr1", 100, "+", 19, 20)])
    tested, stats = differential_sites(a, b, min_coverage=5)
    assert stats["dropped_low_coverage"] == 1
    assert tested.empty


def test_min_diff_filters_small_effects():
    a, b = _counts([("chr1", p, "+", 9, 20) for p in range(100, 120)],
                   [("chr1", p, "+", 11, 20) for p in range(100, 120)])
    loose, _ = differential_sites(a, b, min_coverage=5, min_diff=0.0)
    strict, _ = differential_sites(a, b, min_coverage=5, min_diff=0.5)
    assert strict["significant"].sum() <= loose["significant"].sum()
    assert strict["significant"].sum() == 0


def test_padj_is_never_below_the_raw_pvalue():
    a, b = _counts([("chr1", p, "+", 2, 30) for p in range(100, 140)],
                   [("chr1", p, "+", 20, 30) for p in range(100, 140)])
    tested, _ = differential_sites(a, b, min_coverage=5)
    finite = tested["padj"].notna()
    assert (tested.loc[finite, "padj"] >= tested.loc[finite, "pvalue"] - 1e-12).all()


def test_sites_only_in_one_group_are_still_tested():
    a, b = _counts([("chr1", 100, "+", 1, 20), ("chr1", 200, "+", 1, 20)],
                   [("chr1", 100, "+", 19, 20)])
    tested, stats = differential_sites(a, b, min_coverage=5)
    # the A-only site has zero B coverage, so it fails the per-group filter
    assert stats["sites_union"] == 2
    assert stats["sites_tested"] == 1


# --------------------------------------------------------------------------
# Region merging
# --------------------------------------------------------------------------
def _tested_rows(positions, sig=True):
    rows = []
    for p in positions:
        rows.append({"chrom": "chr1", "pos": p, "strand": "+",
                     "n_meth": 19, "n_total": 20, "m_a": 1, "n_a": 20,
                     "m_b": 19, "n_b": 20, "level_a": 0.05, "level_b": 0.95,
                     "diff": 0.9, "statistic": 5.0, "pvalue": 1e-6,
                     "padj": 1e-5 if sig else 0.9, "significant": sig})
    return pd.DataFrame(rows)


def test_adjacent_sites_merge_into_one_region():
    tested = _tested_rows([100, 110, 120])
    regions, stats = merge_dmr_regions(tested, max_gap=100)
    assert stats["regions"] == 1
    row = regions.iloc[0]
    assert row["n_sites"] == 3
    assert row["start"] == 100 and row["end"] == 121


def test_a_large_gap_splits_regions():
    tested = _tested_rows([100, 110, 500, 510])
    regions, stats = merge_dmr_regions(tested, max_gap=50)
    assert stats["regions"] == 2


def test_region_test_uses_summed_counts_not_the_min_p():
    tested = _tested_rows([100, 110])
    regions, _ = merge_dmr_regions(tested, max_gap=100)
    row = regions.iloc[0]
    assert row["m_a"] == 2 and row["n_a"] == 40
    assert row["m_b"] == 38 and row["n_b"] == 40


def test_region_direction_labels_hyper_and_hypo():
    """Direction follows the region's recomputed effect, so the counts must agree
    with the sign: overwriting ``diff`` alone would leave the fixture contradictory."""
    hyper = _tested_rows([100])
    regions, _ = merge_dmr_regions(hyper)
    assert regions.iloc[0]["direction"] == "hyper"

    hypo = hyper.copy()
    hypo["m_a"], hypo["m_b"] = hyper["m_b"], hyper["m_a"]
    hypo["n_a"], hypo["n_b"] = hyper["n_b"], hyper["n_a"]
    hypo["diff"] = -0.9
    regions2, _ = merge_dmr_regions(hypo)
    assert regions2.iloc[0]["direction"] == "hypo"


def test_min_sites_filters_singleton_regions():
    tested = _tested_rows([100, 500, 510])
    regions, stats = merge_dmr_regions(tested, max_gap=50, min_sites=2)
    assert stats["regions"] == 1
    assert regions.iloc[0]["n_sites"] == 2


def test_no_significant_sites_yields_no_regions():
    tested = _tested_rows([100, 110], sig=False)
    regions, stats = merge_dmr_regions(tested, alpha=0.05)
    assert regions.empty and stats["regions"] == 0


def test_empty_input_is_handled():
    regions, stats = merge_dmr_regions(pd.DataFrame())
    assert regions.empty
    assert stats == {}


# --------------------------------------------------------------------------
# Group handling and CLI
# --------------------------------------------------------------------------
def test_split_bams_expands_commas(tmp_path):
    a = tmp_path / "a.bam"; a.write_bytes(b"x")
    b = tmp_path / "b.bam"; b.write_bytes(b"x")
    assert _split_bams([f"{a},{b}"], "--group-a") == [str(a), str(b)]


def test_split_bams_rejects_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="BAM not found"):
        _split_bams([str(tmp_path / "nope.bam")], "--group-a")


def test_split_bams_rejects_an_empty_group():
    with pytest.raises(ValueError, match="is empty"):
        _split_bams([], "--group-a")


def test_split_bams_rejects_duplicates(tmp_path):
    a = tmp_path / "a.bam"; a.write_bytes(b"x")
    with pytest.raises(ValueError, match="twice"):
        _split_bams([f"{a},{a}"], "--group-a")


def test_group_counts_sums_across_samples(tmp_path):
    b1 = _meth_bam(tmp_path / "s1.bam", [(100, True)])
    b2 = _meth_bam(tmp_path / "s2.bam", [(100, False)])
    counts, stats = group_counts([b1, b2], label="A")
    key = ("chr1", 100, "+")
    assert counts.loc[key, "n_total"] == 2
    assert counts.loc[key, "n_meth"] == 1
    assert stats["sites"] == 1


def test_cli_is_registered():
    from easyepi import cli
    parser = cli.build_parser()
    a = parser.parse_args(["dmr", "--group-a", "a.bam", "--group-b", "b.bam",
                           "--out", "dmr.tsv"])
    assert a.module == "dmr"
    assert a.min_coverage == 5
    assert a.alpha == 0.05


def test_cli_end_to_end_writes_three_files(tmp_path, capsys):
    from easyepi.dmr import run_dmr

    sites_a = [(p, False) for p in range(100, 120)] * 20
    sites_b = [(p, True) for p in range(100, 120)] * 20
    a = _meth_bam(tmp_path / "a.bam", sites_a)
    b = _meth_bam(tmp_path / "b.bam", sites_b)

    class Args:
        group_a = [str(a)]
        group_b = [str(b)]
        label_a = "ctrl"
        label_b = "treat"
        context = "CpG"
        chrom = None
        min_coverage = 5
        min_mapq = 10
        min_diff = 0.0
        alpha = 0.05
        max_gap = 100
        min_sites = 1
        out = str(tmp_path / "out" / "dmr.tsv")

    assert run_dmr(Args) == 0
    assert (tmp_path / "out" / "dmr.tsv").exists()
    assert (tmp_path / "out" / "dmr.sites.tsv").exists()
    summary = json.loads((tmp_path / "out" / "dmr.json").read_text())
    assert summary["groups"]["a"]["label"] == "ctrl"
    assert "not beta-binomial" in summary["limitation"]

    dmr = pd.read_csv(tmp_path / "out" / "dmr.tsv", sep="\t")
    assert not dmr.empty
    assert (dmr["direction"] == "hyper").all()


def test_cli_rejects_one_bam_in_both_groups(tmp_path):
    from easyepi.dmr import run_dmr
    bam = _meth_bam(tmp_path / "same.bam", [(100, True)])

    class Args:
        group_a = [str(bam)]
        group_b = [str(bam)]
        label_a = None
        label_b = None
        context = "CpG"
        chrom = None
        min_coverage = 5
        min_mapq = 10
        min_diff = 0.0
        alpha = 0.05
        max_gap = 100
        min_sites = 1
        out = str(tmp_path / "dmr.tsv")

    with pytest.raises(ValueError, match="same BAM appears in both groups"):
        run_dmr(Args)
