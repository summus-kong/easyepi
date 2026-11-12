"""Genomic-region classification, checked against hand-computed expectations.

`classify_peaks_to_genomic_bins` and `_make_tx_feature_sets` had no test coverage at
all: the only two references to `easyepi.annotate` internals in the whole suite were
`_extract_utr` and `_nearest_tss_distance`. These tests exist so that the ChIPseeker
benchmark has a second, independent check: if the benchmark disagrees, a failure here
says whether the bug is local and obvious or a genuine modelling difference.

Every interval below is chosen so the expected label follows from the priority order
`Promoter > 5'UTR > Exon > Intron > 3'UTR > Intergenic` and can be read off the
coordinates without running the code.

One gene is used, on the plus strand, with these coordinates (0-based, half-open):

    gene      chr1:  1,000 - 20,000
    exon 1    chr1:  1,000 -  2,000   (contains the 5'UTR)
    CDS 1     chr1:  1,800 -  2,000   -> exon 1's 1,000-1,800 is 5'UTR
    intron    chr1:  2,000 -  9,000
    exon 2    chr1:  9,000 - 10,000
    CDS 2     chr1:  9,000 -  9,800   -> exon 2's 9,800-10,000 is 3'UTR
    exon 3    chr1: 12,000 - 13,000
    intron    chr1: 13,000 - 19,000
    exon 4    chr1: 19,000 - 20,000
    TSS       chr1:  1,000 (plus strand)
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.annotate import (GENOMIC_BINS, _make_tx_feature_sets,
                              classify_peaks_to_genomic_bins)

EXONS = [(1000, 2000), (9000, 10000), (12000, 13000), (19000, 20000)]
CDS = [(1800, 2000), (9000, 9800)]


def _tx(chrom="chr1", start=1000, end=20000, strand="+",
        exons=None, cds=None):
    exons = EXONS if exons is None else exons
    return {"gene_id": "G1", "chrom": chrom, "start": start, "end": end,
            "strand": strand, "tss": start if strand == "+" else end - 1,
            "tes": end if strand == "+" else start,
            "exons": list(exons), "cds": list(cds if cds is not None else CDS)}


def _classify(intervals, tx=None, **kwargs):
    """Classify ``[(start, end), ...]`` on chr1 and return the labels."""
    tx = tx or _tx()
    fsets = _make_tx_feature_sets([tx], **kwargs)
    starts = np.array([s for s, _ in intervals], dtype=np.int64)
    ends = np.array([e for _, e in intervals], dtype=np.int64)
    chroms = np.array(["chr1"] * len(intervals), dtype=object)
    labels, _ = classify_peaks_to_genomic_bins(chroms, starts, ends, fsets,
                                               mode="midpoint")
    return list(labels)


# --------------------------------------------------------------------- basics
def test_six_bins_are_declared_in_priority_order():
    """The order is the contract: first match wins.

    UTRs precede Exon, as in ChIPseeker. With Exon first, every UTR peak would be
    captured by the exon interval containing it and the UTR bins would be empty.
    """
    assert GENOMIC_BINS == ["Promoter", "5'UTR", "3'UTR", "Exon", "Intron",
                            "Intergenic"]


def test_midpoint_inside_promoter_is_promoter():
    # TSS at 1000, window 1000-3000 for the default +/-3000 (clamped at 0)
    assert _classify([(1500, 1600)]) == ["Promoter"]


def test_promoter_window_is_parameterised():
    """A 50 bp window must not reach a peak 500 bp from the TSS."""
    assert _classify([(1500, 1600)], promoter_upstream=50,
                     promoter_downstream=50) == ["5'UTR"]


def test_promoter_wins_over_exon_when_both_cover_the_midpoint():
    # 1,500 is inside exon 1 (1,000-2,000) and inside TSS+/-3000.
    assert _classify([(1450, 1550)]) == ["Promoter"]


def test_five_prime_utr_outside_the_promoter_window():
    """A peak in the 5'UTR but beyond a narrow promoter window is 5'UTR.

    Midpoint 1,700 is inside exon 1 (1,000-2,000) and below the CDS start (1,800),
    so it is untranslated; with a +/-50 window around the TSS at 1,000 it is far
    outside the promoter.
    """
    assert _classify([(1650, 1750)], promoter_upstream=50,
                     promoter_downstream=50) == ["5'UTR"]


def test_exonic_region_between_cds_blocks_is_exon():
    # midpoint 9,500 is inside exon 2 and inside CDS 2: exon, not UTR.
    assert _classify([(9450, 9550)], promoter_upstream=0,
                     promoter_downstream=0) == ["Exon"]


def test_three_prime_utr_covers_every_exon_past_the_stop_codon():
    """All exonic bases beyond the stop codon are 3'UTR, not just the first block.

    A 36-exon gene typically has several exons after translation ends; every one of
    them is untranslated. Keeping only the piece abutting the stop codon left 91
    peaks labelled Exon where ChIPseeker says 3'UTR (measured on 43,465 peaks).

    Exon 2 (9,000-10,000) contains CDS 2 (9,000-9,800), so 9,800-10,000 is that
    exon's untranslated tail, and exons 3 and 4 lie wholly past the stop codon.
    """
    fsets = _make_tx_feature_sets([_tx()], promoter_upstream=0,
                                  promoter_downstream=0)
    assert list(zip(fsets["3'UTR"][1].tolist(),
                    fsets["3'UTR"][2].tolist())) == [
        (9800, 10000), (12000, 13000), (19000, 20000)]
    # a peak in a downstream untranslated exon is 3'UTR, not Exon
    assert _classify([(12400, 12600)], promoter_upstream=0,
                     promoter_downstream=0) == ["3'UTR"]


def test_utr_is_exonic_and_non_coding_not_everything_past_the_cds():
    """3'UTR must exclude internal exons that simply carry no CDS.

    Exons 3 (12,000-13,000) and 4 (19,000-20,000) lie between the first and last CDS
    block, so the old global min/max test called them 3'UTR. They are ordinary
    exons: their sequence is transcribed but the gene's coding region is already
    closed. Counting them as 3'UTR inflates that bin on real gene models.
    """
    fsets = _make_tx_feature_sets([_tx()], promoter_upstream=0,
                                  promoter_downstream=0)
    three = list(zip(fsets["3'UTR"][1].tolist(), fsets["3'UTR"][2].tolist()))
    # Exons 3 and 4 carry no CDS and lie past the stop codon, so they are
    # untranslated exonic sequence: 3'UTR. What must NOT happen is the reverse --
    # an exon that *contains* CDS being swallowed into 3'UTR because it sits between
    # the first and last CDS block.
    assert three == [(9800, 10000), (12000, 13000), (19000, 20000)], three
    assert _classify([(12400, 12600)], promoter_upstream=0,
                     promoter_downstream=0) == ["3'UTR"]
    # and a coding exon is still Exon, not UTR
    assert _classify([(9450, 9550)], promoter_upstream=0,
                     promoter_downstream=0) == ["Exon"]


def test_intron_between_exons():
    # midpoint 5,000 is between exon 1 (ends 2,000) and exon 2 (starts 9,000).
    assert _classify([(4950, 5050)], promoter_upstream=0,
                     promoter_downstream=0) == ["Intron"]


def test_far_away_peak_is_intergenic():
    assert _classify([(50000, 50100)], promoter_upstream=0,
                     promoter_downstream=0) == ["Intergenic"]


# ------------------------------------------------------- priority, not pairwise
def test_exon_wins_over_intron_at_a_shared_boundary():
    """Priority decides when two bins both contain the midpoint."""
    # 1,900 is inside exon 1 and inside CDS 1; no intron overlaps it.
    assert _classify([(1850, 1950)], promoter_upstream=0,
                     promoter_downstream=0) == ["Exon"]


def test_five_prime_utr_beats_exon_by_declared_priority():
    """5'UTR is listed before Exon, so a 5'UTR midpoint reports 5'UTR.

    This is a deliberate difference from ChIPseeker, whose default priority is
    Promoter > 5UTR > 3UTR > Exon > Intron > Downstream > Intergenic. The benchmark
    measures the effect; pinning it here stops the value drifting silently.
    """
    # midpoint 1,200: inside exon 1 (1,000-2,000), outside CDS, left of TSS+3000
    # but with the promoter window disabled and the TSS moved far right.
    tx = _tx(start=5000)
    labels = _classify([(1150, 1250)], tx=tx, promoter_upstream=0,
                       promoter_downstream=0)
    # 1,150 is inside exon 1 (1,000-2,000) and inside the 5'UTR piece (1,000-1,800),
    # so the declared priority must return 5'UTR rather than Exon.
    assert labels == ["5'UTR"]
    fsets = _make_tx_feature_sets([tx], promoter_upstream=0, promoter_downstream=0)
    assert fsets["5'UTR"][1].tolist() == [1000]
    assert fsets["5'UTR"][2].tolist() == [1800]


# --------------------------------------------------------------- strand matters
def test_minus_strand_swaps_utr_sides():
    """On the minus strand the leftmost exon block is 3'UTR, not 5'UTR.

    The transcript is reversed so that the gene spans the same coordinates but the
    TSS is at the right-hand end.
    """
    tx = _tx(strand="-")
    # With strand "-", tss = end-1 = 19,999 and the transcript runs right-to-left,
    # so the block at 1,000-1,800 is the 3' end of the mRNA.
    fsets = _make_tx_feature_sets([tx], promoter_upstream=0,
                                  promoter_downstream=0)
    five = list(zip(fsets["5'UTR"][1].tolist(), fsets["5'UTR"][2].tolist()))
    three = list(zip(fsets["3'UTR"][1].tolist(), fsets["3'UTR"][2].tolist()))

    # On the minus strand the transcript runs right-to-left, so translation starts
    # at the right-hand CDS edge (9,800) and ends at the left-hand one (1,800).
    # Everything to the right of 9,800 is 5'UTR, everything left of 1,800 is 3'UTR --
    # including whole exons, on both sides.
    # Getting this backwards is the classic strand bug, which is why it is asserted
    # on the exact intervals rather than only on the labels.
    assert five == [(9800, 10000), (12000, 13000), (19000, 20000)], five
    assert three == [(1000, 1800)], three
    assert _classify([(9850, 9950)], tx=tx, promoter_upstream=0,
                     promoter_downstream=0) == ["5'UTR"]
    assert _classify([(1150, 1250)], tx=tx, promoter_upstream=0,
                     promoter_downstream=0) == ["3'UTR"]


def test_tss_distance_is_signed_by_strand():
    """A peak left of a plus-strand TSS is negative; left of a minus-strand one is
    positive. ChIPseeker reports the same convention, so the benchmark compares it."""
    from easyepi.annotate import _nearest_tss_distance

    # plus strand TSS at 1,000; query midpoint at 800 => -200
    # Argument order is (query_starts, query_ends, query_chrom, tss_pos, tss_chrom,
    # tss_ids, tss_strand) -- tss_pos comes before tss_chrom.
    ids, dists = _nearest_tss_distance(np.array([700]), np.array([900]),
                                       np.array(["chr1"], dtype=object),
                                       np.array([1000]),
                                       np.array(["chr1"], dtype=object),
                                       np.array(["G1"], dtype=object),
                                       np.array(["+"]))
    assert ids[0] == "G1"
    assert dists[0] == pytest.approx(-200)

    # minus strand: the same genomic offset is *downstream*, so the sign flips
    ids, dists = _nearest_tss_distance(np.array([19700]), np.array([19900]),
                                       np.array(["chr1"], dtype=object),
                                       np.array([19999]),
                                       np.array(["chr1"], dtype=object),
                                       np.array(["G1"], dtype=object),
                                       np.array(["-"]))
    assert dists[0] == pytest.approx(199)


# ------------------------------------------------------- robustness / edge cases
def test_peak_on_another_chromosome_is_not_assigned_to_this_gene():
    """A chr2 peak must not pick up chr1 features: the sweep is per chromosome."""
    tx = _tx()
    fsets = _make_tx_feature_sets([tx], promoter_upstream=3000,
                                  promoter_downstream=3000)
    labels, _ = classify_peaks_to_genomic_bins(
        np.array(["chr2"], dtype=object), np.array([1500]), np.array([1600]),
        fsets, mode="midpoint")
    assert list(labels) == ["Intergenic"]


def test_promoter_start_is_clamped_at_zero():
    """A TSS near the contig start must not emit a negative coordinate."""
    tx = _tx(start=100)
    fsets = _make_tx_feature_sets([tx], promoter_upstream=3000,
                                  promoter_downstream=3000)
    starts = fsets["Promoter"][1]
    assert (starts >= 0).all(), starts


def test_a_negative_window_is_rejected_by_the_cli():
    """An inverted window would silently produce an empty Promoter bin."""
    from easyepi.annotate import run_annotate

    class Args:
        peaks = "does-not-matter.bed"
        out = "out"
        gtf = "does-not-matter.gtf"
        backend = "native"
        mode = "midpoint"
        min_overlap = 0.0
        promoter_upstream = -1
        promoter_downstream = 3000

    with pytest.raises((ValueError, FileNotFoundError)):
        run_annotate(Args())


def test_gene_without_cds_produces_no_utr_bins():
    """UTR is defined relative to the CDS; with no CDS there are none."""
    tx = _tx(cds=[])
    fsets = _make_tx_feature_sets([tx])
    assert "5'UTR" not in fsets and "3'UTR" not in fsets


def test_exon_only_gene_still_classifies_exons():
    tx = _tx(cds=[])
    assert _classify([(9450, 9550)], tx=tx, promoter_upstream=0,
                     promoter_downstream=0) == ["Exon"]


# ------------------------------------------------ per-transcript TSS model (0.4.0)
def test_gtf_parser_keys_on_transcript_not_gene(tmp_path):
    """One record per transcript, carrying its gene id.

    Keying on ``gene_id`` gave every gene a single TSS, which is what made the
    Promoter bin 29% instead of 65% against ChIPseeker.
    """
    from easyepi.annotate import _read_gtf_txs

    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        'chr1\tsrc\tgene\t1001\t9000\t.\t+\t.\tgene_id "G1"; gene_name "A";\n'
        'chr1\tsrc\ttranscript\t1001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\tCDS\t1501\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\ttranscript\t7001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
        'chr1\tsrc\texon\t7001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n',
        encoding="utf-8")

    txs = _read_gtf_txs(str(gtf))
    assert len(txs) == 2, [t["transcript_id"] for t in txs]
    assert {t["transcript_id"] for t in txs} == {"T1", "T2"}
    assert {t["gene_id"] for t in txs} == {"G1"}
    by_tx = {t["transcript_id"]: t for t in txs}
    # each transcript's own start, not the gene's
    assert by_tx["T1"]["tss"] == 1000
    assert by_tx["T2"]["tss"] == 7000
    # the gene span is still available for output
    assert by_tx["T1"]["gene_start"] == 1000
    assert by_tx["T1"]["gene_end"] == 9000


def test_multiple_transcripts_widen_the_promoter_bin(tmp_path):
    """Two isoforms starting 6 kb apart must produce two promoter regions."""
    from easyepi.annotate import _read_gtf_txs

    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        'chr1\tsrc\tgene\t1001\t9000\t.\t+\t.\tgene_id "G1";\n'
        'chr1\tsrc\ttranscript\t1001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t1001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\ttranscript\t7001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
        'chr1\tsrc\texon\t7001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n',
        encoding="utf-8")

    fsets = _make_tx_feature_sets(_read_gtf_txs(str(gtf)))
    starts = sorted(fsets["Promoter"][1].tolist())
    assert len(starts) == 2, starts
    # The whole point of the change: a peak at T2's start is a promoter. Under the
    # old one-TSS-per-gene model only the region around 1,000 existed and this peak
    # fell outside it.
    starts_arr = np.array([7050], dtype=np.int64)
    ends_arr = np.array([7150], dtype=np.int64)
    chroms = np.array(["chr1"], dtype=object)
    labels, _ = classify_peaks_to_genomic_bins(chroms, starts_arr, ends_arr, fsets,
                                               mode="midpoint")
    assert labels[0] == "Promoter", labels


def test_minus_strand_transcript_tss_uses_the_right_edge(tmp_path):
    from easyepi.annotate import _read_gtf_txs

    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        'chr1\tsrc\tgene\t1001\t5000\t.\t-\t.\tgene_id "G1";\n'
        'chr1\tsrc\ttranscript\t1001\t5000\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t1001\t2000\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t4001\t5000\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8")

    (tx,) = _read_gtf_txs(str(gtf))
    # Minus strand: the transcript starts at the right-hand end. exons are stored
    # 0-based half-open, so the rightmost exon ends at 5000 and that is the TSS.
    assert tx["tss"] == 5000
    assert tx["tes"] == 1000


def test_gtf_without_transcript_lines_still_parses(tmp_path):
    """GFF3-style exon-only annotation degrades to one unit per gene."""
    from easyepi.annotate import _read_gtf_txs

    gtf = tmp_path / "mini.gff3"
    gtf.write_text(
        'chr1\tsrc\tgene\t1001\t5000\t.\t+\t.\tgene_id "G1";\n'
        'chr1\tsrc\texon\t1001\t2000\t.\t+\t.\tgene_id "G1";\n'
        'chr1\tsrc\texon\t4001\t5000\t.\t+\t.\tgene_id "G1";\n',
        encoding="utf-8")

    txs = _read_gtf_txs(str(gtf))
    assert len(txs) == 1
    assert txs[0]["gene_id"] == "G1"
    assert len(txs[0]["exons"]) == 2


# ------------------------------------------- sweep correctness (the 0.5.1 bug)
def test_overlap_scan_finds_a_feature_that_starts_before_but_ends_after():
    """A feature nested inside a later-starting one must still be found.

    The classifier sorted features by start and then binary-searched the *end* list,
    which is not sorted. The returned position was meaningless, so the scan started
    somewhere arbitrary and skipped overlapping features -- 2,598 peaks were lost on
    a 43,465-peak benchmark (93.78% -> 99.76% once fixed). This case is small enough
    that the wrong start index provably misses the answer.
    """
    # Intervals deliberately sorted by start but with ends out of order: the first
    # one is long, the later ones are short.
    feature_sets = {
        "Intron": (np.array(["chr1", "chr1", "chr1"], dtype=object),
                   np.array([100, 300, 500]),
                   np.array([9000, 400, 600]),   # ends: 9000, 400, 600 (unsorted)
                   np.array(["+", "+", "+"], dtype=object)),
    }
    # A peak at 4,000-4,100 overlaps only the long first interval.
    labels, _ = classify_peaks_to_genomic_bins(
        np.array(["chr1"], dtype=object), np.array([4000]), np.array([4100]),
        feature_sets, mode="overlap")
    assert labels[0] == "Intron", labels


def test_overlap_scan_is_independent_of_feature_order():
    """Shuffling the feature list must not change the answer."""
    rng = np.random.default_rng(7)
    n = 60
    starts = rng.integers(0, 5000, n)
    lengths = rng.integers(50, 900, n)
    ends = starts + lengths
    chroms = np.array(["chr1"] * n, dtype=object)
    strands = np.array(["+"] * n, dtype=object)

    peaks_c = np.array(["chr1"] * 200, dtype=object)
    peaks_s = rng.integers(0, 5500, 200)
    peaks_e = peaks_s + rng.integers(50, 300, 200)

    base, _ = classify_peaks_to_genomic_bins(
        peaks_c, peaks_s, peaks_e,
        {"Intron": (chroms, starts, ends, strands)}, mode="overlap")

    order = rng.permutation(n)
    shuffled, _ = classify_peaks_to_genomic_bins(
        peaks_c, peaks_s, peaks_e,
        {"Intron": (chroms[order], starts[order], ends[order], strands[order])},
        mode="overlap")
    assert (base == shuffled).all(), int((base != shuffled).sum())


def test_no_feature_overlapping_a_peak_leaves_it_unassigned():
    feature_sets = {
        "Intron": (np.array(["chr1"], dtype=object), np.array([100]),
                   np.array([200]), np.array(["+"], dtype=object)),
    }
    labels, _ = classify_peaks_to_genomic_bins(
        np.array(["chr1"], dtype=object), np.array([900]), np.array([950]),
        feature_sets, mode="overlap")
    assert labels[0] == "Intergenic", labels
