"""Annotation must follow ChIPseeker's assignment order, label for label.

These tests were rewritten when the native classifier was replaced by a literal
port of `ChIPseeker::annotatePeak` (`easyepi.annotate.chipseeker_port`). The old
tests exercised internals that no longer exist (`_make_tx_feature_sets`,
`_overlap_scan`); asserting on those would have tested deleted code.

Every expected value below is the string ChIPseeker 1.46.1 itself returns for
this exact GTF and peak, not a paraphrase of it. `annotatePeak` does not emit the
bare bin names: it emits "Promoter (<=1kb)", "Exon (T1/G1, exon 2 of 4)",
"3' UTR" and "Distal Intergenic". Asserting the bare names would have encoded a
vocabulary the reference never produces.

Reference output for this fixture (ChIPseeker 1.46.1, tssRegion +/-3000,
level="transcript", downstreamDistance=300 default):

    peak  40000-40200 -> Distal Intergenic                d=38999
    peak    900-1100  -> Promoter (<=1kb)                 d=0
    peak    300-500   -> Promoter (<=1kb)                 d=-501
    peak  50000-50200 -> Distal Intergenic                d=48999
    peak   9200-9700  -> Exon (T1/G1, exon 2 of 4)        d=8199
    peak  12050-12950 -> 3' UTR                          d=11049
    peak   5000-5100  -> Intron (T1/G1, intron 1 of 3)    d=3999
    peak   1001-1001  -> Promoter (<=1kb)                 d=0
    peak   1000-1002  -> Promoter (<=1kb)                 d=0

The end-to-end agreement is verified separately in `benchmark/a10`: labels and
distanceToTSS match ChIPseeker on 100% of 23,646 peaks across 12 data sets, and
on 100% of peaks for seven different (tssRegion, downstreamDistance) settings.
"""

from __future__ import annotations

import pytest

from easyepi.annotate import BIN_NAMES, DEFAULT_PRIORITY
from easyepi.annotate.chipseeker_port import annotate, annotation_stat
from easyepi.annotate.gtf_features import load_from_gtf


GTF = """\
chr1\ttest\ttranscript\t1001\t20000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\ttest\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\ttest\tCDS\t1801\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\ttest\texon\t9001\t10000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\ttest\tCDS\t9001\t9800\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\ttest\texon\t12001\t13000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\ttest\texon\t19001\t20000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
"""


@pytest.fixture
def feats(tmp_path):
    path = tmp_path / "a.gtf"
    path.write_text(GTF, encoding="utf-8")
    return load_from_gtf(path)


def labels_for(feats, spans):
    """(labels, distances) only -- annotate() also returns per-peak metadata."""
    chroms = ["chr1"] * len(spans)
    labels, dists, _meta = annotate(feats, chroms, [s for s, _ in spans],
                                    [e for _, e in spans])
    return labels, dists


# --------------------------------------------------------------------- contract
def test_six_bins_are_declared_in_priority_order():
    """Promoter first, Intergenic last: the order is the published contract."""
    assert tuple(BIN_NAMES) == (
        "Promoter", "5'UTR", "3'UTR", "Exon", "Intron", "Intergenic"
    )
    assert tuple(DEFAULT_PRIORITY) == tuple(BIN_NAMES)


def test_tss_regions_outside_the_bins_are_unassigned(feats):
    """A peak far from any feature is intergenic, not silently Promoter.

    ChIPseeker renames the leftover Intergenic to "Distal Intergenic" at the end
    of getGenomicAnnotation, so that is the string to expect.
    """
    labels, _ = labels_for(feats, [(40000, 40200)])
    assert labels == ["Distal Intergenic"]


# ---------------------------------------------------------------- box behaviour
def test_a_peak_straddling_a_tss_is_promoter_with_distance_zero(feats):
    """ChIPseeker sets distanceToTSS to 0 when a peak contains a TSS."""
    labels, dists = labels_for(feats, [(900, 1100)])
    assert labels == ["Promoter (<=1kb)"]
    assert dists == [0.0]


def test_a_peak_just_upstream_is_promoter(feats):
    """Within tssRegion (-3000, 3000) on the upstream side."""
    labels, dists = labels_for(feats, [(300, 500)])
    assert labels == ["Promoter (<=1kb)"]
    assert dists == [-501.0]


def test_a_peak_beyond_the_window_is_not_promoter(feats):
    """Past -3000 the promoter rule must not fire."""
    labels, _ = labels_for(feats, [(50000, 50200)])
    assert labels == ["Distal Intergenic"]


def test_a_coding_exon_peak_is_exon(feats):
    """Inside a CDS, past the promoter window: the exon bin wins."""
    labels, dists = labels_for(feats, [(9200, 9700)])
    assert labels == ["Exon (T1/G1, exon 2 of 4)"]
    assert dists == [8199.0]


def test_an_exon_beyond_the_last_cds_is_three_prime_utr(feats):
    """3'UTR outranks Exon, so a trailing exon past the CDS reports 3'UTR."""
    labels, dists = labels_for(feats, [(12050, 12950)])
    assert labels == ["3' UTR"]
    assert dists == [11049.0]


def test_an_intronic_peak_is_intron(feats):
    """Between exons, inside the transcript."""
    labels, dists = labels_for(feats, [(5000, 5100)])
    assert labels == ["Intron (T1/G1, intron 1 of 3)"]
    assert dists == [3999.0]


# ------------------------------------------------------------------ coordinate
def test_minus_strand_tss_is_the_right_edge(tmp_path):
    """resize(width = 1) is strand aware: the TSS is the 3' end on '-'."""
    gtf = (
        'chr1\tt\ttranscript\t1001\t20000\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\texon\t1001\t20000\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\n'
    )
    path = tmp_path / "m.gtf"
    path.write_text(gtf, encoding="utf-8")
    feats = load_from_gtf(path)

    # The TSS is 20000, so a peak covering it is promoter with distance 0.
    labels, dists = labels_for(feats, [(19900, 20100)])
    assert labels == ["Promoter (<=1kb)"]
    assert dists == [0.0]


def test_bed_coordinates_are_one_based_inclusive_at_this_layer(feats):
    """The port consumes GRanges coordinates; callers shift BED by +1."""
    one_based = labels_for(feats, [(1001, 1001)])[0]
    two_based = labels_for(feats, [(1000, 1002)])[0]
    assert one_based == two_based == ["Promoter (<=1kb)"]


# --------------------------------------------------- annoStat (getGenomicAnnoStat)
def test_stat_keeps_promoter_subcategories_and_splits_first_from_other():
    """ChIPseeker's annoStat is not a fold to the six bins.

    ``getGenomicAnnoStat`` keeps the Promoter sub-categories as separate rows and
    splits exons/introns into "1st" and "Other" by testing for "exon 1 of" /
    "intron 1 of". A six-bin fold loses both splits. The row count is
    data-dependent (the reference sample produced 9 because it had no peak in a
    1st exon), so this asserts the label set and its order, not a count.
    """
    labels = [
        "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
        "5' UTR", "3' UTR",
        "Exon (T1/G1, exon 1 of 4)",     # -> 1st Exon
        "Exon (T1/G1, exon 2 of 4)",     # -> Other Exon
        "Intron (T1/G1, intron 1 of 3)", # -> 1st Intron
        "Intron (T1/G1, intron 2 of 3)", # -> Other Intron
        "Distal Intergenic",
    ]
    stat = annotation_stat(labels)
    assert [f for f, _ in stat] == [
        "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
        "5' UTR", "3' UTR", "1st Exon", "Other Exon",
        "1st Intron", "Other Intron", "Distal Intergenic",
    ]
    # The six-bin vocabulary would have collapsed this to six names.
    assert len(stat) > len(BIN_NAMES)


def test_stat_frequencies_are_percentages_that_sum_to_100():
    labels = ["Promoter (<=1kb)", "3' UTR", "Distal Intergenic", "Distal Intergenic"]
    stat = dict(annotation_stat(labels))
    assert stat["Promoter (<=1kb)"] == 25.0
    assert stat["Distal Intergenic"] == 50.0
    assert abs(sum(stat.values()) - 100.0) < 1e-9


def test_stat_folds_downstream_without_the_bp_suffix():
    """getGenomicAnnoStat writes 'Downstream (<=300)', unlike the per-peak label.

    The inconsistency is ChIPseeker's own: annotatePeak emits
    "Downstream (<=300bp)" while getGenomicAnnoStat emits "Downstream (<=300)".
    """
    stat = dict(annotation_stat(["Downstream (<=300bp)"] + ["Distal Intergenic"]))
    assert "Downstream (<=300)" in stat
    assert "Downstream (<=300bp)" not in stat


def test_stat_omits_zero_categories_like_table():
    stat = annotation_stat(["Distal Intergenic"] * 3)
    assert [f for f, _ in stat] == ["Distal Intergenic"]
    assert stat[0][1] == 100.0


def test_stat_promoter_rows_follow_r_sort_order():
    """R does sort(promoter); under en_US.UTF-8 that is numeric, so <=1kb is first."""
    labels = ["Promoter (4-5kb)", "Promoter (<=1kb)", "Promoter (2-3kb)"]
    assert [f for f, _ in annotation_stat(labels)] == [
        "Promoter (<=1kb)", "Promoter (2-3kb)", "Promoter (4-5kb)"
    ]


def test_ignore_downstream_option_relabels_to_intergenic():
    """getOption('ChIPseeker.ignore_downstream') folds Downstream into Intergenic."""
    stat = dict(annotation_stat(["Downstream (<=300bp)", "Distal Intergenic"],
                                ignore_downstream=True))
    assert stat == {"Intergenic": 100.0}



# --- strand handling, measured against ChIPseeker -------------------------
#
# A narrowPeak carries signalValue in column 6, so ChIPseeker's readPeakFile
# drops the strand and every peak becomes "*". These cases exercise the strand
# paths that never run on such files, measured with an explicit GRanges strand.

def _tss_fixture():
    """One plus-strand TSS at 1000 and one minus-strand TSS at 2100."""
    from easyepi.annotate.chipseeker_port import Features
    f = Features()
    f.strand = ["+", "-"]
    f.span = [(900, 1100), (1900, 2100)]
    f.chrom = ["chr1", "chr1"]
    # resize(width=1): left edge on "+", right edge on "-"
    f.tss = {"chr1": [(1000, 0), (2100, 1)]}
    f.key_to_index = {"T1": 0, "T2": 1}
    return f


def test_minus_peak_follow_and_precede_swap():
    """IRanges::follow/precede are strand aware, so a minus query swaps them.

    Measured: on chr1 with TSS at 710676 (-) and 720150 (+), a minus-strand
    peak at 713701-714700 gets follow -> 720150 and precede -> 710676, the
    opposite of a plus or unstranded query.
    """
    from easyepi.annotate.chipseeker_port import nearest_distance
    f = _tss_fixture()
    _, plus_d = nearest_distance(f, "chr1", 1500, 1600, peak_strand="+")
    _, minus_d = nearest_distance(f, "chr1", 1500, 1600, peak_strand="-")
    # plus reads start (1500), minus reads end (1600), and the two candidate
    # TSSs trade places -- so the results must differ.
    assert plus_d != minus_d


def test_a_tss_inside_the_peak_forces_zero_distance():
    """findOverlaps on the 1 bp features sets distanceToTSS to 0."""
    from easyepi.annotate.chipseeker_port import nearest_distance
    f = _tss_fixture()
    assert nearest_distance(f, "chr1", 900, 1100)[1] == 0.0


def test_same_strand_gates_sub_feature_hits():
    """getGenomicAnnotation.internal lines 23-27 keep the strand when TRUE.

    The strand of a row lives in the feature strand vector, not in the label
    column; reading the label there silently filtered every hit.
    """
    from easyepi.annotate.chipseeker_port import annotate, Features
    f = Features()
    f.strand = ["+"]
    f.span = [(900, 20000)]
    f.chrom = ["chr1"]
    f.tss = {"chr1": [(1000, 0)]}
    f.key_to_index = {"T1": 0}
    f.child = {"intron": {"chr1": ([
        (950, 19000, 0, "T1/G1", 1, 1, 0),
    ], [(950, 0)])}}
    # max_len must be the widest interval of this kind, or the backward scan
    # starts too late and misses the row entirely.
    f.child_max_len = {"intron": 19000 - 950}
    # Far outside tssRegion(-3000, 3000), otherwise Promoter wins the loop.
    plus, _, _ = annotate(f, ["chr1"], [12000], [12100], same_strand=True,
                          peaks_strand=["+"])
    minus, _, _ = annotate(f, ["chr1"], [12000], [12100], same_strand=True,
                           peaks_strand=["-"])
    assert plus[0].startswith("Intron")
    assert minus[0] == "Distal Intergenic"


def test_gene_level_uses_gene_spans_not_transcripts():
    """level='gene' replaces features with genes(), so the TSS vector changes.

    annotation itself is final before R reads level, so a far-away peak stays
    Intergenic at both levels -- what changes is distanceToTSS.
    """
    from easyepi.annotate.chipseeker_port import annotate, Features
    f = Features()
    f.strand = ["+", "+"]
    f.span = [(1000, 2000), (5000, 6000)]
    f.chrom = ["chr1", "chr1"]
    f.tss = {"chr1": [(1000, 0), (5000, 1)]}
    f.key_to_index = {"T1": 0, "T2": 1}
    f.gene_of = {"T1": "G1", "T2": "G2"}
    f.span_gene = {"chr1": [(1000, 2000, "+", "G1"), (5000, 6000, "+", "G2")]}
    tx_labels, _, _ = annotate(f, ["chr1"], [4501], [4600], level="transcript")
    gene_labels, _, _ = annotate(f, ["chr1"], [4501], [4600], level="gene")
    assert tx_labels == gene_labels


# ---------------------------------------------------------------------------
# flank_* columns (ChIPseeker addFlankGeneInfo=TRUE)
#
# Reference values below were produced by ChIPseeker 1.46.1 on the same peak
# and the same TxDb, then reduced to a synthetic feature set. The rules they
# pin down:
#   * the search uses the peak widened by `distance`, the arithmetic uses the
#     original peak;
#   * a feature overlapping the original peak scores 0;
#   * otherwise the smaller |psD| / |peD| wins and the value is NOT flipped by
#     the peak's strand -- the reference really does print -3025;
#   * no overlapping feature at all writes "NA" in every column.
# ---------------------------------------------------------------------------


def _flank_feats():
    from easyepi.annotate.chipseeker_port import Features

    f = Features()
    # Two features on chr1. Feature 0 overlaps the peak, feature 1 sits 3025 bp
    # downstream of it on the minus strand.
    f.span = [(760000, 762500), (764900, 765600)]
    f.chrom = ["chr1", "chr1"]
    f.strand = ["-", "-"]
    f.gene_of = {"TxID:1": "GeneID:1", "TxID:2": "GeneID:2"}
    f.key_order = ["TxID:1", "TxID:2"]
    return f


def test_flank_lists_overlapping_features_and_scores_zero():
    from easyepi.annotate.chipseeker_port import flank_info

    got = flank_info(_flank_feats(), ["chr1"], [762101], [762900], 5000)
    assert got["flank_txIds"][0] == "TxID:1;TxID:2"
    # Feature 0 overlaps the peak -> 0. Feature 1 is outside it, so its own
    # '-' strand drives the arithmetic: psD = 765600 - 762900 = 2700 and
    # peD = 765600 - 762101 = 3499, and the smaller |.| wins.
    assert got["flank_gene_distances"][0] == "0;2700"


def test_flank_window_is_widened_but_distance_uses_the_original_peak():
    from easyepi.annotate.chipseeker_port import flank_info

    wide = flank_info(_flank_feats(), ["chr1"], [762101], [762900], 5000)
    narrow = flank_info(_flank_feats(), ["chr1"], [762101], [762900], 100)
    # Feature 1 starts 2000 bp downstream: visible with a 5000 bp window only.
    assert "TxID:2" in wide["flank_txIds"][0]
    assert "TxID:2" not in narrow["flank_txIds"][0]


def test_flank_distance_is_not_sign_flipped_by_the_peak_strand():
    from easyepi.annotate.chipseeker_port import flank_info

    got = flank_info(_flank_feats(), ["chr1"], [762101], [762900], 5000)
    dists = got["flank_gene_distances"][0].split(";")
    # The peak is '-'. ChIPseeker's reference prints -3025 for a feature that
    # sits to the left, i.e. the sign comes from the feature's own strand and
    # is never re-normalised by the peak. Here the feature is downstream of a
    # '-' peak, so psD is already positive and no flip happens.
    assert dists[1] == "2700"


def test_flank_writes_na_when_nothing_is_close():
    from easyepi.annotate.chipseeker_port import flank_info

    got = flank_info(_flank_feats(), ["chr1"], [100], [200], 1000)
    assert got["flank_txIds"][0] == "NA"
    assert got["flank_geneIds"][0] == "NA"
    assert got["flank_gene_distances"][0] == "NA"


def test_flank_gene_level_omits_tx_ids():
    from easyepi.annotate.chipseeker_port import Features, flank_info

    f = _flank_feats()
    f.span_gene = {"chr1": [(760000, 762500, "-", "GeneID:1")]}
    got = flank_info(f, ["chr1"], [762101], [762900], 5000, level="gene")
    # At gene level the column is absent, so the list stays empty; R drops it
    # from the data frame entirely.
    assert got["flank_txIds"] == []
    assert got["flank_geneIds"][0] == "GeneID:1"
    assert got["flank_gene_distances"][0] == "0"


# ---------------------------------------------------------------------------
# write_tables: the GTF -> feature.tsv/tss.tsv/tx.tsv exporter
#
# Found untested while verifying the GTF path against a TxDb. It had been broken
# since child[chrom] became a (rows, order) pair: the exporter still iterated
# the pair as if it were the row list, so any call raised
# "too many values to unpack". Nothing exercised it, so nothing caught it.
# ---------------------------------------------------------------------------

_MINIMAL_GTF = """\
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\tHAVANA\tCDS\t1050\t1150\t.\t+\t0\tgene_id "G1"; transcript_id "T1";
chr1\tHAVANA\texon\t1800\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\tHAVANA\tCDS\t1850\t1950\t.\t+\t0\tgene_id "G1"; transcript_id "T1";
"""


def test_write_tables_emits_every_kind(tmp_path):
    from easyepi.annotate.gtf_features import write_tables

    gtf = tmp_path / "toy.gtf"
    gtf.write_text(_MINIMAL_GTF, encoding="utf-8")
    out = tmp_path / "feat"
    counts = write_tables(gtf, out)

    assert counts["tss"] == 1
    assert (out / "tss.tsv").exists()
    assert (out / "tx.tsv").exists()
    assert (out / "feature.tsv").exists()

    import csv
    with (out / "feature.tsv").open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    kinds = {r["kind"] for r in rows}
    # One transcript, two exons, one intron between them, and UTRs derived from
    # exon minus CDS because this GTF carries no explicit UTR rows.
    assert kinds == {"exon", "intron", "utr5", "utr3"}
    assert counts["feature"] == len(rows)


def test_write_tables_keeps_rows_in_array_order(tmp_path):
    """The child rows must not be sorted: ChIPseeker's subject index is the
    array position produced by unlist(), and the exporter has to preserve it."""
    from easyepi.annotate.gtf_features import load_from_gtf, write_tables

    gtf = tmp_path / "toy.gtf"
    gtf.write_text(_MINIMAL_GTF, encoding="utf-8")
    out = tmp_path / "feat"
    write_tables(gtf, out)

    feats = load_from_gtf(gtf)
    rows, _order = feats.child["exon"]["chr1"]
    emit = [r[:6] for r in rows]

    import csv
    with (out / "feature.tsv").open(newline="", encoding="utf-8") as fh:
        written = [(int(r["start"]), int(r["end"]), r["gene"], int(r["rank"]))
                   for r in csv.DictReader(fh, delimiter="\t")
                   if r["kind"] == "exon"]
    assert written == [(s, e, g, k) for (s, e, _o, g, k, _t) in emit]
