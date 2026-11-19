"""Coordinate contract: every parser and writer must honour its own format.

The plan's rule is that each file format keeps its own standard rather than forcing one
convention everywhere:

    BED / narrowPeak / broadPeak : 0-based, half-open [start, end)
    GTF / GFF3                   : 1-based, closed [start, end]
    BigWig / bedGraph            : 0-based, half-open at the API level

The failure mode this guards against is a *double* conversion -- a parser that already
normalises, followed by a writer that normalises again -- and off-by-one drift at
truncated contig ends. Reading the code cannot prove absence; these tests exercise the
actual round trips on deliberately awkward coordinates.
"""

from __future__ import annotations

import gzip

import numpy as np
import pytest

from easyepi.annotate import read_intervals
from easyepi.annotate.gtf_reader import read_transcripts
from easyepi.gc_background import write_bed
from easyepi.peaks.diffpeaks import read_peaks

# A single-base feature at the very start of the contig is the case where a stray -1
# turns into a negative coordinate.
EDGE_INTERVALS = [
    ("chr1", 0, 1),        # single base at origin
    ("chr1", 0, 100),      # starts at origin
    ("chr1", 99, 100),     # single base, not at origin
    ("chr1", 100, 100),    # zero-length: must be dropped, not silently kept
    ("chr1", 100, 101),    # adjacent single base
]


def _write_bed(path, rows, header=""):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        if header:
            fh.write(header)
        for chrom, s, e in rows:
            fh.write(f"{chrom}\t{s}\t{e}\n")


def test_bed_reading_is_verbatim_zero_based(tmp_path):
    """A BED start is already 0-based; the parser must not shift it."""
    p = tmp_path / "a.bed"
    _write_bed(p, [("chr1", 0, 1), ("chr1", 99, 100), ("chr1", 1000, 1500)],
               header="track name=x\n")
    got = read_peaks(p)
    assert got == [("chr1", 0, 1), ("chr1", 99, 100), ("chr1", 1000, 1500)]


def test_bed_reader_drops_zero_length_but_keeps_single_base(tmp_path):
    """Length 0 is not a feature; length 1 is."""
    p = tmp_path / "a.bed"
    _write_bed(p, [("chr1", 100, 100), ("chr1", 100, 101)])
    assert read_peaks(p) == [("chr1", 100, 101)]


def test_bed_round_trip_preserves_interval_set(tmp_path):
    """write_bed then read_peaks must be the identity on 0-based half-open input."""
    src = [(c, s, e) for c, s, e in EDGE_INTERVALS if e > s]
    out = write_bed(src, tmp_path / "rt.bed")
    assert read_peaks(out) == src


def test_gzipped_bed_matches_plain(tmp_path):
    rows = [("chr1", 0, 1), ("chr1", 500, 750)]
    plain = tmp_path / "p.bed"
    _write_bed(plain, rows)
    gz = tmp_path / "p.bed.gz"
    with gzip.open(gz, "wt", encoding="utf-8", newline="") as fh:
        for c, s, e in rows:
            fh.write(f"{c}\t{s}\t{e}\n")
    assert read_peaks(gz) == read_peaks(plain)


def test_gtf_is_one_based_and_closed_so_start_shifts_by_one(tmp_path):
    """GTF col4 is 1-based; the internal interval must start one base earlier.

    A GTF feature at 1..1 is the first base of the contig. If the parser stored col4
    verbatim it would claim 1..1 in 0-based terms -- the second base -- and a feature
    at the origin would appear to start at 1.
    """
    p = tmp_path / "a.gtf"
    p.write_text(
        'chr1\tsrc\texon\t1\t1\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    txs = read_transcripts(str(p))
    rec = txs[0]
    assert rec["exons"] == [(0, 1), (99, 200)]
    # tss for a plus-strand transcript is the leftmost exon start, already 0-based
    assert rec["tss"] == 0


def test_gff3_and_gtf_agree_on_the_same_model(tmp_path):
    """The two formats describe the same interval; both must yield identical records."""
    gtf = tmp_path / "m.gtf"
    gtf.write_text(
        'chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tsrc\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    gff = tmp_path / "m.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\texon\t101\t200\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tsrc\texon\t301\t400\t.\t+\t.\tID=e2;Parent=t1\n"
        "chr1\tsrc\tmRNA\t101\t400\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tsrc\tgene\t101\t400\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )
    a = read_transcripts(str(gtf))[0]
    b = read_transcripts(str(gff))[0]
    assert a["exons"] == b["exons"] == [(100, 200), (300, 400)]
    assert a["tss"] == b["tss"] == 100
    assert a["strand"] == b["strand"] == "+"


def test_minus_strand_tss_is_the_rightmost_base(tmp_path):
    """On the minus strand the TSS is the highest coordinate, not the lowest."""
    p = tmp_path / "a.gtf"
    p.write_text(
        'chr1\tsrc\texon\t101\t200\t.\t-\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tsrc\texon\t301\t400\t.\t-\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    rec = read_transcripts(str(p))[0]
    assert rec["tss"] == 400, "minus-strand TSS must be the exon end, in 0-based terms"


def test_read_intervals_is_zero_based_half_open(tmp_path):
    """The shared interval reader must agree with the BED convention it reads."""
    p = tmp_path / "a.bed"
    _write_bed(p, [("chr1", 0, 1), ("chr1", 99, 100)])
    chrom, start, end, _strand = read_intervals(p)
    assert list(chrom) == ["chr1", "chr1"]
    assert list(start) == [0, 99]
    assert list(end) == [1, 100]
    assert np.all(np.asarray(end) > np.asarray(start))
