"""Compact input syntax: --peaks short forms and QC --input-dir scanning."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from easyepi.peaks.diffpeaks import parse_peak_specs
from easyepi.qc_inputs import resolve_inputs, scan_folder, split_values


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x")
    return path


# --------------------------------------------------------------------------
# split_values
# --------------------------------------------------------------------------
def test_split_values_flattens_and_strips():
    assert split_values(["a,b", " c "]) == ["a", "b", "c"]
    assert split_values(None) == []


# --------------------------------------------------------------------------
# --peaks: differential-peaks
# --------------------------------------------------------------------------
def test_peaks_explicit_form(tmp_path):
    bam = _touch(tmp_path / "a.bam")
    peak = _touch(tmp_path / "a.bed")
    bam2 = _touch(tmp_path / "b.bam")
    peak2 = _touch(tmp_path / "b.bed")
    got = parse_peak_specs([f"ctrl={bam}:{peak}", f"treat={bam2}:{peak2}"])
    assert got["ctrl"] == {"bam": str(bam), "peak": str(peak)}
    assert got["treat"] == {"bam": str(bam2), "peak": str(peak2)}


def test_peaks_bam_only_finds_peak_in_peak_dir(tmp_path):
    """`sample=BAM` pairs an explicit BAM with the peak file found by sample name."""
    bam = _touch(tmp_path / "bams" / "rep1.bam")
    _touch(tmp_path / "peaks" / "rep1.narrowPeak")
    bam2 = _touch(tmp_path / "bams" / "rep2.bam")
    _touch(tmp_path / "peaks" / "rep2.narrowPeak")
    got = parse_peak_specs([f"rep1={bam}", f"rep2={bam2}"],
                           peak_dir=tmp_path / "peaks")
    assert got["rep1"]["bam"] == str(bam)
    assert got["rep1"]["peak"].endswith("rep1.narrowPeak")
    assert got["rep2"]["peak"].endswith("rep2.narrowPeak")


def test_peaks_bare_labels_resolve_both_dirs(tmp_path):
    for name in ("ctrl1", "ctrl2"):
        _touch(tmp_path / "bams" / f"{name}.bam")
        _touch(tmp_path / "peaks" / f"{name}.narrowPeak")
    got = parse_peak_specs(["ctrl1,ctrl2"],
                           bam_dir=tmp_path / "bams",
                           peak_dir=tmp_path / "peaks")
    assert sorted(got) == ["ctrl1", "ctrl2"]
    assert got["ctrl1"]["bam"].endswith("ctrl1.bam")
    assert got["ctrl1"]["peak"].endswith("ctrl1.narrowPeak")


def test_peaks_missing_bam_reports_the_dir(tmp_path):
    _touch(tmp_path / "peaks" / "x.narrowPeak")
    with pytest.raises(FileNotFoundError, match="could not resolve"):
        parse_peak_specs(["x"], bam_dir=tmp_path / "bams", peak_dir=tmp_path / "peaks")


def test_peaks_missing_peak_reports_the_dir(tmp_path):
    _touch(tmp_path / "bams" / "x.bam")
    with pytest.raises(FileNotFoundError, match="could not resolve"):
        parse_peak_specs(["x"], bam_dir=tmp_path / "bams", peak_dir=tmp_path / "peaks")


def test_peaks_bare_label_without_dirs_is_an_error():
    with pytest.raises(FileNotFoundError, match="could not resolve"):
        parse_peak_specs(["ctrl1"])


def test_peaks_needs_two_samples(tmp_path):
    bam = _touch(tmp_path / "a.bam")
    peak = _touch(tmp_path / "a.bed")
    with pytest.raises(ValueError, match="at least two samples"):
        parse_peak_specs([f"only={bam}:{peak}"])


def test_peaks_rejects_a_missing_explicit_path(tmp_path):
    _touch(tmp_path / "a.bed")
    with pytest.raises(FileNotFoundError, match="BAM not found"):
        parse_peak_specs([f"s=/no/such.bam:{tmp_path / 'a.bed'}"])


def test_peaks_ambiguous_peak_file_is_an_error(tmp_path):
    _touch(tmp_path / "bams" / "x.bam")
    _touch(tmp_path / "peaks" / "x.narrowPeak")
    _touch(tmp_path / "peaks" / "x.bed")
    with pytest.raises(ValueError, match="ambiguous"):
        parse_peak_specs(["x"], bam_dir=tmp_path / "bams", peak_dir=tmp_path / "peaks")


def test_peaks_requires_at_least_one_entry():
    with pytest.raises(ValueError, match="at least one"):
        parse_peak_specs([])


# --------------------------------------------------------------------------
# QC --input-dir
# --------------------------------------------------------------------------
def test_scan_folder_groups_by_sample(tmp_path):
    for name in ("s1", "s2"):
        _touch(tmp_path / f"{name}.bam")
        _touch(tmp_path / f"{name}.bw")
        _touch(tmp_path / f"{name}.narrowPeak")
    got = scan_folder(tmp_path)
    assert sorted(got) == ["s1", "s2"]
    assert set(got["s1"]) == {"bam", "bigwig", "peak"}


def test_scan_folder_strips_bed_gz_completely(tmp_path):
    _touch(tmp_path / "s1.bed.gz")
    got = scan_folder(tmp_path)
    assert list(got) == ["s1"]
    assert got["s1"]["peak"].endswith("s1.bed.gz")


def test_scan_folder_ignores_bai_indexes(tmp_path):
    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.bam.bai")
    got = scan_folder(tmp_path)
    assert list(got) == ["s1"]


def test_scan_folder_rejects_two_files_of_the_same_kind(tmp_path):
    _touch(tmp_path / "s1.bw")
    _touch(tmp_path / "s1.bigWig")
    with pytest.raises(ValueError, match="two bigwig files"):
        scan_folder(tmp_path)


def test_scan_folder_rejects_a_missing_directory(tmp_path):
    with pytest.raises(NotADirectoryError, match="not a directory"):
        scan_folder(tmp_path / "nope")


def test_scan_folder_rejects_an_empty_directory(tmp_path):
    with pytest.raises(ValueError, match="no BAM/BigWig/peak"):
        scan_folder(tmp_path)


def _args(**kw):
    base = dict(bam=None, bigwig=None, peak=None, input_dir=None)
    base.update(kw)
    return SimpleNamespace(**base)


def test_resolve_inputs_accepts_comma_lists(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    got = resolve_inputs(_args(bam=[f"{a},{b}"]))
    assert sorted(got) == ["a", "b"]
    assert got["a"]["bam"] == str(a)


def test_resolve_inputs_pairs_by_sample_name_not_position(tmp_path):
    """The old code paired bigwig[i] with bam[i]; names must decide now."""
    bam_a = _touch(tmp_path / "a.bam")
    bam_b = _touch(tmp_path / "b.bam")
    bw_a = _touch(tmp_path / "a.bw")
    bw_b = _touch(tmp_path / "b.bw")
    # bigwig order reversed relative to bam order
    got = resolve_inputs(_args(bam=[str(bam_a), str(bam_b)],
                               bigwig=[str(bw_b), str(bw_a)]))
    assert got["a"]["bigwig"] == str(bw_a)
    assert got["b"]["bigwig"] == str(bw_b)


def test_resolve_inputs_merges_folder_then_flags(tmp_path):
    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.narrowPeak")
    override = _touch(tmp_path / "other" / "s1.bw")
    got = resolve_inputs(_args(input_dir=tmp_path, bigwig=[str(override)]))
    assert got["s1"]["bam"].endswith("s1.bam")
    assert got["s1"]["bigwig"] == str(override)


def test_resolve_inputs_rejects_missing_flag_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="not found"):
        resolve_inputs(_args(bam=[str(tmp_path / "nope.bam")]))


def test_resolve_inputs_rejects_two_tracks_of_the_same_kind(tmp_path):
    a = _touch(tmp_path / "s.bw")
    b = _touch(tmp_path / "s.bigWig")
    with pytest.raises(ValueError, match="two bigwig files"):
        resolve_inputs(_args(bigwig=[f"{a},{b}"]))


def test_resolve_inputs_requires_something():
    with pytest.raises(ValueError, match="no inputs"):
        resolve_inputs(_args())


def test_qc_cli_exposes_input_dir():
    from easyepi import cli
    parser = cli.build_parser()
    args = parser.parse_args(["qc", "--input-dir", "/data/tracks", "--out", "r.json"])
    assert args.input_dir == "/data/tracks"


def test_differential_peaks_cli_exposes_dirs():
    from easyepi import cli
    parser = cli.build_parser()
    args = parser.parse_args([
        "differential-peaks", "--peaks", "ctrl1,ctrl2",
        "--bam-dir", "/b", "--peak-dir", "/p",
        "--group-a", "ctrl1", "--out", "d.tsv"])
    assert args.bam_dir == "/b" and args.peak_dir == "/p"
