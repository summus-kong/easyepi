"""ChromHMM and HOMER run for real, and their input contracts are pinned.

Both tools are installed in the `easyepi` env, so these tests execute them rather
than mocking. Running them by hand first surfaced four input requirements that
are easy to get wrong and that nothing documented:

ChromHMM `BinarizeBed`
  1. the third positional argument (cellmarkfiletable) is mandatory, and the
     failure modes are unhelpful: with no arguments ChromHMM prints usage and
     exits 0; with some of them it reads `-b` as the lengths file and throws
     FileNotFoundException naming the flag
  2. the file column of that table must be a BARE FILENAME, not a path: ChromHMM
     prefixes it with the bed directory, so a full path yields
     ``<beddir>/<beddir>/<file>.bed`` and a FileNotFoundException
  3. the mark name in column 2 must match the file's contents, else
     "did not find data for ... treating as missing" and no bins are written
  4. BED files need a strand column (6 columns), otherwise
     "strand column expected, but not found"

HOMER
  5. ``run_homer_motifs`` requires an explicit background FASTA in FASTA mode;
     ``findMotifsGenome.pl`` needs an installed genome under
     ``share/homer/data/genomes``, and none is installed here, so FASTA mode is
     the only usable path in this environment.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

CHROMHMM_JAR = sorted(
    Path("/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/share")
    .glob("chromhmm-*/ChromHMM.jar")
)
HOMER_DIR = Path("/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/share/homer")
HUMAN_FAI = Path("/data1/GaoLab/houoweibo/reference/human/genome.fa.fai")
HUMAN_FA = Path("/data1/GaoLab/houoweibo/reference/human/genome.fa")

has_java = shutil.which("java") is not None
requires_chromhmm = pytest.mark.skipif(
    not (CHROMHMM_JAR and has_java and HUMAN_FAI.exists()),
    reason="ChromHMM jar / java / human genome index not available")
requires_homer = pytest.mark.skipif(
    not (HOMER_DIR / "bin/findMotifs.pl").exists(),
    reason="HOMER not installed")


def _run(cmd, cwd=None):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode(errors="replace")


def _write_marks(tmp_path, marks):
    """Six-column BEDs (strand required) plus the ChromHMM table."""
    beddir = tmp_path / "bed"
    beddir.mkdir()
    rows = []
    for i, (name, base) in enumerate(marks):
        path = beddir / ("%s.bed" % name)
        path.write_text("".join(
            "chr1\t%d\t%d\t.\t0\t+\n" % (s, s + 2000)
            for s in (base + k * 10_000 for k in range(200))))
        # bare filename, not the full path
        rows.append("ct%d\t%s\t%s.bed\n" % (i // 2, name, name))
    (tmp_path / "cellmarkfiletable.txt").write_text("".join(rows))
    return beddir


@requires_chromhmm
def test_chromhmm_binarize_writes_binary_tracks(tmp_path):
    """BinarizeBed produces one *_binary.txt per (cell, chromosome)."""
    genome = tmp_path / "genome.txt"
    lines = [ln.split("\t") for ln in HUMAN_FAI.read_text().splitlines()]
    genome.write_text("".join("%s\t%s\n" % (c[0], c[1]) for c in lines))

    beddir = _write_marks(tmp_path, [("H3K4me3", 1_000_000), ("H3K27ac", 5_000_000)])
    out = tmp_path / "bin"
    code, text = _run(["java", "-mx2000M", "-jar", CHROMHMM_JAR[0], "BinarizeBed",
                       "-b", "200", genome, beddir,
                       tmp_path / "cellmarkfiletable.txt", out])
    assert code == 0, text[-1500:]
    # ChromHMM exits 0 on a usage error, so assert on the artifacts instead.
    produced = sorted(p.name for p in out.iterdir() if p.is_file())
    assert produced, "no binary output produced:\n" + text[-1500:]
    assert all(n.endswith("_binary.txt") for n in produced), produced


@requires_chromhmm
def test_chromhmm_ignores_a_full_path_in_the_cellmark_table(tmp_path):
    """A path in column 3 is silently concatenated with the bed dir.

    Pinned because the resulting FileNotFoundException names a doubled path that
    looks like a ChromHMM bug rather than a caller mistake.
    """
    genome = tmp_path / "genome.txt"
    genome.write_text("chr1\t50000000\n")
    beddir = _write_marks(tmp_path, [("H3K4me3", 1_000_000)])
    table = tmp_path / "bad.txt"
    table.write_text("A\tH3K4me3\t%s\n" % (beddir / "H3K4me3.bed"))
    code, text = _run(["java", "-mx2000M", "-jar", CHROMHMM_JAR[0], "BinarizeBed",
                       "-b", "200", genome, beddir, table, tmp_path / "bin"])
    assert code != 0
    assert "FileNotFoundException" in text
    # the doubled prefix is the tell
    assert str(beddir) + str(beddir) in text.replace("\n", ""), text[:600]


@requires_chromhmm
def test_chromhmm_rejects_a_bed_without_a_strand_column(tmp_path):
    genome = tmp_path / "genome.txt"
    genome.write_text("chr1\t50000000\n")
    beddir = tmp_path / "bed"
    beddir.mkdir()
    (beddir / "H3K4me3.bed").write_text("".join(
        "chr1\t%d\t%d\n" % (1_000_000 + k * 10_000, 1_002_000 + k * 10_000)
        for k in range(50)))
    table = tmp_path / "t.txt"
    table.write_text("A\tH3K4me3\tH3K4me3.bed\n")
    code, text = _run(["java", "-mx2000M", "-jar", CHROMHMM_JAR[0], "BinarizeBed",
                       "-b", "200", genome, beddir, table, tmp_path / "bin"])
    assert code != 0
    assert "strand column expected" in text, text[:600]


def test_chromhmm_with_no_arguments_prints_usage_and_exits_zero(tmp_path):
    """Bare `BinarizeBed` shows the usage banner and still exits 0.

    Callers cannot rely on the status code to detect a malformed command line.
    """
    if not (CHROMHMM_JAR and has_java):
        pytest.skip("ChromHMM not available")
    code, text = _run(["java", "-mx2000M", "-jar", CHROMHMM_JAR[0], "BinarizeBed"])
    assert text.lstrip().startswith("usage BinarizeBed"), text[:300]
    assert code == 0, "measured: the usage banner exits 0"


def test_chromhmm_partial_arguments_swallow_the_flag_as_a_filename(tmp_path):
    """Dropping the table makes `-b` be read as the lengths file.

    Pinned because the message (`FileNotFoundException: -b`) points at a flag,
    not at the missing positional argument, which is misleading to diagnose.
    """
    if not (CHROMHMM_JAR and has_java):
        pytest.skip("ChromHMM not available")
    genome = tmp_path / "genome.txt"
    genome.write_text("chr1\t50000000\n")
    beddir = tmp_path / "bed"
    beddir.mkdir()
    code, text = _run(["java", "-mx2000M", "-jar", CHROMHMM_JAR[0], "BinarizeBed",
                       "-b", "200", genome, beddir])
    assert code != 0
    assert "FileNotFoundException: -b" in text, text[:400]
    assert not (tmp_path / "out").exists()


@requires_homer
def test_homer_fasta_mode_requires_a_background(tmp_path):
    """`run_homer_motifs` refuses FASTA mode without an explicit background."""
    import sys
    sys.path.insert(0, "/data1/GaoLab/houoweibo/easyepi_work/easyepiv1_relay")
    from easyepi.motifs.external_tools import run_homer_motifs

    fa = tmp_path / "t.fa"
    fa.write_text(">p1\nACGTACGTACGTACGTACGT\n>p2\nTTTTAAAACCCCGGGGTTTT\n")
    with pytest.raises(ValueError, match="background"):
        run_homer_motifs(str(fa), None, str(tmp_path / "out"),
                         homer_dir=str(HOMER_DIR))


@requires_homer
def test_homer_fasta_mode_runs_with_a_background(tmp_path):
    """With both FASTAs, findMotifs.pl runs and writes its result dirs."""
    import sys
    sys.path.insert(0, "/data1/GaoLab/houoweibo/easyepi_work/easyepiv1_relay")
    from easyepi.motifs.external_tools import run_homer_motifs

    fa = tmp_path / "t.fa"
    bg = tmp_path / "bg.fa"
    fa.write_text("".join(
        ">t%d\nGTTGACCTCAGCTGCTGAGCT\n" % i for i in range(30)))
    bg.write_text("".join(
        ">b%d\nACACACACACACACACACAC\n" % i for i in range(30)))
    out = tmp_path / "out"
    res = run_homer_motifs(str(fa), str(bg), str(out),
                           k=6, homer_dir=str(HOMER_DIR))
    assert res["backend"] == "homer" and res["available"] is True
    assert out.exists()
    # HOMER writes homerResults / knownResults when it got far enough to scan.
    subdirs = {p.name for p in out.iterdir() if p.is_dir()}
    assert subdirs & {"homerResults", "knownResults"}, subdirs
