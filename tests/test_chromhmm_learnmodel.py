"""ChromHMM LearnModel runs end to end, and recovers the planted state structure.

`BinarizeBed` alone was already pinned; this covers the training stage, which the
first attempt could not reach. Two corrections were needed and both match what
`easyepi/specific/histone/chromhmm.py` actually emits:

  * `-peaks` (easyepi passes it; ChIP marks are peak calls, not raw signal)
  * two marks with **disjoint** occupancy, so LearnModel has real state
    structure. With one mark it aborts:

        INFORMATION initialization strategy can only support 1 states

Measured on ChromHMM 1.28: 11 model outputs, 9801 segments over chr1, all three
states used, and the emission matrix separates the two marks cleanly --
state "2" emits H3K27ac at ~1.0 and H3K4me3 at ~1e-222, state "3" the reverse,
state "1" emits neither. That structure is what the assertions below check: the
point is not a golden number, it is that training recovers the signal the fixture
planted.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

ENV = Path("/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi")
JARS = sorted((ENV / "share").glob("chromhmm-*/ChromHMM.jar"))
FAI = Path("/data1/GaoLab/houoweibo/reference/human/genome.fa.fai")

requires_chromhmm = pytest.mark.skipif(
    not (JARS and shutil.which("java") and FAI.exists()),
    reason="ChromHMM jar / java / human genome index not available")

#: Disjoint occupancy across chr1: K4me3 in the first block, K27ac in the second.
MARKS = (("H3K4me3", 1_000_000, 20_000_000),
         ("H3K27ac", 30_000_000, 60_000_000))


def _run(cmd):
    p = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode(errors="replace")


def _fixture(tmp_path):
    genome = tmp_path / "genome.txt"
    rows = [ln.split("\t") for ln in FAI.read_text().splitlines()]
    genome.write_text("".join("%s\t%s\n" % (r[0], r[1]) for r in rows))

    beddir = tmp_path / "input_beds"
    beddir.mkdir()
    for name, lo, hi in MARKS:
        with open(beddir / ("%s.bed" % name), "w") as f:
            s = lo
            while s < hi:
                f.write("chr1\t%d\t%d\t.\t0\t+\n" % (s, s + 2000))
                s += 10_000
    table = tmp_path / "cellmarkfiletable.tsv"
    table.write_text("".join("ct\t%s\t%s.bed\n" % (n, n) for n, _, _ in MARKS))
    return genome, beddir, table


def _read_emissions(path):
    """-> (mark names, {state: {mark: probability}})"""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    header = lines[0].split("\t")
    marks = header[1:]
    table = {}
    for ln in lines[1:]:
        cols = ln.split("\t")
        table[cols[0]] = {m: float(v) for m, v in zip(marks, cols[1:])}
    return marks, table


@requires_chromhmm
def test_learnmodel_recovers_the_two_planted_marks(tmp_path):
    genome, beddir, table = _fixture(tmp_path)

    code, text = _run(["java", "-mx4000M", "-jar", JARS[0], "BinarizeBed",
                       "-peaks", "-b", "200", genome, beddir, table,
                       tmp_path / "binarized"])
    assert code == 0, text[-1200:]
    binary = tmp_path / "binarized"
    assert list(binary.iterdir()), "binarization produced nothing:\n" + text[-1200:]

    model = tmp_path / "model"
    code, text = _run(["java", "-mx4000M", "-jar", JARS[0], "LearnModel",
                       "-b", "200", "-p", "2", "-s", "1",
                       binary, model, "3", "hg38"])
    assert code == 0, text[-1200:]

    # The training artifacts, not just the exit code.
    for name in ("model_3.txt", "emissions_3.txt", "transitions_3.txt",
                 "ct_3_segments.bed"):
        assert (model / name).exists(), "missing %s (got %s)" % (
            name, sorted(p.name for p in model.iterdir()))

    # Every requested state is actually used by the segmentation.
    seg = (model / "ct_3_segments.bed").read_text().splitlines()
    assert seg, "empty segmentation"
    states_used = {ln.split("\t")[3] for ln in seg if len(ln.split("\t")) > 3}
    assert states_used == {"E1", "E2", "E3"}, states_used

    # The emission matrix separates the two marks: one state is dominated by each.
    marks, emis = _read_emissions(model / "emissions_3.txt")
    assert set(marks) == {"H3K4me3", "H3K27ac"}, marks
    k4_max = max(v["H3K4me3"] for v in emis.values())
    k27_max = max(v["H3K27ac"] for v in emis.values())
    assert k4_max > 0.99, "no state emits H3K4me3 strongly: %r" % emis
    assert k27_max > 0.99, "no state emits H3K27ac strongly: %r" % emis
    # ... and they are different states, which is the structure being recovered.
    k4_state = max(emis, key=lambda s: emis[s]["H3K4me3"])
    k27_state = max(emis, key=lambda s: emis[s]["H3K27ac"])
    assert k4_state != k27_state, "one state absorbs both marks: %r" % emis
    # Neither mark-absent state should emit either mark strongly.
    other = [s for s in emis if s not in (k4_state, k27_state)]
    assert other, "expected a background state"
    for s in other:
        assert emis[s]["H3K4me3"] < 0.5 and emis[s]["H3K27ac"] < 0.5, (s, emis[s])


@requires_chromhmm
def test_learnmodel_fails_on_a_single_mark_binarization(tmp_path):
    """Pinned because the message is about states, not about the input.

        ... the INFORMATION initialization strategy can only support 2 states

    A caller who binarizes one mark and asks for 3 states sees a complaint about
    state count; the actual cause is that one mark has only present/absent
    structure to find. Measured with ChromHMM 1.28: the ceiling reported is 2.
    """
    genome = tmp_path / "genome.txt"
    genome.write_text("chr1\t50000000\n")
    beddir = tmp_path / "input_beds"
    beddir.mkdir()
    with open(beddir / "H3K4me3.bed", "w") as f:
        for k in range(100):
            s = 1_000_000 + k * 10_000
            f.write("chr1\t%d\t%d\t.\t0\t+\n" % (s, s + 2000))
    table = tmp_path / "t.tsv"
    table.write_text("ct\tH3K4me3\tH3K4me3.bed\n")

    binary = tmp_path / "binarized"
    code, text = _run(["java", "-mx4000M", "-jar", JARS[0], "BinarizeBed",
                       "-peaks", "-b", "200", genome, beddir, table, binary])
    assert code == 0, text[-800:]

    code, text = _run(["java", "-mx4000M", "-jar", JARS[0], "LearnModel",
                       "-b", "200", "-p", "1", "-s", "1",
                       binary, tmp_path / "model", "3", "hg38"])
    assert code != 0, "expected a single-mark model to fail"
    assert "initialization strategy can only support" in text, text[-800:]
    # The reported ceiling must be lower than requested, whichever number it is.
    import re
    m = re.search(r"can only support (\d+) states", text)
    assert m, "state ceiling not reported: %s" % text[-500:]
    assert int(m.group(1)) < 3, "asked for 3 states, tool allowed %s" % m.group(1)
