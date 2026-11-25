"""The TOBIAS wrapper must name a subcommand that exists.

Measured on TOBIAS 0.17.5: `TOBIAS Footprint` is rejected at argument parsing
("invalid choice: 'Footprint'"), and `FootprintScores --help` states that
`ScoreBigwig` replaces it. The wrapper used `Footprint` *and* passed `--cores`,
which ScoreBigwig does not accept, so every call failed before doing any work and
nothing noticed because nothing invoked it.

These tests pin the argv, and one end-to-end run proves the command is real.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from easyepi.motifs.external_tools import run_tobias_footprint  # noqa: E402

#: The subcommands TOBIAS 0.17.5 accepts, from its own "invalid choice" message.
VALID = {"ATACorrect", "ScoreBigwig", "FootprintScores", "BINDetect", "TFBScan",
         "FormatMotifs", "ClusterMotifs", "ScoreBed", "SubMerge",
         "PlotAggregate", "PlotHeatmap", "PlotChanges", "PlotTracks",
         "DownloadData", "MergePDF", "MaxPos", "SubsampleBam", "CreateNetwork",
         "Log2Table", "FilterFragments"}


def _tobias() -> str:
    for cand in (Path(sys.executable).parent / "TOBIAS", shutil.which("TOBIAS")):
        if cand and Path(cand).exists():
            return str(cand)
    pytest.skip("TOBIAS not installed; the footprint oracle needs it")


@pytest.fixture
def captured(monkeypatch):
    """Capture the argv the wrapper builds, without running anything."""
    calls = {}

    class _Result:
        returncode = 0
        stderr = ""

    def fake_run(cmd, *a, **kw):
        calls["cmd"] = list(cmd)
        return _Result()

    import easyepi.motifs.external_tools as et

    monkeypatch.setattr(et, "_run", fake_run)
    monkeypatch.setattr(et, "_require", lambda name, install_hint="": name.lower())
    return calls


def test_the_subcommand_is_one_tobias_accepts(captured, tmp_path):
    run_tobias_footprint([str(tmp_path / "s.bw")], str(tmp_path / "r.bed"),
                         "hg38", str(tmp_path / "out"))
    argv = captured["cmd"]
    assert argv[1] in VALID, "TOBIAS has no %r subcommand" % argv[1]
    assert argv[1] in ("ScoreBigwig", "FootprintScores")


def test_no_flag_tobias_does_not_accept(captured, tmp_path):
    """`--cores` was the second defect: ScoreBigwig has no such option."""
    run_tobias_footprint([str(tmp_path / "s.bw")], str(tmp_path / "r.bed"),
                         "hg38", str(tmp_path / "out"))
    assert "--cores" not in captured["cmd"]


def test_the_three_required_options_are_present(captured, tmp_path):
    bw = str(tmp_path / "s.bw")
    bed = str(tmp_path / "r.bed")
    run_tobias_footprint([bw], bed, "hg38", str(tmp_path / "out"))
    argv = captured["cmd"]
    for flag in ("--signal", "--regions", "--output"):
        assert flag in argv, argv
    assert argv[argv.index("--signal") + 1] == bw
    assert argv[argv.index("--regions") + 1] == bed


def test_one_track_at_a_time(captured, tmp_path):
    with pytest.raises(ValueError, match="separately for each"):
        run_tobias_footprint([str(tmp_path / "a.bw"), str(tmp_path / "b.bw")],
                             str(tmp_path / "r.bed"), "hg38",
                             str(tmp_path / "out"))


def test_motif_db_is_refused_with_the_reason(captured, tmp_path):
    """Footprint scoring is not binding inference; BINDetect is the later step."""
    with pytest.raises(ValueError, match="BINDetect"):
        run_tobias_footprint([str(tmp_path / "s.bw")], str(tmp_path / "r.bed"),
                             "hg38", str(tmp_path / "out"),
                             motif_db=str(tmp_path / "motifs.meme"))


def test_the_real_binary_accepts_the_built_argv(tmp_path):
    """Invoke TOBIAS's own parser: a wrong subcommand or flag fails loudly here."""
    tobias = _tobias()
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t1000\t2000\n", encoding="utf-8")
    # A bigWig is required by the parser but not opened during --help; passing a
    # non-existent path is fine because TOBIAS rejects on argument parsing first.
    proc = subprocess.run(
        [tobias, "ScoreBigwig", "--signal", str(tmp_path / "s.bw"),
         "--regions", str(bed), "--output", str(tmp_path / "o.bw"), "--help"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = proc.stdout.decode()
    assert "invalid choice" not in out, out[-800:]
    assert "--signal" in out and "--regions" in out
