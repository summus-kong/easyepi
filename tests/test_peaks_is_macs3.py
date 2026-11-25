"""`easyepi peaks` is MACS3, so its output must be MACS3's output.

The command is not a reimplementation: `peaks_cli.run_peaks_command` imports
`MACS3.Commands.callpeak_cmd` and calls its `run()` on a Namespace built from
MACS3's own argument parser. Two consequences are tested here:

  * whatever `--macs-arg` says reaches MACS3, or the call refuses to proceed --
    there is no third outcome where the option is accepted and ignored;
  * a real run produces the same narrowPeak bytes as invoking `macs3 callpeak`
    directly with the same options, because it is executing the same engine.

The second test is the reason the first matters: if a passthrough were dropped,
the two outputs would differ and this test would catch it.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("pysam")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

CHROM = "chr1"
CHROM_LEN = 120_000
READ_LEN = 50


def _has_macs3() -> bool:
    import importlib.util
    return importlib.util.find_spec("MACS3") is not None


def _macs3_binary() -> str | None:
    cand = Path(sys.executable).parent / "macs3"
    if cand.exists():
        return str(cand)
    return shutil.which("macs3")


needs_macs3 = pytest.mark.skipif(not _has_macs3(), reason="MACS3 not installed")


def _fixture(tmp: Path) -> Path:
    """A paired-end BAM with two enriched regions over background."""
    import pysam

    rng = np.random.default_rng(5)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    in_peak = np.zeros(CHROM_LEN, dtype=bool)
    peaks = [(30_000, 800), (80_000, 600)]
    for centre, half in peaks:
        in_peak[centre - half:centre + half] = True

    starts = list(rng.choice(np.flatnonzero(~in_peak), size=4000, replace=False))
    for centre, half in peaks:
        starts.extend(rng.integers(centre - half, centre + half, 3000).tolist())

    recs = []
    for i, start in enumerate(starts):
        frag = 180
        for flag, rs in ((99, start), (147, start + frag - READ_LEN)):
            a = pysam.AlignedSegment()
            a.query_name = "r%d" % i
            a.query_sequence = "A" * READ_LEN
            a.flag = flag
            a.reference_id = 0
            a.reference_start = int(rs)
            a.mapping_quality = 60
            a.cigar = ((0, READ_LEN),)
            a.next_reference_id = 0
            a.next_reference_start = int(start)
            a.template_length = frag if flag == 99 else -frag
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            recs.append(a)

    bam = tmp / "sample.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as fh:
        for r in sorted(recs, key=lambda x: (x.reference_start, x.flag)):
            fh.write(r)
    pysam.index(str(bam))
    return bam


# ---------------------------------------------------------------- passthrough

class _FakeParser:
    """Stands in for MACS3's parser: accepts nothing it was not told about."""

    def __init__(self, accepted):
        self.accepted = set(accepted)

    def parse_args(self, argv):
        ns = Namespace()
        i = 0
        while i < len(argv):
            tok = argv[i]
            if tok.startswith("--"):
                key = tok[2:].replace("-", "_")
                if key not in self.accepted:
                    raise SystemExit(2)
                if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                    setattr(ns, key, argv[i + 1])
                    i += 2
                    continue
                setattr(ns, key, True)
            i += 1
        return ns


def test_passthrough_reaches_the_namespace_when_the_parser_is_available(monkeypatch):
    from easyepi import macs_options

    accepted = {"tfile", "f", "g", "n", "outdir", "q", "keep_dup",
                "call_summits", "bw"}
    monkeypatch.setattr(macs_options, "_call_argparser_factory",
                        lambda argv: (_FakeParser(accepted), None))

    options, prov = macs_options.build_namespace(
        treatment=["a.bam"], out_dir=Path("/tmp/x"), name="s", genome="hs",
        fmt="BAMPE", broad=False, control=None, qvalue=0.01, extsize=None,
        shift=0, keep_dup="auto",
        extra_args=["--call-summits", "--bw 300"])

    assert prov["namespace_from"] == "macs3-argparse"
    assert options.call_summits is True
    assert options.bw == "300"
    assert "--call-summits" in prov["passthrough"]


def test_passthrough_is_refused_when_the_parser_cannot_be_read(monkeypatch):
    """The fallback Namespace cannot honour extra options, so it must not pretend."""
    from easyepi import macs_options

    monkeypatch.setattr(macs_options, "_call_argparser_factory",
                        lambda argv: (None, "no macs3 console script"))

    with pytest.raises(ValueError, match="silently ignored"):
        macs_options.build_namespace(
            treatment=["a.bam"], out_dir=Path("/tmp/x"), name="s", genome="hs",
            fmt="BAMPE", broad=False, control=None, qvalue=0.01, extsize=None,
            shift=0, keep_dup="auto", extra_args=["--call-summits"])


def test_no_passthrough_still_works_on_the_fallback(monkeypatch):
    """Without extra options the fallback is still usable: nothing is dropped."""
    from easyepi import macs_options

    monkeypatch.setattr(macs_options, "_call_argparser_factory",
                        lambda argv: (None, "no macs3 console script"))

    options, prov = macs_options.build_namespace(
        treatment=["a.bam"], out_dir=Path("/tmp/x"), name="s", genome="hs",
        fmt="BAMPE", broad=True, control=None, qvalue=0.01, extsize=200,
        shift=0, keep_dup="auto", extra_args=None)

    assert prov["namespace_from"] == "fallback"
    assert options.format == "BAMPE"
    assert options.broad is True
    assert options.extsize == 200


def test_a_rejected_passthrough_names_the_tokens(monkeypatch):
    from easyepi import macs_options

    monkeypatch.setattr(macs_options, "_call_argparser_factory",
                        lambda argv: (_FakeParser({"tfile"}), None))

    with pytest.raises(ValueError) as exc:
        macs_options.build_namespace(
            treatment=["a.bam"], out_dir=Path("/tmp/x"), name="s", genome="hs",
            fmt="BAMPE", broad=False, control=None, qvalue=0.01, extsize=None,
            shift=0, keep_dup="auto", extra_args=["--not-a-macs-option"])
    msg = str(exc.value)
    assert "not-a-macs-option" in msg
    assert "macs3 callpeak --help" in msg


# ------------------------------------------------------------- end to end

@needs_macs3
def test_the_run_produces_macs3s_own_peak_file(tmp_path):
    """Same engine, same input, same options -> the same peak coordinates."""
    from easyepi.peaks_cli import run_peaks_command

    binary = _macs3_binary()
    if binary is None:
        pytest.skip("macs3 console script not found next to the interpreter")

    bam = _fixture(tmp_path)
    ours = tmp_path / "ours"
    args = Namespace(bam=[str(bam)], out=ours, assay="atac", format="BAMPE",
                     genome="hs", qvalue=0.01, extsize=None, shift_size=None,
                     broad=False, keep_dup_all=False, macs_arg=None,
                     control=None, verbose=False, quiet=True, log_file=None)
    assert run_peaks_command(args) == 0

    produced = sorted(ours.rglob("*_peaks.narrowPeak"))
    assert produced, "easyepi peaks wrote no narrowPeak"
    ours_peaks = produced[0]

    # Now the same thing by hand, with the options easyepi reports it used.
    theirs = tmp_path / "theirs"
    theirs.mkdir()
    cmd = [binary, "callpeak", "-t", str(bam), "-f", "BAMPE", "-g", "hs",
           "-n", "sample", "--outdir", str(theirs), "-q", "0.01",
           "--keep-dup", "auto"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert proc.returncode == 0, proc.stdout.decode()[-2000:]
    theirs_peaks = theirs / "sample_peaks.narrowPeak"
    assert theirs_peaks.exists(), "direct macs3 wrote no narrowPeak"

    def coordinates(path):
        rows = []
        for line in path.read_text().splitlines():
            if not line or line.startswith(("#", "track")):
                continue
            p = line.split("\t")
            rows.append((p[0], int(p[1]), int(p[2])))
        return rows

    assert coordinates(ours_peaks) == coordinates(theirs_peaks), (
        "easyepi peaks and a direct macs3 call disagree; the wrapper is not "
        "passing the same options through")
