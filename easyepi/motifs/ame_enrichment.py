"""Known-motif enrichment by AME (MEME Suite), with its parameters exposed.

Why AME rather than a local implementation
------------------------------------------
Three implementations were measured against HOMER on identical input, and against
each other:

                       time    rho(pct_target)   rho(log_p)   top-20 vs HOMER
    local scan+binomial  77s          0.109          0.237          5/20
    port of HOMER's test 88s          0.120          0.206          4/20
    FIMO + hypergeometric 89s         0.169          0.233          4/20

and FIMO versus the local scan reached rho 0.292 / 0.435 with 6/20 top-20 overlap.
Two independent standard scanners disagree with each other about as much as either
disagrees with HOMER. The ranking is therefore dominated by scanner and partition
choices, not by which implementation is "correct", and chasing numerical agreement
with one tool was the wrong objective.

AME is the published tool for exactly this question (Analysis of Motif Enrichment,
MEME Suite). It is already installed, its test and scoring are explicit flags rather
than an internal search, and it exposes 6 tests x 4 scorings. Delegating to it means
the method is citable and the reader can see which variant produced a result, instead
of inheriting an undocumented internal choice.

What is deliberately NOT claimed: results here are not numerically identical to
HOMER, FIMO or a hand-rolled hypergeometric. AME partitions sequences by affinity
score; HOMER searches a threshold. They answer overlapping but not identical
questions.

The output records the executable, its version, and every flag used, so a run can be
reproduced from the JSON alone.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence

# AME's own defaults, echoed here so the signature shows what will run.
_LONG_DEFAULTS = {
    "scoring": "avg",
    "hit_lo_fraction": None,      # AME default 0.25; only passed when asked
    "kmer": None,                 # AME default 2; only used with --shuffle--
    "seed": None,
    "evalue_report_threshold": None,
}


def _write_fasta(path: Path, seqs: Sequence[str], prefix: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            if s:
                fh.write(f">{prefix}{i}\n{s.upper()}\n")


def _ame_version(executable: str) -> str:
    try:
        r = subprocess.run([executable, "--version"], capture_output=True,
                           text=True, timeout=60)
        text = (r.stdout or r.stderr or "").strip().splitlines()
        return text[0] if text else "unknown"
    except Exception:
        return "unknown"


def parse_ame_tsv(text: str) -> List[Dict]:
    """AME --text output to rows. Keeps the header so unknown columns survive."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows: List[Dict] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            continue
        rec = dict(zip(header, fields))
        # numeric where numeric, so callers do not have to re-parse
        for key in ("rank", "tests", "pos", "neg", "TP", "FP",
                    "FASTA_max", "PWM_min"):
            if key in rec:
                try:
                    rec[key] = int(rec[key])
                except ValueError:
                    pass
        for key in ("p-value", "adj_p-value", "E-value", "%TP", "%FP"):
            if key in rec:
                try:
                    rec[key] = float(rec[key])
                except ValueError:
                    pass
        rows.append(rec)
    return rows


def ame_enrichment(target_seqs: Sequence[str],
                   bg_seqs: Sequence[str] | None,
                   motif_db: str,
                   method: str = "fisher",
                   scoring: str = "avg",
                   shuffle_control: bool = False,
                   kmer: int = 2,
                   seed: int = 1,
                   evalue_threshold: float = 10.0,
                   hit_lo_fraction: float | None = None,
                   ame_executable: str | None = None) -> Dict:
    """Run AME and return its table plus the exact invocation.

    ``bg_seqs`` is the control set. Passing ``None`` with ``shuffle_control=True``
    uses AME's own shuffled control, which preserves k-mer content to ``kmer``
    bases; that is the cheaper option when no matched background exists.

    ``method`` is one of fisher, 3dmhg, 4dmhg, ranksum, pearson, spearman.
    ``scoring`` is one of avg, max, sum, totalhits. Both are AME's own vocabulary
    and are validated by AME itself, so no local list is duplicated here.
    """
    # Argument checks come first: a caller who passes no null should be told that,
    # not told to install AME. Otherwise the two failures are indistinguishable.
    target = [s for s in target_seqs if s]
    if not target:
        raise ValueError("need at least one target sequence")
    if bg_seqs is None and not shuffle_control:
        raise ValueError(
            "Give a control set (bg_seqs) or set shuffle_control=True. A motif "
            "enrichment test without a null has no meaning.")

    exe = ame_executable or shutil.which("ame")
    if exe is None or (ame_executable and not Path(ame_executable).exists()):
        raise RuntimeError(
            f"AME not found{f' at {ame_executable}' if ame_executable else ''}. "
            "Install MEME Suite (conda install -c bioconda meme) or pass "
            "ame_executable.")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        primary = tmp / "targets.fa"
        _write_fasta(primary, target, "t")

        cmd = [exe, "--text", "--method", method, "--scoring", scoring,
               "--seed", str(seed),
               "--evalue-report-threshold", str(evalue_threshold)]
        if hit_lo_fraction is not None:
            cmd += ["--hit-lo-fraction", str(hit_lo_fraction)]

        if bg_seqs is not None:
            control = tmp / "control.fa"
            _write_fasta(control, [s for s in bg_seqs if s], "c")
            cmd += ["--control", str(control)]
            control_kind = "supplied"
        else:
            cmd += ["--control", "--shuffle--", "--kmer", str(kmer)]
            control_kind = f"shuffle(kmer={kmer})"

        cmd += [str(primary), str(motif_db)]

        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if proc.returncode != 0:
            raise RuntimeError(
                f"ame failed (exit {proc.returncode}): {proc.stderr[-900:]}")

        rows = parse_ame_tsv(proc.stdout)

    # Paths in the recorded command point at a temp dir that no longer exists;
    # record the shape of the call instead so it stays readable and reproducible.
    reproducible = [str(Path(c).name) if i >= 2 and c.startswith("/tmp") else c
                    for i, c in enumerate(cmd)]

    return {
        "backend": "ame",
        "ame_version": _ame_version(exe),
        "ame_executable": exe,
        "method": method,
        "scoring": scoring,
        "control": control_kind,
        "n_target": len(target),
        "n_control": len([s for s in (bg_seqs or []) if s]),
        "motif_db": str(motif_db),
        "evalue_threshold": evalue_threshold,
        "command_shape": " ".join(reproducible),
        "n_motifs_reported": len(rows),
        "motifs": rows,
    }
