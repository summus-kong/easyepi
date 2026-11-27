"""Optional external motif-analysis tool backends.

Easyepi does not reimplement de novo motif discovery. An earlier attempt
(``homer_like``) was removed rather than deprecated: on 300-600 real peak
sequences it reported a mean 3.90 bp of shared consensus with HOMER's top ten,
which is what unrelated 8-10 bp motifs look like, and its EM variant claimed
every sequence as a site for every motif it returned. Shipping that meant
producing plausible-looking tables that could not be reproduced by the tool they
claimed to stand in for.

So discovery is delegated, and this module holds the wrappers:

  * **HOMER**   — ``findMotifsGenome.pl`` / ``findMotifs.pl`` (de novo + known)
  * **MEME Suite** — ``meme``/``dreme``/``streme`` (de novo), ``fimo`` (scan),
    ``tomtom`` (compare)
  * **TOBIAS**  — TF footprint (``TOBIAS ScoreBigwig``)

None of these is a hard dependency of the package. Installing and importing
easyepi, calling ``easyepi peaks``/``qc``/``dmr``/``annotate``, scanning a known
PWM and testing known-motif enrichment all work with no external program
present. What needs one says so: a missing executable raises
``ToolNotAvailableError`` naming the install, instead of substituting a weaker
implementation under the same filename.

Usage:
    from easyepi.motifs.external_tools import (tool_available, run_meme_discovery,
                                               run_fimo_scan, run_homer_motifs,
                                               run_tobias_footprint, discover_backends)
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np


class ToolNotAvailableError(RuntimeError):
    """Raised when an optional external tool is not installed."""


# ---------------------------------------------------------------------------
# Tool detection
# ---------------------------------------------------------------------------
def tool_available(name: str, extra_paths: Optional[Sequence[str]] = None) -> bool:
    """Return True if an executable named ``name`` can be resolved.

    Delegates to :func:`which_tool` so the two cannot disagree about what is
    installed; they did, because this function had its own PATH-only copy of the
    search and kept returning False after the shared search was widened.
    """
    return which_tool(name, extra_paths) is not None


def _prefix_bin_dirs() -> List[str]:
    """``bin`` directories belonging to the running interpreter, PATH or not.

    A conda-installed tool (HOMER, MEME, TOBIAS) lands beside the interpreter
    that will import this module. Its directory reaches PATH only when the
    environment is activated, so `python -m easyepi.cli`, a cron job or a
    scheduler submission would otherwise report an installed tool as missing. On
    the analysis host that made `discover_backends()` return ``homer: False``
    while ``<env>/bin/findMotifs.pl`` existed and ran.

    ``sys.prefix`` is the environment root, and ``sys.base_prefix`` is the base
    when a venv is layered over one, so both are checked.
    """
    dirs: List[str] = []
    for prefix in (sys.prefix, getattr(sys, "base_prefix", sys.prefix)):
        if not prefix:
            continue
        cand = os.path.join(prefix, "bin")
        if cand not in dirs and os.path.isdir(cand):
            dirs.append(cand)
    return dirs


def which_tool(name: str, extra_paths: Optional[Sequence[str]] = None) -> Optional[str]:
    """Return the absolute path to an executable, or None if not found.

    Search order: PATH, then the interpreter's own prefix, then ``extra_paths``.
    PATH stays first so an activated environment still wins.
    """
    hit = shutil.which(name)
    if hit:
        return hit
    for directory in _prefix_bin_dirs() + [str(e) for e in (extra_paths or [])]:
        cand = os.path.join(directory, name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def discover_backends(extra_paths: Optional[Sequence[str]] = None) -> Dict[str, bool]:
    """Return a dict of which motif tool backends are currently available."""
    return {
        "homer": tool_available("findMotifsGenome.pl", extra_paths)
                 or tool_available("findMotifs.pl", extra_paths),
        "meme": tool_available("meme", extra_paths),
        "dreme": tool_available("dreme", extra_paths),
        "streme": tool_available("streme", extra_paths),
        "fimo": tool_available("fimo", extra_paths),
        "tomtom": tool_available("tomtom", extra_paths),
        "tobias": tool_available("TOBIAS", extra_paths),
    }


def _require(name: str, extra_paths: Optional[Sequence[str]] = None,
             install_hint: str = "") -> str:
    path = which_tool(name, extra_paths)
    if not path:
        raise ToolNotAvailableError(
            f"External tool '{name}' is not installed. {install_hint}")
    return path


# ---------------------------------------------------------------------------
# Subprocess runner (respects conda paths / PERL5LIB for HOMER)
# ---------------------------------------------------------------------------
def _run(cmd: Sequence[str], env: Optional[Dict[str, str]] = None,
         cwd: Optional[str] = None) -> subprocess.CompletedProcess:
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    return subprocess.run(list(cmd), env=full_env, cwd=cwd,
                          capture_output=True, text=True, check=True)


# ---------------------------------------------------------------------------
# HOMER backend
# ---------------------------------------------------------------------------
def run_homer_motifs(target_fasta: str, bg_fasta: Optional[str], out_dir: str,
                     k: int = 8, homer_dir: Optional[str] = None,
                     tmp_dir: Optional[str] = None) -> dict:
    """Run HOMER's ``findMotifs.pl`` (FASTA mode) through Easyepi.

    Requires the HOMER install; if ``PERL5LIB``/PATH are not set correctly for
    a migrated install, pass ``homer_dir`` (the ``share/homer`` directory) so we
    set ``PERL5LIB`` and prepend ``bin`` to PATH. Returns a summary dict; HOMER
    writes ``homerResults`` / ``knownResults`` under ``out_dir``.
    """
    if homer_dir:
        hp = str(Path(homer_dir))
        env = {"PERL5LIB": os.path.join(hp, "bin"),
               "PATH": os.path.join(hp, "bin") + os.pathsep + os.environ.get("PATH", "")}
        findmotifs = which_tool("findMotifs.pl", [os.path.join(hp, "bin")]) \
                     or which_tool("findMotifs.pl", [os.path.join(hp, "..", "bin")])
    else:
        findmotifs = which_tool("findMotifs.pl")
        env = None

    if not findmotifs:
        raise ToolNotAvailableError(
            "HOMER findMotifs.pl not found. Install HOMER, or pass homer_dir=...")

    # findMotifs.pl calls its helpers (homer2, fasta2tab.pl, findKnownMotifs.pl,
    # compareMotifs.pl, homerTools) as bare commands through the shell. Those sit
    # in the same directory as findMotifs.pl, which is on PATH only when the conda
    # environment is activated -- `python -m easyepi.cli` from an unactivated
    # shell, a cron job or a scheduler therefore left every helper "command not
    # found". HOMER ignores that and still exits 0.
    homer_bin = str(Path(findmotifs).resolve().parent)
    base_path = (env or {}).get("PATH", os.environ.get("PATH", ""))
    env = dict(env or {})
    env.setdefault("PERL5LIB", str(Path(homer_bin).parent / "bin"))
    env["PATH"] = homer_bin + os.pathsep + base_path if homer_bin not in base_path.split(os.pathsep) else base_path

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not bg_fasta:
        raise ValueError("HOMER FASTA mode requires explicit background FASTA")
    cmd = [findmotifs, target_fasta, "fasta", str(out), "-fasta", bg_fasta,
           "-len", str(k), "-nofacts"]
    if tmp_dir:
        env = dict(env or {}) | {"TMPDIR": str(tmp_dir)}
    res = _run(cmd, env=env)

    # HOMER's own exit status does not cover its helpers: with every one of them
    # missing it printed "sh: homer2: command not found" for each step, then
    # "Job finished", and exited 0. Treating that as success wrote `completed`
    # beside a 0-byte homerMotifs.all.motifs, which reads as "no motifs found".
    combined = (res.stdout or "") + (res.stderr or "")
    missing_helpers = sorted({
        line.split(":")[0].strip().removeprefix("sh: ")
        for line in combined.splitlines()
        if "command not found" in line and line.strip().startswith("sh:")
    })
    motifs_file = out / "homerMotifs.all.motifs"
    n_motifs = 0
    if motifs_file.exists() and motifs_file.stat().st_size:
        n_motifs = sum(1 for line in motifs_file.read_text(errors="replace").splitlines()
                       if line.startswith(">"))

    if missing_helpers:
        raise RuntimeError(
            "HOMER ran but its helper programs were not found on PATH: "
            + ", ".join(missing_helpers)
            + ". The exit status is 0 in this case, so this would otherwise be "
              "reported as a successful run producing no motifs. Add HOMER's bin "
              "directory to PATH (or pass homer_dir=...) and re-run.")

    # HOMER also exits 0 when it rejects the input outright. With 3 sequences and
    # a 4 bp motif it printed "!!! Something is wrong... are you sure you chose
    # the right length" and "!!! Filtered out all motifs!!!", then "Job finished",
    # exit 0, and a 0-byte motif file. Reported as success, that reads as "no
    # motifs in this peak set", which is a different and wrong conclusion.
    rejected = [line.strip() for line in combined.splitlines()
                if line.strip().startswith("!!!")]
    if rejected:
        raise RuntimeError(
            "HOMER rejected the input and still exited 0:\n  "
            + "\n  ".join(rejected[:6])
            + f"\nHOMER needs enough sequences for a -len {k} motif (a few "
              "sequences are not enough to estimate a background). Supply more "
              "peak sequences; the run was not a valid negative result.")

    return {
        "backend": "homer", "available": True,
        "command": " ".join(cmd), "returncode": res.returncode,
        "stderr": res.stderr[-2000:] if res.returncode else "",
        "out_dir": str(out),
        "n_denovo_motifs": n_motifs,
    }


# ---------------------------------------------------------------------------
# MEME Suite backends
# ---------------------------------------------------------------------------
def run_meme_discovery(target_fasta: str, out_dir: str,
                       nmotifs: int = 5, minw: int = 8, maxw: int = 20,
                       objfun: str = "classic") -> dict:
    """Run ``meme`` de novo motif discovery (MEME Suite is an optional dep)."""
    executable = _require("meme", install_hint="conda install -c bioconda meme")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cmd = [executable, target_fasta, "-oc", str(out), "-nmotifs", str(nmotifs),
           "-minw", str(minw), "-maxw", str(maxw), "-objfun", objfun, "-dna",
           "-mod", "zoops"]
    res = _run(cmd)
    # MEME writes meme.xml / meme.txt in out_dir
    motif_file = out / "meme.txt"
    return {"backend": "meme", "available": True, "command": " ".join(cmd),
            "returncode": res.returncode, "stderr": res.stderr[-2000:] if res.returncode else "",
            "out_dir": str(out), "motif_file": str(motif_file) if motif_file.exists() else None}


def run_fimo_scan(motif_db: str, target_fasta: str, out_dir: str,
                  thresh: float = 1e-4) -> dict:
    """Run ``fimo`` motif scanning (MEME Suite optional dep)."""
    executable = _require("fimo", install_hint="conda install -c bioconda meme")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cmd = [executable, "--thresh", str(thresh), "--text", motif_db, target_fasta]
    res = _run(cmd)
    hit_file = out / "fimo.tsv"
    hit_file.write_text(res.stdout, encoding="utf-8")
    hits = 0
    for line in res.stdout.splitlines():
        if line.startswith(("#", "motif_id\t", "pattern name\t")) or not line.strip():
            continue
        hits += 1
    return {"backend": "fimo", "available": True, "command": " ".join(cmd),
            "returncode": res.returncode, "n_hits": hits, "hits_file": str(hit_file),
            "stderr": res.stderr[-2000:] if res.returncode else ""}


def run_tomtom_compare(motif_db: str, query_motif_file: str, out_dir: str) -> dict:
    """Compare motifs through Tomtom when the tool is installed."""
    executable = _require("tomtom", install_hint="conda install -c bioconda meme")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cmd = [executable, query_motif_file, motif_db, "-oc", str(out)]
    res = _run(cmd)
    return {"backend": "tomtom", "available": True, "command": " ".join(cmd),
            "returncode": res.returncode, "out_dir": str(out),
            "stderr": res.stderr[-2000:] if res.returncode else ""}


# ---------------------------------------------------------------------------
# TOBIAS backend (TF footprint)
# ---------------------------------------------------------------------------
def run_tobias_footprint(bigwig_paths: Sequence[str], bed_path: str, genome: str,
                         out_dir: str, motif_db: Optional[str] = None) -> dict:
    """Run ``TOBIAS ScoreBigwig`` on an ATAC/CUT&Tag cut-site bigWig.

    TOBIAS renamed this step: 0.17.5 accepts ``ScoreBigwig`` (and the
    legacy ``FootprintScores``); ``Footprint`` was never a subcommand.
    """
    executable = _require("TOBIAS", install_hint="pip install tobias")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if len(bigwig_paths) != 1:
        raise ValueError("Run TOBIAS ScoreBigwig separately for each bias-corrected signal track")
    if motif_db is not None:
        raise ValueError("Footprint scores signal; motif binding inference requires a separate BINDetect step")
    # TOBIAS 0.17.5 has no `Footprint` subcommand; `ScoreBigwig` is its
    # documented replacement for `FootprintScores` and takes no --cores.
    cmd = [executable, "ScoreBigwig", "--signal", str(bigwig_paths[0]),
           "--regions", bed_path, "--output", str(out / "footprint.bw")]
    res = _run(cmd)
    return {"backend": "tobias", "available": True, "command": " ".join(cmd),
            "returncode": res.returncode, "out_dir": str(out),
            "stderr": res.stderr[-2000:] if res.returncode else ""}
