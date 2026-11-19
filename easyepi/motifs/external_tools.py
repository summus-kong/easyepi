"""Optional external motif-analysis tool backends.

Easyepi's native motif modules (``homer_like`` / ``homer_workflow``) are the
zero-dependency path. This module adds *optional* wrappers around the popular
external motif suites, so a user can run them through Easyepi's unified
interface when they are installed:

  * **HOMER**   — ``findMotifsGenome.pl`` / ``findMotifs.pl`` (de novo + known)
  * **MEME Suite** — ``meme``/``dreme``/``streme`` (de novo), ``fimo`` (scan),
    ``tomtom`` (compare)
  * **TOBIAS**  — TF footprint (``TOBIAS Footprint``)

Each backend is only usable if the corresponding executable is present; a
missing tool raises a clear ``ToolNotAvailableError`` telling the user how to
install it, without breaking the core (dependency-free) path.

Usage:
    from easyepi.motifs.external_tools import (tool_available, run_meme_discovery,
                                               run_fimo_scan, run_homer_motifs,
                                               run_tobias_footprint, discover_backends)
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np


class ToolNotAvailableError(RuntimeError):
    """Raised when an optional external tool is not installed."""


# ---------------------------------------------------------------------------
# Tool detection
# ---------------------------------------------------------------------------
def tool_available(name: str, extra_paths: Optional[Sequence[str]] = None) -> bool:
    """Return True if an executable named ``name`` is found on PATH (or in
    ``extra_paths``)."""
    if shutil.which(name):
        return True
    for p in (extra_paths or []):
        cand = os.path.join(p, name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return True
    return False


def which_tool(name: str, extra_paths: Optional[Sequence[str]] = None) -> Optional[str]:
    """Return the absolute path to an executable, or None if not found."""
    hit = shutil.which(name)
    if hit:
        return hit
    for p in (extra_paths or []):
        cand = os.path.join(p, name)
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
        env = None
        findmotifs = which_tool("findMotifs.pl")

    if not findmotifs:
        raise ToolNotAvailableError(
            "HOMER findMotifs.pl not found. Install HOMER, or pass homer_dir=...")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not bg_fasta:
        raise ValueError("HOMER FASTA mode requires explicit background FASTA")
    cmd = [findmotifs, target_fasta, "fasta", str(out), "-fasta", bg_fasta,
           "-len", str(k), "-nofacts"]
    if tmp_dir:
        env = dict(env or {}) | {"TMPDIR": str(tmp_dir)}
    res = _run(cmd, env=env)
    return {
        "backend": "homer", "available": True,
        "command": " ".join(cmd), "returncode": res.returncode,
        "stderr": res.stderr[-2000:] if res.returncode else "",
        "out_dir": str(out),
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
    """Run ``TOBIAS Footprint`` on ATAC/CUT&Tag bigwigs (optional dep)."""
    executable = _require("TOBIAS", install_hint="pip install tobias")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if len(bigwig_paths) != 1:
        raise ValueError("Run TOBIAS Footprint separately for each bias-corrected signal track")
    if motif_db is not None:
        raise ValueError("Footprint scores signal; motif binding inference requires a separate BINDetect step")
    cmd = [executable, "Footprint", "--cores", "4", "--signal", str(bigwig_paths[0]),
           "--regions", bed_path, "--output", str(out / "footprint.bw")]
    res = _run(cmd)
    return {"backend": "tobias", "available": True, "command": " ".join(cmd),
            "returncode": res.returncode, "out_dir": str(out),
            "stderr": res.stderr[-2000:] if res.returncode else ""}
