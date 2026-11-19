"""Motif format layer (L2): HOMER / JASPAR / MEME I/Os and cross-conversion.

Supports reading and writing three widely used Position-Frequency-Matrix (PFM)
formats so Easyepi can consume and produce artifacts compatible with the
popular motif suites:

  * HOMER  ``.motif``  (also the ~/.homer *data/knownTFs* convention)
  * JASPAR ``.pfm`` / ``.jaspar``
  * MEME   ``.meme``

Each PFM is represented as a small dataclass (``Motif``) carrying the name, the
4 x L frequency matrix, an optional motif id and source database.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np

_BASES = "ACGT"
_BASE_TO_I = {b: i for i, b in enumerate(_BASES)}


@dataclass
class Motif:
    """A single position-frequency matrix and its metadata."""

    name: str
    matrix: np.ndarray                      # shape (4, L), ACGT rows, frequencies
    id: Optional[str] = None
    source: Optional[str] = None            # e.g. 'jaspar', 'hocomoco', 'cisbp', 'homer'
    tf: Optional[str] = None                # transcription factor gene symbol
    consensus: Optional[str] = None

    @property
    def length(self) -> int:
        """Motif length in positions."""
        return int(self.matrix.shape[1])

    def to_log_odds(self, bg: Optional[np.ndarray] = None) -> np.ndarray:
        """Convert the frequency matrix to a log-odds (bits) PWM.

        A column whose counts are all zero divides by zero and would turn the whole
        PWM into NaN; every window then scores as non-finite and the motif silently
        yields "no site". Such a column carries no information, so it is replaced by
        the uniform distribution before normalising.
        """
        freq = np.asarray(self.matrix, dtype=np.float64)
        if freq.ndim != 2 or freq.shape[0] != 4:
            raise ValueError("motif matrix must be (4, length)")
        if np.any(freq.sum(axis=0) <= 0):
            freq = np.where(freq.sum(axis=0, keepdims=True) <= 0, 1.0, freq)
        freq = freq / freq.sum(axis=0, keepdims=True)
        if bg is None:
            bg = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
        with np.errstate(divide="ignore"):
            return np.log2(np.maximum(freq, 1e-9) / np.maximum(bg[:, None], 1e-9))

    def consensus_seq(self) -> str:
        """Return the most probable base per position."""
        return "".join(_BASES[int(a)] for a in np.argmax(self.matrix, axis=0))

    @property
    def pfm(self) -> np.ndarray:
        """Raw (counts or frequencies) ACGT x L matrix."""
        return self.matrix


# ---------------------------------------------------------------------------
# HOMER .motif
# ---------------------------------------------------------------------------
def read_homer(path: Union[str, Path]) -> List[Motif]:
    """Parse a HOMER ``.motif`` file (one motif per 4-row block).

    HOMER motif format::
        >MotifName\tconsensus\tweight\t[TF]\t[P-value]\t[...optional...]
        A\t0.2\t0.1\t0.0\t...
        C\t0.2\t0.1\t...
        G\t...
        T\t...
    """
    motifs: List[Motif] = []
    lines = [ln.rstrip("\n") for ln in open(path)]
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        if line.startswith(">"):
            header = line[1:].split("\t")
            name = header[0]
            tf = header[3] if len(header) > 3 else None
            consensus = header[1] if len(header) > 1 else None
            # read the next 4 rows as A,C,G,T
            rows = []
            j = i + 1
            while j < len(lines) and len(rows) < 4:
                p = lines[j].split("\t")
                if not p or p[0] not in _BASE_TO_I:
                    break
                try:
                    vals = [float(x) for x in p[1:]]
                except ValueError:
                    break
                rows.append(vals)
                j += 1
            if len(rows) == 4:
                arr = np.array(rows, dtype=np.float64)
                # rows may be A,C,G,T order already; ensure ACGT
                order = [_BASE_TO_I[p[0]] for p in [lines[x].split("\t") for x in range(i + 1, i + 5)]]
                mat = np.zeros((4, arr.shape[1]), dtype=np.float64)
                for r, base_idx in zip(range(4), order):
                    mat[base_idx] = arr[r]
                motifs.append(Motif(name=name, matrix=mat, id=None,
                                    source="homer", tf=tf, consensus=consensus))
            i = j
        else:
            i += 1
    return motifs


def write_homer(motifs: List[Motif], path: Union[str, Path]) -> None:
    """Write motifs in HOMER ``.motif`` format (one 5-line block each)."""
    with open(path, "w") as fh:
        for m in motifs:
            cons = m.consensus_seq() if m.consensus is None else m.consensus
            tf = m.tf or m.name
            fh.write(f">{m.name}\t{cons}\t1.0\t{tf}\t0.0\n")
            for base_idx, base in enumerate(_BASES):
                fh.write(f"{base}\t" + "\t".join("%.5f" % v for v in m.matrix[base_idx]) + "\n")


def read_homer_motifs(path) -> List[Motif]:
    """Parse a HOMER ``.motifs`` collection file (multiple motifs concatenated).

    Differs from :func:`read_homer`: blocks are ``>consensus<TAB>info...`` then
    exactly 4 lines giving the A, C, G, T position probabilities (no base
    letters on the rows). Used for HOMER's ``known.motifs`` / ``all.motifs``.
    """
    lines = [ln.rstrip("\n") for ln in open(path)]
    motifs: List[Motif] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith(">"):
            header = line[1:].split("\t")
            consensus = header[0]
            src = header[1] if len(header) > 1 else None
            name = (src or consensus).split("/")[0]
            tf = (src.split("/")[0] if src else None)
            # Parse ALL consecutive numeric rows (A,C,G,T probabilities), one
            # row per motif position, i.e. L rows for a length-L motif.
            rows = []
            j = i + 1
            while j < len(lines):
                p = lines[j].split()
                try:
                    vals = [float(x) for x in p]
                except ValueError:
                    break
                if len(vals) < 4:
                    break
                rows.append(vals[:4])
                j += 1
            if len(rows) >= 4:
                mat = np.array(rows, dtype=np.float64)  # (L, 4) A,C,G,T
                # HOMER rows are ACGT order; transpose to (4, L) for the Motif
                mat = mat.T
                motifs.append(Motif(name=name, matrix=mat, id=None, source="homer",
                                    tf=tf, consensus=consensus))
                i = j
            else:
                i += 1
        else:
            i += 1
    return motifs


# ---------------------------------------------------------------------------
# JASPAR .pfm / .jaspar
# ---------------------------------------------------------------------------
def read_jaspar(path: Union[str, Path]) -> List[Motif]:
    """Parse a JASPAR ``.pfm`` file (4 lines per matrix: A, C, G, T counts).

    Supports the plain .pfm (4 lines of whitespace/comma-separated numbers) and
    the matrix with an optional ID (e.g. '>MA0004.1  ARNT').
    """
    motifs: List[Motif] = []
    motifs_raw = []
    cur_name = None
    cur_id = None
    cur_rows: List[List[float]] = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_rows:
                motifs_raw.append((cur_name, cur_id, np.array(cur_rows)))
            parts = line[1:].split()
            cur_name = parts[1] if len(parts) > 1 else parts[0]
            cur_id = parts[0] if parts else None
            cur_rows = []
        else:
            nums = [float(x) for x in line.replace(",", " ").split()]
            if nums:
                cur_rows.append(nums)
    if cur_rows:
        motifs_raw.append((cur_name, cur_id, np.array(cur_rows)))
    for name, mid, mat in motifs_raw:
        if mat.shape[0] != 4:
            continue
        motifs.append(Motif(name=name, matrix=mat, id=mid, source="jaspar"))
    return motifs


def write_jaspar(motifs: List[Motif], path: Union[str, Path]) -> None:
    """Write motifs in JASPAR ``.pfm`` format (A,C,G,T rows)."""
    with open(path, "w") as fh:
        for m in motifs:
            mid = m.id or m.name
            fh.write(f">{mid}\t{m.name}\n")
            for base_idx in range(4):
                fh.write(" ".join("%.4f" % v for v in m.matrix[base_idx]) + "\n")


# ---------------------------------------------------------------------------
# MEME .meme
# ---------------------------------------------------------------------------
def read_meme(path: Union[str, Path]) -> List[Motif]:
    """Parse a MEME ``.meme`` motif file (``motif <id>`` blocks with a 4xL matrix)."""
    motifs: List[Motif] = []
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("MOTIF"):
            parts = ln.split()
            mid = parts[1] if len(parts) > 1 else None
            name = " ".join(parts[2:]) or mid
            # find 'letter-probability matrix' on following lines
            j = i + 1
            while j < len(lines) and "letter-probability" not in lines[j]:
                j += 1
            if j < len(lines) and "alength" in lines[j]:
                k = int(lines[j].split("alength=")[-1].split()[0])
                mat = np.zeros((4, k), dtype=np.float64)
                r = 0
                j += 1
                while j < len(lines) and r < 4:
                    nums = [float(x) for x in lines[j].split()]
                    if nums:
                        mat[r] = nums[:k]
                        r += 1
                    j += 1
                if r == 4:
                    motifs.append(Motif(name=name, matrix=mat, id=mid, source="meme"))
                i = j
            else:
                i += 1
        else:
            i += 1
    return motifs


def write_meme(motifs: List[Motif], path: Union[str, Path]) -> None:
    """Write motifs in MEME ``.meme`` format."""
    with open(path, "w") as fh:
        fh.write("MEME version 5.5.0\n\n")
        fh.write("ALPHABET= ACGT\n\n")
        fh.write("strands: + -\n\n")
        fh.write("Background letter frequencies (from uniform):\n")
        fh.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for m in motifs:
            mid = m.id or m.name
            fh.write(f"MOTIF {mid} {m.name}\n\n")
            k = m.length
            fh.write(f"letter-probability matrix: alength= {k} w= {k}\n")
            freq = m.matrix / m.matrix.sum(axis=0, keepdims=True)
            for base_idx in range(4):
                fh.write("  " + " ".join("%.5f" % v for v in freq[base_idx]) + "\n")
            fh.write("\n")


# ---------------------------------------------------------------------------
# Generic reader + cross-conversion + dedup merge
# ---------------------------------------------------------------------------
READERS = {
    ".homer": read_homer,
    ".motif": read_homer,
    ".motifs": read_homer_motifs,
    ".pfm": read_jaspar,
    ".jaspar": read_jaspar,
    ".meme": read_meme,
}


def read_motifs(path: Union[str, Path]) -> List[Motif]:
    """Dispatch to the appropriate reader based on file extension."""
    ext = Path(path).suffix.lower()
    if ext not in READERS:
        raise ValueError(f"unknown motif format '{ext}' (supported: {list(READERS)})")
    return READERS[ext](path)


def write_motifs(motifs: List[Motif], path: Union[str, Path]) -> None:
    """Write motifs in the format implied by the output extension."""
    ext = Path(path).suffix.lower()
    if ext in (".homer", ".motif"):
        write_homer(motifs, path)
    elif ext in (".pfm", ".jaspar"):
        write_jaspar(motifs, path)
    elif ext == ".meme":
        write_meme(motifs, path)
    else:
        raise ValueError(f"unknown motif format '{ext}'")


def dedup_merge(motifs: List[Motif], threshold: float = 0.9) -> List[Motif]:
    """Remove near-duplicate motifs (by reversed-complement-aware similarity).

    Uses a simple, fast heuristic: same length and a high (>= ``threshold``)
    Pearson correlation between one motif's column frequencies and the other
    motif (or its reverse complement). Keeps the first occurrence of each
    cluster.
    """
    out: List[Motif] = []
    for m in motifs:
        dup = False
        for seen in out:
            if _similar(m, seen, threshold):
                dup = True
                break
        if not dup:
            out.append(m)
    return out


def _rc(mat: np.ndarray) -> np.ndarray:
    """Reverse complement of a 4 x L frequency matrix."""
    # swap A<->T (0<->3) and C<->G (1<->2), then reverse columns
    swap = np.array([3, 2, 1, 0], dtype=int)
    return mat[swap, ::-1]


def _similar(a: Motif, b: Motif, threshold: float) -> bool:
    if a.length != b.length:
        return False
    fa = a.matrix / a.matrix.sum(axis=0, keepdims=True)
    fb = b.matrix / b.matrix.sum(axis=0, keepdims=True)
    cands = [fb, _rc(fb)]
    for cb in cands:
        corr = float(np.corrcoef(fa.ravel(), cb.ravel())[0, 1])
        if np.isnan(corr):
            continue
        if corr >= threshold:
            return True
    return False


def merge_databases(sources: List[List[Motif]], threshold: float = 0.9) -> List[Motif]:
    """Merge several motif collections and deduplicate across sources."""
    flat = [m for group in sources for m in group]
    return dedup_merge(flat, threshold=threshold)
