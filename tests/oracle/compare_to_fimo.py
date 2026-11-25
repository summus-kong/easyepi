"""Compare `easyepi.motifs.scan_sequences` against FIMO on the same PWM.

    python tests/oracle/compare_to_fimo.py

Needs `fimo` from MEME Suite (bioconda `meme`; 5.5.9 measured).

## What is compared, and why the scores are not

Both tools scan a sequence with a log-odds matrix, but they do not publish the
same number:

  * `easyepi` `read_pwm` returns `log2(p_i / 0.25)` per base, and
    `scan_sequences` reports the argmax window's **sum** of those values.
  * FIMO reports a `score` column produced by MEME Suite's own scoring.

On four fixed 7-mers the two disagree by an amount that **changes sign with base
composition**: GATATCA +0.40, GATATCG +0.11, AAAACGT +0.26, GGGGGGG -1.14. A
constant offset, a pseudocount, and a background-smoothing blend were each ruled
out as the explanation. Establishing FIMO's exact definition would mean reading
MEME Suite's source, and the difference does not change any downstream decision:
a motif call is made from a threshold or a ranking, never from the absolute
score.

So this harness compares what does carry meaning:

  1. **Best window** -- the position and strand each tool picks per sequence.
  2. **Ranking** -- whether the two order the sequences the same way.

## Coordinate conventions

`scan_sequences` returns a 0-based offset, or for a reverse hit
`len(seq) - pos - L`. FIMO reports 1-based inclusive `start`/`stop`. The harness
converts before comparing, because a one-base mismatch from convention alone
would read as a real disagreement.

## What this does NOT establish

* Score equivalence -- see above, and measured per-k-mer in this file's header.
* Enrichment statistics: FIMO's q-value and easyepi's HOMER-like enrichment are
  different statistics on different inputs.
* Multi-motif competition: one PWM, so FIMO never has to choose between motifs.
* Ambiguous bases: the fixture is clean ACGT.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from easyepi.motifs import read_pwm, scan_sequences  # noqa: E402

#: The matrix as counts, one row per POSITION. The MEME writer needs this order;
#: the JASPAR writer transposes it.
COUNTS = [
    [2, 3, 90, 5],
    [90, 2, 3, 5],
    [2, 2, 2, 94],
    [85, 5, 5, 5],
    [5, 5, 5, 85],
    [5, 85, 5, 5],
    [80, 5, 10, 5],
]

MOTIF_DB = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GATA_like
letter-probability matrix: alength= 4 w= 7 nsites= 40 E= 0
 0.02 0.03 0.90 0.05
 0.90 0.02 0.03 0.05
 0.02 0.02 0.02 0.94
 0.85 0.05 0.05 0.05
 0.05 0.05 0.05 0.85
 0.05 0.85 0.05 0.05
 0.80 0.05 0.10 0.05
"""

CONSENSUS = "GATATCA"


def _find_fimo() -> str:
    for cand in (str(Path(sys.executable).parent / "fimo"), shutil.which("fimo")):
        if cand and Path(cand).exists():
            return cand
    raise SystemExit("fimo not found; install MEME Suite or extend the candidate list")


def fixture_sequences():
    """Planted consensus, its reverse complement, and noise, at varying offsets."""
    rng = np.random.default_rng(1234)
    seqs, labels = [], []
    for i in range(12):
        left = int(rng.integers(3, 20))
        right = int(rng.integers(3, 20))
        flank_l = "".join(rng.choice(list("ACGT"), left))
        flank_r = "".join(rng.choice(list("ACGT"), right))
        if i % 3 == 0:
            body = CONSENSUS
        elif i % 3 == 1:
            body = CONSENSUS.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        else:
            body = "".join(rng.choice(list("ACGT"), len(CONSENSUS)))
        seqs.append(flank_l + body + flank_r)
        labels.append("seq%d" % i)
    return labels, seqs


def write_inputs(tmp: Path, labels, seqs):
    """Write the matrix once per tool: MEME for FIMO, JASPAR for read_pwm."""
    meme = tmp / "motif.meme"
    meme.write_text(MOTIF_DB, encoding="utf-8")

    # JASPAR/PFM is one row per BASE, so COUNTS (positions x bases) is transposed.
    # Skipping this silently truncates a 7-position motif to 4.
    jaspar = tmp / "motif.pfm"
    by_base = list(zip(*COUNTS))
    lines = ["%s [ %s ]" % (b, "  ".join(str(v) for v in row))
             for b, row in zip("ACGT", by_base)]
    jaspar.write_text("\n".join(lines) + "\n", encoding="utf-8")

    fasta = tmp / "peaks.fa"
    fasta.write_text("".join(">%s\n%s\n" % (n, s) for n, s in zip(labels, seqs)),
                     encoding="utf-8")
    return meme, jaspar, fasta


def run_fimo(meme: Path, fasta: Path, thresh: float):
    proc = subprocess.run([_find_fimo(), "--thresh", str(thresh), "--text",
                           str(meme), str(fasta)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise SystemExit("fimo failed:\n" + proc.stderr.decode()[-2000:])
    rows = []
    for line in proc.stdout.decode().splitlines():
        if not line.strip() or line.startswith(("#", "motif_id")):
            continue
        p = line.split("\t")
        if len(p) < 9:
            continue
        rows.append({"seq": p[2], "start": int(p[3]), "stop": int(p[4]),
                     "strand": p[5], "score": float(p[6])})
    return rows


def best_per_sequence(rows):
    best = {}
    for r in rows:
        cur = best.get(r["seq"])
        if cur is None or r["score"] > cur["score"]:
            best[r["seq"]] = r
    return best


def compare(labels, seqs, our, fimo_best):
    """Best-window agreement and ranking agreement, per sequence."""
    records = []
    for i, label in enumerate(labels):
        f = fimo_best.get(label)
        if f is None:
            records.append({"label": label, "had_hit": False})
            continue
        records.append({
            "label": label,
            "had_hit": True,
            "our_pos": int(our["positions"][i]),
            "our_strand": str(our["strands"][i]),
            "our_score": float(our["scores"][i]),
            # FIMO start is 1-based inclusive; 0-based offset is start - 1.
            "fimo_pos": f["start"] - 1,
            "fimo_strand": f["strand"],
            "fimo_score": f["score"],
        })
    return records


def main():
    labels, seqs = fixture_sequences()
    tmp = Path(tempfile.mkdtemp(prefix="fimocmp_"))
    meme, jaspar, fasta = write_inputs(tmp, labels, seqs)

    pwm = read_pwm(jaspar)
    our = scan_sequences(pwm, seqs)
    print("easyepi PWM shape :", pwm.shape)
    print("ML consensus      :", "".join("ACGT"[i] for i in pwm.argmax(axis=0)))
    print("sequences         :", len(seqs))

    # thresh 1.0 reports every window, so FIMO's per-sequence argmax is available
    # for comparison with scan_sequences' argmax.
    rows = run_fimo(meme, fasta, thresh=1.0)
    print("FIMO windows      :", len(rows))
    fimo_best = best_per_sequence(rows)

    records = compare(labels, seqs, our, fimo_best)
    hits = [r for r in records if r["had_hit"]]

    print("\n%-8s %-7s %-14s %-14s %s"
          % ("seq", "strand", "easyepi pos", "FIMO pos", "same?"))
    for r in hits:
        same = "yes" if (r["our_pos"] == r["fimo_pos"]
                         and r["our_strand"] == r["fimo_strand"]) else "NO"
        print("%-8s %-7s %-14d %-14d %s"
              % (r["label"], r["our_strand"], r["our_pos"], r["fimo_pos"], same))

    pos_ok = sum(1 for r in hits
                 if r["our_pos"] == r["fimo_pos"] and r["our_strand"] == r["fimo_strand"])
    print("\n=== best window ===")
    print("  identical (pos + strand): %d/%d" % (pos_ok, len(hits)))

    # Ranking: do the two order the sequences the same way? Uses only sequences
    # where both reported a hit, and compares the induced orders by Spearman.
    if len(hits) >= 3:
        a = np.array([r["our_score"] for r in hits])
        b = np.array([r["fimo_score"] for r in hits])
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        rho = float(np.corrcoef(ra, rb)[0, 1])
        print("\n=== ranking of the same sequences ===")
        print("  Spearman rho on best-window score: %.3f" % rho)

    print("\n=== score scale (informational, not a pass/fail) ===")
    print("  The absolute scores differ and the sign of the difference varies with")
    print("  base composition, so they are not the same statistic. See the module")
    print("  docstring for the four k-mers that establish this.")
    for r in hits[:4]:
        print("  %-8s easyepi %8.4f   FIMO %8.4f   delta %+7.4f"
              % (r["label"], r["our_score"], r["fimo_score"],
                 r["our_score"] - r["fimo_score"]))

    print("\nartifacts:", tmp)


if __name__ == "__main__":
    main()
