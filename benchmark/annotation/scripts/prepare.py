"""Prepare peak sets for the annotation agreement benchmark.

Why this task, and why these peaks
----------------------------------
The evaluation names the ChIPseeker agreement (43,465 H3K4me3 peaks, ~99.9954%) as the
strongest existing evidence, then warns it was reached "after multiple rounds of
correction on the same dataset" and that the version should be "frozen and tested on
independent data".

These peaks qualify as independent: they are MACS3 calls on ENCODE H3K4me3 alignments
downloaded for the previous task, and never took part in the ChIPseeker alignment work.

One peak file per cell line is used -- the deepest library of each -- so the peak set is
not the shallow, threshold-limited one.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_annotation")
PEAKS = Path("/data1/GaoLab/houoweibo/easyepi_work/pilot_h3k4me3/peaks")

PICKS = {"HepG2": "ENCFF223UKR", "K562": "ENCFF855ZMQ"}


def find_peaks(acc: str) -> Path:
    for cand in PEAKS.rglob(f"*{acc}*.narrowPeak"):
        return cand
    raise SystemExit(f"no narrowPeak for {acc}")


def main() -> int:
    ROOT.mkdir(parents=True, exist_ok=True)
    manifest = {}

    for label, acc in PICKS.items():
        path = find_peaks(acc)
        bed = ROOT / f"{label}.bed"
        n = 0
        with path.open(encoding="utf-8") as fh, bed.open("w", encoding="utf-8") as out:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.split("\t")
                if len(p) >= 3:
                    out.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")
                    n += 1
        manifest[label] = {
            "accession": acc, "source": str(path), "bed": str(bed), "peaks": n,
        }
        print(f"{label}: {n:,} peaks from {acc} -> {bed}")

    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nwritten -> {ROOT / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
