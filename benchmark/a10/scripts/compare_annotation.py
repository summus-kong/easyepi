"""A10 core: per-peak annotation agreement between EasyEpi and ChIPseeker.

Both annotators receive byte-identical peak sets (called once by EasyEpi in
prepare.py) and annotations derived from the same TxDb, so a difference can
only come from the annotation rules themselves.

Reported per dataset, with no thresholded "pass":
  * exact label match           - same bin label on the same peak
  * confusion by label          - where the disagreements live
  * per-label recall            - how often each bin is recovered

H3K27me3 is called but excluded from the agreement table: it is a broad mark,
a narrowPeak set of 15 intervals cannot represent it, and comparing two
annotators on 15 peaks would be noise presented as a result. The exclusion is
explicit in the output rather than silent.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")
PY = "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/bin/python"
RSCRIPT = "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/seurat5/bin/Rscript"

# Broad marks whose narrowPeak representation is not meaningful.
EXCLUDED = {"H3K27me3"}


def sh(cmd: list[str], timeout: int = 1800) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def read_tsv(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    try:
        idx = header.index("annotation")
    except ValueError:
        idx = len(header) - 1
    return [ln.split("\t")[idx] for ln in lines[1:] if ln.strip()]


def main() -> int:
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    gtf = manifest["gtf"]
    results = {}

    for name, info in manifest["peaks"].items():
        if name in EXCLUDED:
            results[name] = {"excluded": "broad mark; narrowPeak not meaningful"}
            print(f"{name}: excluded (broad mark)")
            continue
        peaks = Path(info["peaks"])
        outdir = ROOT / "annot" / name
        outdir.mkdir(parents=True, exist_ok=True)

        # EasyEpi
        ee_tsv = outdir / "easyepi.tsv"
        r = sh([PY, "-m", "easyepi.cli", "annotate", "--gtf", gtf,
                "--peaks", str(peaks), "--out", str(outdir / "ee")])
        cands = sorted((outdir / "ee").rglob("annotated.tsv"))
        if not cands:
            results[name] = {"error": f"easyepi annotate failed: "
                                      f"{(r.stderr or r.stdout)[-300:]}"}
            print(f"{name}: EasyEpi FAILED")
            continue
        ee_tsv = cands[0]

        # ChIPseeker
        cs_tsv = outdir / "chipseeker.tsv"
        r2 = sh([RSCRIPT, "--vanilla",
                 str(ROOT / "scripts" / "annotate_chipseeker.R"),
                 str(peaks), str(GTF := gtf), str(cs_tsv)])
        if not cs_tsv.exists():
            results[name] = {"error": f"ChIPseeker failed: "
                                      f"{(r2.stderr or r2.stdout)[-300:]}"}
            print(f"{name}: ChIPseeker FAILED")
            continue

        a = read_tsv(ee_tsv)
        b = read_tsv(cs_tsv)
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
        agree = sum(1 for x, y in zip(a, b) if x == y)
        conf = Counter((x, y) for x, y in zip(a, b) if x != y)

        per_label = {}
        for want in sorted(set(a) | set(b)):
            tot = sum(1 for x in a if x == want)
            hit = sum(1 for x, y in zip(a, b) if x == want and y == want)
            if tot:
                per_label[want] = {"n_easyepi": tot, "matched": hit,
                                   "recall": round(hit / tot, 4)}

        results[name] = {
            "n_peaks": n,
            "exact_match": agree,
            "exact_match_frac": round(agree / n, 6) if n else None,
            "per_label": per_label,
            "top_disagreements": [
                {"easyepi": x, "chipseeker": y, "n": c}
                for (x, y), c in conf.most_common(8)
            ],
        }
        print(f"{name}: {agree:,}/{n:,} exact ({100 * agree / n:.4f}%)")
        for want, d in per_label.items():
            print(f"    {want:12s} n={d['n_easyepi']:6,d} recall={d['recall']:.4f}")

    out = ROOT / "results" / "annotation_agreement.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
