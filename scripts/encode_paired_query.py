"""Locate real paired-end ENCODE BAMs and report their mate structure.

What the API actually provides (measured, not assumed)
-----------------------------------------------------
* `paired_end` is NOT a filterable field on `type=File`: every query carrying it
  returns HTTP 404 with an empty `@graph`. A File document fetched directly also
  reports `paired_end = None`.
* `type=Experiment` search works (463 ATAC-seq GRCh38 experiments).
* An experiment document's `files` key holds **full embedded file objects**, not path
  strings, so alignment files can be read without a second request per file.

So end-type must be established from the BAM itself (`flag & 0x1`), which is the only
reliable signal. This script collects candidates and records whatever metadata exists;
`s3_verify_paired.py` proves mate structure from the file.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

BASE = "https://www.encodeproject.org"
OUT = Path("/data1/GaoLab/houoweibo/easyepi_work/pilot_paired")


def get(url: str) -> dict:
    proc = subprocess.run(
        ["curl", "-sS", "-L", "--fail", "--max-time", "180",
         "-H", "Accept: application/json", url],
        capture_output=True, text=True, check=True,
    )
    return json.loads(proc.stdout)


def experiments(assay: str, limit: int = 12) -> list[dict]:
    url = (f"{BASE}/search/?type=Experiment&format=json&limit={limit}"
           f"&assay_title={assay.replace(' ', '%20')}&status=released")
    return get(url).get("@graph", [])


def bam_files(acc: str) -> list[dict]:
    exp = get(f"{BASE}/experiments/{acc}/?format=json")
    out = []
    for f in (exp.get("files") or []):
        if not isinstance(f, dict):
            continue
        if f.get("file_format") != "bam":
            continue
        if f.get("output_type") not in ("alignments", "unfiltered alignments"):
            continue
        out.append({
            "accession": f.get("accession"),
            "output_type": f.get("output_type"),
            "size_mb": round(int(f.get("file_size") or 0) / 1e6, 1),
            "assembly": f.get("assembly"),
            "paired_end": f.get("paired_end"),
            "run_type": f.get("run_type"),
            "href": f.get("href"),
            "biological_replicates": f.get("biological_replicates"),
        })
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report = {}

    for assay in ("ATAC-seq", "Histone ChIP-seq"):
        print(f"\n=== {assay} ===")
        exps = experiments(assay)
        print(f"  experiments returned: {len(exps)}")
        picked = []
        for e in exps:
            acc = e.get("accession")
            if not acc:
                continue
            try:
                bams = bam_files(acc)
            except subprocess.CalledProcessError:
                continue
            if not bams:
                continue
            picked.append({"experiment": acc, "bams": bams})
            sizes = [b["size_mb"] for b in bams]
            print(f"  {acc}: {len(bams)} bam(s)  sizes {sizes[:4]} MB  "
                  f"paired_end fields {[b['paired_end'] for b in bams][:4]}")
            if len(picked) >= 5:
                break
        report[assay] = picked

    (OUT / "paired_candidates.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    total = sum(len(v) for v in report.values())
    print(f"\n{total} experiments with BAMs written -> {OUT / 'paired_candidates.json'}")
    print("\nNOTE: end-type is NOT determined here. `paired_end` is absent from these "
          "documents, so it must be read from each BAM's FLAG (0x1).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
