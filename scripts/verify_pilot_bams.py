"""Verify the downloaded pilot BAMs before any analysis touches them.

A benchmark whose inputs are silently wrong is worse than no benchmark, so each
file is checked for the properties the pilot actually depends on:

* **Readable and non-empty** -- a truncated transfer resumes as a valid-looking
  file if only its size was checked, and pysam is what will read it later.
* **Coordinate-sorted** -- peak calling and window counting assume it.
* **Same assembly** -- K562 and HepG2 must both be GRCh38. ENCODE ships hg19 for
  some older experiments and a mixed set produces zero overlaps that look like
  biology. Chromosome names are the cheap discriminator: GRCh38 ENCODE files use
  `chr1`.
* **Replicate depth** -- an uneven pair is a broken two-group design, so the read
  counts are reported per file rather than assumed comparable.
* **Index present** -- created here if missing, since `fetch` on an unindexed BAM
  returns nothing rather than failing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pysam

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/data1/GaoLab/houoweibo/easyepi_work/pilot_h3k4me3")

REPORT = []
print(f"{'group':7s} {'accession':12s} {'seqs':>6s} {'sorted':>7s} "
      f"{'chr1':>5s} {'mapped(M)':>10s} {'index':>6s}  notes")
print("-" * 92)

for bam_path in sorted(ROOT.rglob("*.bam")):
    group = bam_path.parent.name
    row = {"group": group, "accession": bam_path.stem, "path": str(bam_path)}
    notes = []
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            refs = list(bam.references)
            row["n_references"] = len(refs)
            row["has_chr1"] = "chr1" in refs
            row["has_plain_1"] = "1" in refs
            row["sorted"] = bool(bam.header.to_dict()
                                 .get("HD", {}).get("SO") == "coordinate")
            if row["has_plain_1"] and not row["has_chr1"]:
                notes.append("hg19-style names")
            elif not row["has_chr1"]:
                notes.append("no chr1!")

            mapped = total = 0
            for rec in bam.fetch(until_eof=True):
                total += 1
                if not rec.is_unmapped:
                    mapped += 1
                if total >= 400_000:      # enough to detect a truncated file
                    break
            row["sampled_reads"] = total
            row["sampled_mapped"] = mapped
            row["index"] = bam.has_index()
    except Exception as exc:  # noqa: BLE001
        row["error"] = f"{type(exc).__name__}: {exc}"
        notes.append("UNREADABLE")
        print(f"{group:7s} {bam_path.stem:12s}  ---  FAILED: {exc}")
        REPORT.append(row)
        continue

    if not row.get("index"):
        try:
            pysam.index(str(bam_path))
            row["index"] = True
            notes.append("index built")
        except Exception as exc:  # noqa: BLE001
            notes.append(f"index failed: {exc}")

    print(f"{group:7s} {row['accession']:12s} {row['n_references']:>6d} "
          f"{str(row['sorted']):>7s} {str(row['has_chr1']):>5s} "
          f"{row['sampled_mapped'] / 1e6:>10.2f} {str(row['index']):>6s}  "
          f"{', '.join(notes)}")
    REPORT.append(row)

print()
# The two properties that would silently ruin the comparison.
assemblies = {}
for row in REPORT:
    if "error" in row:
        continue
    assemblies.setdefault("GRCh38" if row.get("has_chr1") else "hg19", []).append(
        row["accession"])
print("assembly groups:", {k: len(v) for k, v in assemblies.items()})
if len(assemblies) > 1:
    print("  MIXED ASSEMBLIES -> do not compare across groups:",
          {k: v for k, v in assemblies.items()})

per_group = {}
for row in REPORT:
    if "error" not in row:
        per_group.setdefault(row["group"], []).append(row)
print()
for group, rows in sorted(per_group.items()):
    m = [r["sampled_mapped"] for r in rows]
    print(f"{group}: {len(rows)} file(s), sampled mapped reads "
          f"{min(m) / 1e6:.2f}-{max(m) / 1e6:.2f}M")
    if len(rows) < 2:
        print(f"  WARNING: {group} has fewer than 2 replicates")

bad = [r for r in REPORT if "error" in r]
(ROOT / "bam_check.json").write_text(json.dumps(REPORT, indent=2), encoding="utf-8")
print(f"\n{'OK' if not bad else 'PROBLEMS'}: "
      f"{len(REPORT) - len(bad)}/{len(REPORT)} readable")
print(f"report -> {ROOT / 'bam_check.json'}")
