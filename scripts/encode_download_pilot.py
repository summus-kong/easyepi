"""Download the H3K4me3 pilot BAMs for a K562 vs HepG2 shape comparison.

Design properties that matter for a benchmark that must be reproducible:

* **Idempotent.** A file whose size already matches the manifest is skipped, so a
  re-run after an interruption resumes instead of restarting 6.7 GB.
* **Resumable.** ``curl -C -`` continues a partial transfer. ENCODE's download
  endpoint honours range requests, and on this link a 1 GB file does get
  interrupted.
* **Verified, not assumed.** Every file's byte size is checked against the manifest
  after transfer; a mismatch is recorded as a failure rather than left to surface
  as a truncated BAM deep inside peak calling.
* **Grouped by cell line**, because the pilot's whole premise is a two-group
  contrast and a flat directory loses which BAM belongs to which group.

The manifest is the output of ``scripts/encode_manifest_pilot.py``; nothing here
queries ENCODE for metadata, so the file list is fixed before the transfer starts.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/data1/GaoLab/houoweibo/easyepi_work/pilot_h3k4me3")
MANIFEST = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/tmp/encode_pilot_manifest.json")


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.2f} {unit}"
        n /= 1024
    return f"{n:.2f} PB"


def fetch(url: str, target: Path, expected_gb: float, attempts: int = 5) -> dict:
    """Download one file, resuming across attempts.

    Returns a status dict rather than raising: one unreachable file must not lose
    the other seven, and the caller records which ones need a retry.
    """
    expected_bytes = int(expected_gb * 1e9)
    # ENCODE reports file_size in bytes; the manifest rounded to 3 decimals of GB,
    # so compare with a tolerance of 1 MB rather than demanding byte equality.
    tolerance = 1_000_000

    if target.exists():
        have = target.stat().st_size
        if abs(have - expected_bytes) <= tolerance:
            return {"file": target.name, "status": "skipped",
                    "bytes": have, "note": "already complete"}

    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(target.suffix + ".part")

    for attempt in range(1, attempts + 1):
        cmd = ["curl", "-L", "--fail", "--retry", "3", "--retry-delay", "5",
               "--connect-timeout", "60", "--max-time", "7200",
               "-C", "-", "-o", str(part), url]
        print(f"    attempt {attempt}/{attempts}: {target.name}", flush=True)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        size = part.stat().st_size if part.exists() else 0
        if proc.returncode == 0 and abs(size - expected_bytes) <= tolerance:
            part.replace(target)
            return {"file": target.name, "status": "downloaded", "bytes": size}
        if size and abs(size - expected_bytes) <= tolerance:
            # curl can report a non-zero code after completing a resumed transfer.
            part.replace(target)
            return {"file": target.name, "status": "downloaded", "bytes": size}
        print(f"      curl exit {proc.returncode}, have {human(size)} "
              f"of {human(expected_bytes)}", flush=True)
        if proc.stderr:
            print("      " + proc.stderr.strip().splitlines()[-1][:200], flush=True)
        time.sleep(min(60, 5 * attempt))

    return {"file": target.name, "status": "failed",
            "bytes": part.stat().st_size if part.exists() else 0,
            "error": f"after {attempts} attempts"}


def main() -> int:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    print(f"manifest: {len(manifest)} files, "
          f"{sum(m['size_gb'] for m in manifest):.2f} GB")
    print(f"destination: {DEST}\n")

    DEST.mkdir(parents=True, exist_ok=True)
    results = []
    for i, entry in enumerate(manifest, 1):
        group = entry["group"]
        acc = entry["accession"]
        target = DEST / group / f"{acc}.bam"
        print(f"[{i}/{len(manifest)}] {group}/{acc} "
              f"({entry['size_gb']:.2f} GB)", flush=True)
        res = fetch(entry["url"], target, entry["size_gb"])
        res.update({"group": group, "experiment": entry["experiment"],
                    "accession": acc, "expected_gb": entry["size_gb"]})
        results.append(res)
        print(f"    -> {res['status']} ({human(res.get('bytes', 0))})", flush=True)

    # Recording the outcome next to the data is what makes a later benchmark run
    # auditable: it says exactly which files were used and whether any were partial.
    status_path = DEST / "download_status.json"
    status_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    ok = [r for r in results if r["status"] in ("downloaded", "skipped")]
    bad = [r for r in results if r["status"] == "failed"]
    print(f"\n=== {len(ok)}/{len(results)} complete, {len(bad)} failed ===")
    for r in bad:
        print(f"  FAILED {r['group']}/{r['accession']} "
              f"{human(r.get('bytes', 0))} of {r['expected_gb']:.2f} GB")
    print(f"status written to {status_path}")

    # Per-group tallies, because an uneven group is a broken two-group design.
    for group in sorted({r["group"] for r in results}):
        rows = [r for r in results if r["group"] == group]
        total = sum(r.get("bytes", 0) for r in rows)
        print(f"  {group}: {len(rows)} file(s), {human(total)}")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
