"""Produce a verified download manifest for a 2-group x 2-replicate H3K4me3 pilot.

What the earlier query got wrong
--------------------------------
The first summary counted `(biosample.accession, biological_replicate)` from the
*search* projection, where `biological_replicate` comes back null. Every experiment
therefore looked like one replicate, and the whole assembly looked unusable.

ENCODE's convention for histone ChIP-seq is two biological replicates per assay.
Where both libraries point at the *same* biosample accession they are resequencing
of one sample, not replication -- which is exactly what ENCSR339TAQ turned out to
be. Distinguishing the two needs the replicate field fetched explicitly, which is
what this script does.

Output: for each candidate biosample, the experiments with a true 2-replicate
design, plus the alignment files to download. Nothing is downloaded here.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://www.encodeproject.org"
HEADERS = {"Accept": "application/json",
           "User-Agent": "easyepi-benchmark-probe/1.0 (research)"}


def get(url: str, timeout: int = 120, attempts: int = 3):
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout * (i + 1)) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            print(f"    retry {i+1}/{attempts}: {type(exc).__name__}",
                  file=sys.stderr)
            time.sleep(2)
    raise RuntimeError(f"{url}: {last}")


def experiment_replicates(acc: str) -> dict:
    """Resolve an experiment's replication structure from its full document.

    Returns the number of *biological* replicates, meaning distinct biosample
    accessions. Two libraries on one biosample are technical/sequencing repeats and
    cannot separate biological from technical variance, which is the entire point
    of the pilot.
    """
    doc = get(f"{BASE}/experiments/{acc}/?format=json&field=replicates")
    reps = [r for r in (doc.get("replicates") or []) if isinstance(r, dict)]
    samples = []
    for rep in reps:
        lib = rep.get("library") or {}
        bs = lib.get("biosample") or {}
        acc_bs = bs.get("accession") if isinstance(bs, dict) else bs
        if acc_bs:
            samples.append(acc_bs)
    distinct = sorted(set(samples))
    return {"accession": acc, "n_libraries": len(reps),
            "n_biological": len(distinct), "biosamples": distinct}


def alignment_files(acc: str, assembly: str = "GRCh38") -> list[dict]:
    """Released GRCh38 alignment files, plus a matching index if published."""
    doc = get(f"{BASE}/experiments/{acc}/?format=json&field=files")
    out = []
    for f in doc.get("files") or []:
        if not isinstance(f, dict):
            continue
        if f.get("status") != "released":
            continue
        if assembly not in (f.get("assembly") or ""):
            continue
        if f.get("output_type") != "alignments":
            continue
        if not str(f.get("file_format", "")).lower() == "bam":
            continue
        out.append({"accession": f.get("accession"),
                    "size_gb": round((f.get("file_size") or 0) / 1e9, 2),
                    "href": f.get("href"),
                    "output_type": f.get("output_type")})
    return out


def main() -> int:
    target = sys.argv[1] if len(sys.argv) > 1 else "H3K4me3"
    assembly = sys.argv[2] if len(sys.argv) > 2 else "GRCh38"

    params = {
        "type": "Experiment",
        "assay_title": "Histone ChIP-seq",
        "target.label": target,
        "assembly": assembly,
        "status": "released",
        "format": "json",
        "limit": "300",
        "field": "accession",
    }
    url = f"{BASE}/search/?{urllib.parse.urlencode(params)}"
    print(f"searching {target} / {assembly} ...", file=sys.stderr)
    graph = get(url).get("@graph", [])
    accs = [e.get("accession") for e in graph if e.get("accession")]
    print(f"{len(accs)} released experiments\n", file=sys.stderr)

    # Only experiments that resolve to >=2 distinct biosamples are candidates:
    # anything else cannot supply replicate-level variance.
    usable = []
    for acc in accs:
        info = experiment_replicates(acc)
        if info["n_biological"] >= 2:
            usable.append(info)
            print(f"  {acc}: {info['n_biological']} bio reps", file=sys.stderr)

    print(f"\n=== {len(usable)} experiments with >= 2 biological replicates ===")
    for info in usable:
        print(f"{info['accession']}  bio_reps={info['n_biological']}  "
              f"libs={info['n_libraries']}  biosamples={','.join(info['biosamples'])}")

    if len(usable) < 2:
        print("\nNOT enough for a two-group design.")
        return 1

    print("\n=== alignment files for the first two candidates ===")
    manifest = []
    for info in usable[:2]:
        files = alignment_files(info["accession"], assembly)
        total = sum(f["size_gb"] for f in files)
        print(f"\n{info['accession']}  ({len(files)} BAM, {total:.2f} GB)")
        for f in files:
            print(f"  {f['accession']}  {f['size_gb']:.2f} GB")
            manifest.append({"group": info["accession"], **f})
    print(f"\ntotal download: {sum(m['size_gb'] for m in manifest):.2f} GB "
          f"({len(manifest)} files)")

    out = sys.argv[3] if len(sys.argv) > 3 else "encode_manifest.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"manifest written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
