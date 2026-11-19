"""Aim an ENCODE download for a two-cell-line H3K4me3 shape pilot.

Design decision
---------------
The pilot investigates shape/breadth change in H3K4me3 between two conditions. The
two groups are two cell lines (K562 vs HepG2): the biological difference is large
and well documented, both are ENCODE production cell lines with stable data
quality, and both have the 2-biological-replicate design ENCODE ships for histone
ChIP-seq. A cross-cell-line contrast is not the same question as a treatment
contrast -- it also moves peak positions -- which is recorded here rather than
discovered later: the shape test has to show it adds something over a breadth-only
and a position-only baseline, or it is measuring the wrong thing.

Why the previous attempt timed out
----------------------------------
It resolved ~300 experiments one at a time (one full document each, several seconds
per request on this link). The search endpoint can project the replicate biosample
accessions directly, so one request per cell line is enough to decide.

Biological vs technical replication
-----------------------------------
ENCODE histone experiments carry two *libraries*. Where both point at the same
biosample accession they are resequencing of one sample and provide no
between-sample variance -- measured case: ENCSR339TAQ, 2 libraries, 1 biosample.
This script counts distinct biosample accessions, not library entries.
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

#: Groups for the pilot. Two cell lines that ENCODE has profiled deeply.
GROUPS = ("K562", "HepG2")


def get(url: str, timeout: int = 180, attempts: int = 3):
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout * (i + 1)) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            print(f"    retry {i + 1}/{attempts}: {type(exc).__name__}",
                  file=sys.stderr)
            time.sleep(2)
    raise RuntimeError(f"{url}: {last}")


def search(target: str, cell_line: str, assembly: str, limit: int = 40):
    """One request per cell line, with the replicate projection inline."""
    params = [
        ("type", "Experiment"),
        ("assay_title", "Histone ChIP-seq"),
        ("target.label", target),
        ("assembly", assembly),
        ("status", "released"),
        ("biosample_ontology.term_name", cell_line),
        ("format", "json"),
        ("limit", str(limit)),
        # Projections keep the payload small enough to finish on this link.
        ("field", "accession"),
        ("field", "replicates.library.biosample.accession"),
        ("field", "replicates.library.biosample.term_name"),
        ("field", "date_released"),
    ]
    url = f"{BASE}/search/?{urllib.parse.urlencode(params)}"
    print(f"  querying {cell_line} / {target} ...", file=sys.stderr)
    return get(url).get("@graph", [])


def biosamples_of(doc: dict) -> list[str]:
    """Distinct biosample accessions reached through the replicate entries."""
    out = []
    for rep in doc.get("replicates") or []:
        if not isinstance(rep, dict):
            continue
        lib = rep.get("library") or {}
        bs = lib.get("biosample") or {}
        acc = bs.get("accession") if isinstance(bs, dict) else bs
        if acc:
            out.append(acc)
    return sorted(set(out))


def alignment_files(acc: str, assembly: str) -> list[dict]:
    """Released GRCh38 BAMs for one experiment."""
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
        if str(f.get("file_format", "")).lower() != "bam":
            continue
        out.append({"accession": f.get("accession"),
                    "size_gb": round((f.get("file_size") or 0) / 1e9, 3),
                    "href": f.get("href")})
    return out


def main() -> int:
    target = sys.argv[1] if len(sys.argv) > 1 else "H3K4me3"
    assembly = sys.argv[2] if len(sys.argv) > 2 else "GRCh38"
    out_path = sys.argv[3] if len(sys.argv) > 3 else "/tmp/encode_manifest.json"

    chosen = {}
    for cell_line in GROUPS:
        graph = search(target, cell_line, assembly)
        scored = []
        for doc in graph:
            acc = doc.get("accession")
            if not acc:
                continue
            bios = biosamples_of(doc)
            scored.append((len(bios), doc.get("date_released") or "", acc, bios))
        scored.sort(reverse=True)
        print(f"  {cell_line}: {len(scored)} experiments; "
              f"max bio reps = {scored[0][0] if scored else 0}", file=sys.stderr)
        for n, date, acc, bios in scored[:3]:
            print(f"    {acc}  bio_reps={n}  released={date[:10]}  "
                  f"biosamples={','.join(bios)}", file=sys.stderr)
        if scored and scored[0][0] >= 2:
            chosen[cell_line] = scored[0]

    print()
    if len(chosen) < 2:
        print("RESULT: could not find >=2 biological replicates in both cell lines.")
        for cell_line, row in chosen.items():
            print(f"  {cell_line}: {row[2]} with {row[0]} bio reps")
        return 1

    manifest = []
    print("=== download plan ===")
    for cell_line, (n_bio, date, acc, bios) in chosen.items():
        files = alignment_files(acc, assembly)
        total = sum(f["size_gb"] for f in files)
        print(f"\n{cell_line}  {acc}  bio_reps={n_bio}  "
              f"{len(files)} BAM  {total:.2f} GB  (released {date[:10]})")
        for f in files:
            print(f"  {f['accession']}  {f['size_gb']:.2f} GB")
            manifest.append({"group": cell_line, "experiment": acc,
                             "accession": f["accession"],
                             "size_gb": f["size_gb"],
                             "url": f"https://www.encodeproject.org{f['href']}"})

    total = sum(m["size_gb"] for m in manifest)
    print(f"\ntotal: {len(manifest)} files, {total:.2f} GB")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"manifest -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
