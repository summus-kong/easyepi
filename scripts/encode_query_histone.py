"""Query ENCODE for H3K4me3 experiments usable as a two-group, multi-replicate set.

Why this exists
---------------
The shape-vs-breadth differential prototype needs what the current host data does
not have: two groups with **>= 3 biological replicates each**, on the **same
assembly**. The four BAMs already on the host are one replicate per mark, from four
different experiments, and one of them (H3K4me1) is hg19 while the rest are
GRCh38 -- mixing them would silently produce zero overlaps.

This script only *queries and reports*; nothing is downloaded here. It prints a
candidate table so the download can be aimed rather than guessed.

ENCODE's API is slow from this host: a bare `limit=1` search already streams for
several seconds. Requests therefore ask for a narrow `field=` projection and retry
once on timeout.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict

BASE = "https://www.encodeproject.org"
HEADERS = {
    "Accept": "application/json",
    "User-Agent": "easyepi-benchmark-probe/1.0 (research; contact: local)",
}

#: Fields requested explicitly. `field=` keeps the payload small, which matters
#: because the full experiment documents are megabytes each on this link.
EXP_FIELDS = (
    "accession",
    "assay_title",
    "target.label",
    "biosample_ontology.term_name",
    "biosample_ontology.classification",
    "assembly",
    "replicates.library.biosample.accession",
    "replicates.library.biosample.biological_replicate",
    "replicates.library.biosample.donor",
)


def fetch(url: str, *, attempts: int = 3, timeout: int = 90):
    """GET JSON with a bounded retry.

    A single timeout is common here because the link streams slowly; retrying the
    same URL with a longer timeout is what actually succeeds, so the timeout grows
    per attempt rather than the request changing.
    """
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout * (i + 1)) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            print(f"    retry {i + 1}/{attempts} ({type(exc).__name__}: {exc})",
                  file=sys.stderr)
            time.sleep(2)
    raise RuntimeError(f"could not fetch {url}: {last}")


def search_experiments(target: str, assembly: str, limit: int = 200):
    params = {
        "type": "Experiment",
        "assay_title": "Histone ChIP-seq",
        "target.label": target,
        "assembly": assembly,
        "status": "released",
        "format": "json",
        "limit": str(limit),
    }
    query = "&".join(f"{k}={urllib.parse.quote(str(v))}" for k, v in params.items())
    for f in EXP_FIELDS:
        query += f"&field={urllib.parse.quote(f)}"
    url = f"{BASE}/search/?{query}"
    print(f"querying {target} / {assembly} ...", file=sys.stderr)
    data = fetch(url)
    return data.get("@graph", [])


def main() -> int:
    target = sys.argv[1] if len(sys.argv) > 1 else "H3K4me3"
    assembly = sys.argv[2] if len(sys.argv) > 2 else "GRCh38"

    experiments = search_experiments(target, assembly)
    print(f"{len(experiments)} released {target} experiments on {assembly}\n")

    # Group by biosample term so a "group" is a real biological unit rather than
    # an arbitrary split of one cell type's replicates.
    by_term = defaultdict(list)
    for exp in experiments:
        term = (exp.get("biosample_ontology") or {}).get("term_name") or "?"
        reps = exp.get("replicates") or []
        # Count distinct biological replicates, not library entries: one biosample
        # can be sequenced twice and that is technical, not biological, replication.
        bio = set()
        for rep in reps:
            sample = (rep.get("library") or {}).get("biosample") or {}
            key = (sample.get("accession"), sample.get("biological_replicate"))
            if sample.get("accession"):
                bio.add(key)
        by_term[term].append((exp.get("accession"), len(bio)))

    rows = []
    for term, items in by_term.items():
        best = max(items, key=lambda x: x[1])
        rows.append((best[1], term, best[0], len(items)))
    rows.sort(reverse=True)

    print(f"{'bio_reps':>8}  {'experiments':>11}  biosample")
    print("-" * 78)
    for reps, term, acc, n in rows:
        flag = "  <-- usable (>=3)" if reps >= 3 else ""
        print(f"{reps:>8}  {n:>11}  {term}   [{acc}]{flag}")

    usable = [r for r in rows if r[0] >= 3]
    print(f"\n{len(usable)} biosample group(s) with >= 3 biological replicates")
    if len(usable) >= 2:
        print("candidate two-group set (pick by biology, not by this list):")
        for reps, term, acc, _n in usable[:6]:
            print(f"  {term:34s} {reps} reps   {acc}")
    else:
        print("NOT enough groups for a two-group pilot on this assembly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
