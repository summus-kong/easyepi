"""Fetch small ENCODE BAMs with a retry loop that survives per-attempt failures.

Why not curl in a shell loop
----------------------------
The link to `encode-public.s3.amazonaws.com` answers 200 when healthy but drops
connections intermittently ("Empty reply from server"). A bash `for` loop run through
the job wrapper aborts on the first non-zero exit, so the retries never happened.

A Python loop controls its own error handling: each attempt is wrapped in try/except,
so one dropped connection costs one attempt rather than the whole run. The per-attempt
timeout is short on purpose -- a stalled socket is cheaper to abandon and retry than to
wait out.

Nothing here decides end-type. That is read from the BAM FLAG afterwards, because the
ENCODE metadata does not expose `paired_end`.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://www.encodeproject.org"
OUT = Path("/data1/GaoLab/houoweibo/easyepi_work/pilot_paired/small")

# Small alignment BAMs found by sorting `type=File&output_type=alignments` by size.
# Mixed assemblies on purpose: the goal is a workable paired-end file, not a specific
# genome.
CANDIDATES = [
    "ENCFF015YFD", "ENCFF503DMU", "ENCFF948IFA", "ENCFF011ITM",
    "ENCFF234DZO", "ENCFF637NOS", "ENCFF291HAS", "ENCFF268UNA",
    "ENCFF050UNA", "ENCFF567AVK", "ENCFF080LJI", "ENCFF397QTR",
]

ATTEMPTS = 8
TIMEOUT = 60


def resolve(acc: str) -> str | None:
    """The portal 307s to a signed S3 URL; return that URL."""
    url = f"{BASE}/files/{acc}/@@download/{acc}.bam"
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return resp.geturl()
    except Exception:
        return None


def fetch(acc: str, dest: Path) -> int:
    """Download one file, retrying; returns bytes written (0 on failure)."""
    for attempt in range(1, ATTEMPTS + 1):
        url = resolve(acc)
        if not url:
            time.sleep(1.5)
            continue
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp, \
                    open(dest, "wb") as fh:
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    fh.write(chunk)
            if dest.stat().st_size > 0:
                return dest.stat().st_size
        except Exception:
            pass
        time.sleep(1.5)
    return 0


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    got = {}
    for acc in CANDIDATES:
        dest = OUT / f"{acc}.bam"
        if dest.exists() and dest.stat().st_size > 0:
            got[acc] = dest.stat().st_size
            print(f"{acc} already present {got[acc]} bytes")
            continue
        n = fetch(acc, dest)
        if n:
            got[acc] = n
            print(f"{acc} OK {n} bytes")
        else:
            print(f"{acc} FAILED after {ATTEMPTS} attempts")

    (OUT / "downloaded.json").write_text(json.dumps(got, indent=2), encoding="utf-8")
    print(f"\n{len(got)} file(s) downloaded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
