"""Download ENCODE peak files with bounded retries that never abort the batch.

Why this is Python rather than a shell loop: the remote Run wrapper stops the
whole script as soon as one command exits non-zero, so a `for` loop over curl
never gets past its first failure. Driving curl from Python keeps every failure
local, and a 1.2 MB peak file is small enough that a patient retry is cheap even
on the degraded link we measured (ENCODE ~8 KB/s, GEO ~855 KB/s).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

OUT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10/enc_peaks")

# (experiment, file accession) for GRCh38 pseudoreplicated H3K27ac peaks.
FILES = [
    ("ENCSR864OOO", "ENCFF837ELM"),
    ("ENCSR954JMZ", "ENCFF149QUM"),
    ("ENCSR837DVF", "ENCFF299FWI"),
    ("ENCSR779ISJ", "ENCFF271HJI"),
]

ATTEMPTS = 8


def fetch(acc: str, dest: Path) -> tuple[bool, str]:
    url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
    last = ""
    for attempt in range(1, ATTEMPTS + 1):
        proc = subprocess.run(
            [
                "curl", "-sS", "-L",
                "--connect-timeout", "15", "--max-time", "180",
                "-o", str(dest), "-w", "%{http_code}:%{size_download}",
                url,
            ],
            capture_output=True,
            text=True,
        )
        code = (proc.stdout or "").strip() or f"exit{proc.returncode}"
        last = code
        if code.startswith("200:") and dest.exists() and dest.stat().st_size > 0:
            return True, f"attempt {attempt}: {code}"
        # A dangling 307 means the portal redirected but the S3 hop dropped.
        if dest.exists():
            dest.unlink()
    return False, last


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ok = 0
    for experiment, acc in FILES:
        dest = OUT / f"{acc}.bed.gz"
        if dest.exists() and dest.stat().st_size > 0:
            print(f"{experiment}/{acc}: already present")
            ok += 1
            continue
        good, detail = fetch(acc, dest)
        size = dest.stat().st_size if dest.exists() else 0
        print(f"{experiment}/{acc}: {'OK' if good else 'FAILED'} "
              f"({detail}, {size} bytes)")
        ok += int(good)

    print(f"\ndownloaded {ok}/{len(FILES)}")
    for p in sorted(OUT.glob("*.bed.gz")):
        print(f"  {p.name} {p.stat().st_size} bytes")
    return 0 if ok == len(FILES) else 1


if __name__ == "__main__":
    raise SystemExit(main())
