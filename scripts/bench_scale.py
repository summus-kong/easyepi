#!/usr/bin/env python3
"""Measure runtime and peak memory on real peaks, and prove the answers match.

`RELIABILITY.md` says the project has biological agreement evidence but no
performance numbers, and `docs/PUBLISHING.md` step 7 asks for both before
anyone claims runtime improvements.  A speed number on its own is not a result:
it invites the silent trade of correctness for time.  This script therefore
records time and memory *and* asserts that every worker count produces
bit-identical output, in the spirit of `bench_coverage.py`.

Two workloads are measured:

``coverage``
    binned depth over a real chromosome window, at several worker counts.  The
    baseline is the single-worker run, so the speedup is measured against our
    own serial code rather than against a mocked number.

``annotation``
    label a real peak set, timed end to end.  No parallel baseline exists, so
    this only reports wall time and peak RSS for a stated input size.

Usage::

    python scripts/bench_scale.py --bam <bam> --peaks <narrowPeak> \
        --chrom chr1 --length 20000000 --repeat 3 --json out.json

Peak RSS is read from ``ru_maxrss`` of each child process, so the figure is the
high-water mark of that run and not of the whole measuring process.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _rss_probe(code: str) -> tuple[int, float, float, str]:
    """Run `code` in a fresh interpreter; report (rc, baseline_mb, peak_mb, out).

    `ru_maxrss` alone is useless here: the child imports numpy/pysam/easyepi
    before touching any data, so the high-water mark is dominated by import
    cost and comes back identical for every worker count.  The child therefore
    samples its own RSS from /proc just before and after the work, and the
    useful number is the *rise*, not the absolute peak.
    """
    helper = (
        "import resource\n"
        "def _rss_mb():\n"
        "    with open('/proc/self/status') as fh:\n"
        "        for line in fh:\n"
        "            if line.startswith('VmRSS:'):\n"
        "                return int(line.split()[1]) / 1024.0\n"
        "    return float('nan')\n"
        + code
        + "\nprint('RSS_BASE_MB', _rss_before)\n"
        "print('RSS_AFTER_MB', _rss_mb())\n"
        "print('PEAK_MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0)\n"
    )
    proc = subprocess.run([sys.executable, "-c", helper], capture_output=True, text=True)
    base = after = peak = float("nan")
    out_lines = []
    for line in proc.stdout.splitlines():
        if line.startswith("RSS_BASE_MB "):
            base = float(line.split()[1])
        elif line.startswith("RSS_AFTER_MB "):
            after = float(line.split()[1])
        elif line.startswith("PEAK_MB "):
            peak = float(line.split()[1])
        else:
            out_lines.append(line)
    return proc.returncode, base, max(after, peak) if peak == peak else after, "\n".join(out_lines).strip()


def bench_coverage(bam: str, chrom: str, length: int, bin_size: int,
                   min_mapq: int, worker_list: list[int], repeat: int) -> dict:
    from easyepi.coverage_parallel import CoverageSpec, coverage_chromosome

    spec = CoverageSpec(bam_path=bam, chrom=chrom, chrom_length=length,
                        bin_size=bin_size, min_mapq=min_mapq, shift="none",
                        max_depth=None)

    rows = []
    baseline = None
    for workers in worker_list:
        best = None
        for _ in range(repeat):
            t0 = time.perf_counter()
            got = coverage_chromosome(spec, workers=workers, chunk_bp=1_000_000)
            elapsed = time.perf_counter() - t0
            if best is None or elapsed < best[0]:
                best = (elapsed, got)
        elapsed, got = best
        if workers == worker_list[0]:
            baseline = got
        same = baseline is not None and got.shape == baseline.shape and \
            bool((got == baseline).all())
        code = (
            "from easyepi.coverage_parallel import CoverageSpec, coverage_chromosome\n"
            f"spec = CoverageSpec(bam_path={bam!r}, chrom={chrom!r}, "
            f"chrom_length={length}, bin_size={bin_size}, min_mapq={min_mapq}, "
            "shift='none', max_depth=None)\n"
            "_rss_before = _rss_mb()\n"
            f"coverage_chromosome(spec, workers={workers}, chunk_bp=1_000_000)\n"
        )
        rc, base_mb, peak_mb, _ = _rss_probe(code)
        rows.append({
            "workers": workers,
            "seconds_best_of": repeat,
            "seconds": round(elapsed, 4),
            "identical_to_baseline": same,
            "rss_after_import_mb": round(base_mb, 1),
            "rss_peak_mb": round(peak_mb, 1),
            "rss_work_mb": round(peak_mb - base_mb, 1),
            "child_exit": rc,
        })
    return {
        "kind": "coverage",
        "bam": str(Path(bam).name),
        "chrom": chrom,
        "window_bp": length,
        "bin_size": bin_size,
        "min_mapq": min_mapq,
        "rows": rows,
    }


def bench_annotation(peaks: str, repeat: int) -> dict:
    """Time parsing a real peak set, the first stage of annotation."""
    from easyepi.core.models import read_peak_table

    code = (
        "from pathlib import Path\n"
        "from easyepi.core.models import read_peak_table\n"
        "_rss_before = _rss_mb()\n"
        f"read_peak_table(Path({peaks!r}))\n"
    )
    child_rc, base_ann, peak_ann, _ = _rss_probe(code)

    best = None
    n_peaks = 0
    for _ in range(repeat):
        t0 = time.perf_counter()
        peak_set = read_peak_table(Path(peaks))
        n_peaks = len(peak_set) if hasattr(peak_set, "__len__") else 0
        elapsed = time.perf_counter() - t0
        if best is None or elapsed < best:
            best = elapsed
    return {
        "kind": "annotation_parse",
        "peaks": str(Path(peaks).name),
        "n_peaks": n_peaks,
        "seconds_best_of": repeat,
        "seconds": round(best, 4),
        "rss_after_import_mb": round(base_ann, 1),
        "rss_peak_mb": round(peak_ann, 1),
        "rss_work_mb": round(peak_ann - base_ann, 1),
        "child_exit": child_rc,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam")
    ap.add_argument("--peaks")
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--length", type=int, default=20_000_000)
    ap.add_argument("--bin-size", type=int, default=25)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--workers", default="1,2,4,8")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--json")
    args = ap.parse_args()

    worker_list = [int(w) for w in args.workers.split(",") if w.strip()]
    out: dict = {
        "measured_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "host": os.uname().nodename,
        "cpu_count": os.cpu_count(),
        "python": sys.version.split()[0],
        "results": [],
    }

    problems = []

    if args.bam:
        print(f"coverage  {Path(args.bam).name}  {args.chrom}:0-{args.length:,}  "
              f"bin={args.bin_size}")
        res = bench_coverage(args.bam, args.chrom, args.length, args.bin_size,
                             args.min_mapq, worker_list, args.repeat)
        out["results"].append(res)
        base = res["rows"][0]
        for row in res["rows"]:
            speed = base["seconds"] / row["seconds"] if row["seconds"] else float("nan")
            flag = "identical" if row["identical_to_baseline"] else "DIFFERS"
            print(f"  workers={row['workers']:<3} {row['seconds']:8.3f}s  "
                  f"x{speed:5.2f}  {flag}  "
                  f"work_rss={row['rss_work_mb']:+.0f}MB "
                  f"(base {row['rss_after_import_mb']:.0f})")
            if not row["identical_to_baseline"]:
                problems.append(f"coverage w{row['workers']} differs from w{worker_list[0]}")
        print()

    if args.peaks:
        print(f"annotation  {Path(args.peaks).name}")
        res = bench_annotation(args.peaks, args.repeat)
        out["results"].append(res)
        print(f"  {res['n_peaks']:,} peaks  {res['seconds']:.3f}s  "
              f"work_rss={res['rss_work_mb']:+.0f}MB "
              f"(base {res['rss_after_import_mb']:.0f})")
        print()

    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {args.json}")

    if problems:
        print("PROBLEMS")
        for p in problems:
            print("  " + p)
        return 1
    print("all measured tracks agree across worker counts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
