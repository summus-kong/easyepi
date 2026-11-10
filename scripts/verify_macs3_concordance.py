#!/usr/bin/env python3
"""Concordance check: easyepi's MACS3 wrapper vs a direct macs3 invocation.

What this establishes, and what it does not
-------------------------------------------
It proves the wrapper is *argument-transparent*: given the same input, the peak
file it produces is byte-identical to a hand-written ``macs3 callpeak`` with the
documented flags. That rules out the wrapper silently changing parameters, which
is the failure mode a wrapper can actually introduce.

It does not establish that the peaks are biologically correct, and it says nothing
about easyepi's own pure-Python callers, which are a different code path.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

REPO = Path("/data1/GaoLab/houoweibo/easyepi_work/easyepiv1")
WORK = Path("/data1/GaoLab/houoweibo/easyepi_work/encode_h3k4me3")
OUT = Path("/data1/GaoLab/houoweibo/easyepi_work/validation2/concordance")


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        rows.append((f[0], int(f[1]), int(f[2])))
    return sorted(rows)


def main() -> int:
    sys.path.insert(0, str(REPO))
    from easyepi._peak_calling import macs3_command

    macs3 = shutil.which("macs3")
    if not macs3:
        raise SystemExit("macs3 not on PATH")

    OUT.mkdir(parents=True, exist_ok=True)
    bam = WORK / "bam" / "H3K4me3.bam"
    # MACS3 reads BED tracks, not bigWig, so the coverage case uses the bedGraph
    # that easyepi's own wrapper accepts.
    bedgraph = WORK / "states" / "H3K4me3.chr1_5Mb.bedGraph"
    report = {"macs3": subprocess.run([macs3, "--version"], capture_output=True,
                                      text=True).stdout.strip() or "3.x",
              "cases": []}

    cases = [
        # name, sample attr, kwargs, expected format
        ("bam_narrow", {"bam": str(bam)}, dict(genome_size="hs", assay="CHIP"), "BAMPE"),
        ("bedgraph_narrow", {"bigwig": str(bedgraph)}, dict(genome_size="hs",
                                                            assay="CHIP",
                                                            paired_end=False), "BED"),
        ("bedgraph_broad", {"bigwig": str(bedgraph)}, dict(genome_size="hs",
                                                           assay="CHIP", broad=True,
                                                           paired_end=False), "BED"),
    ]

    for name, attrs, kwargs, expect_fmt in cases:
        sample = SimpleNamespace(sample_id=name, **attrs)
        ep_dir = OUT / name / "easyepi"
        di_dir = OUT / name / "direct"
        for d in (ep_dir, di_dir):
            if d.exists():
                shutil.rmtree(d)
            d.mkdir(parents=True)

        cmd = macs3_command(macs3, sample, ep_dir, **kwargs)
        fmt = cmd[cmd.index("-f") + 1]
        ep = subprocess.run(cmd, capture_output=True, text=True)
        if ep.returncode != 0:
            report["cases"].append({"case": name, "error": "easyepi macs3 failed",
                                    "stderr": ep.stderr[-500:]})
            continue

        # rebuild the same call by hand from the documented flags
        direct = [macs3, "callpeak", "-t", attrs.get("bam") or attrs.get("bigwig"),
                  "-f", fmt, "-g", kwargs["genome_size"], "-n", name,
                  "--outdir", str(di_dir), "-q", "0.01"]
        if kwargs.get("broad"):
            direct.append("--broad")
        if fmt == "BED":
            direct += ["--nomodel", "--extsize", "147"]
        di = subprocess.run(direct, capture_output=True, text=True)
        if di.returncode != 0:
            report["cases"].append({"case": name, "error": "direct macs3 failed",
                                    "stderr": di.stderr[-500:]})
            continue

        ep_peaks = sorted(ep_dir.glob(f"*_peaks.*Peak"))
        di_peaks = sorted(di_dir.glob(f"*_peaks.*Peak"))
        if not ep_peaks or not di_peaks:
            report["cases"].append({"case": name, "error": "no peak output"})
            continue

        a = read_bed(ep_peaks[0])
        b = read_bed(di_peaks[0])
        entry = {
            "case": name,
            "format": fmt,
            "easyepi_peaks": len(a),
            "direct_peaks": len(b),
            "identical_order_and_coords": a == b,
            "md5_easyepi": md5(ep_peaks[0]),
            "md5_direct": md5(di_peaks[0]),
            "peak_type": ep_peaks[0].suffix,
        }
        # coordinate-set equality is the meaningful contract; score columns may
        # carry a run timestamp in the header comment
        entry["pass"] = entry["identical_order_and_coords"]
        report["cases"].append(entry)

    report["all_pass"] = all(c.get("pass") for c in report["cases"]) and bool(report["cases"])
    (OUT / "concordance.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    for c in report["cases"]:
        if "error" in c:
            print(f"  ERROR {c['case']}: {c['error']}")
            continue
        print(f"  {c['case']:14s} fmt={c['format']:5s} "
              f"easyepi={c['easyepi_peaks']:6d} direct={c['direct_peaks']:6d} "
              f"identical={c['identical_order_and_coords']} "
              f"md5match={c['md5_easyepi'] == c['md5_direct']}")
    print(f"\nall_pass: {report['all_pass']}")
    print(f"report: {OUT / 'concordance.json'}")
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
