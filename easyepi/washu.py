"""WashU Epigenome Browser tracks visualization (optional backend).

Generates a WashU Epigenome Browser compatible track configuration (a JSON
track-hub-like description) plus the clickable browser URL, so Easyepi samples
can be visualized interactively in the WashU browser. Tracks use the native
WashU formats: ``bigWig`` (signal) and ``bigBed``/``bed`` (features).

The module is a lightweight, optional integration: it never requires the WashU
browser itself, only that the referenced data files are readable over HTTP
(hosted by the user / a file server), since the WashU browser fetches remote
files by URL.
"""

from __future__ import annotations

import json
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence



@dataclass(frozen=True)
class Sample:
    """Sample helper."""
    name: str
    bigwig: str
    group: str = "default"
    peak_file: str | None = None

    @property
    def stem(self) -> str:
        """Filename stem without directories or a trailing suffix."""
        return Path(self.name).stem


def read_manifest(path: str) -> List[Sample]:
    """Read sample_id/name, bigwig and group columns for WashU track export."""
    rows = [line for line in Path(path).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if not rows:
        return []
    header = rows[0].split("\t")
    samples: List[Sample] = []
    if "sample_id" in header or "name" in header:
        for row in csv.DictReader(rows, delimiter="\t"):
            name = row.get("sample_id") or row.get("name") or ""
            bigwig = row.get("bigwig") or row.get("track") or ""
            group = row.get("group") or row.get("condition") or "default"
            if not name or not bigwig:
                raise ValueError("manifest requires sample_id/name and bigwig/track")
            samples.append(Sample(name, bigwig, group, row.get("peak_file") or row.get("peak") or row.get("peaks")))
    else:
        for line in rows:
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError("manifest rows require name and bigwig")
            samples.append(Sample(fields[0], fields[1], fields[2] if len(fields) > 2 else "default",
                                  fields[3] if len(fields) > 3 and fields[3] else None))
    if len({sample.name for sample in samples}) != len(samples):
        raise ValueError("manifest sample names must be unique")
    return samples

# WashU genome name -> our genome label. The browser knows hg38 by this name.
_SUPPORTED_GENOMES = {"hg38": "hg38", "hg19": "hg19", "mm10": "mm10"}


def build_washu_hub(samples: Sequence[Sample], out_dir: str,
                    genome: str = "hg38",
                    base_url: Optional[str] = None) -> Dict:
    """Build a WashU browser track-hub description for a set of samples.

    Args:
        samples: list of :class:`Sample` (name, bigwig, group, optional peak_file).
        out_dir: where to write the hub JSON (and the data links).
        genome: reference genome WashU name (hg38 default).
        base_url: optional base URL (http/https) under which the track data are
            served. If None, tracks are emitted as file paths (for a local
            arrangement), and a note is produced that the browser needs URLs.
    Returns:
        dict with hub_path, browser_url, n_samples, n_tracks, genome.
    """
    if genome not in _SUPPORTED_GENOMES:
        raise ValueError(f"unsupported genome '{genome}' (supported: {list(_SUPPORTED_GENOMES)})")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    tracks = []
    for smp in samples:
        label = smp.stem
        track = {
            "name": label,
            "type": "bigWig",
            "url": _to_url(smp.bigwig, base_url),
            "color": smp.group,
            "height": 50,
            "autoScale": True,
        }
        tracks.append(track)
        if getattr(smp, "peak_file", None) and Path(smp.peak_file).exists():
            tracks.append({
                "name": f"{label}-peaks",
                "type": "bigBed",
                "url": _to_url(smp.peak_file, base_url),
                "color": "#d62728",
                "height": 40,
            })

    hub = {
        "genome": genome,
        "sample_mode": "grouped",
        "tracks": tracks,
    }
    hub_path = out / f"washu_hub_{genome}.json"
    hub_path.write_text(json.dumps(hub, indent=2), encoding="utf-8")

    browser_url = make_washu_url(samples, genome=genome, base_url=base_url)
    return {"hub_path": str(hub_path), "browser_url": browser_url,
            "n_samples": len(samples), "n_tracks": len(tracks), "genome": genome,
            "notes": ("WashU browser fetches tracks by URL. If data files are not "
                      "publicly served, host them first or use the track hub file with "
                      "a local WashU instance.")}


def make_washu_url(samples: Sequence[Sample], genome: str = "hg38",
                   region: str = "chr1:1,000,000-2,000,000",
                   base_url: Optional[str] = None) -> str:
    """Build a WashU Epigenome Browser URL that loads the given tracks.

    URL form (WashU browser accepts ``bam``/``bigwig``/``bigbed`` as comma
    separated track URLs):
        https://epigenomegateway.wustl.edu/browser/?genome=hg38&bigwig=<url1>,<url2>...
    """
    track_urls = [_to_url(s.bigwig, base_url) for s in samples]
    param = ",".join(track_urls)
    return (f"https://epigenomegateway.wustl.edu/browser/?genome={genome}"
            f"&bigwig={param}&region={region}")


def _to_url(path: str, base_url: Optional[str]) -> str:
    """Return the browser-fetchable URL for a data file."""
    p = Path(path)
    name = p.name
    if base_url:
        return f"{base_url.rstrip('/')}/{name}"
    # No base -> emit the absolute local path so a local/self-hosted WashU file
    # server or the user can map it to a URL. Mark it clearly.
    return str(p)


def write_washu_readme(out_dir: str, samples: Sequence[Sample],
                       genome: str = "hg38") -> str:
    """Write a short README explaining how to view the tracks in WashU."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lines = [
        "# WashU Epigenome Browser tracks",
        "",
        "Easyepi exported WashU-compatible tracks for the following samples:",
        "",
        "| sample | signal (bigWig) | peak (bigBed) |",
        "|---|---|---|",
    ]
    for s in samples:
        pk = Path(s.peak_file).name if getattr(s, "peak_file", None) else "-"
        lines.append(f"| {s.stem} | {Path(s.bigwig).name} | {pk} |")
    lines += [
        "",
        "## How to view",
        "1. Open the WashU Epigenome Browser:",
        "   `https://epigenomegateway.wustl.edu/browser/`",
        "2. Set the genome (hg38) and paste the track file URLs into the browser.",
        "3. The browser fetches signal (bigWig) and peak (bigBed) tracks by URL.",
        "",
        "> Note: WashU fetches tracks over HTTP(S). If your files are on the "
        "cluster (not publicly served), host them on a web/file server first, or "
        "run a self-hosted WashU instance with a file server mapping these paths.",
    ]
    readme = out / "WASHU_README.md"
    readme.write_text("\n".join(lines), encoding="utf-8")
    return str(readme)
