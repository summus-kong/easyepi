"""CLI driver for the pure-Python chromatin-state model."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np


def _read_signal_per_bin(path: str, chrom: str, n_bins: int, resolution: int,
                       origin: int = 0) -> np.ndarray:
    """Fold a bedGraph into per-bin signal for one chromosome.

    A bedGraph line can straddle several bins, so the overlap length is used to
    weight the value; simply taking the line start would bias signal toward the
    left edge of every bin.

    ``origin`` is the genomic coordinate of bin 0. It was previously assumed to be
    0, so a run with a non-zero ``--start`` folded the genome from coordinate 0
    into a window that was declared to begin at ``--start``, shifting every bin.
    """
    out = np.zeros(n_bins, dtype=float)
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] != chrom:
                continue
            try:
                start, end, value = int(f[1]), int(f[2]), float(f[3])
            except ValueError:
                continue
            if end <= start:
                continue
            lo = max(0, (start - origin) // resolution)
            hi = min(n_bins - 1, (end - 1 - origin) // resolution)
            if hi < 0 or lo >= n_bins:
                continue
            for b in range(lo, hi + 1):
                b_start = origin + b * resolution
                b_end = b_start + resolution
                overlap = min(end, b_end) - max(start, b_start)
                if overlap > 0:
                    out[b] += value * overlap / (b_end - b_start)
    return out


def _version() -> str | None:
    try:
        import easyepi
        return getattr(easyepi, "__version__", None)
    except Exception:  # noqa: BLE001 - reporting must never break a run
        return None


def _module_version(name: str) -> str | None:
    try:
        from importlib import metadata
        return metadata.version(name)
    except Exception:  # noqa: BLE001
        return None


def run_histone_states(args) -> int:
    """Binarize -> LearnModel -> MakeSegmentation -> Annotate."""
    from easyepi import histone_states as hs

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    chrom = args.chrom
    resolution = int(args.resolution)
    n_bins = int(args.end - args.start) // resolution
    if n_bins <= 0:
        raise ValueError("--end must exceed --start by at least one bin")

    from easyepi.viz.panels import expand_csv_values

    marks: dict[str, str] = {}
    for spec in expand_csv_values(args.signal):
        if "=" not in spec:
            raise ValueError(f"--signal expects NAME=PATH, got {spec!r}")
        name, path = spec.split("=", 1)
        marks[name] = path
    if not marks:
        raise ValueError("at least one --signal NAME=PATH is required")

    # ---- 1. Binarize -------------------------------------------------------
    binarised = []
    binarise_info = {}
    for name, path in marks.items():
        values = _read_signal_per_bin(path, chrom, n_bins, resolution,
                                       origin=int(args.start))
        if args.binarize == "peaks":
            track = (values > 0).astype(np.int8)
            binarise_info[name] = {"mode": "peaks", "n_on": int(track.sum())}
        else:
            track, info = hs.poisson_threshold(values, cutoffs=tuple(args.poisson_cutoff))
            binarise_info[name] = {"mode": "poisson", **info, "n_on": int(track.sum())}
        binarised.append(track)
    matrix = np.column_stack(binarised)
    mark_names = list(marks)

    on_per_mark = matrix.sum(axis=0)
    if np.all(on_per_mark == 0):
        raise ValueError("every mark is all-zero after binarisation; "
                         "raise --poisson-cutoff or use --binarize peaks")

    # ---- 2. LearnModel ----------------------------------------------------
    state = hs.learn_states(matrix, int(args.states), seed=int(args.seed),
                            max_iter=int(args.max_iter),
                            n_init=int(getattr(args, "n_init", 8) or 8))
    emissions = hs.emission_probabilities(state)

    # ---- 3. MakeSegmentation ---------------------------------------------
    paths = hs.segment(matrix, state)
    names = hs.default_state_names(int(args.states))
    bed = hs.states_to_bed(paths, chrom, resolution, names, origin=int(args.start))
    bed_path = out / "segmentation.bed"
    with bed_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for row in bed:
            w.writerow(row)

    # ---- 4. Annotate ------------------------------------------------------
    summary = hs.state_summary(paths, emissions, names, resolution)

    # Provenance: which files produced this model. Without it a re-run cannot be
    # shown to have used the same input, and two runs that disagree cannot be
    # traced to a data change rather than a seed change. The hash is over file
    # bytes only; it is not a substitute for archiving the inputs.
    def _digest(path: str) -> dict:
        import hashlib
        target = Path(path)
        if not target.exists():
            return {"path": str(path), "sha256": None, "bytes": None,
                    "note": "input not present at run time"}
        h = hashlib.sha256()
        with target.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return {"path": str(target.resolve()), "sha256": h.hexdigest(),
                "bytes": target.stat().st_size}

    provenance = {
        "easyepi_version": _version(),
        "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "hmmlearn_version": _module_version("hmmlearn"),
        "inputs": {name: _digest(path) for name, path in marks},
    }

    (out / "states.json").write_text(json.dumps({
        "chrom": chrom, "start": int(args.start), "end": int(args.end),
        "provenance": provenance,
        "resolution": resolution, "n_bins": n_bins,
        "marks": mark_names,
        "binarisation": binarise_info,
        "n_states": int(args.states),
        "fitted_on_bins": state["fitted_on"],
        "n_init": state.get("n_init"),
        "loglik": state.get("loglik"),
        "transition": np.round(state["model"].transmat_, 6).tolist(),
        "states": summary,
    }, indent=2), encoding="utf-8")

    tsv = out / "state_summary.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["state", "n_bins", "coverage_bp", "fraction"] + mark_names)
        for row in summary:
            w.writerow([row["state"], row["n_bins"], row["coverage_bp"],
                        row["fraction"]] + row["emissions"])

    print(f"[easyepi] histone-states: {len(names)} states over {n_bins} bins "
          f"({len(marks)} marks) -> {out}")
    for row in summary[:5]:
        print(f"    {row['state']:5s} {row['fraction'] * 100:5.1f}%  emissions={row['emissions']}")
    return 0
