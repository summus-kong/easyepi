#!/usr/bin/env python3
"""End-to-end real-data validation for tf-annotate and footprint.

Runs entirely from already-present local data:
  * peaks      : ENCODE H3K4me3 narrowPeak (GRCh38, chr-named)
  * genome     : indexed GRCh38 FASTA (genome.fa + genome.fa.fai)
  * motifs     : JASPAR PFMs fetched once with the stdlib (cached under --cache)
  * signal     : ENCODE H3K4me3 bigWig for the aggregate profile

Everything is reported as measured numbers plus explicit pass/fail gates, so a
silent degradation (e.g. every motif skipped, or zero sites found) fails loudly.
"""
from __future__ import annotations

import json
import ssl
import sys
import urllib.request
from pathlib import Path

import numpy as np

WORK = Path("/data1/GaoLab/houoweibo/easyepi_work/encode_h3k4me3")
GENOME = Path("/data1/GaoLab/houoweibo/reference/human/genome.fa")
REPO = Path("/data1/GaoLab/houoweibo/easyepi_work/easyepiv1")
CACHE = WORK / "motifs"
OUT = WORK / "validation"

# JASPAR core PFMs for factors expected at active promoters / enhancers.
JASPAR_IDS = {
    "CTCF": "MA0139.1",
    "GATA1": "MA0035.4",
    "GATA2": "MA0036.3",
    "NFYA": "MA0060.3",
    "SP1": "MA0079.5",
    "YY1": "MA0095.2",
    "E2F1": "MA0024.4",
    "MYC": "MA0147.3",
    "MAX": "MA0058.3",
    "NRF1": "MA0506.1",
    "EKLF": "MA0012.2",
    "TBP": "MA0108.2",
    "ELK1": "MA0028.2",
    "CREB1": "MA0018.4",
    "ATF2": "MA0090.2",
}


def fetch(url: str, timeout: int = 60) -> bytes:
    ctx = ssl.create_default_context()
    try:
        import certifi
        ctx.load_verify_locations(certifi.where())
    except Exception:
        pass
    req = urllib.request.Request(url, headers={"User-Agent": "easyepi-validation"})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        return r.read()


def build_motif_file() -> Path:
    """Download each JASPAR PFM once and concatenate into one .pfm."""
    CACHE.mkdir(parents=True, exist_ok=True)
    out = CACHE / "jaspar_subset.pfm"
    if out.exists() and out.stat().st_size > 1000:
        print(f"  cached: {out}")
        return out
    blocks = []
    for tf, mid in JASPAR_IDS.items():
        try:
            text = fetch(f"https://jaspar.elixir.no/api/v1/matrix/{mid}.pfm").decode()
        except Exception as exc:
            print(f"  ! {tf} ({mid}) fetch failed: {exc}")
            continue
        body = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith(">")]
        if len(body) != 4:
            print(f"  ! {tf} ({mid}) unexpected PFM shape ({len(body)} rows)")
            continue
        blocks.append(f">{mid} {tf}\n" + "\n".join(body) + "\n")
    if not blocks:
        raise SystemExit("no JASPAR motifs could be fetched")
    out.write_text("".join(blocks), encoding="utf-8")
    print(f"  wrote {out} ({len(blocks)} motifs)")
    return out


def main() -> int:
    sys.path.insert(0, str(REPO))
    OUT.mkdir(parents=True, exist_ok=True)

    peaks = WORK / "peaks" / "H3K4me3_peaks.narrowPeak"
    bigwig = WORK / "data" / "H3K4me3.bw"
    report = {"inputs": {"peaks": str(peaks), "genome": str(GENOME),
                         "bigwig": str(bigwig)}}

    print("== 1. real JASPAR motif set ==")
    motif_file = build_motif_file()

    from easyepi.tf_annotate import (annotate_sequences, fetch_sequences,
                                     load_motif_database, read_peak_intervals)

    print("\n== 2. tf-annotate on real H3K4me3 peaks ==")
    sub = OUT / "figures"
    sub.mkdir(parents=True, exist_ok=True)

    all_records = read_peak_intervals(peaks)
    records = all_records[:1000]
    motifs = load_motif_database(motif_file)
    print(f"  peaks in file: {len(all_records)}, using: {len(records)}")
    print(f"  motifs loaded: {len(motifs)}")

    seqs = fetch_sequences(records, GENOME)
    lengths = [len(s) for s in seqs]
    print(f"  sequences fetched: {len(seqs)}, length min/median/max = "
          f"{min(lengths)}/{sorted(lengths)[len(lengths)//2]}/{max(lengths)}")
    assert len(seqs) == len(records) and min(lengths) > 0

    ann = annotate_sequences(seqs, motifs, fpr=1e-4, top=5)
    usable = ann["n_motifs"] - len(ann["skipped_motifs"])
    n_hit = len({r["seq_index"] for r in ann["per_peak"]})
    print(f"  motifs usable: {usable} / skipped: {len(ann['skipped_motifs'])}")
    if ann["skipped_motifs"]:
        for name, rec in ann["skipped_motifs"].items():
            print(f"    skipped {name}: best fpr {rec['best_achievable_fpr']:.2e}")
    print(f"  peaks with >=1 site: {n_hit}/{len(seqs)} "
          f"({100 * n_hit / len(seqs):.1f}%)")
    print("  top motifs by peak coverage:")
    for row in ann["summary"][:8]:
        print(f"    {row['tf']:8s} {row['pct_peaks_with_site']:5.1f}%  "
              f"thr={row['score_threshold']:7.3f}  best={row['best_score']}")

    report["tf_annotate"] = {
        "n_peaks": len(records),
        "n_motifs": len(motifs),
        "n_motifs_usable": usable,
        "skipped": ann["skipped_motifs"],
        "n_peaks_with_site": n_hit,
        "pct_peaks_with_site": round(100 * n_hit / len(seqs), 2),
        "top_motifs": ann["summary"][:8],
    }

    print("\n== 3. aggregate profile on real H3K4me3 bigWig ==")
    from easyepi.footprint import (extract_profiles, fig_footprint,
                                   fig_tf_annotation, footprint_profile,
                                   footprint_significance, local_signal)

    chrom, start, end = "chr1", 0, 5_000_000
    signal = local_signal(str(bigwig), chrom, start, end)
    print(f"  signal bins read: {len(signal):,} "
          f"(mean={signal.mean():.4f}, max={signal.max():.2f})")

    # Centre the profile on each peak's midpoint. Using the interval start would
    # place the reference at the left edge of every peak, which shows up as a
    # monotonic ramp instead of the summit that H3K4me3 must have.
    loci = [r for r in records
            if r["chrom"] == chrom and start <= (r["start"] + r["end"]) // 2 < end]
    offsets = [(r["start"] + r["end"]) // 2 - start for r in loci]
    print(f"  target loci in {chrom}:{start}-{end}: {len(offsets)}")
    assert offsets, "no target loci inside the probe window"

    flank = 100
    bg = []
    for pos in offsets:
        for shift in (5_000, 25_000, 100_000):
            for cand in (pos + shift, pos - shift):
                if 0 <= cand - flank and cand + flank < len(signal):
                    bg.append(cand)
                    break
    print(f"  shifted background loci: {len(bg)}")

    prof = footprint_profile(signal, offsets, flank=flank, background_sites=bg)
    sig = footprint_significance(extract_profiles(signal, offsets, 0, flank),
                                 extract_profiles(signal, bg, 0, flank))
    print(f"  centre-vs-shoulder score: {prof['score']:.4f} (n={prof['n_sites']})")
    print(f"  target mean dip: {sig['target_mean_dip']:.4f}, "
          f"background: {sig.get('background_mean_dip')}")
    print(f"  Mann-Whitney p: {sig.get('p_mannwhitney')}")
    print(f"  Wilcoxon p (dip>0): {sig.get('p_wilcoxon')}")

    t = np.asarray(prof["target_profile"], dtype=float)
    b = np.asarray(prof["background_profile"], dtype=float)
    mid = len(t) // 2
    centre_vs_bg = float(t[mid] - b[mid]) if b.size else float("nan")
    print(f"  centre: target={t[mid]:.3f} background={b[mid]:.3f} "
          f"diff={centre_vs_bg:+.3f}")

    report["footprint"] = {
        "region": f"{chrom}:{start}-{end}",
        "n_target": len(offsets), "n_background": len(bg),
        "score": prof["score"], "target_mean_dip": sig["target_mean_dip"],
        "background_mean_dip": sig.get("background_mean_dip"),
        "p_mannwhitney": sig.get("p_mannwhitney"),
        "p_wilcoxon": sig.get("p_wilcoxon"),
        "centre_vs_background_centre": round(centre_vs_bg, 4),
    }

    fig_footprint(prof, title=f"H3K4me3 aggregate profile {chrom}:{start}-{end}",
                  save=str(sub / "footprint"), formats=["png"])
    fig_tf_annotation(ann["summary"], n_peaks=len(seqs),
                      title=f"Known-motif annotation, {len(seqs)} H3K4me3 peaks",
                      save=str(sub / "tf_annotation"), formats=["png"])
    figs = sorted(p.name for p in sub.glob("*.png"))
    print(f"  figures written: {figs}")
    report["figures"] = figs

    (OUT / "validation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    # The informative direction for a histone mark at its own peaks is centre vs
    # background centre: H3K4me3 should be enriched at the peak summit, not
    # depleted. ``footprint_profile``'s ``score`` is centre relative to the target's
    # OWN shoulders, a different (and for a broad mark, uninformative) contrast, so
    # the gate tests the direction that the biology predicts.
    print("\n== 4. gates ==")
    gates = {
        "sequences_fetched": len(seqs) == len(records),
        "some_motif_usable": usable > 0,
        "site_rate_plausible": 50 <= report["tf_annotate"]["pct_peaks_with_site"] <= 100,
        "signal_nonzero": float(signal.max()) > 0,
        "centre_enriched_over_background": centre_vs_bg > 0,
        "figures_written": len(figs) == 2,
    }
    for name, ok in gates.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\nreport: {OUT / 'validation_report.json'}")
    return 0 if all(gates.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
