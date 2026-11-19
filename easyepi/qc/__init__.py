"""QC summary entry point shared by the generic ``easyepi qc`` command.

When ``--assay`` is supplied this produces the *same* metric payload and
``qc_verdict`` block as the assay-specific commands (``atac-qc``, ``chip-qc`` …),
because both call :func:`easyepi.qc_payload.build_sample_payload`. Without
``--assay`` it stays the original lightweight BAM summary.
"""

from __future__ import annotations

import json
from pathlib import Path


def bam_qc_summary(bam_path) -> dict:
    """Return measured BAM QC; missing readers raise an explicit dependency error."""
    from ..core.io import bam_qc

    q = bam_qc(bam_path)
    summary = {"file": str(bam_path), **q, "mapping_rate": round(q.get("mapping_rate", 0.0), 4)}
    # Preserve provenance so a filtered BAM's trivially-1.0 ratio is never
    # presented as an alignment success rate.
    summary["mapping_rate_source"] = q.get("mapping_rate_source")
    return summary


def run_qc(args) -> int:
    """Generic QC entry point; assay-aware when ``--assay`` is given."""
    from ..qc_payload import build_sample_payload
    from ..qc_thresholds import load_rules, summarize_verdicts

    from ..qc_inputs import resolve_inputs

    # Resolve by sample name, not list position: positional pairing
    # mismatched tracks as soon as one sample lacked a BigWig or the flags
    # were written in a different order.
    inputs = resolve_inputs(args)
    assay = getattr(args, "assay", None)
    if assay:
        # Validate and canonicalise before anything reads it. A typo used to be
        # passed straight to the payload, which returned an empty assay_qc and a
        # report that looked successful: "chp" silently produced no FRiP. The
        # canonical label is also what the metric table is keyed by, and what the
        # input check uses to say which of the supplied tracks it will ignore.
        from ..qc_assay import canonical_assay, check_inputs
        assay = canonical_assay(assay)
        for note in check_inputs(assay, {
                sample: {kind: path for kind, path in tracks.items() if path}
                for sample, tracks in inputs.items()}):
            print(f"[easyepi] {note}")
    rules = load_rules(getattr(args, "qc_thresholds", None))

    result: dict[str, object] = {"samples": {}}
    if assay:
        result["assay"] = assay
    rows: list[dict] = []

    # Loaded once for the whole run: a TSS table is identical for every sample, and
    # re-parsing it per sample would dominate the runtime of a small cohort.
    #
    # ``--annotation`` is the documented flag and the one in the QC command page,
    # but TSS enrichment used to ignore it entirely and accept only ``--tss``, so
    # the obvious invocation silently produced no score. It is now read as the
    # TSS source too. ``--tss`` still wins when both are given, because a
    # hand-curated table is a deliberate statement and a GTF is a convenience.
    _tss_rows = None
    _tss_path = getattr(args, "tss", None) or getattr(args, "annotation", None)
    if _tss_path:
        from ..tss_enrichment import read_tss_table

        try:
            _tss_rows = read_tss_table(_tss_path)
        except (ValueError, OSError) as exc:
            # A bad annotation must not abort a QC run whose other metrics are
            # fine; say why TSS enrichment is missing and continue.
            _tss_rows = None
            print(f"[easyepi] TSS enrichment skipped: {exc}")
        else:
            print(f"[easyepi] TSS enrichment: {len(_tss_rows)} start sites "
                  f"from {_tss_path}")

    for sample_id, tracks in inputs.items():
        bam = tracks.get("bam")
        bigwig = tracks.get("bigwig")
        peak = tracks.get("peak")
        summary = bam_qc_summary(bam) if bam else {}
        metrics: dict[str, object] = {"sample_id": sample_id, "bam": summary}

        if assay:
            metrics["assay"] = assay
            built = build_sample_payload(
                assay,
                bam=bam,
                bigwig=bigwig,
                peak=peak,
                annotation=getattr(args, "annotation", None),
                blacklist=getattr(args, "blacklist", None),
                bam_qc=summary,
                max_reads=getattr(args, "max_reads", 200000),
                rules=rules,
                tss_rows=_tss_rows,
                tss_flank=int(getattr(args, "tss_flank", 2000)),
            )
            metrics["assay_qc"] = built["assay_qc"]
            metrics["qc_verdict"] = built["qc_verdict"]
            rows.append({"sample_id": sample_id, "qc_verdict": built["qc_verdict"]})

        result["samples"][sample_id] = metrics

    if assay:
        result["verdict_summary"] = summarize_verdicts(rows)

    # A report with no samples is a silent failure: it looks successful but
    # describes nothing. Fail instead of writing an empty file.
    if not result["samples"]:
        raise ValueError(
            "no samples were resolved; check --bam/--bigwig/--peak or "
            "--input-dir (no BAM, BigWig or peak file was recognised)")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[easyepi] QC summary written to {out}")

    # Figures are rendered only after the JSON is on disk. Every plot helper in
    # `easyepi._plots` reads the written payload (`json.load(path)`), so rendering
    # first found no data and produced nothing -- the reason a QC report carried
    # no images while a second, default-format pass wrote PDFs on its own.
    figures_made: list[str] = []
    try:
        from .._plots import render_module_figures, selected_formats

        formats = selected_formats(args)
        # The report embeds figures and a PDF cannot be embedded, so guarantee a
        # raster format by testing what the resolved list actually contains.
        want_html = not getattr(args, "no_html", False)
        raster = {"png", "svg", "jpg", "jpeg", "gif"}
        if formats and want_html and not any(
                str(f).lower().lstrip(".") in raster for f in formats):
            formats = tuple(formats) + ("png",)
        if formats:
            figures_made = render_module_figures("qc", out, formats) or []
    except Exception as exc:  # noqa: BLE001 - figures are secondary to the numbers
        print(f"[easyepi] figure rendering skipped: {type(exc).__name__}: {exc}")

    # The HTML report sits beside the JSON and reads it, rather than being a second
    # source of truth: whichever is opened, the numbers are the same ones.
    if not getattr(args, "no_html", False):
        from ..qc_html import render_qc_html

        html_path = out.with_suffix(".html")
        # Only the figures this run just wrote. Scanning the whole `figures/`
        # directory instead pulled in plots from any other run that shared the
        # output folder, so a MNase report embedded the ATAC run's charts.
        try:
            render_qc_html(result, html_path, figure_paths=figures_made)
        except Exception as exc:  # a broken report must not fail a finished run
            print(f"[easyepi] HTML report skipped: {type(exc).__name__}: {exc}")
        else:
            print(f"[easyepi] QC report written to {html_path}")
    return 0
