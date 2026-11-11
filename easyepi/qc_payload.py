"""Single source of truth for per-sample QC payload assembly.

Three call sites used to duplicate this logic and had drifted apart:

* ``assay_commands.run_assay_qc``  (``chip-qc`` / ``atac-qc`` / …)
* ``qc.run_qc``                    (generic ``qc`` command)
* ``project_commands.make_qc_report_command``

Each re-implemented "resolve paths -> compute assay metrics -> compute extended
metrics -> merge -> evaluate thresholds" slightly differently, which is exactly
how the batch report ended up disagreeing with the assay report for the same
sample. ``build_sample_payload`` is now the only place that logic lives; the call
sites keep their own extra outputs (``bigwig_qc``, ``peak_qc``, logging,
``resolved_config.json``) but share the metric + verdict core.

Every step is defensive: a metric that cannot be computed records an ``error``
entry instead of aborting, because one bad sample must not sink a batch.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .assay_metrics import compute_assay_metrics
from .qc_extended import compute_extended_metrics
from .qc_thresholds import evaluate_thresholds


def _existing(path: Any) -> str | None:
    """Return ``str(path)`` when it exists on disk, else None."""
    if not path:
        return None
    try:
        text = str(path)
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return None
    try:
        return text if Path(text).exists() else None
    except OSError:
        return None


def methylation_beta_values(table: Any) -> list[float] | None:
    """Extract the beta column from a parsed methylation table, if present."""
    if table is None:
        return None
    columns = getattr(table, "columns", None)
    if columns is None or "beta" not in columns:
        return None
    try:
        return [float(v) for v in table["beta"].dropna().tolist()]
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return None


def _beta_from_bam(bam_path: str, payload: dict[str, Any],
                   min_coverage: int = 5) -> list[float] | None:
    """Per-CpG beta values read from a bisulfite BAM, or None when it has no tags.

    The BAM is the documented WGBS input, so this is the normal path. Failure is
    recorded in ``payload`` rather than raised: a BAM without MM/ML or XM tags is a
    user-input problem, and one unreadable sample must not abort a batch QC run.
    The returned values are the per-CpG-dyad betas above ``min_coverage``, which is
    the same unit `beta_distribution` summarises for a table input.
    """
    try:
        from .methylation_bam import (calls_to_dataframe,
                                      cpg_level_methylation,
                                      read_methylation_calls)
        calls, stats = read_methylation_calls(bam_path)
        frame = calls_to_dataframe(calls, {})
        sites, site_stats = cpg_level_methylation(frame, min_coverage=min_coverage)
    except Exception as exc:
        # Includes "no methylation calls found in this BAM", which is the message a
        # user needs to see: the tag is missing, not the methylation.
        payload["methylation_error"] = f"{type(exc).__name__}: {exc}"
        return None
    payload["methylation_bam"] = {**stats, **site_stats}
    if sites.empty:
        return None
    return [float(v) for v in sites["beta"].dropna().tolist()]


def _peak_covers_bam(peak_path: str | None, bam_path: str | None) -> bool:
    """True when the peak set touches most of the BAM's aligned reference span.

    FRiP is only comparable to genome-wide conventions when peaks could have
    captured the reads. A peak set confined to one chromosome (a slice, a panel)
    against a whole-genome BAM cannot, so the caller must be told.
    """
    if not peak_path or not bam_path:
        return True
    try:
        import pysam

        peak_chroms = set()
        with open(peak_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                peak_chroms.add(line.split("\t", 1)[0])
                if len(peak_chroms) > 3:
                    return True
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            bam_chroms = {r for r, n in zip(bam.references, bam.lengths) if n > 1_000_000}
        if not bam_chroms:
            return True
        covered = len(peak_chroms & bam_chroms)
        return covered >= max(1, len(bam_chroms) // 2)
    except Exception:
        return True


def build_sample_payload(assay: str | None, *,
                         bam: Any = None,
                         bigwig: Any = None,
                         peak: Any = None,
                         methylation: Any = None,
                         annotation: Any = None,
                         blacklist: Any = None,
                         bam_qc: Mapping[str, Any] | None = None,
                         max_reads: int = 200_000,
                         rules: dict[str, dict[str, Any]] | None = None,
                         ) -> dict[str, Any]:
    """Assemble the shared QC payload for one sample.

    Returns a dict with:
      * ``assay_qc``   — merged assay-specific + extended metrics (may be empty),
      * ``qc_verdict`` — threshold verdict, or None when no assay is given.

    ``bam_qc`` is folded into the verdict payload under ``bam_qc`` so a low
    mapping rate is judged alongside FRiP, TSS enrichment and the rest.
    """
    payload: dict[str, Any] = {}
    if not assay:
        return {"assay_qc": payload, "qc_verdict": None}

    bam_path = _existing(bam)
    bigwig_path = _existing(bigwig)
    peak_path = _existing(peak)
    annotation_path = _existing(annotation)
    blacklist_path = _existing(blacklist)

    beta_values = None
    if assay.lower().startswith("methyl"):
        # Two sources, in the order the documentation promises. A bisulfite BAM
        # carries the calls on the record (MM/ML or Bismark XM), so `--bam` alone
        # is the normal case; a pre-called `.cov`/`.bedMethyl` table is accepted
        # for data that arrives that way. Without this, `qc --assay methylation
        # --bam x.bam` -- the exact command on the WGBS page -- reported generic
        # read flags and no beta distribution at all.
        methylation_path = _existing(methylation)
        if methylation_path:
            try:
                from .methylation import read_methylation_table
                beta_values = methylation_beta_values(read_methylation_table(methylation_path))
            except Exception as exc:
                payload["methylation_error"] = f"{type(exc).__name__}: {exc}"
        elif bam_path:
            beta_values = _beta_from_bam(bam_path, payload)

    try:
        payload.update(compute_assay_metrics(
            assay,
            bam=bam_path,
            bigwig=bigwig_path,
            peak=peak_path,
            methylation_values=beta_values,
            max_reads=max_reads,
        ))
    except Exception as exc:
        payload["assay_error"] = f"{type(exc).__name__}: {exc}"

    try:
        payload.update(compute_extended_metrics(
            assay,
            bam=bam_path,
            bigwig=bigwig_path,
            annotation=annotation_path,
            peak=peak_path,
            blacklist=blacklist_path,
            max_reads=max_reads,
        ))
    except Exception as exc:
        payload["extended_error"] = f"{type(exc).__name__}: {exc}"

    verdict_payload = dict(payload)
    if bam_qc:
        verdict_payload.setdefault("bam_qc", dict(bam_qc))
    # Tell the threshold layer when FRiP cannot be judged against genome-wide norms.
    if peak_path and bam_path and not _peak_covers_bam(peak_path, bam_path):
        verdict_payload["_partial_peak_coverage"] = True
    return {"assay_qc": payload,
            "qc_verdict": evaluate_thresholds(verdict_payload, rules)}


def merge_into_row(row: dict[str, Any], result: dict[str, Any]) -> dict[str, Any]:
    """Merge a :func:`build_sample_payload` result into an existing row dict."""
    payload = result.get("assay_qc") or {}
    if payload:
        merged = dict(row.get("assay_qc") or {})
        merged.update(payload)
        row["assay_qc"] = merged
    if result.get("qc_verdict") is not None:
        row["qc_verdict"] = result["qc_verdict"]
    return row
