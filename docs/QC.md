# QC commands, thresholds, and batch roll-up

Easyepi exposes QC through one CLI. Every path below shares a single metric and
verdict implementation (`easyepi.qc_payload.build_sample_payload`), so the same
sample gets the same numbers and the same verdict regardless of which command
produced the report.

## Commands

| Command | Assay | Extra outputs |
|---|---|---|
| `chip-qc` | ChIP-seq / CUT&Tag / CUT&RUN | `bam_qc`, `peak_qc`, `bigwig_qc`, `input_qc` |
| `atac-qc` | ATAC-seq | same + fragment-length structure |
| `histone-qc` | histone marks | same as `chip-qc` |
| `mnase-qc` | MNase-seq | same + nucleosome periodicity |
| `methylation-qc` | bisulfite / WGBS / RRBS | same + `methylation_qc` beta summary |
| `qc` | any, via `--assay` | generic; legacy BAM-only summary when `--assay` is omitted |
| `make-qc-report` | any, from a project | batch-level `qc_report.json` |
| `merge-qc-reports` | — | merges several `qc_report.json` |

All commands accept the same QC flags:

```
--project/--bam/--bigwig/--peak   inputs
--annotation GTF/GFF              enables TSS enrichment
--blacklist BED                   enables blacklist_fraction
--qc-thresholds JSON              override warn thresholds
--max-reads N                     cap reads inspected per metric (default 200000)
--figure-format png|pdf|svg       write figures
```

Assay names are normalised, so `ChIP`, `chip`, `ChIP-seq`, `CUT&Tag`, `CUT&RUN`,
`ATAC`, `ATAC-seq`, `histone`, `MNase`, `MNase-seq`, `WGBS`, `RRBS` and
`bisulfite` all resolve to the right metric profile.

## Metrics

| Group | Metrics |
|---|---|
| Alignment | `bam_total_reads`, `bam_mapping_rate` (+`mapping_rate_source`), `duplicate_rate` (+`duplicate_source`), `mitochondrial_fraction`, `median_read_length`, `mapq_ge_30`, `mapq_lt_10` |
| Signal | `frip`, `reads_in_peaks`, `signal_enrichment`, `mean_signal_in_peaks`, `tss_enrichment` |
| Peak support | `peaks_with_zero_reads`, `fraction_peaks_zero`, `median_reads_per_peak`, `max_reads_per_peak`, `top_decile_read_share` |
| Fragment structure | `median_fragment_length`, `fraction_nfr`, `fraction_mononucleosome`, `fraction_dinucleosome`, `mono_to_nfr_ratio` |
| Methylation | `mean_beta`, `is_bimodal` |
| Contamination | `blacklist_fraction` |

`mapping_rate_source` is `uninformative_no_unaligned_records` for BAMs that were
already filtered upstream; the ratio is then trivially 1.0 and must not be read
as an alignment success rate.

## Thresholds are advisory

Every rule breach produces `warn`. **`overall` is never `fail`**, because these
thresholds encode conventions, not biology — one metric must not reject a library
on its own. `failures` remains in the JSON as an always-empty list for backwards
compatibility.

```json
{
  "overall": "warn",
  "warnings": ["frip", "blacklist_fraction"],
  "warning_details": [
    {"metric": "frip", "status": "breached", "value": 0.1136,
     "threshold": 0.2, "op": "lt", "severity": "warn", "message": "low FRiP (<0.20)"}
  ],
  "checks": [ ... every evaluated rule, including "not_evaluated" ... ]
}
```

Metrics that were never computed are reported as `not_evaluated`, never as pass.

### Default rules

| Metric | Breach when | Note |
|---|---|---|
| `bam_mapping_rate` | < 0.80 | |
| `duplicate_rate` | > 0.50 | |
| `mitochondrial_fraction` | > 0.20 | lysis / cytoplasmic carry-over |
| `mapq_lt_10` | > 0.20 | |
| `frip` | < 0.20 | |
| `signal_enrichment` | < 2.0 | peak vs flank |
| `tss_enrichment` | < 2.0 | ENCODE-style |
| `fraction_peaks_zero` | > 0.50 | |
| `top_decile_read_share` | > 0.90 | signal carried by few peaks |
| `blacklist_fraction` | > 0.01 | |
| `fraction_nfr` | < 0.05 | ATAC tagmentation |
| `mono_to_nfr_ratio` | > 20.0 | under-digested |

### Overriding rules

Precedence (later wins, merged per metric — unspecified metrics keep defaults):

```
built-in defaults  <  --qc-thresholds rules.json  <  project [qc_thresholds]
```

```bash
easyepi chip-qc --project project.toml --out results --qc-thresholds rules.json
```

```json
{"frip": {"op": "lt", "value": 0.35, "message": "custom relaxed FRiP"},
 "tss_enrichment": {"op": "lt", "value": 3.0}}
```

```toml
[qc_thresholds.frip]
op = "lt"
value = 0.35
```

Malformed entries (unknown `op`, non-numeric `value`) and unreadable files are
skipped silently so one typo cannot break a batch.

## Migrating away from `fail`

If a downstream script gated on `verdict["failures"]`, it now never triggers.
Pick one:

```python
# Preferred: decide what is critical in your own code.
CRITICAL = {"bam_mapping_rate", "frip"}
breached = {c["metric"] for c in verdict["checks"] if c["status"] == "breached"}
if breached & CRITICAL:
    sys.exit(1)
```

```json
// Or: tighten thresholds so only real problems breach, then gate on warnings.
{"frip": {"op": "lt", "value": 0.10}}
```

## Multi-batch layout

One batch = one directory = one project. `merge-qc-reports` recurses for
`qc_report.json`, so the canonical name matters.

```
study/
├── qc_thresholds.json              # shared rules (optional)
├── batch_2025_01/
│   ├── project.toml                # set project_id here
│   └── qc_report.json              # <- from make-qc-report
├── batch_2025_03/
│   ├── project.toml
│   └── qc_report.json
└── qc_merged.json                  # output
```

```bash
for b in study/batch_*; do
  easyepi make-qc-report --project "$b/project.toml" --out "$b" \
      --qc-thresholds study/qc_thresholds.json
done

easyepi merge-qc-reports "study/batch_*/qc_report.json" --out study/qc_merged.json
```

Directory mode finds `qc_report.json`; to merge assay-command output, name the
files explicitly:

```bash
easyepi merge-qc-reports "batch_*/chip_seq_qc.json" --out merged.json
```

Merged output:

```json
{
  "module": "merged-qc-report",
  "n_reports": 2,
  "n_samples": 4,
  "batches": {"proj_A": 2, "proj_B": 2},
  "verdict_summary": {
    "samples": {"pass": 2, "warn": 2, "not_evaluated": 0},
    "warning_metrics": {"frip": 2}
  },
  "problems": []
}
```

Behaviour worth knowing:

- **Batch name** comes from `project_id`, then `module`, then the directory name.
  Set `project_id` in `project.toml` or several batches will collide as `qc-report`.
- **Duplicate sample ids** across batches are namespaced `batch_b::s1`; nothing is
  silently overwritten.
- **Broken reports** land in `problems` and do not abort the merge.
- **Both payload shapes** are accepted: `make-qc-report` writes `samples` as an
  object, assay commands write it as a list.
- Exit code is 1 only when nothing matched at all.
