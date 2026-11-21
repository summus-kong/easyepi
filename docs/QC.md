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
--figure-format pdf|svg|png       write figures (pdf is the default)
--no-figures                      do not write figures
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
| Library | `duplicate_fraction` (+`duplicate_status`, `duplicate_note`), `library_complexity`, `mitochondrial_fraction` (+`mitochondrial_reads`, `mitochondrial_denominator`, `mitochondrial_mapped_records`, `mitochondrial_status`, `mitochondrial_pattern`) |
| Methylation | `mean_beta`, `is_bimodal` |
| Contamination | `blacklist_fraction` |

`mapping_rate_source` is `uninformative_no_unaligned_records` for BAMs that were
already filtered upstream; the ratio is then trivially 1.0 and must not be read
as an alignment success rate.

The two library fractions can both be *undefined*, and each says so rather than
reporting zero:

- `duplicate_status` is `measured` only when some record carries the duplicate flag.
  A BAM straight out of an aligner has it clear on every record, so a bare `0.0`
  would claim a clean library that was never assessed; that case is reported as
  `not_marked` with an explanatory `duplicate_note`. `library_complexity` (distinct
  templates over reads inspected) is always defined and is the fallback signal.
- `mitochondrial_fraction` requires `--mito-pattern`. Contig names differ between
  references (`chrM`, `MT`, `chrMT`, accessions), so none is assumed; without the
  flag the value is `null` with `mitochondrial_status: not_applicable`. When no
  contig matches the pattern the status is `no_contig_matched`, again `null` — never
  `0.0`, because a zero here reads as "free of mitochondrial reads" when the
  question was never answered.

`mitochondrial_fraction` is counted through the BAM index over the **whole file**,
not from the scanned sample: chrM is the last contig of a coordinate-sorted BAM, so a
scan bounded by `--max-reads` never reaches it. The reported triple is auditable —
`mitochondrial_reads / mitochondrial_mapped_records`, with
`mitochondrial_denominator: index_mapped_records` naming the population. That
denominator counts every mapped record in the file, including secondary alignments,
so it is not the same population as `library_complexity`'s scan; do not divide one by
the other.

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

## Which reads are discarded

Every metric below is computed after the same baseline filter, so a count that looks
low is usually a filter rather than a data problem. Each stage reports what it
skipped (`skipped`, `records_seen`, `used`) so the loss is visible instead of
implied.

### Baseline filter, applied everywhere

| Test | Records dropped | Why |
|---|---|---|
| `is_unmapped` | yes | no position to contribute |
| `is_secondary` | yes | one primary alignment per read; counting all placements would inflate coverage by repeat copy number |
| `is_supplementary` | yes | chimeric segments are not independent observations |
| `is_duplicate` | yes, when the BAM carries the flag | flagged duplicates are PCR copies. A BAM that was never duplicate-marked has the flag clear on every record and nothing is dropped; `library.duplicate_status` reports that case as `not_marked` rather than as zero duplicates |
| `mapping_quality < --min-mapq` | yes | see the defaults below |

`is_qcfail` is **not** tested. That is deliberate: the flag is set by aligner-specific
conventions and is frequently absent, so filtering on it would silently behave
differently across input files.

### `--min-mapq` defaults, and why they differ

| Command | Default | Note |
|---|---|---|
| `bam2bw` | 30 | |
| `differential-peaks` | 30 | |
| `nucleosome-positioning` | 30 | |
| `dmr` | **10** | bisulfite aligners assign systematically lower MAPQ; 30 would discard a large share of valid calls |

The `dmr` exception is pinned by `tests/test_min_mapq_invariants.py`: raising it to 30
fails three assertions, including one that requires the help text to keep stating the
reason.

### Fragment-length window (`nucleosome-positioning` only)

`--min-length 120` / `--max-length 200` select single-nucleosome fragments: 147 bp of
protected DNA plus linker. A fragment below 120 bp is usually a sub-nucleosome or a
dyad-less fragment and would place a false dyad; above 200 bp it is increasingly a
di-nucleosome, whose midpoint is the linker between two nucleosomes rather than
either dyad. Both bounds are flags because digestion strength varies by protocol.

### One record per fragment

Paired-end input is deduplicated to one observation per fragment: the read-1 record
alone is kept. The fragment and cut-site paths both do this, so a paired library is
not counted twice, and `library_complexity` (distinct templates over reads) remains
comparable between single-end and paired-end input.

### Blacklist

`--blacklist BED` drops signal in the listed intervals and reports the fraction lost
as `blacklist_fraction`. It is opt-in: no blacklist is assumed, because the ENCODE
lists are genome-build specific and applying the wrong build removes real signal.

### Where the two coverage paths are checked against each other

`bam2bw` selects a path from the assay: cut sites for `--assay atac`
(`coverage_readends`), accumulated fragment depth otherwise (`coverage_parallel`).
The two disagreed on one filter until recently — the fragment path dropped flagged
duplicates and the cut-site path did not — so a duplicate-marked BAM produced a
different track depending only on which assay chose the path.
`tests/test_coverage_path_parity.py` now fails if they drift apart again.

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
