# easyepi qc

```bash
easyepi qc --assay atac --bam sample.bam --out qc.json
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--bam` | optional |  | BAM file(s); repeatable and/or comma-separated |
| `--bigwig` | optional |  | BigWig signal(s); repeatable and/or comma-separated. Enables assay-level metrics |
| `--peak` | optional |  | BED/narrowPeak(s); repeatable and/or comma-separated. Enables FRiP and peak metrics |
| `--input-dir` | optional |  | folder scanned for *.bam / *.bw / *.bigWig and peak files; each sample becomes <stem>. Saves listing every file by hand |
| `--out` | required |  | output report path |
| `--assay` | optional |  | assay label: atac, chip, histone, mnase or methylation; CUT&Tag/CUT&RUN and common synonyms are accepted. Enables the same verdict block as the assay-specific commands |
| `--project` | optional |  | Project TOML or TSV sample manifest |
| `--annotation` | optional |  | GTF/GFF annotation for TSS enrichment |
| `--blacklist` | optional |  | ENCODE-style blacklist BED |
| `--qc-thresholds` | optional |  | JSON file overriding QC warn thresholds |
| `--max-reads` | optional | `200000` | cap on reads inspected per metric (default: 200000) |
| `--verbose` | optional |  |  |
| `--quiet` | optional |  |  |
| `--figure-format` | optional |  | Optional export format; repeat to export multiple formats |

Generated from the live parser by `scripts/gen_docs.py`.
