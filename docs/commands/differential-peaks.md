# easyepi differential-peaks

```bash
easyepi differential-peaks --peaks ctrl1,ctrl2 --bam-dir bams --peak-dir peaks --group-a ctrl1 --out diff.tsv
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--peaks` | required |  | repeatable and/or comma-separated. Three forms: SAMPLE=BAM:PEAK (explicit), SAMPLE=BAM (peak file found via --peak-dir), or SAMPLE alone (both found via --bam-dir and --peak-dir), e.g. --peaks ctrl1,ctrl2 |
| `--bam-dir` | optional |  | folder holding <sample>.bam for bare-label --peaks entries |
| `--peak-dir` | optional |  | folder holding <sample>.narrowPeak / .broadPeak / .bed |
| `--group-a` | optional |  | comma-separated sample names forming group A |
| `--design` | optional |  | design TSV (sample_id + group) instead of --group-a |
| `--min-overlap` | optional | `1` | bp overlap needed to merge peaks into one region (default 1) |
| `--min-mapq` | optional | `30` |  |
| `--out` | required |  | output TSV file path |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
