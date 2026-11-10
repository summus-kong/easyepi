# easyepi bam2bw

```bash
easyepi bam2bw sample.bam -o tracks/ --assay atac
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| bam | required |  | coordinate-sorted BAM file(s) |
| `-o`, `--out` | required |  | output directory for one track per sample |
| `--output-file` | optional |  | explicit BigWig path; only valid with a single BAM |
| `--assay` | optional |  | sets sensible defaults for shift/bin-size/normalize |
| `--shift` | optional |  | tn5 applies the +4/-5 cut-site shift (ATAC/CUT&Tag) |
| `--bin-size` | optional |  | bin width in bp |
| `--normalize` | optional |  | value scaling; cpm/rpkm make samples comparable |
| `--chrom` | optional |  | restrict to one chromosome (default: all) |
| `--min-mapq` | optional | `30` |  |
| `--max-depth` | optional |  | cap per-base depth before binning |
| `--name` | optional |  | override a sample name (default: BAM file stem) |
| `--jobs`, `-j` | optional |  | worker processes for coverage accumulation; default: every core the scheduler allows, 1 = in-process |
| `--manifest` | optional |  | write a JSON manifest of the produced tracks |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
