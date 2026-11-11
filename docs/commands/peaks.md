# easyepi peaks

```bash
easyepi peaks tracks/sample.bw -o peaks/ --assay atac
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| bam | required |  | treatment BAM(s), coordinate-sorted |
| `-c`, `--control` | optional |  | control BAM(s); repeatable or comma-separated |
| `-o`, `--out` | required |  | output directory |
| `--assay` | optional |  | sets format and broad/narrow defaults |
| `--format` | optional |  | MACS3 format: BAM, BAMPE, BED, BEDPE (default from --assay) |
| `--genome` | optional | `hs` | MACS3 genome size: hs, mm, dm or a number (default hs) |
| `--qvalue` | optional | `0.01` | q-value cutoff (default 0.01) |
| `--extsize` | optional |  |  |
| `--shift-size` | optional |  |  |
| `--broad` | optional |  | broad peaks (default from --assay) |
| `--macs-arg` | optional |  | extra MACS3 options passed through verbatim, e.g. --macs-arg="--nomodel --extsize 200"; repeatable |
| `--keep-dup-all` | optional |  | keep all duplicate reads (MACS3 default is auto) |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
