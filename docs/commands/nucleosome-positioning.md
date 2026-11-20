# easyepi nucleosome-positioning

```bash
easyepi nucleosome-positioning --bam paired.bam --chrom chr1 --start 0 --end 1000000 --out np/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--bam` | required |  | coordinate-sorted paired-end BAM (MNase-seq or ATAC-seq) |
| `--chrom` | required |  |  |
| `--start` | required |  |  |
| `--end` | required |  |  |
| `--bin-size` | optional | `10` |  |
| `--min-length` | optional | `120` | minimum fragment length kept (default 120 bp) |
| `--max-length` | optional | `200` | maximum fragment length kept (default 200 bp) |
| `--smooth-bins` | optional | `3` |  |
| `--mad-k` | optional | `2.0` | call threshold in scaled MADs above the median |
| `--min-spacing` | optional | `100` | minimum bp between accepted calls (default 100) |
| `--min-height` | optional |  |  |
| `--min-relative-height` | optional | `0.05` | fraction of the strongest peak a call must reach (default 0.05); guards against noise when MAD is 0 |
| `--fuzziness-fraction` | optional | `0.5` |  |
| `--fuzziness-radius` | optional | `150` |  |
| `--max-fragments` | optional | `2000000` |  |
| `--min-mapq` | optional | `30` |  |
| `--out` | required |  | output directory for positioning results |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
