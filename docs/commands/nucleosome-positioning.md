# easyepi nucleosome-positioning

```bash
easyepi nucleosome-positioning --bam paired.bam --chrom chr1 --start 0 --end 1000000 --out np/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--bam` | required |  | coordinate-sorted paired-end BAM (MNase-seq or ATAC-seq). Single-end input is refused: one cut site gives no fragment span, so the dyad cannot be located. Use `easyepi nucleosome` for a measured occupancy signal instead |
| `--chrom` | required |  | contig name as written in the BAM header, e.g. chr1 |
| `--start` | required |  | 0-based, inclusive start of the region to call |
| `--end` | required |  | 0-based, exclusive end of the region to call |
| `--bin-size` | optional | `10` | dyad density bin width in bp (default 10); spacing and fuzziness are reported in bp, not bins |
| `--min-length` | optional | `120` | minimum fragment length kept (default 120 bp) |
| `--max-length` | optional | `200` | maximum fragment length kept (default 200 bp) |
| `--smooth-bins` | optional | `3` | running-mean window in bins applied to the dyad density before peak calling (default 3) |
| `--mad-k` | optional | `2.0` | call threshold in scaled MADs above the median |
| `--min-spacing` | optional | `100` | minimum bp between accepted calls (default 100) |
| `--min-height` | optional |  | absolute floor on the smoothed dyad density a call must reach, in density units (default 0 = no absolute floor; the relative floor below still applies) |
| `--min-relative-height` | optional | `0.05` | fraction of the strongest peak a call must reach (default 0.05); guards against noise when MAD is 0 |
| `--fuzziness-fraction` | optional | `0.5` | fraction of a dyad cluster's signal the reported span must cover, measured outward from the peak (default 0.5) |
| `--fuzziness-radius` | optional | `150` | how far either side of a peak the span search may reach, in bp (default 150, about one nucleosome plus linker). Values near the neighbour spacing push the span to the cap |
| `--max-fragments` | optional | `2000000` | cap on fragments read from the region before density is built (default 2000000); raise it for deep sequencing |
| `--min-mapq` | optional | `30` | minimum mapping quality for a fragment to be used (default 30) |
| `--out` | required |  | output directory for positioning results |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
