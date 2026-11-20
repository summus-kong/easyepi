# easyepi nucleosome

```bash
easyepi nucleosome --bigwig mnase.bw --out nuc/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--bigwig` | required |  | occupancy track (bigWig); for dyad calls from paired-end BAM use nucleosome-positioning |
| `--chrom` | required |  | contig name as written in the bigWig, e.g. chr1 |
| `--start` | optional |  | 0-based, inclusive start (default 0) |
| `--end` | required |  | 0-based, exclusive end of the region to summarise |
| `--bin-size` | optional | `10` | bin width in bp (default 10) |
| `--out` | required |  | output directory for nucleosome results |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
