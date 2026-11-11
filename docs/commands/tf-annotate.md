# easyepi tf-annotate

```bash
easyepi tf-annotate --peaks p.narrowPeak --genome-fasta g.fa --motifs jaspar.pfm --out motifs/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--peaks` | required |  | BED/narrowPeak/broadPeak to annotate |
| `--genome-fasta` | required |  | indexed genome FASTA (.fai required); peak sequences are fetched from it |
| `--motifs` | required |  | motif database: JASPAR .pfm/.jaspar, HOMER .motif or MEME .meme; comma-separate to merge several |
| `--fpr` | optional | `0.0001` | per-motif per-sequence false-positive rate used for the exact threshold (default 1e-4) |
| `--background` | optional |  | background-sequences FASTA; enables the Fisher enrichment test and BH correction |
| `--top` | optional | `5` | hits reported per peak (default 5) |
| `--out` | required |  | output directory for TF annotation tables |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
