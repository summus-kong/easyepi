# Install

## Requirements

- Python 3.9 or newer (tested on 3.11)
- Linux, macOS or Windows

## Core install

```bash
pip install easyepi
```

This installs `numpy`, `scipy` and `pandas`, which every analysis command needs.

## With BAM and BigWig support

```bash
pip install "easyepi[bam]"
```

Adds `pysam` (reading alignments) and `pyBigWig` (reading and writing signal). Almost
every command needs these two; install them unless you only intend to use the
table-level tools.

## With figures

```bash
pip install "easyepi[viz]"
```

Adds `matplotlib`. Only needed when `--figure-format` is used.

## Everything

```bash
pip install "easyepi[bam,viz,states]"
```

`states` adds `hmmlearn` for `easyepi histone-states`.

## Check the installation

```bash
easyepi --version
python -c "import easyepi, pysam, pyBigWig; print('ok')"
```

If an optional dependency is missing, the first command that needs it fails with the
exact install line for the interpreter that is running — not a bare
`ModuleNotFoundError`.

## External programs (optional)

| Program | Used by | Without it |
| --- | --- | --- |
| `macs3` | `easyepi peaks` | falls back to the built-in Poisson caller |
| `samtools` | convenience only | not required |
| `bedtools` | convenience only | not required |

`bam2bw`, `qc`, `dmr`, `tf-annotate`, `annotate`, `footprint` and the figure commands
need no external program at all.

## From source

```bash
git clone https://github.com/summus-kong/easyepi
cd easyepi
pip install -e ".[bam,viz]"
pytest tests -q
```

## Reference files

Several commands need a genome annotation or a reference sequence. They are not
bundled; point the flags at your own copies.

| File | Flags | Needed by |
| --- | --- | --- |
| Genome FASTA with a `.fai` index | `--genome-fasta`, `references.genome` | `tf-annotate`, `motifs` |
| GTF/GFF annotation (gzip is fine) | `--gtf`, `references.annotation` | `annotate`, `gene-panel` |
| Motif database (JASPAR `.pfm`, HOMER `.motif`, MEME `.meme`) | `--motifs` | `tf-annotate`, `motifs` |

A `.fai` can be made with `samtools faidx genome.fa`.
