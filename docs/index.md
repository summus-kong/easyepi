# easyepi

Bulk epigenomics downstream analysis in pure Python. Every tool takes an aligned
**BAM** and writes ordinary files — BigWig, BED, TSV, JSON, PNG. Peak calling imports
MACS3, so that package must be installed; no other external program is
needed, and no step depends on another having been run first.

```bash
pip install "easyepi[bam,viz]"
easyepi bam2bw sample.bam -o tracks/ --assay atac
easyepi peaks sample.bam -o peaks/ --assay atac
```

## What it covers

| Assay | Status |
| --- | --- |
| [ATAC-seq](assays/atac.md) | full — Tn5 shift, FRiP, peaks, footprint, nucleosome |
| [ChIP-seq](assays/chip.md) | full — narrow/broad peaks, FRiP, chromatin states |
| [CUT&Tag / CUT&RUN](assays/cuttag.md) | full — same Tn5 path as ATAC |
| [MNase-seq](assays/mnase.md) | partial — occupancy and periodicity; dyads need paired-end |
| [WGBS / bisulfite](assays/wgbs.md) | full — per-cytosine calls from the BAM, two-group DMR |

## Design commitments

**BAM is the only required input format.** BigWig, peaks and methylation calls are
things this suite produces, so a study that starts from alignments never has to
assemble a second toolkit.

**One tool, one job.** 28 subcommands, each reading files and writing files. Nothing
is hidden behind a "run everything" wrapper, and a partial failure is reported per
step.

**Limits are printed, not buried.** A motif too short to reach the requested FPR is
listed as skipped rather than shown as "no site". A FRiP computed over a regional
peak set is marked `interpretable: false`. A methylation test that is binomial rather
than beta-binomial says so in its own output.

## Where to go next

- [Install](install.md) — dependencies, optional extras, external programs
- [Quick start](quickstart.md) — one command per analysis step
- [End-to-end examples](EXAMPLES.md) — a complete sequence per assay
- [Commands](commands/index.md) — every subcommand
- [API reference](api/bam2bigwig.md) — the Python interface
- [Reliability and limitations](RELIABILITY.md) — read before interpreting results

## Status

Release candidate. The test suite establishes software behaviour, not biological
validity: `pytest tests -q` passes on Linux, `maintenance/verify_linux.py` checks the
CLI surface, and `scripts/verify_macs3_concordance.py` proves the MACS3 wrapper
produces byte-identical peaks to a direct `macs3 callpeak`. The legacy pure-Python
helpers exported from `easyepi` have no equivalent external concordance study and are
not reachable from `easyepi peaks`.

License: BSD-3-Clause.
