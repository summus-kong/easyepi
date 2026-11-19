# Scope and evidence for future evaluation

The current contribution is an independent-module Python interface with shared sample metadata, provenance, explicit optional dependencies and reusable plotting. It combines established methods rather than proposing a new statistical method.

## What has actually been measured

A manuscript would need evidence beyond interface organization. The evidence that now exists, and its exact limits:

| Claim | Evidence | Limit |
| --- | --- | --- |
| The MACS3 wrapper does not alter parameters | Byte-identical peak files vs a direct `macs3 callpeak` in 3 configurations (BAM narrow, bedGraph narrow, bedGraph broad), `scripts/verify_macs3_concordance.py` | Says nothing about biological correctness, and nothing about the legacy pure-Python helpers |
| Known-motif annotation behaves sensibly on real data | 375/622 ENCODE H3K4me3 peaks carry a JASPAR site at FPR 1e-4; top factors SP1 28.5%, CTCF 19.0%, NRF1 17.4% | Not compared against FIMO/HOMER output |
| Peak-set similarity is width-robust | Base-pair Jaccard reproduces the expected mark ordering on ENCODE data (H3K4me3/H3K27ac highest, H3K27me3 lowest) | Single 5 Mb region, no genome-wide benchmark |
| The Linux runtime is reproducible here | 48 verifier checks, 0 failures, `maintenance/verify_linux.py` | One host; the GitHub Actions workflow has still never executed |

## What remains unsupported

No accuracy comparison against an established caller exists for the pure-Python peak caller, the chromatin-state learner or the nucleosome positioning module. No speed benchmark has been run at genome scale. No demonstration yet shows that multi-track integration improves a specific biological analysis. Historical benchmark images are not validation evidence for this release candidate.

The prioritized development work and limitations are in [RELIABILITY.md](RELIABILITY.md). Method references are in [REFERENCES.md](REFERENCES.md).
