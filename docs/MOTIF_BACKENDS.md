# Motif and regulator methods

| Task | Recommended entry | Dependency | Interpretation |
| --- | --- | --- | --- |
| Matched control regions | `easyepi background` | NumPy, SciPy | GC- and length-matched controls plus a report of how well they matched. Accessibility is **not** matched |
| Single best PWM match | `easyepi.motifs.scan_sequences`, `universal.motif_scan` | NumPy | Best sequence score and orientation |
| All sites above a threshold | `easyepi.motifs.scan_sites` | NumPy | Per-site coordinates, strand, score and p-value, as offsets into the supplied sequence |
| Known-motif sequence-set enrichment | `easyepi.motifs.ame_enrichment.ame_enrichment` (what `easyepi motifs --db` calls) | MEME Suite (`ame`); non-MEME databases converted automatically | AME partitions sequences by affinity and reports its own p, adj-p and E-value. The test (`--ame-method`, 6 choices) and scoring (`--ame-scoring`, 4 choices) are recorded in the result because they do not agree with each other |
| De novo discovery | `easyepi motifs --backend homer` | HOMER | HOMER's own output. There is no built-in engine and no fallback |
| De novo via MEME or FIMO | `universal.motif_analysis --backend meme/fimo/homer` | MEME Suite or HOMER | Official tool output, written under the run directory |
| Regulator activity | `universal.tf_activity` | decoupler 1.9 | ULM t-statistics from a measured matrix and a declared network. Not chromVAR |
| Footprint score wrapper | `external_tools.run_tobias_footprint` | TOBIAS | One bias-corrected track per call; not automatic bias correction or BINDetect |

The PFM parser accepts one plain or labelled JASPAR matrix. Uniform-background PWM thresholds from Biopython are model-specific; the per-strand false positive rate is not a genome-wide FDR. Matching sequence length, GC, accessibility, mappability and the tested universe remains the user's responsibility. Genomic annotation and motif database labels answer different questions.

ULM is the published decoupler method. A motif-peak network may estimate motif-associated accessibility, but it does not implement chromVAR background deviations or establish direct binding. Supply measured normalized features and a biologically justified network. TFs with indistinguishable motifs cannot be resolved using motif identity alone.

### Why enrichment runs through AME

Four implementations were measured on the same 300 peaks and 300 background
sequences, and against each other:

| implementation | wall time | rank correlation with HOMER |
| --- | --- | --- |
| local scan + binomial | 77 s | rho(hit rate) 0.109, rho(log P) 0.237 |
| port of HOMER's own hypergeometric and threshold search | 88 s | 0.120 / 0.206 |
| FIMO scan + the same hypergeometric | 89 s | 0.169 / 0.233 |
| **AME 5.5.9 (fisher / avg)** | **5 s** | 0.396 (n=144) |

FIMO and the local scan, compared directly to each other, reach only rho 0.292 on
hit rates and 0.435 on log P. Two independent standard scanners therefore disagree
with each other about as much as either disagrees with HOMER: the ranking is set by
scanner and partition choices, not by which implementation is "correct". Chasing
numerical agreement with one tool was the wrong objective, so the method was
delegated to the published tool for the question instead of reimplemented.

What this does **not** claim: results are not numerically identical to HOMER. AME
partitions sequences by affinity score, HOMER searches a threshold. They answer
overlapping but different questions, which is why the test and scoring are written
into the output rather than left implicit.

Native de novo helpers were **removed** on 2026-11-27; they are not kept as
compatibility functions. `easyepi/motifs/homer_like.py`, `de_novo_homer` and
`discover_motifs` were measured against HOMER on real peak sequences and did not
reproduce it: the k-mer counter shared a mean 3.90 bp of consensus with HOMER's top
ten, and the EM variant returned `AAAAAAAA` and `CCAGCCTG` as its best two motifs
while reporting all 600 sequences as sites.

`easyepi motifs --backend auto` therefore either runs HOMER or reports that HOMER is
required. `--backend homer` raises when HOMER is absent, and a non-zero HOMER exit
raises as well — including the specific case where HOMER exits 0 with a 0-byte motif
file because its helper programs were missing. Known-motif enrichment and PWM
scanning remain and need no external program, so HOMER is a requirement of one
analysis, not of the package.

See [usage](USAGE.md), [reliability](RELIABILITY.md) and [primary references](REFERENCES.md).
