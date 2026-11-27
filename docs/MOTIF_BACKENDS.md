# Motif and regulator methods

| Task | Recommended entry | Dependency | Interpretation |
| --- | --- | --- | --- |
| Matched control regions | `easyepi background` | NumPy, SciPy | GC- and length-matched controls plus a report of how well they matched. Accessibility is **not** matched |
| Single best PWM match | `easyepi.motifs.scan_sequences`, `universal.motif_scan` | NumPy | Best sequence score and orientation |
| All sites above a threshold | `easyepi.motifs.scan_sites` | NumPy | Per-site coordinates, strand, score and p-value, as offsets into the supplied sequence |
| Known-motif sequence-set enrichment | `easyepi.motifs.homer_workflow.known_motif_enrichment_homer`, `universal.motif_enrichment` | NumPy, SciPy (Biopython optional) | Thresholded sequence presence; binomial or Fisher exact + BH. The sequence is the unit, and the per-motif threshold is not a cross-motif FDR |
| De novo discovery | `easyepi motifs --backend homer` | HOMER | HOMER's own output. There is no built-in engine and no fallback |
| De novo via MEME or FIMO | `universal.motif_analysis --backend meme/fimo/homer` | MEME Suite or HOMER | Official tool output, written under the run directory |
| Regulator activity | `universal.tf_activity` | decoupler 1.9 | ULM t-statistics from a measured matrix and a declared network. Not chromVAR |
| Footprint score wrapper | `external_tools.run_tobias_footprint` | TOBIAS | One bias-corrected track per call; not automatic bias correction or BINDetect |

The PFM parser accepts one plain or labelled JASPAR matrix. Uniform-background PWM thresholds from Biopython are model-specific; the per-strand false positive rate is not a genome-wide FDR. Matching sequence length, GC, accessibility, mappability and the tested universe remains the user's responsibility. Genomic annotation and motif database labels answer different questions.

ULM is the published decoupler method. A motif-peak network may estimate motif-associated accessibility, but it does not implement chromVAR background deviations or establish direct binding. Supply measured normalized features and a biologically justified network. TFs with indistinguishable motifs cannot be resolved using motif identity alone.

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
