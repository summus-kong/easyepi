# Motif and regulator methods

| Task | Recommended entry | Dependency | Interpretation |
| --- | --- | --- | --- |
| Single best PWM match | `universal.motif_scan` | NumPy | Best sequence score and orientation |
| Known-motif sequence-set enrichment | `universal.motif_enrichment` | Biopython, SciPy | Thresholded sequence presence; Fisher exact + BH |
| Complete calibrated site calls | `universal.motif_analysis --backend fimo` | MEME Suite FIMO | Site table retained as fimo.tsv |
| Established de novo discovery | `universal.motif_analysis --backend meme` | MEME Suite | Official MEME output |
| Optional HOMER analysis | `universal.motif_analysis --backend homer` | HOMER/Perl | Official HOMER output, no fallback |
| Regulator activity | `universal.tf_activity` | decoupler 1.9 | ULM t-statistics from measured matrix and declared network |
| Footprint score wrapper | Python `external_tools.run_tobias_footprint` | TOBIAS | One bias-corrected track per call; not automatic bias correction or BINDetect |

The PFM parser accepts one plain or labelled JASPAR matrix. Uniform-background PWM thresholds from Biopython are model-specific; the per-strand false positive rate is not a genome-wide FDR. Matching sequence length, GC, accessibility, mappability and the tested universe remains the user's responsibility. Genomic annotation and motif database labels answer different questions.

ULM is the published decoupler method. A motif-peak network may estimate motif-associated accessibility, but it does not implement chromVAR background deviations or establish direct binding. Supply measured normalized features and a biologically justified network. TFs with indistinguishable motifs cannot be resolved using motif identity alone.

Native de novo helpers and HOMER-style exploratory utilities remain as low-level compatibility functions. They have not been validated for equivalence to HOMER/MEME; no exact-replication claim is made. `auto` in the historical HOMER API now requires HOMER and never silently changes methods.

See [usage](USAGE.md), [reliability](RELIABILITY.md) and [primary references](REFERENCES.md).
