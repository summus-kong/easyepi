# Python API

Install from the project root with `python -m pip install -e .`. Optional libraries load at feature use, not at package import.

## Input contracts

```python
from easyepi import load_project, validate_project

project = load_project("project.toml")
validate_project(project, check_files=True)
for sample in project.samples:
    print(sample.sample_id, sample.assay, sample.time)
```

## Methylation

```python
from easyepi import read_methylation_table, summarize_methylation

sites = read_methylation_table("methylation.tsv")
summary = summarize_methylation(sites, min_coverage=5)
```

Coordinates are zero-based half-open. Beta must be a fraction, not a percentage. Counts determine coverage; beta-only input has unknown coverage, so filtering then requires `min_coverage=0`.

## Explicit normalization

```python
import numpy as np
from easyepi import cpm_normalize

matrix = np.array([[100, 200], [40, 80]], dtype=float)
normalized = cpm_normalize(matrix, np.array([1_000_000, 2_000_000]))
```

Library sizes must be measured and use the same read/fragment units and filtering as the matrix. RPGC additionally requires positive mean fragment length and coverage input.

## Replicated exploratory comparisons

```python
from easyepi import differential_features

result = differential_features(
    np.array([[1, 2, 3, 8, 9, 10]], dtype=float),
    np.array([True, True, True, False, False, False]),
    method="ttest",  # Welch by default
)
```

For low-level paired Wilcoxon calls, pair order is the caller's responsibility. For design-aware ID matching use `easyepi diff --pair subject`. These tests do not fit a count dispersion model or adjust covariates.

## Motif parsing and scanning

```python
from pathlib import Path
from easyepi import read_pwm, scan_sequences

pwm = read_pwm(Path("TF.pfm"))
hits = scan_sequences(pwm, ["ACGTACGT", "TTTTACGT"], reverse_complement=True)
```

One best score per sequence is returned; positions use forward-sequence zero-based starts, with strand reported separately. No score threshold or binding conclusion is implied. For thresholded sequence-set enrichment use `easyepi.universal.motif_enrichment`, which requires the `motif` extra. Use FIMO for full site tables.

## TF activity from measured features

```python
import pandas as pd
from easyepi import estimate_tf_activity

matrix = pd.read_csv("normalized_accessibility.tsv", sep="\t", index_col=0)
network = pd.read_csv("tf_peak_network.tsv", sep="\t")
activity = estimate_tf_activity(matrix, network, min_targets=5)
```

Network columns are `source`, `target`, `weight`; targets match the matrix row index. Requires `easyepi[activity]`. This calls decoupler ULM and does not reproduce chromVAR or directly predict binding.

## Plot functions

```python
from easyepi.plot import plot_volcano, plot_distribution, plot_heatmap, plot_profile

fig = plot_volcano(pd.DataFrame({"log2fc": result["log2fc"], "padj": result["padj"]}))
fig.savefig("volcano.svg")
profile = plot_profile(np.arange(30))
# Display profile in a notebook, or choose your own output path and format.
```

Shared generic functions accept pandas tables, except `plot_profile`, which accepts a one-dimensional vector. Missing p-values remain missing. `save_figure` is optional and closes the figure after explicit export; call `Figure.savefig` directly when retaining the Figure is useful.

## Independent script calls from Python

```python
from easyepi.specific.atac.qc import main

exit_code = main(["--project", "project.toml", "--out", "results/qc", "--verbose"])
```

All maintained independent entries accept an argv list. Parser errors use `SystemExit`, as in argparse; shared project `--dry-run` also returns through `SystemExit(0)` after writing metadata. Keep distinct output directories per analysis to avoid overwriting shared `status.json` and `metrics.json` files.

## Integration

`read_gene_bed`, `gene_signal_matrix`, `cluster_tracks`, `link_peaks_to_genes`, `read_rna_results` and `run_integration` are importable from `easyepi`. They build descriptive gene-centric tables. Clustering requires scikit-learn; track extraction requires pyBigWig. Linkage is nearest TSS, not causal inference.

## Legacy functions

`call_peaks_from_bam` and `poisson_peaks` remain available for compatibility and inspection; neither is a replacement equivalent to MACS3, and `easyepi peaks` wraps MACS3 instead.

**Motif discovery has no built-in engine.** `discover_motifs`, `homer_like` and `de_novo_homer` were removed in 0.7.1 rather than deprecated. Both were measured against HOMER on real peak sequences and neither reproduced it: the k-mer counter shared a mean 3.90 bp of consensus with HOMER's top ten, which is the overlap of unrelated 8-10 bp motifs, and the EM variant returned `AAAAAAAA` and `CCAGCCTG` as its best two motifs while claiming all 600 sequences as sites. Keeping them meant producing plausible tables that the tool they stood in for would not reproduce.

What does run without any external program is in `easyepi.motifs`: `read_pwm` and `scan_sequences` (numpy sliding-window PWM scanning, compared against FIMO in `docs/dev/TOOL_COMPARISONS.md`), and `known_motif_enrichment_homer` (PWM scan plus a binomial tail and BH correction), which is what `easyepi tf-annotate` and `easyepi motifs --pwm/--db` use. De novo discovery delegates to HOMER through `easyepi.motifs.external_tools.run_homer_motifs`, and a missing HOMER raises `ToolNotAvailableError` naming the install. HOMER is not a dependency of the package.
