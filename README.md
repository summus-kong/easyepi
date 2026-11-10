# Easyepi

A Python toolkit of independent bulk epigenomics analysis modules. Run only the step you need, on one sample or multiple samples. Group labels and controls are optional until a method explicitly requires them. There is no fixed pipeline or standard/discovery mode.

Application code is Python. Optional established tools execute their own published methods: MACS3 for peaks, ChromHMM for chromatin states, MEME/FIMO or HOMER for motif analysis, and decoupler for regulator activity. No R runtime or R backend is used.

## Install

The supported runtime is Linux with Python 3.11 or newer. Use a dedicated environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m pip install -e ".[bam,bigwig,motif,stats,activity,dev]"
```

Install only the extras needed for your analysis. MACS3, Java/ChromHMM, MEME Suite and HOMER are separate optional installations. Development targets Linux only; native Windows support is outside the scope.

A reproducible Linux environment specification is provided in [environment.yml](environment.yml). Linux CI runs the API suite, all script help entries, BAM/BigWig reader integration and ULM smoke checks. CI results are pending until the workflow runs.

## Start with an independent module

```bash
python -m easyepi.universal.validate_project --project examples/project.toml --out results/validation --no-check-files
python -m easyepi.specific.atac.qc --project project.toml --out results/atac_qc --verbose
python -m easyepi.specific.atac.peaks --project project.toml --genome-size hs --out results/peaks
python -m easyepi.universal.motif_enrichment --target peaks.fa --background matched.fa --pwm TF.pfm --out results/motifs
```

The example project references placeholder files; schema-only validation works without downloading them. Replace the paths before analysis. Each module accepts `--help`; Python functions can be imported directly. Existing `easyepi` CLI commands remain independent conveniences. Historical `peaks` and `nucleosome` metadata-only commands now reject execution and point to working modules.

## Documentation

- [Complete usage and input contracts](docs/USAGE.md)
- [QC commands, thresholds and batch roll-up](docs/QC.md)
- [Python API](docs/PYTHON_API.md)
- [Module and file map](docs/DEVELOPMENT.md)
- [Reliability audit and remaining work](docs/RELIABILITY.md)
- [Motif methods and interpretation](docs/MOTIF_BACKENDS.md)
- [Release verification](docs/PUBLISHING.md)
- [Actual verification results](docs/VERIFICATION.md)
- [Method references](docs/REFERENCES.md)

## QC

Assay-specific QC commands share one metric and verdict implementation, so the
same sample gets identical numbers whichever command produced the report:

```bash
easyepi chip-qc    --project project.toml --annotation genes.gtf --out results/chip_qc
easyepi atac-qc    --project project.toml --blacklist blacklist.bed --out results/atac_qc
easyepi mnase-qc   --project project.toml --out results/mnase_qc --figure-format png
easyepi methylation-qc --project project.toml --out results/meth_qc
```

`qc` runs the same analysis for any assay via `--assay`; without it, `qc` keeps
its original BAM-only summary. `make-qc-report` produces a batch-level
`qc_report.json`, and `merge-qc-reports` combines several of those into one
study-wide verdict summary:

```bash
easyepi make-qc-report    --project batch_01/project.toml --out batch_01
easyepi merge-qc-reports  "batch_*/qc_report.json" --out study/qc_merged.json
```

Thresholds are **advisory**: a breach always yields `warn`, never `fail`, and
`overall` is only ever `pass`, `warn` or `not_evaluated`. Override the rule table
per metric with `--qc-thresholds rules.json` or a `[qc_thresholds]` table in the
project TOML (later wins; unspecified metrics keep defaults). Metrics that were
never computed are reported as `not_evaluated`, never as pass.

See [QC commands, thresholds and batch roll-up](docs/QC.md) for the full metric
list, default thresholds, the migration path off `failures`, and the recommended
multi-batch directory layout.

## Plotting

Analysis does not export figures by default. Functions return Matplotlib Figures:

```python
import pandas as pd
from easyepi.plot import plot_volcano

fig = plot_volcano(pd.read_csv("differential.tsv", sep="\t"))
fig.savefig("volcano.pdf", bbox_inches="tight")  # Format and destination are your choice.
```

Modules with optional figure rendering accept repeated `--figure-format` arguments. Exporting a PDF does not also create PNG files.

## Status

This is a release candidate undergoing method and input-contract validation. Unit tests establish software behavior, not biological validity. Exploratory signal comparisons are not count-model differential analysis; methylation site comparisons are not DMR calling; motif-associated accessibility is not direct TF binding. See the reliability audit before interpreting results.

License: BSD-3-Clause. Cite the underlying methods actually used and the software version; there is no assigned Easyepi paper or software DOI.
