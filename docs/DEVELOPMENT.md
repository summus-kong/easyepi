# Development guide and source map

## Design decisions

The implementation follows one independently callable function module per analysis task, with universal functions separated from assay-specific ones. The optional project manifest describes stable sample/reference inputs; per-run options remain explicit CLI arguments. No fixed standard/discovery pipeline is provided.

The application is developed in English and Python for Linux. Native Windows compatibility is outside the support scope. Optional external backends retain their original implementations. New statistical backends must call established, cited methods and declare their input assumptions, units and experimental-design restrictions. Preserve underlying legacy functions for compatibility, but label exploratory behavior and reject unsupported apparent analyses.

## Refactoring sequence

1. Audited the actual live files separately from old build outputs; old cached package contents were not treated as current source.
2. Moved universal/specific trees out of the redundant python directory; moved dispatch, runtime and plotting helpers to package root.
3. Updated imports/tests to the new paths and constrained package discovery to maintained packages.
4. Examined assay adapters, external commands, input parsing, missingness, pairing and plotting behavior.
5. Connected MACS3, corrected ChromHMM orchestration, added decoupler ULM and Biopython motif thresholds, and replaced hand-written k-means with scikit-learn.
6. Corrected data-integrity errors and disabled placeholder or unrelated implementations that could produce misleading success.
7. Reworked dependency extras, English documentation, release manifests and regression tests.
8. Ran source/import checks, tests, a real decoupler smoke example and isolated installation-content verification. See VERIFICATION.md for recorded results and unresolved limitations.

## Maintained directory layout

```text
easyepiv1/
  easyepi/                  Importable Python package
    universal/              Independent shared analysis scripts
    specific/
      atac/                 ATAC QC, cuts, peaks, consensus, matrices, comparison
      chip/                 ChIP/CUT&Tag/CUT&RUN entries
      histone/              Histone QC and official ChromHMM stages
      methylation/          Methylation QC, profile, annotation, site comparison
      mnase/                MNase summaries; positioning explicitly unsupported
    core/                   Models, readers and statistical primitives
    annotate/               Genomic annotation helpers
    differential/           SciPy comparisons and design matching
    motifs/                 Motif formats, scanning and optional external tools
    normalize/              Explicit matrix normalization
    nucleosome/             Descriptive occupancy/periodicity helpers
    peaks/                  Legacy exploratory low-level utilities
    qc/                     BAM QC interface
    viz/                    Compatibility plotting wrappers
    _runtime.py             Project-aware module setup
    _dispatch.py            Shared CLI dispatch
    _plots.py               Shared plotting implementations
    _peak_calling.py         Official MACS3 command adapter
    activity.py             decoupler ULM adapter
    plot.py                 Public plotting exports
  docs/                     English user/developer/reliability documentation
  examples/                 Manifest and Python API examples
  tests/                    API and data-integrity regression tests
  environment.yml           Linux analysis environment specification
  .github/workflows/        Linux CI configuration
  maintenance/              Linux verifier and host cleanup helper
  pyproject.toml            Build, dependencies and entry-point metadata
  MANIFEST.in               Source-distribution inclusion/exclusion policy
  LICENSE                   BSD-3-Clause license
  CITATION.cff              Software citation metadata
  CHANGELOG.md              Release-candidate changes
  README.md                 Entry point to documentation
```

Historical build, cache and output directories have been removed from the working tree with `maintenance/clean_legacy.ps1 -Apply`. They were never part of the maintained layout and are excluded from distributions by `pyproject.toml` and `MANIFEST.in`. The sibling project outside easyepiv1 was intentionally not modified. The only retained local-only tree is `.verification`, which holds project-local validation dependencies and results; it is ignored and excluded from distributions.

## File responsibilities

| File | Responsibility |
| --- | --- |
| `easyepi/__init__.py` | Independent Python tools for bulk epigenomics data and signal summaries. |
| `easyepi/_dispatch.py` | Dispatch independent script entry points to shared CLI functions. |
| `easyepi/_peak_calling.py` | Assay-aware command construction for the official MACS3 executable. |
| `easyepi/_plots.py` | Shared plotting functions and optional module-output rendering. Plotting reads completed analysis tables without altering statistical results. |
| `easyepi/_runtime.py` | Shared project parsing and run initialization for independent scripts. |
| `easyepi/activity.py` | Regulator activity estimates using the published decoupler ULM implementation. |
| `easyepi/annotate/__init__.py` | Peak annotation and functional enrichment. |
| `easyepi/assay_commands.py` | Independent assay input QC using the shared project and logging contract. These summaries do not replace dedicated analysis backends. |
| `easyepi/cli.py` | Easyepi command line interface. |
| `easyepi/core/__init__.py` | Core package: data model, IO and statistics primitives. |
| `easyepi/core/io.py` | Genomic IO: BAM (pysam), BigWig/BedGraph (pyBigWig), BED. |
| `easyepi/core/models.py` | Core data model: unified in-memory representation for epigenomics downstream analysis. |
| `easyepi/core/stats.py` | Core statistics primitives: descriptive stats, significance tests, correction. |
| `easyepi/differential/__init__.py` | Differential analysis and statistical modeling (native). |
| `easyepi/differential_matrix.py` | Descriptive differential analysis of integrated BigWig signal. |
| `easyepi/footprint.py` | Aggregate motif-site signal profiles for descriptive visualization. |
| `easyepi/integration.py` | Gene-centric integration of normalized BigWig tracks, peaks and RNA result tables. Nearest-TSS links are candidates, not causal regulation. |
| `easyepi/logging_utils.py` | Shared console panels, optional detailed logs, status records and command provenance. |
| `easyepi/methylation.py` | Strict zero-based methylation tables, coverage summaries and exploratory two-group site comparisons. This module does not call DMRs. |
| `easyepi/motifs/__init__.py` | PWM parsing and scanning, optional external motif tools and legacy utilities. |
| `easyepi/motifs/external_tools.py` | Optional external motif-analysis tool backends. |
| `easyepi/motifs/formats.py` | Motif format layer (L2): HOMER / JASPAR / MEME I/Os and cross-conversion. |
| `easyepi/motifs/homer_like.py` | Legacy exploratory motif utilities retained for compatibility. |
| `easyepi/motifs/homer_workflow.py` | Optional official HOMER delegation and legacy exploratory motif utilities. |
| `easyepi/motifs/logo.py` | Motif figure helpers: sequence logo + enrichment bar chart. |
| `easyepi/normalize/__init__.py` | Signal / count normalization: library size, spike-in, RPGC, CPM. |
| `easyepi/nucleosome/__init__.py` | Nucleosome occupancy / positioning analysis (native, dependency-free). |
| `easyepi/peaks/__init__.py` | Peak calling / event identification (native). |
| `easyepi/plot.py` | Shared public plotting API. Functions return Matplotlib Figures; callers choose whether and how to export them. |
| `easyepi/project.py` | Project and sample input contracts. Explicit CLI reference overrides take precedence over project values; resolved inputs are recorded for each run. |
| `easyepi/project_commands.py` | Independent project validation and QC-only reporting. |
| `easyepi/pubfig.py` | Easyepi publication figure helpers (Nature-style, Python/matplotlib). |
| `easyepi/qc/__init__.py` | QC and statistical summary module. |
| `easyepi/specific/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/specific/atac/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/specific/atac/consensus.py` | Construct consensus intervals across samples. |
| `easyepi/specific/atac/count.py` | Count filtered alignment overlaps or integrate BigWig signal over fixed regions. |
| `easyepi/specific/atac/differential.py` | Exploratory replicated two-group comparison. |
| `easyepi/specific/atac/peaks.py` | ATAC narrow peak calling with MACS3. |
| `easyepi/specific/atac/qc.py` | ATAC-seq QC。 |
| `easyepi/specific/atac/tn5_shift.py` | Independent tn5 shift analysis entry. |
| `easyepi/specific/chip/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/specific/chip/broad_peaks.py` | CHIP broad peak calling with MACS3. |
| `easyepi/specific/chip/consensus.py` | Construct consensus intervals across samples. |
| `easyepi/specific/chip/count.py` | Independent count analysis entry. |
| `easyepi/specific/chip/differential.py` | Exploratory replicated two-group comparison. |
| `easyepi/specific/chip/narrow_peaks.py` | CHIP narrow peak calling with MACS3. |
| `easyepi/specific/chip/qc.py` | ChIP-seq/CUT&Tag/CUT&RUN QC。 |
| `easyepi/specific/histone/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/specific/histone/chromhmm.py` | Prepare inputs and run the official ChromHMM BinarizeBed/LearnModel commands. |
| `easyepi/specific/histone/chromhmm_annotate.py` | Summarize ChromHMM state coverage and optional user-reviewed labels. |
| `easyepi/specific/histone/chromhmm_binarize.py` | Standalone ChromHMM binarize entry. |
| `easyepi/specific/histone/chromhmm_prepare.py` | Prepare headerless ChromHMM CellMark inputs. |
| `easyepi/specific/histone/chromhmm_train.py` | Standalone ChromHMM train entry. |
| `easyepi/specific/histone/qc.py` | Independent assay input QC. |
| `easyepi/specific/methylation/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/specific/methylation/annotate.py` | Annotate methylation intervals with candidate nearest TSS links. |
| `easyepi/specific/methylation/differential.py` | Exploratory replicated two-group comparison. |
| `easyepi/specific/methylation/profile.py` | Assay-specific signal or methylation profile summaries. |
| `easyepi/specific/methylation/qc.py` | Independent assay input QC. |
| `easyepi/specific/mnase/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/specific/mnase/nucleosome.py` | MNase signal summaries; no de novo positioning is inferred. |
| `easyepi/specific/mnase/positioning.py` | Reserved entry for an established MNase nucleosome-positioning backend. |
| `easyepi/specific/mnase/profile.py` | Assay-specific signal or methylation profile summaries. |
| `easyepi/specific/mnase/qc.py` | MNase-seq QC。 |
| `easyepi/universal/__init__.py` | Independent analysis scripts grouped by assay or shared scope. |
| `easyepi/universal/bam_qc.py` | Independent BAM input QC. |
| `easyepi/universal/make_qc_report.py` | Summarize input QC without generating downstream biological conclusions. |
| `easyepi/universal/motif_analysis.py` | Run established external motif discovery or scanning without heuristic fallback. |
| `easyepi/universal/motif_enrichment.py` | Known-motif enrichment using Biopython PWM thresholds, Fisher exact tests and BH correction. |
| `easyepi/universal/motif_scan.py` | Scan known JASPAR/PFM motifs and return best sequence scores. |
| `easyepi/universal/multiomics_integrate.py` | Integrate normalized epigenomic tracks and RNA result tables by gene. |
| `easyepi/universal/normalize_signal.py` | Normalize a measured matrix with explicit factors and run provenance. |
| `easyepi/universal/peak_annotation.py` | Annotate peak genomic features and distance to transcription start sites. |
| `easyepi/universal/rna_prepare.py` | Prepare RNA result identifiers with explicit logging and no implicit model fitting. |
| `easyepi/universal/tf_activity.py` | Estimate TF activity from measured normalized features using decoupler ULM. |
| `easyepi/universal/tf_activity_differential.py` | Compare externally supported TF activity estimates between two replicated groups. |
| `easyepi/universal/validate_project.py` | Validate project inputs and record resolved references. |
| `easyepi/viz/__init__.py` | Visualization (native matplotlib/seaborn). |
| `easyepi/washu.py` | WashU Epigenome Browser tracks visualization (optional backend). |

Additional files: `examples/project.toml` illustrates sample/reference metadata; `examples/demo_import.py` demonstrates import use. `tests/test_release_contract.py` protects corrected data/method boundaries; other tests retain existing API and motif/integration coverage. Historical PNG files under docs are not cited as release validation and are excluded by the explicit source manifest.

## Add a feature

Expose an importable function first; keep command parsing in the independent module. Use RunRecorder/setup_logging for status and optional detailed logs. Add no implicit group split, control, coverage, library size or statistical fallback. Return Figures from plotting helpers and export only on explicit request. For mature external tools, check nonzero exits and preserve outputs; record the exact command. Add a targeted test for a real failure mode rather than mirroring the implementation.

Use unique output directories per module. Shared files such as status.json are not a workflow scheduler or concurrent-run registry. Threads are passed only to supported backends; the generic --threads metadata flag does not parallelize every Python function. Paths in input data may contain Unicode, but maintained source prose and CLI output are English.

## Where commits are made and pushed

The analysis host is the development machine: it runs the tests, and its interpreter
is the one `docs/VERIFICATION.md` records. It is not where a push happens.

**The local workstation is the only portal to GitHub.** Commits are made and
`git push` runs on the local workstation, never on the analysis host. A change
developed and verified remotely therefore travels back first:

```bash
# 1. on the analysis host: verify, then copy the changed files to a staging path
python -m pytest tests -q
tar -C "$STAGE" -czf /tmp/out.tar.gz .

# 2. pull the artifacts to the local workstation, then unpack into the local tree
#    and confirm byte identity before committing
md5sum <files>            # compared on both sides

# 3. on the local workstation only
git add -A && git commit && git push origin main
```

Two properties this ordering protects:

* **The pushed tree is exactly the verified tree.** Files are compared by
  checksum after the round trip; a silent edit on either side shows up as a
  mismatch rather than as a commit that was never tested.
* **Credentials and history stay in one place.** The analysis host holds no push
  credentials, so a compromised or misconfigured shared host cannot publish.

The analysis host can reach `github.com` (a `git ls-remote` succeeds from it), so
this is a working rule, not a network limitation. Do not substitute a remote
`git push` for the round trip just because it happens to work.

## Linux verification commands

```bash
python -m pip install -e ".[dev,stats,motif]"
python -m pytest -q
python -m easyepi.cli --help
python -m easyepi.specific.atac.peaks --help
python maintenance/verify_linux.py
```

External backends need separate integration checks on real data. Keep environment and method validation separate; do not interpret a green unit suite as proof of scientific validity. See RELIABILITY.md and REFERENCES.md.

Linux verification is implemented in `maintenance/verify_linux.py` and `.github/workflows/linux.yml`. Optional-reader tests use synthetic real BAM/BigWig files when their extras are installed; ULM tests call the real decoupler implementation. The available host has no running Linux environment, so CI execution is pending. The PowerShell cleanup helper is for removing pre-existing artifacts on the current development host only; it is not part of the Linux analysis runtime.
