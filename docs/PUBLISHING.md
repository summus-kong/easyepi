# Local release verification

This project has not been uploaded or published by the refactor. Build/install checks are local verification only; the user did not request a distributable archive.

## Before release

1. Resolve the historical cleanup limitation documented in RELIABILITY.md. Review `maintenance/clean_legacy.ps1` in preview mode first.
2. Check author/project URLs, version and citation metadata. Unverified repository/documentation URLs were removed; add confirmed public URLs only.
3. Install the intended extras in a clean Python >=3.11 environment and run the full test suite.
4. Run each independent module's help; import from outside the source checkout to detect shadowing.
   `scripts/verify_installed.py` does this against an installed wheel: it refuses to run
   from inside the checkout (where the answer is meaningless), then checks every
   subcommand's `--help`, that no legacy module shipped, and that all 135 modules import.
   Point it at any interpreter that has the wheel installed and run it from another
   directory, for example `cd /tmp && <venv>/bin/python <repo>/scripts/verify_installed.py`.
5. Build from a clean source snapshot and inspect wheel/sdist contents: no R files, legacy python nesting, cached analyses, datasets, private paths or old package modules.
6. Run real MACS3/ChromHMM/MEME and BAM/BigWig examples on a supported Linux environment, recording tool/reference versions and expected results.
7. Complete biological and performance benchmarks before claiming method equivalence, runtime improvements or publication readiness.

A standard local build is `python -m build`, after cleaning old build/egg metadata. Old build directories have been removed from the working tree, so `python -m build` can be run directly from a clean checkout. See VERIFICATION.md for actual checks; do not reuse stale build directories as release evidence.

## Dependency policy

The base package uses NumPy, SciPy, pandas and Matplotlib. Specialist extras are activity, bam, bigwig, deseq2, dev, motif, states, stats, tobias, viz, and `all` (which unions them). `enrich` and `plot` were listed here historically but never existed in `pyproject.toml`; they are gone from this list rather than silently mapping to something else. Tool discovery does not require an activated environment: executable lookup
searches `PATH`, then the interpreter's own prefix, then any explicit
`search_paths`. HOMER's helper programs are on PATH for the child process even
when the parent shell never activated the environment, and a run in which HOMER
reports a missing helper raises rather than reporting an empty result.

External executables and databases are not bundled; `pyproject.toml` declares
them under `[tool.easyepi.external]` so the requirement is readable rather than
inferred from a failure. The analyses themselves are not all pure Python: peak
calling imports `macs3`, `motifs` delegates de novo discovery to HOMER, and
`chromhmm` drives ChromHMM's jar.

No runtime code path launches R. `easyepi/dmr_r.py` does launch it, but it is a
test bridge only: its single importer is `tests/test_beta_binomial_matches_dss.py`,
which uses it to check the beta-binomial port against the real DSS package.
`tools.dmr_r.legacy_dss_bridge` records this so the distinction is machine
readable, and `tests/test_release_metadata.py` asserts both the classification and
that the module carries no hardcoded workstation path. The activity adapter pins decoupler to the verified 1.9 API; review and test API changes before broadening that requirement.

Do not publish a wheel from a dirty prior build directory. This document does not grant permission to upload to PyPI or create a remote release.

## Linux validation target

Create `environment.yml` or install the CI dependencies, then run `python maintenance/verify_linux.py`. The verifier writes command logs and dependency versions under `.verification/linux`. The GitHub Actions workflow runs on Ubuntu with Python 3.11 and 3.12. It is supplied as configuration and has not been executed remotely during this session. Native Windows support is not a release target.
