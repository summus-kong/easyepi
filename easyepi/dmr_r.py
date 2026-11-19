"""Bridge to DSS through Rscript, with no silent fallback.

Why a subprocess

DSS is an R package with no Python equivalent, so the test runs in R. The bridge
deliberately fails loudly when R or DSS is missing: a methylation test that
silently degrades to a less valid method is worse than one that refuses, because
the caller cannot then tell which model produced the numbers.

API notes (verified against DSS 2.58.0)

* ``DMLtest(BSobj, group1, group2)`` is the two-group entry point. There is no
  ``DMLfit`` export, and ``DMLtest``'s arguments are ``group1``/``group2``, not
  ``sample1``/``sample2``.
* ``BSseq()`` takes ``M`` and ``Cov`` as loci x samples, so the samples x loci
  matrices from Python must be transposed.
* The returned frame carries ``diff``, ``pval`` and ``fdr``.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np

#: Rscript locations tried before falling back to PATH.
_CANDIDATE_ENVS = (
    "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/seurat5",
    "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi",
)

_R_SCRIPT = r"""
data_path <- commandArgs(trailingOnly = TRUE)[1]
con <- file(data_path, "r")
on.exit(close(con), add = TRUE)

n_sites   <- as.integer(readLines(con, n = 1))
n_samples <- as.integer(readLines(con, n = 1))
meth  <- matrix(scan(con, what = integer(), n = n_sites * n_samples, quiet = TRUE),
                nrow = n_samples, byrow = TRUE)
total <- matrix(scan(con, what = integer(), n = n_sites * n_samples, quiet = TRUE),
                nrow = n_samples, byrow = TRUE)
design <- scan(con, what = integer(), n = n_samples, quiet = TRUE)

suppressPackageStartupMessages({
  library(DSS)
  library(bsseq)
  library(GenomicRanges)
  library(IRanges)
})

loci <- GRanges(seqnames = "chr1",
                ranges = IRanges(start = seq_len(n_sites), width = 1))
sample_names <- paste0("s", seq_len(n_samples))
# M and Cov must be loci x samples; Python sends samples x loci.
bs <- BSseq(M = t(meth), Cov = t(total), gr = loci,
            sampleNames = sample_names, rmZeroCov = FALSE)

group1 <- sample_names[design == 0]
group2 <- sample_names[design == 1]

# DSS parallelises with mcfork(); on a loaded or memory-constrained host that
# fails with "unable to fork, Cannot allocate memory". One core is plenty for a
# per-site test and removes the failure mode entirely.
fit <- tryCatch(
  DMLtest(bs, group1 = group1, group2 = group2, smoothing = TRUE, ncores = 1),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (is.null(fit)) quit(status = 1)

# Report the dispersion actually used, so a run can be traced to a model.
cat("DISPERSION\t", mean(c(fit$phi1, fit$phi2), na.rm = TRUE), "\n", sep = "")
out <- data.frame(diff = fit$diff, pval = fit$pval)
write.table(out, stdout(), sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
"""


@lru_cache(maxsize=8)
def find_rscript() -> str | None:
    """Locate an Rscript that can load DSS, or None.

    Cached: probing imports R and DSS, which costs seconds, and callers ask
    repeatedly.
    """
    candidates: list[str] = []
    for env in _CANDIDATE_ENVS:
        cand = Path(env) / "bin" / "Rscript"
        if cand.exists():
            candidates.append(str(cand))
    on_path = shutil.which("Rscript")
    if on_path:
        candidates.append(on_path)

    for cand in candidates:
        try:
            probe = subprocess.run(
                [cand, "--vanilla", "-e",
                 'quit(status = if (requireNamespace("DSS", quietly = TRUE)) 0 else 1)'],
                capture_output=True, timeout=180,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        if probe.returncode == 0:
            return cand
    return None


def run_dss_test(meth: np.ndarray, total: np.ndarray, design: np.ndarray,
                 timeout: int = 3600) -> dict:
    """Run DSS::DMLtest on (n_samples, n_sites) counts.

    ``design`` is 0 for group A and 1 for group B. Raises RuntimeError when R or
    DSS is unavailable -- it never falls back to a different model.
    """
    rscript = find_rscript()
    if rscript is None:
        raise RuntimeError(
            "DSS is not available: no Rscript with the DSS package was found. "
            "Install with BiocManager::install('DSS'), or run the pooled test "
            "explicitly and report that between-replicate variance is unmodelled.")

    meth = np.asarray(meth, dtype=np.int64)
    total = np.asarray(total, dtype=np.int64)
    if meth.ndim != 2:
        raise ValueError("meth must be 2-D (n_samples, n_sites)")
    n_samples, n_sites = meth.shape
    if total.shape != meth.shape:
        raise ValueError("meth and total must have the same shape")
    if design.shape != (n_samples,):
        raise ValueError(f"design must have shape ({n_samples},)")

    payload = "\n".join([
        str(n_sites),
        str(n_samples),
        " ".join(str(int(v)) for v in meth.ravel()),
        " ".join(str(int(v)) for v in total.ravel()),
        " ".join(str(int(v)) for v in design),
    ]) + "\n"

    with tempfile.TemporaryDirectory(prefix="easyepi-dss-") as tmp:
        rsrc = Path(tmp) / "dss_bridge.R"
        data = Path(tmp) / "counts.txt"
        rsrc.write_text(_R_SCRIPT, encoding="utf-8")
        data.write_text(payload, encoding="utf-8")
        proc = subprocess.run(
            [rscript, "--vanilla", str(rsrc), str(data)],
            capture_output=True, text=True, timeout=timeout,
        )

    if proc.returncode != 0:
        message = (proc.stderr.strip() or proc.stdout.strip())[:2000]
        raise RuntimeError("DSS failed: " + message)

    dispersion = None
    numbers: list[tuple[float, float]] = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("DISPERSION"):
            tail = line.split("\t")[-1].strip()
            dispersion = None if tail in ("NA", "") else float(tail)
            continue
        parts = line.split("\t")
        if len(parts) == 2:
            try:
                numbers.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue

    if len(numbers) != n_sites:
        raise RuntimeError(
            f"DSS returned {len(numbers)} results for {n_sites} sites; "
            "refusing to align them by position")

    arr = np.asarray(numbers, dtype=float)
    return {"diff": arr[:, 0], "pvalue": arr[:, 1], "dispersion": dispersion}
