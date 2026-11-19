"""Compare `easyepi.dmr_beta_binomial` against DSS, intermediate by intermediate.

This is the harness behind the numbers in `docs/dev/DSS_REPLACEMENT_PLAN.md`. It
is not a pytest module (that is `tests/test_beta_binomial_matches_dss.py`); it
prints the full error table so a mismatch can be localised to one step.

    Rscript --vanilla tests/oracle/dss_oracle.R <counts.txt> <out.tsv>   # driven here
    python tests/oracle/compare_to_dss.py

Needs `Rscript` with DSS for the reference side and nothing else for ours. Run
it from the repository root.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from easyepi.dmr_beta_binomial import beta_binomial_test  # noqa: E402

#: Same fixture as the test: two chromosomes, irregular spacing, scatter inside
#: each group. Any of those three missing makes the comparison vacuous.
PER_CHR = 150
N_A = 3
N_B = 3
N_DIFF = 40

RSCRIPT_CANDIDATES = (
    "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/seurat5/bin/Rscript",
    "Rscript",
)


def fixture():
    rng = np.random.default_rng(11)
    n_sites = 2 * PER_CHR
    chrom = np.array(["chr1"] * PER_CHR + ["chr2"] * PER_CHR)
    gaps = rng.integers(1, 120, n_sites)
    pos = np.concatenate([
        np.cumsum(gaps[:PER_CHR]) + 1000,
        np.cumsum(gaps[PER_CHR:]) + 5000,
    ]).astype(np.int64)
    total = rng.integers(12, 45, (N_A + N_B, n_sites)).astype(np.int64)
    base = np.empty((N_A + N_B, n_sites), dtype=float)
    base[:] = 0.30
    base[N_A:, :N_DIFF] = 0.65
    base[N_A:, PER_CHR:PER_CHR + N_DIFF] = 0.65
    level = np.clip(rng.normal(base, 0.06), 0.01, 0.99)
    meth = rng.binomial(total, level).astype(np.int64)
    design = np.concatenate([np.zeros(N_A, dtype=int), np.ones(N_B, dtype=int)])
    return meth, total, design, chrom, pos


def write_payload(meth, total, design, chrom, pos, path: Path) -> None:
    n_sites, n_samples = meth.shape[1], meth.shape[0]
    path.write_text("\n".join([
        str(n_sites), str(n_samples),
        " ".join(str(c) for c in chrom),
        " ".join(str(int(v)) for v in pos),
        " ".join(str(int(v)) for v in meth.ravel()),
        " ".join(str(int(v)) for v in total.ravel()),
        " ".join(str(int(v)) for v in design),
    ]) + "\n", encoding="utf-8")


def find_rscript() -> str | None:
    import shutil
    for cand in RSCRIPT_CANDIDATES:
        found = shutil.which(cand) if "/" not in cand else (
            cand if Path(cand).exists() else None)
        if found:
            return found
    return None


def main() -> int:
    rscript = find_rscript()
    if rscript is None:
        print("no Rscript with DSS found; nothing to compare against")
        return 2

    meth, total, design, chrom, pos = fixture()
    oracle = REPO / "tests" / "oracle" / "dss_oracle.R"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        write_payload(meth, total, design, chrom, pos, tmp / "counts.txt")
        proc = subprocess.run([rscript, "--vanilla", str(oracle),
                               str(tmp / "counts.txt"), str(tmp / "dss.tsv")],
                              capture_output=True, text=True, timeout=3600)
        for line in proc.stdout.splitlines():
            if not line.startswith("  |"):        # DSS's progress bar
                print(line)
        if proc.returncode != 0:
            print(proc.stderr[-4000:])
            return 1
        ref = pd.read_csv(tmp / "dss.tsv", sep="\t", comment="#")

    mine = beta_binomial_test(meth, total, design, smoothing_span=500,
                              chrom=chrom, pos=pos)
    tested = np.asarray(mine["tested"], bool)
    print(f"\ntested={tested.sum()}  dss={len(ref)}  kernels=numpy")
    if tested.sum() != len(ref):
        print("SITE FILTER MISMATCH")
        return 1

    ref = ref.sort_values(["chr", "pos"]).reset_index(drop=True)
    order = np.lexsort((pos[tested], chrom[tested].astype(str)))
    idx = np.flatnonzero(tested)[order]
    cols = {
        "e1": mine["level1"], "e2": mine["level2"],
        "e1w": mine["level1_weighted"], "e2w": mine["level2_weighted"],
        "phi1": mine["phi1"], "phi2": mine["phi2"],
        "var1": mine["var1"], "var2": mine["var2"],
        "diff": mine["diff"], "se": mine["se"], "stat": mine["stat"],
        "pval": mine["pvalue"],
    }
    print(f"{'column':6s} {'max|err|':>12s} {'max|rel|':>12s} {'worst site':>10s}")
    for key, arr in cols.items():
        got = np.asarray(arr, float)[idx]
        want = ref[key].to_numpy(float)
        err = np.abs(got - want)
        rel = np.where(want != 0, err / np.abs(want), err)
        print(f"{key:6s} {np.nanmax(err):12.3e} {np.nanmax(rel):12.3e} "
              f"{int(np.nanargmax(err)):10d}")

    print("\npriors:")
    print("  mine", mine["prior1"], mine["prior2"])
    print("pval<0.05 agreement:",
          bool(((np.asarray(mine["pvalue"], float)[idx] < 0.05)
                == (ref["pval"].to_numpy(float) < 0.05)).all()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
