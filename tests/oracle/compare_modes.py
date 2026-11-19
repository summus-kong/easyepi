"""Compare every DMLtest mode against the reference, not just the default one.

Four modes: smoothing x equal_disp, plus the no-smoothing single-replicate case
that DSS resolves by sharing one dispersion on its own. Each is a separate code
path in both implementations, and only the default one was verified before.

    python tests/oracle/compare_modes.py
"""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from easyepi.dmr_beta_binomial import beta_binomial_test  # noqa: E402

RSCRIPT_CANDIDATES = (
    "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/seurat5/bin/Rscript",
    "Rscript",
)

PER_CHR = 150
N_DIFF = 40


def fixture(n_a=3, n_b=3, seed=11):
    rng = np.random.default_rng(seed)
    n_sites = 2 * PER_CHR
    chrom = np.array(["chr1"] * PER_CHR + ["chr2"] * PER_CHR)
    gaps = rng.integers(1, 120, n_sites)
    pos = np.concatenate([
        np.cumsum(gaps[:PER_CHR]) + 1000,
        np.cumsum(gaps[PER_CHR:]) + 5000,
    ]).astype(np.int64)
    total = rng.integers(12, 45, (n_a + n_b, n_sites)).astype(np.int64)
    base = np.empty((n_a + n_b, n_sites), dtype=float)
    base[:] = 0.30
    base[n_a:, :N_DIFF] = 0.65
    base[n_a:, PER_CHR:PER_CHR + N_DIFF] = 0.65
    level = np.clip(rng.normal(base, 0.06), 0.01, 0.99)
    meth = rng.binomial(total, level).astype(np.int64)
    design = np.concatenate([np.zeros(n_a, dtype=int), np.ones(n_b, dtype=int)])
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


def find_rscript():
    """First Rscript that exists; an absolute path is trusted, a name is looked up."""
    import shutil
    for cand in RSCRIPT_CANDIDATES:
        if "/" in cand:
            if Path(cand).exists():
                return cand
        else:
            found = shutil.which(cand)
            if found:
                return found
    return None


def run_mode(name, rscript, oracle, smoothing, equal_disp, direct=False,
             **fixture_args):
    meth, total, design, chrom, pos = fixture(**fixture_args)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        write_payload(meth, total, design, chrom, pos, tmp / "counts.txt")
        proc = subprocess.run(
            [rscript, "--vanilla", str(oracle), str(tmp / "counts.txt"),
             str(tmp / "dss.tsv"), "TRUE" if smoothing else "FALSE",
             "TRUE" if equal_disp else "FALSE", "TRUE" if direct else "FALSE"],
            capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            print(f"\n{name}: ORACLE FAILED")
            print(proc.stderr[-2500:])
            return False
        ref = pd.read_csv(tmp / "dss.tsv", sep="\t", comment="#")

    mine = beta_binomial_test(meth, total, design, smoothing_span=500,
                              chrom=chrom, pos=pos,
                              equal_disp=equal_disp, smoothing=smoothing)
    tested = np.asarray(mine["tested"], bool)
    if tested.sum() != len(ref):
        print(f"\n{name}: SITE FILTER MISMATCH {tested.sum()} vs {len(ref)}")
        return False

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
    tol = {"phi1": 2e-3, "phi2": 2e-3}
    print(f"\n=== {name}  (smoothing={smoothing}, equal_disp={equal_disp}, "
          f"reps={fixture_args.get('n_a', 3)}v{fixture_args.get('n_b', 3)}) ===")
    print(f"{'column':6s} {'max|err|':>12s} {'max|rel|':>12s}   tolerance")
    ok = True
    for key, arr in cols.items():
        got = np.asarray(arr, float)[idx]
        want = ref[key].to_numpy(float)
        err = np.abs(got - want)
        rel = np.where(want != 0, err / np.abs(want), err)
        rtol = tol.get(key, 1e-6)
        bad = np.nanmax(rel) > rtol and np.nanmax(err) > 1e-9
        ok = ok and not bad
        print(f"{key:6s} {np.nanmax(err):12.3e} {np.nanmax(rel):12.3e}"
              f"   rtol={rtol:g}{'   <-- FAIL' if bad else ''}")

    same_calls = bool(((np.asarray(mine["pvalue"], float)[idx] < 0.05)
                       == (ref["pval"].to_numpy(float) < 0.05)).all())
    print(f"p<0.05 agreement: {same_calls}")
    print(f"equal_disp reported: {mine['equal_disp']}  "
          f"smoothing reported: {mine['smoothing']}")
    print(f"reference: {mine['reference']}")
    return ok and same_calls


def main() -> int:
    rscript = find_rscript()
    if rscript is None:
        print("no Rscript with DSS found")
        return 2
    oracle = REPO / "tests" / "oracle" / "dss_oracle.R"

    results = {
        "sm=TRUE  eq=FALSE": run_mode("default (smoothing, per-group disp)",
                                      rscript, oracle, True, False),
        "sm=TRUE  eq=TRUE": run_mode("smoothing, equal dispersion",
                                     rscript, oracle, True, True),
        "sm=FALSE eq=FALSE": run_mode("no smoothing, per-group disp",
                                      rscript, oracle, False, False),
        "sm=FALSE eq=TRUE": run_mode("no smoothing, equal dispersion",
                                     rscript, oracle, False, True),
        # DMLtest refuses smoothing=FALSE with a single replicate per group
        # ("There is no biological replicates in at least one condition"), even
        # though DSS:::DMLtest.noSmooth handles it by sharing one dispersion.
        # Reaching that branch therefore means calling the internal function.
        "sm=FALSE 1v3 (direct)": run_mode(
            "no smoothing, one replicate, via DSS:::DMLtest.noSmooth",
            rscript, oracle, False, False, direct=True, n_a=1, n_b=3),
    }
    print("\n" + "=" * 66)
    for k, v in results.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
