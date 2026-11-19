# DSS replacement in pure Python — final state

`easyepi dmr` no longer calls R. Differential methylation on CpG sites is
computed in numpy/scipy and reproduces `DSS::DMLtest(smoothing = TRUE)` —
including the shrunken dispersion — to the limits of DSS's own numerics.

DSS remains the reference and the test oracle; it is not a runtime dependency.
The package needs numpy and scipy, both already required by every other module,
and nothing else: there is no C extension and no compile step at install or at
run time.

---

## 1. Why this is not a small change

DSS's test is not a textbook beta-binomial test. Reproducing it means
reproducing five coupled steps, in this order, per group:

```
1. estprob <- compute.mean.Smooth(x, n, allchr, allpos, span)   smoothed level
2. phi     <- est.dispersion.BSseq(x, n, estprob, ncores)       logN-shrunk dispersion
3. w       <- 1/(1+(n-1)*phi);  w <- w/mean(w)                  site weights
4. estprob <- compute.mean.Smooth(x*w, n*w, allchr, allpos, span) level, recomputed
5. v       <- compute.var.Smooth(estprob[,1], n*w, phi, span)   variance
--------------------------------------------------------------------
diff = estprob1 - estprob2
vv   = var1 + var2   (with DSS's two-sided guard)
se   = sqrt(vv);  p = 2*pnorm(-|diff/se|)
```

All five matter. Reproducing only 1, 2 and 4 gives a `diff` that looks right
while `se` and `p` are wrong by orders of magnitude — the failure mode that
took longest to notice.

### The six details that actually break a reimplementation

| # | DSS behaviour | What a plausible rewrite does instead |
|---|---|---|
| 1 | `est.prior.BSseq.logN` uses `IQR/1.39`, no lower bound | `1.4826*IQR`, sometimes with a `sqrt(2/n)` floor |
| 2 | `compute.mean.Smooth` adds `N+1` and `X+const` before smoothing, with `const = mean(X/N)` | smooths the bare ratio — loses the pseudocount entirely |
| 3 | `compute.var.Smooth` is given `n*wt`, then adds `const = 1/nreps` | given the raw depth, or no `1/nreps` |
| 4 | `phi` is recycled down the columns of a `(sites, reps)` matrix | transposed (`phi[None, :]`), which raises or silently misaligns |
| 5 | the variance window is `[ipos-left, ipos+right)` — **right edge exclusive** | copies `filter.c`'s inclusive convention, so every variance drifts |
| 6 | the diagonal term is added once **per site in the window** | added once per site, leaving `se` ~13 % too small |

Detail 5 is worth stating twice: `filter.c` and `compute_var_smooth.c` use the
same half-widths with *different* edge conventions. The off-by-one is in DSS
itself, and copying either convention into the other place is wrong.

---

## 2. What is verified, and against what

Fixture: 300 CpG sites over **two chromosomes**, **irregular spacing**
(1–120 bp gaps), 3 vs 3 replicates, coverage 12–45, effect on 40 sites per
chromosome, and scatter *within* each group so the dispersion is identifiable.

The oracle (`tests/oracle/dss_oracle.R`) computes every intermediate by calling
DSS's own functions in DSS's own order, and then also runs `DMLtest`. It
asserts internally that the two agree to 1e-12 before writing anything, so a
green comparison below is a comparison with **DSS**, not with the test's idea
of DSS.

Max absolute error over all 300 sites:

| quantity | max abs error | max rel error |
|---|---|---|
| `e1`, `e2` (first-pass level) | 5.0e-16 | 1.8e-15 |
| `e1w`, `e2w` (weighted level) | 7.5e-12 | 2.0e-11 |
| `phi1`, `phi2` (shrunken dispersion) | 1.8e-11 | 4.0e-10 |
| `var1`, `var2` (smoothed variance) | 1.3e-13 | 1.6e-10 |
| `diff` | 7.5e-12 | 2.2e-08 |
| `se` | 1.5e-12 | 3.8e-11 |
| `stat` | 5.2e-10 | 2.2e-08 |
| `pval` | 1.1e-10 | 7.2e-09 |
| sites called at p < 0.05 | identical set | — |

The prior matches exactly:

```
mine  (-3.2037106937234414, 1.0422753711192465)
DSS   (-3.2037106937234414, 1.0422753711192462)
```

The residual is the last digit of `IQR/1.39`, i.e. R's `median`/`IQR`
arithmetic — not a formula difference.

`phi` is the loosest quantity because it is a bounded 1-D optimisation and DSS
fixes `tol = 0.001` on `log(phi)`. Two correct implementations may differ by
that much and no more; the test pins `phi` at `rtol=2e-3` and everything
downstream of it at `rtol=1e-6`.

### The numpy kernels are the C kernels

DSS's smoothing lives in C. `easyepi/kernels/` reimplements the three functions
in numpy so there is no compiler to depend on. `tests/test_dss_kernels.py`
checks each one against a literal transcription of the C loops (24 tests) across
window widths 1–500, uniform / irregular / duplicate positions, and window
sizes down to a single site. The full DSS comparison above was run with the
numpy kernels, so the table is the numpy path, not a compiled stand-in.

---

## 3. What changed in the product path

`easyepi/dmr_replicate_aware.py` calls `beta_binomial_test` with the **real**
`chrom`/`pos`. The old R bridge fabricated `chr1` with positions `1..n_sites`,
which put every smoothing window on a one-base-per-site axis and could smooth
across chromosome boundaries — `hasCoverage` and `nitem_bin` are both genomic,
so the fabricated axis changed which sites were tested and how.

`beta_binomial_test` now applies DSS's own coverage filter (`hasCoverage`,
`ws = 200`, threshold `2 * n_samples`, in both groups) *before* any statistic is
computed, because the dropped sites move every later window. Dropped sites come
back as `NaN` with `tested == False`, and the caller filters them out before BH
correction — a `NaN` p-value would otherwise poison `p.adjust` and the region
merge. `stats["dropped_by_dss_coverage_filter"]` reports the count.

End-to-end check on a two-chromosome fixture (240 sites, 3 vs 3, effect on 80
sites): 240 rows, no `NaN`, `padj` inside [0, 1], both chromosomes preserved,
85 significant of which 82 are on true-difference sites, median |diff| on
significant sites 0.389.

### The unsmoothed variant is reachable from the CLI

`easyepi dmr --no-smoothing` selects `DMLtest(smoothing = FALSE)`. It is not the
default, and the run summary distinguishes the two DSS variants so a table can
never be attributed to the wrong one:

```
model      dss_beta_binomial                 | dss_beta_binomial_no_smoothing
method     "... on smoothed levels ..."      | "... on unsmoothed per-site levels ..."
stats      smoothing = dss_compute_mean_smooth | dss_compute_mean_noSmooth
```

`--no-smoothing` is ignored on the pooled path (fewer than two replicates per
group, or `--pooled`), and the summary then still says `pooled_binomial`: the
flag reaches only the DSS branch, so asking for it cannot relabel a pooled run.

---

## 4. Files

```
easyepi/dmr_beta_binomial.py     rewritten; steps 1-5 of section 1
easyepi/dmr_dispersion.py        rewritten; est.prior.BSseq.logN + dispersion.shrinkage.BSseq
easyepi/kernels/__init__.py      rewritten; the three DSS kernels in numpy
easyepi/dmr_replicate_aware.py   uses beta_binomial_test with real coordinates
easyepi/dmr.py                   summary fields name the model that actually ran
easyepi/cli.py                   --no-smoothing selects DMLtest(smoothing = FALSE)
tests/test_beta_binomial_matches_dss.py   one comparison per DMLtest mode (5)
tests/test_dss_kernels.py                 numpy kernels vs literal C loops (24)
tests/test_dmr_summary_is_honest.py       the summary cannot misdescribe itself
tests/test_dmr_replicate_aware.py         method name and R-skip updated
tests/oracle/dss_oracle.R                 DSS reference values, all modes (test-only)
tests/oracle/compare_modes.py             prints every mode's error table
```

`easyepi/dmr_r.py` is unchanged and still correct. It is no longer imported from
the product path, so it is now an optional route rather than the only one.

Removed: `easyepi/kernels/dss_kernels.c` and the ctypes/`BACKEND` build probe.
The package no longer compiles anything.

---

## 5. Test status

Server, `easyepi` conda environment (Python 3.11.16, numpy 2.2.6, scipy 1.17.1),
with DSS 2.58.0 installed so the oracle tests run rather than skip:

```
$ python -m pytest tests -q
722 passed, 6 warnings in 127.02s
```

The 2025 round added the mode coverage: `test_beta_binomial_matches_dss.py` now
holds one comparison per `DMLtest` mode (5 oracle tests), and
`test_dss_kernels.py` holds the 24 kernel-vs-C-loop tests. The oracle tests skip
cleanly when R or DSS is absent, which is the ordinary case for an installed
package, so the suite is meaningful on a machine without DSS and *stronger* on
one with it.

---

## 6. Provenance

- DSS 2.58.0 (R library path on the analysis server:
  `…/envs/seurat5/lib/R/library/DSS`), all of:
  `DMLtest`, `DMLtest.Smooth`, `waldTest.DML`, `compute.mean.Smooth`,
  `compute.var.Smooth`, `compute.waldStat.Smooth`, `smooth.chr`, `hasCoverage`,
  `est.prior.BSseq.logN`, `est.dispersion.BSseq`, `dispersion.shrinkage.BSseq`,
  `dbb`
- C behaviour transcribed from `DSS/src/nitem.c`, `src/filter.c`,
  `src/compute_var_smooth.c`, `DSS.h`
- `BINSIZE = 50`, `rho = 0.8`, coverage window `ws = 200` with threshold
  `2 * ncol`, variance guard `vv<1e-5 -> 1e-4` and `vv>1-1e-5 -> 1-1e-4`
  applied to the **sum**

## 7. Every `DMLtest` mode, and what is still not covered

`DMLtest` has four combinations of `smoothing` x `equal.disp`, plus a
one-replicate case that the no-smoothing path resolves by itself. All five are
implemented and **all five are compared against DSS numerically**
(`tests/test_beta_binomial_matches_dss.py`, one test each):

| mode | how DSS branches | verified |
|---|---|---|
| `smoothing=TRUE, equal.disp=FALSE` | `DMLtest.Smooth`, per-group `est.dispersion.BSseq` | yes |
| `smoothing=TRUE, equal.disp=TRUE` | `cbind` the groups, one dispersion | yes |
| `smoothing=FALSE, equal.disp=FALSE` | `DMLtest.noSmooth`, per-group | yes |
| `smoothing=FALSE, equal.disp=TRUE` | `cbind`, one dispersion | yes |
| `smoothing=FALSE`, one replicate per group | `DMLtest.noSmooth` shares a dispersion | yes |

### The no-smoothing path is not "span = 0"

It is a separate pair of functions with its own arithmetic and, importantly, its
own variance guard:

```
compute.waldStat.noSmooth : vv[vv < 1e-05] <- 1e-05                  (one-sided)
compute.waldStat.Smooth   : vv[vv < 1e-05] <- 1e-04
                            vv[vv > 1 - 1e-05] <- 1 - 1e-04          (two-sided)
```

Two guards, deliberately not factored into one helper. Reusing the smoothed
floor here would change every no-smoothing p-value at low coverage.

`compute.mean.noSmooth` keeps the pseudocount as well: `(rowSums(X) + const) /
(rowSums(N) + 1)` with the same `const = mean(X/N)`. It drops the window, not the
correction.

### One DSS behaviour worth knowing before calling the internal path

`DMLtest` refuses `smoothing = FALSE` with a single replicate per group
("There is no biological replicates in at least one condition"), yet
`DSS:::DMLtest.noSmooth` handles that combination by sharing one dispersion. The
public entry point and the internal one therefore disagree about which inputs
are legal. The oracle takes a `direct` route to reach the internal branch, and
Python reports `equal_disp = True` in that case because that is what was *used*,
not what the caller asked for.

### Still not covered

- **Performance at genome scale** (measured, see below). The variance kernel is
  O(m^2) per site where m is the sites per window, i.e. the same algorithmic
  cost as DSS's own C double loop, so this is not asymptotically worse than the
  reference -- but the constant is Python, not C.

### Performance: measured, then cut by 2.7x

`tests/oracle/bench_kernels.py` times each stage separately at 1k / 5k / 20k
sites on one chromosome, so a slow stage cannot hide behind a fast one. At 20k
sites with ~101 CpG per 500 bp window, before optimisation:

| stage | time | share |
|---|---|---|
| variance kernel x3 replicates | 7.17 s | 94.5% |
| dispersion (per-site optimiser) | 0.39 s | 5.2% |
| smoothing | 0.02 s | 0.3% |

123.6 us/site, i.e. 5.9 h extrapolated to a 28M-CpG genome, two groups.

Profiling the kernel body (`tests/oracle/profile_kernel.py`) split that further:

| per-site component | time |
|---|---|
| building the lag matrix | 31.3 us |
| `rho ** lag` | 96.4 us |
| `u^T W u` (BLAS) | 7.9 us |

So the transcendental power call was 60% of the whole kernel, and the exponent
only ever takes the integer values `0 .. span//BINSIZE + 1` (12 of them at
span 500). Replacing `rho ** lag.astype(float)` with a table lookup is exact --
`rho ** np.arange(n)` and `rho ** lag` agree bitwise -- and
`test_rho_powers_table_is_bit_identical` asserts that with `array_equal`, not
`allclose`.

`nitem_bin` was also still a Python two-pointer loop, though it does not show in
the table above because it runs once per call rather than once per replicate.
Since positions are sorted, the window is two binary searches; 20.6 ms dropped
to 1.0 ms.

After both changes, at the same size:

| | before | after | ratio |
|---|---|---|---|
| variance kernel x3 | 7.17 s | 2.56 s | 2.8x |
| `nitem_bin` | 20.6 ms | 1.0 ms | 21x |
| total, per site | 379.4 us | 138.8 us | **2.7x** |
| genome, 28M CpG, 2 groups | 5.9 h | **2.2 h** | |

The DSS oracle tests were re-run after the change and every intermediate is
unchanged, which is the point: the optimisation touches how a constant is
looked up, not what is computed.

The kernel remains the dominant cost (92%), and it is now genuinely O(m^2)
indexing and BLAS work per site rather than a mis-spent power call. Going
further would mean changing the algorithm -- e.g. exploiting the Toeplitz
structure of the lag weights to get O(m log m) per window -- which would no
longer be a line-by-line reproduction of DSS. That is a deliberate stopping
point, not an oversight.

Extrapolations above are linear because the per-site stages are O(1) in the
number of sites; only the window size (not the site count) sets the inner cost.
They are extrapolations from 20k sites, not a full genome run.
- Region-level DSS functions (`callDMR`, `DMR.plot`, the count-based variants)
  are out of scope. This port covers the site-level two-group test only: both
  smoothing modes, both dispersion rules, and the replicate filter.
