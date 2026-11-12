"""Negative-binomial count model for differential analysis, via PyDESeq2.

Why this exists
---------------
`easyepi.differential.differential_test` offers rank-sum and t-tests on whatever matrix
it is given. Those are signal comparisons: they cannot model counts, cannot absorb a
covariate, and have no concept of dispersion or independent filtering. A "significant"
call from them is exploratory.

This module does the inference properly. PyDESeq2 is the maintained Python
reimplementation of DESeq2's negative-binomial GLM, so the model is the published one
rather than something reimplemented here.

What is deliberately not hidden
-------------------------------
* **Counts must be integers.** DESeq2's likelihood is defined on counts. Feeding it
  normalised signal silently changes the variance model, so a non-integer matrix is
  rejected with the reason rather than rounded.
* **Zero-count genes are filtered, and the number is reported.** Independent filtering
  removes low-count features before testing; a caller that does not say how many were
  removed makes the FDR look better than it is.
* **The dispersion fit is reported.** If the parametric trend fails and PyDESeq2 falls
  back to a local fit, that changes the p-values and belongs in the output, not only in
  a log line.
* **The PyDESeq2 version is recorded.** A model is only reproducible with the version
  that produced it.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = ["deseq2_available", "run_deseq2"]


def deseq2_available() -> bool:
    """Whether PyDESeq2 is importable, without importing it as a side effect of --help."""
    import importlib.util

    return importlib.util.find_spec("pydeseq2") is not None


def _require() -> Dict[str, Any]:
    """Import PyDESeq2 or explain the install line."""
    try:
        import pydeseq2
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "deseq2 differential analysis needs pydeseq2; install it with "
            "`micromamba install -c conda-forge pydeseq2` or `pip install pydeseq2`"
        ) from exc
    version = getattr(pydeseq2, "__version__", None)
    if version is None:
        import importlib.metadata as md

        version = md.version("pydeseq2")
    return {"DeseqDataSet": DeseqDataSet, "DeseqStats": DeseqStats,
            "version": version}


def run_deseq2(
    region_ids: Sequence[str],
    matrix: np.ndarray,
    sample_ids: Sequence[str],
    design: Dict[str, Dict[str, str]],
    mask_a: np.ndarray,
    levels: Tuple[str, str],
    *,
    design_factors: Optional[Sequence[str]] = None,
    ref_level: Optional[Tuple[str, str]] = None,
    lfc_shrink: bool = False,
    independent_filter: bool = True,
    quiet: bool = True,
) -> Dict[str, Any]:
    """Fit a negative-binomial model and return per-feature statistics.

    ``matrix`` is regions x samples, matching the layout
    :func:`easyepi.differential.read_count_matrix` produces, and must hold integer
    counts. ``mask_a`` marks the samples in the numerator group; the denominator is
    the complement, and the contrast is therefore explicit rather than inferred from
    column order.

    ``design_factors`` names extra design columns (batch, sex, ...) that must already
    be present in ``design``. They are modelled as additive terms; an interaction
    would need an explicit formula and is not guessed here.
    """
    mod = _require()

    values = np.asarray(matrix, dtype=float)
    if values.size == 0:
        raise ValueError("count matrix is empty")
    if not np.all(np.isfinite(values)):
        raise ValueError("count matrix contains NaN or inf; DESeq2 needs complete counts")
    # DESeq2's NB likelihood is defined on counts. Rounding silently would change the
    # variance model without saying so, so refuse and name the likely cause.
    if not np.allclose(values, np.round(values), atol=1e-8):
        worst = float(np.max(np.abs(values - np.round(values))))
        raise ValueError(
            "deseq2 needs integer counts but the matrix has non-integer values "
            f"(largest deviation {worst:.3g}); this looks like normalised signal. "
            "Use --method mannwhitney for signal, or supply raw counts.")
    if np.any(values < 0):
        raise ValueError("count matrix has negative values; DESeq2 needs non-negative counts")

    n_regions, n_samples = values.shape
    if n_samples < 4:
        raise ValueError(
            f"deseq2 needs at least 4 samples to estimate dispersion (got {n_samples}); "
            "use --method mannwhitney for a small design")

    mask_a = np.asarray(mask_a, dtype=bool)
    if mask_a.sum() < 2 or (~mask_a).sum() < 2:
        raise ValueError(
            f"deseq2 needs at least 2 replicates per group "
            f"(got {int(mask_a.sum())} vs {int((~mask_a).sum())})")

    # Build the sample metadata. The contrast column is always present; extra factors
    # are validated against the design rather than discovered by name-matching.
    group_col = "easyepi_group"
    metadata: Dict[str, List[str]] = {group_col: []}
    extra = list(design_factors or [])
    for name in extra:
        metadata[name] = []
    for sid, in_a in zip(sample_ids, mask_a):
        metadata[group_col].append(levels[0] if in_a else levels[1])
        for name in extra:
            value = design.get(sid, {}).get(name, "")
            if not value:
                raise ValueError(
                    f"design factor {name!r} is missing for sample {sid!r}; "
                    "every sample in the analysis needs a value for every factor")
            metadata[name].append(value)

    for name in extra:
        if len(set(metadata[name])) < 2:
            raise ValueError(
                f"design factor {name!r} has a single level "
                f"({sorted(set(metadata[name]))}); it cannot be modelled")

    import pandas as pd

    counts_df = pd.DataFrame(
        np.round(values).astype(int).T,      # PyDESeq2 wants samples x features
        index=list(sample_ids), columns=list(region_ids))
    meta_df = pd.DataFrame(metadata, index=list(sample_ids))

    factors = [group_col] + extra
    # DESeq2 parameterises relative to a reference level. The contrast below asks
    # for `levels[0]` over `levels[1]`, so `levels[1]` must be the reference: with
    # the reference the other way round PyDESeq2 reports the reciprocal ratio,
    # which flips the sign of every log2FoldChange while leaving p-values
    # untouched -- an error that reads as real biology rather than a bug.
    # PyDESeq2 0.5 deprecates `ref_level`: the contrast names both operands, so the
    # reference only affects the GLM parameterisation and not the reported
    # direction. Passing it produced a DeprecationWarning and changed nothing, so
    # the adapter neither accepts nor emits it -- a parameter that is silently
    # ignored is worse than one that does not exist.
    if ref_level:
        if list(ref_level)[0] not in levels:
            raise ValueError(
                f"reference level {list(ref_level)[0]!r} is not one of the "
                f"contrast levels {levels}")
    dds = mod["DeseqDataSet"](
        counts=counts_df, metadata=meta_df, design="~" + " + ".join(factors),
        quiet=quiet)
    dds.deseq2()

    # PyDESeq2 reports log2(second / first) for contrast=[factor, first, second],
    # i.e. the *opposite* of R DESeq2's resultsNames convention. Measured, not
    # assumed: on data where group B is doubled and group A is halved, the
    # contrast [g, "A", "B"] returns -0.98 for the doubled features and +1.04 for
    # the halved ones, and swapping the operands flips both signs. The adapter's
    # contract is log2(levels[0] / levels[1]), so levels[1] goes first.
    stats = mod["DeseqStats"](dds, contrast=[group_col, levels[1], levels[0]],
                              independent_filter=independent_filter, quiet=quiet)
    use_shrink = lfc_shrink
    stats.summary(lfc_shrink=use_shrink)
    res = stats.results_df

    # Report what the model actually did, so the numbers can be read with the right
    # caveats.
    dispersion_fit = getattr(dds, "disp_function_type", None)
    fitted = getattr(dds, "varm", {}).get("dispersions", None) if hasattr(dds, "varm") else None
    tested = int(np.isfinite(res["pvalue"].to_numpy(dtype=float)).sum())

    ordered = res.reindex(list(region_ids))
    return {
        "baseMean": ordered["baseMean"].to_numpy(dtype=float),
        "log2FoldChange": ordered["log2FoldChange"].to_numpy(dtype=float),
        "lfcSE": ordered["lfcSE"].to_numpy(dtype=float) if "lfcSE" in ordered else None,
        "stat": ordered["stat"].to_numpy(dtype=float) if "stat" in ordered else None,
        "pvalue": ordered["pvalue"].to_numpy(dtype=float),
        "padj": ordered["padj"].to_numpy(dtype=float),
        "method": "deseq2",
        "info": {
            "pydeseq2_version": mod["version"],
            "design": "~" + " + ".join(factors),
            "contrast": [group_col, levels[0], levels[1]],
            "contrast_reported_as": f"log2({levels[0]} / {levels[1]})",
            # The reference no longer parameterises the model (PyDESeq2 0.5 dropped
            # it); the contrast operands are what determine the direction.
            "reference_level": levels[1],
            "design_factors": extra,
            "lfc_shrink": bool(use_shrink),
            "independent_filter": bool(independent_filter),
            "treated_as_na_before_filtering": int(n_regions - tested),
            "replicates_a": int(mask_a.sum()),
            "replicates_b": int((~mask_a).sum()),
            "dispersion_fit": dispersion_fit,
            "dispersion_fitted": bool(fitted is not None),
        },
    }
