"""Regulator activity estimates using the published decoupler ULM implementation."""
from __future__ import annotations
from importlib.metadata import version
import numpy as np
import pandas as pd
from easyepi.core.stats import benjamini_hochberg


def estimate_tf_activity(matrix: pd.DataFrame, network: pd.DataFrame, min_targets: int = 5) -> pd.DataFrame:
    """Fit decoupler ULM to a normalized features-by-samples measurement matrix.

    The network requires source (TF), target (feature ID), weight. Motif-derived
    networks yield motif-associated accessibility estimates, not direct binding
    probabilities or chromVAR scores. Prepare and normalize inputs explicitly.
    """
    if matrix.empty or not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("Activity matrix must be nonempty with unique feature and sample IDs")
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all() or matrix.shape[0] < 3:
        raise ValueError("Activity matrix requires finite measured values and at least three features")
    if not {"source","target","weight"}.issubset(network.columns):
        raise ValueError("Network requires source, target and weight columns")
    net = network[["source","target","weight"]].copy()
    if net.isna().any().any() or net.duplicated(["source","target"]).any():
        raise ValueError("Network edges must be complete and unique")
    net["weight"] = pd.to_numeric(net["weight"],errors="raise")
    if not np.isfinite(net.weight).all() or min_targets < 2:
        raise ValueError("Network weights must be finite and min_targets must be >=2")
    net = net[net.target.isin(matrix.index)]
    keep = net.groupby("source")["target"].nunique()
    net = net[net.source.isin(keep[keep >= min_targets].index)]
    if net.empty:
        raise ValueError("No regulator has sufficient measured target features")
    for source, group in net.groupby("source"):
        weights = group.set_index("target").weight.reindex(matrix.index,fill_value=0).to_numpy()
        if np.std(weights) == 0:
            raise ValueError(f"{source}: constant network weights cannot identify ULM activity")
    if any(matrix[col].std() == 0 for col in matrix.columns):
        raise ValueError("Constant sample profiles cannot identify ULM activity")
    estimates, pvalues, pvalues_are_adjusted = _run_ulm(matrix, net, min_targets)
    rows = [{"sample_id":sample,"tf":tf,"activity":float(estimates.loc[sample,tf]),
             "pvalue":float(pvalues.loc[sample,tf]),"method":"decoupler_ULM"}
            for sample in estimates.index for tf in estimates.columns]
    result = pd.DataFrame(rows)
    if pvalues_are_adjusted:
        # decoupler 2.x already applies Benjamini-Hochberg along the sources of
        # each observation, which is the same axis and the same correction this
        # adapter applies for 1.9.x. Correcting twice would shrink `padj` without
        # any statistical basis, so the values are taken as given.
        result["padj"] = result["pvalue"].to_numpy()
    else:
        result["padj"] = result.groupby("sample_id")["pvalue"].transform(
            lambda x: benjamini_hochberg(x.to_numpy()))
    return result


def _run_ulm(matrix: pd.DataFrame, net: pd.DataFrame, min_targets: int):
    """Call decoupler ULM on either API generation; return (es, pv, pv_adjusted).

    decoupler moved the ULM entry point between major versions:

    * 1.9.x -- ``dc.run_ulm(mat=..., net=..., source=, target=, weight=, min_n=)``,
      where ``mat`` is features-by-samples and the returned p-values are raw.
    * 2.x  -- ``dc.mt.ulm(data, net, tmin=)``, where ``data`` is a DataFrame of
      samples-by-features (verified in ``decoupler/pp/data.py:_extract``, which
      reads ``data.index`` as observations and ``data.columns`` as features),
      ``net`` keeps the same source/target/weight columns, and the p-values
      come back already Benjamini-Hochberg corrected by ``mt/_run.py``.

    Both return DataFrames indexed by sample with one column per source, so the
    caller does not need to know which ran. The version is read from the package
    metadata rather than from a probe call, so a broken install reports itself
    instead of being mistaken for the other generation.
    """
    import decoupler as dc

    installed = version("decoupler")
    major = int(installed.split(".")[0])
    if major >= 2:
        estimates, pvalues = dc.mt.ulm(matrix.T, net, tmin=min_targets,
                                       verbose=False)
        return estimates, pvalues, True
    if installed.startswith("1.9."):
        estimates, pvalues = dc.run_ulm(mat=matrix.T, net=net, source="source",
                                        target="target", weight="weight",
                                        min_n=min_targets, verbose=False,
                                        use_raw=False)
        return estimates, pvalues, False
    raise ImportError(
        f"decoupler {installed} is not supported by this adapter; it handles the "
        "1.9.x and 2.x ULM entry points. Install 'easyepi[activity]' for a tested "
        "version.")
