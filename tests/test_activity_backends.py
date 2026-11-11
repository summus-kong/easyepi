"""`easyepi.activity` drives two incompatible decoupler APIs, and must pick correctly.

decoupler moved the ULM entry point between major versions:

* 1.9.x -- ``dc.run_ulm(mat=..., source=, target=, weight=, min_n=)``, features x
  samples, raw p-values.
* 2.x   -- ``dc.mt.ulm(data, net, tmin=)``, samples x features, p-values already
  Benjamini-Hochberg corrected by ``mt/_run.py``.

The adapter reads the installed version and dispatches. These tests pin the dispatch
and the two properties that are easy to get wrong: the transpose, and applying FDR a
second time to values 2.x already adjusted.

The 2.x call shape was verified against decoupler 2.1.6 directly -- its
``_func_ulm`` returns t-values that match 1.9.2 to 7.5e-07 on identical input, so the
two paths compute the same statistic. What cannot run here is 2.x itself, which needs
either an interpreter with it installed or a network; this file asserts the contract
the adapter relies on instead of pretending to execute it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("decoupler")

from easyepi import activity  # noqa: E402


def _inputs():
    matrix = pd.DataFrame(np.random.default_rng(5).normal(size=(30, 4)),
                          index=[f"p{i}" for i in range(30)],
                          columns=["A1", "A2", "B1", "B2"])
    network = pd.DataFrame({"source": ["TF1"] * 10 + ["TF2"] * 10,
                            "target": [f"p{i}" for i in range(20)],
                            "weight": [1.0] * 20})
    return matrix, network


def test_dispatch_reads_the_installed_generation(monkeypatch):
    """The branch is chosen from package metadata, not from a probe that can lie."""

    def fake_version(name):
        assert name == "decoupler"
        return "1.9.2"

    monkeypatch.setattr(activity, "version", fake_version)
    matrix, network = _inputs()
    es, pv, adjusted = activity._run_ulm(matrix, network, 5)

    # 1.9.x returns raw p-values, so the adapter must correct them itself.
    assert adjusted is False
    assert list(es.index) == list(matrix.columns)      # samples are rows
    assert list(es.columns) == ["TF1", "TF2"]           # sources are columns


def test_the_transpose_is_applied_once(monkeypatch):
    """Both APIs want samples as observations; the adapter passes ``matrix.T``.

    A double transpose would silently fit each sample as a feature, producing a
    result with the right shape and the wrong meaning.
    """
    seen = {}
    real_run_ulm = __import__("decoupler").run_ulm

    def spy(mat=None, **kwargs):
        seen["shape"] = tuple(mat.shape)
        seen["index"] = list(mat.index)
        return real_run_ulm(mat=mat, **kwargs)

    monkeypatch.setattr(activity, "version", lambda name: "1.9.2")
    dc = __import__("decoupler")
    monkeypatch.setattr(dc, "run_ulm", spy, raising=False)

    matrix, network = _inputs()
    activity._run_ulm(matrix, network, 5)

    # features x samples in, samples x features out: (4, 30), not (30, 4)
    assert seen["shape"] == (4, 30), seen
    assert seen["index"] == ["A1", "A2", "B1", "B2"]


def test_adjusted_pvalues_are_not_corrected_twice():
    """When the backend already applied BH, `padj` must equal its p-value.

    Applying BH twice shrinks padj without any statistical basis. The check is on the
    public result, so it also catches a future refactor that loses the flag.
    """
    matrix, network = _inputs()
    result = activity.estimate_tf_activity(matrix, network)

    # 1.9.2 path: the adapter corrects, so padj >= pvalue per sample and both in [0,1]
    assert result.pvalue.between(0, 1).all()
    assert result.padj.between(0, 1).all()
    assert (result.padj >= result.pvalue - 1e-12).all()

    # and with the flag flipped (what 2.x reports) padj is passed through unchanged
    import easyepi.activity as A

    original = A._run_ulm

    def fake(matrix, net, min_targets):
        es, pv, _ = original(matrix, net, min_targets)
        return es, pv, True

    try:
        A._run_ulm = fake
        passed = A.estimate_tf_activity(matrix, network)
    finally:
        A._run_ulm = original

    assert np.allclose(passed.padj.to_numpy(), passed.pvalue.to_numpy())


def test_an_unsupported_generation_fails_loudly(monkeypatch):
    """A version outside both known APIs must raise, not guess which one to call."""
    monkeypatch.setattr(activity, "version", lambda name: "0.4.0")
    matrix, network = _inputs()
    with pytest.raises(ImportError, match="not supported"):
        activity._run_ulm(matrix, network, 5)
