"""`easyepi.activity` drives two incompatible decoupler APIs, and must pick correctly.

decoupler moved the ULM entry point between major versions:

* 1.9.x -- ``dc.run_ulm(mat=..., source=, target=, weight=, min_n=)``, features x
  samples, raw p-values.
* 2.x   -- ``dc.mt.ulm(data, net, tmin=)``, samples x features, p-values already
  Benjamini-Hochberg corrected by ``mt/_run.py``.

The adapter reads the installed version and dispatches. These tests pin the dispatch
and the two properties that are easy to get wrong: the transpose, and applying FDR a
second time to values 2.x already adjusted.

Both generations are faked rather than imported: CI installs 1.9.x on 3.11 and 2.x on
3.12, so a test that calls whichever one happens to be present would assert different
things on different rows of the matrix. The fake records what the adapter asked for;
the real ULM arithmetic is covered by
``tests/test_optional_readers.py::test_real_decoupler_ulm_returns_measured_activity``,
which runs against the installed backend on both interpreters.
"""
from __future__ import annotations

import types

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


def _fake_generation(calls: dict, *, adjusted: bool):
    """Build a decoupler stand-in exposing exactly one generation's entry point."""

    def result(obs, sources):
        # 2 sources, named so the caller can tell which axis was used
        es = pd.DataFrame(np.zeros((len(obs), len(sources))),
                          index=list(obs), columns=list(sources))
        pv = pd.DataFrame(np.full((len(obs), len(sources)), 0.5),
                          index=list(obs), columns=list(sources))
        return es, pv

    def run_ulm(mat=None, net=None, source=None, target=None, weight=None,
                min_n=None, verbose=False, use_raw=False):
        calls["backend"] = "1.9"
        calls["matrix_shape"] = tuple(mat.shape)
        calls["index"] = list(mat.index)
        calls["sources"] = sorted(set(net[source]))
        return result(mat.index, calls["sources"])

    def mt_ulm(data, net, tmin=5, raw=False, empty=True, bsize=250_000,
               verbose=False):
        calls["backend"] = "2"
        calls["matrix_shape"] = tuple(data.shape)
        calls["index"] = list(data.index)
        calls["sources"] = sorted(set(net["source"]))
        return result(data.index, calls["sources"])

    module = types.SimpleNamespace(run_ulm=run_ulm)
    module.mt = types.SimpleNamespace(ulm=mt_ulm)
    return module


def test_dispatches_to_1_9_when_that_is_installed(monkeypatch):
    calls: dict = {}
    monkeypatch.setattr(activity, "version", lambda name: "1.9.2")
    monkeypatch.setitem(__import__("sys").modules, "decoupler",
                        _fake_generation(calls, adjusted=False))
    matrix, network = _inputs()
    es, pv, adjusted = activity._run_ulm(matrix, network, 5)

    assert calls["backend"] == "1.9"
    # features x samples in, samples x features out: (4, 30), not (30, 4)
    assert calls["matrix_shape"] == (4, 30)
    assert calls["index"] == ["A1", "A2", "B1", "B2"]
    # 1.9.x returns raw p-values, so the adapter must correct them itself
    assert adjusted is False
    assert list(es.columns) == ["TF1", "TF2"]


def test_dispatches_to_2_x_when_that_is_installed(monkeypatch):
    calls: dict = {}
    monkeypatch.setattr(activity, "version", lambda name: "2.3.1")
    monkeypatch.setitem(__import__("sys").modules, "decoupler",
                        _fake_generation(calls, adjusted=True))
    matrix, network = _inputs()
    es, pv, adjusted = activity._run_ulm(matrix, network, 5)

    assert calls["backend"] == "2"
    # 2.x also wants observations on the index, so the same single transpose
    assert calls["matrix_shape"] == (4, 30)
    assert calls["index"] == ["A1", "A2", "B1", "B2"]
    # 2.x already applied BH inside mt/_run.py
    assert adjusted is True
    assert list(es.columns) == ["TF1", "TF2"]


def test_the_transpose_is_applied_exactly_once(monkeypatch):
    """A double transpose yields the right shape with the wrong meaning."""
    calls: dict = {}
    monkeypatch.setattr(activity, "version", lambda name: "1.9.2")
    monkeypatch.setitem(__import__("sys").modules, "decoupler",
                        _fake_generation(calls, adjusted=False))
    matrix, network = _inputs()
    activity._run_ulm(matrix, network, 5)
    # 30 features x 4 samples must reach the backend as 4 x 30
    assert calls["matrix_shape"] == (4, 30)
    assert calls["matrix_shape"] != matrix.shape


@pytest.mark.parametrize("already_adjusted", [False, True])
def test_padj_is_not_corrected_twice(monkeypatch, already_adjusted):
    """`padj` must equal the backend's value when it says it already adjusted.

    Applying BH twice shrinks padj with no statistical basis. Checked on the public
    result so a refactor that loses the flag fails here.
    """
    monkeypatch.setattr(activity, "version",
                        lambda name: "2.3.1" if already_adjusted else "1.9.2")
    monkeypatch.setitem(__import__("sys").modules, "decoupler",
                        _fake_generation({}, adjusted=already_adjusted))
    matrix, network = _inputs()
    result = activity.estimate_tf_activity(matrix, network)

    assert result.pvalue.between(0, 1).all()
    assert result.padj.between(0, 1).all()
    if already_adjusted:
        assert np.allclose(result.padj.to_numpy(), result.pvalue.to_numpy())
    else:
        assert (result.padj >= result.pvalue - 1e-12).all()


def test_an_unsupported_generation_fails_loudly(monkeypatch):
    """A version outside both known APIs must raise, not guess which one to call."""
    monkeypatch.setattr(activity, "version", lambda name: "0.4.0")
    monkeypatch.setitem(__import__("sys").modules, "decoupler",
                        _fake_generation({}, adjusted=False))
    matrix, network = _inputs()
    with pytest.raises(ImportError, match="not supported"):
        activity._run_ulm(matrix, network, 5)
