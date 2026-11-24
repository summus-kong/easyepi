"""A pairing column the chosen test cannot use must be refused, not dropped.

`--pair` is read only by the wilcoxon branch. Combined with `mannwhitney` or
`ttest` it used to be accepted and ignored, so the run produced an *unpaired*
p-value the caller had every reason to read as paired. These tests pin the
refusal, and also pin that wilcoxon still validates the pairing itself.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from easyepi.differential import run_diff

COLS = ["s1", "s2", "s3", "s4"]


class _Args:
    """Minimal stand-in for the parsed CLI namespace."""

    def __init__(self, **kw):
        self.counts = kw["counts"]
        self.design = kw["design"]
        self.out = kw["out"]
        self.method = kw.get("method", "mannwhitney")
        self.pair = kw.get("pair")
        self.design_factors = kw.get("design_factors")
        self.alpha = 0.05
        self.min_count = 10
        self.fold_change = 1.0


@pytest.fixture()
def fixture(tmp_path):
    rng = np.random.default_rng(11)
    n = 40
    header = "region\t" + "\t".join(COLS)
    lines = [header]
    for i in range(n):
        base = int(rng.integers(80, 120))
        vals = [base, base + 4, base + 45, base + 48]
        lines.append("r%02d\t%s" % (i, "\t".join(str(v) for v in vals)))
    counts = tmp_path / "counts.tsv"
    counts.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return tmp_path, counts


def _design(tmp_path, text):
    path = tmp_path / "design.tsv"
    path.write_text(text, encoding="utf-8")
    return path


BALANCED = (
    "sample_id\tgroup\tsubject\n"
    "s1\tcontrol\ts_a\n"
    "s2\ttreated\ts_b\n"
    "s3\ttreated\ts_a\n"
    "s4\tcontrol\ts_b\n"
)

ONE_SIDED = (
    "sample_id\tgroup\tsubject\n"
    "s1\tcontrol\ts_a\n"
    "s2\tcontrol\ts_b\n"
    "s3\ttreated\ts_c\n"
    "s4\ttreated\ts_d\n"
)


@pytest.mark.parametrize("method", ["mannwhitney", "ttest"])
def test_pair_is_refused_for_unpaired_methods(fixture, method):
    tmp, counts = fixture
    design = _design(tmp, BALANCED)
    out = tmp / "out.tsv"

    with pytest.raises(ValueError) as err:
        run_diff(_Args(counts=counts, design=design, out=out,
                       method=method, pair="subject"))

    msg = str(err.value)
    assert "--pair" in msg, msg
    assert "wilcoxon" in msg, msg
    assert method in msg, "the refusal should name the method that was chosen"
    assert not out.exists(), "a refused run must not leave a table behind"


def test_wilcoxon_still_rejects_a_pairing_that_carries_no_information(fixture):
    tmp, counts = fixture
    design = _design(tmp, ONE_SIDED)
    out = tmp / "out.tsv"

    with pytest.raises(ValueError) as err:
        run_diff(_Args(counts=counts, design=design, out=out,
                       method="wilcoxon", pair="subject"))

    assert "one sample per group" in str(err.value)
    assert not out.exists()


def test_wilcoxon_with_a_balanced_pairing_runs(fixture):
    tmp, counts = fixture
    design = _design(tmp, BALANCED)
    out = tmp / "out.tsv"

    assert run_diff(_Args(counts=counts, design=design, out=out,
                          method="wilcoxon", pair="subject")) == 0
    assert out.exists()


def test_plain_unpaired_run_is_untouched(fixture):
    tmp, counts = fixture
    design = _design(tmp, BALANCED)
    out = tmp / "out.tsv"

    assert run_diff(_Args(counts=counts, design=design, out=out)) == 0
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert "pvalue" in header
