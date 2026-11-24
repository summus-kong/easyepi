"""A covariate the chosen test cannot use must be refused, not dropped.

`--design-factors` is read only by the deseq2 adapter. With `mannwhitney` and
the rest it used to be accepted and ignored: exit 0, no warning, and a p-value
that a reader would take as batch-adjusted. These tests pin the refusal.
"""

from __future__ import annotations

import pathlib

import pytest

from easyepi.differential import run_diff


def _write_inputs(tmp_path, extra_design=""):
    counts = tmp_path / "counts.tsv"
    rows = [["region", "s1", "s2", "s3", "s4"]]
    for i in range(30):
        base = 100 + i
        rows.append(["r%02d" % i, str(base), str(base + 3),
                     str(base + 40), str(base + 44)])
    counts.write_text("\n".join("\t".join(r) for r in rows) + "\n", encoding="utf-8")

    design = tmp_path / "design.tsv"
    design.write_text(
        "sample_id\tgroup\tbatch\n"
        "s1\tcontrol\tb1\n"
        "s2\tcontrol\tb2\n"
        "s3\ttreated\tb1\n"
        "s4\ttreated\tb2\n",
        encoding="utf-8",
    )
    return counts, design


class _Args:
    def __init__(self, **kw):
        self.counts = kw["counts"]
        self.design = kw["design"]
        self.out = kw["out"]
        self.method = kw["method"]
        self.design_factors = kw.get("design_factors")
        self.pair = kw.get("pair")
        self.ref_level = None
        self.lfc_shrink = False
        self.no_independent_filter = False
        self.equal_var = False
        self.figures = None
        self.figure_formats = None


@pytest.mark.parametrize("method", ["mannwhitney", "ttest", "wilcoxon"])
def test_a_covariate_is_refused_when_the_test_cannot_use_it(tmp_path, method):
    counts, design = _write_inputs(tmp_path)
    args = _Args(
        counts=str(counts), design=str(design), out=str(tmp_path / "out.tsv"),
        method=method, design_factors="batch",
        pair="batch" if method == "wilcoxon" else None,
    )

    with pytest.raises(ValueError) as exc:
        run_diff(args)

    message = str(exc.value)
    assert "--design-factors" in message
    assert "deseq2" in message
    assert method in message
    assert not (tmp_path / "out.tsv").exists(), (
        "refusing must not leave a result table behind"
    )


def test_without_the_flag_the_same_input_still_runs(tmp_path):
    """The guard must not break the ordinary path."""
    counts, design = _write_inputs(tmp_path)
    out = tmp_path / "out.tsv"
    args = _Args(counts=str(counts), design=str(design), out=str(out),
                 method="mannwhitney")

    assert run_diff(args) == 0
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert "pvalue" in header
