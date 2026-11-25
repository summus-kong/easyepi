"""Pin the wide table's "one column per stage" promise.

The three-column reading that prompted this test came from a stale fixture with
only three time points, not from the writer. The invariant is worth pinning
anyway: a header that silently drops a stage is indistinguishable from a design
that never had it.
"""

from __future__ import annotations

import csv

import pytest

np = pytest.importorskip("numpy")

from easyepi import timecourse as tc
from easyepi.timecourse_cli import run_timecourse


class _Args:
    def __init__(self, counts, design, out):
        self.counts = counts
        self.design = design
        self.out = out


def _write_case(tmp_path, stages, reps=2):
    """Counts and design for len(stages) time points, each with ``reps`` samples."""
    sample_ids = ["S_%s%d" % (s, r + 1) for s in stages for r in range(reps)]
    rows = []
    for i, s in enumerate(stages):
        # Feature i peaks at stage i and is near-zero everywhere else.
        values = []
        for j in range(len(stages)):
            level = 950.0 if j == i else 20.0
            values.extend([level + k for k in range(reps)])
        rows.append(("peak_%s" % s, values))

    counts = tmp_path / "counts.tsv"
    with counts.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["feature"] + sample_ids)
        for name, values in rows:
            w.writerow([name] + values)

    design = tmp_path / "design.tsv"
    with design.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "time"])
        for s in stages:
            for r in range(reps):
                w.writerow(["S_%s%d" % (s, r + 1), s])
    return counts, design


def _header(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return next(csv.reader(fh, delimiter="\t"))


@pytest.mark.parametrize("n_stages", [3, 4, 5, 7])
def test_every_stage_gets_its_own_column(tmp_path, n_stages):
    stages = ["ABCDEFG"[i] for i in range(n_stages)]
    counts, design = _write_case(tmp_path, stages)
    out = tmp_path / "out.tsv"
    run_timecourse(_Args(str(counts), str(design), str(out)))

    header = _header(out)
    spec = [c for c in header if c.startswith("spec_")]
    assert spec == ["spec_%s" % s for s in stages], header


def test_a_stage_column_is_blank_not_zero_when_unmeasurable(tmp_path):
    """An all-zero feature must not look like a measured zero anywhere."""
    stages = ["A", "B", "C"]
    counts, design = _write_case(tmp_path, stages)
    body = counts.read_text(encoding="utf-8")
    counts.write_text(
        body + "silent\t0\t0\t0\t0\t0\t0\n", encoding="utf-8")

    out = tmp_path / "out.tsv"
    run_timecourse(_Args(str(counts), str(design), str(out)))

    with open(out, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    silent = [r for r in rows if r["feature"] == "silent"]
    assert len(silent) == 1
    for stage in stages:
        assert silent[0]["spec_%s" % stage] == "", silent[0]


def test_the_stage_scores_reproduce_tau(tmp_path):
    """max(spec_s) == tau, so the two tables cannot disagree."""
    stages = ["A", "B", "C", "D", "E"]
    counts, design = _write_case(tmp_path, stages)
    out = tmp_path / "out.tsv"
    run_timecourse(_Args(str(counts), str(design), str(out)))

    with open(out, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            scores = [float(row["spec_%s" % s]) for s in stages]
            assert max(scores) == pytest.approx(float(row["tau"]), abs=1e-5)
