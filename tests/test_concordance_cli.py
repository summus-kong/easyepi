"""`easyepi concordance`: the input contract, the output contract, and the traps.

The two traps worth a test each are the ones that would pass silently:

* NaN in the JSON. ``json.dumps`` writes bare ``NaN``, so a report can be
  produced that every strict parser rejects. The suite therefore parses with
  ``parse_constant`` set to raise, which is the only way the defect shows up.
* "not computed" collapsing into "zero". An empty Jaccard written as 0.0 reads
  as "these peak sets share nothing", which is a biological claim, when the
  truth is that nothing was measured.
"""
from __future__ import annotations

import csv
import json

import numpy as np
import pytest

from easyepi.cli import build_parser
from easyepi.concordance import replicate_concordance
from easyepi.concordance_cli import run_concordance


def _strict_loads(text: str):
    """Parse JSON the way a strict consumer does: bare NaN/inf is an error."""
    def reject(token):
        raise ValueError(f"non-JSON numeric constant {token!r}")

    return json.loads(text, parse_constant=reject)


def _design(path, groups):
    lines = ["sample_id\tgroup"]
    lines += [f"{s}\t{g}" for s, g in groups.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _matrix(path, ids, *, seed=0, group_offsets=None):
    """Features x samples table honouring the `diff --counts` layout."""
    rng = np.random.default_rng(seed)
    common = rng.normal(size=(60, 1))
    cols = []
    for i, name in enumerate(ids):
        offset = 0.0 if group_offsets is None else group_offsets[name]
        cols.append(common[:, 0] + offset + rng.normal(scale=0.1, size=60))
    mat = np.column_stack(cols)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["feature"] + list(ids))
        for r in range(mat.shape[0]):
            w.writerow([f"f{r}"] + [f"{v:.6f}" for v in mat[r]])
    return path


def _peaks(path, regions):
    with path.open("w", encoding="utf-8") as fh:
        for chrom, start, end in regions:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    return path


def _cli(argv):
    """Run the real command through the parser, as a user would."""
    from easyepi.cli import main

    return main(argv)


def _args(**kw):
    base = dict(peaks=None, matrix=None, design=None, method="spearman",
                min_samples=2, out=None, figure_formats=None, no_figures=True,
                verbose=False, log_file=None, quiet=True)
    base.update(kw)
    return type("Args", (), base)()


# --------------------------------------------------------------------- parser
def test_parser_accepts_both_input_modes():
    parser = build_parser()
    a = parser.parse_args(["concordance", "--design", "d.tsv",
                           "--matrix", "m.tsv", "--out", "o"])
    assert a.matrix == "m.tsv" and a.peaks is None
    b = parser.parse_args(["concordance", "--design", "d.tsv",
                           "--peaks", "A=a.bed", "B=b.bed", "--out", "o"])
    assert b.peaks == [["A=a.bed", "B=b.bed"]] and b.matrix is None


def test_parser_requires_design_and_out():
    parser = build_parser()
    for missing in (["--out", "o"], ["--design", "d.tsv"]):
        argv = ["concordance", "--matrix", "m.tsv"] + missing
        with pytest.raises(SystemExit):
            parser.parse_args(argv)


# --------------------------------------------------------------- matrix mode
def test_matrix_mode_writes_json_tsv_and_has_no_jaccard(tmp_path):
    ids = ["A1", "A2", "B1", "B2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    matrix = _matrix(tmp_path / "m.tsv", ids,
                     group_offsets={"A1": -3, "A2": -3, "B1": 3, "B2": 3})
    out = tmp_path / "out"

    assert run_concordance(_args(matrix=str(matrix), design=str(design),
                                 out=str(out))) == 0
    assert (out / "concordance.json").exists()
    assert (out / "concordance_pairs.tsv").exists()

    payload = _strict_loads((out / "concordance.json").read_text())
    assert payload["level"] == "matrix"
    assert payload["n_samples"] == 4
    # A count matrix has no coordinates, so every Jaccard is null, never 0.0.
    assert all(p["jaccard"] is None for p in payload["pairs"])
    assert any("no genomic coordinates" in n for n in payload["notes"])


def test_matrix_mode_correlation_matches_the_api(tmp_path):
    """The CLI must not compute its own numbers."""
    ids = ["A1", "A2", "A3", "B1", "B2", "B3"]
    groups = {s: s[0] for s in ids}
    design = _design(tmp_path / "d.tsv", groups)
    matrix = _matrix(tmp_path / "m.tsv", ids, seed=5)
    out = tmp_path / "out"
    run_concordance(_args(matrix=str(matrix), design=str(design), out=str(out)))

    payload = _strict_loads((out / "concordance.json").read_text())

    from easyepi.differential import read_count_matrix

    _r, mat, sids = read_count_matrix(matrix)
    expected = replicate_concordance(sids, mat, groups=groups)["summary"]
    assert payload["summary"]["within_group_mean_correlation"] == pytest.approx(
        expected["within_group_mean_correlation"])
    assert payload["summary"]["grouping_supported"] == \
        expected["grouping_supported"]


# ---------------------------------------------------------------- peaks mode
def test_peaks_mode_computes_jaccard_and_records_groups(tmp_path):
    ids = ["A1", "A2", "B1"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    shared = [("chr1", 100, 200), ("chr1", 300, 400)]
    a1 = _peaks(tmp_path / "A1.bed", shared)
    a2 = _peaks(tmp_path / "A2.bed", shared)
    b1 = _peaks(tmp_path / "B1.bed", [("chr1", 900, 950)])
    out = tmp_path / "out"

    run_concordance(_args(peaks=[f"{a1}", f"{a2}", f"{b1}"],
                          design=str(design), out=str(out)))
    payload = _strict_loads((out / "concordance.json").read_text())

    assert payload["level"] == "peaks"
    by_pair = {(p["sample_a"], p["sample_b"]): p for p in payload["pairs"]}
    assert by_pair[("A1", "A2")]["jaccard"] == pytest.approx(1.0)
    assert by_pair[("A1", "B1")]["jaccard"] == pytest.approx(0.0)
    assert payload["groups"] == {"A1": "A", "A2": "A", "B1": "B"}


def test_nested_peaks_values_from_the_parser_are_flattened(tmp_path):
    """argparse gives a list of lists; the loader must accept both shapes.

    ``action="append"`` with ``nargs="+"`` produces one inner list per flag
    occurrence. Flattening with ``str()`` turned an inner list into a single
    bogus path, so the command failed with a not-found message quoting a Python
    list. This pins both the nested and the flat form.
    """
    from easyepi.concordance import read_region_sets

    ids = ["A1", "A2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    p1 = _peaks(tmp_path / "A1.bed", [("chr1", 1, 9)])
    p2 = _peaks(tmp_path / "A2.bed", [("chr1", 5, 12)])

    nested = read_region_sets([[str(p1), str(p2)]])
    flat = read_region_sets([str(p1), str(p2)])
    assert set(nested) == set(flat) == {"A1", "A2"}
    assert nested["A1"] == flat["A1"]

    # And through the real parser, end to end.
    out = tmp_path / "out"
    assert _cli(["concordance", "--peaks", str(p1), str(p2),
                 "--design", str(design), "--out", str(out),
                 "--no-figures"]) == 0
    payload = _strict_loads((out / "concordance.json").read_text())
    assert payload["n_samples"] == 2


def test_empty_region_sets_are_null_not_zero(tmp_path):
    """Two empty peak files have no Jaccard; that is not a Jaccard of zero."""
    ids = ["A1", "A2", "B1"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    empty1 = _peaks(tmp_path / "A1.bed", [])
    empty2 = _peaks(tmp_path / "A2.bed", [])
    one = _peaks(tmp_path / "B1.bed", [("chr1", 10, 20)])
    out = tmp_path / "out"

    run_concordance(_args(peaks=[str(empty1), str(empty2), str(one)],
                          design=str(design), out=str(out)))
    payload = _strict_loads((out / "concordance.json").read_text())
    pair = next(p for p in payload["pairs"]
                if {p["sample_a"], p["sample_b"]} == {"A1", "A2"})
    assert pair["jaccard"] is None


def test_all_empty_peak_files_are_refused(tmp_path):
    ids = ["A1", "A2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    e1 = _peaks(tmp_path / "A1.bed", [])
    e2 = _peaks(tmp_path / "A2.bed", [])
    with pytest.raises(ValueError, match="nothing to compare"):
        run_concordance(_args(peaks=[str(e1), str(e2)], design=str(design),
                              out=str(tmp_path / "o")))


# ------------------------------------------------------------ output contract
def test_json_is_strict_with_nan_never_emitted(tmp_path):
    """json.dumps writes bare NaN; the report must contain none of it."""
    ids = ["A1", "A2", "A3"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    e = _peaks(tmp_path / "A1.bed", [])
    e2 = _peaks(tmp_path / "A2.bed", [])
    b = _peaks(tmp_path / "A3.bed", [("chr1", 5, 9)])
    out = tmp_path / "out"
    run_concordance(_args(peaks=[str(e), str(e2), str(b)],
                          design=str(design), out=str(out)))
    text = (out / "concordance.json").read_text()
    assert "NaN" not in text and "Infinity" not in text
    _strict_loads(text)  # would raise on a non-JSON constant


def test_pairs_tsv_writes_blank_for_uncomputed(tmp_path):
    ids = ["A1", "A2", "B1"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    m = _matrix(tmp_path / "m.tsv", ids)
    out = tmp_path / "out"
    run_concordance(_args(matrix=str(m), design=str(design), out=str(out)))
    with (out / "concordance_pairs.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows and all(r["jaccard"] == "" for r in rows)
    # And no literal "nan" string, which some table tools read as zero.
    assert all("nan" not in str(v).lower() for r in rows for v in r.values())


def test_figures_are_written_when_requested(tmp_path):
    ids = ["A1", "A2", "B1"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    m = _matrix(tmp_path / "m.tsv", ids)
    out = tmp_path / "out"
    # Through the entry point, not run_concordance: figure rendering lives in
    # the CLI wrapper, so calling the module function would test neither.
    assert _cli(["concordance", "--matrix", str(m), "--design", str(design),
                 "--out", str(out)]) == 0
    assert (out / "figures" / "correlation_heatmap.pdf").exists()


def test_no_figures_writes_no_figures(tmp_path):
    ids = ["A1", "A2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    m = _matrix(tmp_path / "m.tsv", ids)
    out = tmp_path / "out"
    assert _cli(["concordance", "--matrix", str(m), "--design", str(design),
                 "--out", str(out), "--no-figures"]) == 0
    assert not (out / "figures").exists()


def test_one_sample_per_group_prints_na_instead_of_crashing(tmp_path):
    """One sample per group: the within-group mean is absent, and is said to be.

    Two samples in different groups form one between-group pair and no
    within-group pair at all, so `within_group_mean_correlation` is None while
    the between-group mean is a real number. Formatting that None with ``:.4f``
    raised TypeError before this case existed; the report must print it as
    absent and must carry the module's own note rather than a clean-looking
    summary that implies replicates were compared.
    """
    ids = ["A1", "B1"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    m = _matrix(tmp_path / "m.tsv", ids)
    out = tmp_path / "out"
    assert run_concordance(_args(matrix=str(m), design=str(design),
                                 out=str(out))) == 0
    payload = _strict_loads((out / "concordance.json").read_text())
    assert payload["summary"]["within_group_mean_correlation"] is None
    assert payload["summary"]["between_group_mean_correlation"] is not None
    # Nothing to compare within a group means the grouping cannot be judged.
    assert payload["summary"]["grouping_supported"] is None
    assert any("at least two samples per group" in n for n in payload["notes"])


# --------------------------------------------------------------- error paths
def test_single_sample_is_refused_with_a_named_reason(tmp_path):
    design = _design(tmp_path / "d.tsv", {"A1": "A"})
    m = _matrix(tmp_path / "m.tsv", ["A1"])
    with pytest.raises(ValueError, match="at least 2 samples"):
        run_concordance(_args(matrix=str(m), design=str(design),
                              out=str(tmp_path / "o")))


def test_min_samples_raises_the_floor(tmp_path):
    ids = ["A1", "A2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    m = _matrix(tmp_path / "m.tsv", ids)
    with pytest.raises(ValueError, match="at least 4 samples"):
        run_concordance(_args(matrix=str(m), design=str(design),
                              min_samples=4, out=str(tmp_path / "o")))


def test_design_and_input_mismatch_names_both_sides(tmp_path):
    design = _design(tmp_path / "d.tsv", {"A1": "A", "A9": "A"})
    m = _matrix(tmp_path / "m.tsv", ["A1", "A2"])
    with pytest.raises(ValueError, match="missing=\\['A2'\\].*extra=\\['A9'\\]"):
        run_concordance(_args(matrix=str(m), design=str(design),
                              out=str(tmp_path / "o")))


def test_design_without_group_column_is_refused(tmp_path):
    d = tmp_path / "d.tsv"
    d.write_text("sample_id\tbatch\nA1\tx\n", encoding="utf-8")
    m = _matrix(tmp_path / "m.tsv", ["A1"])
    with pytest.raises(ValueError, match="group or condition"):
        run_concordance(_args(matrix=str(m), design=str(d),
                              out=str(tmp_path / "o")))


def test_missing_design_is_refused(tmp_path):
    m = _matrix(tmp_path / "m.tsv", ["A1", "A2"])
    with pytest.raises(ValueError, match="--design is required"):
        run_concordance(_args(matrix=str(m), design=None,
                              out=str(tmp_path / "o")))


def test_both_inputs_at_once_is_refused(tmp_path):
    ids = ["A1", "A2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    m = _matrix(tmp_path / "m.tsv", ids)
    p = _peaks(tmp_path / "A1.bed", [("chr1", 1, 2)])
    with pytest.raises(ValueError, match="mutually exclusive"):
        run_concordance(_args(matrix=str(m), peaks=[str(p)],
                              design=str(design), out=str(tmp_path / "o")))


def test_no_input_is_refused(tmp_path):
    design = _design(tmp_path / "d.tsv", {"A1": "A"})
    with pytest.raises(ValueError, match="no input"):
        run_concordance(_args(design=str(design), out=str(tmp_path / "o")))


def test_duplicate_sample_in_one_peaks_list_is_refused(tmp_path):
    ids = ["A1", "A2"]
    design = _design(tmp_path / "d.tsv", {s: s[0] for s in ids})
    p = _peaks(tmp_path / "A1.bed", [("chr1", 1, 2)])
    q = _peaks(tmp_path / "A2.bed", [("chr1", 3, 4)])
    with pytest.raises(ValueError, match="more than once"):
        run_concordance(_args(peaks=[str(p), str(p), str(q)],
                              design=str(design), out=str(tmp_path / "o")))


# ------------------------------------------------------------------ ordering
def test_samples_are_reordered_to_the_design(tmp_path):
    """The heatmap reads group by group, so the output order follows the design."""
    design_groups = {"B1": "B", "A1": "A", "B2": "B", "A2": "A"}
    design = _design(tmp_path / "d.tsv", design_groups)
    m = _matrix(tmp_path / "m.tsv", ["A1", "A2", "B1", "B2"])
    out = tmp_path / "out"
    run_concordance(_args(matrix=str(m), design=str(design), out=str(out)))
    payload = _strict_loads((out / "concordance.json").read_text())
    assert payload["samples"] == ["B1", "A1", "B2", "A2"]
    corr = np.asarray(payload["correlation_matrix"], dtype=float)
    # The matrix must be permuted consistently, so it stays symmetric.
    assert np.allclose(corr, corr.T, equal_nan=True)
    assert np.allclose(np.diag(corr), 1.0)
