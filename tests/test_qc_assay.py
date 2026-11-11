"""`qc_assay` must accept every spelling it advertises, and key its tables safely.

`resolve_assay` accepts 21 spellings while the tables are keyed by five canonical
names. An earlier version returned the canonical name from `resolve_assay` but then
indexed `ASSAY_INPUTS[assay]` directly, so `check_inputs("ATAC", ...)` raised
KeyError -- a crash one call after the value was accepted.
"""
from __future__ import annotations

import pytest

from easyepi.qc_assay import (ASSAY_ALIASES, ASSAY_INPUTS, AssayError,
                              canonical_assay, check_inputs, describe,
                              resolve_assay)

AVAILABLE = {"S1": {"bam": "a.bam"}}


@pytest.mark.parametrize("spelling", sorted(ASSAY_ALIASES))
def test_every_advertised_spelling_resolves_and_keys_tables(spelling):
    canonical = resolve_assay(spelling, AVAILABLE)
    assert canonical == ASSAY_ALIASES[spelling.strip().lower()]
    assert canonical in ASSAY_INPUTS
    # the same value must be usable by the table-indexing function
    assert canonical_assay(spelling) == canonical
    assert isinstance(check_inputs(spelling, AVAILABLE), list)


def test_check_inputs_accepts_the_spellings_resolve_assay_returns():
    for spelling in ("ATAC", "atac-seq", "cut&tag", "cuttag", "wgbs", "bisulfite"):
        assert check_inputs(spelling, AVAILABLE) == check_inputs(
            canonical_assay(spelling), AVAILABLE)


def test_a_bam_is_a_usable_methylation_input():
    # WGBS.md documents `qc --assay methylation --bam ctrl1.bam`; the note must not
    # call that BAM unusable, which it used to ("methylation QC does not use: bam").
    notes = check_inputs("methylation", {"S1": {"bam": "ctrl1.bam"}})
    assert not any("does not use" in n for n in notes), notes
    assert not any("beta metrics are skipped" in n for n in notes), notes


def test_absent_methylation_input_names_both_sources():
    notes = check_inputs("methylation", {"S1": {"bigwig": "x.bw"}})
    joined = " ".join(notes)
    assert "MM/ML" in joined and "XM" in joined


def test_an_unknown_assay_is_still_an_error():
    with pytest.raises(AssayError):
        canonical_assay("rna-seq")
    with pytest.raises(AssayError):
        resolve_assay("rna-seq", AVAILABLE)


def test_inference_refuses_to_guess_between_bam_assays():
    # BAM alone cannot separate ATAC from ChIP from histone, so it must not pick one.
    with pytest.raises(AssayError, match="cannot infer"):
        resolve_assay(None, {"S1": {"bam": "a.bam"}})


def test_describe_lists_every_canonical_assay():
    table = describe()
    for name in set(ASSAY_ALIASES.values()):
        assert name in table
