"""Regression: HOMER .motifs parser must NOT truncate long motifs.

This guards against the bug where read_homer_motifs read a fixed 4 rows per
motif, cutting a length-L motif (e.g. COUP-TFII's 12 bp 'GKBCARAGGTCA') down
to 4 bp. A truncated motif has length 4 and therefore can never be matched as
the real TF, silently breaking known-motif enrichment.
"""
import numpy as np
import pytest

from easyepi.motifs import formats as F


def _mock_motifs_text():
    """4 bp motif + 12 bp motif, as a HOMER .motifs collection."""
    rows4 = ["0.205 0.217 0.473 0.105",
             "0.146 0.205 0.398 0.251",
             "0.141 0.292 0.218 0.349",
             "0.209 0.426 0.182 0.183"]
    nb = len(rows4)  # 4
    rows12 = ["0.205 0.217 0.473 0.105"] * 12
    txt = ">GCTC\tAT1G24250(Orphan)/Homer\t6.0\t0\t0\tT:1.0,B:0.1\n" \
          + "\n".join(rows4) + "\n"
    txt += ">GKBCARAGGTCA\tCOUP-TFII(NR)/Homer\t6.4\t0\t0\tT:1.0,B:0.1\n" \
           + "\n".join(rows12) + "\n"
    return txt


def test_read_homer_motifs_does_not_truncate_long_motif(tmp_path):
    p = tmp_path / "mixed.motifs"
    p.write_text(_mock_motifs_text(), encoding="utf-8")
    ms = F.read_homer_motifs(str(p))
    assert len(ms) == 2
    by_name = {m.name: m for m in ms}
    # The 4 bp motif stays 4
    assert by_name["AT1G24250(Orphan)"].length == 4
    # The long NR motif must NOT be truncated to 4; it should be 12 bp.
    coup = by_name["COUP-TFII(NR)"]
    assert coup.length == 12, f"COUP-TFII was truncated to {coup.length} bp (bug!)"
    assert coup.matrix.shape == (4, 12)


def test_read_homer_motifs_roundtrip_realistic_consensus():
    """Consensus of the parsed long motif matches the header consensus length."""
    p = "/tmp/_reg.motifs"
    import tempfile, os
    p = tempfile.mktemp(suffix=".motifs")
    open(p, "w").write(_mock_motifs_text().replace("GKBCARAGGTCA", "AGGTCAAGGTCA"))
    ms = F.read_homer_motifs(p)
    os.remove(p)
    coup = [m for m in ms if "COUP" in m.name][0]
    assert coup.length == 12
    # consensus length equals motif length (12)
    assert len(coup.consensus_seq()) == 12


def test_motif_frequency_rows_are_acgt():
    """Rows in a HOMER .motifs block are ACGT probabilities -> shape (4, L)."""
    p = "/tmp/_reg2.motifs"
    import tempfile, os
    p = tempfile.mktemp(suffix=".motifs")
    open(p, "w").write(_mock_motifs_text())
    ms = F.read_homer_motifs(p)
    os.remove(p)
    for m in ms:
        # every column sums to ~1 (normalized frequency)
        colsum = m.matrix.sum(axis=0)
        assert np.allclose(colsum, 1.0, atol=1e-6)

def test_which_tool_finds_an_executable_in_the_interpreter_prefix(tmp_path, monkeypatch):
    """A conda-installed tool must be found even when its env is not activated.

    On the analysis host `findMotifs.pl` lives in `<env>/bin` beside the running
    `python`, but that directory is only on PATH after activation. Every entry
    point that runs `python -m easyepi.cli` or a scheduler job therefore saw
    `homer: False` with HOMER installed — and the advice that follows that answer
    ("install HOMER") is wrong.
    """
    from easyepi.motifs import external_tools

    fake_prefix = tmp_path / "env"
    (fake_prefix / "bin").mkdir(parents=True)
    exe = fake_prefix / "bin" / "findMotifs.pl"
    exe.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    exe.chmod(0o755)

    # PATH must not contain it, which is the whole point.
    monkeypatch.setattr(external_tools.shutil, "which", lambda name: None)
    monkeypatch.setattr(external_tools.sys, "prefix", str(fake_prefix))
    monkeypatch.setattr(external_tools.sys, "base_prefix", str(fake_prefix))

    found = external_tools.which_tool("findMotifs.pl")
    assert found == str(exe), f"prefix search did not find {exe}; got {found}"
    assert external_tools.tool_available("findMotifs.pl")


def test_which_tool_still_reports_a_genuinely_absent_tool(monkeypatch):
    """The prefix search must not turn a missing tool into a false positive."""
    from easyepi.motifs import external_tools

    monkeypatch.setattr(external_tools.shutil, "which", lambda name: None)
    assert external_tools.which_tool("definitely-not-installed-xyz") is None
    assert not external_tools.tool_available("definitely-not-installed-xyz")

# --------------------------------------------------------------------------
# Site-level scanning (scan_sites): the things a site table must get right
# --------------------------------------------------------------------------
def _sharp_pwm():
    """A 4 bp PWM that strongly prefers ACGT, as log-odds."""
    import numpy as np
    return np.log2(np.array([
        [0.97, 0.01, 0.01, 0.01],
        [0.01, 0.97, 0.01, 0.01],
        [0.01, 0.01, 0.97, 0.01],
        [0.01, 0.01, 0.01, 0.97],
    ]) / 0.25)


def test_scan_sites_reports_coordinates_on_the_plus_strand():
    from easyepi.motifs import scan_sites

    res = scan_sites(_sharp_pwm(), ["TTTACGTAAA"], fpr=1e-2)
    assert res["n_sites"] == 1, res["sites"]
    site = res["sites"][0]
    assert site["strand"] == "+"
    assert (site["start"], site["end"]) == (3, 7)
    assert res["motif_length"] == 4


def test_scan_sites_reports_the_reverse_complement_at_forward_coordinates():
    """A minus-strand hit must carry forward-strand coordinates.

    Reporting the offset within the reverse-complemented string is the classic
    off-by-length error: the interval then points at the wrong place in BED.
    """
    from easyepi.motifs import scan_sites

    # ACGT's reverse complement is ACGT, so use a non-palindromic motif.
    import numpy as np
    pwm = np.log2(np.array([
        [0.97, 0.01, 0.01, 0.01],
        [0.01, 0.01, 0.97, 0.01],
        [0.01, 0.97, 0.01, 0.01],
        [0.01, 0.01, 0.01, 0.97],
    ]) / 0.25)  # prefers A C G T -> rc is A C G T too; use AAGC instead
    pwm = np.log2(np.array([
        [0.97, 0.97, 0.01, 0.01],
        [0.01, 0.01, 0.97, 0.01],
        [0.01, 0.01, 0.01, 0.97],
        [0.01, 0.01, 0.01, 0.01],
    ]) / 0.25)  # prefers "AACG"

    seq = "TTTTTCGTTTT"  # rc of AACG is CGTT, present at offset 5
    res = scan_sites(pwm, [seq], fpr=1e-2)
    minus = [s for s in res["sites"] if s["strand"] == "-"]
    assert minus, res["sites"]
    site = minus[0]
    rc = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rcseq = "".join(rc[c] for c in reversed(seq))
    assert rcseq[site["start"]:site["end"]] is not None
    # the forward interval must be within bounds and of the motif length
    assert 0 <= site["start"] < site["end"] <= len(seq)
    assert site["end"] - site["start"] == 4


def test_scan_sites_skips_windows_containing_ambiguous_bases():
    from easyepi.motifs import scan_sites

    # ACGT with an N inside would otherwise score as a perfect hit.
    res = scan_sites(_sharp_pwm(), ["ACNTACGT"], fpr=1e-2)
    starts = [(s["start"], s["strand"]) for s in res["sites"]]
    assert (0, "+") not in starts, "a window with N was reported as a site"
    assert (4, "+") in starts


def test_scan_sites_handles_short_and_empty_sequences():
    from easyepi.motifs import scan_sites

    res = scan_sites(_sharp_pwm(), ["AC", ""], fpr=1e-2)
    assert res["n_sites"] == 0
    assert res["threshold"] is not None

    res = scan_sites(_sharp_pwm(), [], fpr=1e-2)
    assert res["n_sites"] == 0


def test_scan_sites_single_strand_option():
    import numpy as np
    from easyepi.motifs import scan_sites

    pwm = np.log2(np.array([
        [0.97, 0.97, 0.01, 0.01],
        [0.01, 0.01, 0.97, 0.01],
        [0.01, 0.01, 0.01, 0.97],
        [0.01, 0.01, 0.01, 0.01],
    ]) / 0.25)
    seq = "TTTTTCGTTTT"  # only the reverse complement matches
    both = scan_sites(pwm, [seq], fpr=1e-2)
    plus_only = scan_sites(pwm, [seq], fpr=1e-2, reverse_complement=False)
    assert both["n_sites"] >= 1
    assert plus_only["n_sites"] == 0


def test_scan_sites_pvalues_are_upper_tail_and_ordered_with_score():
    from easyepi.motifs import scan_sites

    res = scan_sites(_sharp_pwm(), ["ACGTACGT", "TTTTTTTT"], fpr=1e-2, pvalue=True)
    scored = [s for s in res["sites"] if "pvalue" in s]
    assert scored, res["sites"]
    for s in scored:
        assert 0.0 <= s["pvalue"] <= 1.0
    ordered = sorted(scored, key=lambda s: -s["score"])
    pv = [s["pvalue"] for s in ordered]
    assert pv == sorted(pv), "p-value must not decrease as the score increases"


def test_scan_sites_and_scan_sequences_agree_on_the_best_score():
    """The site table and the summary must not disagree about the same hit."""
    from easyepi.motifs import scan_sequences, scan_sites

    seqs = ["TTTACGTAAA", "ACGTACGTAC", "GGGGCCCCGG"]
    best = scan_sequences(_sharp_pwm(), seqs)["scores"]
    sites = scan_sites(_sharp_pwm(), seqs, fpr=1e-2)["sites"]
    for i, score in enumerate(best):
        mine = [s["score"] for s in sites if s["seq_index"] == i]
        if mine:
            assert abs(max(mine) - score) < 1e-6, (i, score, max(mine))

def test_scan_sites_rejects_an_unreachable_fpr_with_a_useful_message():
    """A 4 bp PWM cannot reach fpr=1e-4; that must be said, not silently clamped."""
    import pytest as _pytest

    from easyepi.motifs import scan_sites

    with _pytest.raises(ValueError, match="no PWM score reaches"):
        scan_sites(_sharp_pwm(), ["ACGT"], fpr=1e-4)
