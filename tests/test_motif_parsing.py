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
