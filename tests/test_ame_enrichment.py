"""The AME wrapper: parsing, parameter plumbing, and refusal of useless tests.

The parsing tests matter because a column-index mistake in the FIMO path earlier in
this project silently produced zero hits instead of failing. Anything that reads a
tool's tabular output gets a test that asserts the values land in the right fields.
"""
import shutil
import textwrap

import pytest

from easyepi.motifs.ame_enrichment import (
    ame_enrichment,
    parse_ame_tsv,
)

AME = shutil.which("ame")

TSV = textwrap.dedent("""\
    rank\tmotif_DB\tmotif_ID\tmotif_alt_ID\tconsensus\tp-value\tadj_p-value\tE-value\ttests\tFASTA_max\tpos\tneg\tPWM_min\tTP\t%TP\tFP\t%FP
    1\tlib.meme\tHEB(bHLH)\tHEB(bHLH)\tVCAGSTGBNN\t1.06e-26\t3.33e-24\t1.50e-21\t313\t300\t300\t300\t1\t221\t73.67\t92\t30.67
    2\tlib.meme\tKLF5(Zf)\tKLF5(Zf)\tDGGGYGKGGC\t1.91e-21\t3.65e-19\t1.64e-16\t191\t300\t300\t300\t1.68\t124\t41.33\t26\t8.67
""")


class TestParse:
    def test_fields_land_in_the_right_columns(self):
        rows = parse_ame_tsv(TSV)
        assert len(rows) == 2
        first = rows[0]
        # sequence name / motif id confusion is exactly the bug to guard against
        assert first["motif_ID"] == "HEB(bHLH)"
        assert first["consensus"] == "VCAGSTGBNN"
        assert first["TP"] == 221
        assert first["FP"] == 92
        assert abs(first["%TP"] - 73.67) < 1e-9
        assert abs(first["%FP"] - 30.67) < 1e-9
        assert abs(first["p-value"] - 1.06e-26) < 1e-30
        assert first["rank"] == 1

    def test_rank_is_int_not_string(self):
        rows = parse_ame_tsv(TSV)
        assert isinstance(rows[0]["rank"], int)
        assert isinstance(rows[0]["tests"], int)

    def test_short_rows_are_skipped_not_misaligned(self):
        broken = TSV + "3\tlib.meme\tTRUNCATED\n"
        rows = parse_ame_tsv(broken)
        assert len(rows) == 2, "a truncated row must not become a dict with holes"

    def test_empty_input_gives_no_rows(self):
        assert parse_ame_tsv("") == []
        assert parse_ame_tsv("\n\n") == []

    def test_header_only_gives_no_rows(self):
        header = TSV.splitlines()[0]
        assert parse_ame_tsv(header) == []


class TestRefusals:
    def test_no_control_and_no_shuffle_is_refused(self):
        with pytest.raises(ValueError, match="null"):
            ame_enrichment(["ACGT" * 10], None, "db.meme", shuffle_control=False)

    def test_empty_targets_refused(self):
        with pytest.raises(ValueError, match="target"):
            ame_enrichment([], ["ACGT" * 10], "db.meme")

    def test_missing_executable_names_the_install(self):
        with pytest.raises(RuntimeError, match="MEME Suite"):
            ame_enrichment(["ACGT" * 10], ["ACGT" * 10], "db.meme",
                           ame_executable="/nonexistent/ame")


@pytest.mark.skipif(AME is None, reason="AME not installed")
class TestAgainstAme:
    @pytest.fixture(scope="class")
    def motif_db(self, tmp_path_factory):
        from easyepi.motifs.formats import Motif, write_meme
        import numpy as np
        d = tmp_path_factory.mktemp("db")
        motifs = []
        for name, core in (("M1", "ACGTACGT"), ("M2", "TTTTAAAA")):
            p = np.full((4, len(core)), 0.1)
            for j, ch in enumerate(core):
                p["ACGT".index(ch), j] = 0.7
            m = Motif(name=name, matrix=p)
            motifs.append(m)
        out = d / "lib.meme"
        write_meme(motifs, out)
        return str(out)

    def test_runs_and_reports_its_invocation(self, motif_db):
        target = ["ACGTACGT" + "AAAA" * 20 + "ACGTACGT",
                  "ACGTACGT" + "CCCC" * 20 + "ACGTACGT",
                  "GGGG" * 10 + "ACGTACGT" + "TTTT" * 10]
        control = ["GGGGCCCC" * 15, "CCCCGGGG" * 15, "ATATATAT" * 15]
        res = ame_enrichment(target, control, motif_db)
        assert res["backend"] == "ame"
        assert res["control"] == "supplied"
        assert res["n_target"] == 3
        assert res["n_control"] == 3
        assert "ame" in res["command_shape"]
        assert res["ame_version"]
        assert isinstance(res["motifs"], list)

    def test_shuffle_control_needs_no_background(self, motif_db):
        target = ["ACGTACGT" + "AAAA" * 20 + "ACGTACGT",
                  "ACGTACGT" + "CCCC" * 20 + "ACGTACGT"]
        res = ame_enrichment(target, None, motif_db, shuffle_control=True, kmer=2)
        assert res["control"].startswith("shuffle")
        assert res["n_control"] == 0

    def test_method_and_scoring_are_forwarded(self, motif_db):
        target = ["ACGTACGT" + "AAAA" * 20, "ACGTACGT" + "CCCC" * 20]
        control = ["GGGGCCCC" * 15, "CCCCGGGG" * 15]
        res = ame_enrichment(target, control, motif_db,
                             method="ranksum", scoring="totalhits")
        assert res["method"] == "ranksum"
        assert res["scoring"] == "totalhits"
        assert "--method ranksum" in res["command_shape"]
        assert "--scoring totalhits" in res["command_shape"]
