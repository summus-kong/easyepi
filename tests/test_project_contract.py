from pathlib import Path

from easyepi.cli import build_parser
from easyepi.project import load_project, validate_project
from easyepi.project_commands import make_qc_report_command, validate_project_command
from easyepi.assay_commands import run_assay_qc


def test_project_toml_resolution_and_cli_entries(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    project = tmp_path / "project.toml"
    project.write_text(
        "[project]\nproject_id='p'\n[references]\ngenome='ref/genome.fa'\n"
        "[[samples]]\nsample_id='S1'\nassay='ATAC-seq'\nbigwig='s1.bw'\n",
        encoding="utf-8",
    )
    cfg = load_project(project)
    assert cfg.genome == str((tmp_path / "ref" / "genome.fa").resolve())
    assert cfg.samples[0].sample_id == "S1"
    assert validate_project(cfg, check_files=False) == []
    parser = build_parser()
    assert parser.parse_args(["validate-project", "--project", str(project), "--out", "out"]).module == "validate-project"
    assert parser.parse_args(["make-qc-report", "--project", str(project), "--out", "out"]).module == "make-qc-report"


def test_validate_and_qc_write_run_metadata(tmp_path):
    project = tmp_path / "project.toml"
    project.write_text("[[samples]]\nsample_id='S1'\nassay='ATAC-seq'\n", encoding="utf-8")
    out = tmp_path / "out"
    args = build_parser().parse_args(["validate-project", "--project", str(project), "--out", str(out), "--no-check-files", "--verbose"])
    assert validate_project_command(args) == 0
    assert (out / "resolved_config.json").exists()
    assert (out / "inputs.tsv").exists()
    assert (out / "status.json").exists()
    qc_out = tmp_path / "qc"
    args = build_parser().parse_args(["make-qc-report", "--project", str(project), "--out", str(qc_out), "--no-check-files"])
    assert make_qc_report_command(args) == 0
    assert (qc_out / "qc_report.json").exists()
    assert '"status": "completed"' in (qc_out / "status.json").read_text(encoding="utf-8")


def test_assay_qc_works_without_groups_or_controls(tmp_path):
    peak = tmp_path / "sample.bed"
    peak.write_text("chr1\t10\t20\nchr1\t30\t50\n", encoding="utf-8")
    out = tmp_path / "atac"
    args = build_parser().parse_args(["atac-qc", "--peak", str(peak), "--out", str(out)])
    assert run_assay_qc(args, "ATAC-seq") == 0
    assert (out / "atac_seq_qc.tsv").exists()
    assert (out / "atac_seq_qc.json").exists()


def test_python_script_tree_is_split_by_scope():
    from easyepi.universal import validate_project as universal_validate
    from easyepi.specific.atac import qc as atac_qc
    from easyepi.specific.methylation import qc as methylation_qc
    assert callable(universal_validate.main)
    assert callable(atac_qc.main)
    assert callable(methylation_qc.main)


def test_methylation_python_api_and_script_contract(tmp_path):
    table = tmp_path / "methylation.tsv"
    table.write_text("chrom\tstart\tbeta\tcoverage\nchr1\t10\t0.2\t10\nchr1\t20\t0.8\t12\n", encoding="utf-8")
    from easyepi.methylation import read_methylation_table, summarize_methylation
    frame = read_methylation_table(table)
    summary = summarize_methylation(frame, min_coverage=5)
    assert summary["n_retained"] == 2
    from easyepi.specific.methylation import profile
    project = tmp_path / "project.toml"
    project.write_text(f"[[samples]]\nsample_id='M1'\nassay='methylation'\nmethylation='{table.name}'\n", encoding="utf-8")
    assert profile.main(["--project", str(project), "--out", str(tmp_path / "profile"), "--no-check-files"]) == 0
    assert (tmp_path / "profile" / "M1.methylation.tsv").exists()

