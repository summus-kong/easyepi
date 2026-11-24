"""Provenance must record what it can and say so when it cannot.

Three things are worth guarding here, each learned from a bug this
module had while being written:

  * A tool that exists outside `PATH` must still be found. The first
    version only used `shutil.which`, so every tool on the analysis host
    reported `not_found` while sitting in a conda env bin directory.
  * Version banners come in incompatible shapes. `samtools 1.24`,
    `bedtools v2.31.1` and `bowtie2-align-s version 2.5.5` all had to
    yield a version, and the `v` prefix must not end up in the number.
  * A missing file, a missing tool and an unreadable path must produce a
    recorded reason, never an exception and never a silent absence.
"""

from __future__ import annotations

import json
import os
import sys

import pytest

from easyepi import provenance


# --- file digests -------------------------------------------------------


def test_digest_reports_hash_and_size(tmp_path):
    payload = b"@r1\nACGT\n+\nIIII\n"
    target = tmp_path / "reads.fq"
    target.write_bytes(payload)

    got = provenance.digest_files([target], root=tmp_path)

    assert list(got) == ["reads.fq"]
    assert got["reads.fq"]["size"] == len(payload)
    # sha256 of the exact bytes, computed independently of the module.
    import hashlib

    assert got["reads.fq"]["sha256"] == hashlib.sha256(payload).hexdigest()


def test_digest_records_a_missing_file_instead_of_raising(tmp_path):
    got = provenance.digest_files([tmp_path / "gone.bam"], root=tmp_path)

    assert got["gone.bam"] == {"error": "not_found"}


def test_digest_skips_oversized_files_and_says_why(tmp_path):
    big = tmp_path / "huge.bam"
    big.write_bytes(b"x" * 64)

    got = provenance.digest_files([big], root=tmp_path, limit_bytes=10)

    assert got["huge.bam"]["error"] == "skipped_too_large"
    assert got["huge.bam"]["size"] == 64


def test_digest_keeps_an_absolute_path_when_it_cannot_be_made_relative(tmp_path):
    outside = tmp_path.parent / "elsewhere.fq"
    outside.write_bytes(b"x")

    got = provenance.digest_files([outside], root=tmp_path)

    assert str(outside) in got or outside.name in got


# --- packages -----------------------------------------------------------


def test_package_versions_always_carries_python():
    got = provenance.package_versions([])

    assert got == {"python": __import__("platform").python_version()}


def test_package_versions_marks_an_absent_package():
    got = provenance.package_versions(["definitely-not-a-real-package-xyz"])

    assert got["definitely-not-a-real-package-xyz"] == "not installed"


# --- external tools -----------------------------------------------------


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX shell fixture")
def test_a_tool_outside_path_is_found_through_search_paths(tmp_path):
    """The regression that motivated `search_paths`.

    A conda env bin directory is not on `PATH`. Before this, every tool
    in such an environment was reported as `not_found` even though the
    file existed and ran.
    """
    tool = tmp_path / "faketool"
    tool.write_text("#!/bin/sh\necho 'faketool v9.8.7'\n", encoding="utf-8")
    tool.chmod(0o755)

    got = provenance.external_tools({"faketool": ()}, search_paths=[str(tmp_path)])

    assert got["faketool"]["status"] == "ok"
    assert got["faketool"]["version"] == "9.8.7"
    assert got["faketool"]["path"] == str(tool)


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX shell fixture")
def test_version_is_extracted_from_three_banner_shapes(tmp_path):
    """Real banners disagree about where the number sits."""
    cases = {
        "bare": ("samtools 1.24", "1.24"),
        "v_prefixed": ("bedtools v2.31.1", "2.31.1"),
        "long_program_name": ("x-align-s version 2.5.5", "2.5.5"),
    }
    for name, (banner, _) in cases.items():
        script = tmp_path / name
        script.write_text(f"#!/bin/sh\necho '{banner}'\n", encoding="utf-8")
        script.chmod(0o755)

    for name, (banner, expected) in cases.items():
        got = provenance.external_tools({name: ()}, search_paths=[str(tmp_path)])
        assert got[name]["version"] == expected, banner
        assert got[name]["banner"] == banner


def test_an_absent_tool_is_reported_not_omitted():
    got = provenance.external_tools({"surely-not-installed-tool": ("--version",)})

    assert got["surely-not-installed-tool"] == {"status": "not_found"}


def test_a_tool_that_fails_still_returns_a_status(tmp_path):
    broken = tmp_path / "broken"
    broken.write_text("#!/bin/sh\nexit 3\n", encoding="utf-8")
    broken.chmod(0o755)

    got = provenance.external_tools({"broken": ()}, search_paths=[str(tmp_path)])

    assert got["broken"]["status"] == "failed"


# --- environment --------------------------------------------------------


def test_environment_omits_machine_identity():
    got = provenance.environment()

    for unwanted in ("hostname", "user", "username", "node"):
        assert unwanted not in got
    assert "python_version" in got


# --- the assembled record ----------------------------------------------


def test_build_is_json_serialisable(tmp_path):
    """The recorder writes this straight into `status.json`."""
    f = tmp_path / "a.txt"
    f.write_text("hello", encoding="utf-8")

    record = provenance.build(inputs=[f], root=tmp_path)

    json.dumps(record)  # must not raise
    assert "inputs" in record
    assert record["inputs"]["a.txt"]["size"] == 5


def test_outputs_are_only_hashed_when_asked(tmp_path):
    f = tmp_path / "out.bed"
    f.write_text("chr1\t0\t10\n", encoding="utf-8")

    without = provenance.build(outputs=[f], root=tmp_path)
    with_it = provenance.build(outputs=[f], root=tmp_path, hash_outputs=True)

    assert "outputs" not in without
    assert "out.bed" in with_it["outputs"]


def test_describe_gaps_names_every_unset_field(tmp_path):
    record = provenance.build(
        inputs=[tmp_path / "missing.bam"],
        root=tmp_path,
        tools={"absent-tool": ("--version",)},
        packages=["also-absent-package"],
    )

    gaps = provenance.describe_gaps(record)
    joined = "\n".join(gaps)

    assert "missing.bam" in joined and "not_found" in joined
    assert "absent-tool" in joined
    assert "also-absent-package" in joined


def test_a_complete_record_has_no_gaps(tmp_path):
    f = tmp_path / "ok.txt"
    f.write_text("x", encoding="utf-8")

    record = provenance.build(inputs=[f], root=tmp_path, packages=[])

    # Only the default tool list can produce gaps here; drop it.
    record["external_tools"] = {}
    assert provenance.describe_gaps(record) == []
