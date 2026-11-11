"""Release metadata must stay consistent, and the wheel must contain the package.

`c27cd27` removed a `License :: OSI Approved :: BSD License` classifier because it
made every install fail, and nothing in CI would have noticed: the editable install
tests the source tree, which never reads the packaged metadata. These tests cover the
parts that can be checked without building a wheel; `.github/workflows/linux.yml`
builds and installs one for the rest.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

tomllib = pytest.importorskip("tomllib", reason="Python 3.11+")

REPO = Path(__file__).resolve().parents[1]


def _pyproject() -> dict:
    return tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))


def test_version_agrees_across_pyproject_and_package():
    declared = _pyproject()["project"]["version"]
    src = (REPO / "easyepi" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', src, re.M)
    assert match, "__version__ not found in easyepi/__init__.py"
    assert match.group(1) == declared, (
        f"__version__ {match.group(1)!r} != pyproject {declared!r}; "
        "`easyepi --version` would report something other than the installed one")


def test_no_legacy_license_classifier():
    """The classifier that broke every install must not come back.

    With PEP 639 a project declares `license = "BSD-3-Clause"`, and an added
    `License :: OSI Approved :: ...` classifier conflicts with the license
    expression closely enough that setuptools refuses to build.
    """
    project = _pyproject()["project"]
    classifiers = project.get("classifiers", [])
    legacy = [c for c in classifiers if c.startswith("License ::")]
    assert not legacy, f"legacy License classifiers break installs: {legacy}"
    assert project.get("license"), "project.license is not declared"


def test_required_modules_are_importable_from_the_package():
    """Cheap guard against a packaging `exclude` silently dropping a real module."""
    import importlib

    for name in ("easyepi", "easyepi.cli", "easyepi.assays", "easyepi._plots",
                 "easyepi.qc_assay", "easyepi.peaks_cli", "easyepi.assay_metrics"):
        assert importlib.import_module(name) is not None


def test_packages_find_does_not_exclude_a_real_subpackage():
    """`exclude` names three subpackages that no longer exist; a stale entry that
    matches real code would drop it from the wheel without any test noticing."""
    cfg = _pyproject()
    find = cfg.get("tool", {}).get("setuptools", {}).get("packages", {}).get("find", {})
    excludes = [e for e in find.get("exclude", []) if e.startswith("easyepi")]
    existing = {p.name for p in (REPO / "easyepi").iterdir() if p.is_dir()}
    for pattern in excludes:
        # `easyepi.<x>*` matches nothing today; assert it still matches nothing
        stem = pattern.replace("easyepi.", "").rstrip("*")
        assert stem not in existing, (
            f"packages.find excludes {pattern} but easyepi/{stem} exists; "
            "that subpackage would be missing from the wheel")


def test_the_changelog_has_an_entry_for_the_current_version():
    declared = _pyproject()["project"]["version"]
    changelog = (REPO / "CHANGELOG.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", changelog, re.M)
    assert any(h.startswith(declared) for h in headings), (
        f"no '## {declared}' heading in CHANGELOG.md; headings are {headings[:4]}")
