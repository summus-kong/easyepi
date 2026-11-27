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


def test_dmr_r_is_declared_a_test_bridge():
    """Product code and test scaffolding are separated by a recorded fact.

    `dmr_r.py` launches R, which the project promises no runtime path does. The
    distinction lived only in prose; declaring it in metadata makes a review able
    to check it, and makes it obvious if the module ever gains a runtime importer.
    """
    cfg = _pyproject()
    declared = cfg["tool"]["easyepi"]["classification"]["test_bridge"]
    assert "easyepi.dmr_r" in declared, (
        "dmr_r is not classified as a test bridge; PUBLISHING documents it as one")

    # No runtime module may import it. Match the import statement, not every
    # mention: `dmr_replicate_aware` contains the substring and is real product
    # code, so a naive search reported two false positives.
    import re as _re
    pattern = _re.compile(
        r"^\s*(?:from|import)\s+easyepi\.dmr_r\b|"
        r"^\s*from\s+\.\s*import\s+dmr_r\b",
        _re.M)
    runtime_importers = []
    for path in (REPO / "easyepi").rglob("*.py"):
        if path.name == "dmr_r.py":
            continue
        if pattern.search(path.read_text(encoding="utf-8")):
            runtime_importers.append(str(path.relative_to(REPO)))
    assert not runtime_importers, (
        f"{runtime_importers} import dmr_r at runtime; it is a DSS test bridge, "
        "not a supported analysis path")


def test_all_extra_unions_every_other_extra():
    """`all` promised to union every extra and did not include `hmmlearn`.

    README and PUBLISHING both describe `all` as the union, so a user who
    installed `easyepi[all]` and then ran `easyepi histone-states` got a
    ModuleNotFoundError: the one extra the docs singled out was the one missing
    from the union. Compare requirement *names*, not strings, so a version bump
    is not a failure.
    """
    extras = _pyproject()["project"]["optional-dependencies"]
    assert "all" in extras, "no `all` extra declared"

    def names(requirements):
        return {
            re.split(r"[<>=!~\[; ]", req, maxsplit=1)[0].strip().lower()
            for req in requirements
        }

    union = set()
    for extra, requirements in extras.items():
        if extra in {"all", "dev"}:
            continue
        union |= names(requirements)
    missing = sorted(union - names(extras["all"]))
    assert not missing, (
        f"`all` does not include {missing}; docs describe it as the union of every "
        "extra, so an install that follows them silently lacks a dependency")


def test_ci_installs_every_optional_extra():
    """The workflow claims to install every extra, so it must actually do so.

    `tobias` was added to pyproject in 0.7.0 and never reached EXTRAS. The TOBIAS
    port tests guard on their imports, so they skipped, and the workflow's own
    "no test may be silently skipped" step is what caught it — after the comment
    beside the list had already asserted the opposite.

    The comparison is by *package name*, not by extra label: `bigwig` and `stats`
    are aliases that `bam` and `viz` already satisfy (pyBigWig, scikit-learn), so
    requiring the label would demand a redundant install while proving nothing
    about what the runner actually has.
    """
    workflow = (REPO / ".github" / "workflows" / "linux.yml").read_text(encoding="utf-8")
    match = re.search(r'EXTRAS:\s*"([^"]+)"', workflow)
    assert match, "no EXTRAS list found in linux.yml"

    extras = _pyproject()["project"]["optional-dependencies"]

    def names(requirements):
        return {
            re.split(r"[<>=!~\[; ]", req, maxsplit=1)[0].strip().lower()
            for req in requirements
        }

    installed_labels = {name.strip() for name in match.group(1).split(",") if name.strip()}
    installed = set()
    for label in installed_labels:
        assert label in extras, f"linux.yml installs unknown extra {label!r}"
        installed |= names(extras[label])

    # `all` is a union, so it would satisfy everything and hide a gap.
    required = set()
    for label, requirements in extras.items():
        if label == "all":
            continue
        required |= names(requirements)

    missing = sorted(required - installed)
    assert not missing, (
        f"linux.yml installs no extra providing {missing}; tests that import them "
        "would skip, and the workflow asserts that no test skips")


def test_dmr_r_has_no_workstation_path():
    """A packaged module must not hardcode one lab's filesystem.

    `easyepi/dmr_r.py` is a test bridge to DSS, not a product code path, and it
    searched two absolute paths under a specific account before falling back to
    PATH. Those are host-specific guesses compiled into a wheel that anyone may
    install. The module now reads `EASYEPI_RSCRIPT` and falls back to PATH only.
    """
    src = (REPO / "easyepi" / "dmr_r.py").read_text(encoding="utf-8")
    offenders = [
        line.strip()
        for line in src.splitlines()
        if ("/data1/" in line or "/home/" in line or "GaoLab" in line)
    ]
    assert not offenders, f"workstation-specific paths in dmr_r.py: {offenders}"


def test_the_changelog_has_an_entry_for_the_current_version():
    declared = _pyproject()["project"]["version"]
    changelog = (REPO / "CHANGELOG.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", changelog, re.M)
    assert any(h.startswith(declared) for h in headings), (
        f"no '## {declared}' heading in CHANGELOG.md; headings are {headings[:4]}")
