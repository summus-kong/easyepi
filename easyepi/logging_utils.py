"""Shared console panels, optional detailed logs, status records and command provenance."""

from __future__ import annotations

import json
import logging
import shlex
import subprocess
import sys
import traceback
import time
import platform
from importlib.metadata import version, PackageNotFoundError

from easyepi import provenance
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


class PanelFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        """Format one record as ``[EASYEPI][MODULE][LEVEL] message``."""
        module = getattr(record, "module_name", record.name)
        return f"[Easyepi][{module}][{record.levelname}] {record.getMessage()}"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def setup_logging(module: str, out_dir: str | Path, *, verbose: bool = False,
                  quiet: bool = False, log_file: str | Path | None = None) -> logging.Logger:
    """Create console panel logs and optional detailed file logs."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"easyepi.{module.lower().replace(' ', '_')}")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    if not quiet:
        panel = logging.StreamHandler(sys.stdout)
        panel.setLevel(logging.DEBUG if verbose else logging.INFO)
        panel.setFormatter(PanelFormatter())
        logger.addHandler(panel)
    if verbose or log_file:
        detailed = Path(log_file) if log_file else out / "logs" / "run.log"
        detailed.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(detailed, encoding="utf-8")
        handler.setLevel(logging.DEBUG)
        logging.Formatter.converter = time.gmtime
        handler.setFormatter(logging.Formatter("%(asctime)sZ\t%(levelname)s\t%(name)s\t%(message)s", "%Y-%m-%dT%H:%M:%S"))
        logger.addHandler(handler)
    return logging.LoggerAdapter(logger, {"module_name": module})


class RunRecorder:
    """Record run state, metrics, software versions and commands."""

    def __init__(self, module: str, out_dir: str | Path, logger: logging.LoggerAdapter,
                 *, figure_formats: Sequence[str] = (), render_figures: bool = True):
        """Store the output directory, module label and figure formats for this run."""
        self.module = module
        self.out = Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.logger = logger
        # `run_external` receives this same adapter; the back-reference is how
        # an external tool version reaches the manifest without threading it
        # through every call site.
        try:
            setattr(logger, "recorder", self)
        except AttributeError:  # a bare Logger refuses attributes; acceptable
            pass
        self.figure_formats = tuple(figure_formats)
        self.render_figures = render_figures
        self.started = _now()
        self.metrics: dict[str, Any] = {}
        self.status: dict[str, Any] = {"module": module, "status": "running", "started_at": self.started}
        # Filled on demand by `note_external`; see that method.
        self.external_tools: dict[str, dict[str, Any]] = {}
        self._hash_outputs = False

    def write_command(self, argv: Sequence[str] | None = None) -> Path:
        """Record the exact command line for reproducibility."""
        path = self.out / "command.txt"
        path.write_text(" ".join(shlex.quote(str(x)) for x in (argv or sys.argv)), encoding="utf-8")
        return path

    def note_external(self, executable: str, version_args: Sequence[str] = ("--version",)) -> None:
        """Record one external tool's version banner the first time it runs.

        Called from `run_external`, so only tools a run really used appear in the
        manifest. Repeat calls for the same executable are ignored: the first
        banner is the one the run actually ran under.
        """
        name = str(executable)
        if name in self.external_tools:
            return
        self.external_tools[name] = provenance.external_tools({name: version_args})[name]

    def hash_outputs(self) -> None:
        """Ask `finish` to hash whatever `--out` contains.

        Off by default because a caller hashing a multi-gigabyte BAM pays real
        IO for a field it may never read. Opt in at the call site.
        """
        self._hash_outputs = True

    def set_metrics(self, **values: Any) -> None:
        """Attach summary metrics written to ``metrics.json`` on finish."""
        self.metrics.update(values)

    def finish(self, *, status: str = "completed", error: BaseException | None = None) -> None:
        """Write ``status.json``, ``metrics.json`` and any requested figures."""
        if error is None and self.render_figures and self.figure_formats:
            try:
                from easyepi._plots import render_module_figures
                figures = render_module_figures(self.module, self.out, self.figure_formats)
                if figures:
                    self.metrics.setdefault("figures", figures)
                    self.logger.info("Wrote %d figure file(s)", len(figures))
            except Exception as exc:  # Plot export failures must not hide completed analysis results.
                self.logger.warning("Figure rendering skipped: %s", exc)
        record = provenance.build(hash_inputs=False)
        if self.external_tools:
            record["external_tools"] = self.external_tools
        if self._hash_outputs:
            produced = sorted(
                p for p in self.out.rglob("*") if p.is_file() and p.name != "status.json"
            )
            record["outputs"] = provenance.digest_files(produced, root=self.out)
        self.status["provenance"] = record
        gaps = provenance.describe_gaps(record)
        for gap in gaps:
            self.logger.warning("provenance gap: %s", gap)
        self.status.update({"status": status, "finished_at": _now(), "metrics": self.metrics})
        if error:
            self.status["error"] = {"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()}
        (self.out / "status.json").write_text(json.dumps(self.status, indent=2, ensure_ascii=False), encoding="utf-8")
        (self.out / "metrics.json").write_text(json.dumps(self.metrics, indent=2, ensure_ascii=False), encoding="utf-8")

    def __enter__(self) -> "RunRecorder":
        self.write_command()
        self.status["versions"] = provenance.package_versions()
        (self.out / "status.json").write_text(json.dumps(self.status, indent=2), encoding="utf-8")
        self.logger.info("Starting %s", self.module)
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc is None:
            self.finish()
            self.logger.info("DONE")
        else:
            self.finish(status="failed", error=exc)
            self.logger.error("FAILED: %s", exc)
        return False


def run_external(command: Sequence[str], logger: logging.LoggerAdapter, *, cwd: str | Path | None = None,
                 check: bool = True) -> subprocess.CompletedProcess[str]:
    """Execute an external command and record stdout, stderr and failure status."""
    logger.info("External command: %s", shlex.join(map(str, command)))
    # `RunRecorder` attaches itself to its LoggerAdapter, so a run that
    # shells out can record the tool without every call site passing it.
    recorder = getattr(logger, "recorder", None)
    if recorder is not None and command:
        try:
            recorder.note_external(str(command[0]))
        except Exception:  # a version probe must never break the real command
            pass
    result = subprocess.run(list(map(str, command)), cwd=cwd, text=True, capture_output=True)
    if result.stdout:
        logger.debug("stdout: %s", result.stdout.rstrip())
    if result.stderr:
        logger.debug("stderr: %s", result.stderr.rstrip())
    logger.info("External command exit code: %s", result.returncode)
    if check and result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command, result.stdout, result.stderr)
    return result

