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
        self.figure_formats = tuple(figure_formats)
        self.render_figures = render_figures
        self.started = _now()
        self.metrics: dict[str, Any] = {}
        self.status: dict[str, Any] = {"module": module, "status": "running", "started_at": self.started}

    def write_command(self, argv: Sequence[str] | None = None) -> Path:
        """Record the exact command line for reproducibility."""
        path = self.out / "command.txt"
        path.write_text(" ".join(shlex.quote(str(x)) for x in (argv or sys.argv)), encoding="utf-8")
        return path

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
        self.status.update({"status": status, "finished_at": _now(), "metrics": self.metrics})
        if error:
            self.status["error"] = {"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()}
        (self.out / "status.json").write_text(json.dumps(self.status, indent=2, ensure_ascii=False), encoding="utf-8")
        (self.out / "metrics.json").write_text(json.dumps(self.metrics, indent=2, ensure_ascii=False), encoding="utf-8")

    def __enter__(self) -> "RunRecorder":
        self.write_command()
        versions = {"python": platform.python_version()}
        for name in ("easyepi", "numpy", "scipy", "pandas", "pysam", "pyBigWig"):
            try:
                versions[name] = version(name)
            except PackageNotFoundError:
                versions[name] = "not installed"
        self.status["versions"] = versions
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
    result = subprocess.run(list(map(str, command)), cwd=cwd, text=True, capture_output=True)
    if result.stdout:
        logger.debug("stdout: %s", result.stdout.rstrip())
    if result.stderr:
        logger.debug("stderr: %s", result.stderr.rstrip())
    logger.info("External command exit code: %s", result.returncode)
    if check and result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command, result.stdout, result.stderr)
    return result

