"""Shared pytest configuration.

Matplotlib defaults to an interactive backend, which fails on a headless host
("invalid command name / tcl_findLibrary") and would make the plotting tests
depend on the machine rather than on the code. Selecting Agg before any test
imports pyplot keeps the suite runnable on CI and on a server without a display.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
