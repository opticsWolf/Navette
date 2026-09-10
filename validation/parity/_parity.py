# -*- coding: utf-8 -*-
"""Shared helpers for the parity scripts under ``validation/parity/``.

These files are dual-shaped on purpose: they run as plain scripts
(``python validation/parity/smatrix/test_w_function.py``) and are collected
by pytest. The helpers here make both behave sensibly, and exist because of
two defects found when the scripts were finally collected (R2.4):

``console_utf8()``
    The scripts print ``Delta``/``mu``/``->`` and died with
    ``UnicodeEncodeError`` on a cp1252 console (the same §17 defect fixed for
    the benches in R2.2).

``require_reference()``
    The numba reference was imported from a ``loom/`` directory that does not
    exist, so it never loaded -- and the scripts responded by scoring every
    comparison ``"PASS (rust-only)"``. A parity test that compares nothing and
    passes is worse than no test: it reports coverage it does not have. A
    missing reference is now a **skip**, never a pass.

``report()``
    The scripts printed ``OUTPUT_STATUS FAIL`` and exited **0**, so any caller
    -- a CI step, a shell loop -- saw success. Failure now fails.
"""

from __future__ import annotations

import sys


def console_utf8() -> None:
    """Make stdout/stderr UTF-8-safe. Idempotent; never raises."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            pass


def require_reference(available: bool, what: str, how: str = "") -> None:
    """Skip the module unless the comparison reference is loaded.

    Under pytest this is a module-level skip with a reason; run as a script it
    exits 0 with a SKIP line (nothing failed -- the environment is incomplete).
    Either way the module does not proceed to compare Rust against itself.
    """
    if available:
        return
    reason = f"{what} unavailable, nothing to compare against"
    if how:
        reason += f" ({how})"
    if "pytest" in sys.modules:
        import pytest

        pytest.skip(reason, allow_module_level=True)
    print(f"SKIP: {reason}")
    raise SystemExit(0)


def report(unit: str, all_pass: bool) -> bool:
    """Print the verdict line. Returns `all_pass` so callers can exit on it."""
    print(f"\n{unit}: {'ALL CORRECTNESS TESTS PASSED' if all_pass else 'SOME CORRECTNESS TESTS FAILED'}")
    print(f"OUTPUT_STATUS {'PASS' if all_pass else 'FAIL'}")
    return all_pass
