# -*- coding: utf-8 -*-
"""Every file in ``examples/`` must still run (R4.3).

``examples/spectralweave_example.py`` failed on a fresh clone. It called the
raw native ``OpticalWeaver`` with a key tuple rather than the documented
``SimulationWeaver`` / ``OpticalFragment`` API, and handed ``unweave`` a
``linspace`` over the same range as the stored frame -- which shares one
wavelength with it, not a hundred. The engine reported ``Length mismatch``
from four call levels down.

Nothing ran the examples, so the API could drift away from them indefinitely.
This runs each one in a subprocess: pytest's process state, matplotlib
backends and any module-level side effects stay out of the suite, and the
example is exercised the way a reader would actually invoke it.
"""

import subprocess
import sys
from pathlib import Path

import pytest

_EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
_SCRIPTS = sorted(_EXAMPLES.glob("*.py"))


def test_the_examples_directory_is_not_empty():
    """A glob that silently matches nothing would make this file a no-op."""
    assert _SCRIPTS, f"no *.py under {_EXAMPLES}"


@pytest.mark.parametrize("script", _SCRIPTS, ids=lambda p: p.name)
def test_example_runs(script):
    """Run it as a script, from the repo root, and require a clean exit."""
    proc = subprocess.run(
        [sys.executable, str(script)],
        cwd=str(_EXAMPLES.parent),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env={**__import__("os").environ, "PYTHONIOENCODING": "utf-8",
             "MPLBACKEND": "Agg"},
        timeout=120,
    )
    assert proc.returncode == 0, (
        f"examples/{script.name} exited {proc.returncode}\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")
    assert proc.stdout.strip(), (
        f"examples/{script.name} printed nothing; an example that shows the "
        f"reader no output cannot have demonstrated anything")
