# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Shared preamble for every bench under ``validation/benches``.

Two jobs, both of which exist because of a failure that actually happened:

``setup_bench()``
    Makes ``stdout``/``stderr`` UTF-8 and puts the repo's ``src/`` on
    ``sys.path``. Several benches print ``x``, ``->`` or emoji and died with
    ``UnicodeEncodeError`` on a cp1252 console; others did
    ``sys.path.insert(0, "src")``, which only works when the bench is launched
    from the repo root.

``require_release()``
    Refuses to benchmark a debug build. The dev venv once shipped a plain
    ``maturin develop`` (debug) extension for long enough that a full round of
    committed timings, and the review conclusions drawn from them, had to be
    discarded -- a debug build is several times slower but looks perfectly
    healthy from Python. Nothing detects that by eye, so the benches assert it.

Usage, at the very top of a bench, before importing ``navette``::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from _bench_common import require_release, setup_bench

    setup_bench()
    require_release()

The check belongs here and not in the pytest suite on purpose: a contributor
running the tests against a debug build is fine, a contributor *timing* one is
not.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable, Optional

_REPO_ROOT = Path(__file__).resolve().parents[2]

_BUILD_HINT = "maturin develop --release   # from the repo root"


def setup_bench() -> None:
    """Make console output UTF-8-safe and ``navette`` importable.

    Idempotent, and never raises: a stream that cannot be reconfigured
    (already-detached, or a plain pipe on some hosts) is left alone rather
    than taking the bench down with it.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            pass

    src = str(_REPO_ROOT / "src")
    if src not in sys.path:
        sys.path.insert(0, src)


def build_profile() -> str:
    """Cargo profile of the installed extension: ``"release"`` or ``"debug"``.

    Raises ``ImportError`` if the extension is not built.
    """
    from navette._navette import build_profile as _native_build_profile

    return str(_native_build_profile())


def require_release(profile_fn: Optional[Callable[[], str]] = None) -> str:
    """Exit the process unless the native extension is an optimized build.

    Returns the profile string on success, so a bench can stamp it into its
    results JSON (see ``bench_provenance``).

    ``profile_fn`` is an injection point for tests; benches call this with no
    arguments.
    """
    try:
        profile = (profile_fn or build_profile)()
    except ImportError as exc:
        sys.exit(
            f"REFUSING TO BENCHMARK: navette's native extension is not "
            f"importable ({exc}). Build it with:\n    {_BUILD_HINT}"
        )
    except AttributeError:
        sys.exit(
            "REFUSING TO BENCHMARK: the installed extension predates "
            "build_profile() and its Cargo profile cannot be determined. "
            f"Rebuild with:\n    {_BUILD_HINT}"
        )

    if profile != "release":
        sys.exit(
            f"REFUSING TO BENCHMARK: extension built with '{profile}' "
            f"profile -- every timing from it is meaningless. Rebuild "
            f"with:\n    {_BUILD_HINT}"
        )
    return profile


def bench_provenance() -> dict:
    """Provenance block to embed in committed bench results.

    Committed JSON from before this existed is of *unknown* profile and should
    not be compared against new runs.
    """
    import platform

    from navette import __version__

    return {
        "build_profile": build_profile(),
        "navette_version": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
