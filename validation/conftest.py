# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Pytest collection rules for validation/.

Collected:
  smoke/                 - package import surface, no build quirks assumed
  goldens/spectralweave/ - pinned numeric regression tests
  regression/            - differential/bit-identity suites
  parity/                - numba-vs-Rust parity (R2.4; see below)

Ignored:
  benches/               - timing scripts; they hard-exit on a debug build,
                           which is correct for a bench and wrong for a test
                           run. Run them explicitly:
                             python validation/benches/<area>/<script>.py
  parity/**/refs/        - the numba reference implementations themselves
  parity/**/gen_*.py     - golden regenerators (they WRITE files)

Why parity/ is collected now (R2.4, review §15)
-----------------------------------------------
It was blanket-ignored as "standalone scripts", and that hid the strongest
oracle in the repo — plus three defects that only surfaced on collection:

  * the numba reference was imported from a `loom/` directory that does not
    exist, so it never loaded and every comparison scored
    "PASS (rust-only)" — passing while comparing nothing;
  * the scripts printed `OUTPUT_STATUS FAIL` and exited 0, so failure was
    invisible to any caller;
  * two files called `sys.exit(1)` at import, which turns collection into a
    pytest INTERNALERROR.

The dual-shaped scripts keep working as `python <file>`; each now also
exposes a `test_parity()` that asserts the verdict. See `parity/_parity.py`.
"""

collect_ignore = ["benches"]
collect_ignore_glob = [
  "parity/**/refs/*",
  "parity/**/gen_*.py",
  "parity/*/golden_mirror.py",
]

import pytest  # noqa: E402  (fixture support below)


@pytest.fixture
def rng_for():
  """Injectable flip-proof RNG factory (see `_rng_for` below)."""
  return _rng_for


def _rng_for(seed):
  """Flip-proof RNG handle: `Generator` pre-flip, raw seed post-flip.

  All error-path tests take their randomness through this helper so the
  twin files run byte-identical against the Python implementation (which
  consumes `Generator`s) and the bound classes (which consume seeds).
  Post-flip this returns the seed itself; until then, a `Generator`.
  """
  import numpy as np
  try:
    from navette.structure import Layer
    from navette._structure import Layer as RsLayer
    if Layer is RsLayer:
      return seed
  except ImportError:
    pass
  return np.random.default_rng(seed)
