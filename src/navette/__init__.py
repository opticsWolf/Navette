# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Navette: Weaving the mathematics of light in thin film systems.

Unified Python package. Rust-accelerated subpackages are thin wrappers
over private native submodules of the single aggregated extension
``navette._navette`` (built with ``maturin develop --release`` from the
repo root)::

==================  ========================
public wrapper      native submodule
==================  ========================
``navette.color``         ``navette._color``
``navette.interpolate``   ``navette._interpolate``
``navette.smatrix``       ``navette._smatrix``
``navette.spectralweave`` ``navette._spectralweave``
``navette.materials``     ``navette._materials``
==================  ========================

Pure-Python subpackages (no build needed):

- ``navette.structure`` — layer stacks, architect, solver arrays
- ``navette.config`` — YAML/JSON material libraries and stack configs
- ``navette.data`` — bundled CIE reference spectra

Build the Rust extension with::

    maturin develop --release

Use ``--release``. Plain ``maturin develop`` produces a *debug* build: it
imports and computes correctly but runs several times slower, which silently
invalidates every timing taken against it. ``build_profile()`` reports which
one is installed, and the benches refuse to run on ``"debug"``.

Wrappers raise a helpful ``ImportError`` (with the exact maturin command)
when their native module is missing, so ``import navette`` always works.
"""

from __future__ import annotations

from .__about__ import (
  __title__,
  __version__,
  __description__,
  __author__,
  __license__,
  __copyright__,
  metadata_summary,
)


def build_profile() -> str:
  """Returns the Cargo profile of the installed extension.

  Either ``"release"`` or ``"debug"``. Imported lazily so that
  ``import navette`` still succeeds without a built extension.

  A ``"debug"`` build is functionally correct but several times slower, so
  any benchmark run against one is meaningless -- see
  ``validation/benches/_bench_common.py::require_release``, which exits
  rather than produce such numbers.

  Raises:
    ImportError: If the native extension is not built.
  """
  from navette._navette import build_profile as _native_build_profile

  return str(_native_build_profile())


__all__ = [
  "build_profile",
  "__title__",
  "__version__",
  "__description__",
  "__author__",
  "__license__",
  "__copyright__",
  "metadata_summary",
]
