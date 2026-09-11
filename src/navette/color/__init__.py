# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""
``navette.color`` — colorimetry (sRGB/XYZ/Lab/LCH/LUV/Oklab, ΔE, photometry).

Thin wrapper over the private native extension ``navette._color`` built
from ``rust/navette/src/color``::

    maturin develop --release  # from the repo root
"""

from __future__ import annotations

try:
  from navette._color import *  # noqa: F401,F403
  from navette._color import __doc__ as _native_doc  # noqa: F401
except ImportError as exc:  # pragma: no cover - native not built yet
  raise ImportError(
    "Could not import the compiled `navette._color` extension. "
    "Build the Rust crate so it is importable, then retry:\n"
    "    maturin develop --release  # from the repo root"
  ) from exc
