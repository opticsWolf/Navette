# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The quick-start snippet shown on the site, kept runnable.

The page quotes this file between the two markers; CI runs it, so the
snippet on the site cannot drift away from the API it advertises.
"""
# --- snippet ---
import numpy as np
from navette.smatrix.smatrix import Request, ScatterMatrix

wl = np.linspace(380.0, 780.0, 401)            # nm
H, L, lam0 = 2.35, 1.46, 560.0                 # TiO2 / SiO2, design wavelength

n = [1.0] + [H, L] * 8 + [1.52]                # ambient, 16 films, glass
d = [0.0] + [lam0 / (4 * H), lam0 / (4 * L)] * 8 + [0.0]
nk = np.array(n, dtype=complex)[:, None] * np.ones_like(wl)

mirror = ScatterMatrix(nk, d, wavelengths=wl, angles=[0.0, 45.0])
out = mirror.compute(Request.RS | Request.RP)  # arrays of shape (angles, wavelengths)
print(f"peak Rs at normal incidence: {out['Rs'][0].max():.4f}")
# --- end snippet ---

assert out["Rs"].shape == (2, wl.size)
assert out["Rs"][0].max() > 0.99
