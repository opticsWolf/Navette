#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""``core_engine`` against the numba kernel it replaced, across grid sizes.

R2.4a put the two engines on the same inputs for the first time and found the
rewrite running at 0.5-0.65x the kernel it replaced. That number came out of
the parity scripts' speed section, which times one 50x10 grid -- 500 points,
small enough that fixed per-call cost dominates and a 20% swing between runs is
normal. R5.3 needed the shape of the curve, not one point on it, and the shape
is where the answer was: the cost was in *building* the solver, so it scaled
with the grid and was invisible in a ratio taken at a single size.

This bench is that curve, kept so the ratio is tracked rather than
rediscovered. It sweeps the grid from 500 to 60 000 points at two request
breadths:

``photometry``
    Rs/Rp/Ts/Tp -- exactly what the numba ``core_engine_photometry_only``
    kernel computes, so the comparison is like for like.

``rigorous``
    the twelve-channel ellipsometric request, against
    ``core_engine_rigorous_ellipsometry``. More emitted channels per solved
    point, which is what separates the two columns.

Usage:  python validation/benches/smatrix/bench_core_engine_scaling.py [--out FILE]
Writes JSON {unit: {points: {numba_ms, rust_ms, ratio}}} to stdout (and FILE).
Without numba installed it still runs and reports the Rust timings alone.
"""

# --- bench preamble (R2.1/R2.2) ------------------------------------------
# UTF-8 console, repo `src/` on sys.path, and a hard gate against timing a
# debug build. Must precede any `navette` import.
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from _bench_common import bench_provenance, require_release, setup_bench  # noqa: E402

setup_bench()
require_release()
# -------------------------------------------------------------------------

import json
import os
import sys
import time

import numpy as np

import navette._smatrix as rust_mod
from navette.smatrix.smatrix import CoherenceMode, Request

# The numba reference lives beside the parity scripts that own it.
_REFS = _Path(__file__).resolve().parents[2] / "parity" / "smatrix" / "refs"
if _REFS.is_dir():
    sys.path.insert(0, str(_REFS))
try:
    from loom_matrix import core_engine_photometry_only as numba_photometry
    from loom_matrix import core_engine_rigorous_ellipsometry as numba_rigorous

    HAVE_NUMBA = True
except Exception as exc:  # pragma: no cover - optional reference
    numba_photometry = numba_rigorous = None
    HAVE_NUMBA = False
    _NUMBA_ERR = exc

N_LAYERS = 6
THICK = np.array([50.0, 100.0, 200.0, 150.0, 80.0, 500.0])
ZERO_I = np.zeros(N_LAYERS, dtype=np.int32)
ZERO_F = np.zeros(N_LAYERS)

# One angle, so "points" is the wavelength count and the sweep is a clean
# 1-D scaling curve. The per-point work does not depend on which axis the
# points came from.
GRIDS = [500, 2_000, 5_000, 20_000, 60_000]

PHOTOMETRY = Request.RS | Request.RP | Request.TS | Request.TP
RIGOROUS = (Request.RS | Request.RP | Request.TS | Request.TP
            | Request.R_AVG | Request.T_AVG
            | Request.PSI_R | Request.DELTA_R | Request.DOP_R
            | Request.PSI_T | Request.DELTA_T | Request.DOP_T)

WARMUP = 3
REPS = 25
COOLDOWN_S = 1.0

def cooldown():
    """Let the other engine's thread pool go quiet before timing this one.

    numba's default threading layer leaves its workers spin-waiting after a
    call returns, and a spinning worker holds a core while the *other* engine
    is being timed. Measured at 20 000 points: the same Rust call is 1.52 ms
    run on its own, 2.00 ms immediately after a block of numba calls, and
    1.58 ms again after a one-second pause. Forcing numba onto the
    non-spinning `workqueue` layer also removes it, but at the cost of making
    numba itself 70% slower here -- which is not a fair comparison either.
    A pause leaves both engines in their native configuration (R5.3).
    """
    time.sleep(COOLDOWN_S)



def make_case(n_points):
    """One stack on an `n_points` grid, in both engines' input layouts."""
    wavls = np.linspace(400.0, 800.0, n_points)
    sin_theta = np.sin(np.radians(np.array([12.0])))
    nr = [1.0, 1.6, 2.2, 1.9, 2.6, 1.0]
    ni = [0.0, 0.005, 0.03, 0.01, 0.05, 0.0]
    ramp = ((wavls - 400.0) / 400.0)[:, None]
    n_complex = (np.array(nr)[None, :] + 0.05 * ramp
                 + 1j * np.array(ni)[None, :]).astype(np.complex128)
    n_flat = np.empty(n_points * N_LAYERS * 2)
    n_flat[0::2] = n_complex.real.ravel()
    n_flat[1::2] = n_complex.imag.ravel()
    return dict(wavls=wavls, sin_theta=sin_theta, n_complex=n_complex,
                n_flat=n_flat)


def call_rust(case, requested):
    return rust_mod.core_engine(
        case["wavls"], case["sin_theta"], N_LAYERS, case["n_flat"], THICK,
        ZERO_I, ZERO_I, ZERO_F, int(CoherenceMode.FRONT_BLOCK), int(requested))


def call_numba(case, func, tail):
    """`tail` is the kernel's trailing flags: the photometry kernel takes
    ``calc_s, calc_p``, the rigorous one a single ``debug_flag``."""
    return func(case["wavls"], case["sin_theta"], N_LAYERS, case["n_complex"],
                THICK, ZERO_I, ZERO_I, ZERO_F, *tail)


def timeit(fn):
    for _ in range(WARMUP):
        fn()
    ts = []
    for _ in range(REPS):
        t0 = time.perf_counter()
        fn()
        ts.append((time.perf_counter() - t0) * 1e3)
    return float(np.median(ts))


UNITS = [
    ("photometry", PHOTOMETRY, "numba_photometry", (np.int32(1), np.int32(1))),
    ("rigorous", RIGOROUS, "numba_rigorous", (np.int32(0),)),
]


def main():
    res = {}
    print(f"core_engine scaling ({N_LAYERS} layers, 1 angle; median of {REPS})")
    if not HAVE_NUMBA:
        print(f"  numba reference unavailable ({_NUMBA_ERR}); Rust timings only")
    for unit, mask, ref_name, tail in UNITS:
        ref = globals()[ref_name]
        rows = {}
        print(f"\n  {unit}")
        print(f"    {'points':>8s} {'numba (ms)':>12s} {'rust (ms)':>11s} {'ratio':>8s}")
        for n in GRIDS:
            case = make_case(n)
            cooldown()
            rust_ms = timeit(lambda: call_rust(case, mask))
            row = {"rust_ms": rust_ms}
            line = f"    {n:>8,d} {'-':>12s} {rust_ms:>11.3f} {'-':>8s}"
            if ref is not None:
                cooldown()
                numba_ms = timeit(lambda: call_numba(case, ref, tail))
                row["numba_ms"] = numba_ms
                row["ratio"] = numba_ms / rust_ms if rust_ms > 0 else float("inf")
                line = (f"    {n:>8,d} {numba_ms:>12.3f} {rust_ms:>11.3f} "
                        f"{row['ratio']:>7.2f}x")
            print(line)
            rows[str(n)] = row
        res[unit] = rows
    # Stamp the build the numbers came from -- results without it are of
    # unknown profile and must not be compared against these (R2.1).
    res["_provenance"] = bench_provenance()
    res["_provenance"]["numba"] = HAVE_NUMBA
    text = json.dumps(res, indent=1)
    if len(sys.argv) == 3 and sys.argv[1] == "--out":
        with open(sys.argv[2], "w") as f:
            f.write(text)
    elif os.environ.get("NAVETTE_BENCH_JSON"):
        print(text)


if __name__ == "__main__":
    main()
