"""Garbage-in behaviour of the public solver entry points (review harness).

Originally this script *documented* what the library did with malformed
input, and what it documented was bad: negative and NaN thicknesses both
returned the numbers for a stack with that layer deleted, 120 deg aliased
onto 60 deg, a duplicated wavelength produced NaN dispersion channels. None
of it raised.

R3.1 turned the ScatterMatrix cases into construction-time errors and R3.2
the needle ones, so the script now asserts an expected verdict per case and
exits non-zero when one drifts.

Verdicts: `raises` | `silent-clean` | `NaN-in-output`.
"""

import warnings
import numpy as np
from navette.smatrix.smatrix import ScatterMatrix, Request
from navette.smatrix.needle import NeedleRequest, needle_gradient

FAILURES = []


def run(name, expect, fn, note=""):
    """Run `fn`, classify the outcome, and compare with `expect`."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            out = fn()
            verdict, detail = "silent-clean", ""
            first = None
            for v in (out.values() if isinstance(out, dict) else [out]):
                a = np.asarray(v)
                if first is None:
                    first = a
                flat = a.view(float) if np.issubdtype(a.dtype, np.complexfloating) else a.astype(float)
                if flat.size and np.isnan(flat).any():
                    verdict, detail = "NaN-in-output", ""
                    break
            if verdict == "silent-clean":
                detail = f"sample={np.asarray(first).ravel()[:2]}"
            detail += f" (warnings={len(w)})"
        except Exception as e:
            verdict = "raises"
            detail = f"{type(e).__name__}: {str(e)[:90]}"
    mark = "ok " if verdict == expect else "DRIFT"
    if verdict != expect:
        FAILURES.append(f"{name}: expected {expect}, got {verdict}")
    tail = f"  [{note}]" if note else ""
    print(f"  {mark} {name}: {verdict} {detail}{tail}")

N = np.array([1.0+0j, 2.35+0j, 1.46+0j, 2.10+0j, 1.52+0j])
D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
WLS = np.linspace(450.0, 750.0, 31)

print("=== ScatterMatrix garbage-in ===")
run("NaN index in one layer", "raises",
    lambda: ScatterMatrix(np.array([1.0+0j, np.nan, 1.46+0j, 2.10+0j, 1.52+0j]),
                          D, wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("inf (lossless-style) index", "raises",
    lambda: ScatterMatrix(np.array([1.0+0j, 1e308+0j, 1.46+0j, 2.10+0j, 1.52+0j]),
                          D, wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("NaN wavelength", "raises",
    lambda: ScatterMatrix(N, D, wavelengths=np.where(WLS == WLS[5], np.nan, WLS),
                          angles=[30.0]).compute(Request.RS))
run("angle 90 deg (grazing, sin=1)", "silent-clean",
    lambda: ScatterMatrix(N, D, wavelengths=WLS, angles=[90.0]).compute(Request.RS), note="grazing is physical: R = 1")
run("angle 95 deg (sin>1, TIR-invalid)", "raises",
    lambda: ScatterMatrix(N, D, wavelengths=WLS, angles=[95.0]).compute(Request.RS))
run("angle 120 deg (retrograde)", "raises",
    lambda: ScatterMatrix(N, D, wavelengths=WLS, angles=[120.0]).compute(Request.RS))
run("single wavelength grid", "silent-clean",
    lambda: ScatterMatrix(N, D, wavelengths=np.array([550.0]), angles=[30.0]).compute(
        Request.RS | Request.DISP_R_S), note="no spacing to differentiate; dispersion is 0, not NaN")
run("two identical wavelengths", "raises",
    lambda: ScatterMatrix(N, D, wavelengths=np.array([550.0, 550.0]), angles=[30.0]).compute(
        Request.RS | Request.DISP_R_S))
run("descending wavelengths", "raises",
    lambda: ScatterMatrix(N, D, wavelengths=WLS[::-1].copy(), angles=[30.0]).compute(Request.RS))
run("negative thickness (already known, S19)", "raises",
    lambda: ScatterMatrix(N, np.array([0.0, 120.0, -50.0, 80.0, 0.0]),
                          wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("NaN thickness", "raises",
    lambda: ScatterMatrix(N, np.array([0.0, 120.0, np.nan, 80.0, 0.0]),
                          wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("huge thickness 1e9 nm", "silent-clean",
    lambda: ScatterMatrix(N, np.array([0.0, 120.0, 1e9, 80.0, 0.0]),
                          wavelengths=WLS, angles=[30.0]).compute(Request.RS), note="no upper cap by design")

print("=== needle_gradient garbage-in ===")
st = ScatterMatrix(N, D, wavelengths=WLS, angles=[30.0])
NN = np.full(WLS.size, 1.8+0.05j, dtype=np.complex128)
run("NaN needle index", "raises",
    lambda: needle_gradient(st, np.where(WLS == WLS[3], np.nan+0j, NN), [220.0],
                            NeedleRequest.P, pol="s"), note="R3.2 -- fixed")
run("needle index below total-internal floor", "silent-clean",
    lambda: needle_gradient(st, np.full(WLS.size, 0.3+0j), [220.0],
                            NeedleRequest.P, pol="s"), note="n < 1 is a real metallic index, not garbage")
run("z outside the stack (z=1000)", "raises",
    lambda: needle_gradient(st, NN, [1000.0], NeedleRequest.P, pol="s"), note="R3.2 -- fixed")
run("negative z", "raises",
    lambda: needle_gradient(st, NN, [-10.0], NeedleRequest.P, pol="s"), note="R3.2 -- fixed")

print()
if FAILURES:
    print(f"DRIFT: {len(FAILURES)} case(s) no longer behave as recorded:")
    for f in FAILURES:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL AS EXPECTED")
