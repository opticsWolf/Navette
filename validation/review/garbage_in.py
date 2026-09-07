import warnings
import numpy as np
from navette.smatrix.smatrix import ScatterMatrix, Request
from navette.smatrix.needle import NeedleRequest, needle_gradient

OK = True
def report(name, verdict, detail=""):
    # verdict: "silent-garbage" | "error" | "clean" | "nan"
    global OK
    print(f"  {name}: {verdict} {detail}")

def run(name, fn):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            out = fn()
            first = None
            for v in (out.values() if isinstance(out, dict) else [out]):
                a = np.asarray(v, float)
                if first is None:
                    first = a
                if np.asarray(a).size and (np.isnan(np.asarray(a, float)).any()
                                           if not np.issubdtype(np.asarray(a).dtype, np.complexfloating)
                                           else np.isnan(np.asarray(a).view(float)).any()):
                    report(name, "NaN-in-output", f"(warnings={len(w)})")
                    return
            report(name, "silent-clean", f"sample={np.asarray(first).ravel()[:2]} (warnings={len(w)})")
        except Exception as e:
            report(name, "raises", f"{type(e).__name__}: {str(e)[:70]}")

N = np.array([1.0+0j, 2.35+0j, 1.46+0j, 2.10+0j, 1.52+0j])
D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
WLS = np.linspace(450.0, 750.0, 31)

print("=== ScatterMatrix garbage-in ===")
run("NaN index in one layer",
    lambda: ScatterMatrix(np.array([1.0+0j, np.nan, 1.46+0j, 2.10+0j, 1.52+0j]),
                          D, wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("inf (lossless-style) index",
    lambda: ScatterMatrix(np.array([1.0+0j, 1e308+0j, 1.46+0j, 2.10+0j, 1.52+0j]),
                          D, wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("NaN wavelength",
    lambda: ScatterMatrix(N, D, wavelengths=np.where(WLS == WLS[5], np.nan, WLS),
                          angles=[30.0]).compute(Request.RS))
run("angle 90 deg (grazing, sin=1)",
    lambda: ScatterMatrix(N, D, wavelengths=WLS, angles=[90.0]).compute(Request.RS))
run("angle 95 deg (sin>1, TIR-invalid)",
    lambda: ScatterMatrix(N, D, wavelengths=WLS, angles=[95.0]).compute(Request.RS))
run("angle 120 deg (retrograde)",
    lambda: ScatterMatrix(N, D, wavelengths=WLS, angles=[120.0]).compute(Request.RS))
run("single wavelength grid",
    lambda: ScatterMatrix(N, D, wavelengths=np.array([550.0]), angles=[30.0]).compute(
        Request.RS | Request.DISP_R_S))
run("two identical wavelengths",
    lambda: ScatterMatrix(N, D, wavelengths=np.array([550.0, 550.0]), angles=[30.0]).compute(
        Request.RS | Request.DISP_R_S))
run("descending wavelengths",
    lambda: ScatterMatrix(N, D, wavelengths=WLS[::-1].copy(), angles=[30.0]).compute(Request.RS))
run("negative thickness (already known, S19)",
    lambda: ScatterMatrix(N, np.array([0.0, 120.0, -50.0, 80.0, 0.0]),
                          wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("NaN thickness",
    lambda: ScatterMatrix(N, np.array([0.0, 120.0, np.nan, 80.0, 0.0]),
                          wavelengths=WLS, angles=[30.0]).compute(Request.RS))
run("huge thickness 1e9 nm",
    lambda: ScatterMatrix(N, np.array([0.0, 120.0, 1e9, 80.0, 0.0]),
                          wavelengths=WLS, angles=[30.0]).compute(Request.RS))

print("=== needle_gradient garbage-in ===")
st = ScatterMatrix(N, D, wavelengths=WLS, angles=[30.0])
NN = np.full(WLS.size, 1.8+0.05j, dtype=np.complex128)
run("NaN needle index",
    lambda: needle_gradient(st, np.where(WLS == WLS[3], np.nan+0j, NN), [220.0],
                            NeedleRequest.P, pol="s"))
run("needle index below total-internal floor",
    lambda: needle_gradient(st, np.full(WLS.size, 0.3+0j), [220.0],
                            NeedleRequest.P, pol="s"))
run("z outside the stack (z=1000)",
    lambda: needle_gradient(st, NN, [1000.0], NeedleRequest.P, pol="s"))
run("negative z",
    lambda: needle_gradient(st, NN, [-10.0], NeedleRequest.P, pol="s"))
