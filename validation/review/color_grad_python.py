# -*- coding: utf-8 -*-
"""R4.2: the color gradient on the *Python* needle path.

The documented Python flow is `build_needle_targets` -> `needle_gradient`.
For pointwise demands it worked. For a color demand -- Lab/DE2000, White,
Yellow, dominant wavelength -- it returned zero, with no error:

  * a color demand integrates the whole spectrum into ONE residual, so it
    has no per-point (target, weight) pair; the native fold instead emits
    `grad_r`/`grad_t`, the analytic dF/dcurve per solver point;
  * `build_needle_targets`' Python dict dropped those two arrays; and
  * `needle_gradient` had no argument that could have accepted them.

So the folded `r` pair was all zeros, `needle_gradient` dutifully returned
a gradient of zero, and a caller assembling their own needle cycle in
Python optimized a color demand against nothing. `run_design` was never
affected -- it folds and deposits entirely inside Rust.

This harness drives the fixed Python path end to end and checks it four
ways: the fold's own numbers against finite differences of the merit, the
deposit against a needle call it can be re-derived from (1e-12, the two are
the same kernel and must agree to arithmetic noise), the whole chain
against a layer-thickness FD of the merit, and the contract that a color
bucket handed to a channel that is not being computed is an error rather
than a silent drop -- which is the original bug wearing a new hat.
"""

import json

import numpy as np

from navette._smatrix import SimCurves, build_needle_targets, compile_merit_spec
from navette.smatrix.needle import NeedleRequest, needle_gradient
from navette.smatrix.smatrix import Request, ScatterMatrix

FAILURES = []


def report(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'} {name}" + (f"  {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


# --------------------------------------------------------------------------
# The same 3-layer coating `color_merit_check.py` uses, so the two harnesses
# can be read against each other.
# --------------------------------------------------------------------------
N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
WLS = np.linspace(450.0, 750.0, 31)
ANGLES = np.array([30.0])
WEIGHT = 1.7
REF = [55.0, 12.0, -18.0]
NW = WLS.size

# Needle host depths, absolute from the top of layer 1, and the host index.
LAY_Z = {1: 60.0, 2: 220.0, 3: 360.0}
LAY_N = {1: 2.35 + 0j, 2: 1.46 + 0j, 3: 2.10 + 0j}

st0 = ScatterMatrix(N, D, wavelengths=WLS, angles=ANGLES)
out0 = st0.compute(Request.RS | Request.TS)
RS = np.atleast_2d(np.asarray(out0["Rs"], float))
TS = np.atleast_2d(np.asarray(out0["Ts"], float))


def make_spec(curve="Rs", quantity="Lab", distance="DeltaE2000", reference=REF):
    doc = {
        "spectral": [], "angular": [],
        "color": [{
            "curve": curve, "angle": 30.0,
            "illuminant": "D65", "observer": "1931_2deg",
            "quantity": quantity, "reference": reference,
            "distance": distance, "weight": WEIGHT,
        }],
    }
    return compile_merit_spec(json.dumps(doc))


def sim_of(rs=None, ts=None):
    sim = SimCurves(ANGLES, WLS, float(D[1:-1].sum()), 1.0, 1.52)
    sim.set_curve("Rs", np.ascontiguousarray(
        (RS if rs is None else rs).ravel()))
    sim.set_curve("Ts", np.ascontiguousarray(
        (TS if ts is None else ts).ravel()))
    return sim


SPEC = make_spec()
FOLD = build_needle_targets(SPEC, ANGLES, WLS, sim_of())


# --------------------------------------------------------------------------
# 1. The fold emits the buckets at all, and they are dF/dR
# --------------------------------------------------------------------------
print("=== 1. build_needle_targets emits the color buckets ===")

report("grads_r / grads_t present in the fold dict",
       "grads_r" in FOLD and "grads_t" in FOLD,
       f"keys={sorted(FOLD)}")

gr = np.asarray(FOLD["grads_r"], float)
gt = np.asarray(FOLD["grads_t"], float)
report("grads_r has one entry per solver point", gr.shape == (NW,), f"{gr.shape}")
report("grads_r is not identically zero", np.any(gr != 0.0),
       f"max|g|={np.max(np.abs(gr)):.6g}")
report("grads_t is zero for an Rs-only color demand", np.all(gt == 0.0))

# The pointwise buckets must stay empty: a color demand has no per-point
# target, and a non-zero weight here would double-count it.
report("the r (target, weight) pair stays all-zero",
       np.all(np.asarray(FOLD["r"]["weights"], float) == 0.0)
       and np.all(np.asarray(FOLD["r"]["targets"], float) == 0.0))


def merit_of_rs(rs):
    return SPEC.merit(sim_of(rs=rs), 1e6)


# dF/dR_i by central difference on the sim row: the fold claims to have
# computed exactly this analytically.
F0 = merit_of_rs(RS)
fd_dfdr = np.zeros(NW)
h = 1e-6
for i in range(NW):
    rp = RS.copy(); rp[0, i] += h
    rm = RS.copy(); rm[0, i] -= h
    fd_dfdr[i] = (merit_of_rs(rp) - merit_of_rs(rm)) / (2 * h)

rel = np.max(np.abs(gr - fd_dfdr)) / max(np.max(np.abs(fd_dfdr)), 1e-300)
report("grads_r == dF/dR by finite difference (<= 1e-6)", rel < 1e-6,
       f"max rel dev {rel:.3e}")


# --------------------------------------------------------------------------
# 2. The regression itself: without the buckets the gradient is zero
# --------------------------------------------------------------------------
print("=== 2. the silent zero, before and after ===")


def grad_at(layer, *, with_color, request=NeedleRequest.P, **kw):
    return needle_gradient(
        st0, np.full(NW, LAY_N[layer]), [LAY_Z[layer]], request,
        targets_r=np.asarray(FOLD["r"]["targets"], float),
        weights_r=np.asarray(FOLD["r"]["weights"], float),
        grads_r=gr if with_color else None,
        **kw)


old = grad_at(2, with_color=False)["P_s"]
new = grad_at(2, with_color=True)["P_s"]
report("omitting the buckets reproduces the zero gradient",
       np.all(old == 0.0), f"max|P|={np.max(np.abs(old)):.3e}")
report("the fixed path returns a non-zero gradient",
       np.max(np.abs(new)) > 0.0, f"max|P|={np.max(np.abs(new)):.3e}")


# --------------------------------------------------------------------------
# 3. The deposit against a needle call it can be re-derived from (1e-12)
# --------------------------------------------------------------------------
print("=== 3. deposit vs the same kernel driven through the R channel ===")
# Both kernels multiply the SAME Re{conj(r) dr/ddz} by a scalar: the
# pointwise one by 2*w*(R - t), the color one by g. Driving the pointwise
# channel with targets_r = 0, weights_r = 1 therefore yields
# P_ref = 2*R*Re{...}, so the deposit must be g * P_ref / (2*R) point for
# point -- to arithmetic noise, not to FD accuracy.
for lay in (1, 2, 3):
    ref = needle_gradient(st0, np.full(NW, LAY_N[lay]), [LAY_Z[lay]],
                          NeedleRequest.P, targets_r=np.zeros(NW),
                          weights_r=np.ones(NW))["P_s"].ravel()
    got = grad_at(lay, with_color=True)["P_s"].ravel()
    want = 0.5 * gr * ref / RS.ravel()
    dev = np.max(np.abs(got - want)) / max(np.max(np.abs(want)), 1e-300)
    report(f"layer {lay}: deposit == g * P_ref / (2*R) (<= 1e-12)", dev < 1e-12,
           f"max rel dev {dev:.3e}")


# --------------------------------------------------------------------------
# 4. Additivity: a mixed spec is the sum of its parts
# --------------------------------------------------------------------------
print("=== 4. pointwise and color superpose on the same channel ===")
# A color bucket riding the R channel must not disturb the pointwise term
# already there; the native pass adds both into one accumulator and so
# must this one.
tg = np.full(NW, 0.05)
wg = np.full(NW, 3.0)
kw = dict(targets_r=tg, weights_r=wg)
only_pt = needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                          NeedleRequest.P, **kw)["P_s"]
only_col = needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                           NeedleRequest.P, targets_r=np.zeros(NW),
                           weights_r=np.zeros(NW), grads_r=gr)["P_s"]
both = needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                       NeedleRequest.P, grads_r=gr, **kw)["P_s"]
dev = np.max(np.abs(both - (only_pt + only_col))) / max(
    np.max(np.abs(both)), 1e-300)
report("P(pointwise + color) == P(pointwise) + P(color) (<= 1e-12)",
       dev < 1e-12, f"max rel dev {dev:.3e}")


# --------------------------------------------------------------------------
# 5. End to end: the chain rule, and a thickness FD of the merit
# --------------------------------------------------------------------------
print("=== 5. dF/dthickness through the Python path ===")
# P is the HALF gradient (docstring convention), so dF/ddelta = 2*sum_k P_k;
# a needle of the host's own index inside that host is layer growth.


def merit_of_thickness(th):
    st = ScatterMatrix(N, np.asarray(th, float), wavelengths=WLS, angles=ANGLES)
    o = st.compute(Request.RS)
    return SPEC.merit(sim_of(rs=np.atleast_2d(np.asarray(o["Rs"], float))), 1e6)


# The hand-assembled chain rule of color_merit_check.py Part C, rebuilt
# here on the fold's own g rather than on FD residual derivatives.
h_th = 0.01
for lay in (1, 2, 3):
    th = D.copy()
    step = np.eye(D.size)[lay] * h_th
    fd = (merit_of_thickness(th + step) - merit_of_thickness(th - step)) / (2 * h_th)

    p = grad_at(lay, with_color=True)["P_s"].ravel()
    analytic = 2.0 * float(np.sum(p))

    ref = needle_gradient(st0, np.full(NW, LAY_N[lay]), [LAY_Z[lay]],
                          NeedleRequest.P, targets_r=np.zeros(NW),
                          weights_r=np.ones(NW))["P_s"].ravel()
    chain = float(np.sum(gr * (ref / RS.ravel())))

    rel_fd = abs(analytic - fd) / max(abs(fd), abs(analytic), 1e-300)
    rel_ch = abs(analytic - chain) / max(abs(chain), 1e-300)
    report(f"layer {lay}: vs merit FD (<= 2e-3)", rel_fd < 2e-3,
           f"python={analytic:.6g} fd={fd:.6g} rel={rel_fd:.2e}")
    report(f"layer {lay}: vs the hand-assembled chain rule (<= 1e-12)",
           rel_ch < 1e-12, f"rel={rel_ch:.2e}")


# --------------------------------------------------------------------------
# 6. The T channel carries its own bucket
# --------------------------------------------------------------------------
print("=== 6. a transmission color demand rides P_T, not P ===")
spec_t = make_spec(curve="Ts")
fold_t = build_needle_targets(spec_t, ANGLES, WLS, sim_of())
gr_t = np.asarray(fold_t["grads_r"], float)
gt_t = np.asarray(fold_t["grads_t"], float)
report("grads_t is populated and grads_r is empty",
       np.any(gt_t != 0.0) and np.all(gr_t == 0.0),
       f"max|gt|={np.max(np.abs(gt_t)):.6g}")

o = needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                    NeedleRequest.P | NeedleRequest.P_T,
                    targets_r=np.zeros(NW), weights_r=np.zeros(NW),
                    targets_t=np.zeros(NW), weights_t=np.zeros(NW),
                    grads_r=gr_t, grads_t=gt_t)
report("P_T_s is non-zero", np.max(np.abs(o["P_T_s"])) > 0.0,
       f"max|P_T|={np.max(np.abs(o['P_T_s'])):.3e}")
report("P_s stays zero (nothing leaked across channels)",
       np.all(o["P_s"] == 0.0))

ref_t = needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                        NeedleRequest.P_T, targets_t=np.zeros(NW),
                        weights_t=np.ones(NW))["P_T_s"].ravel()
want = 0.5 * gt_t * ref_t / TS.ravel()
dev = np.max(np.abs(o["P_T_s"].ravel() - want)) / max(np.max(np.abs(want)), 1e-300)
report("T deposit == g * P_T_ref / (2*T) (<= 1e-12)", dev < 1e-12,
       f"max rel dev {dev:.3e}")


# --------------------------------------------------------------------------
# 7. The contract: a bucket with nowhere to land is an error
# --------------------------------------------------------------------------
print("=== 7. a dropped bucket is reported, not swallowed ===")


def raises(fn, fragment):
    try:
        fn()
    except ValueError as exc:
        return fragment in str(exc), str(exc)
    except Exception as exc:  # noqa: BLE001 - any other type is also a failure
        return False, f"{type(exc).__name__}: {exc}"
    return False, "no exception"


ok, msg = raises(
    lambda: needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                            NeedleRequest.P_T, grads_r=gr),
    "NREQ_P was not requested")
report("grads_r without NREQ_P raises", ok, msg)

ok, msg = raises(
    lambda: needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                            NeedleRequest.P, grads_t=gt_t),
    "NREQ_P_T was not requested")
report("grads_t without NREQ_P_T raises", ok, msg)

bad = gr.copy(); bad[7] = np.nan
ok, msg = raises(
    lambda: needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                            NeedleRequest.P, grads_r=bad),
    "grads_r[7]")
report("a NaN bucket names its index", ok, msg)

# An all-zero array is not a demand: the fold hands both arrays through
# unfiltered, so a color-free spec must not be forced to filter them.
try:
    needle_gradient(st0, np.full(NW, LAY_N[2]), [LAY_Z[2]],
                    NeedleRequest.P, grads_t=np.zeros(NW))
    report("an all-zero unused bucket is accepted", True)
except Exception as exc:  # noqa: BLE001
    report("an all-zero unused bucket is accepted", False, str(exc))


print()
if FAILURES:
    print(f"DRIFT in {len(FAILURES)} check(s): {FAILURES}")
    raise SystemExit(1)
print("ALL OK")
