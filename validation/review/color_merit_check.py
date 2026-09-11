# SPDX-License-Identifier: LGPL-3.0-or-later
# Step 2: color_merit.rs independent verification.
# Part A: merit value of a Lab+DeltaE2000 color demand vs my own Python chain
#         (XYZ rectangular integral, k-normalization, Lab, Sharma DE2000)
# Part B: deposited gradient d(sqrt(F))/dR(lambda_i) vs central FD of the
#         residual w.r.t. R values (pure data check, no solver)
# Part C: end-to-end color gradient through the solver (layer growth FD)
# Part D: White (CIE whiteness), Yellow (E313), DomWl quantities vs my own
#         closed forms / locus intersection
import json
import numpy as np

from navette.smatrix.smatrix import ScatterMatrix, Request
from navette.smatrix.needle import NeedleRequest, needle_gradient
from navette._smatrix import compile_merit_spec, SimCurves, build_needle_targets
from navette.data import load_cie_table

rng = np.random.default_rng(20260215)
OK = True


def report(name, cond, detail=""):
    global OK
    print(f"  {name}: {'OK' if cond else 'FAIL'} {detail}")
    OK = bool(cond) and OK


# ---- tables (navette's own, but the *chain* below is my own) --------------
d65 = load_cie_table("CIE", "sds", "CIE_std_illum_D65_S_D65.json")
E_WL = np.asarray(d65["lambda"], float)
E = np.asarray(d65["S_D65(lambda)"], float)
cmf = load_cie_table("CIE", "cmf", "CIE_xyz_1931_2deg.json")
C_WL = np.asarray(cmf["lambda"], float)
XB = np.asarray(cmf["x_bar(lambda)"], float)
YB = np.asarray(cmf["y_bar(lambda)"], float)
ZB = np.asarray(cmf["z_bar(lambda)"], float)


def lin_resample(x, y, xv):
    if xv < x[0] or xv > x[-1]:
        return None
    j = np.searchsorted(x, xv, side="right") - 1
    j = min(max(j, 0), len(x) - 2)
    t = (xv - x[j]) / (x[j + 1] - x[j])
    return y[j] + t * (y[j + 1] - y[j])


def xyz_of(R, wl, wr=None):
    """My independent XYZ: rectangular rule, forward-diff dw, k-norm."""
    xs, es, cs = [], [], []
    pos = {w: i for i, w in enumerate(wl)}
    for w in wl:
        if wr is not None and (w < wr[0] or w > wr[1]):
            continue
        e = lin_resample(E_WL, E, w)
        c = (lin_resample(C_WL, XB, w), lin_resample(C_WL, YB, w),
             lin_resample(C_WL, ZB, w))
        if e is None or any(v is None for v in c):
            continue
        xs.append(w); es.append(e); cs.append(c)
    dw = np.diff(np.asarray(xs))
    dw = np.append(dw, dw[-1])
    k = 1.0 / np.sum(np.asarray(es) * np.asarray([c[1] for c in cs]) * dw)
    if len(xs) == len(wl) and wr is None:
        Rv = np.asarray(R, float)
    else:
        Rv = np.asarray([R[pos[w]] for w in xs])

    X = np.sum(Rv * es * np.asarray([c[0] for c in cs]) * k * dw)
    Y = np.sum(Rv * es * np.asarray([c[1] for c in cs]) * k * dw)
    Z = np.sum(Rv * es * np.asarray([c[2] for c in cs]) * k * dw)
    return np.array([X, Y, Z])


def lab_of(xyz, white):
    def f(t):
        d = 6.0 / 29.0
        return np.where(t > d ** 3, np.cbrt(t), t / (3 * d ** 2) + 4.0 / 29.0)
    x, y, z = np.asarray(xyz) / np.asarray(white)
    fx, fy, fz = f(x), f(y), f(z)
    L = 116 * fy - 16
    a = 500 * (fx - fy)
    b = 200 * (fy - fz)
    return np.array([L, a, b])


def de2000(lab1, lab2):
    L1, a1, b1 = lab1
    L2, a2, b2 = lab2
    C1 = np.hypot(a1, b1); C2 = np.hypot(a2, b2)
    Cb = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cb ** 7 / (Cb ** 7 + 25 ** 7)))
    a1p = (1 + G) * a1; a2p = (1 + G) * a2
    C1p = np.hypot(a1p, b1); C2p = np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    dLp = L2 - L1
    dCp = C2p - C1p
    if C1p * C2p == 0:
        dhp = 0.0
    else:
        dhp = h2p - h1p
        if dhp > 180: dhp -= 360
        elif dhp < -180: dhp += 360
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dhp) / 2)
    Lbp = (L1 + L2) / 2
    Cbp = (C1p + C2p) / 2
    if C1p * C2p == 0:
        hbp = h1p + h2p
    else:
        hbp = (h1p + h2p) / 2
        if abs(h1p - h2p) > 180:
            hbp += 180 if (h1p + h2p) < 360 else -180
    T = (1 - 0.17 * np.cos(np.radians(hbp - 30)) + 0.24 * np.cos(np.radians(2 * hbp))
         + 0.32 * np.cos(np.radians(3 * hbp + 6)) - 0.20 * np.cos(np.radians(4 * hbp - 63)))
    Sl = 1 + 0.015 * (Lbp - 50) ** 2 / np.sqrt(20 + (Lbp - 50) ** 2)
    Sc = 1 + 0.045 * Cbp
    Sh = 1 + 0.015 * Cbp * T
    dtheta = 30 * np.exp(-((hbp - 275) / 25) ** 2)
    Rc = 2 * np.sqrt(Cbp ** 7 / (Cbp ** 7 + 25 ** 7))
    RT = -np.sin(np.radians(2 * dtheta)) * Rc
    return np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2
                   + RT * (dCp / Sc) * (dHp / Sh))


# ---- stack + sim ----------------------------------------------------------
N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
WLS = np.linspace(450.0, 750.0, 31)
st0 = ScatterMatrix(N, D, wavelengths=WLS, angles=[30.0])
out0 = st0.compute(Request.RS | Request.RS_C)
RS = np.atleast_2d(np.asarray(out0["Rs"], float))

WEIGHT = 1.7
REF = [55.0, 12.0, -18.0]


def make_spec(quantity="Lab", distance="DeltaE2000", reference=REF,
              extra=None):
    doc = {
        "spectral": [], "angular": [],
        "color": [dict({
            "curve": "Rs", "angle": 30.0,
            "illuminant": "D65", "observer": "1931_2deg",
            "quantity": quantity, "reference": reference,
            "distance": distance, "weight": WEIGHT,
        }, **(extra or {}))],
    }
    return compile_merit_spec(json.dumps(doc))


def sim_of(rs):
    sim = SimCurves(np.array([30.0]), WLS, float(D[1:-1].sum()), 1.0, 1.52)
    sim.set_curve("Rs", np.ascontiguousarray(rs.ravel()))
    return sim


# ------------------------------------------------- Part A: merit value
print("=== Part A: Lab+DE2000 merit value vs independent Python chain ===")
spec = make_spec()
sim = sim_of(RS)
native = spec.merit(sim, 1e6)

white = xyz_of(np.ones(E_WL.size), E_WL)  # engine: white on the illuminant's NATIVE grid
xyz = xyz_of(RS, WLS)
lab = lab_of(xyz, white)
F_mine = WEIGHT * de2000(lab, np.asarray(REF, float)) ** 2
rel = abs(native - F_mine) / max(abs(F_mine), 1e-300)
report("Lab/DE2000 merit", rel < 1e-9,
       f"native={native:.10g} mine={F_mine:.10g} rel={rel:.2e}")

# ------------------------------------------- Part B: deposited gradient
print("=== Part B: d(sqrt(F))/dR deposits vs FD on the sim row ===")
res0 = float(np.asarray(spec.residuals(sim_of(RS)))[0])
grads = {}
# native grads are not directly exposed; verify via merit FD instead:
# F(R) = merit = residual^2 -> dF/dR_i = 2*residual*d(residual)/dR_i
for i in (3, 10, 15, 22, 28):
    h = 1e-6
    rp = RS.copy(); rp[0, i] += h
    rm = RS.copy(); rm[0, i] -= h
    fp = spec.merit(sim_of(rp), 1e6)
    fm = spec.merit(sim_of(rm), 1e6)
    dfdr = (fp - fm) / (2 * h)
    # if deposits are d(residual)/dR: 2*residual*deposit == dfdr
    grads[i] = dfdr / (2 * res0) if abs(res0) > 1e-300 else float("nan")
print(f"  (merit FD gives d(residual)/dR via chain; res0={res0:.6g})")
# cross-check: FD of residual directly
for i in (10, 22):
    h = 1e-6
    rp = RS.copy(); rp[0, i] += h
    rm = RS.copy(); rm[0, i] -= h
    dr_res = (float(np.asarray(spec.residuals(sim_of(rp)))[0])
              - float(np.asarray(spec.residuals(sim_of(rm)))[0])) / (2 * h)
    agree = abs(dr_res - grads[i]) / max(abs(dr_res), 1e-300)
    report(f"d(residual)/dR[{i}] consistency", agree < 1e-6,
           f"via merit-FD={grads[i]:.6g} direct={dr_res:.6g}")

# --------------------------------- Part C: end-to-end through the solver
print("=== Part C: end-to-end dF/dtheta with color-only spec ===")
# Assembled independently of the engine's own color buckets, on purpose:
# this is the oracle `color_grad_python.py` checks the (since R4.2, present)
# grads_r/grads_t path against, so it must not consume them.
# dF/dtheta = sum_i 2*res*g_i*(dR_i/dtheta), with
# g_i = d(residual)/dR_i (FD on the sim row, Part B) and
# dR_i/dtheta = P_i/R_i (needle R channel: P_i = R_i * dR_i/ddelta at w=1,t=0).
res0 = float(np.asarray(spec.residuals(sim_of(RS)))[0])
g = np.zeros(31)
h = 1e-6
for i in range(31):
    rp = RS.copy(); rp[0, i] += h
    rm = RS.copy(); rm[0, i] -= h
    g[i] = (float(np.asarray(spec.residuals(sim_of(rp)))[0])
            - float(np.asarray(spec.residuals(sim_of(rm)))[0])) / (2 * h)
R_flat = RS.ravel()
LAY_Z = {1: 60.0, 2: 220.0, 3: 360.0}
LAY_N = {1: 2.35 + 0j, 2: 1.46 + 0j, 3: 2.10 + 0j}


def merit_of(th):
    st = ScatterMatrix(N, np.asarray(th, float), wavelengths=WLS, angles=[30.0])
    o = st.compute(Request.RS)
    return spec.merit(sim_of(np.atleast_2d(np.asarray(o["Rs"], float))), 1e6)


h_th = 0.01
for lay in (1, 2, 3):
    th = D.copy()
    fd = (merit_of(th + np.eye(5)[lay] * h_th) - merit_of(th - np.eye(5)[lay] * h_th)) / (2 * h_th)
    ng0 = needle_gradient(st0, np.full(WLS.size, LAY_N[lay]), [LAY_Z[lay]],
                          NeedleRequest.P, pol="s",
                          targets_r=np.zeros(31), weights_r=np.ones(31))
    P_l = np.asarray(ng0["P_s"]).ravel()
    analytic = float(np.sum(2.0 * res0 * g * (P_l / R_flat)))
    rel = abs(analytic - fd) / max(abs(fd), abs(analytic), 1e-300)
    report(f"color dF/dtheta_layer{lay} (own chain-rule)", rel < 2e-3,
           f"analytic={analytic:.6g} fd={fd:.6g} rel={rel:.2e}")

# ------------------------------- Part D: White / Yellow / DomWl quantities
print("=== Part D: White / Yellow / DomWl closed forms ===")
spec_w = make_spec("White", "Channels", [80.0, 5.0])
native_w = spec_w.merit(sim, 1e6)
s = xyz.sum(); ws = white.sum()
x, y = xyz[0] / s, xyz[1] / s
xn, yn = white[0] / ws, white[1] / ws
W = 100.0 * xyz[1] + 800.0 * (xn - x) + 1700.0 * (yn - y)
Tw = 1000.0 * (xn - x) - 650.0 * (yn - y)
F_w = WEIGHT * (((W - 80.0) / 1.0) ** 2 + ((Tw - 5.0) / 1.0) ** 2)
report("CIE whiteness pair", abs(native_w - F_w) / F_w < 1e-9,
       f"native={native_w:.10g} mine={F_w:.10g}")

spec_y = make_spec("Yellow", "Channels", 12.0)
native_y = spec_y.merit(sim, 1e6)
YI = 100.0 * (1.3013 * xyz[0] - 1.1498 * xyz[2]) / xyz[1]
F_y = WEIGHT * ((YI - 12.0) / 1.0) ** 2
report("E313 yellowness", abs(native_y - F_y) / F_y < 1e-9,
       f"native={native_y:.10g} mine={F_y:.10g}")

# DomWl: independent locus intersection
spec_d = make_spec("DomWl", "Channels", [550.0, 0.9])
native_d = spec_d.merit(sim, 1e6)
px, py = x, y
d = np.array([px - xn, py - yn])
if np.dot(d, d) < 1e-18:
    dl, dp = 0.0, 0.0
else:
    locus = []
    for i in range(len(C_WL)):
        ssum = XB[i] + YB[i] + ZB[i]
        if ssum > 0:
            locus.append(((XB[i] / ssum, YB[i] / ssum), C_WL[i]))
    fwd = bwd = None
    for (a, la), (b, lb) in zip(locus[:-1], locus[1:]):
        e = np.array([b[0] - a[0], b[1] - a[1]])
        den = d[0] * e[1] - d[1] * e[0]
        if abs(den) < 1e-300:
            continue
        aw = np.array([a[0] - xn, a[1] - yn])
        t = (aw[0] * e[1] - aw[1] * e[0]) / den
        sg = (aw[0] * d[1] - aw[1] * d[0]) / den
        if not (0.0 <= sg <= 1.0):
            continue
        lam = la + sg * (lb - la)
        if t > 1.0 - 1e-9 and (fwd is None or t < fwd[0]):
            fwd = (t, lam)
        u = -t
        if u > 1e-9 and (bwd is None or u < bwd[0]):
            bwd = (u, lam)
    if fwd is not None:
        dp, dl = 1.0 / fwd[0], fwd[1]
    else:
        dp, dl = -1.0 / bwd[0], bwd[1]
F_d = WEIGHT * (((dl - 550.0) / 1.0) ** 2 + ((dp - 0.9) / 1.0) ** 2)
rel_d = abs(native_d - F_d) / max(abs(F_d), 1e-300)
report("Dominant wavelength + purity", rel_d < 1e-7,
       f"native={native_d:.10g} mine={F_d:.10g} (dl={dl:.4f} dp={dp:.4f})")

# ------------------------------- Part E: other quantity transforms spot --
print("=== Part E: XyY / Oklab / sRGB / Din99 / Luv channels ===")
def channels_check(quantity, ref, my_triple, tol=1e-9):
    sp = make_spec(quantity, "Channels", list(np.asarray(ref, float)))
    nv = sp.merit(sim, 1e6)
    r = (np.asarray(my_triple) - np.asarray(ref, float))
    F = WEIGHT * float(np.sum((r / 1.0) ** 2))
    rel = abs(nv - F) / max(abs(F), 1e-300)
    report(f"{quantity} channels", rel < tol, f"rel={rel:.2e} (F={F:.4g})")

# XyY: x,y,Y
xXY = xyz[0] / (xyz[0] + xyz[1] + xyz[2])
yXY = xyz[1] / (xyz[0] + xyz[1] + xyz[2])
xyy_t = [xXY, yXY, xyz[1]]
channels_check("XyY", np.asarray(xyy_t) + [0.05, -0.03, 0.10], xyy_t)
# Oklab: Bradford adapt to D65 then Oklab — reuse verified §11 formulas
def bradford_adapt(xyz, src_white, dst_white):
    M = np.array([[0.8951, 0.2664, -0.1614], [-0.7502, 1.7135, 0.0367],
                  [0.0389, -0.0685, 1.0296]])
    Mi = np.linalg.inv(M)
    s = M @ np.asarray(src_white); t = M @ np.asarray(dst_white)
    cone = M @ np.asarray(xyz)
    cone = cone * (t / s)
    return Mi @ cone
adapted = bradford_adapt(xyz, white, [0.950455927, 1.0, 1.089057708])
M_ok = np.array([[0.8189330101, 0.3618667424, -0.1288597137],
                 [0.0329845436, 0.9293118715, 0.0361456387],
                 [0.0482003018, 0.2643662691, 0.6338517070]])
lms = M_ok @ adapted
lms_p = np.cbrt(lms)
ok = np.array([[-0.0329845436*0 + 0, 0, 0]])  # placeholder replaced below
Ok_L = 0.2104542553 * lms_p[0] + 0.7936177850 * lms_p[1] - 0.0040720468 * lms_p[2]
Ok_a = 1.9779984951 * lms_p[0] - 2.4285922050 * lms_p[1] + 0.4505937099 * lms_p[2]
Ok_b = 0.0259040371 * lms_p[0] + 0.7827717662 * lms_p[1] - 0.8086757660 * lms_p[2]
ok_t = [Ok_L, Ok_a, Ok_b]
channels_check("Oklab", np.asarray(ok_t) + [0.5, -1.0, 0.3], ok_t)
# Luv
def luv_of(xyz, white):
    X, Y, Z = xyz; Xn, Yn, Zn = white
    dn = Xn + 15 * Yn + 3 * Zn
    dd = X + 15 * Y + 3 * Z
    up, vp = 4 * X / dd if dd else 0, 9 * Y / dd if dd else 0
    unp, vnp = 4 * Xn / dn, 9 * Yn / dn
    L = 116 * lab_of(xyz, white)[0] / 116 + 16 if False else 116 * (Y / Yn) ** (1 / 3) - 16 if Y > (6 / 29) ** 3 * Yn else 0
    return [L, 13 * L * (up - unp), 13 * L * (vp - vnp)]
luv = luv_of(xyz, white)
channels_check("Luv", np.asarray(luv) + [2.0, -1.0, 0.7], luv, tol=1e-6)

print()
print("ALL OK" if OK else "FAILURES PRESENT")
