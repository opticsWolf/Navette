# Independent (non-loom) verification of the synthesis gradient chain.
# Part A: Python reimplementation of the merit fold vs native MeritSpec.merit
# Part B: end-to-end dF/dtheta_k: FD through the real solver+fold vs
#         analytic from build_needle_targets + needle_gradient
# Part C: dispersion channels (dphi/dgd/dgdd/dtod/dfod) vs central-difference
#         spectral derivatives across an inserted needle slab
import json
import numpy as np

from navette.smatrix.smatrix import ScatterMatrix, Request
from navette.smatrix.needle import NeedleRequest, needle_gradient
from navette._smatrix import MeritSpec, SimCurves, build_needle_targets

rng = np.random.default_rng(20260214)
OK = True


def report(name, cond, detail=""):
    global OK
    print(f"  {name}: {'OK' if cond else 'FAIL'} {detail}")
    OK = bool(cond) and OK


# ---------------------------------------------------------------- stack
# air | TiO2 120 | SiO2 200 | Ta2O5 80 | glass   (roughness-free)
N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
WLS = np.linspace(450.0, 750.0, 31)
ANGLES = [0.0, 30.0]


def stack_with(d=None):
    dd = D.copy() if d is None else np.asarray(d, float)
    return ScatterMatrix(N, dd, wavelengths=WLS, angles=ANGLES)


def solver_curves(st):
    req = (Request.RS | Request.TS | Request.PHI_RS | Request.RS_C
           | Request.DISP_R_S)
    out = st.compute(req)
    rs = np.asarray(out["Rs"], float)   # (n_angles, n_wavs)
    ts = np.asarray(out["Ts"], float)
    rs_c = np.asarray(out["rs_c"])
    return out


def sim_from(st, total_d):
    out = solver_curves(st)
    sim = SimCurves(np.array(ANGLES, float), WLS, total_d, 1.0, 1.52)
    sim.set_curve("Rs", np.ascontiguousarray(rs.ravel()))
    sim.set_curve("Ts", np.ascontiguousarray(ts.ravel()))
    sim.set_complex("Rs", np.ascontiguousarray(rs_c.ravel().astype(np.complex128)))
    return sim, out


def build_spec():
    spec = MeritSpec()
    wl_e = np.ascontiguousarray(WLS[::5])   # exact target on coarse grid
    key = spec.add_key(30.0, "Rs")
    spec.add_target(key, wl_e, 0.32 * np.ones(wl_e.size),
                    0.05 * np.ones(wl_e.size), "e", "linear", 1.0, None,
                    False, None, 1.7, None, False)
    # above-target on T at 30 deg, log transform
    keyT = spec.add_key(30.0, "Ts")
    spec.add_target(keyT, wl_e, 0.80 * np.ones(wl_e.size),
                    0.10 * np.ones(wl_e.size), "a", "log", 1.0, None,
                    False, None, 1.3, None, False)
    # range band on R at 0 deg
    keyR0 = spec.add_key(0.0, "Rs")
    wls_c = np.ascontiguousarray(WLS)
    spec.add_target(keyR0, wls_c, 0.10 * np.ones(WLS.size),
                    0.02 * np.ones(WLS.size), "r", "linear", 1.0,
                    0.03 * np.ones(WLS.size), False, None, 0.9, None, False)
    # phase (absolute) on r at 30 deg
    spec.add_target(key, wl_e, 0.5 * np.ones(wl_e.size),
                    0.20 * np.ones(wl_e.size), "e", "phase", 1.0, None,
                    True, None, 1.1, None, False)
    return spec


# ------------------------------------------------------- Part A: fold value
print("=== Part A: independent Python fold vs native MeritSpec.merit ===")
TAU = 2.0 * np.pi


def kind_residual(kind, sd, tol, bw):
    if kind == "e":
        return sd / tol
    if kind == "a":
        return sd / tol if sd < 0.0 else 0.0
    if kind == "b":
        return sd / tol if sd > 0.0 else 0.0
    if kind == "r":
        bwe = tol if bw <= 0.0 else bw
        ad = abs(sd)
        return 0.0 if ad <= bwe else (ad - bwe) / tol
    if kind == "c":
        if bw <= 0.0:
            return sd / tol
        ad = abs(sd)
        return sd / bw if ad <= bw else np.sqrt(((ad - bw) / tol) ** 2 + 1.0)
    raise ValueError(kind)


def my_merit(sim_data, spec_desc):
    """spec_desc: list of (angle, curve, wl, targets, tol, kind, transform,
    nf, band, phase, weight, count_norm, integral)."""
    total = 0.0
    for (ang, cur, twl, tgt, tol, kind, tr, nf, band, phase, w, cn, integ) in spec_desc:
        row = sim_data[cur][ANGLES.index(ang)]
        crow = sim_data[cur + "_c"][ANGLES.index(ang)] if cur + "_c" in sim_data else None
        rscale = np.sqrt(w / (cn if cn else 1.0))
        acc_d = acc_tol = acc_bw = 0.0
        res = []
        for i, wv in enumerate(twl):
            k = np.argmin(np.abs(WLS - wv))
            raw = float(row[k])
            if phase:
                raw = float(np.angle(crow[k]))
                sd = (raw - tgt[i] + np.pi) % TAU - np.pi
            elif tr == "log":
                sd = np.log10(max(raw, 1e-12)) * nf - tgt[i]
            else:
                sd = raw * nf - tgt[i]
            b = band[i] if band is not None else 0.0
            if integ:
                acc_d += sd; acc_tol += tol[i]; acc_bw += b
            else:
                res.append(kind_residual(kind, sd, tol[i], b) * rscale)
        if integ:
            n = len(twl)
            res.append(kind_residual(kind, acc_d / n, max(acc_tol / n, 1e-300),
                                     acc_bw / n) * rscale)
        total += float(np.sum(np.asarray(res) ** 2))
    return total


st0 = stack_with()
out0 = solver_curves(st0)
sim0_native = SimCurves(np.array(ANGLES, float), WLS, float(D[1:-1].sum()), 1.0, 1.52)
sim0_native.set_curve("Rs", np.ascontiguousarray(np.asarray(out0["Rs"], float).ravel()))
sim0_native.set_curve("Ts", np.ascontiguousarray(np.asarray(out0["Ts"], float).ravel()))
sim0_native.set_complex("Rs", np.ascontiguousarray(np.asarray(out0["rs_c"]).ravel().astype(np.complex128)))
spec = build_spec()

sim_data = {
    "Rs": np.asarray(out0["Rs"], float),
    "Ts": np.asarray(out0["Ts"], float),
    "Rs_c": np.asarray(out0["rs_c"]),
}
spec_desc = [
    (30.0, "Rs", WLS[::5], 0.32 * np.ones(7), 0.05 * np.ones(7), "e", "lin", 1.0, None, False, 1.7, None, False),
    (30.0, "Ts", WLS[::5], 0.80 * np.ones(7), 0.10 * np.ones(7), "a", "log", 1.0, None, False, 1.3, None, False),
    (0.0, "Rs", WLS, 0.10 * np.ones(31), 0.02 * np.ones(31), "r", "lin", 1.0, 0.03 * np.ones(31), False, 0.9, None, False),
    (30.0, "Rs", WLS[::5], 0.5 * np.ones(7), 0.20 * np.ones(7), "e", "phase", 1.0, None, True, 1.1, None, False),
]
mine = my_merit(sim_data, spec_desc)
native = spec.merit(sim0_native, 1e6)
err = abs(mine - native) / max(abs(native), 1e-300)
report("fold merit value", err < 1e-12, f"mine={mine:.12g} native={native:.12g} rel={err:.2e}")

# --------------------------------------------- Part B: end-to-end FD gradient
print("=== Part B: dF/dtheta (FD via solver) vs analytic needle fold ===")


def build_spec_lin():
    """All-linear spec (needle-foldable): e on Rs@30, a on Ts@30, r on Rs@0."""
    spec = MeritSpec()
    wl_e = np.ascontiguousarray(WLS[::5])
    key30 = spec.add_key(30.0, "Rs")
    spec.add_target(key30, wl_e, 0.32 * np.ones(wl_e.size),
                    0.05 * np.ones(wl_e.size), "e", "linear", 1.0, None,
                    False, None, 1.7, None, False)
    keyT = spec.add_key(30.0, "Ts")
    spec.add_target(keyT, wl_e, 0.80 * np.ones(wl_e.size),
                    0.10 * np.ones(wl_e.size), "a", "linear", 1.0, None,
                    False, None, 1.3, None, False)
    key0 = spec.add_key(0.0, "Rs")
    wls_c = np.ascontiguousarray(WLS)
    spec.add_target(key0, wls_c, 0.10 * np.ones(WLS.size),
                    0.02 * np.ones(WLS.size), "r", "linear", 1.0,
                    0.03 * np.ones(WLS.size), False, None, 0.9, None, False)
    return spec


spec_lin = build_spec_lin()
native_lin = spec_lin.merit(sim0_native, 1e6)
fold = build_needle_targets(spec_lin, np.array(ANGLES, float), WLS, sim0_native)
print(f"  native merit (linear spec) = {native_lin:.10g}")
for ch in ("r", "t", "a", "phi0"):
    tt = np.asarray(fold[ch]["targets"]).ravel()
    ww = np.asarray(fold[ch]["weights"]).ravel()
    act = int((np.abs(ww) > 0).sum())
    print(f"    fold[{ch}]: shape={tt.shape} active_w={act} "
          f"w-range=[{ww.min():.3g},{ww.max():.3g}] t-range=[{tt.min():.3g},{tt.max():.3g}]")


def merit_of_stack(th):
    st = stack_with(th)
    out = st.compute(Request.RS | Request.TS | Request.RS_C)
    sim = SimCurves(np.array(ANGLES, float), WLS, float(np.sum(th[1:-1])), 1.0, 1.52)
    sim.set_curve("Rs", np.ascontiguousarray(np.asarray(out["Rs"], float).ravel()))
    sim.set_curve("Ts", np.ascontiguousarray(np.asarray(out["Ts"], float).ravel()))
    sim.set_complex("Rs", np.ascontiguousarray(np.asarray(out["rs_c"]).ravel().astype(np.complex128)))
    return spec_lin.merit(sim, 1e6)


h = 1e-2
LAYERS = {1: 120.0, 2: 200.0, 3: 80.0}
N_HOST = {1: 2.35 + 0j, 2: 1.46 + 0j, 3: 2.10 + 0j}
Z_HOST = {1: 60.0, 2: 220.0, 3: 360.0}   # abs depths (layer1:0-120, 2:120-320, 3:320-400)

for lay in (1, 2, 3):
    th = D.copy()
    f_plus = merit_of_stack(th + np.eye(5)[lay] * h)
    f_minus = merit_of_stack(th - np.eye(5)[lay] * h)
    fd = (f_plus - f_minus) / (2 * h)

    nn = np.full(WLS.size, N_HOST[lay], dtype=np.complex128)
    req = NeedleRequest.P | NeedleRequest.P_T | NeedleRequest.P_A | NeedleRequest.P_PHI
    ng = needle_gradient(st0, nn, [Z_HOST[lay]], req, pol="s",
                         targets_r=np.asarray(fold["r"]["targets"]).ravel(),
                         weights_r=np.asarray(fold["r"]["weights"]).ravel(),
                         targets_t=np.asarray(fold["t"]["targets"]).ravel(),
                         weights_t=np.asarray(fold["t"]["weights"]).ravel(),
                         targets_a=np.asarray(fold["a"]["targets"]).ravel(),
                         weights_a=np.asarray(fold["a"]["weights"]).ravel(),
                         targets_phi=np.asarray(fold["phi0"]["targets"]).ravel(),
                         weights_phi=np.asarray(fold["phi0"]["weights"]).ravel(),
                         gain_shift_phi=float(np.asarray(fold["phi0"]["gain_shift"]).ravel()[0]),
                         channel=0)
    p_r = np.asarray(ng["P_s"]).ravel()
    p_t = np.asarray(ng["P_T_s"]).ravel()
    p_a = np.asarray(ng["P_A_s"]).ravel()
    p_phi = np.asarray(ng["P_PHI_s"]).ravel()
    analytic = 2.0 * (p_r + p_t + p_a).sum() + p_phi.sum()
    rel = abs(analytic - fd) / max(abs(fd), abs(analytic), 1e-300)
    report(f"dF/dtheta_layer{lay}", rel < 1e-6,
           f"analytic={analytic:.10g} fd={fd:.10g} rel={rel:.2e}")


# --------------------------------------------- Part C: dispersion channels
print("=== Part C: dispersion ladder, independent verification ===")
# Method: (1) dphi vs physical FD of the phase (unwrap-free robust formula,
# positive-only one-sided 2nd-order stencil — negative thicknesses are
# silently clamped by the solver, see review S21); (2) dgd..dfod vs
# independent np.gradient spectral derivatives of the FD-verified dphi.
# (FD of the engine's own GD outputs across delta-stacks is unreliable:
# the internal phase unwrap can flip branches between stacks.)


def stack_insert(xi, n_prime, delta):
    idx = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, complex(n_prime),
                    1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
    thick = np.array([0.0, 120.0, xi, delta, 200.0 - xi, 80.0, 0.0])
    return ScatterMatrix(idx, thick, wavelengths=WLS, angles=[30.0])


def amp_of(st):
    return np.asarray(st.compute(Request.RS_C)["rs_c"])


def fd_dphi(h, nprime=1.8 + 0.05j, xi=100.0, pol_ang=30.0):
    st = ScatterMatrix(N, D, wavelengths=WLS, angles=[pol_ang])
    a0 = amp_of(st)
    ap = amp_of(stack_insert(xi, nprime, h))
    a2 = amp_of(stack_insert(xi, nprime, 2 * h))
    return ((np.conj(a0) * (-3 * a0 + 4 * ap - a2) / (2 * h)).imag
            / np.abs(a0) ** 2)


nn = np.full(WLS.size, 1.8 + 0.05j, dtype=np.complex128)
ng = needle_gradient(st0, nn, [220.0], NeedleRequest.DISPERSION,
                     channel=0, pol="s")
eng = {k: np.asarray(v).reshape(2, WLS.size)[1] for k, v in ng.items()}  # 30deg row

# (1) dphi convergence series
print("  dphi convergence (rel err vs engine):")
prev = None
for h in (8.0, 4.0, 2.0, 1.0, 0.5):
    fd = fd_dphi(h).ravel()
    rel = np.abs(fd - eng["dphi_s"]).max() / max(np.abs(eng["dphi_s"]).max(), 1e-300)
    rate = f"  rate={prev / rel:.1f}" if prev else ""
    print(f"    h={h:5.1f}nm rel={rel:.3e}{rate}")
    prev = rel
fd_final = fd_dphi(0.5).ravel()
scale = np.abs(eng["dphi_s"]).max()
report("dispersion/dphi (h=0.5nm)", np.abs(fd_final - eng["dphi_s"]).max() / scale < 1e-3,
       f"rel={np.abs(fd_final - eng['dphi_s']).max() / scale:.2e}")

# (2) spectral chain: independent np.gradient in omega
C_NM_PER_FS = 299.792458
omega = 2 * np.pi * C_NM_PER_FS / WLS
chain = [("dgd", "dphi_s", "dgd_s"), ("dgdd", "dgd_s", "dgdd_s"),
         ("dtod", "dgdd_s", "dtod_s"), ("dfod", "dtod_s", "dfod_s")]
base = eng["dphi_s"]
for name, base_key, target_key in chain:
    mine = np.gradient(base, omega, axis=-1)
    rel = np.abs(mine - eng[target_key]).max() / max(np.abs(eng[target_key]).max(), 1e-300)
    report(f"dispersion/{name} = d/domega(previous)", rel < 1e-12, f"rel={rel:.2e}")
    base = mine

# (3) p-polarization spot check (dphi only)
ng_p = needle_gradient(st0, nn, [220.0], NeedleRequest.DPHI,
                       channel=0, pol="p")
eng_p = np.asarray(ng_p["dphi_p"]).reshape(2, WLS.size)[1]
fd_p = fd_dphi(1.0)  # amp_of uses RS_C = s-pol; rebuild for p
def fd_dphi_p(h):
    st = ScatterMatrix(N, D, wavelengths=WLS, angles=[30.0])
    a0 = np.asarray(st.compute(Request.RP_C)["rp_c"])
    def amp_p(delta):
        stx = stack_insert(100.0, 1.8 + 0.05j, delta)
        return np.asarray(stx.compute(Request.RP_C)["rp_c"])
    ap, a2 = amp_p(h), amp_p(2 * h)
    return ((np.conj(a0) * (-3 * a0 + 4 * ap - a2) / (2 * h)).imag
            / np.abs(a0) ** 2)
fd_p = fd_dphi_p(1.0).ravel()
rel_p = np.abs(fd_p - eng_p).max() / max(np.abs(eng_p).max(), 1e-300)
report("dispersion/dphi p-pol (h=1nm)", rel_p < 5e-3, f"rel={rel_p:.2e}")

print()
print("ALL OK" if OK else "FAILURES PRESENT")
