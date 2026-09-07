# Step 5: Tauc-Lorentz golden from first principles (break the reference oracle).
# See docs/code_review.md section 23. kk_validate.py provides the
# pair-sampled principal-value Kramers-Kronig scheme (validated against an
# analytic Lorentz oscillator first — see that file's main block).
import numpy as np
from navette.materials import MaterialSpec, evaluate
from kk_validate import eps1_pair

OK = True
def report(name, cond, detail=""):
    global OK
    print(f"  {name}: {'OK' if cond else 'FAIL'} {detail}")
    OK = bool(cond) and OK

HC_EV_NM = 1239.841984

# ---- my own TL eps2 (Jellison-Modine 1996 closed form) --------------------
def eps2_tl(E, Eg, oscs):
    E = np.asarray(E, float)
    out = np.zeros_like(E)
    for (A, E0, C) in oscs:
        mask = E > Eg
        out[mask] += (A * E0 * C * (E[mask] - Eg) ** 2
                      / (((E[mask] ** 2 - E0 ** 2) ** 2 + C ** 2 * E[mask] ** 2) * E[mask]))
    return out

# ---- my own Lab + Sharma DeltaE2000 (independent of func_04/func_16) ------
def lab_of(xyz, white):
    def f(t):
        d = 6.0 / 29.0
        return np.where(t > d ** 3, np.cbrt(t), t / (3 * d ** 2) + 4.0 / 29.0)
    x, y, z = np.asarray(xyz) / np.asarray(white)
    fx, fy, fz = f(x), f(y), f(z)
    return np.array([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)])


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


def eps1_kk_vec(E_t, Eg, oscs, eps_inf, E_hi=80.0, h=2e-3):
    """Independent KK via the corrected pair-sampled PV scheme."""
    return np.array([eps1_pair(eps2_tl, (Eg, oscs), E, eps_inf=eps_inf,
                               E_hi=E_hi, h=h) for E in np.atleast_1d(E_t)])


# ---- compare against navette ----------------------------------------------
spec = MaterialSpec("TaucLorentz", {
    "Eg": 3.0,
    "osc": [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)],
    "epsilon_inf": 1.0,
})
wls = np.array([300.0, 400.0, 550.0, 700.0, 900.0, 1200.0])
nk = evaluate(spec, wls)
E = HC_EV_NM / wls

# eps2 check (exact closed form; 550nm+ sit below the 3 eV gap -> eps2 == 0)
e2_mine = eps2_tl(E, 3.0, [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)])
e2_engine = 2.0 * np.real(nk) * np.imag(nk)
err2 = np.abs(e2_engine - e2_mine)
scale2 = np.maximum(np.abs(e2_mine), 1e-3)   # abs-tol floor for below-gap pts
rel2 = err2 / scale2
# NOTE: 1e-6 not 1e-12 — HC_EV_NM here is 1239.841984 while the engine's
# `energy_ev` carries more digits; the ~3e-10 energy difference shows up as
# a ~1e-8 relative eps2 difference exactly on the steep E0=4 eV resonance.
report("eps2 closed form (Jellison-Modine)", rel2.max() < 1e-6,
       f"max scaled err={rel2.max():.2e}  (below-gap pts e2=0 exact both; "
       "residual = hc-constant digit truncation on the sharp resonance)")

# eps1 check via independent PV-KK quadrature
eps1_engine = np.real(nk) ** 2 - np.imag(nk) ** 2
eps1_mine = eps1_kk_vec(E, 3.0, [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)], 1.0)
rel1 = np.abs(eps1_engine - eps1_mine) / np.maximum(np.abs(eps1_mine), 1e-300)
report("eps1 via independent pair-sampled PV-KK", rel1.max() < 2e-2,
       f"max rel={rel1.max():.2e} (residual = documented FFT-KK grid accuracy "
       "near the E0=4eV resonance; my quadrature is refinement-stable)")
print("    eps1 engine:", np.round(eps1_engine, 5))
print("    eps1 mine  :", np.round(eps1_mine, 5))

# below-gap: eps2 must be exactly 0
wls_bg = np.array([1300.0, 1500.0])
nk_bg = evaluate(spec, wls_bg)
e2_bg = 2.0 * np.real(nk_bg) * np.imag(nk_bg)
report("eps2 == 0 below the gap", np.abs(e2_bg).max() < 1e-12,
       f"max|e2|={np.abs(e2_bg).max():.2e}")

print()
print("ALL OK" if OK else "FAILURES PRESENT")
