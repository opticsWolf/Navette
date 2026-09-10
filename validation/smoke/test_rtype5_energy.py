"""Roughness energy-conservation regression (R1.1, review §3.2).

The Névot-Croce path (roughness type 5) historically applied the reflection
Debye-Waller factor ``f = exp(-2 kz1 kz2 sigma^2)`` to the transmission
amplitudes as well, silently destroying up to ~7.5% of transmitted energy at a
single near-matched interface. The fix applies the wavevector-difference factor
``ga = exp(+((kz1-kz2) sigma)^2/2)`` to transmission instead.

Because ``ga >= 1`` compensates the reflection damping, the corrected build
gives ``R + T >= 1`` for every rtype-5 case (a tiny, one-sided overshoot that is
the perturbative NC amplitude-factor model's residual, largest at high-contrast
grazing p-pol where the factoring approximation is weakest). The historical bug
gave ``R + T < 1`` (a loss) — so the direction of the residual alone is a strong
regression guard.

Types 1-4 are graded-profile approximations that intentionally lose a little
energy; they are checked separately as a bounded, one-sided loss.
"""

import numpy as np
import pytest

from navette.smatrix.smatrix import ScatterMatrix, Request

# Single interface: ambient (n=1) -> substrate (real n, so lossless: energy has
# nowhere to go and R+T must equal 1 up to the roughness model's fidelity).
_LAM = np.array([550.0])
_REQ = Request.RS | Request.RP | Request.TS | Request.TP
_CONTRASTS = (1.5001, 2.35, 4.28)   # near-matched -> high-contrast
_SIGMAS = (0.0, 3.0, 6.0, 10.0)
_ANGLES = (0.0, 30.0, 60.0)


def _energy(nsub, sigma, rtype, theta):
    """Return (R_s+T_s, R_p+T_p) for a single rough interface."""
    n = np.array([1.0 + 0j, complex(nsub, 0.0)])
    d = np.array([0.0, 0.0])
    rt = np.array([0, int(rtype)], dtype=np.int32)
    rv = np.array([0.0, float(sigma)])
    out = ScatterMatrix(
        n, d, wavelengths=_LAM, angles=[theta],
        roughness_types=rt, roughness_values=rv,
    ).compute(_REQ)
    rs, rp, ts, tp = (float(out[k][0]) for k in ("Rs", "Rp", "Ts", "Tp"))
    return rs + ts, rp + tp


@pytest.mark.parametrize("nsub", _CONTRASTS)
@pytest.mark.parametrize("sigma", _SIGMAS)
@pytest.mark.parametrize("theta", _ANGLES)
def test_rtype5_never_damps_transmission(nsub, sigma, theta):
    """The primary regression guard: transmission is no longer damped.

    Fixed build: R+T >= 1 (tiny one-sided overshoot, <= ~15% only at the most
    extreme high-contrast grazing p-pol corner where the NC factor model is
    weakest). The pre-fix bug gave R+T as low as 0.9247 -- restoring the old
    ``t12*f`` line fails the ``>= 1 - 1e-6`` bound by a wide margin.
    """
    es, ep = _energy(nsub, sigma, 5, theta)
    for e in (es, ep):
        assert e >= 1.0 - 1e-6, (
            f"rtype-5 transmission damping regressed: R+T={e} "
            f"(nsub={nsub}, sig={sigma}, th={theta})"
        )
        assert e <= 1.16, f"rtype-5 R+T unreasonably large: {e}"


@pytest.mark.parametrize("nsub", _CONTRASTS)
@pytest.mark.parametrize("sigma", _SIGMAS)
def test_rtype5_normal_incidence_near_conservation(nsub, sigma):
    """At normal incidence the NC factoring is well-behaved: |1-(R+T)| small,
    scaling with contrast and sigma (never the pre-fix ~7.5% loss)."""
    es, ep = _energy(nsub, sigma, 5, 0.0)
    # tolerance grows with contrast (perturbation ~ kz*sigma)
    tol = {1.5001: 5e-4, 2.35: 3e-3, 4.28: 2e-2}[nsub]
    assert abs(es - 1.0) < tol and abs(ep - 1.0) < tol, (
        f"nsub={nsub} sig={sigma}: R+T=({es},{ep}), tol={tol}"
    )


def test_rtype5_worked_example():
    """§3.2 worked example: n 1.0->1.5001, sigma=10 nm, normal incidence.

    Reflection factor damps, transmission factor compensates -> R+T ~= 1.0001.
    The old bug gave 0.924678 here (fails abs<3e-4 by >250x).
    """
    es, ep = _energy(1.5001, 10.0, 5, 0.0)
    assert abs(es - 1.0) < 3e-4, f"R+T (s) = {es}"
    assert abs(ep - 1.0) < 3e-4, f"R+T (p) = {ep}"


@pytest.mark.parametrize("nsub", _CONTRASTS)
@pytest.mark.parametrize("sigma", (3.0, 6.0, 10.0))
def test_rtype5_matches_nevot_croce_closed_form(nsub, sigma):
    """Independent literature anchor: engine R and T (not just their sum) match
    the closed-form Névot-Croce single-interface prediction to machine
    precision.

    Névot & Croce, Rev. Phys. Appl. 15, 761 (1980); de Boer, Phys. Rev. B 44,
    498 (1991); Stearns, J. Appl. Phys. 65, 491 (1989):

        r~ = r0 * exp(-2 kz1 kz2 sigma^2)        =>  R = R0 * exp(-4 kz1 kz2 s^2)
        t~ = t0 * exp(+(kz1-kz2)^2 sigma^2 / 2)  =>  T = T0 * exp(+(kz1-kz2)^2 s^2)

    with ideal Fresnel R0 = ((n1-n2)/(n1+n2))^2, T0 = 1 - R0 (normal incidence,
    lossless). This checks the fix against physics, not against its own output.
    """
    n1, lam = 1.0, float(_LAM[0])
    kz1 = 2.0 * np.pi * n1 / lam
    kz2 = 2.0 * np.pi * nsub / lam
    r0 = (n1 - nsub) / (n1 + nsub)
    R0 = r0 * r0
    T0 = 1.0 - R0
    R_lit = R0 * np.exp(-4.0 * kz1 * kz2 * sigma * sigma)
    T_lit = T0 * np.exp((kz1 - kz2) ** 2 * sigma * sigma)

    n = np.array([1.0 + 0j, complex(nsub, 0.0)])
    d = np.array([0.0, 0.0])
    out = ScatterMatrix(
        n, d, wavelengths=_LAM, angles=[0.0],
        roughness_types=np.array([0, 5], dtype=np.int32),
        roughness_values=np.array([0.0, sigma]),
    ).compute(Request.RS | Request.TS)
    R_eng, T_eng = float(out["Rs"][0]), float(out["Ts"][0])
    assert abs(R_eng - R_lit) < 1e-12, f"R engine={R_eng} lit={R_lit}"
    assert abs(T_eng - T_lit) < 1e-12, f"T engine={T_eng} lit={T_lit}"


def test_rtype0_strict_conservation():
    """No roughness -> exact energy conservation."""
    for nsub in _CONTRASTS:
        for th in _ANGLES:
            es, ep = _energy(nsub, 10.0, 0, th)
            assert abs(es - 1.0) < 1e-9 and abs(ep - 1.0) < 1e-9


@pytest.mark.parametrize("rtype", (1, 2, 3, 4))
def test_types_1_4_bounded_graded_loss(rtype):
    """Graded-profile roughness (types 1-4) loses a little energy by design:
    the loss is one-sided (R+T <= 1) and bounded -- never a gain and never the
    catastrophic NC-transmission-damping regression."""
    for nsub in _CONTRASTS:
        for sig in _SIGMAS:
            for th in _ANGLES:
                es, ep = _energy(nsub, sig, rtype, th)
                for e in (es, ep):
                    assert e <= 1.0 + 1e-9, f"type {rtype} gained energy: {e}"
                    assert e >= 0.84, f"type {rtype} lost too much energy: {e}"


def test_type4_near_matched_documented_value():
    """Type-4 (Gaussian w-function) near-matched sigma=10 stays ~0.9948 (§3.2)."""
    es, _ = _energy(1.5001, 10.0, 4, 0.0)
    assert abs(es - 0.994837) < 1e-3, f"type-4 R+T = {es}"
