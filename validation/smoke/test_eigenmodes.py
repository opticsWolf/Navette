# -*- coding: utf-8 -*-
"""Eigenmode search: the pinned SPP, and the runaway that used to pass (R3.3).

``char_func`` is ``|1/r(n_eff)|^2``, so a pole of the reflection coefficient
drives it to zero and the minimizer's job is to find one. The trouble is that
``|1/r|^2`` *also* goes to zero as ``|n_eff|`` runs away, and the search was
unbounded. Seeded on this stack where no s-polarized mode exists, it sprinted
to ``n_eff = -2.0e8 + 8.4e6j`` and reported a characteristic value of
``1.3e-217`` -- twelve orders "better" than the real surface-plasmon pole it
never found, for a trial index two hundred million times the largest index in
the stack. Nothing said anything.

Two changes, both tested here: the minimizer is confined to a physical box
(``3 x max |n|``), and ``refine_mode`` refuses to call something a mode when
its characteristic value says otherwise.

Nothing pinned the eigenmode path before this file, which is why the runaway
survived a full review cycle.
"""

import numpy as np
import pytest

from navette.smatrix.smatrix import Pol, ScatterMatrix

# Air / 50 nm metal / glass -- the textbook Kretschmann geometry. eps = -12+1j
# is a lossy Drude-like metal in the red; at 632.8 nm the glass-side interface
# carries a surface-plasmon pole and there is no s-polarized mode at all.
_EPS_METAL = -12.0 + 1.0j
_N = np.array([1.0 + 0j, np.sqrt(_EPS_METAL), 1.52 + 0j])
_D = np.array([0.0, 50.0, 0.0])
_LAM = 632.8
_FILM_BOTTOM = 50.0

# 3 x max|n| -- the box the minimizer is confined to (see optimizer::n_eff_bound).
_BOUND = 3.0 * max(abs(complex(n)) for n in _N)

# The pole, as this code converges to it from any nearby seed.
_SPP = 1.7139416834 + 0.0226069147j


def _stack():
    return ScatterMatrix(_N, _D, wavelengths=np.array([_LAM]), angles=[0.0])


# --------------------------------------------------------------------------
# 1. The mode that is really there
# --------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [1.70 + 0.02j, 1.75 + 0.03j, 1.72 + 0.025j],
                         ids=["below", "above", "near"])
def test_glass_side_spp_is_found_from_any_nearby_seed(seed):
    """A true pole polishes to a characteristic value near machine zero.

    ``1e-17`` here against the ``1e-6`` rejection threshold is eleven orders of
    margin, which is what makes the threshold safe to leave strict by default:
    a real mode is not close to it.
    """
    n_eff, val = _stack().refine_mode(seed, pol=Pol.P)
    assert abs(n_eff - _SPP) < 1e-6, f"converged to {n_eff}, expected {_SPP}"
    assert val < 1e-10, f"characteristic value {val:.3e} is not a pole"


def test_the_spp_field_peaks_at_the_metal_glass_interface():
    """The physics check the characteristic value cannot give you.

    A surface plasmon is bound *to an interface*: the field decays away from
    it on both sides. Converging to the right number for the wrong reason
    would leave a profile that peaks somewhere else, or does not decay at all.
    """
    st = _stack()
    n_eff, _ = st.refine_mode(1.70 + 0.02j, pol=Pol.P)
    prof = st.field_profile(n_eff, pol=Pol.P, wavelength=_LAM,
                            points_per_layer=200)
    z = np.asarray(prof["z"], float)
    e = np.asarray(prof["E"], float)

    assert z[int(np.argmax(e))] == pytest.approx(_FILM_BOTTOM, abs=1e-9)
    assert np.all(np.diff(e) > 0), "field is not monotonic across the film"
    # Evanescent through the metal: an order-of-magnitude drop, not a ripple.
    assert e[0] / e[-1] < 0.3


def test_find_eigenmodes_recovers_the_same_pole():
    """The landscape -> refine workflow must agree with the manual seed.

    ``find_eigenmodes`` refines inside Rust and does not go through
    ``refine_mode``'s residual contract, so this also checks the box did not
    change what the documented path returns.
    """
    modes = _stack().find_eigenmodes((1.0, 2.0), (0.0, 0.3),
                                     resolution=(240, 120), pol=Pol.P)
    assert any(abs(m - _SPP) < 1e-6 for m in modes), f"SPP missing from {modes}"


# --------------------------------------------------------------------------
# 2. The mode that is not there
# --------------------------------------------------------------------------

_NO_MODE_SEEDS = [1.1 + 0.05j, 1.3 + 0.1j, 1.45 + 0.2j]


@pytest.mark.parametrize("seed", _NO_MODE_SEEDS, ids=lambda s: f"{s.real:g}")
def test_s_polarization_has_no_mode_and_says_so(seed):
    """There is no s-polarized mode on this stack; asking for one must fail."""
    with pytest.raises(ValueError, match="no eigenmode near seed"):
        _stack().refine_mode(seed, pol=Pol.S)


@pytest.mark.parametrize("seed", _NO_MODE_SEEDS + [1e6 + 1e6j, -1e8 + 1e7j],
                         ids=lambda s: f"{s.real:g}")
def test_the_search_never_leaves_the_physical_box(seed):
    """The runaway regression, including from seeds already outside the box.

    ``max_residual=None`` is what the caller used to get unconditionally.
    Whatever comes back must still be a physically conceivable trial index:
    the previous behaviour returned ``-2.0e8 + 8.4e6j`` on the first of these
    seeds, which is 2e7 times the bound.
    """
    n_eff, _ = _stack().refine_mode(seed, pol=Pol.S, max_residual=None)
    assert abs(n_eff.real) <= _BOUND + 1e-9, f"{n_eff} escaped the box"
    assert abs(n_eff.imag) <= _BOUND + 1e-9, f"{n_eff} escaped the box"


def test_the_no_mode_result_is_reported_as_a_bad_residual_not_a_good_one():
    """Boxing alone is not enough: the answer must also look wrong.

    Confined to the box, the search settles on the boundary with a
    characteristic value around 1e-3 -- nine orders worse than the real pole,
    and above the threshold. If the box had merely capped a *good*-looking
    value, the contract would still pass it through.
    """
    _, val = _stack().refine_mode(1.1 + 0.05j, pol=Pol.S, max_residual=None)
    assert val > 1e-6, f"a non-mode scored {val:.3e}, which reads as a pole"


def test_s_polarized_landscape_is_flat():
    """No minimum to seed from, and ``find_eigenmodes`` reports none.

    The documented workflow was never the broken one: a bounded scan finds
    nothing and says nothing was found. It is the manual seed that lied.
    """
    st = _stack()
    land = st.eigenmode_landscape((1.0, 2.0), (0.0, 0.3), resolution=(60, 60),
                                  pol=Pol.S)
    assert float(land.values.min()) == pytest.approx(1.0, abs=1e-6)
    assert st.find_eigenmodes((1.0, 2.0), (0.0, 0.3), resolution=(60, 60),
                              pol=Pol.S) == []


# --------------------------------------------------------------------------
# 3. The contract itself
# --------------------------------------------------------------------------

def test_max_residual_none_returns_the_raw_result():
    n_eff, val = _stack().refine_mode(1.1 + 0.05j, pol=Pol.S, max_residual=None)
    assert np.isfinite(val) and np.isfinite(n_eff.real)


def test_a_loose_threshold_accepts_what_the_strict_one_rejects():
    """The threshold is a knob, not a hidden constant."""
    st = _stack()
    with pytest.raises(ValueError):
        st.refine_mode(1.1 + 0.05j, pol=Pol.S)
    n_eff, val = st.refine_mode(1.1 + 0.05j, pol=Pol.S, max_residual=1.0)
    assert val <= 1.0


def test_the_rejection_message_carries_the_evidence():
    """It must say where the search went and how bad it was, not just "failed"."""
    with pytest.raises(ValueError) as excinfo:
        _stack().refine_mode(1.1 + 0.05j, pol=Pol.S)
    message = str(excinfo.value)
    for fragment in ("n_eff", "characteristic value", "max_residual",
                     "eigenmode_landscape", "max_residual=None"):
        assert fragment in message, f"message omits {fragment!r}: {message}"


# --------------------------------------------------------------------------
# 4. The landscape is deliberately NOT boxed
# --------------------------------------------------------------------------

def test_the_landscape_is_not_walled_off_outside_the_box():
    """The box belongs to the minimizer, not to the scanner.

    The remediation plan put the guard in ``char_func``, which both share.
    That would paint ``1e30`` across any part of a user-requested scan range
    lying outside the box -- corrupting a diagnostic the caller explicitly
    asked for, at whatever range they asked for it. The guard lives in
    ``char_func_xy`` instead: the minimizer picks its own points and needs
    walls; the scanner is already bounded by the caller.
    """
    land = _stack().eigenmode_landscape((5.0, 20.0), (0.0, 1.0),
                                        resolution=(40, 10), pol=Pol.P)
    values = np.asarray(land.values, float)
    assert land.n_real.max() > _BOUND, "scan did not reach outside the box"
    assert np.all(np.isfinite(values))
    assert values.max() < 1e6, f"landscape shows a wall: max {values.max():.3e}"
