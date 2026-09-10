# -*- coding: utf-8 -*-
"""Cross-path consistency for Nevot-Croce roughness (type 5).

The type-5 interface factor is built in FOUR places in the Rust core:

  * ``coherent_block::solve_pol_specialized``          (single-pol solve)
  * ``coherent_block::solve_coherent_block_fields_dual`` (dual-pol solve)
  * ``solver::field_prof``                             (field profile)
  * ``needle_operator::interface_matrix``              (needle sensitivity)

R1.1 originally fixed only the first two -- the remediation plan's site
inventory said "two verbatim sites". The result was worse than the original
bug in one respect: the synthesis *merit* (via ``coherent_block``) and the
synthesis *gradient* (via ``needle_operator``) started using different
interface physics, so the needle P-function was no longer a derivative of the
function it was optimizing.

These tests make a partial fix impossible to land again, three ways:

1. ``test_index_matched_interface_is_invisible`` -- a physical invariant every
   path must satisfy, checked bit-exactly (see its docstring).
2. ``test_needle_gradient_matches_finite_difference`` -- ties the needle
   operator to the coherent-block solver numerically on a genuinely rough
   stack.
3. ``test_all_rtype5_sites_use_shared_factor`` -- a source-level guard: every
   ``rtype == 5`` branch must call the single shared helper.
"""

import re
from pathlib import Path

import numpy as np
import pytest

from navette.smatrix.smatrix import ScatterMatrix, Request
from navette.smatrix.needle import NeedleRequest, needle_gradient

_WLS = np.linspace(450.0, 750.0, 11)
_ANGLES = [0.0, 30.0, 60.0]
_SIGMAS = (1.0, 5.0, 20.0)


# --------------------------------------------------------------------------
# 1. Index-matched interface must be invisible -- on every path
# --------------------------------------------------------------------------

# air | TiO2 120 | TiO2 200 | glass. Interface 1|2 separates two IDENTICAL
# media, and only that interface carries roughness.
_MATCHED_N = np.array([1.0 + 0j, 2.35 + 0j, 2.35 + 0j, 1.52 + 0j])
_MATCHED_D = np.array([0.0, 120.0, 200.0, 0.0])


def _matched_stack(rtype, sigma):
    rt = np.zeros(4, dtype=np.int32)
    rv = np.zeros(4)
    rt[2] = int(rtype)      # the index-matched interface, and only it
    rv[2] = float(sigma)
    return ScatterMatrix(_MATCHED_N, _MATCHED_D, wavelengths=_WLS,
                         angles=_ANGLES, roughness_types=rt,
                         roughness_values=rv)


def _worst_diff(a, b):
    return max(float(np.max(np.abs(np.asarray(a[k]) - np.asarray(b[k]))))
               for k in a)


@pytest.mark.parametrize("sigma", _SIGMAS)
def test_index_matched_interface_is_invisible_compute(sigma):
    """A rough interface between two identical media must not exist optically.

    When n1 == n2 the Fresnel coefficients are r = 0, t = 1, and the
    Nevot-Croce factors collapse to f = exp(-2 kz^2 sigma^2) (multiplying a
    zero reflection, so irrelevant) and ga = exp(0) = 1 exactly. So a type-5
    interface here must be bit-identical to no interface at all, for any
    sigma.

    The historical bug applied the *reflection* factor f to transmission,
    which made this fictitious interface absorb: at sigma = 20 nm it damps
    transmission by a large factor. This is therefore a zero-tolerance probe
    that fires on any site still using the wrong factor.
    """
    req = Request.PHOTOMETRY | Request.ABSORPTION
    ref = _matched_stack(0, 0.0).compute(req, squeeze=False)
    rough = _matched_stack(5, sigma).compute(req, squeeze=False)
    assert _worst_diff(ref, rough) == 0.0


@pytest.mark.parametrize("sigma", _SIGMAS)
def test_index_matched_interface_is_invisible_field_profile(sigma):
    """Same invariant, through ``solver::field_prof``."""
    kw = dict(pol=0, wavelength=600.0, points_per_layer=60)
    ref = _matched_stack(0, 0.0).field_profile(1.7 + 0.01j, **kw)["E"]
    rough = _matched_stack(5, sigma).field_profile(1.7 + 0.01j, **kw)["E"]
    assert np.array_equal(np.asarray(ref), np.asarray(rough))


@pytest.mark.parametrize("sigma", _SIGMAS)
def test_index_matched_interface_is_invisible_needle(sigma):
    """Same invariant, through ``needle_operator::interface_matrix``."""
    nn = np.full(_WLS.size, 1.46 + 0j, dtype=np.complex128)
    kw = dict(pol="s", targets_r=0.0, weights_r=1.0,
              targets_t=0.0, weights_t=1.0)
    req = NeedleRequest.P | NeedleRequest.P_T
    ref = needle_gradient(_matched_stack(0, 0.0), nn, [60.0, 220.0], req, **kw)
    rough = needle_gradient(_matched_stack(5, sigma), nn, [60.0, 220.0], req, **kw)
    assert _worst_diff(ref, rough) == 0.0


# --------------------------------------------------------------------------
# 2. Needle operator vs coherent-block solver, on a genuinely rough stack
# --------------------------------------------------------------------------

_FD_N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
_FD_D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
_FD_HOSTS = ((1, 60.0), (2, 220.0), (3, 360.0))   # (layer, absolute depth)


def _fd_stack(d, rtype, sigma):
    rt = np.full(5, int(rtype), dtype=np.int32)
    rv = np.full(5, float(sigma))
    rt[0] = 0
    rv[0] = 0.0
    return ScatterMatrix(_FD_N, np.asarray(d, float), wavelengths=_WLS,
                         angles=[0.0, 30.0], roughness_types=rt,
                         roughness_values=rv)


def _merit_t(d, rtype, sigma):
    """sum(T^2) -- the merit whose needle sensitivity is ``P_T`` at
    ``targets_t=0, weights_t=1``."""
    out = _fd_stack(d, rtype, sigma).compute(Request.TS, squeeze=False)
    t = np.asarray(out["Ts"], float).ravel()
    return float(np.sum(t * t))


@pytest.mark.parametrize("rtype,sigma", [(0, 0.0), (5, 2.0), (5, 4.0), (5, 8.0)])
def test_needle_gradient_matches_finite_difference(rtype, sigma):
    """``needle_operator``'s P-function must be the derivative of the merit
    that ``coherent_block`` evaluates.

    Inserting a needle of the host's own index at depth z is equivalent to
    growing the host layer, so the analytic sensitivity ``2 * P_T`` must equal
    a central difference of sum(T^2) through the ordinary solver. The two go
    through *different* interface builders, so any disagreement in the type-5
    factor shows up here directly.

    The residual is set by the finite-difference step (~1e-7 relative), not by
    the physics. Before the four-site fix this test failed at rtype 5 by ~20%.
    """
    st = _fd_stack(_FD_D, rtype, sigma)
    h = 1e-3
    for lay, z in _FD_HOSTS:
        d_plus = _FD_D.copy()
        d_plus[lay] += h
        d_minus = _FD_D.copy()
        d_minus[lay] -= h
        fd = (_merit_t(d_plus, rtype, sigma)
              - _merit_t(d_minus, rtype, sigma)) / (2.0 * h)

        nn = np.full(_WLS.size, _FD_N[lay], dtype=np.complex128)
        ng = needle_gradient(st, nn, [z], NeedleRequest.P_T, pol="s",
                             targets_t=0.0, weights_t=1.0)
        analytic = 2.0 * float(np.asarray(ng["P_T_s"]).ravel().sum())

        rel = abs(analytic - fd) / max(abs(fd), 1e-30)
        assert rel < 1e-6, (
            f"layer {lay}: analytic={analytic:.12g} fd={fd:.12g} rel={rel:.3e}")


# --------------------------------------------------------------------------
# 3. Source-level guard: one implementation, no open-coded copies
# --------------------------------------------------------------------------

_RUST_SRC = Path(__file__).resolve().parents[2] / "rust" / "navette" / "src"
_SHARED = "nevot_croce_factors"
_OWNER = "optics_core.rs"


def _rust_files():
    if not _RUST_SRC.is_dir():
        pytest.skip("Rust sources not present (installed wheel, not a checkout)")
    return sorted(_RUST_SRC.rglob("*.rs"))


def test_all_rtype5_sites_use_shared_factor():
    """Every ``rtype == 5`` branch must delegate to the shared helper.

    This is the guard that R1.1 lacked: the fix landed in two of four sites
    because nothing tied them together. Now they all call one function, and
    this test fails if a fifth site appears -- or an existing one drifts back
    to open-coding the exponentials.
    """
    offenders = []
    for path in _rust_files():
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if re.search(r"rtype\s*==\s*5", line):
                window = "\n".join(lines[i:i + 12])
                if _SHARED not in window:
                    offenders.append(f"{path.name}:{i + 1}")
    assert not offenders, (
        f"rtype == 5 branch(es) not using {_SHARED}(): {offenders}")


def test_nevot_croce_factor_has_a_single_implementation():
    """The reflection exponential may be written out in exactly one file."""
    pattern = re.compile(r"-\s*2\.0\s*\*\s*kz1\s*\*\s*kz2")
    owners = [p.name for p in _rust_files() if pattern.search(
        p.read_text(encoding="utf-8"))]
    assert owners == [_OWNER], (
        f"Nevot-Croce reflection factor open-coded outside {_OWNER}: {owners}")
