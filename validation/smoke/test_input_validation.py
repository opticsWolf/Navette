# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Input validation for the public solver entry points (R3.1, R3.2, R3.4).

The engine is permissive on purpose -- it is also the optimizer's inner loop.
Before R3.1 that permissiveness reached the user unfiltered, and
``validation/review/garbage_in.py`` recorded what came back:

  * a negative thickness and a NaN thickness both produced *the same numbers
    as deleting the layer* -- a plausible answer to a question nobody asked;
  * 120 deg aliased onto 60 deg, because only ``sin(theta)`` reaches the
    engine;
  * a duplicated wavelength made every dispersion channel NaN (the GD/GDD
    kernels divide by the grid spacing);
  * a NaN index turned the whole output NaN with nothing naming the layer.

None raised. None warned. The needle operator had the same hole one layer
down, guarded by a ``debug_assert!`` that no release build compiles in (R3.2,
bottom section).

These tests hold the line on both sides of it: the garbage must raise, with
the offending index and value in the message, and the legitimately odd inputs
must still go through untouched.
"""

import numpy as np
import pytest

from navette.smatrix.smatrix import CoherenceMode, Request, ScatterMatrix

_N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 1.52 + 0j])
_D = np.array([0.0, 120.0, 200.0, 0.0])
_WLS = np.linspace(450.0, 750.0, 31)
_ANG = [30.0]


def _sm(**kw):
    args = dict(layer_indices=_N, thicknesses=_D, wavelengths=_WLS, angles=_ANG)
    args.update(kw)
    layer_indices = args.pop("layer_indices")
    thicknesses = args.pop("thicknesses")
    return ScatterMatrix(layer_indices, thicknesses, **args)


def _with(arr, idx, value):
    out = np.array(arr, copy=True)
    out[idx] = value
    return out


# --------------------------------------------------------------------------
# Rejected inputs -- one row per silent-wrong-answer path
# --------------------------------------------------------------------------

_REJECT = [
    # (id, kwargs, fragments the message must contain)
    ("nan_index", dict(layer_indices=_with(_N, 1, np.nan)), ["layer_indices", "layer 1", "nan"]),
    ("inf_index", dict(layer_indices=_with(_N, 2, np.inf)), ["layer_indices", "layer 2", "inf"]),
    ("nan_imag_index", dict(layer_indices=_with(_N, 1, 2.0 + np.nan * 1j)), ["layer_indices", "layer 1"]),
    ("overflow_index", dict(layer_indices=_with(_N, 1, 1e308 + 0j)), ["layer_indices", "layer 1", "magnitude"]),
    ("nan_thickness", dict(thicknesses=_with(_D, 2, np.nan)), ["thicknesses", "layer 2", "nan"]),
    ("inf_thickness", dict(thicknesses=_with(_D, 2, np.inf)), ["thicknesses", "layer 2", "inf"]),
    ("negative_thickness", dict(thicknesses=_with(_D, 2, -50.0)), ["thicknesses", "layer 2", "-50.0"]),
    ("nan_wavelength", dict(wavelengths=_with(_WLS, 5, np.nan)), ["wavelengths", "index 5", "nan"]),
    ("zero_wavelength", dict(wavelengths=_with(_WLS, 0, 0.0)), ["wavelengths", "index 0", "0.0"]),
    ("negative_wavelength", dict(wavelengths=_with(_WLS, 3, -1.0)), ["wavelengths", "index 3", "-1.0"]),
    ("duplicate_wavelength", dict(wavelengths=np.array([500.0, 550.0, 550.0, 600.0])),
     ["wavelengths", "strictly increasing", "index 2", "550.0"]),
    ("descending_wavelengths", dict(wavelengths=_WLS[::-1].copy()),
     ["wavelengths", "strictly increasing", "index 1"]),
    ("nan_angle", dict(angles=[30.0, np.nan]), ["angles", "index 1", "nan"]),
    ("angle_above_90", dict(angles=[95.0]), ["angles", "index 0", "95.0"]),
    ("angle_retrograde", dict(angles=[120.0]), ["angles", "index 0", "120.0"]),
    ("angle_negative", dict(angles=[-30.0]), ["angles", "index 0", "-30.0"]),
    ("angle_radians_out_of_range", dict(angles=[2.0], angles_in_radians=True),
     ["angles", "index 0", "2.0", "rad"]),
]


@pytest.mark.parametrize("kwargs,fragments", [r[1:] for r in _REJECT],
                         ids=[r[0] for r in _REJECT])
def test_bad_input_raises_and_names_the_offender(kwargs, fragments):
    """Every rejected input must raise ``ValueError`` naming index and value.

    The index matters as much as the rejection: "a thickness is negative" sends
    the caller hunting through a 40-layer stack, "layer 2: -50.0" does not.
    """
    with pytest.raises(ValueError) as excinfo:
        _sm(**kwargs)
    message = str(excinfo.value)
    missing = [f for f in fragments if f not in message]
    assert not missing, f"message does not name {missing}: {message!r}"


# --------------------------------------------------------------------------
# Accepted inputs -- the checks must not overreach
# --------------------------------------------------------------------------

_ACCEPT = [
    ("grazing_90_deg", dict(angles=[90.0])),
    ("normal_0_deg", dict(angles=[0.0])),
    ("radians_pi_over_2", dict(angles=[np.pi / 2.0], angles_in_radians=True)),
    ("single_wavelength", dict(wavelengths=np.array([550.0]))),
    ("zero_thickness_everywhere", dict(thicknesses=np.zeros(4))),
    ("huge_thickness_1e9", dict(thicknesses=_with(_D, 2, 1e9))),
    ("strongly_absorbing_index", dict(layer_indices=_with(_N, 1, 1.2 + 7.5j))),
    ("index_below_one", dict(layer_indices=_with(_N, 1, 0.2 + 3.0j))),
    ("two_dim_index_array", dict(layer_indices=np.repeat(_N[:, None], _WLS.size, axis=1))),
    # R3.4 singles out the *incident* side. Absorption anywhere else is
    # ordinary physics and must stay untouched.
    ("absorbing_substrate", dict(layer_indices=_with(_N, 3, 1.52 + 0.05j))),
    ("strongly_absorbing_substrate", dict(layer_indices=_with(_N, 3, 0.9 + 6.5j))),
    # 0.6.26: an absorbing ambient is corrected and warned about, not refused,
    # so it belongs here. That it still *solves* is the point -- the numbers it
    # produces are pinned below.
    ("absorbing_incident_medium", dict(layer_indices=_with(_N, 0, 1.52 + 0.05j))),
    ("faintly_absorbing_incident_medium",
     dict(layer_indices=_with(_N, 0, 1.0 + 1e-14j))),
]


@pytest.mark.parametrize("kwargs", [a[1] for a in _ACCEPT], ids=[a[0] for a in _ACCEPT])
def test_legitimate_input_still_constructs_and_solves(kwargs):
    """Odd but physical inputs must survive.

    Grazing incidence, a one-point grid, a metallic index with n < 1, a
    millimetre-thick layer: all legal, all previously working, and a
    validation layer that rejects any of them has broken the library to fix a
    bug. The solve is included because construction alone would not catch a
    check that mangles the arrays on the way through.
    """
    out = _sm(**kwargs).compute(Request.RS | Request.TS)
    assert set(out) == {"Rs", "Ts"}
    assert np.all(np.isfinite(np.asarray(out["Rs"], float)))


# --------------------------------------------------------------------------
# Positive controls: rejecting is not the same as repairing
# --------------------------------------------------------------------------

def test_ascending_grid_is_stored_verbatim():
    """Reject, never reorder.

    Silently sorting a descending grid would desynchronize it from the
    caller's own wavelength-indexed arrays -- their n(lambda) table would no
    longer line up with the grid the engine used, and nothing would say so.
    This asserts the grid comes out exactly as it went in.
    """
    st = _sm()
    assert np.array_equal(st.wavls, _WLS)
    assert st.wavls[0] < st.wavls[-1]


def test_a_rejected_construction_does_not_mutate_the_caller_arrays():
    """The validator reads; it must not write.

    A check that repaired its input in place would leave the caller holding a
    silently different array after catching the error.
    """
    wls = _WLS[::-1].copy()
    thick = _with(_D, 2, -50.0)
    idx = _with(_N, 1, np.nan)
    wls_before, thick_before = wls.copy(), thick.copy()
    idx_before = idx.copy()

    for kwargs in (dict(wavelengths=wls), dict(thicknesses=thick),
                   dict(layer_indices=idx)):
        with pytest.raises(ValueError):
            _sm(**kwargs)

    assert np.array_equal(wls, wls_before)
    assert np.array_equal(thick, thick_before)
    assert np.array_equal(idx[~np.isnan(idx)], idx_before[~np.isnan(idx_before)])
    assert np.isnan(idx[1])


def test_the_first_offender_is_the_one_reported():
    """With several bad entries, the message names the first, not an arbitrary one."""
    bad = _with(_with(_D, 1, -1.0), 3, -2.0)
    with pytest.raises(ValueError, match=r"layer 1: -1\.0"):
        _sm(thicknesses=bad)


def test_validation_runs_before_the_native_solver_is_built():
    """Bad input must not reach the engine at all.

    If the native Solver were constructed first, a rejected stack would still
    have paid for (and possibly warned from) the native setup, and a future
    native panic on garbage input would surface instead of the clear error.
    """
    import navette._smatrix as native

    calls = []
    real = native.Solver

    class _Spy:
        def __new__(cls, *a, **kw):
            calls.append(1)
            return real(*a, **kw)

    import navette.smatrix.smatrix as mod
    original = mod._NativeSolver
    mod._NativeSolver = _Spy
    try:
        with pytest.raises(ValueError):
            _sm(thicknesses=_with(_D, 2, -50.0))
        assert calls == [], "native Solver was constructed for a rejected stack"
        _sm()
        assert calls == [1], "native Solver was not constructed for a valid stack"
    finally:
        mod._NativeSolver = original


# --------------------------------------------------------------------------
# Needle depth and index range (R3.2)
# --------------------------------------------------------------------------
# `needle_slopes4_ddz` guarded its host-layer invariant with `debug_assert!`,
# which is compiled out of every release build -- the one everybody actually
# runs. `locate_depth_in` does not reject an out-of-range z either: it falls
# through to the last layer of the range and returns a depth past that layer's
# thickness, so the kernel computed a gradient for a needle that is not where
# the caller put it, and returned it without complaint.
#
# The contract now lives once in `solver::needle_gradient`, which both PyO3
# entry points funnel through. The kernel keeps its `debug_assert!`s as Rust-
# caller invariants: a release `assert!` in an O(1) hot kernel would trade a
# silent wrong answer for a panic, which is worse on a library path.

_NEEDLE_N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
_NEEDLE_D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])   # host span = 400 nm
_NEEDLE_SPAN = 400.0


def _needle_stack(**kw):
    return ScatterMatrix(_NEEDLE_N, _NEEDLE_D, wavelengths=_WLS, angles=[30.0], **kw)


def _needle(z, stack=None, n_prime=None, **kw):
    from navette.smatrix.needle import NeedleRequest, needle_gradient

    st = stack if stack is not None else _needle_stack()
    nn = n_prime if n_prime is not None else np.full(_WLS.size, 1.8 + 0j,
                                                     dtype=np.complex128)
    kw.setdefault("pol", "s")
    return needle_gradient(st, nn, np.atleast_1d(z), NeedleRequest.P, **kw)


@pytest.mark.parametrize("z,fragments", [
    (1000.0, ["z_grid[0]", "1000", "[0, 400]", "coherent block"]),
    (-10.0, ["z_grid[0]", "-10", "[0, 400]"]),
    (_NEEDLE_SPAN + 1e-3, ["z_grid[0]", "[0, 400]"]),
    (np.nan, ["z_grid[0]", "not finite"]),
    (np.inf, ["z_grid[0]", "not finite"]),
], ids=["past_the_stack", "negative", "just_past_the_span", "nan", "inf"])
def test_needle_depth_outside_the_span_raises(z, fragments):
    with pytest.raises(ValueError) as excinfo:
        _needle(z)
    message = str(excinfo.value)
    missing = [f for f in fragments if f not in message]
    assert not missing, f"message does not name {missing}: {message!r}"


@pytest.mark.parametrize("z", [0.0, 120.0, _NEEDLE_SPAN, _NEEDLE_SPAN + 1e-12,
                               -1e-12],
                         ids=["top", "interior", "bottom", "within_tol_high",
                              "within_tol_low"])
def test_legal_needle_depth_is_accepted(z):
    """Both endpoints are legal, and the 1e-9 tolerance absorbs round-off.

    A caller that builds its z grid as ``np.linspace(0, sum(d), k)`` lands on
    the bottom endpoint with whatever error the sum accumulated. Rejecting
    that would make the obvious way to write the call fail intermittently.
    """
    out = _needle(z)
    assert np.all(np.isfinite(np.asarray(out["P_s"], float)))


def test_the_span_follows_the_block_not_the_stack():
    """``end_idx`` narrows the eligible span, and the error says by how much.

    The coherent kernels are confined to ``[start_idx, end_idx]``; a z that is
    inside the stack but outside the block has no host, and used to be clamped
    into the block's last layer.
    """
    _needle(50.0, end_idx=2)                       # inside layer 1: fine
    with pytest.raises(ValueError, match=r"\[0, 120\].*coherent block"):
        _needle(200.0, end_idx=2)


def test_the_multiblock_span_is_not_the_coherent_one():
    """A multiblock request walks the whole stack, so it keeps the wider span.

    Checking one span for both paths would reject legitimate multiblock depths
    (or wave through illegitimate coherent ones). Each is checked only when the
    request actually reaches that path.
    """
    from navette.smatrix.needle import NeedleRequest, needle_gradient

    flags = np.array([0, 0, 1, 0, 0], dtype=np.int32)
    st = _needle_stack(incoherent_flags=flags)
    nn = np.full(_WLS.size, 1.8 + 0j, dtype=np.complex128)

    # 350 nm sits in layer 3 -- past the default coherent block, inside the
    # cascade. The multiblock request must accept it.
    needle_gradient(st, nn, [350.0], NeedleRequest.P_MB, pol="s")
    with pytest.raises(ValueError, match="multiblock cascade"):
        needle_gradient(st, nn, [500.0], NeedleRequest.P_MB, pol="s")


def test_non_finite_needle_index_raises():
    """A NaN needle index made the entire gradient NaN with nothing to point at."""
    nn = np.full(_WLS.size, 1.8 + 0j, dtype=np.complex128)
    nn[3] = np.nan
    with pytest.raises(ValueError, match=r"needle_n_per_wav\[3\]"):
        _needle(200.0, n_prime=nn)


# --------------------------------------------------------------------------- #
# The raw core_engine binding (R5.3)
# --------------------------------------------------------------------------- #
# The binding hands its index cache to the solver in the wav-major layout it
# arrives in, instead of transposing it into a layer-major Vec first. The old
# transpose loop indexed that buffer directly, so a cache of the wrong length
# was an out-of-bounds panic -- and a panic here aborts the interpreter rather
# than raising. The length is now a checked precondition of the constructor.


def _core_engine_args(n_wavs=8, n_layers=3):
    """The raw binding's argument tuple, in its own flat layout."""
    import navette._smatrix as native
    from navette.smatrix.smatrix import CoherenceMode

    wavls = np.linspace(450.0, 750.0, n_wavs)
    n_flat = np.tile(np.array([1.0, 0.0, 2.35, 0.01, 1.52, 0.0]), n_wavs)
    return native.core_engine, (
        wavls, np.array([0.0]), n_layers, n_flat,
        np.array([0.0, 120.0, 0.0]), np.zeros(n_layers, dtype=np.int32),
        np.zeros(n_layers, dtype=np.int32), np.zeros(n_layers),
        int(CoherenceMode.FRONT_BLOCK), int(Request.RS))


def test_a_short_index_cache_raises_instead_of_panicking():
    engine, args = _core_engine_args()
    assert engine(*args)["Rs"].shape == (1, 8)
    short = list(args)
    short[3] = args[3][:-2]
    with pytest.raises(ValueError, match="index cache length"):
        engine(*short)


def test_a_long_index_cache_is_refused_too():
    """Too much data is as wrong as too little -- silently ignoring the tail
    would solve a different stack than the caller described."""
    engine, args = _core_engine_args()
    long = list(args)
    long[3] = np.concatenate([args[3], np.zeros(2)])
    with pytest.raises(ValueError, match="index cache length"):
        engine(*long)


# ---------------------------------------------------------------------------
# The Sellmeier domain (R6.2)
# ---------------------------------------------------------------------------

# Schott BK7: the first resonance is at sqrt(C1) um = 77.46 nm.
_BK7 = dict(B1=1.03961212, C1=0.00600069867, B2=0.231792344,
            C2=0.0200179144, B3=1.01046945, C3=103.560653)


def test_a_sellmeier_grid_inside_the_fit_evaluates():
    from navette.materials import MaterialSpec, evaluate
    n = evaluate(MaterialSpec("Sellmeier", dict(_BK7)), np.linspace(400.0, 800.0, 9))
    assert np.all(np.isfinite(n.real))
    assert abs(n[np.searchsorted(np.linspace(400.0, 800.0, 9), 550.0)].real - 1.5185) < 5e-3


def test_a_sellmeier_grid_past_a_resonance_raises_instead_of_returning_nan():
    """A NaN index is the worst kind of wrong answer: it survives every layer
    and arrives at the user as a NaN spectrum with nothing pointing back at
    the material that produced it."""
    from navette.materials import MaterialSpec, evaluate
    spec = MaterialSpec("Sellmeier", dict(_BK7))
    with pytest.raises(ValueError, match="Sellmeier"):
        evaluate(spec, np.array([550.0, 70.0]))
    try:
        evaluate(spec, np.array([70.0]))
    except ValueError as e:
        msg = str(e)
    assert "77.4" in msg, msg
    assert msg.isascii(), "the message must survive a cp1252 console"


def test_the_urbach_variant_is_guarded_too():
    from navette.materials import MaterialSpec, evaluate
    p = dict(_BK7, alpha0=1e5, Eu=0.06, lambda_g=380.0)
    with pytest.raises(ValueError, match="Sellmeier"):
        evaluate(MaterialSpec("SellmeierUrbach", p), np.array([70.0]))


def test_dop_r_never_exceeds_one():
    """Reflected DOP was unclamped while transmitted was not (review 3.3);
    round-off alone was measured at 1.0000000000000004."""
    rng = np.random.default_rng(20260911)
    worst = 0.0
    for _ in range(40):
        nl = int(rng.integers(3, 7))
        n = np.array([1.0 + 0j]
                     + [complex(rng.uniform(1.3, 3.5), rng.uniform(0.0, 0.4))
                        for _ in range(nl - 2)]
                     + [complex(1.515, 0.001)])
        d = np.array([0.0] + [rng.uniform(5.0, 400.0) for _ in range(nl - 2)] + [0.0])
        sm = ScatterMatrix(n, d, wavelengths=np.linspace(300.0, 1200.0, 301),
                           angles=[float(rng.uniform(0.0, 85.0))],
                           coherence_mode=CoherenceMode.COHERENCY_MATRIX)
        out = sm.compute(Request.DOP_R | Request.DOP_T)
        worst = max(worst, float(np.max(np.asarray(out["DOP_R"]))))
        assert float(np.max(np.asarray(out["DOP_T"]))) <= 1.0
    assert worst <= 1.0
    # And the clamp is not hiding a channel that never gets near the bound.
    assert worst > 0.99


# --------------------------------------------------------------------------
# R3.4 -- the incident medium (refused in 0.6.21; corrected + warned since 0.6.26)
# --------------------------------------------------------------------------

def test_absorbing_incident_medium_warns_and_explains_itself():
    """The correction has to announce itself and carry its justification.

    A silent fix is the failure mode this whole item exists to end: before
    0.6.21 the engine quietly returned R = 417 at 10 degrees. The warning has
    to say three things -- that layer 0's absorption was dropped, why the
    absorbing-ambient problem has no reflectance to return, and what the
    caller is getting instead -- or it invites the reader to assume the number
    means something it does not.
    """
    with pytest.warns(UserWarning) as record:
        _sm(layer_indices=_with(_N, 0, 1.52 + 0.05j))
    msgs = [str(w.message) for w in record]
    hits = [m for m in msgs if "incident medium" in m]
    assert len(hits) == 1, f"expected exactly one incident-medium warning, got {msgs!r}"
    msg = hits[0]
    assert msg.isascii(), "the message must survive a cp1252 console"
    for fragment in ("incident medium", "dropped", "transparent ambient",
                     "energy ratio", "substrate", "0.05"):
        assert fragment in msg, f"{fragment!r} missing from {msg!r}"


def test_absorbing_ambient_is_solved_as_its_transparent_twin():
    """The correction is exact, not approximate: same stack, Im(n[0]) = 0.

    Bit-for-bit rather than ``allclose`` -- "we dropped k" has one correct
    implementation and any drift from it means something else was touched. The
    energy check is the payoff: R + T was 1.0096 at k = 0.1 and 1.44 at k = 1
    before this item, and 417 at 10 degrees.
    """
    req = Request.RS | Request.TS | Request.RP | Request.TP
    with pytest.warns(UserWarning):
        absorbing = ScatterMatrix(_with(_N, 0, 1.0 + 0.3j), _D,
                                  wavelengths=_WLS, angles=_ANG).compute(req, squeeze=False)
    transparent = ScatterMatrix(_with(_N, 0, 1.0 + 0.0j), _D,
                                wavelengths=_WLS, angles=_ANG).compute(req, squeeze=False)
    for key in transparent:
        assert np.array_equal(absorbing[key], transparent[key]), f"{key} drifted"

    # Lossless stack: the residual is now round-off, not 44 percent.
    for r, t in (("Rs", "Ts"), ("Rp", "Tp")):
        resid = np.abs(np.asarray(absorbing[r], float)
                       + np.asarray(absorbing[t], float) - 1.0)
        assert float(np.max(resid)) < 1e-12, f"{r}+{t} off by {float(np.max(resid))}"


def test_only_the_incident_row_is_touched_in_a_2d_index_array():
    """A per-wavelength index grid is judged and corrected on row 0 only.

    A dispersive absorbing *layer* three rows down shares the array and must
    survive untouched -- and so must the caller's array itself, which reaches
    the constructor without being copied when it is already complex128.
    """
    grid = np.repeat(_N[:, None], _WLS.size, axis=1)
    grid[1, :] = 2.35 + 0.4j          # absorbing interior layer: fine
    grid[3, :] = 1.52 + 0.02j         # absorbing substrate: fine
    out = ScatterMatrix(grid, _D, wavelengths=_WLS, angles=_ANG).compute(Request.RS)
    assert np.all(np.isfinite(np.asarray(out["Rs"], float)))

    # One absorbing wavelength in the ambient row is enough to warn, and the
    # warning counts them rather than claiming the whole row.
    grid[0, 7] = 1.0 + 1e-9j
    before = grid.copy()
    with pytest.warns(UserWarning, match=r"absorbing at 1 of "):
        sm = ScatterMatrix(grid, _D, wavelengths=_WLS, angles=_ANG)
    assert np.array_equal(grid, before), "the caller's array was mutated"

    # The interior and substrate rows still absorb; only layer 0 was flattened.
    got = sm.compute(Request.RS | Request.A_S, squeeze=False)
    assert float(np.max(np.asarray(got["A_s"], float))) > 1e-3


def test_the_native_solver_is_still_permissive():
    """The escape hatch R3.1 promised has to actually be open.

    The wrapper corrects; the engine does not. A caller who knows what an
    absorbing ambient means for their amplitudes can still reach the raw ones,
    which is what makes the correction a convenience rather than a lost
    capability.
    """
    from navette._smatrix import Solver, solver_rt_request

    n = _with(_N, 0, 1.52 + 0.05j)
    idx2d = np.repeat(np.asarray(n, dtype=np.complex128)[:, None], _WLS.size, axis=1)
    solver = Solver(_WLS, np.asarray(_ANG, dtype=float),
                    np.ascontiguousarray(idx2d).ravel(), _N.size,
                    thicknesses=_D)
    out = solver.solve(solver_rt_request(pol="s"))
    assert out  # it answers; interpreting the answer is the caller's problem
