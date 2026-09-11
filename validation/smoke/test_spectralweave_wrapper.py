# -*- coding: utf-8 -*-
"""``SimulationWeaver``'s edges, as the broken example found them (R4.3).

Three defects, all reachable from a first-time reader's transcript:

  * ``unweave`` reported ``Length mismatch`` from deep inside the native
    distribution plan when the supplied target grid did not *contain* each
    frame's wavelengths. Distribution is by exact match (1e-12), not by
    interpolation, and nothing said so;
  * a key tuple handed in where an ``OpticalFragment`` belongs raised
    ``AttributeError: 'tuple' object has no attribute '_rust_key'``; and
  * ``OpticalFragment`` was a frozen dataclass holding two numpy arrays, so
    the synthesised ``__hash__`` raised -- which made ``unweave_batch``, whose
    signature is ``dict[OpticalFragment, np.ndarray]``, impossible to call.
"""

import numpy as np
import pytest

from navette.spectralweave import OpticalFragment, SimulationWeaver

_WL_A = np.linspace(400.0, 499.0, 100)
_WL_B = np.linspace(500.0, 600.0, 201)
_KEY = dict(base_wavelength=550.0, data_type="R", polarization="s")


@pytest.fixture
def woven():
    w = SimulationWeaver(cache_size=16)
    for wl in (_WL_A, _WL_B):
        w.add_fragment(OpticalFragment(wavelengths=wl,
                                       values=np.sin(wl / 40.0), **_KEY))
    return w


def _template(wl, values):
    return OpticalFragment(wavelengths=wl, values=values, **_KEY)


# --------------------------------------------------------------------------
# 1. The round trip the example demonstrates
# --------------------------------------------------------------------------

def test_weave_then_unweave_is_a_round_trip(woven):
    wl, values = woven.get_continuous_curve(**_KEY)
    assert wl.size == _WL_A.size + _WL_B.size

    target = values * 1.5
    assert woven.unweave(_template(wl, target), wl, target) == 2

    wl_after, after = woven.get_continuous_curve(**_KEY)
    assert np.array_equal(wl_after, wl)
    assert np.allclose(after, target)


# --------------------------------------------------------------------------
# 2. The opaque native error
# --------------------------------------------------------------------------

def test_a_grid_that_only_spans_the_range_is_rejected_by_name(woven):
    """The example's own failure, with the message it should have had.

    ``linspace(400, 600, 500)`` covers the same interval as the two frames
    and shares almost no wavelength with either. The message must say how
    many points were covered, which frame, and that matching is exact.
    """
    stranger = np.linspace(400.0, 600.0, 500)
    with pytest.raises(ValueError) as excinfo:
        woven.unweave(_template(stranger, np.cos(stranger / 100.0)),
                      stranger, np.cos(stranger / 100.0))
    msg = str(excinfo.value)
    for fragment in ("unweave", "exact wavelength match", "get_weaved", "100"):
        assert fragment in msg, f"message omits {fragment!r}: {msg}"


def test_a_grid_covering_one_frame_exactly_still_works(woven):
    """Partial coverage is fine as long as each *touched* frame is complete.

    Guards the check against over-rejecting: a caller updating only the blue
    half must not be told their grid is wrong.
    """
    vals = np.full(_WL_A.size, 0.25)
    assert woven.unweave(_template(_WL_A, vals), _WL_A, vals) == 1
    _, after = woven.get_continuous_curve(**_KEY)
    assert np.allclose(after[:_WL_A.size], 0.25)
    # The untouched frame kept its own values.
    assert not np.allclose(after[_WL_A.size:], 0.25)


# --------------------------------------------------------------------------
# 3. The wrapper guards
# --------------------------------------------------------------------------

def test_a_key_tuple_where_a_fragment_belongs_is_a_TypeError(woven):
    with pytest.raises(TypeError, match="needs an OpticalFragment"):
        woven.unweave((550.0, "R", "s"), _WL_A, np.zeros(_WL_A.size))


def test_the_type_error_names_what_it_got(woven):
    with pytest.raises(TypeError, match="got NoneType"):
        woven.unweave(None, _WL_A, np.zeros(_WL_A.size))


def test_mismatched_curve_lengths_are_caught_before_the_backend(woven):
    with pytest.raises(ValueError, match="one value per wavelength"):
        woven.unweave(_template(_WL_A, np.zeros(_WL_A.size)),
                      _WL_A, np.zeros(_WL_A.size + 1))


# --------------------------------------------------------------------------
# 4. OpticalFragment as a dict key
# --------------------------------------------------------------------------

def test_a_fragment_can_be_a_dict_key(woven):
    """``unweave_batch``'s whole signature depends on this."""
    frag = _template(_WL_A, np.zeros(_WL_A.size))
    assert {frag: 1}[frag] == 1


def test_fragments_are_distinguished_by_identity(woven):
    """Equal arrays, two objects, two keys -- and neither comparison raises."""
    a = _template(_WL_A, np.zeros(_WL_A.size))
    b = _template(_WL_A.copy(), np.zeros(_WL_A.size))
    assert a != b and a == a
    assert len({a: 0, b: 0}) == 2


def test_unweave_batch_writes_every_key(woven):
    wl, values = woven.get_continuous_curve(**_KEY)
    batch = {
        OpticalFragment(wavelengths=wl, values=values, base_wavelength=550.0,
                        data_type="R", polarization=pol): values * scale
        for pol, scale in (("s", 1.0), ("p", 0.5))
    }
    # Two keys x two frames.
    assert woven.unweave_batch(wl, batch) == 4
    _, got_p = woven.get_continuous_curve(550.0, "R", "p")
    assert np.allclose(got_p, values * 0.5)


def test_a_fragment_still_rejects_mismatched_arrays():
    with pytest.raises(ValueError, match="same shape"):
        OpticalFragment(wavelengths=_WL_A, values=np.zeros(3), **_KEY)


# --------------------------------------------------------------------------
# 5. The batch path's own guards (R5.1)
# --------------------------------------------------------------------------

def test_a_short_curve_in_a_batch_is_an_error_not_a_crash(woven):
    """``unweave`` checks this in the wrapper; ``unweave_batch`` never did.

    The distribution plan addresses the source curve by position, so a short one
    indexed straight out of the buffer -- a Rust panic surfacing as
    ``pyo3_runtime.PanicException``, which does not even derive from
    ``Exception``. It is an error message now.
    """
    wl, values = woven.get_continuous_curve(**_KEY)
    frag = OpticalFragment(wavelengths=wl, values=values, **_KEY)
    with pytest.raises(ValueError, match="one value per wavelength"):
        woven.unweave_batch(wl, {frag: values[:-1]})


def test_the_batch_length_error_names_the_offending_key(woven):
    wl, values = woven.get_continuous_curve(**_KEY)
    frag = OpticalFragment(wavelengths=wl, values=values,
                           base_wavelength=550.0, data_type="R",
                           polarization="p")
    with pytest.raises(ValueError, match=r"\(550, R, p\)"):
        woven.unweave_batch(wl, {frag: values[:-1]})


def test_a_rejected_batch_writes_nothing(woven):
    """All-or-nothing: one bad curve must not leave the others half-applied."""
    wl, values = woven.get_continuous_curve(**_KEY)
    good = OpticalFragment(wavelengths=wl, values=values, base_wavelength=550.0,
                           data_type="R", polarization="good")
    bad = OpticalFragment(wavelengths=wl, values=values, base_wavelength=550.0,
                          data_type="R", polarization="bad")
    with pytest.raises(ValueError):
        woven.unweave_batch(wl, {good: values * 2.0, bad: values[:-1]})
    with pytest.raises(Exception):
        woven.get_continuous_curve(550.0, "R", "good")


def test_a_batch_writes_what_the_single_calls_write(woven):
    """The batch path shares one plan and one materialised source across keys;
    ``unweave`` does neither. They must still agree (R5.1)."""
    wl, values = woven.get_continuous_curve(**_KEY)
    scales = {"s": 1.0, "p": 0.5, "u": -3.25}
    frags = {
        pol: OpticalFragment(wavelengths=wl, values=values,
                             base_wavelength=550.0, data_type="R",
                             polarization=pol)
        for pol in scales
    }
    woven.unweave_batch(wl, {f: values * scales[pol] for pol, f in frags.items()})

    single = SimulationWeaver(cache_size=16)
    for sub in (_WL_A, _WL_B):
        single.add_fragment(OpticalFragment(wavelengths=sub,
                                            values=np.sin(sub / 40.0), **_KEY))
    for pol, f in frags.items():
        single.unweave(f, wl, values * scales[pol])

    for pol in scales:
        _, from_batch = woven.get_continuous_curve(550.0, "R", pol)
        _, from_single = single.get_continuous_curve(550.0, "R", pol)
        assert np.array_equal(from_batch, from_single)
        assert np.array_equal(from_batch, values * scales[pol])
