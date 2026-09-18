# SPDX-License-Identifier: LGPL-3.0-or-later
"""0.7.12 (C3, C5): a flag that cannot mean what it says, and one that does nothing.

Two warnings, neither an error, because in both cases the engine's behaviour
is well defined and the caller may have meant it:

* **C3 — thin flagged layer.** The flag asserts the layer destroys the phase
  relation between its two surfaces, which needs the path-length spread
  across it to exceed the source coherence length ``L_C = lambda^2 /
  delta_lambda``. Nothing checked that. The threshold here is *derived*
  rather than picked: the warning quotes the source bandwidth that layer
  would need, so "too thin" is a number the caller can argue with. Round 1
  measured why it matters — a flagged 200 nm air gap past the critical angle
  returns ``R = 1.000000`` exactly, the thick limit unconditionally, where
  the coherent stack gives ``R = 0.763`` — and that the split alone moves the
  answer even at zero thickness (``Rs 0.193432 -> 0.109790``).
* **C5 — half-space flag.** Rows 0 and last are half-spaces. The block sweep
  scans ``current_idx + 1`` up to ``idx_n`` and applies the attenuation
  element only below ``idx_n``, so their flags are never read: bit-identical
  output, no diagnostic. The constructor docstring made this worse by naming
  a "thick substrate" as the example, which is exactly the row a caller would
  flag to no effect.

Neither refuses. A half-space has no second surface to lose coherence
against, so a flag there is a no-op rather than a mistake; and the thin
incoherent limit is a legitimate model, as long as it is the one the caller
meant. C2's refusal deliberately uses the same interior-only predicate, so
these three findings cannot contradict each other.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from navette.smatrix.smatrix import Request, ScatterMatrix

WLS = np.array([500.0, 600.0])
_THIN = dict(
    layer_indices=np.array([1.0, 1.45, 2.35, 1.52]),
    thicknesses=[0.0, 60.0, 100.0, 0.0],
    wavelengths=WLS,
    angles=[0.0],
)
_THICK = dict(
    layer_indices=np.array([1.0, 1.45, 1.52]),
    thicknesses=[0.0, 50_000.0, 0.0],
    wavelengths=WLS,
    angles=[0.0],
)


def _warns(**kwargs):
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        ScatterMatrix(**kwargs)
    return [str(w.message) for w in rec]


# ---- C3 --------------------------------------------------------------------
def test_a_thin_flagged_layer_warns_and_still_builds():
    msgs = _warns(**_THIN, incoherent_flags=[0, 1, 0, 0])
    assert len(msgs) == 1, msgs
    msg = msgs[0]
    assert "row 1" in msg
    assert "87.0000" in msg, "the optical thickness n*d, not the thickness"
    assert "0.174 wavelengths" in msg
    assert msg.isascii(), "the message must survive a cp1252 console"


def test_the_threshold_is_derived_not_chosen():
    """The warning quotes the bandwidth the layer would need to be incoherent.

    ``delta_lambda > lambda**2 / (2*n*d)``. That is the honest form of "too
    thin": a number the caller can check against their own source instead of
    a constant somebody picked.
    """
    msg = _warns(**_THIN, incoherent_flags=[0, 1, 0, 0])[0]
    nd = 1.45 * 60.0
    lam = 500.0
    need = lam * lam / (2.0 * nd)
    assert f"{need:.4f}" in msg, f"expected delta_lambda > {need:.4f} in: {msg}"
    assert f"{100.0 * need / lam:.0f}%" in msg
    assert "L_C = lambda^2 / delta_lambda" in msg


def test_a_genuinely_thick_flagged_layer_says_nothing():
    """50 um is the substrate regime the physics actually supports."""
    assert _warns(**_THICK, incoherent_flags=[0, 1, 0]) == []


def test_the_flag_is_still_honoured_after_the_warning():
    """A warning, not a refusal: the numbers must change, or it is a lie.

    This is also the C3 measurement in miniature — the flagged answer is the
    incoherent one, which is the whole reason a thin flag is dangerous rather
    than merely pointless.
    """
    req = Request.RS | Request.TS
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flagged = ScatterMatrix(**_THIN, incoherent_flags=[0, 1, 0, 0]).compute(
            req, squeeze=False)
    plain = ScatterMatrix(**_THIN).compute(req, squeeze=False)
    assert not np.array_equal(
        np.asarray(flagged["Rs"]), np.asarray(plain["Rs"])
    ), "the flag must still partition the stack"


def test_a_zero_thickness_flagged_layer_does_not_warn_about_thickness():
    """It has no optical thickness to be wrong about.

    It still decoheres — the join breaks the p-s phase relation whatever the
    thickness, which is round 1's measurement — but a bandwidth figure for a
    zero-thickness layer would be infinite and meaningless, so the warning
    stays quiet rather than printing one.
    """
    msgs = _warns(
        layer_indices=np.array([1.0, 1.45, 2.35, 1.52]),
        thicknesses=[0.0, 0.0, 100.0, 0.0],
        wavelengths=WLS, angles=[0.0], incoherent_flags=[0, 1, 0, 0],
    )
    assert [m for m in msgs if "optical thickness" in m] == []


# ---- C5 --------------------------------------------------------------------
@pytest.mark.parametrize("flags", [[1, 0, 0], [0, 0, 1], [1, 0, 1]])
def test_a_half_space_flag_warns_that_it_does_nothing(flags):
    msgs = [m for m in _warns(**_THICK, incoherent_flags=flags)
            if "no effect" in m]
    assert len(msgs) == 1, msgs
    assert "half-spaces" in msgs[0]
    assert "interior layer" in msgs[0], "must say what to do instead"
    assert msgs[0].isascii()


@pytest.mark.parametrize("flags", [[1, 0, 0], [0, 0, 1], [1, 0, 1]])
def test_a_half_space_flag_really_is_bit_identical(flags):
    """The claim in the warning, measured rather than asserted.

    If this ever stops holding, the warning is wrong and C2's interior-only
    predicate is too.
    """
    req = Request.RS | Request.TS
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        got = ScatterMatrix(**_THICK, incoherent_flags=flags).compute(
            req, squeeze=False)
    plain = ScatterMatrix(**_THICK).compute(req, squeeze=False)
    for key in plain:
        assert np.array_equal(np.asarray(got[key]), np.asarray(plain[key])), key


def test_an_interior_flag_is_not_a_half_space_flag():
    """The control that keeps the C5 test above from being vacuous."""
    req = Request.RS | Request.TS
    flagged = ScatterMatrix(**_THICK, incoherent_flags=[0, 1, 0]).compute(
        req, squeeze=False)
    plain = ScatterMatrix(**_THICK).compute(req, squeeze=False)
    assert not np.array_equal(
        np.asarray(flagged["Rs"]), np.asarray(plain["Rs"])
    )


# ---- one rule, two doors ---------------------------------------------------
def test_both_doors_explain_thin_flags_the_same_way():
    """``ScatterMatrix`` (Python) and ``solve_arrays`` (Rust) carry one text.

    Same contract as the ambient, cross-channel and gain rules: separate
    implementations so the Python one can point ``stacklevel`` at the caller,
    one shared explanation so they cannot drift.
    """
    from navette.structure import Layer, Navette_Structure, solve_structure
    from navette.structure.materials import DictMaterialProvider

    py_msg = _warns(**_THIN, incoherent_flags=[0, 1, 0, 0])[0]
    shared = py_msg[py_msg.index("An incoherent flag says"):]
    assert len(shared) > 600, f"parity slice looks truncated: {shared!r}"

    mats = {"glass": np.full(WLS.size, 1.52 + 0j),
            "thin": np.full(WLS.size, 1.45 + 0j)}
    st = Navette_Structure(
        [Layer(0.0, "glass"), Layer(60.0, "thin", coherent=False),
         Layer(0.0, "glass")],
        {}, DictMaterialProvider(dict(mats), wavelength=WLS))
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        solve_structure(st, WLS, 0.0, request=int(Request.RS))
    rust_msgs = [str(w.message) for w in rec if "incoherent layer at row" in str(w.message)]
    assert rust_msgs, f"the Rust door said nothing: {[str(w.message) for w in rec]}"
    assert shared in rust_msgs[0], (
        f"the two doors have drifted apart / python: {shared!r} / "
        f"rust: {rust_msgs[0]!r}"
    )
    assert rust_msgs[0].isascii()


def test_the_two_findings_do_not_contradict_each_other():
    """C3, C5 and C2 all use the same interior-only rule.

    A stack flagged only on its half-spaces gets C5's no-op warning, never
    C3's thickness warning and never C2's refusal — those two are about a
    layer the sweep actually partitions on, and this is not one.
    """
    msgs = _warns(**_THICK, incoherent_flags=[1, 0, 1])
    assert all("optical thickness" not in m for m in msgs), msgs
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # C2 would refuse a cross observable here if the predicate included
        # half-space rows. It must not.
        ScatterMatrix(**_THICK, incoherent_flags=[1, 0, 1]).compute(
            Request.DOP_R, squeeze=False)
