# SPDX-License-Identifier: LGPL-3.0-or-later
"""0.7.11 (C4): ``Im(n) < 0`` is refused at every door, with one explanation.

Optical gain has no representation in either solver path, and the two paths
mangle it *differently*:

* the coherent kernels conjugate the propagation phase back to decay
  (``coherent_block.rs:146``/``:345``, mirrored at
  ``needle_operator.rs:248``), so a gain layer returns the LOSS layer's
  answer -- on a symmetric stack bit-for-bit identical to ``+k`` -- and
  reports a POSITIVE absorptance for a medium that amplifies;
* the incoherent cascade clamps it to zero instead (``core_engine.rs:503``/
  ``:665``, ``spacer_tau`` at ``needle_operator.rs:1290``), so the layer goes
  transparent and the energy books open by ``A = -4.3e-5``.

Neither is any physical system, and neither is recoverable from the output:
it reads as an ordinary absorbing stack. Round 1 measured both (physics
review C4); round 2 sharpened the first into "bit-identical for gain and
loss", i.e. the sign is erased outright rather than merely mis-booked.

WHERE THE REFUSAL SITS, and why that differs from C2. C2's refusal had to
stay at the doors because Mode A's numbers are a deliberate legacy port whose
parity suite enters through the raw ``core_engine`` pyfunction, below them.
Gain has no such constituency -- nothing below the doors computes anything
worth preserving -- so the check goes at ``Solver::assemble``, the one place
every constructor funnels through, and covers ``ScatterMatrix``,
``solve_arrays``/``solve_structure`` and the raw FFI at once. The free
``needle_gradient`` takes a flat cache and never builds a ``Solver``, so it
carries its own copy, the same gap C8's flag canonicalization had to close.

SCOPE. This gates the INDEX array only. ``test_layer_gate.py`` pins that the
flat-array *roughness* surface stays permissive; that decision is untouched.
Indices were already gated here for non-finite values and for ``|n|**2``
overflow -- both "no correction exists, and the damage is unrecoverable from
the output", which is exactly what gain is.
"""
from __future__ import annotations

import numpy as np
import pytest

import navette._smatrix as _native
from navette.smatrix.smatrix import Request, ScatterMatrix

WLS = np.array([500.0, 550.0, 600.0])
THICK = [0.0, 2000.0, 0.0]
_LOSS = np.array([1.0, 2.35 + 0.01j, 1.52])
_GAIN = np.array([1.0, 2.35 - 0.01j, 1.52])


def _raw_core_engine(k: float):
    """The raw pyfunction, which enters `Solver::from_wav_major_flat`."""
    cache = []
    for _ in WLS:
        for n in (complex(1.0, 0.0), complex(2.35, k), complex(1.52, 0.0)):
            cache += [n.real, n.imag]
    return _native.core_engine(
        WLS, np.array([0.0]), 3, np.array(cache),
        np.array(THICK, float), np.zeros(3, np.int32),
        np.zeros(3, np.int32), np.zeros(3), 0,
        int(Request.RS | Request.TS),
    )


# ---- the refusal -----------------------------------------------------------
def test_the_python_door_refuses_gain():
    with pytest.raises(ValueError) as exc:
        ScatterMatrix(_GAIN, THICK, wavelengths=WLS, angles=[0.0])
    msg = str(exc.value)
    assert "`layer_indices`" in msg, "must name the argument the caller passed"
    assert "layer 1" in msg, "must name the row, not just the fact"
    assert "wavelength index 0" in msg
    assert "3 of the index grid's values" in msg, (
        "a count separates one stray wavelength from a whole bad material"
    )
    assert "optical gain" in msg
    assert msg.isascii(), "the message must survive a cp1252 console"


def test_the_engine_refuses_gain_too_unlike_the_cross_channel():
    """The deliberate difference from C2, pinned so it stays deliberate.

    ``test_the_engine_still_computes_what_the_doors_refuse`` pins the opposite
    for C2: the raw pyfunction must keep computing Mode A's cross channel,
    because a legacy parity port depends on those exact numbers. Nothing
    depends on a gain layer's output -- the two paths do not even agree with
    each other -- so here the raw pyfunction refuses as well.
    """
    assert _raw_core_engine(+0.01), "loss is the control and must still solve"
    with pytest.raises(ValueError) as exc:
        _raw_core_engine(-0.01)
    msg = str(exc.value)
    assert "Solver" in msg and "layer 1" in msg, msg
    assert "optical gain" in msg, msg


def test_gain_anywhere_in_the_stack_refuses_including_the_ambient():
    """The ambient is checked too, and the order matters.

    ``Im(n[0]) != 0`` is otherwise *dropped* with a warning (R3.1), and that
    rule counts a negative imaginary part as absorption. Left alone it would
    have swallowed a sign error silently, reporting that it had dropped
    absorption that was never there. The gain gate runs first.
    """
    for row in (0, 1, 2):
        n = np.array([1.0, 2.35, 1.52], dtype=complex)
        n[row] -= 0.01j
        with pytest.raises(ValueError, match="optical gain"):
            ScatterMatrix(n, THICK, wavelengths=WLS, angles=[0.0])


def test_the_needle_door_refuses_a_gain_needle_material():
    """The needle material arrives separately and would otherwise miss the gate.

    The host stack cannot carry gain past ``ScatterMatrix``, but the needle's
    own index is handed to ``needle_gradient`` directly. It goes into the same
    conjugating kernels (``needle_operator.rs:248``) and its spacer tau takes
    the same clamp (``:1290``), so a gain needle is mangled exactly like a
    gain layer -- and it would be a strange rule that refused the stack being
    inserted into but not the thing being inserted.
    """
    from navette.smatrix.needle import NeedleRequest, needle_gradient

    idx = np.array([1.0 + 0j, 2.35 + 0j, 1.45 + 0j, 1.52 + 0j])
    thick = np.array([0.0, 40.0, 80.0, 0.0])
    st = ScatterMatrix(idx, thick, wavelengths=WLS, angles=[0.0])
    z = [40.0 + 30.0]

    ok = needle_gradient(st, np.full(WLS.size, 1.9 + 0j), z, NeedleRequest.P)
    assert ok, "the control must produce a gradient"

    with pytest.raises(ValueError) as exc:
        needle_gradient(st, np.full(WLS.size, 1.9 - 0.02j), z, NeedleRequest.P)
    msg = str(exc.value)
    assert "needle material" in msg, msg
    assert "optical gain" in msg, msg
    assert msg.isascii()


# ---- what must NOT refuse --------------------------------------------------
def test_loss_is_the_control_and_still_solves():
    out = ScatterMatrix(_LOSS, THICK, wavelengths=WLS, angles=[0.0]).compute(
        Request.RS | Request.TS | Request.A_S, squeeze=False
    )
    a = np.asarray(out["A_s"], float)
    assert np.all(a > 0.0), f"an absorbing stack must absorb: {a}"


def test_negative_zero_is_not_gain():
    """``-0.0 < 0`` is False in IEEE, and the engine's own branch rules agree.

    ``forward_branch`` and ``sanitize_incident_index`` both already treat a
    signed zero as zero, and a provider that writes ``-0.0`` for a transparent
    material is not describing an amplifier. Pinned because the obvious "tidy
    up the sign" edit would start refusing real grids.
    """
    n_neg = np.array([1.0, complex(2.35, -0.0), 1.52])
    n_pos = np.array([1.0, complex(2.35, 0.0), 1.52])
    req = Request.RS | Request.TS
    a = ScatterMatrix(n_neg, THICK, wavelengths=WLS, angles=[0.0]).compute(
        req, squeeze=False)
    b = ScatterMatrix(n_pos, THICK, wavelengths=WLS, angles=[0.0]).compute(
        req, squeeze=False)
    for key in a:
        assert np.array_equal(np.asarray(a[key]), np.asarray(b[key])), key


def test_a_transparent_stack_still_builds():
    ScatterMatrix(np.array([1.0, 2.35, 1.52]), THICK,
                  wavelengths=WLS, angles=[0.0])


# ---- one rule, three doors -------------------------------------------------
def test_both_doors_explain_gain_the_same_way():
    """Python ``ScatterMatrix``, the Rust ``Solver`` and the Structure door.

    Three separate implementations on purpose: the Python one can name the
    ``layer_indices`` argument the caller actually wrote, the Solver one
    reaches every flat-array entry including the raw FFI, and the Structure
    one fires during nominal expansion, before a solver array exists. The
    explanation is the thing that has to stay identical -- whichever is
    edited, this fails until the others follow. Same contract as
    ``test_both_doors_explain_it_the_same_way`` (ambient) and
    ``test_both_doors_explain_the_cross_channel_the_same_way`` (C2).

    The Structure door is the one that was already guarded when round 1 was
    written: it refused ``k < 0`` while the two flat-array doors did not.
    Its message said "check provider data" and nothing about what the solver
    would otherwise have done, so it now carries this text too.
    """
    from navette.structure import Layer, Navette_Structure, solve_structure
    from navette.structure.materials import DictMaterialProvider

    with pytest.raises(ValueError) as pexc:
        ScatterMatrix(_GAIN, THICK, wavelengths=WLS, angles=[0.0])
    py_msg = str(pexc.value)
    shared = py_msg[py_msg.index("Im(n) < 0 is optical gain"):]
    assert len(shared) > 600, f"parity slice looks truncated: {shared!r}"

    with pytest.raises(ValueError) as rexc:
        _raw_core_engine(-0.01)
    rust_msg = str(rexc.value)
    assert shared in rust_msg, (
        f"Solver door drifted / python: {shared!r} / rust: {rust_msg!r}"
    )

    mats = {"glass": np.full(WLS.size, 1.52 + 0j),
            "amplifier": np.full(WLS.size, 2.35 - 0.01j)}
    st = Navette_Structure(
        [Layer(0.0, "glass"), Layer(2000.0, "amplifier"), Layer(0.0, "glass")],
        {}, DictMaterialProvider(dict(mats), wavelength=WLS))
    with pytest.raises(ValueError) as sexc:
        solve_structure(st, WLS, 0.0, request=int(Request.RS))
    st_msg = str(sexc.value)
    assert "Nominal expansion produced k < 0" in st_msg, st_msg
    assert shared in st_msg, (
        f"Structure door drifted / python: {shared!r} / structure: {st_msg!r}"
    )
    assert rust_msg.isascii() and st_msg.isascii()


def test_the_explanation_says_what_each_path_does_and_how_to_fix_it():
    """A refusal that does not say why is a wall, not a message."""
    with pytest.raises(ValueError) as exc:
        ScatterMatrix(_GAIN, THICK, wavelengths=WLS, angles=[0.0])
    msg = str(exc.value)
    for fragment in (
        "conjugates",          # what the coherent path does
        "clamps",              # what the incoherent path does
        "-4.3e-5",             # the measured size of the incoherent damage
        "POSITIVE",            # the sign that makes it unrecoverable
        "exp(+i*w*t)",         # the way out for a convention mismatch
        "-0.0 is not gain",    # the edge the caller will hit next
    ):
        assert fragment in msg, f"missing {fragment!r} from: {msg}"
