# SPDX-License-Identifier: LGPL-3.0-or-later
"""C2/C9 -- the front-block cross channel is refused at the doors, not fixed.

Mode A takes the p-s cross channel from the FIRST coherent block while its
intensities are totals over every incoherent echo. Asking for Delta, DOP,
S2/S3, the retardance or the raw cross terms on a stack with an interior
incoherent flag therefore hands back a ratio of two different stacks: DOP_R
arrives as ``|rs_c|^2/Rs`` instead of 1 for a stack that does not depolarize.

Mode A's numbers are a deliberate bit-for-bit legacy port, so they are NOT
changed. The refusal sits at the two doors ABOVE the port's entry -- the
`ScatterMatrix` class and `solve_arrays` -- while the port itself drives the
raw ``core_engine`` pyfunction straight into ``Solver::solve`` and still
computes. `test_the_engine_still_computes_what_the_doors_refuse` pins that
the door check did not leak inward.
"""
import numpy as np
import pytest

from navette.smatrix.smatrix import CoherenceMode, Request, ScatterMatrix

WLS = np.array([550.0])
ANG = np.array([0.0])

# air / H / thick spacer / L / glass -- row 2 is the flagged one.
IDX = np.array([[1.0 + 0j], [2.35 + 0j], [1.5 + 0j], [1.46 + 0j], [1.52 + 0j]])
THICK = np.array([0.0, 95.0, 50_000.0, 110.0, 0.0])
INTERIOR = np.array([0, 0, 1, 0, 0], dtype=np.int32)
NONE = np.zeros(5, dtype=np.int32)
HALF_SPACES = np.array([1, 0, 0, 0, 1], dtype=np.int32)


def _sm(flags, mode=CoherenceMode.FRONT_BLOCK):
    return ScatterMatrix(IDX, THICK, wavelengths=WLS, angles=ANG,
                         incoherent_flags=flags, coherence_mode=mode)


# ---- the refusal ------------------------------------------------------------
@pytest.mark.parametrize("bits", [
    Request.DOP_R, Request.DOP_T, Request.DELTA_R, Request.DELTA_T,
    Request.S2_R, Request.S3_R, Request.S2_T, Request.S3_T,
    Request.CROSS_R, Request.CROSS_T, Request.RETARD_R, Request.RETARD_T,
])
def test_every_cross_bit_is_refused_under_mode_a_with_an_interior_flag(bits):
    """All twelve NEEDS_CROSS bits, one rule -- no 'raw cross is fine' split.

    The raw channel is the same defective object, and letting it through
    would let a caller rebuild the broken DOP by hand.
    """
    with pytest.raises(ValueError, match="cross-channel observable"):
        _sm(INTERIOR).compute(bits)


def test_the_refusal_names_the_row_and_the_way_out():
    with pytest.raises(ValueError) as exc:
        _sm(INTERIOR).ellipsometry()
    msg = str(exc.value)
    assert "row(s) 2" in msg, msg
    assert "coherence_mode=1" in msg, msg
    assert "front_block" in msg, msg
    assert msg.isascii(), "the message must survive a cp1252 console"


# ---- the controls: exactly the affected population, nothing else ------------
def test_mode_b_is_allowed_and_is_the_physical_answer():
    """Mode B cascades the cross channel with the echoes.

    On a stack of isotropic layers at normal incidence nothing depolarizes,
    so the true DOP is exactly 1 -- which is the number Mode B returns and
    the one Mode A cannot.
    """
    out = _sm(INTERIOR, CoherenceMode.COHERENCY_MATRIX).ellipsometry()
    assert float(np.ravel(out["DOP_R"])[0]) == pytest.approx(1.0, abs=1e-12)


def test_mode_c_is_allowed():
    """FULLY_COHERENT is one block over the whole stack; its cross is correct."""
    _sm(INTERIOR, CoherenceMode.FULLY_COHERENT).ellipsometry()


def test_no_flag_is_allowed():
    """With nothing flagged, A and B are bit-identical -- no hazard, no refusal."""
    a = _sm(NONE).ellipsometry()
    b = _sm(NONE, CoherenceMode.COHERENCY_MATRIX).ellipsometry()
    for k in ("DOP_R", "Delta_R", "Rs"):
        assert np.array_equal(np.asarray(a[k]), np.asarray(b[k])), k


def test_half_space_flags_alone_are_allowed():
    """C5: rows 0 and last are half-spaces the sweep never consults.

    Flagging only those is a no-op on every door, so refusing it would
    contradict C5 rather than implement C2.
    """
    _sm(HALF_SPACES).ellipsometry()


def test_intensities_are_never_refused():
    """R/T/A are totals and correct under Mode A even with a flag."""
    out = _sm(INTERIOR).reflectance_transmittance()
    assert np.isfinite(np.asarray(out["Rs"])).all()


# ---- the engine underneath is untouched ------------------------------------
def test_the_engine_still_computes_what_the_doors_refuse():
    """The legacy port enters below the doors and must keep working.

    `test_core_engine_rigorous_ellipsometry` drives Mode A + interior flags +
    cross observables through the raw ``core_engine`` pyfunction, which builds
    a Solver and calls ``Solver::solve`` directly. If the door check ever
    leaks into the engine, that bit-for-bit port dies -- this fails first.
    """
    import navette._smatrix as native

    n_layers, n_wavs = 5, 1
    cache = np.empty(n_layers * n_wavs * 2, dtype=np.float64)
    for row in range(n_layers):
        cache[row * 2] = IDX[row, 0].real
        cache[row * 2 + 1] = IDX[row, 0].imag
    req = int(Request.DOP_R | Request.CROSS_R | Request.CROSS_T | Request.RS)
    out = native.core_engine(
        WLS, np.sin(np.radians(ANG)), n_layers, cache, THICK,
        INTERIOR, np.zeros(n_layers, dtype=np.int32),
        np.zeros(n_layers, dtype=np.float64),
        int(CoherenceMode.FRONT_BLOCK), req,
    )
    assert np.isfinite(np.asarray(out["DOP_R"])).all()


def test_mode_a_and_mode_b_cross_terms_are_different_objects():
    """C9. cross_R is the front block; cross_T is a third thing again.

    Mode A's transmitted cross term is a product of per-block
    ``t_p*conj(t_s)`` across the joins with no bounce series -- not the front
    block and not the Mode B cascade. Both are real numbers the engine will
    hand back; the point of the record is that they are not the same object.
    """
    import navette._smatrix as native

    n_layers, n_wavs = 5, 1
    cache = np.empty(n_layers * n_wavs * 2, dtype=np.float64)
    for row in range(n_layers):
        cache[row * 2] = IDX[row, 0].real
        cache[row * 2 + 1] = IDX[row, 0].imag
    req = int(Request.CROSS_R | Request.CROSS_T | Request.RS | Request.TS)

    def run(mode):
        return native.core_engine(
            WLS, np.sin(np.radians(ANG)), n_layers, cache, THICK,
            INTERIOR, np.zeros(n_layers, dtype=np.int32),
            np.zeros(n_layers, dtype=np.float64), int(mode), req,
        )

    a = run(CoherenceMode.FRONT_BLOCK)
    b = run(CoherenceMode.COHERENCY_MATRIX)
    assert not np.array_equal(np.asarray(a["cross_T"]), np.asarray(b["cross_T"])), (
        "Mode A's cross_T is supposed to be a different object from Mode B's"
    )
    # The photometric identity the remediation plan pins: A and B agree
    # bit-for-bit on intensities, and only the cross channel differs.
    for k in ("Rs", "Ts"):
        assert np.array_equal(np.asarray(a[k]), np.asarray(b[k])), k


# ---- one rule, two doors ---------------------------------------------------
def test_both_doors_explain_the_cross_channel_the_same_way():
    """`ScatterMatrix` (Python) and `solve_arrays` (Rust) carry one text.

    They are separate implementations on purpose -- the Python one can point
    ``stacklevel`` at the caller and a Rust one cannot -- so the explanation
    is the thing that has to stay identical. Whichever is edited, this fails
    until the other follows. Same contract as
    `test_both_doors_explain_it_the_same_way` for the ambient rule.
    """
    from navette.structure import Layer, Navette_Structure, solve_structure
    from navette.structure.materials import DictMaterialProvider

    with pytest.raises(ValueError) as exc:
        _sm(INTERIOR).ellipsometry()
    py_msg = str(exc.value)
    shared = py_msg[py_msg.index("Mode A (front_block)"):]
    assert len(shared) > 300, f"parity slice looks truncated: {shared!r}"

    mats = {"glass": np.full(1, 1.52 + 0j), "spacer": np.full(1, 1.5 + 0j)}
    st = Navette_Structure(
        [Layer(0.0, "glass"),
         Layer(50_000.0, "spacer", coherent=False),
         Layer(0.0, "glass")],
        {}, DictMaterialProvider(dict(mats), wavelength=WLS))
    with pytest.raises(ValueError) as rexc:
        solve_structure(st, WLS, 0.0, request=int(Request.DOP_R))
    rust_msg = str(rexc.value)
    assert shared in rust_msg, (
        f"the two doors have drifted apart / python: {shared!r} / "
        f"rust: {rust_msg!r}"
    )
    assert rust_msg.isascii(), "the message must survive a cp1252 console"
