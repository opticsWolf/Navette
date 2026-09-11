# -*- coding: utf-8 -*-
"""Cross-language sync guards for the request bits and the schema versions.

Three constant tables are written twice, once per language, and nothing tied
them together:

  * 49 ``REQ_*`` solver bits -- ``core_engine.rs`` vs ``Request`` in
    ``smatrix.py``. Hand-copied, and not exported by the extension.
  * 18 ``NREQ_*`` needle bits -- ``needle_engine.rs``, re-exported by
    ``navette-py``, and bound (not copied) into ``NeedleRequest``.
  * two schema versions -- ``PROGRAM_SCHEMA_VERSION`` (``config.rs`` vs
    ``config/program.py``) and ``SCHEMA_VERSION`` (``structure/version.rs``
    vs ``structure/types.py``).

A drift in any of them is silent: the caller asks for one observable and is
handed another, with the right shape and the wrong physics. The tests here
close that off from both ends --

  * **source-level**, parsing the Rust constants and comparing name-for-name
    and bit-for-bit, which catches a constant added on one side only; and
  * **behavioural**, driving each bit through the real engine and asserting
    the returned channel keys, which catches a renumbering that both source
    tables happen to agree on being wrong about.

The source-level half skips when the Rust tree is absent (installed wheel);
the behavioural half always runs.
"""

import re
from pathlib import Path

import numpy as np
import pytest

from navette.smatrix import smatrix as _sm
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.smatrix.needle import NeedleRequest

_RUST_SRC = Path(__file__).resolve().parents[2] / "rust"
_CORE_ENGINE = _RUST_SRC / "navette" / "src" / "smatrix" / "core_engine.rs"
_NEEDLE_ENGINE = _RUST_SRC / "navette" / "src" / "smatrix" / "needle_engine.rs"

_CONST_RE = re.compile(
    r"^pub const (?P<name>[A-Z0-9_]+): u64 = 1 << (?P<bit>\d+);", re.M)


def _rust_bits(path, prefix):
    """``{NAME_TAIL: value}`` for every ``pub const <prefix>_X: u64 = 1 << n``."""
    if not path.is_file():
        pytest.skip(f"{path.name} not present (installed wheel, not a checkout)")
    text = path.read_text(encoding="utf-8")
    out = {}
    for m in _CONST_RE.finditer(text):
        name = m.group("name")
        if not name.startswith(prefix + "_"):
            continue
        out[name[len(prefix) + 1:]] = 1 << int(m.group("bit"))
    assert out, f"no {prefix}_* constants parsed from {path.name} -- regex rotted"
    return out


def _single_bits(flag_cls):
    """Members that name exactly one bit (i.e. not the convenience bundles)."""
    return {m.name: m.value for m in flag_cls
            if m.value != 0 and m.value & (m.value - 1) == 0}


# --------------------------------------------------------------------------
# 1. REQ_* -- source-level, exact
# --------------------------------------------------------------------------

def test_request_bits_match_the_rust_constants():
    """``Request`` must be name-for-name and bit-for-bit ``core_engine.rs``.

    A bit added in Rust and forgotten in Python is unreachable; a bit added in
    Python and missing in Rust is silently ignored by the engine, which then
    returns a dict without that channel. Both directions fail here.
    """
    rust = _rust_bits(_CORE_ENGINE, "REQ")
    py = _single_bits(Request)

    assert set(rust) - set(py) == set(), (
        f"REQ_* in core_engine.rs with no Request member: "
        f"{sorted(set(rust) - set(py))}")
    assert set(py) - set(rust) == set(), (
        f"Request members with no REQ_* in core_engine.rs: "
        f"{sorted(set(py) - set(rust))}")
    mismatched = {k: (py[k], rust[k]) for k in rust if py[k] != rust[k]}
    assert not mismatched, f"bit value drift (python, rust): {mismatched}"


def test_request_bits_are_dense_and_unique():
    """No duplicate and no skipped bit positions.

    Two members sharing a value is the failure mode that survives every
    name-based check: the engine cannot tell the two requests apart.
    """
    py = _single_bits(Request)
    positions = sorted(v.bit_length() - 1 for v in py.values())
    assert len(set(positions)) == len(positions), (
        f"duplicate bit positions in Request: {positions}")
    assert positions == list(range(len(positions))), (
        f"gap in Request bit positions: {positions}")


def test_every_request_bit_has_a_declared_output_key():
    """The bit -> channel-key map must cover every bit, with no strays."""
    py = set(_single_bits(Request))
    mapped = {m.name for m in _sm._SCALAR_KEYS} | {m.name for m in _sm._DISP_KEYS}
    assert py - mapped == set(), f"Request bits with no output key: {sorted(py - mapped)}"
    assert mapped - py == set(), f"output keys for unknown bits: {sorted(mapped - py)}"

    keys = list(_sm._SCALAR_KEYS.values())
    keys += [k for v in _sm._DISP_KEYS.values() for k in v]
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    assert not dupes, f"two request bits claim the same output key: {dupes}"


# --------------------------------------------------------------------------
# 2. REQ_* -- behavioural, through the engine
# --------------------------------------------------------------------------

_N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0.001j, 1.52 + 0j])
_D = np.array([0.0, 120.0, 200.0, 0.0])
_WLS = np.linspace(400.0, 800.0, 9)   # >= 5 points: the DISP_* bits differentiate
_ANGLES = [0.0, 45.0]


def _stack():
    return ScatterMatrix(_N, _D, wavelengths=_WLS, angles=_ANGLES)


_SINGLE = sorted(
    (m for m in Request if m.value and m.value & (m.value - 1) == 0),
    key=lambda m: m.value)


@pytest.mark.parametrize("bit", _SINGLE, ids=lambda m: m.name)
def test_single_bit_returns_exactly_its_own_channels(bit):
    """Ask for one bit; get exactly the channels that bit declares.

    This is the half that a source-level comparison cannot do. If Python's
    ``Request.RS`` and Rust's ``REQ_RS`` were renumbered in step -- both
    tables self-consistent, both wrong relative to the compute kernel's
    dispatch -- the source test passes and this one fails, because the engine
    would hand back the neighbouring channel under the requested name.
    """
    got = set(_stack().compute(bit, squeeze=False))
    expected = set(_sm.expected_keys(bit))
    assert got == expected, (
        f"{bit.name} (1 << {bit.value.bit_length() - 1}): "
        f"missing={sorted(expected - got)} unexpected={sorted(got - expected)}")


@pytest.mark.parametrize("bundle", [
    Request.PHOTOMETRY, Request.ELLIPSOMETRY, Request.ABSORPTION,
    Request.BACKSIDE, Request.STOKES_R, Request.STOKES_T,
], ids=lambda m: m.name or "bundle")
def test_convenience_bundles_are_unions_of_their_bits(bundle):
    """A bundle must be exactly the union of the bits it ORs together."""
    members = [m for m in _SINGLE if bundle & m]
    expected = set()
    for m in members:
        expected |= set(_sm.expected_keys(m))
    assert set(_stack().compute(bundle, squeeze=False)) == expected


def test_all_bits_at_once_returns_every_channel():
    """The union of all 49 bits emits the union of all declared keys.

    Requesting everything at once exercises the shared intermediates the
    single-bit calls skip, so a bit that only works in isolation fails here.
    """
    every = Request(0)
    for m in _SINGLE:
        every |= m
    expected = set()
    for m in _SINGLE:
        expected |= set(_sm.expected_keys(m))
    got = set(_stack().compute(every, squeeze=False))
    assert got == expected, (
        f"missing={sorted(expected - got)} unexpected={sorted(got - expected)}")


# --------------------------------------------------------------------------
# 3. NREQ_* -- three-way: Rust source, native export, Python flag
# --------------------------------------------------------------------------

def _native_nreq():
    import navette._smatrix as native
    return {n[len("NREQ_"):]: getattr(native, n)
            for n in dir(native) if n.startswith("NREQ_")}


def test_needle_bits_match_rust_native_and_python():
    """``needle_engine.rs`` -> extension exports -> ``NeedleRequest``.

    ``NeedleRequest`` *binds* the native values rather than copying them, so
    the values cannot drift. What can drift is the membership: a constant
    added in Rust and exported but never surfaced in the IntFlag is
    unreachable from Python, and nothing else notices.
    """
    rust = _rust_bits(_NEEDLE_ENGINE, "NREQ")
    native = _native_nreq()
    py = _single_bits(NeedleRequest)

    assert set(native) == set(rust), (
        f"extension exports vs needle_engine.rs: "
        f"only-rust={sorted(set(rust) - set(native))} "
        f"only-native={sorted(set(native) - set(rust))}")
    assert set(py) == set(rust), (
        f"NeedleRequest vs needle_engine.rs: "
        f"only-rust={sorted(set(rust) - set(py))} "
        f"only-python={sorted(set(py) - set(rust))}")
    bad = {k: (py[k], native[k], rust[k]) for k in rust
           if not (py[k] == native[k] == rust[k])}
    assert not bad, f"value drift (python, native, rust): {bad}"


def test_needle_bits_are_dense_and_unique():
    py = _single_bits(NeedleRequest)
    positions = sorted(v.bit_length() - 1 for v in py.values())
    assert len(set(positions)) == len(positions), f"duplicates: {positions}"
    assert positions == list(range(len(positions))), f"gap: {positions}"


# --------------------------------------------------------------------------
# 4. Schema versions -- both gates live natively, so probe them
# --------------------------------------------------------------------------

def _accepted_version(gate, lo=0, hi=8):
    """The one version a native gate accepts, discovered by probing it."""
    ok = []
    for v in range(lo, hi + 1):
        try:
            gate(v)
        except Exception:
            continue
        ok.append(v)
    assert len(ok) == 1, f"gate accepts {ok}, expected exactly one version"
    return ok[0]


def test_state_schema_version_matches_the_native_gate():
    """``structure.types.SCHEMA_VERSION`` vs the version Rust actually accepts.

    ``check_schema_version`` is a thin shim over ``structure/version.rs``, so
    probing it reads the Rust constant at runtime -- no source parsing, and it
    works from an installed wheel.
    """
    from navette.structure.types import SCHEMA_VERSION, check_schema_version

    native = _accepted_version(
        lambda v: check_schema_version({"schema_version": v}, "sync-probe"))
    assert SCHEMA_VERSION == native, (
        f"python SCHEMA_VERSION={SCHEMA_VERSION}, rust accepts {native}")


def test_state_schema_gate_refuses_an_untagged_state():
    """An untagged state is malformed, not legacy -- v1 has no past."""
    from navette.structure.types import check_schema_version
    with pytest.raises(ValueError, match="missing schema_version"):
        check_schema_version({}, "sync-probe")


def test_program_schema_version_matches_the_native_gate(tmp_path):
    """``config.program.PROGRAM_SCHEMA_VERSION`` vs ``config.rs``'s gate.

    The gate used to test ``Some(1)`` as a literal while quoting the constant
    in its error message, so a bump would have rejected exactly the version it
    claimed to read. It compares against the constant now, and this test is
    what holds that.
    """
    import json

    from navette.config import program

    doc = tmp_path / "materials.json"

    def gate(v):
        doc.write_text(json.dumps(
            {"kind": "materials", "schema_version": v, "materials": []}),
            encoding="utf-8")
        return program.load_document(doc)

    native = _accepted_version(gate)
    assert program.PROGRAM_SCHEMA_VERSION == native, (
        f"python PROGRAM_SCHEMA_VERSION={program.PROGRAM_SCHEMA_VERSION}, "
        f"rust accepts {native}")


def test_program_schema_gate_refuses_an_untagged_document(tmp_path):
    """No envelope version means no envelope contract."""
    import json

    from navette.config import program

    doc = tmp_path / "materials.json"
    doc.write_text(json.dumps({"kind": "materials", "materials": []}),
                   encoding="utf-8")
    with pytest.raises(ValueError, match="schema_version"):
        program.load_document(doc)
