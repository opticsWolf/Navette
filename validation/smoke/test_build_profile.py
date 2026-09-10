# -*- coding: utf-8 -*-
"""The build-profile probe and the bench release gate (R2.1).

Deliberately does NOT assert that the installed extension is a release build:
a contributor running the tests against a debug build is fine, a contributor
*timing* one is not, so the profile assertion belongs at bench time. What is
tested here is that the mechanism works -- the probe exists and reports a real
profile, and the gate actually refuses a debug one.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

import navette

_BENCHES = Path(__file__).resolve().parents[1] / "benches"


def _bench_common():
    spec = importlib.util.spec_from_file_location(
        "_bench_common_under_test", _BENCHES / "_bench_common.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_build_profile_reports_a_real_profile():
    assert navette.build_profile() in ("release", "debug")


def test_build_profile_matches_the_native_symbol():
    """The package-level accessor must not drift from the extension."""
    from navette._navette import build_profile as native

    assert navette.build_profile() == native()


def test_build_profile_is_exported():
    assert "build_profile" in navette.__all__


def test_require_release_accepts_release():
    assert _bench_common().require_release(lambda: "release") == "release"


@pytest.mark.parametrize("profile", ["debug", "dev", ""])
def test_require_release_refuses_anything_else(profile):
    """The gate must exit, not warn: a bench that prints numbers from a debug
    build is worse than one that prints nothing."""
    with pytest.raises(SystemExit) as exc:
        _bench_common().require_release(lambda: profile)
    assert "REFUSING TO BENCHMARK" in str(exc.value)
    assert "--release" in str(exc.value)


def test_require_release_refuses_an_unbuilt_extension():
    def missing():
        raise ImportError("No module named 'navette._navette'")

    with pytest.raises(SystemExit) as exc:
        _bench_common().require_release(missing)
    assert "REFUSING TO BENCHMARK" in str(exc.value)


def test_every_bench_gates_on_the_release_build():
    """No bench may skip the preamble -- that is how the debug build got in.

    A new bench added without ``require_release()`` fails here rather than
    quietly producing numbers nobody can trust.
    """
    benches = [p for p in sorted(_BENCHES.rglob("*bench*.py"))
               if p.name != "_bench_common.py"
               and "refs" not in p.parts
               and "gen" not in p.parts]
    # Guard the guard: an empty scan would pass vacuously.
    assert len(benches) >= 7, f"bench discovery broke: found {benches}"

    offenders = [p.relative_to(_BENCHES).as_posix() for p in benches
                 if "require_release()" not in p.read_text(encoding="utf-8")]
    assert not offenders, f"benches missing require_release(): {offenders}"


def test_bench_provenance_stamps_the_profile():
    prov = _bench_common().bench_provenance()
    assert prov["build_profile"] == navette.build_profile()
    assert prov["navette_version"] == navette.__version__


def test_setup_bench_is_idempotent_and_puts_src_on_path():
    """A bench must import ``navette`` from any cwd, and repeated calls must
    not grow ``sys.path``.

    Counted as a delta, not an absolute: under pytest, ``src`` is typically
    already present twice (``pytest.ini``'s ``pythonpath = src`` plus the
    editable install's ``.pth``), which says nothing about this function.
    """
    mod = _bench_common()
    src = str(Path(__file__).resolve().parents[2] / "src")
    saved = list(sys.path)
    try:
        before = sys.path.count(src)
        mod.setup_bench()
        once = sys.path.count(src)
        mod.setup_bench()
        twice = sys.path.count(src)
    finally:
        sys.path[:] = saved

    assert once >= 1, "setup_bench() did not make src importable"
    assert once - before <= 1, "setup_bench() duplicated an existing entry"
    assert twice == once, "setup_bench() is not idempotent"


def test_setup_bench_survives_an_unreconfigurable_stream(monkeypatch):
    """Never take a bench down over console encoding -- warn-free no-op."""

    class Stubborn:
        def reconfigure(self, **kw):
            raise ValueError("underlying buffer has been detached")

    monkeypatch.setattr(sys, "stdout", Stubborn())
    monkeypatch.setattr(sys, "stderr", Stubborn())
    _bench_common().setup_bench()
