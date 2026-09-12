# SPDX-License-Identifier: LGPL-3.0-or-later
"""Toolchain lint: is the local clippy as new as the one CI runs?

`cargo clippy` only reports the lints its own version knows about. A local run
can be genuinely clean and CI still red, with nothing in the local battery able
to show it.

That is not hypothetical. Between 0.6.13 and 0.6.30 the local toolchain was
rustc 1.97.1 while CI installed `dtolnay/rust-toolchain@stable` (1.98.0). The
lint `chunks_exact_to_as_chunks` shipped in 1.98 and fired on one line of
`solver.rs`. CI was red for **17 consecutive pushes** and every local
verification passed. Because the `rust` job aborts at the clippy step, every
step after it -- the feature-gated builds, and the rustfmt gate added in
0.6.30 -- had never run on CI at all.

Usage: python tools/check_toolchain.py
Exit 1 when the active toolchain is older than the newest one installed, or
older than what the CI workflow declares.
"""
import re
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
CI = ROOT / ".github" / "workflows" / "ci.yml"


def run(*args: str) -> str:
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError) as exc:
        return f"<unavailable: {exc}>"
    return (p.stdout or p.stderr).strip()


def version_of(text: str) -> tuple[int, ...] | None:
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", text)
    return tuple(int(g) for g in m.groups()) if m else None


def fmt(v: tuple[int, ...] | None) -> str:
    return ".".join(str(x) for x in v) if v else "unknown"


def main() -> int:
    active_rustc = version_of(run("rustc", "--version"))
    active_clippy = version_of(run("cargo", "clippy", "--version"))

    # Newest toolchain rustup has on disk -- what the developer *could* be
    # linting with. `rustup toolchain list` prints e.g.
    # "1.98.0-x86_64-pc-windows-msvc" and "stable-... (active, default)".
    installed = []
    for line in run("rustup", "toolchain", "list").splitlines():
        v = version_of(line)
        if v:
            installed.append(v)
    newest = max(installed, default=None)

    print(f"active rustc   : {fmt(active_rustc)}")
    print(f"active clippy  : {fmt(active_clippy)}")
    print(f"newest local   : {fmt(newest)}")

    # What CI declares. `dtolnay/rust-toolchain@stable` means "whatever stable
    # is on the day the job runs", which is a moving target this script cannot
    # resolve offline -- so it only reports the channel and leans on the
    # installed-toolchain comparison for the hard check.
    channel = "unknown"
    if CI.is_file():
        text = CI.read_text(encoding="utf-8")
        m = re.search(r"dtolnay/rust-toolchain@(\S+)", text)
        if m:
            channel = m.group(1)
    print(f"CI toolchain   : dtolnay/rust-toolchain@{channel}")

    if active_rustc is None:
        print("\nFAIL: no working rustc on PATH.")
        return 1

    if newest and active_rustc < newest:
        print(
            f"\nFAIL: linting with {fmt(active_rustc)} while {fmt(newest)} is "
            f"installed.\n"
            f"      Clippy reports only the lints its own version knows, so a "
            f"clean run here\n"
            f"      proves nothing about CI. Either:\n"
            f"        rustup update stable\n"
            f"      or run the gate explicitly against the newer toolchain:\n"
            f"        cargo +{fmt(newest)} clippy --workspace --all-targets -- -D warnings\n"
            f"        cargo +{fmt(newest)} fmt --all --check"
        )
        return 1

    if channel == "stable":
        print(
            "\nOK. Note: CI tracks `stable`, which moves. Run `rustup update "
            "stable`\nbefore trusting a green clippy near a Rust release."
        )
    else:
        print("\nOK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
