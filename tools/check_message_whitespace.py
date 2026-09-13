# SPDX-License-Identifier: LGPL-3.0-or-later
"""Message hygiene scanner: lost line-continuations and forbidden non-ASCII.

`format!` message literals are user-facing error text. When a literal is
line-wrapped in source, a lost `\\` continuation leaves the indent
whitespace of the next line INSIDE the string: the user sees a run of
spaces mid-sentence. That happened for real - two refusals shipped at
0.6.34 (F0.2) with a 22-space run, one at 0.6.36 (F1.1) with 18, and two
older sites predate the series (amendment 3, findings 3 + V1). The exact
-format! twins added at C3 pin the five known messages; this tool guards
the class, because the next one will be a message nobody pinned.

Ground rule 7 (ASCII error messages) exists because of cp1252 consoles:
a byte the console cannot encode is where the rule actually bites. So
the non-ASCII check is two-staged (C3/V2):

- Space runs of 3+ inside string literals: **blocking**, allowlist
  empty. That class is unambiguous - a run of three spaces inside a
  message is always a lost `\\`.
- Non-ASCII: blocking for characters cp1252 cannot encode (e.g. `>=`,
  lambda), allowlisted advisory for the encodable ones (`em dash`, `x`,
  `cdot`) - F3.1's exposure re-audit reads every one of those messages
  anyway and retires the list there, at which point the tool goes fully
  blocking.

Scope: the trailing `#[cfg(test)] mod ... {` region and `rust/*/tests/`
are skipped (test assert messages are not user-facing; their wrapping is
nobody's problem), and `//`/`///`/`/* */` comments are never scanned (a
quoted string inside a trailing comment is not a message). The cutoff is
the trailing INLINE test module only (`_test_region_start`) - `cfg(test)`
attributes on helpers inside production impl blocks and out-of-line test
modules are not cutoffs; see the helper's docstring for the measured
counts.

Usage: python tools/check_message_whitespace.py
Exit 1 on any blocking finding.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Literal-aware scan: file sets that are test code in their entirety.
TEST_DIRS = ("rust/navette/tests", "rust/navette-py/tests")

# cp1252-encodable non-ASCII still outstanding in production messages.
# Each entry names the FILE, and the granularity is the whole file: a
# NEWLY INTRODUCED encodable non-ASCII character anywhere in an
# allowlisted file passes unremarked - acceptable only while the list
# shrinks toward empty, which is why F3.1 retires it. Measured at C3
# (review C/E4 corrected the count): SEVEN files, ten literals -
# tables.rs and solver.rs carry an em dash and a multiplication sign,
# color_merit.rs and synthesis_pipeline.rs em dashes, jacobian.rs middle
# dots and an em dash, optimizer.rs and trf.rs em dashes. All
# cp1252-encodable, all advisory until F3.1's re-audit.
NON_ASCII_ALLOWLIST = {
    "rust/navette/src/color/tables.rs",
    "rust/navette/src/smatrix/solver.rs",
    "rust/navette/src/smatrix/synthesis/color_merit.rs",
    "rust/navette/src/smatrix/synthesis/jacobian.rs",
    "rust/navette/src/smatrix/synthesis/optimizer.rs",
    "rust/navette/src/smatrix/synthesis/trf.rs",
    "rust/navette-py/src/synthesis_pipeline.rs",
}


def is_encodable(char: str) -> bool:
    try:
        char.encode("cp1252")
        return True
    except UnicodeEncodeError:
        return False


def _test_region_start(lines: list[str]) -> int | None:
    """1-based line of the trailing `#[cfg(test)] mod ... {` block, or None.

    NOT the first `#[cfg(test)]` in the file: that attribute also sits on
    test-only helpers inside production impl blocks (structure.rs:1185,
    solver.rs:831, providers.rs:321) and on out-of-line test modules
    (color/mod.rs:83, `mod parity;`). Cutting at the first one hid 63-75%
    of those files from the gate (review C, finding E2). Only an INLINE
    `mod ... {` at column 0 runs to EOF and is safe to cut at; a file
    with no such block is scanned whole.
    """
    for i in range(len(lines) - 1, -1, -1):
        if not lines[i].startswith("#[cfg(test)]"):
            continue
        j = i + 1
        while j < len(lines) and lines[j].startswith("#["):  # further attrs
            j += 1
        if (
            j < len(lines)
            and lines[j].startswith(("mod ", "pub mod "))
            and lines[j].rstrip().endswith("{")
        ):
            return i + 1
    return None


def scan_file(path: pathlib.Path):
    """Yield (line_no, kind, detail) findings for one .rs file.

    UNFILTERED: every string literal in the file is scanned, including
    test regions and comments' contents-adjacent literals - the caller
    applies the test-region rule, so the cutoff lives in one place.
    """
    text = path.read_text(encoding="utf-8", errors="replace")

    in_str = False
    in_line_comment = False
    cont_skip_ws = False  # `\` + newline: rustc strips the next line's indent
    content: list[str] = []
    start_line = 0
    findings = []

    def close_string(lineno: int):
        nonlocal content
        s = "".join(content)
        if "   " in s:
            findings.append((start_line, "space-run", repr(s)))
        non_ascii = sorted({c for c in s if ord(c) > 127})
        for c in non_ascii:
            kind = "non-ascii-nonencodable" if not is_encodable(c) else "non-ascii"
            findings.append((start_line, kind, f"{c!r} in {s[:60]!r}"))
        content = []

    i = 0
    lineno = 1
    raw = text
    while i < len(raw):
        c = raw[i]
        if c == "\n":
            lineno += 1
            if in_str:
                # A `\` at line end inside a literal is the continuation:
                # rustc strips the newline and the next line's indent.
                if raw[i - 1] == "\\":
                    cont_skip_ws = True
                else:
                    # A bare newline inside a normal string literal is a
                    # compile error, so reaching here means we mis-tracked
                    # (e.g. a quote inside a comment we did not model) -
                    # bail out of the string rather than cascade.
                    in_str = False
                    content = []
            in_line_comment = False
            i += 1
            continue
        if in_line_comment:
            i += 1
            continue
        if in_str:
            if cont_skip_ws and c == " ":
                i += 1
                continue
            cont_skip_ws = False
            if c == "\\":
                # skip the escaped char (may be the newline handled above)
                if i + 1 < len(raw):
                    if raw[i + 1] != "\n":
                        i += 1  # escaped quote/backslash/... not content
                i += 1
                continue
            if c == '"':
                close_string(lineno)
                in_str = False
                i += 1
                continue
            content.append(c)
            i += 1
            continue
        # not in string, not in line comment
        if c == "/" and i + 1 < len(raw):
            if raw[i + 1] == "/":
                in_line_comment = True
                i += 2
                continue
            if raw[i + 1] == "*":
                # skip to the closing */
                j = raw.find("*/", i + 2)
                if j == -1:
                    # unterminated: skip to EOF
                    break
                skipped_newlines = raw.count("\n", i, j + 2)
                lineno += skipped_newlines
                i = j + 2
                continue
        if c == '"':
            in_str = True
            content = []
            start_line = lineno
            i += 1
            continue
        i += 1

    return findings


def main() -> int:
    blocking = []
    advisory = []
    for base in ("rust",):
        for path in sorted((ROOT / base).rglob("*.rs")):
            rel = path.relative_to(ROOT).as_posix()
            if any(rel.startswith(d) for d in TEST_DIRS):
                continue
            findings = scan_file(path)
            if not findings:
                continue
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            test_from = _test_region_start(lines)
            for (ln, kind, detail) in findings:
                if test_from is not None and ln >= test_from:
                    continue  # test region: not user-facing
                if kind == "non-ascii":
                    if rel in NON_ASCII_ALLOWLIST:
                        advisory.append(f"{rel}:{ln} {kind}: {detail}")
                    else:
                        blocking.append(f"{rel}:{ln} {kind}: {detail}")
                else:
                    blocking.append(f"{rel}:{ln} {kind}: {detail}")

    if advisory:
        print("advisory (cp1252-encodable; retire at F3.1):")
        for a in advisory:
            print(f"  {ascii(a)}")
    if blocking:
        print("BLOCKING findings - fix the message, not the checker:")
        for b in blocking:
            print(f"  {ascii(b)}")
        print(
            "\nA 3+ space run inside a string literal is a lost "
            "`\\` line-continuation (wrap the SOURCE line, not the string). "
            "A non-encodable non-ASCII char breaks cp1252 consoles "
            "(ground rule 7)."
        )
        return 1
    print("message hygiene OK (space runs: 0; non-encodable non-ASCII: 0)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
