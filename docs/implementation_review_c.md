# Navette — C-series implementation review (amendment 3)

Review of the execution of `docs/implementation_amendment_3.md` (revision 2,
the `C*` worklist with its `V*` verification ledger) against the tree at
`823b7ab` (`dev_feature`, 0.6.44), 2026-09-13.

Commits under review:

| Item | Feature commit | Docs commit | CI run | Version |
|---|---|---|---|---|
| C1 | `e751852` | `d998381` | [34771251228](https://github.com/opticsWolf/Navette/actions/runs/34771251228) ✅ | none (gate repair) |
| C2 | `db3d5e1` | `0f679c1` | [34771889430](https://github.com/opticsWolf/Navette/actions/runs/34771889430) ✅ | 0.6.43 |
| C3 | `a91fc3d` | `823b7ab` | [34772770616](https://github.com/opticsWolf/Navette/actions/runs/34772770616) ✅ | 0.6.44 |

Method: read each item's specification in amendment 3, read the diff that
claims to implement it, ran the full ground-rule-3 battery locally, read
the CI logs for all three pushes, and settled the two behaviour claims
with independent probes rather than by reading the twins that ship with
the fix (a twin written by the same hand as the fix is evidence, not
proof). Both probes were built, run and deleted; the tree is unchanged
apart from this file.

Findings here use the prefix **E\*** — disjoint from the plan's `A*`/`N*`,
the amendment's `B*`/`C*`, the plan's `U*`/`R*`, the re-check pass's `V*`,
and the plan's `§D` design sections.

---

## 1. Verdict

**All three items are implemented as specified, including the nine `V*`
corrections, and the blocking finding is genuinely resolved.** The needle
pin now runs and **passes on both platforms in CI** — it is not silently
taking the new skip path, which was the obvious way this repair could have
looked green while gating nothing. The behaviour fix does what it claims:
the review's measured defect reproduces as fixed under an independent
probe, and the plain-film control is bitwise unmoved.

The whole battery is green here, and — for the first time in this series —
green in CI on a per-item push rather than on one batch. Amendment 3 §7's
push-per-item rule was honoured: `origin/dev_feature` and `HEAD` are equal,
and each of the three pushes has its own passing run.

**Four defects and two nits, all in the bookkeeping and tooling half of
C3.** None is in C1's re-specification or C2's behaviour change. The two
that matter are **E1** (the bookkeeping fix broke a neighbouring table row —
a fresh instance of the very rule it was fixing) and **E2** (the new
blocking gate's coverage is a fraction of what its docstring claims, in
the very file that produced the finding).

Nothing here blocks F2.1. E1 and E2 should land first, as one bookkeeping
commit, because both make the *next* audit harder rather than the current
code wrong.

---

## 2. What each item delivered

### 2.1 C1 — the needle pin, re-specified per platform

Specified and delivered, with the probe outcome recorded honestly.

| Required | In the tree |
|---|---|
| Platform-keyed `RECORDED_DIGESTS` | `test_needle_pin.py:145`, keyed by `platform.system()` |
| Provenance per entry | File header: tree commit, runner image, toolchain, build profile, probe run id |
| Loud skip on an unrecorded platform | `pytest.skip` naming the recording procedure and the amendment item |
| Failure message names the platform | `needle fingerprint moved on {platform}: recorded …, got …` |
| Runner images pinned | `python` matrix → `ubuntu-24.04` / `windows-2025`; `dependency floors` → `ubuntu-24.04` |
| `rust` job left alone | Yes — same-process comparisons, no recorded cross-machine constant |
| Ground rule 5 gains the sentence | `implementation_plan.md:100` |
| No version bump | Correct (R5's test-only precedent) |

**Outcome A confirmed, and confirmed the right way.** The recorded Linux
digest `5c309813…f547da` is the value the probe run printed, and it matches
the tip's Linux digest from the red run — so the series moved neither
platform and no bisect was needed. The claim is checkable after the fact
because the probe run id is written down.

**The skip path is not masking anything.** This was the failure mode worth
checking, and CI log for run 34772770616 settles it:

```
python (ubuntu-24.04)  test_needle_run_is_deterministic_and_matches_the_recorded_digest PASSED
python (windows-2025)  test_needle_run_is_deterministic_and_matches_the_recorded_digest PASSED
dependency floors      validation/regression/synthesis/test_needle_pin.py ..
```

Both platforms execute the digest comparison. The only skip in the whole
suite is the pre-existing `opt-minpack-lm` feature gate.

One deviation from the written spec, benign: the amendment specified a
`@pytest.mark.skipif` decorator; the implementation uses an in-body
`pytest.skip` with the same message. Behaviour is identical (the decorator
would report the skip at collection rather than at call). Not a finding.

### 2.2 C2 — the interface-carrying film clamps up

Specified and delivered, and the design is sounder than the amendment's
own snippet.

The predicate swap is at `structure.rs:1118` (`clamp_up_rows &&
sp.is_singleton_bulk()`), and the slice-aware arithmetic is exactly the
specified three lines. **The implementer caught something the amendment's
snippet did not cover:** the old branch hard-coded `end: p + 1` and
`p += 1`, which is correct only for a one-row span. Both were generalised
to `e = p + thin.len()` / `p = e`. Without that the span table would have
been corrupted for every clamped interface film — the snippet as written
in the amendment would have shipped a bug.

**The two branches are provably disjoint**, which is what makes the swap
safe and is worth recording because it is not obvious from the call site:

- `Span::is_singleton_bulk()` is `end - bulk_start == 1` — exactly one bulk row.
- `span_is_scalable_rows` (`structure.rs:1395`) requires `bulk >= 2`.

So a slice-carrying film can never enter F1.6's scale branch and have its
interface slice scaled with the bulk, which would have violated B1's
ownership rule. No ordering dependency between the two `if`s is load-bearing.

**Independent verification.** I rebuilt the review's original measured
scenario as a throwaway integration test, written against the public API
and without reference to the committed twins — 100 nm SiO2 lead, 0.8 nm
TiO2 carrying a 0.2 nm interface, 2.0 nm floor, `ClampUpFinal`:

```
before: total=100.8 rows=3
after:  total=102.0 rows=3 removed=[]
plain control: total=102.0
```

The film is kept, the span total lands on the floor, the slice is bitwise
0.2, the bulk takes 1.8, and the plain (no-interface) control is unmoved at
102.0 — the old behaviour is bitwise the `slice_rows == 0` case of the new
code, as the design claimed. The defect I measured in the status review is
gone.

The edge argument holds as written: `total < min_nm` implies
`min_nm − slice_total > bulk >= 0`, so the assigned bulk is positive by
construction; `bulk_start` is reconstructed as `p + slice_rows`, which
keeps `is_singleton_bulk`'s own `debug_assert` satisfied.

V3 and V4 were both honoured rather than quietly dropped: the bound twin
asserts `bulk == clamp_min_nm` (not `total == clamp_min_nm`), and the
conservative-bound sentence is in `ThinLayerPolicy`'s doc comment where
V3 asked for it. No `debug_assert!(recipe.is_none())` appears next to the
branch, which is what V4 existed to prevent.

### 2.3 C3 — the sweep

All five literals fixed by wrapping the **source** rather than the string,
plus the two cp1252-unencodable characters V2 identified
(`≥` → `>=` at `config.rs:191`, `λ` → `wl=` at `inflate.rs:44`). The two
existing refusal twins are upgraded from substring asserts to exact
full-message asserts, and `test_roundtrip.py` gains the Python-side twin
pinning V1's `MixRule` wording.

**The scanner is a real tool, not a stub.** I exercised
`tools/check_message_whitespace.py` against a synthetic file carrying one
of each class:

| Injected | Verdict |
|---|---|
| 22-space run inside a literal | caught (`space-run`) |
| quoted string inside a trailing `//` comment | correctly ignored |
| literal wrapped with a proper `\` continuation | correctly ignored |
| `≥` in a production literal | caught as `non-ascii-nonencodable` (blocking) |
| `—` in a production literal | caught as `non-ascii` (advisory) |
| space run inside `#[cfg(test)]` | found, then filtered by the test-region rule |

The staging V2 asked for is implemented correctly, including the property
that matters: a **non-encodable** character stays blocking even in an
allowlisted file (the allowlist is consulted only for `kind ==
"non-ascii"`). The rest of C3 — the `Param` comment rewritten to present
tense, `sub_layer_count`'s gradient exception documented, ground rule 7's
stale paragraph retired, the six §8 decisions stamped, `uv.lock` tracked —
all landed as written.

---

## 3. The battery, re-run here

Every gate in ground rule 3, run against `823b7ab`:

| Gate | Result |
|---|---|
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo test --workspace` | **546 passed**, 0 failed |
| `--features opt-minpack-lm` | 551 passed |
| `--features opt-argmin` | 555 passed |
| both features | 560 passed |
| `pytest` (validation) | **764 passed, 1 skipped** (pre-existing feature gate) |
| both fingerprints | hold — random-stack and needle pin green |
| `tools/check_*.py` (now **five**) | all OK |
| `validation/review/*.py` (ten) | all OK |
| Seven version sites | all `0.6.44`; README `:139`/`:170` untouched (N6 respected) |
| `uv.lock` | tracked; working tree clean |
| Push state | `origin/dev_feature` == `HEAD`, 0 ahead / 0 behind |

The counts in C2's and C3's DONE markers (546 / 560 / 764) reproduce
exactly.

---

## 4. Findings

### E1 — the F1.4 strike swallowed the F1.5 row (medium)

`implementation_plan.md:159` is **two table rows joined into one physical
line**. The newline between F1.4's row and F1.5's was lost when the strike
was applied:

```
| ~~F1.4~~ **DONE (0.6.41)** | … | §D5 + correction §1.2 || ~~F1.5~~ **DONE (0.6.42)** | … | §D5 |
```

The master table has 7 columns; this line carries **15 cells**. Markdown
renders the first 7 and discards the rest, so **F1.5's row no longer
appears in the master table at all** — the ladder now shows 14 of 15
items, and the `||` leaves F1.4's row with a stray empty cell.

Why it matters beyond cosmetics: the item that introduced this was
C3's bookkeeping, whose stated job was fixing ground rule 1's *only*
outstanding miss. It fixed one and created another, in the same table, in
the same edit. The next audit that greps for unstruck DONE rows will find
none — because the row it should check is no longer a row.

**Proposal.** Split at the `||`:

```markdown
| ~~F1.4~~ **DONE (0.6.41)** | Schema v2 + a readable-version **range**, not a point | 0.6.41 | **P0** | M (every state file reads through this gate) | M | §D5 + correction §1.2 |
| ~~F1.5~~ **DONE (0.6.42)** | `design_config` rows + Python `Layer.gradient` surface | 0.6.42 | P1 | S | M | §D5 |
```

Worth adding to the ritual, since this is the second table-edit slip in the
series: a cell-count check over the plan's tables is three lines of Python
and would have caught it. Candidate for the C3 scanner's tool directory at
F3.1 rather than now.

### E2 — the scanner's test-region cutoff hides 63–75% of two production files (medium-high)

`scan_file` treats **everything after the first `#[cfg(test)]` line in a
file** as test code:

```python
for i, line in enumerate(lines, 1):
    if "#[cfg(test)]" in line:
        test_from = i
        break
```

That attribute is not only on the trailing test module. It also sits on
**test-only helpers inside production `impl` blocks**, and on out-of-line
test module declarations. Measured over the tree:

| File | First `#[cfg(test)]` | What it is | Real test `mod` | Production lines skipped |
|---|---|---|---|---|
| `smatrix/solver.rs` | 831 | `fn none()` inside a macro | 2561 | 2499 — **75% of the file** |
| `synthesis/structure.rs` | 1185 | `set_row_optimize_for_test` | 1550 | 2031 — **63% of the file** |
| `structure/providers.rs` | 321 | test helper | 412 | 146 — 31% |
| `color/mod.rs` | 83 | `#[path] mod parity;` (another file) | — | tail of the file |

The tool's own docstring states the assumption and asks for it to be
enforced — *"This tree keeps every test module at the file end under that
attribute; assert differently if that ever stops being true"* — but it is
**already not true**, in four files, and nothing asserts it.

The blind spot is **latent, not live**: I re-scanned the hidden windows and
there are no findings in them today. But this is a new *blocking* gate
whose real coverage is far smaller than its documentation claims, and the
worst-affected file is `structure.rs` — the file that produced the original
finding, and the file this series edits most. A message added below line
1185 there is invisible to the gate that exists to catch exactly that.
`solver.rs:181`'s `×` was caught only because it happens to sit at line 181,
above the cutoff.

**Proposal.** Cut at the **last** column-0 `#[cfg(test)]` that introduces an
**inline** `mod … {` — the only form that runs to EOF. An out-of-line
`mod parity;` must *not* be treated as a cutoff, because production code
follows it:

```python
def _test_region_start(lines: list[str]) -> int | None:
    """1-based line of the trailing `#[cfg(test)] mod ... {` block, or None.

    NOT the first `#[cfg(test)]` in the file: that attribute also sits on
    test-only helpers inside production impl blocks (structure.rs:1185,
    solver.rs:831, providers.rs:321) and on out-of-line test modules
    (color/mod.rs:83, `mod parity;`). Cutting at the first one hid 63-75%
    of those files from the gate. Only an INLINE `mod ... {` at column 0
    runs to EOF and is safe to cut at; a file with no such block is
    scanned whole.
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
```

**Verified before proposing:** under this rule the scanned region widens by
**+1730 lines** (`solver.rs`), **+365** (`structure.rs`) and **+91**
(`providers.rs`), `color/mod.rs` becomes scanned whole — and **no new
findings appear, so the gate stays green**. It is a pure coverage gain with
no message work attached.

### E3 — the ceiling refusal still keys on a row count, unmarked (low-medium)

C2 swapped the clamp-up predicate to `is_singleton_bulk()` at
`structure.rs:1118`. Sixty lines above it, the pre-mutation ceiling refusal
still reads:

```rust
if sp.end - sp.start <= 1 || self.span_is_scalable(sp) {
    continue;
}
```

`structure.rs:1056` — a **row count**, which is precisely the construct the
status review flagged as wrong at the other site. The divergence *is*
deliberate: amendment 3 §3 scopes the ceiling side out explicitly, and
CHANGELOG `[0.6.43]` records it as a known asymmetry. But neither record is
in the code. A reader who notices the two predicates disagree has nothing
at the call site telling them the disagreement was chosen, and the two most
likely reactions — "fix" it, or re-file C2 — are both wrong.

**Proposal.** One comment inside the refusal loop, above the `if` at
`structure.rs:1056` — no behaviour change:

```rust
// The row count here is DELIBERATE and diverges from the clamp-up
// branch's `is_singleton_bulk()` below (C2): the floor side clamps a
// plain interface film, the CEILING side keeps refusing it. Shrinking a
// plain film to the ceiling would be well-defined, but changing F0.2's
// refusal semantics is a separate decision - amendment 3 section 3
// records the asymmetry. Revisit only with a real user case.
```

### E4 — the advisory allowlist's accounting disagrees with itself (low)

Four statements about the same list, three of them wrong:

| Source | Says | Actual |
|---|---|---|
| `NON_ASCII_ALLOWLIST` | — | **7 entries** |
| The set's own explanatory comment | names 5 files | omits `optimizer.rs`, `trf.rs` |
| CHANGELOG `[0.6.44]` | "five production files" | 7 |
| `implementation_plan.md` (C3 entry) | "ten entries" | 10 is the *finding* count, not the entry count |

This matters more than a typo because the list is a scheduled deliverable:
F3.1's exposure re-audit is supposed to retire it entry by entry, and an
audit that starts from a list which misstates its own contents will leave
two files behind.

Second, smaller point on the same structure: the allowlist is **file-granular**
(`"any literal in it"`), so once a file is on the list a *newly introduced*
encodable non-ASCII character anywhere in it passes unremarked. That is a
deliberate simplification and acceptable while the list is shrinking toward
empty — but it is worth one sentence in the tool's comment, because
"allowlisted" currently reads as "these known sites", not "this whole file".

**Proposal.** Correct the comment to enumerate all seven, fix the CHANGELOG
line to "seven production files (ten literals)", correct the plan's "ten
entries" to "ten literals across seven files", and add the file-granularity
sentence. Alternatively — and better if F3.1 is near — make the allowlist
`(file, line)` pairs so the retirement is a checklist rather than a
judgement call.

### E5 — dead locals in `scan_file` (trivial)

Three assigned-never-read locals in the new tool:

- `test_from` (lines 76–79) — computed in `scan_file`, never used there;
  `main` recomputes it at 189–195. Misleading, since it looks like the
  function filters test regions itself (it does not).
- `in_block_comment` (line 84) — never read; block comments are handled by
  the `raw.find("*/")` skip instead.
- `col` (lines 103, 123) — never read.

**Proposal.** Delete all three. Removing `test_from` from `scan_file` also
removes the false impression that the cutoff lives in two places — relevant
once E2 moves it into a helper.

### E6 — amendment §4.4 undercounts the open decisions (trivial)

The bookkeeping item says *"Decisions 1–5 and 11 stay open — they are
Phase B's to resolve at their items"*. Decisions **6** (per-layer EMA host
selection) and **7** (`ProfileShape` beyond `Linear`) are also unstamped.
Nothing is lost — both are self-documented in §8 as deferred until a real
case, which is a different status from "Phase B's to resolve" — but the
sentence as written implies §8 has no other open items.

**Proposal.** Amend to "Decisions 1–5 and 11 stay open for Phase B;
6 and 7 remain deferred-until-a-real-case, as §8 already records."

---

## 5. What the implementation got right

Worth recording, because these are the things that usually slip:

- **The probe was actually run, and its id written down.** C1's outcome is
  checkable by anyone later, which is the whole difference between
  "Outcome A" as a finding and as an assumption.
- **The V-corrections were honoured, not absorbed.** V3's wording change
  reached the assertion (`bulk == clamp_min_nm`, not `total ==`), V4's
  warning reached the absence of a `debug_assert`, V1's third site was
  fixed and given the Python-side twin it asked for, V7's corrected
  rationale is the one recorded in the commit message.
- **V2's measurement was re-measured, and the amendment corrected back.**
  C3's DONE marker says the ledger was wrong in two directions and names
  both — two files misplaced, two more missed. An item that corrects its
  own specification while executing it is the ritual working.
- **The snippet was improved rather than pasted.** The `end`/`p`
  generalisation in C2 (§2.2) is a bug the amendment would have shipped.
- **The scanner's own output is ASCII-safe**, and the commit says why: its
  first draft crashed printing findings on a cp1252 console — the exact
  class it polices. That is the kind of detail that normally goes
  unrecorded.
- **Per-item pushes, each with a green run.** The failure mode that
  produced C1 in the first place — a whole series pushed as one batch, so
  CI never evaluated an individual release — has been retired in practice,
  not just in the rule.

---

## 6. Recommended order of work

1. **E1 + E2 together, one bookkeeping commit, no version bump** (docs +
   tooling only — the R5/C1 precedent). These are the two that degrade the
   *next* audit: one hides a ladder row, the other hides two thirds of two
   files from a blocking gate. E2's fix is verified green, so the commit
   cannot turn CI red.
2. **E3, E4, E5 in the same commit** if convenient — all three are comments
   and deletions, none touches behaviour.
3. **E6** wherever amendment 3 is next edited.
4. **Then F2.1 (0.6.45).** Nothing in this review blocks it. The Phase B
   re-versioning (0.6.45–0.6.49) is already stamped in the plan's ladder
   and in amendment 3 §5.

One item to carry forward into F3.1's scope, beyond what §6 of the
amendment already lists: **the allowlist retirement needs E4's corrected
list**, and is easier if E4 takes the `(file, line)` option.

---

## 7. Disposition summary

| ID | Finding | Severity | Fix |
|---|---|---|---|
| E1 | F1.4 strike joined two master-table rows; F1.5's row no longer renders | medium | split line 159 |
| E2 | Scanner cuts at the *first* `#[cfg(test)]`; 63–75% of two production files unscanned | medium-high | `_test_region_start` helper (verified: +2186 lines scanned, gate stays green) |
| E3 | Ceiling refusal keeps a row-count predicate 60 lines from the one C2 replaced, unmarked | low-medium | one comment inside the refusal loop at `structure.rs:1056` |
| E4 | Allowlist size stated three different ways (7 / "five" / "ten"); file-granularity undocumented | low | correct all three; consider `(file, line)` entries |
| E5 | `test_from`, `in_block_comment`, `col` assigned and never read | trivial | delete |
| E6 | §4.4 undercounts §8's open decisions (6 and 7 omitted) | trivial | one sentence |

No finding requires a behaviour change, a version bump, or a fingerprint
re-record.
