# Navette Feature Plan — Implementation Status

Review of `docs/implementation_plan.md` (revision 4) against the tree at
`31b8f6d` (`dev_feature`, 0.6.42), 2026-09-13.

Method: read the plan end to end, checked every **DONE** item's design
claims against the code that implements it, ran the full verification
battery from ground rule 3 locally, and read the CI history for the
branch. One behaviour claim was settled by a throwaway Rust probe
(built, run, deleted — the tree is unchanged apart from this file).

---

## 1. Verdict

**Phase 0 and Phase A are complete and the code matches the plan.** Ten of
fifteen items have landed (`0.6.33 → 0.6.42`), each with its CHANGELOG
entry, its CORRECTIONS block, and the named gate twins actually present in
the tree. The plan's ritual has been followed closely enough that spot
checks kept confirming it rather than catching it out — including the
easy-to-miss ones (N6's README line-not-file rule; R5's fixture committed
one commit *before* the gate it tests).

**But the branch is red.** The needle-run fingerprint — the pin F0.1
created precisely so that "fingerprint unmoved" would name a real test —
**fails on Linux in CI** while passing on Windows and locally. That is the
one finding that blocks the series, and it is the failure mode ground rule
4 was written about: a local battery that is honestly green while CI is
not.

Phase B (multi-environment, F2.1–F2.4) and Phase C (F3.1, docs/release)
have not been started; no environment code exists in the tree, which is
consistent with the progress log.

---

## 2. Where the plan stands

| ID | Version | Plan status | Verified in tree |
|---|---|---|---|
| F0.1 | 0.6.33 | done | `DesignStack::spans`, `Span::is_singleton_bulk`, `bulk_start`, `assert_spans_partition`, `needle_host_refusal`, `emit_entry` extraction, needle pin — all present; the four named twins exist (`n1_interface_film_bulk_is_still_a_host`, `is_singleton_bulk_four_cases`, `needle_refusal_at_the_mutator_names_span_and_range`, `merge_bookkeeping_intra_and_cross_span`) |
| F0.2 | 0.6.34 | done | All four rules are span quantities: floor and cap read the span's slice-inclusive total (`clamp_all_policy`), the budget counts `spans().len()`, inflate filters on `is_singleton_bulk()`. `ClampReport` + `final_clamp_report` present; refusal checked before mutation |
| F0.3 | 0.6.35 | done | `ThinLayerPolicy::{Remove, ClampUpFinal, ClampUpAlways}`, parsed, defaulted to `Remove`, wired to the LM lower bound; twins present. **One gap — see finding 2** |
| F1.1 | 0.6.36 | done | `GradientSpec`/`GradientMode::FixedSpan`/`MixRule` reuse, `gradient_sub_layer_count`, homogenize path, `_norm_gradient` on the Python design door |
| F1.2 | 0.6.37 | done | `GradientMode::RateCapped` with the two-sided semantics the CORRECTIONS block records |
| F1.3 | 0.6.38 | done | `InhMode::RateCapped`, `delta_layer()` as the single source for all six readers incl. the PyO3 getter (B6) |
| F1.6 | 0.6.39 | done | `Param::{Row, Span}` in `build_params`/`apply_params`; `span_is_scalable_rows`; scalable spans capped as a bound instead of refused |
| F1.7 | 0.6.40 | done | `SpanRecipe` + `emitted_total` staleness marker; `refresh_profiles` called at **exactly three** production sites (`from_design`, macro-cycle top before the budget check, after the final clamp) and nowhere else — U4 holds as written |
| F1.4 | 0.6.41 | done | `SCHEMA_VERSION = 2`, `MIN_READABLE_SCHEMA_VERSION = 1`, range gate in `version.rs`, shared `_state_version_ok` in `config/models.py`, v1 fixture + oracle committed at `ab7e27b` (before the gate changed) |
| F1.5 | 0.6.42 | done | `LayerRow.gradient`, `GradientSpec` with `deny_unknown_fields` on the nested spec, Python `Layer(gradient=…)` + getter/setter, `.pyi` in sync |
| F2.1–F2.4 | 0.6.43–46 | not started | confirmed: no `environments` / `EnvironmentSegment` / `residuals_multi` anywhere in `rust/` or `src/`; no `bench_eval.py` |
| F3.1 | 0.6.47 | not started | confirmed: no user docs or examples for any Phase A feature (see finding 8) |

Version sites: all seven carry `0.6.42` and agree
(`pyproject.toml:14`, `Cargo.toml:9,29`, `__about__.py:12`,
`Cargo.lock:465,492`, `README.md:190`). N6 honored — `README.md:139` and
`:170` still say 0.6.32 as historical facts.

Size of the series: 50 files, +12 680 / −385 since `72a2d4d`
(rust +7 252, validation +1 091, python +161, docs/CHANGELOG +4 157).
`structure.rs` alone went 853 → 3 060 lines.

---

## 3. Verification battery — local vs CI

### Local (this machine, release build, `build_profile() == "release"`)

| Gate | Result |
|---|---|
| `pytest validation -q` | **763 passed, 1 skipped** (the skip needs the `opt-minpack-lm` feature) |
| `cargo test --workspace` | **543 + 22 + 15 doc-tests passed**, 0 failed |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| `check_exposure` / `check_cie_sync` / `check_pyi_sync` / `check_toolchain` | all OK |
| the ten `validation/review/*.py` harnesses | all exit 0 |
| `--features opt-minpack-lm` / `opt-argmin` / both | 548 / 552 / 557 passed |

Every gate the plan names is green here.

### CI — red

[Run 34749211804](https://github.com/opticsWolf/Navette/actions/runs/34749211804)
on `31b8f6d`: **failure**.

- `python (ubuntu-latest)` — **FAILED** `test_needle_pin.py::test_needle_run_is_deterministic_and_matches_the_recorded_digest`
- `dependency floors` (ubuntu) — **FAILED**, same test
- `python (windows-latest)` — passed
- `rust (test + warnings)` — passed

```
needle fingerprint moved:
  recorded cf753a910c7bf1d2e6c5ea5b9f5c696601ab675abaf497db28bac5340a63044d
  got      5c30981320d71b9228a2dbc3b1942599a59a0bcdecd1c921c8b49cc8dbf547da
```

CI history for the branch: `072c5cd` (the last docs-only commit, still a
0.6.32 tree) was green; the next run is the tip. **All twelve commits of
the series were pushed as one batch, so CI has never evaluated any
individual release of the ladder** — including `8a0ea5c`, the commit that
recorded the pin.

---

## 4. Findings

Ordered by what they cost. Everything below was verified against the tree
or measured, not inferred from the plan.

### 1. The needle pin is not portable, and nothing in the plan says it has to be — blocking

`validation/regression/synthesis/test_needle_pin.py:122` hard-compares a
SHA-256 of hex-rendered floats against a constant recorded on a **Windows
release build** at 0.6.32. On Linux the run is deterministic (the test's
own `d1 == d2` assertion passes) but hashes to a different value, and both
Linux jobs fail.

Two things make this worth treating as the headline item:

- It is the repo's **only** cross-machine recorded float constant. The
  other two fingerprints are not exposed to this: the random-stack pin
  (`test_differential.py:95`) compares the Rust engine against the Python
  engine *in the same process*, and `test_state_fingerprint` hashes a key
  set, not floats. The pin R1 asked for introduced a portability
  requirement that no prior gate had.
- **Whether the Linux digest also differed at 0.6.32 is currently
  unknown**, because the pin has never run on Linux CI. So the failure is
  one of two very different things: a platform float divergence (libm
  `exp`/`cos` differences feeding an LM trajectory — a pin-design defect,
  harmless to the features) or a genuine behaviour change that Windows
  happens not to show (a real regression somewhere in F0.1–F1.5).

**Settle it before anything else, and it is cheap:** push `8a0ea5c` (the
pin at a 0.6.32 tree) to a throwaway branch and read the Linux job. Green
there and red at the tip means a regression; red in both means the pin
needs a platform-keyed digest or a tolerance-based oracle, and the plan's
ground rule 5 needs a sentence about what "the fingerprint" means across
platforms. (`ci.yml` has no `workflow_dispatch`, so a branch push is the
only trigger.)

Until this is resolved, no later item can honestly claim "fingerprint
unmoved" — the gate itself is in an unknown state on half the CI matrix.

### 2. Under a clamp-up policy, an interface-carrying film is deleted where the same film without an interface is kept — measured

`clamp_all_policy` selects the clamp-up branch with `sp.end - sp.start == 1`
([structure.rs:1107](rust/navette/src/smatrix/synthesis/structure.rs:1107))
— a **row count** — while the repo's own predicate for "one physical
layer" is `is_singleton_bulk()`, which N1 introduced precisely because an
interface-carrying plain film is *two* rows (slice + bulk) and still one
layer.

Measured on a throwaway probe, identical 0.8 nm TiO₂ film under a 2.0 nm
floor with `clamp_up = true`:

| Film | Rows | Result |
|---|---|---|
| plain | 1 | clamped up to 2.0 nm, kept |
| `interface = true`, `interface_thickness = 0.2` | 2 | **span removed whole**, total 100.8 → 100.0 nm |

The deletion *is* reported (`spans_removed: ["TiO2 (0.8 nm)"]`), so it is
not silent — but it is the opposite of what the user asked for by
selecting the policy whose whole purpose (U1) is "set that layer to the
minimum you can actually deposit instead of deleting it."

The doc comment does describe the mechanism honestly ("sets a surviving
ONE-ROW span"), and F0.3's own deferral sentence justified it — but that
justification was written for *graded* spans, on the grounds that clamping
a span up is F1.6's scale operation, and **F1.6 has since landed**. For a
2-row interface span the operation is not even a scale: it is
`bulk = min_nm − t_slice`, well-defined, and it is the case B1 already
established the slice belongs to the authored thickness for.

No twin covers it: `f03_thin_graded_span_is_removed_whole_under_clamp_up_too`
uses a 4-row graded span, and the F1.6 clamp twins all use scalable spans.

Suggested home: F3.1 is too late for a behaviour fix; this wants either a
tail item on the Phase A ladder or an explicit "known and deliberate"
sentence in the `ThinLayerPolicy` docs plus a twin that pins the choice.

### 3. Two user-facing refusal messages carry a 22-space run — introduced at F0.2

[structure.rs:1057](rust/navette/src/smatrix/synthesis/structure.rs:1057)
and [pipeline.rs:104](rust/navette/src/smatrix/synthesis/pipeline.rs:104),
both from `e7a0e61`:

```
clamp_all: span 'TiO2' is 1000.0 nm thick, above the 300.0 nm ceiling -                      refusing rather than rescaling a profile
```

A `\` line continuation was lost when the literal was wrapped. Both are
ASCII (ground rule 7 holds) and both name both sides (ground rule 6
holds), so only the rendering is wrong — but these are the two messages a
user meets when the licence item 2 refusal fires. The existing twin
(`f02_new_refuses_pinned_span_above_the_ceiling`) asserts substrings
(`'TiO2'`, `1000.0`, `300.0`), which is why it passes. Two pre-existing
messages have the same shape (`interpolate/mod.rs:210`,
`thick_opt.rs:323`, the latter also carrying a non-ASCII `×`); those are
not this series' doing and the plan's "do not fix as drive-by" rule
applies to them.

### 4. A doc comment that F1.7 should have retired

[evaluator.rs:425](rust/navette/src/smatrix/synthesis/evaluator.rs:425),
written at F1.6 and untouched by F1.7:

> the RateCapped modes depend on absolute depth and stay homogenized
> until F1.7.

F1.7 shipped. `span_is_scalable_rows` reads optimize flags only, so rate
spans *are* scalable now and `refresh_profiles` is what makes that sound.
The comment sits on the `Param` enum — the first thing a reader of the
parameter machinery meets.

### 5. Plan bookkeeping drift — two spots

- **F1.4's master-table row is not struck.**
  [implementation_plan.md:144](docs/implementation_plan.md:144) still reads
  `| F1.4 | Schema v2 … |` while §2002 says **DONE (0.6.41)** and §9 lists
  it done with its commit. Ground rule 1 requires the row struck. Every
  other done item's row carries `~~Fx.y~~ **DONE (…)**`.
- **§8 decisions 8–14 shipped without being stamped.** Only decision 15
  carries a `RESOLVED` block. Decision 12 in particular ("how a span
  declares itself scalable", marked *decide before F1.4, because F1.4
  freezes the Phase A field set") was in fact decided as recommended —
  `GradientSpec` has no `scalable` field and `span_is_scalable_rows` reads
  `optimize` — but the section still reads as open, and F1.4 has already
  frozen the schema around that choice. The Phase B decisions (1–5, 11)
  are genuinely still open.

### 6. `Layer.sub_layer_count` reports 1 for a gradient carrier

`Layer::sub_layer_count()` ([layer.rs:144](rust/navette/src/structure/layer.rs:144))
is gated on `self.inhomogen`, so a gradient-only layer returns 1 while it
expands to ≥3 rows. The getter's doc calls it the "solver sub-layer
count", and it is public API (`navette.structure.Layer.sub_layer_count`).

This is a design constraint rather than an oversight — the gradient rule
(`gradient_sub_layer_count`) needs the wavelength grid and both endpoint
spectra, which a bare `Layer` does not have — but the getter says nothing
about it, and B6 established this reader as one that must not diverge from
the emitted count. A doc note naming the gradient exception (or returning
an honest sentinel) belongs with F3.1's exposure re-audit.

### 7. `uv.lock` is neither tracked nor ignored

It shows up as untracked in every `git status` on the branch. The repo's
stated stance (`.gitignore:4`, `:23`) is that lockfiles pinning
reproducible builds stay tracked. Track it or ignore it; the third state
is noise in a workflow whose ritual reads `git status` fifteen times.

### 8. Phase A is shipped and undocumented — expected, but worth stating

No README section, no example, and no line in `docs/materials-*.md`
covers gradients, `ThinLayerPolicy`, scalable spans or profile refresh.
That is F3.1's scope by design, but it means five releases of user-visible
features currently have their only documentation in CHANGELOG entries and
in the plan. If anything ships from this branch before F3.1, that gap
ships with it.

---

## 5. What the implementation got right

Worth recording, because a findings list reads as if nothing did.

- **The gates named in the plan exist as tests with the names the plan
  uses.** `assert_spans_partition`, `is_singleton_bulk_four_cases`, the N1
  twin, the needle-refusal twin, the merge twins, the F1.7 refresh-count
  twin, the F1.4 accepted-range probe. This is unusual and it is what
  made this review cheap.
- **U4 holds exactly.** `refresh_profiles` has three production call sites
  and they are the three the plan lists. The staleness marker
  (`emitted_total`, the *measured* sum at last emission) is a better
  answer than the plan's own "re-derive at current D", and the CORRECTIONS
  block explains why rather than quietly substituting it.
- **The corrections are honest about where reality bit.** The one-ulp knee
  at F1.2, the unreachable scale-invariance numbers at F1.6, the
  "`emit_entry` gained two parameters" admission against §8.15's own
  contract at F1.1 — each names the plan sentence it invalidates.
- **The bit-exactness discipline is real where it can be checked.** F0.2's
  pin diagnosis (strip `clamp_report`, and the 0.6.32 digest holds) is the
  strongest single piece of evidence in the series, and it is exactly the
  kind of argument the plan asked for.

---

## 6. Recommended order of work

1. **Settle finding 1.** Push `8a0ea5c` to a throwaway branch, read the
   Linux job, and then either fix the regression or re-specify the pin.
   Nothing else on the ladder means much until the fingerprint gate is
   trustworthy on the whole CI matrix.
2. **Decide finding 2** — fix the interface clamp-up, or pin the current
   behaviour as deliberate with a twin and a doc sentence.
3. Sweep findings 3–7 (one commit; none of them touch behaviour except 2).
4. Strike F1.4's row, stamp §8 decisions 8–14, and record decision 12's
   resolution where F1.4 froze it.
5. Then F2.1, with §8 decisions 1–5 and 11 resolved at the item as the
   plan prescribes — and consider pushing per item from here on, so CI
   sees each release rather than the end state.
