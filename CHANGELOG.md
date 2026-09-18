# Changelog

All notable changes to Navette are recorded here. Work items reference
`docs/remediation_plan.md` (Rx.y), `docs/code_review.md` (§),
`docs/implementation_plan.md` (Fx.y), and `docs/implementation_plan_pd.md`
(PD1–PD4).

## [0.7.12] - C3/C5/C6/C10: what an incoherent flag can and cannot mean

### Added
- **C3.** A flagged layer thinner than five wavelengths of optical
  thickness now warns at construction, at both engine doors. The flag
  asserts the layer destroys the phase relation between its two surfaces,
  which needs the path-length spread across it to exceed the source
  coherence length `L_C = lambda^2 / delta_lambda`; nothing checked that.
  The threshold is **derived, not chosen**: the warning quotes the source
  bandwidth that layer would need — `delta_lambda > lambda^2 / (2*n*d)`,
  and what percentage of the wavelength that is — so "too thin" is a
  number the caller can argue with rather than a constant somebody picked.
- **C5.** A flag on row 0 or the last row now warns that it does nothing.
  The block sweep scans `current_idx + 1` up to `idx_n` and applies the
  attenuation element only below `idx_n`, so half-space flags are never
  consulted: bit-identical output, previously with no diagnostic at all.
  The warning says what to do instead — a thick substrate is an *interior*
  layer with a real thickness, flagged, between the film stack and the
  exit medium.
- `optics_core::THIN_FLAGGED_LAYER_EXPLANATION` and
  `THIN_FLAG_WAVELENGTHS`, shared by both doors in the
  `AMBIENT_DROP_EXPLANATION` shape, with
  `test_both_doors_explain_thin_flags_the_same_way` to keep them from
  drifting.
- `validation/smoke/test_incoherent_flag_sanity.py` (14 tests), including
  the measurement behind each warning: a half-space flag really is
  bit-identical, an interior one really is not, and the flag is still
  honoured after the thin-layer warning.

### Changed
- The `incoherent_flags` constructor docstring no longer offers "thick
  substrate" as its example. That was the exact row a caller would flag to
  no effect (C5). It now states the interior-only rule, the `L_C` criterion
  with the 50–100 µm practical substrate threshold, and that the partition
  changes the answer by itself — a **zero**-thickness flagged layer still
  decoheres, because the join breaks the p-s phase relation whatever the
  thickness (`Rs 0.193432 -> 0.109790`).
- **C6.** `field_profile` gained the front-block/single-block caveat,
  completing the sweep begun in 0.7.10 (`ellipsometry`, `stokes`,
  `complex_amplitudes`).
- **C10.** `eigenmode_landscape`, `find_eigenmodes`, `refine_mode` and
  `field_profile` now say that the guided-mode kernels solve `[0, n-1]` as
  one coherent block and never consult `incoherent_flags`, and
  `smatrix/optimizer.rs` says it at the module level. This is correct by
  intent — a guided mode is a coherent-stack concept and there is no
  eigenmode to find across a partition — but every other surface on the
  same `Solver` honours the flags, so a caller had no way to tell these do
  not.

### Unchanged on purpose
- Both are warnings, not refusals. A half-space has no second surface to
  lose coherence against, so a flag there is a no-op rather than a
  mistake; and the thin incoherent limit is a legitimate thing to model,
  as long as it is what the caller meant. The engine honours the array
  either way.
- All three coherence findings now share one interior-only predicate, so
  C2's refusal, C3's thickness warning and C5's no-op warning cannot
  contradict each other. A stack flagged only on its half-spaces gets C5's
  warning and neither of the others — pinned by
  `test_the_two_findings_do_not_contradict_each_other`.
- Three existing tests now emit the new warnings and still pass
  (`test_needle_t_a_phi`, `test_the_multiblock_span_is_not_the_coherent_one`,
  `test_half_space_flags_alone_are_allowed`). All three use thin or
  half-space flags deliberately, as FD and partition geometry rather than
  as physics; the warning saying so is the warning working.

## [0.7.11] - C4: optical gain is refused at every door

### Fixed
- `Im(n) < 0` is now refused instead of silently mangled (physics review
  rounds 1 and 2, C4). Gain has no representation in either solver path,
  and the two paths destroyed it *differently*: the coherent kernels
  conjugate the propagation phase back to decay, so a gain layer came back
  wearing the loss layer's answer — bit-identical to `+k` on a symmetric
  stack, reporting a **positive** absorptance for a medium that
  amplifies — while the incoherent cascade clamps the same quantity to
  zero, turning the layer transparent and opening the energy books by
  `A = -4.3e-5`. Neither is any physical system, and neither is
  recoverable from the output: it reads as an ordinary absorbing stack.
- The check sits at `Solver::assemble`, the one place every constructor
  funnels through, so `ScatterMatrix`, `solve_arrays` / `solve_structure`
  and the raw `core_engine` FFI are covered by one rule rather than three
  copies. This is deliberately **unlike** C2, whose refusal had to stay at
  the doors because Mode A's numbers are a legacy port reached from below;
  nothing below the doors computes anything about gain worth preserving.
- The free `needle_gradient` takes a flat cache and never builds a
  `Solver`, so it carries the check itself — the same gap C8's flag
  canonicalization had to close separately. It also checks the needle
  material, which arrives as its own argument: the needle's index goes
  into the same conjugating kernels and its spacer `tau` takes the same
  clamp.
- The `Structure` door, which was the one door already refusing `k < 0`,
  now carries the same explanation as the other two. Its message said
  "check provider data" and nothing about what the solver would otherwise
  have done with the number.

### Added
- `optics_core::GAIN_MEDIUM_EXPLANATION`, the shared text behind all three
  refusals, in the `AMBIENT_DROP_EXPLANATION` shape, plus
  `gain_medium_message` and `scan_for_gain` so every Rust site reports the
  same grid position the same way. `test_both_doors_explain_gain_the_same_way`
  fails if any one door is edited without the others.
- `validation/smoke/test_gain_refusal.py` (9 tests): each door refuses and
  names the layer, the wavelength index and how many grid values are
  negative; the explanation states what each path does, the measured size
  of the damage and the way out for a time-convention mismatch; and the
  controls pin what must **not** refuse.

### Unchanged on purpose
- `-0.0` is not gain. `-0.0 < 0.0` is false in IEEE, `forward_branch` and
  `sanitize_incident_index` already treat a signed zero as zero, and a
  provider that writes `-0.0` for a transparent material is not describing
  an amplifier. Pinned on both sides, because the obvious "tidy up the
  sign" edit would start refusing real grids.
- The gate is on the **index** array only. The flat-array roughness
  surface stays permissive exactly as `test_layer_gate.py` pins it; the
  dividing line is whether a correction exists, not which surface the
  value arrived on. Indices were already gated there for non-finite values
  and for `|n|^2` overflow, and gain is the same kind of problem.
- The ambient is checked like any other row, and the gain gate runs before
  the absorbing-ambient drop. That drop counts a negative imaginary part
  as absorption, so left alone it would have swallowed a sign error and
  reported dropping absorption that was never there.

## [0.7.10] - C2/C9: the front-block cross channel is refused, not mixed

### Fixed
- Cross-channel observables under the default `FRONT_BLOCK` mode on a stack
  with an interior incoherent layer now refuse instead of returning a
  vector built from two different stacks (physics review rounds 1 and 2,
  C2). Mode A takes the p-s cross channel from the FIRST coherent block
  while its intensities are totals over every incoherent echo, so `DOP_R`
  came back as `|rs_c|^2/Rs` — measured 0.763904485 where the physical
  answer for a non-depolarizing stack at normal incidence is exactly 1,
  and `Delta` up to 22.2 degrees out at 70 degrees.
- The refusal covers all twelve `NEEDS_CROSS` bits (Delta, DOP, S2/S3,
  retardance, and the raw `cross_R`/`cross_T`) as one rule. There is no
  "raw cross is fine" split: the raw channel is the same defective object,
  and letting it through would let a caller rebuild the broken DOP by hand.
- The predicate is deliberately narrow — `FRONT_BLOCK` **and** a cross bit
  **and** a non-zero flag on an INTERIOR row. Rows 0 and last are
  half-spaces whose flags the engine never consults (C5), so flagging only
  those does not refuse. With nothing flagged, Modes A and B are
  bit-identical, so nothing refuses there either. Mode C is one block over
  the whole stack and its cross channel is correct.
- Both doors carry one shared explanation, `FRONT_BLOCK_CROSS_EXPLANATION`,
  in the `AMBIENT_DROP_EXPLANATION` shape:
  `ScatterMatrix.compute` (Python, where `stacklevel` can point at the
  caller) and `solver::solve_arrays` (Rust, which is what `solve_structure`
  and the raw FFI reach). `test_both_doors_explain_the_cross_channel_the_same_way`
  fails if either is edited without the other.
- `NEEDS_CROSS` is now exported to Python and imported by the door rather
  than re-declared there, so the door and the engine cannot test different
  bits — the `NREQ_*` pattern. `test_needs_cross_is_bound_not_copied` pins
  the composition against the Rust constants and pins Psi's absence from
  the mask (it is an amplitude ratio, not one of the mixed objects).

### Unchanged on purpose
- **The engine is untouched.** `solve_point`, `solve_point_intensity` and
  `Solver::solve` compute exactly what they computed before. The legacy
  parity port drives Mode A with interior flags and cross observables
  through the raw `core_engine` pyfunction, which enters below both doors,
  so it still runs and still pins its numbers bit-for-bit.
  `test_the_engine_still_computes_what_the_doors_refuse` fails first if
  the door check ever leaks inward.
- Mode A's `cross_T` is recorded rather than changed (C9). It is a third
  object — a product of per-block `t_p*conj(t_s)` across the joins with no
  multiple-bounce series — neither the front block nor the Mode B cascade.
  `test_mode_a_and_mode_b_cross_terms_are_different_objects` pins that,
  alongside the A/B photometric identity.
- The Python default stays `FRONT_BLOCK`. It matches the legacy port by a
  remediation-plan decision, and flipping it would silently re-litigate
  that for stored results. With the refusal in place the default is no
  longer dangerous, only loud; a flip would be its own bump.
- The R6.2 DOP clamp is unchanged. Its comment now states the scope
  condition it always had — the reflected identity holds "for a single
  coherent block" — and records that Mode A on a flagged stack produces a
  DEFICIT the `.min(1.0)` clamp cannot see.

## [0.7.9] - C1: synthesis refuses the flag it cannot honor

### Fixed
- A film marked `coherent: false` now refuses at the synthesis door
  instead of being silently ignored (physics review rounds 1 and 2, C1).
  `DesignStack::solver_arrays` materialized the flag faithfully
  (`incoherent_flags[slot] = i32::from(!layer.coherent)`) and
  `simulate_inner` then solved `[0, nl-1)` as ONE coherent block without
  ever reading it; the needle pass built its stack fields over the same
  single block. Every consumer of that funnel — the design door, the
  multi-environment door, the LM step, the FD jacobian, the needle scan —
  returned the fully-coherent answer for a stack the caller had marked
  otherwise. Measured: merit bit-identical with the flag on and off,
  while the engine door moves on the same stack.
- The refusal sits in `DesignStack::from_parts`, the one internal
  constructor every `DesignStack` passes through, so no route into
  synthesis can carry the flag — `with_films`, `from_design`,
  `insert_needle_seed`, `merge_adjacent` and the rest reach it. It names
  the film by index and material, and it names the doors that DO honor
  the flag (`solve_structure`, `ScatterMatrix`) rather than only saying
  no.

### Unchanged on purpose
- Ambient and substrate are not checked. They are half-spaces at rows 0
  and last, where the engine never consults the flag either (C5), so a
  flag there is the same no-op on every door and refusing it would be a
  different finding.
- This is the refusal half of C1 only. `p_function_multiblock` still
  exists, validated and unwired, and wiring the synthesis funnel to the
  multiblock path stays available as a later phase — unblocked by this,
  which only converts a silent wrong answer into a loud one.

## [0.7.8] - C8: one coherence flag, one meaning

### Fixed
- A coherence flag of any non-zero value other than 1 no longer deletes
  the flagged layer's absorption (physics review round 2, C8). The
  published contract is "non-zero where a layer breaks phase coherence"
  and the block sweep agrees — it extends a coherent run only while
  `flag == 0` — but the attenuation element was gated on `flag == 1`.
  A flag of 2 therefore partitioned the stack and then skipped `tau`:
  the layer decohered but never absorbed. Measured on an air / 2 um
  `n = 1.5 + 0.01i` slab / air stack at 550 nm: `flag = 1` gives
  `Ts = 0.583945129` with `A = 0.361243521`, while `flag = 2`, `-1` and
  `7` all gave `Ts = 0.923089546` with `A = -0.000042666` — the exact
  lossless transmittance, and a negative absorptance from the open
  energy books.
- `incoherent_flags` is now canonicalized to 0/1 in `Solver::assemble`,
  the shared tail every constructor funnels through (`new`,
  `from_wav_major_flat`, `from_raw`, and so `solve_arrays`, `ScatterMatrix`
  and `core_engine` with them), and in the free `needle_gradient`, whose
  flags arrive from the caller rather than from `self`. Both gates now
  read the same array, and a future gate is honest whichever comparison
  it picks.

### Unchanged on purpose
- `coherence_mode` needed no work: `Solver::validate` already refuses
  anything outside {0, 1, 2} with a named message, every door reaches it
  through `new`, and `validation_refuses` has pinned it since before this
  review. Round 2 read the raw FFI door as unvalidated; it is not.
- Flags of 0 and 1 are bit-for-bit unchanged, so no fingerprint moves.
  The canonicalization can only alter a stack that was already getting a
  physically impossible answer.

## [0.7.7] - M7: the merit spec remembers its roster

### Fixed
- A pre-built `MeritSpec` whose environment roster is in a different
  order from the design request's now refuses at the run door instead of
  scoring every demand against the wrong surroundings (review PB, M7).
  A demand's environment is a POSITION in the roster — `resolve_env`
  turns each `environment=` tag into an index at compile time and the
  names were then discarded, leaving `n_envs` (a COUNT) as the only
  thing binding a spec to a design. A typo always refused, early, at
  `build_merit_spec`; a permutation has no shape error at either gate,
  and the run completed: measured on the 0.7.5 wheel at
  `8523.676912295265` where the correct pairing gives
  `19876.13664037235`, two environments that differ by two orders of
  magnitude in cover thickness with asymmetric demands.
- `MeritSpec` carries `env_names` with `env_names()` / `set_env_roster()`;
  `compile_merit_spec` records the roster **only when the target set
  named one**, and `set_n_envs` refuses a count that contradicts a
  recorded roster. `run_environments` compares names position by
  position when it has them and refuses with both lists and the way out.

### Unchanged on purpose
- A spec whose target set named no environments records nothing and
  keeps the count check alone. The alternative — recording the synthetic
  `["default"]` roster — would newly refuse every untagged
  single-environment request whose design environment is called anything
  else, which is most of them. Twinned: a matching roster is bit-equal
  to the unrostered spec in Rust and to the `TargetCollection` path in
  Python, and `env_names()` stays empty for a set that named nothing.

### Changed
- `run_needle`'s `environments` docstring and
  `docs/spectralweave-target-kinds.md` described the positional binding
  as a live hazard with "pass the collection, not a pre-built spec";
  both now describe the refusal and say exactly how narrow it is (name
  your environments in both places and a mismatch refuses; name them in
  neither and nothing changes).

## [0.7.6] - M1: a per-film override cannot free a surrounding

### Fixed
- `per_film_flags[<surrounding's material code>] = {"optimize": true}`
  (or `"needle"`) now refuses at compile instead of turning a fixed
  surrounding into a free variable (review PB, M1). Flag application
  runs global map -> row -> per-film override, and the override is keyed
  by MATERIAL CODE, which a surrounding row carries just as a design row
  does — so it landed after `FixedLayerRow::to_row` had forced both
  flags false and undid it. Under K > 1 the LM then moved environment
  0's copy of that surrounding against environment-0-only residuals
  while every other environment kept the compiled thickness, and the
  shared design film was driven to its clamp floor: measured on the
  0.7.5 wheel at 500.0 -> 534.0298 nm, no error. The refusal names the
  environment, the surrounding row (`{env}.fixed[{seg}][{i}]`), the flag,
  the material code, and the two ways out.
- The refusal lives one stage after `to_row`'s because that is the first
  stage that can see the override: `compile_env_rows` now reports which
  rows came from a fixed segment (films only — a half-space row in a
  fixed segment never reaches the override) and `RowAssembly` carries
  that set. Reported by the compile rather than recovered downstream
  from the auto-naming pattern: a name is presentation, a forced flag is
  a contract. The flat door passes an empty set, so it is unchanged —
  there every row is a `LayerRow` with both flags defaulting true and
  there is no invariant to protect.
- `film_flags` (the global map) is deliberately NOT guarded: it is
  applied before the row, so the forced false already wins. Measured as
  a control, and twinned on both sides, because refusing it would break
  every segmented run that asks to optimize its design.

### Changed
- `run_needle`'s `environments` docstring and
  `docs/spectralweave-target-kinds.md` documented M1 as an open gap with
  a "keep surroundings out of `per_film_flags`" instruction; both now
  document the refusal, and the target-kinds section spells out why the
  global map is exempt.

## [0.7.5] - F3.1: docs, worked examples, release

### Added
- `docs/graded-media.md`: spans versus solver rows (16 authored spans
  reach the solver as 80), the Nevot-Croce validity caveat for rtype 5 on
  a graded span (the factor applies at every sublayer boundary, and the
  model's sigma << thickness condition fails first because a sublayer is
  thin by construction), and the thin-gradient floor (a 3 nm gradient is
  three 1 nm rows; exempt from the thin-layer sweep, still below what a
  staircase approximation means).
- `docs/spectralweave-target-kinds.md`: an **Environments** section --
  the `environment=` tag and where it resolves, absent-means-first, one
  merit over all K, parameter identity by film name, the joint needle,
  the xK cost (worse when a surrounding is graded), and what a joint run
  gives up at the thin-layer floor. Plus a coverage-matrix row for the
  per-environment fold.
- `examples/multi_environment_ar.py`: one AR coating bare and laminated,
  optimized jointly, with each finished design scored alone in each
  surrounding (the two runs' own merits carry different residual counts
  and are not comparable), and the same design rebuilt as a program
  document to show identical merit bits.
- `examples/rugate_gradient_vs_discrete.py`: a triangular `FixedSpan`
  rugate against the discrete quarter-wave equivalent -- row count,
  merit and peak R, with the Nevot-Croce and thin-gradient notes inline.

### Changed
- README: a **One Design, Several Surroundings** feature entry, and the
  graded-media entry now points at `docs/graded-media.md`.

### Review PB (no bump; `docs/implementation_review_pb.md`)
- `docs/spectralweave-target-kinds.md` and `run_needle`'s docstring now
  state two open gaps instead of promising they cannot happen. (1) A
  `per_film_flags` entry is applied after the row and keyed by material
  code, so it can re-enable `optimize`/`needle` on a fixed surrounding;
  under K > 1 the LM then moves environment 0's copy against
  environment-0-only residuals while the others keep the compiled
  thickness (M1). (2) After tag resolution a `MeritSpec` keeps only
  `n_envs`, so the run door compares environment COUNTS, not names -- a
  roster written in a different order in the two places is accepted and
  scores every demand against the wrong surroundings (M7). Passing the
  `TargetCollection` to `run_needle` instead of a pre-built spec avoids
  the second entirely. Both are refusal holes and both are open; the
  fixes take their own bump.
- `environments.rs`: four texts that pointed at F2.3 as future work
  reworded, two of which had become false -- `expand`'s docstring listed
  needle insertion among what the driver refuses under K > 1 (F2.3 gave
  it `insert_seed`), and the module header said nothing there was
  reachable from Python (F2.4 made it so). `check_alignment`'s refusal
  message now names `insert_seed` as the only structural move, and its
  twin asserts on that word (M4). Message text only; nothing computed
  changes.
- Measured, not changed: gate 4.6 (c)'s wall-time half, which had no
  number anywhere. With the LM pinned to one iteration, K=2 costs 1.99x
  / 2.14x K=1 across two sessions and each further environment about
  0.74x of a first one. Left to converge the same problem reads 11x,
  which is an iteration count rather than a per-eval cost -- recorded in
  the plan's F2.2 corrections with that caveat (M5).

## [0.7.4] - F2.4: the Python surface and the program sections

### Added
- `run_needle(design={...}, environments=[...])`: named design segments
  defined once and shared by every environment that references them.
  `layers` and `design` are exclusive and exactly one is required. The
  keyword call shapes the same `DesignRequest` a program document is, so
  a run described in a file and the same run described in Python compile
  through one compiler and land on the same merit bits.
- `design_from_program`: a loaded program's `design:` / `environments:`
  sections as `run_needle` kwargs, each row's `material_code` resolved
  against the program's own materials. The material SPEC crosses, not the
  provider's evaluated curve, so it is evaluated once — on the run's grid.
- `run_design_environments` (native): the multi-environment door,
  `DesignRequest` JSON in, result dict out.
- `environment=` on `SpectralTarget`, `AngularTarget` and `ColorTarget`,
  and `build_merit_spec(environments=[...])` to declare the roster. The
  tag is emitted only when set; absent means the first environment, so
  every existing target set keeps its meaning and its JSON.
- Program documents gain `design:` and `environments:` sections (inside
  `sections`, like every other payload). `LoadedProgram` carries both,
  verbatim: the loader checks their shape, the compiler checks their
  meaning.

### Changed
- Program envelope schema is now a readable RANGE: `PROGRAM_SCHEMA_VERSION`
  is 2 and `MIN_READABLE_PROGRAM_SCHEMA_VERSION` is 1, in both homes. A
  document from a newer build refuses saying so; one below the floor
  refuses as stale. A v1 program assembles bit-identically under the v2
  build.
- `run_needle`'s `layers` parameter is optional (it was a required
  positional). Passing neither `layers` nor `design`, or both, is a
  `ValueError` naming which.

### Fixed
- An unknown program SECTION name now refuses, naming it. Every section
  lookup is a bare `get`, so an unrecognised name was dropped in silence —
  which is exactly how a v2 program would have lost its environments on a
  v1 build, and how the next section added would have been lost again.
  This is the defect that makes the version bump necessary; the whitelist
  is what stops it recurring.
- `program_to_dict` dropped the two new sections on the way into Python,
  and the `context=` branch of `load_program` never read them. A section
  that parses and then evaporates is the same failure the bump prevents,
  one layer up.
- `_load_program_native` stamped its rebuilt envelope with a literal
  `schema_version: 1`. Harmless under a point gate; under the range gate
  it would have read a v2 program's sections under a v1 label.

## [0.7.3] - F2.3: the joint needle

### Added
- The needle sweep runs over every environment at once. Each environment
  scans its own assembly with its own fold; every site that lands in the
  shared design is routed by design parameter into a shared bucket and
  summed, so the candidate the sweep picks is the one that helps the
  joint merit rather than the one that helps environment 0.
- `CompiledEnvironments::slot_of_row`, `host_row` and `insert_seed`: the
  two halves of the locus translation, and the split itself. One
  insertion edits the shared design object and every environment's
  template at the same time, at each environment's own row.
- `build_needle_targets_env`: the fold for one environment of a joint
  spec. `build_needle_targets` keeps its signature, becomes the
  single-environment door, and refuses a multi-environment spec rather
  than folding environment 0 and dropping the rest.
- `DesignContext::environments` / `environments_mut`, defaulted to
  `None` — the flat path answers "no environments" and nothing else
  about it changes.

### Changed
- `SpectralInputs::fold` is now `folds`, one per environment in roster
  order. A fold activates one-sided and banded kinds at the operating
  point, and environment 1's operating point is not environment 0's.
- `run_environments` no longer refuses `needles_per_cycle > 0`.

### Refused
- `thin_layer_policy = 'remove'` while K > 1. Elimination is the inverse
  of the insertion F2.3 added, and the compile does not have it: a span
  deleted from the shared design would leave K templates still carrying
  it. Joint runs take `clamp_up_final` — the documented default for
  needle runs anyway — and the context pins the floor as an LM bound for
  the duration, so the span layout moves only through `insert_seed`.
- `enable_cleanup` and `enable_inflate` while K > 1, as in 0.7.2.

### Known limitation
- A seed the optimizer wants to reject parks at the floor instead of
  disappearing, because nothing is eliminated while K > 1. The needle
  loop still stops on its convergence test, and `needles_per_cycle` and
  `max_macro_cycles` still bound the growth, but a joint run can end
  carrying floor-thickness layers a flat run would have dropped.
- The thickness LM still differences its Jacobian under K > 1 (0.7.2's
  decline, unchanged — §4.4 scopes the thickness LM out of this item).
- Still not reachable from Python; that is F2.4.

## [0.7.2] - F2.2: K assemblies, K solves, joint merit

Second rung of Phase B. `[0.7.1]` compiled environments; this release
**evaluates** them. A joint run now assembles the shared design once per
environment, solves each, and scores every demand against the curves of
the environment it was tagged with. Still not reachable from Python - the
`environments=` surface is `[0.7.4]`.

### Added

- **`merit_multi` / `residuals_multi`** over one simulation per
  environment. Residuals concatenate environment-major, insertion-minor;
  the missing-curve penalty is charged per `(environment, key)` PAIR, and
  a pair no demand belongs to is skipped free - it has no curve to miss.
  `merit` and `residuals` keep their signatures and are the one-element
  case, not a separate path.
- **The K-expansion.** The shared design's current thicknesses are
  re-expressed as K stacks through the routing table, per ROW rather than
  per span total, so a graded design film's profile survives exactly.
  Surroundings are never touched.
- **`run_environments`** - the end-to-end joint run.

### Unchanged

A single-environment run takes the pre-0.7.2 call sequence op for op: the
branch is one `Option` on the solver context, read once per eval, outside
every loop. Both bit-exactness fingerprints unmoved, and the recorded
`eval_baseline.json` merit reproduces bit for bit. Every directly measured
bench phase is at or below the `[0.7.1]` baseline.

### Refused

A merit spec whose environment count disagrees with the design's. One
simulation handed to a K-environment spec (it would silently drop every
demand tagged with another environment). Needle insertion, cleanup and
inflate while K > 1 - a structural move changes the shared design's span
layout, and translating it back into the shared object so it propagates
everywhere is `[0.7.3]`. A design stack that no longer matches the span
layout the compile recorded.

### Known limitation

While K > 1 the analytic Jacobian declines and the thickness optimizer
runs on finite differences: the sensitivity pass walks one simulation and
knows nothing about environments. Routing deposits by design-film name is
`[0.7.3]`. Single-environment runs keep the analytic path.

## [0.7.1] - F2.1: environment segment schema and compile

First rung of Phase B (multi-environment optimization,
`docs/plans/multi_environment_plan.md`). **Compile only**: this release
adds the schema, the validation and the routing table. It does not
evaluate anything, and none of it is reachable from Python yet - the
driver loop is `[0.7.2]`, needle/LM routing `[0.7.3]`, the Python
surface `[0.7.4]`.

### Added

- **Named design segments and per-environment segment lists** in the
  native design request. A design segment is defined once and shared:
  every environment references it, so a thickness step or a needle
  insertion propagates to all of them by construction. An environment
  is an ordered list of segments, each either inline fixed layers (its
  own cover glass, substrate, housing) or a reference to a design
  segment.
- **The `design_slot -> (env, span)` routing table.** Parameter
  identity across environments is the film NAME, not its position:
  assemblies of different length cannot share positional identity. Span
  -keyed rather than row-keyed, because after Phase A one design film
  can be a span of up to 64 rows.
- **`environment=` on all three demand kinds**, plus the roster on
  `TargetSet`. An untagged demand is the first environment, which is
  what keeps every pre-0.7.1 target set meaningful without an edit.
- **`validation/benches/synthesis/bench_eval.py`** and its recorded
  `results/eval_baseline.json`: assemble / simulate / merit timed
  separately on one fixed design. Nothing under `validation/benches`
  timed the design pipeline's inner eval, and `[0.7.2]`'s "the K=1
  branch costs nothing" gate needs a baseline that predates the change.

### Refused, each naming what it found

An environment that omits a design segment or references it twice; an
unknown segment reference; a segment carrying both `layers` and
`design` or neither; `optimize`/`needle` set true inside fixed
surroundings (named down to `<env>.fixed[<seg>][<i>]`); one film name
defined by two design segments; a demand tagging an environment that
does not exist (named alongside the ones that do); flat films and
environments in the same request; design segments with no environment
to reach them, or environments with no design segments to share.

### Unchanged

A request with no `environments` takes the flat path unchanged - it is
delegated to, not re-implemented, so there is no second assembler to
drift. Both bit-exactness fingerprints unmoved; the solver, `SimCurves`,
the merit formulas and the fold arms learn nothing about environments.

## [0.7.0] — the differential-phase series, released

The minor marker, taken at the close of the PD series rather than at
F3.1 (main plan §0.3). **This is the first release since `0.6.44`** —
`0.6.45` through `0.6.49` were development rungs that never reached
PyPI or crates.io — so upgrading from `0.6.44` picks up all of the
following at once. Per-item detail is in the sections below; this entry
is the orientation.

### Fixed — the one that matters if you use Δφ

- **The differential-phase reference index now follows the
  wavelength.** Through `0.6.44`, Δφ's reference was one scalar index
  applied at every λ: a dispersive incidence medium was frozen at its
  centre-λ value, which made Δφ wrong by
  `err(λ) = 2π·D·cosθ·[n(λ) − n(λ_centre)]/λ`. Measured on a mildly
  dispersive ambient, that is **0.226 rad** — a silent, systematic
  error in a published release, which is why the series ran before
  Phase B rather than after. A non-dispersive medium is unaffected,
  bit for bit. See `[0.6.45]`.

### Added

- **`PDts`/`PDtp` as first-class `compute()` observables** (`PD_TS =
  1<<49`, `PD_TP = 1<<50`), plus a `differential_phase(*, s_pol=True,
  p_pol=True)` view. Δφ was reachable through the synthesis merit and
  the numpy rotation recipe but not from the simulation surface
  itself. The merit op point and the compute key agree **bitwise**.
  Coherent stacks only, as with the dispersion keys. See `[0.6.48]`.
- **Per-λ reference arrays at the FFI and Python doors**, with a
  permanent guard on the scalar door: a differential demand with a
  scalar reference index warns once, naming the remedy. See `[0.6.46]`
  and `[0.6.49]` for the guard's final shape.
- **A decided, pinned GD/GDD convention over Δφ**: group delay taken
  over the corrected Δφ carries the reference's own dispersion. The
  absolute `GD`/`GDD` keys are unchanged and carry no reference term.
  See `[0.6.47]`.

### Changed

- `MeritSpec::merit` and `::residuals` now **panic** on a malformed
  reference-row length rather than mis-sampling silently. Unreachable
  from Python — the FFI constructor refuses it at build time — and
  documented under `# Panics` for the hand-built Rust caller.
  `curve_sensitivity` deliberately has no such guard: the reference is
  additive and independent of the curve, so it drops out of
  `d(residual)/d(curve)` and that path never reads the reference rows.

### Reviewed

The series was reviewed twice against `docs/implementation_plan_pd.md`,
and both rounds are recorded in `docs/implementation_review_pd.md`
rather than summarized away:

- **Round 1** (§1–§6) found five issues. G1 — the scalar guard warned
  about `n_back`, a side no label can read, so a caller who supplied
  exactly the per-λ front array the guard asked for was still warned —
  was fixed before release, along with G2 (two `total_d` derivations,
  one claiming to be the other).
- **Round 2** (§7) reviewed the fixes against a fresh build and found
  H1: the G2 fix was pinned by nothing. Reverting it left the entire
  battery green, because every test stack sets both half-spaces to
  zero — exactly where the two expressions agree. The twin now exists
  and was watched fail (5.03 rad, exit 1) before being kept.

One verification gap is recorded and still open: the
`nondispersive-bitwise` literals have not been reproduced from outside
the repo (review §5). The claim they carry is held independently by
§1.3's route.

**Closed after the tag** (docs only, no bump): the literals were
re-derived from a fresh `v0.6.44` build in an isolated worktree, using
check 11's own construction, and all four reproduce bit-exactly on both
the scalar-door and engine paths. Review §8 carries the comparison.

## [0.6.49] - PD review applied: the guard's sides (G1) + one `total_d` (G2)

`docs/implementation_review_pd.md` reviewed the whole PD series against
the tree and found five issues; the two that touch behaviour land here,
the three bookkeeping items in the docs commit that follows.

### Fixed

- **G1 (P0): the PD2 scalar guard warned about a side no demand can
  read.** The guard gated on a spec-wide `uses_differential()` and then
  reported *every* length-1 reference row — but `n_back_re` is read
  only under `key.curve.is_back()`, and no differential label maps to
  a back curve (both `PDts`/`PDtp` are front), so the back half of the
  warning was unconditionally a false positive: a caller who did
  exactly what the guard asks (per-λ `n_front`, default `n_back`) was
  still warned, about an index nothing reads. `MeritSpec` now exposes
  `demanded_reference_sides()` (walk the targets, collect
  `key.curve.is_back()` per differential demand); the guard reports
  only demanded sides. Front-only spec + per-λ front + default back is
  now **silent**; a scalar front is reported alone. Twin added: the
  silent case is the one that was wrong, and no check covered it.
  The engine-fill path is unaffected (full-length rows both sides).
- **G2 (P1): the synthesis evaluator's coating thickness `D` was the
  plain sum of all layer thicknesses while the engine's PD keys sum
  the interior `[1..nl-1]`** — the same quantity only as long as the
  half-spaces are zero, which `DesignStack`'s direct door does not
  enforce (`LayerSpec("amb", nk, 100.0)` was accepted). The evaluator
  now uses the interior sum too, so the merit's differential op point
  and `compute(PDts)` agree on *any* stack, not just well-formed ones;
  the comment that claimed the two expressions were one rule is now
  true. Bitwise identical for every legal stack (zero half-spaces —
  adding 0.0 terms cannot move a partial sum).

### Twins

- `guard-sides` (in the PD2 door twin): front-only demand + per-λ
  `n_front` + default `n_back` emits zero warnings; the scalar warning
  names only the demanded side.
- Rust: `demanded_reference_sides_tracks_the_labels` (front/back/none
  from the label table, including the hypothetical back-curve case).
- `compute-observable` (the PD4 twin) gains the half-space case: the
  `PDts`/`PDtp` merit op point is **bitwise** unmoved at ambient/
  substrate thicknesses 999/777 nm. Added in the review's second round,
  which found the G2 fix unpinned — reverting the interior sum left the
  whole battery green; it now fails this check by 5.03 rad.

## [0.6.48] - PD4: PDts/PDtp as first-class compute() observables

The differential phase was reachable through the synthesis merit and the
numpy rotation recipe, but not from the simulation surface itself:
`ScatterMatrix.compute()` had no way to ask for it.

### Added

- **Two request bits**: `PD_TS = 1 << 49`, `PD_TP = 1 << 50`
  (`REQ_PD_TS`/`REQ_PD_TP` in `core_engine.rs`, the `Request` IntFlag in
  `smatrix.py`), emitting keys `PDts` / `PDtp`, wired into
  `expected_keys`. The request-bit smoke test sweeps them end to end.
- **The derivation**: the complex forward-t rows the engine already
  computes for `TS_C`/`TP_C`, minus `reference_phase(λ, n_front_re[λ],
  θ, D, 1)` — `D` the sum of the interior thicknesses (ambient and
  substrate carry zero, the same rule the synthesis evaluator relies
  on) and `n_front_re` layer 0's real index per wavelength (PD1's
  column, at a second call site). No new user plumbing.
- **A convenience view**: `differential_phase(*, s_pol=True,
  p_pol=True)` alongside `complex_amplitudes()` and `dispersion()`.
- **Decisions stated** (per the plan, not discovered): the keys emit
  the **wrapped principal value** in `(−π, π]`, consistent with
  `phi_ts`/`phi_tp` — unwrapping across the grid is the caller's job;
  and Δφ inherits `dispersion()`'s coherent-stacks caveat verbatim.
- `.pyi` stubs for the new request builder (`check_pyi_sync.py`
  blocking); `check_exposure.py` allowlists the row-derivation helper
  with its rationale (the Python surface is the view and the keys, not
  the kernel).

### Twins

- `compute(PD_TS|PD_TP)` equals `apply_reference_rotation` on
  `compute(TS_C|TP_C)` followed by `np.angle`, to 1e-12 (dispersive
  ambient, both polarizations) — the rotation oracle, now spanning the
  new surface.
- `compute(PD_TS)` equals the native merit's op-point Δφ for the same
  stack and grid (bitwise, 0.0e+00) — the two surfaces agree rather
  than merely both existing.
- `expected_keys` round-trips the new bits.

## [0.6.47] - PD3: GD/GDD over Δφ carry the reference's dispersion (decision pinned)

PD1 made the differential reference per-λ. That silently changed what
"GD/GDD over Δφ" means: a frozen scalar index contributed a constant
group delay and exactly zero GDD (which is what made the old doc remark
"finite differences kill the reference anyway" true); a per-λ index
contributes a λ-dependent GD shift and a genuine GDD term.

### Decided

**Option 1 — report GD/GDD over the corrected Δφ.** The series exists to
make the reference correct; carving the dispersion orders out would
reintroduce, one level up, exactly the inconsistency being removed.
Option 3 (both sets of keys) is deferred, not refused: additive keys can
be added later without revisiting the decision.

### Changed

- `docs/spectralweave-target-kinds.md`: the stale "finite differences
  kill the reference anyway" sentence replaced with the corrected
  recipe and the analytic corrections (`GD_ref = (n + ω·dn/dω)·D·cosθ/c`,
  the reference's GDD term `2·dn/dω + ω·d²n/dω²`).

### Numeric consequence (dispersive ambients)

For a dispersive incidence index the GD over Δφ acquires the
λ-dependent shift `−(n_inc + ω·dn_inc/dω)·D·cosθ_inc/c` and the GDD over
Δφ acquires `−(2·dn_inc/dω + ω·d²n_inc/dω²)·D·cosθ_inc/c` relative to
the pre-0.6.45 numbers. For a constant ambient index: **no change at
all** — GD shifts by the constant `−n_inc·D·cosθ/c` only in the
differential keys (as always) and GDD is untouched.

### Added

- Parity check `gd-gdd-convention` (13th check in the PD file): a
  two-point hand case through the library `dispersion()` is
  bitwise-exact (film index = substrate index kills the Fabry-Perot
  denominator, so `arg(t) = n·ω·D/c` exactly); GDD over Δφ is non-zero
  for a dispersive ambient and matches the analytic
  `−3·B·ω·D/(2π²c³)` at the uniform-grid centre to 1e-12, and is zero
  (1e-12) for a constant ambient; the absolute GD/GDD keys carry no
  reference term (the 0.6.46 behaviour, pinned bitwise for the linear
  row).

## [0.6.46] - PD2: the doors accept a per-λ reference; the scalar door gets its permanent guard

After PD1 the native path is correct by construction (it reads the
stack). The remaining way to hand a frozen scalar index to a dispersive
medium is a hand-assembled `SimCurves` — and every door that takes a
reference index still took a scalar.

### Changed

- **`SimCurves.__init__`** (`n_front`/`n_back`): accept `float | FloatArray`
  — a float becomes a length-1 broadcast row, an array passes through for
  the native length rule (length 1 or `len(wavelengths)`; anything else
  refuses via `reference_length_issue`, naming both numbers).
- **`reference_rotation`**: `n_inc` accepts `float | FloatArray`, length 1
  or `len(wavelengths)`; a wrong length refuses naming both numbers. The
  core kernel is per-λ (`n_inc: &[f64]`, broadcast) — the per-λ arithmetic
  is bitwise the pre-PD2 scalar kernel when the row is constant.
- **`apply_reference_rotation` / `sim_curves_from_arrays`**: pass-through
  (thin, as before — the validation stays native).
- `.pyi` stubs updated (`check_pyi_sync.py` green).

### Added

- **The permanent guard (the part that outlives PD1):** when a scalar
  (length-1) index meets an active differential demand, the door emits a
  `UserWarning` naming the quantity, the supplied value, and the remedy
  — pass the per-λ array (ground rule 6: announced, never silent;
  ground rule 7: ASCII). A warning, not a refusal: air is a legitimate
  scalar and by far the common case, and the door cannot know the stack
  it was not given. Python's warning registry dedupes per call site, so
  an LM loop sees it once. Wired at `reference_rotation` (covering
  `apply_reference_rotation`'s numpy path) and at the merit/residuals/
  `build_needle_targets` FFI doors (a hand-built sim + a differential
  spec). The engine fill (per-λ rows) never warns.
- Parity check `dispersive-rotation-door`: scalar == length-1 array
  bitwise at both doors; the dispersive oracle (native differential ==
  per-λ-rotated absolute) holds at 1e-12 — what keeps the numpy path a
  real oracle instead of a co-drifting copy; the guard warns once with
  the remedy and is ASCII; the length mismatch refuses naming both
  numbers. Rust unit tests for the kernel's broadcast rule and refusals
  (551 lib tests).

## [0.6.45] - PD1: the differential-phase reference index follows the wavelength

The differential-phase reference was a **scalar frozen at the centre
wavelength** (`SimCurves.n_front_re`/`n_back_re: f64`, picked at
`nw / 2` from the per-λ stack cache). A dispersive incidence medium
silently yielded a wrong Δφ by
`err(λ) = 2π·D·cosθ·[n(λ) − n(λ_centre)]/λ` — up to 145° (N-SF11, 3 µm,
380–780 nm), and live in a published release (0.6.44 on PyPI,
2026-09-14). Two comments called dispersive ambients "pathological";
there was no runtime warning, no refusal, and no test pinning the
approximation.

### Changed

- **`SimCurves`'s reference indices are per-λ rows** (`Arc<[f64]>`,
  length `wavelengths.len()` or 1). Length 1 broadcasts — the
  load-bearing decision: `Default` stays `[1.0]` so `total_d = 0`
  reproduces absolute phase bit-for-bit, the air case stays cheap, and
  the scalar FFI door keeps working without a shim. Any other length is
  a refusal (`SimCurves::reference_length_issue`) naming both numbers,
  wired into `merit`/`residuals`/the needle fold.
- **The fill collects the whole column** (`evaluator.rs`): layer 0's
  real index per wavelength (and the substrate's for the exit
  reference) instead of the centre-λ element.
- **The readers sample the row at the demand wavelength** — merit's
  inner loop (same two-pointer interpolation class as the complex rows,
  sharing its bracket state), the needle op-point sample, and both
  gain-shift sites (`interp_n_row`: length-1 broadcasts, longer rows
  interpolate against `sim.wavelengths`). The demand grid and the sim
  grid need not coincide; the FD `dM/dD` twin is the gate that pins
  this.
- `needle_pass.rs`'s `n_inc.unwrap_or(1.0)` is now an `expect` with the
  reachability proof in its message: `sample_op_value` runs only under
  `(Some(sim), Some(rows))` and `n_inc` derives from that same sim, so
  the `None` arm was unreachable all along (plan §7 decision 2,
  resolved during PD1 as the plan directs).
- `synthesis_merit.rs`'s PyO3 ctor wraps the scalar doors into length-1
  rows (compile-level; the array door is PD2).

### Fixed

- Δφ, merit, residuals, LM trajectories and the needle gain shift for a
  spec with a differential demand **and** a wavelength-varying incidence
  (or exit) index — the only licensed numeric change (plan §2); a
  non-dispersive medium is **bitwise unchanged**, pinned by the new
  `nondispersive-bitwise` parity check (merit and residuals are the
  recorded 0.6.44 literals through both doors).

### Added

- Parity checks `dispersive-hand` (dispersive ambient through the engine
  fill: hand targets embed the per-λ reference — merit ≈ 0 post-fix,
  ≈ 300 with the frozen index, the check the defect would have failed)
  and `nondispersive-bitwise`; Rust twins for the gain-shift
  interpolation at the demand wavelength and for needle-site invariance
  under a differential demand (the shift is uniform in z, so the argmax
  never moves); `SimCurves` broadcast/refusal unit tests. 550 lib
  tests.

## [0.6.44] - C3: the sweep — message rendering, stale comment, doc note, bookkeeping

Everything the status review found that does not change behaviour except
message text (findings 3–7), plus the tool that guards the class.

### Changed

- **Five message literals fixed** (a lost `\` line-continuation leaves
  the next line's indent INSIDE the string): the two F0.2 refusals the
  user meets when the ceiling licence fires
  (`clamp_all: ... ceiling - refusing rather than rescaling a profile`;
  the `NeedlePipeline:` variant), the F1.1 `MixRule` parse refusal (V1 —
  one typo away on every gradient a user authors), the pre-existing
  `extrap='error'` range message, and the pre-existing
  `levenberg_marquardt: analytic jacobian` message — which also carried
  a non-ASCII `×`, now ASCII (ground rule 7). The two existing refusal
  twins upgraded from substring asserts to the EXACT full message; a
  new Python-side twin pins the `MixRule` wording (the only one of the
  five reachable without constructing a stack).
- **Two non-encodable characters left production messages** (ground
  rule 7's actual bite — cp1252 consoles cannot encode them):
  `stagnation_window must be >= 2.` (was `≥`) and
  `material has non-positive n={n_real} at wl={...} nm` (was `λ`).
- **New tool: `tools/check_message_whitespace.py`**, wired into CI and
  the battery, staged exactly as amendment 3/V2 specifies:
  space runs of 3+ inside string literals are BLOCKING with an empty
  allowlist (the class is unambiguous); cp1252-unencodable characters
  are BLOCKING; the encodable non-ASCII still outstanding (—, ×, · in
  seven production files, ten literals) is written out as an explicit
  advisory allowlist that F3.1's exposure re-audit retires. Test
  regions (`#[cfg(test)]`, `rust/*/tests/`) and comments are skipped — the two
  space-run sites the review's naive scan would have flagged are test
  assert messages, and one quoted string sits inside a trailing
  comment.
- **A stale doc comment retired:** the `Param` enum still said the rate
  modes "stay homogenized until F1.7" — F1.7 shipped; the comment now
  states the refresh reconciliation in present tense.
- **`Layer::sub_layer_count` documents its gradient exception:** a
  gradient-only layer reports 1 (the count rule needs the wavelength
  grid and both endpoint spectra, which a bare `Layer` does not have);
  the emitted count is the expansion's `gradient_sub_layer_count`.
  Whether the getter gets an honest sentinel is queued for F3.1.

### Bookkeeping

- F1.4's master-table row struck (ground rule 1's only miss).
- §8 decisions 8, 9, 10, 12, 13, 14 stamped RESOLVED with their
  resolving commits (decisions 1–5 and 11 stay open for Phase B).
- Ground rule 7's reality-check paragraph retired (V8): the em-dash
  exception it carved out no longer exists; the paragraph now points at
  the scanner tool and the allowlist.
- `uv.lock` tracked (V7: no gate consumes it — it records what the dev
  environment resolved to when digests were recorded, which C1 made
  load-bearing provenance; the alternative was `.gitignore`, not the
  status quo).

## [0.6.43] - C2: an interface-carrying film clamps up instead of being deleted

The review's measured case (finding 2): under a clamp-up policy, a plain
film that carries an interface slice was deleted whole where the same
film without an interface was clamped up and kept - because the clamp-up
branch keyed on a ROW COUNT (`end - start == 1`) instead of the repo's
own predicate for "one physical layer" (`is_singleton_bulk`, N1), and an
interface film is two rows and still one layer. Deleting is the opposite
of what the policy's whole purpose (U1) asks for.

### Changed

- **The clamp-up branch keys on `is_singleton_bulk()`** - one bulk row,
  interface slice optional - so a plain interface-carrying film under
  `ClampUpFinal`/`ClampUpAlways` is CLAMPED: the slice stays bitwise
  (B1: it belongs to the authored thickness) and the bulk row takes
  `clamp_min_nm - t_slice`, landing the slice-inclusive total exactly on
  the floor. The one-row plain film is the `slice == 0` case of the same
  arithmetic - bitwise the old behaviour, which is why every existing
  twin holds unchanged.
- **The graded/span deferral is fully retired.** F0.3's sentence ("a
  multi-row span is removed whole ... clamping a span up is F1.6's scale
  operation and does not exist yet") predates F1.6; the scalable branch
  (F1.6) already scales scalable spans to the floor, and C2 completes
  the singleton-bulk case. What still is removed whole, deliberately:
  every other multi-row span - a plain multi-row span has no principled
  way to choose which row grows, and a non-scalable profiled span has
  neither (the rate modes read absolute depth; scaling rows the user
  did not offer for scaling is not a clamp). Removal stays reported.
- **`ClampUpAlways`'s LM bound is stated as conservative, not exact**
  (doc comment): the floor binds the BULK row (the LM parameter); an
  interface-carrying film therefore settles one slice thickness above
  the floor (`bulk >= clamp_min_nm`, total `>= clamp_min_nm + t_slice`).
  Conservative by construction - a floor is a minimum, not a target -
  and pinned by the interface bound twin.
- **Scoped out:** the CEILING side keeps refusing a plain
  interface-carrying span above `clamp_max_nm` (the refusal is loud,
  names both sides, and is safe); shrinking a plain film to the ceiling
  would also be well-defined, but changing F0.2's refusal semantics is
  not what the review found. Recorded as a known asymmetry.

### Twins

- The review's probe scenario, promoted: 0.8 nm TiO2 with a 0.2 nm
  interface under a 2.0 nm floor - `ClampUpFinal` and `ClampUpAlways`
  variants keep it (slice bitwise 0.2, bulk `2.0 - 0.2`, total within
  the established post-scale tolerance, empty report); `Remove` still
  eliminates it whole (2 rows, reported).
- Above the floor: untouched (no clamp branch reachable, slice rows
  never move through clamp).
- The LM bound's interface variant: a thin interface-carrying film under
  `ClampUpAlways` binds the BULK row at `clamp_min_nm` (total
  `clamp_min_nm + slice`), merit history monotone over five sweeps,
  never deleted; under `Remove` the carrier goes whole and the lead
  survives.

## [0.6.42] - F1.5: config rows and the Python gradient surface

The last transport gap: gradients were constructible only through the
synthesis film-dict door (which carries Python-evaluated spectra). Now
they ride the config documents and the native `Layer` itself.

### Added

- **`LayerRow.gradient`** - the named SPEC (`GradientSpec`: material_a,
  material_b, ema, mode, shape, sublayers), not the `GradientJson`
  transport: both native paths resolve endpoint spectra themselves (the
  config path through the provider at expansion, the design path
  through the library's nk table), and a config document holding
  evaluated spectra would freeze a library snapshot into a hand-edited
  file. The `nk_b`-carrying `GradientJson` stays the Python film-dict
  door's transport (A5), where no library exists.
- **The native `Layer` constructor takes `gradient=`** (the named spec
  dict), with a getter returning the full serde shape and a setter
  gated on the probe (a rejected spec leaves the layer as it was).
  `sublayers` and `shape` are optional on the ctor (F1.1 ships
  `Linear` only).
- **`builders.py` passes `gradient` through**, so config documents
  reach gradients end to end: `LayerConfig` -> `layer_from_config` ->
  `Layer` -> `Structure::expand` (endpoints resolved through the
  provider).
- **Validation is native everywhere** (no Python pre-checks): the
  spec's own rule surface runs at the config door (`LayerRow::validate`),
  at the Layer constructor's gate, and at expansion; the
  two-profile-engines conflict (gradient + inhomogen) now runs in the
  LAYER's rule surface, where both flags live, instead of surfacing
  first at expansion.

### Changed

- **The exact endpoints of every EMA kernel are bitwise the endpoint
  spectra.** The structure path's nominal-expansion gate refused the
  f=1 mixture row as `k < 0` - the Newton root stops on
  `|d_eps| <= tol` and can leave a `-4e-22` artifact in a lossless
  endpoint's k. There is no mixture at f = 0 or f = 1: the kernel now
  short-circuits to the endpoint itself, in the SHARED kernel, so both
  the synthesis path and the `evaluate()` oracle change together and
  every existing bitwise twin still holds.
- **The nested spec refuses unknown fields too** (`deny_unknown_fields`
  on `GradientSpec` and `MixRule`, extending N9): a typo'd key inside
  `gradient` must not vanish. The state path is safe - its version
  gate runs before any field parses.
- **`build_design` gates a gradient row's rule surface** (the same
  `ValidationIssue::gate` the ArrayFilm driver runs) and its
  background rule counts gradient carriers as profiled films (the A4
  mirror).

### Twins

- Rust: a gradient row through `build_design` expands bitwise the
  direct EMA oracle (count rule included); bad specs refuse at
  `LayerRow::validate`; a typo'd nested key refuses loudly.
- Python: the ctor/getter/setter/state round-trip; native validation
  wordings at the door (duplicate materials, out-of-range fractions,
  two profile engines, unknown kernels); the config door end to end
  (`LayerConfig` with gradient -> validate -> expand); and the
  structure-path expansion - a gradient `Layer` in a `Structure`
  against a gridded provider, rows bitwise the direct `evaluate()`
  oracle (the F1.1 branch, Python-reachable for the first time).

## [0.6.41] - F1.4: schema v2 and a readable version range

Every state file on disk reads through this gate; getting it wrong is a
data-loss-shaped bug, so the item is P0 and landed before the Python
surface.

### Added

- **A readable RANGE, not a point** (F1.4's whole point): the state gate
  now accepts `[MIN_READABLE_SCHEMA_VERSION, SCHEMA_VERSION]` =
  `[1, 2]` and refuses outside it with a reason that names direction -
  below the range is *stale* ("refusing a stale state"), above the
  range is *newer* ("refusing a state written by a newer build -
  upgrade navette to read it"). Rationale: v1 states written by
  pre-F1.4 builds are still truthfully readable (the v2 keys ride
  additively and reconstruct on read), so refusing them would discard
  user data for no misread; the top of the range stays closed because a
  newer writer's keys/meanings cannot be known here.
- **`MIN_READABLE_SCHEMA_VERSION`** at both halves (Rust `version.rs`,
  Python `structure/types.py`) and in the sync probe.
- **`_accepted_range`** in `test_request_bits.py`: probes the gate and
  asserts the pair (oldest readable, current writer) plus CONTIGUITY
  (a hole - a version refused while both neighbors are accepted -
  would mean a refusal no schema policy justifies). The state test
  asserts `(1, 2)`. The old single-version helper stays for the
  program gate, which is still a point until F2.4.

### Changed

- **`SCHEMA_VERSION` = 2 at every half in the same commit**: the Rust
  constant, the Python constant, and BOTH policy comments (the
  pre-F1.4 text said additive keys were safe without a bump; F1.4's
  newer-writer hazard justifies bumping anyway, and the rewritten
  comments now say so - including the fingerprint test's own header).
- **`gradient` rides the state additively, only when `Some`** (the
  same ethic as F1.3's `inh_mode`): a stack with no gradient
  serializes byte-identically at v1 and v2 apart from the version
  tag. The nested `GradientSpec` key set has its own fingerprint
  entry - a nested object would otherwise have weaker protection
  than every top-level key.
- **`Layer`'s hand-written serde** (no derives to hang defaults on -
  the plan's N4 correction): the field is read with an explicit
  `None` fallback matching the surrounding style, and the map length
  became `13 + rate_capped + has_gradient` instead of a bare 13.
- **The config-path state wrappers** (`StructureState` / `ArchitectState`
  in `config/models.py`, the Python half A3 understated further)
  pre-check the same range through one shared helper.

### Inverted deliberately (the amendment's test list)

- `test_stale_schema_versions_refused`: v1 now loads (the pre-F1.4
  `SCHEMA_VERSION - 1` case), v0 refuses as stale, v+999 refuses with
  the newer-build reason, untagged still refuses as malformed.
- The v1 fixture (committed at `ab7e27b`, BEFORE the gate changed -
  R5) loads and expands bit-identically to its live-built equivalent;
  the loaded state re-tags itself at the current version on write, so
  the oracle is stable across the bump.
- A hand-built v2 gradient layer state round-trips v2 -> v2 through
  `from_state` (the constructor surface for gradient is F1.5, but the
  STATE surface exists at F1.4 and is pinned at F1.4).

## [0.6.40] - F1.7: profile refresh for the rate modes, at construction only

A rate-type profile is a function of absolute thickness, so scaling the
span makes its stored nk stale. The rate modes now keep their profiles
like the fixed modes do, and the profile is re-derived at the
construction points.

### Added

- **`SpanRecipe`** - everything one iteration of `expand`'s emission
  loop reads, captured at the door that built the span (B2: the
  extraction was done at F0.1, so refresh and construction are the
  SAME code, not two functions that agree). Rides
  `DesignStack::recipes`, aligned index-for-index with `spans` (never a
  map keyed by `logical` - a second index that renumbers on
  `remove_film` is the failure mode F0.1 exists to prevent). `None`
  for every plain film. Any mutator that alters a span's row set drops
  its recipe; spans that only shift keep theirs; the merge's intra-span
  case (the RateCapped saturated tail merging into itself) keeps the
  recipe because refresh rebuilds from the carrier at the new total.
- **`refresh_profiles(&mut stack)`** - rebuilds every rate-mode span's
  profile from its current total. Called at exactly THREE places:
  `from_design`, the top of each macro cycle (before the budget check
  and needle scan), and after the final clamp pass. NOWHERE else -
  never inside an LM round, never inside the Jacobian, never inside a
  needle scan or a cleanup trial (U4). The refresh-count twin pins the
  counter at exactly `1 + N cycles + 1` (5 for a 3-cycle pipeline); a
  fourth call site later makes it fail loudly.
- **A span whose rows still sum (bitwise) to the measured total of its
  last emission is skipped** - the marker is exact staleness detection,
  not a heuristic: refresh at construction is a bitwise no-op for
  freshly built spans, and the marker re-measures from the new rows
  after every move (one-step convergence, no ulp chasing).
- **Rate spans are now LM parameters too** (completing F1.6's licence):
  the from_design keep-profile branch widens from the two scale-free
  modes to ALL profiled modes. `optimize = true` never means
  "homogenize me" for any profiled film.
- **B3's refusal:** a recipe carrying `apply_errors: true` is refused
  by `refresh_profiles` (refresh would re-roll the error ensemble);
  unreachable through `from_design` today, proved shut by the twin
  before anything can open it.

### Changed (the completion of F1.6's one non-additive part)

- **`optimize = true` on a rate-profiled film keeps the profile** (was:
  homogenize with a warning, for the two RateCapped modes). The
  homogenize warning still fires for the posture that always
  homogenized (`optimize = false, needle = true`) and for pinned
  backgrounds. Within one macro cycle the rate span's profile is the
  one built at that cycle's top (the plan's staleness bound: a 5 nm
  excursion moves the end fraction by `0.05 * rate`; the drift twin
  measures the merit correction at under half the excursion's own
  response and under 3% of the merit); every merit number shown to the
  user is evaluated after a refresh (refresh point 3).

### The Jacobian stays analytic

With the profile frozen for the solve, the stack LM minimizes genuinely
is the frozen-profile stack, so F1.6's exact contraction is the exact
derivative of the actual objective. U4 buys the analytic Jacobian for
rate modes twice over.

### Twins (all landed)

- **Refresh-count (the item's most important test, written per the
  plan's sequencing):** a moving mock over 3 macro cycles - the counter
  reads exactly 5.
- **Refresh-correctness, gradient engine, full bitwise:** 100 -> 200 nm
  (count 6 -> 12, the ceil re-derived) equals the same span expanded
  from scratch at 200 nm, films bitwise; with an interface slice, the
  slice is re-carved (not scaled) bitwise.
- **Refresh-correctness, legacy engine, full bitwise:** the scenario is
  chosen to make that honest - at 100 nm with rate 0.25 the count is 16
  (binary-exact 6.25 nm rows), so the hand-doubled total is exactly
  200.0 and both sides run the same code at the same extent, the delta
  clamped on both sides. A generic (t, rate) pair only reaches 1-ulp
  agreement (N7: the legacy uniform split sums to the carrier to float
  precision) - the plan's flat "bitwise" holds where the arithmetic is
  exact and is documented where it cannot.
- **Fixed-mode no-op:** bitwise identity over 300 seeded randomized
  stacks (mixed Fixed graded / FixedSpan gradient / plain films). The
  plan's citation points at `test_differential.py`, which is the
  structure engine's py-rs comparison - the synthesis-level randomized
  differential runs Rust-side, same seeded style (plan correction).
- **Row-count-frozen:** a real solve moves 221 -> 260 nm (across the
  12/13 boundary): the count stays 12 through the solve and the refresh
  re-derives 14. Frozen-then-healed is the whole U4 contract.
- **Drift-bound + reported-merit:** the staleness correction (6.9 on a
  26.9 excursion response, 2.7% of the merit); `final_mf`'s snapshot is
  bitwise the returned (refreshed) stack.

## [0.6.39] - F1.6: one thickness parameter per graded span

A profiled layer is one physical layer, so its total thickness is one
number the optimizer moves. The LM parameter list stops being a row
list.

### Added

- **`Param` (`Row` | `Span { rows, fractions }`)** replaces the
  parameter row list at the heart of `optimize_thicknesses`. A scalable
  span - at least two bulk rows, every one optimize-flagged, the shape
  `from_design` gives a profiled (`inhomogen` or `gradient`) carrier
  with `optimize = true` whose mode scales exactly (`InhMode::Fixed`,
  `GradientMode::FixedSpan`) - contributes ONE parameter: its total
  thickness `D`, written back as `phi_r * D` per bulk row with the
  fractions frozen at build time. The interface slice is never scaled
  (an interface property, not part of the thickness). A stack with no
  scalable span produces only `Row` params, in the same order as the
  old row list - bit-identical, asserted.
- **Bounds are span quantities:** `lb`/`ub` apply to `D` (the LM sees
  one bound per physical layer). The manufacturing ceiling that never
  fired on a graded film now fires correctly, as a bound.
- **The F0.2 construction refusal narrows to PINNED spans**
  (`NeedlePipeline::new` and `clamp_all`'s pre-scan): a scalable span
  above the ceiling is CAPPED instead (fractions preserved,
  `spans_capped` reported).
- **F0.3's span branch lands:** an under-thickness scalable span under
  an active clamp-up posture is SCALED to the floor (fractions
  preserved) instead of removed whole; under `Remove`, and for pinned
  spans, F0.3's rules stand unchanged. Clamp-ups remain non-Report.
- **The analytic Jacobian stays analytic, exact, cheap:**
  `assemble_jacobian_mapped` contracts the per-row deposit columns into
  per-parameter columns - a `Row` param is a one-element group with
  weight 1.0 (bitwise the old assembly, asserted), a `Span` param is
  `sum_r phi_r * col_r` in row order (bitwise the closed form,
  asserted). The deposits cost O(1) per row, so a 57-row span costs 57
  O(1) terms, not 57 simulates. The map is positional (deposits laid
  out in parameter order) and length-checked instead of indexed blind
  - the Python merge/optimizer differential caught the first version
  indexing by film row.

### Changed (the one non-additive part)

- **`optimize = true` on a profiled film no longer means "homogenize
  me"** (licence item 1): the film keeps its full profile as one
  scalable span and its thickness moves during optimization. The
  homogenize warning stops firing for these films (item 2); it still
  fires for the posture that always homogenized (`optimize = false,
  needle = true`) and for the RateCapped modes, which stay homogenized
  until F1.7 (they depend on absolute depth, not fractional position).
  `optimize = false, needle = false` (pinned background) is untouched.
- Needle host selection unchanged: a scalable span is still not a host
  (scaling preserves a profile; splitting a foreign row into it does
  not).

### Twins (all landed)

- **Scale-invariance (written before the feature):** a FixedSpan
  gradient at 128/256 nm with the count pinned at 4 - nk rows bitwise
  equal, thicknesses exactly doubled (binary-exact split); the legacy
  engine has no count override and no doubling preserves its count
  (t^0.4 grows 1.32x per doubling), so its twin uses the same-count
  pair 100/120 nm (both ceil(t^0.4) = 7): nk rows bitwise equal, every
  row thickness bitwise `t/sub`.
- **Jacobian, two ways:** the span column against a central difference
  on D to 1e-6 relative; and against `sum_r phi_r *` (the row columns)
  BITWISE. Plus the no-span assertion: the mapped assembly of an all-Row
  stack is bitwise `assemble_jacobian`.
- **Profile preservation:** after a real LM solve, every row is bitwise
  the frozen-fraction writeback of the returned D and every nk row is
  untouched; the re-summed fractions agree to 1e-12 (the plan's literal
  "d_r / D unchanged to the last bit" is not well-posed - D is
  re-summed from the rows - so the writeback contract is the bitwise
  oracle).
- **Interface:** a scalable span's slice row thickness is unchanged by
  the scale; the span partition holds.
- **Bound:** the F0.2 cap case (1000 nm graded, ceiling 300) marked
  optimize=true constructs with no warning, `NeedlePipeline::new`
  builds, the LM never returns D above 300, the clamp caps the span
  with fractions preserved.

## [0.6.38] - F1.3: `InhMode::RateCapped` - thickness-relative single-material drift

The legacy scaling drift gains an application mode: the grading
strength can grow with thickness (a deposition drift) instead of being
a fixed property. The frozen legacy arithmetic is untouched - the mode
only replaces the NUMBER that flows into it.

### Added

- **`InhMode`** (`Fixed` default | `RateCapped { rate, ref_thickness,
  cap }`): `delta_layer = min(rate * thickness / ref_thickness, cap)`
  (D2's formula) as a pure function of thickness, exposed as
  `Layer::delta_layer` - the ONE source every reader consumes
  (F1.3/B6): the emission, both row-count predictions (`Structure` and
  `Architect`), the advisory message and the PyO3
  `Layer.sub_layer_count` getter all read it, so a predicted count and
  an emitted count cannot diverge between modes. Pinned by
  prediction-vs-emission twins on every Rust reader and a getter twin
  in Python.
- **The frozen combination order, unchanged:** `delta_nominal =
  (delta_layer + group.inh_delta_summand) * 0.5`; in RateCapped only,
  `delta_nominal` is clamped BEFORE the stochastic error draw - the
  cap binds the nominal, the noise is allowed to exceed it (clamping
  draws would bias Monte-Carlo statistics). Pinned by a statistical
  twin over 200 seeds: the mean drawn delta sits AT the clamped
  nominal and individual draws exceed it.
- **Validation** in `property_issues`: rate finite, `ref_thickness >
  0`, `cap in [0, 1]`; an advisory when the mode is set while
  `inhomogen` is false (the mode is inert there).
- **Python surface:** the native `Layer` constructor accepts
  `inh_mode` (`"fixed"` default, or `{'RateCapped': {'rate':...,
  'ref_thickness': 100.0, 'cap': 0.3}}`); the `inh_mode` getter
  returns the tagged form (`None` = Fixed).

### Notes

- **The nominal clamp is two-sided (`[-cap, cap]`), not the plan's
  literal `[0, cap]`:** the plan's own validation says "rate finite
  (sign free - negative inverts the drift direction)", which a one-
  sided clamp would nullify (a negative rate's drift would clamp to
  zero and vanish). The two-sided clamp keeps the saturation semantics
  (positive rates saturate at `+cap` exactly as the twins test) AND
  the inverted-direction reading (the frozen ramp runs inverted).
- **`inh_mode` rides the state ADDITIVELY** (serialized only when not
  `Fixed`): the PyO3 constructor makes the mode user-reachable at this
  version, and the repo's own additive-key policy (types.py /
  test_roundtrip.py's fingerprint comment) covers old readers - every
  pre-F1.3 state parses unchanged and a Fixed layer's state is
  byte-identical. F1.4's version bump then covers the newer-writer
  hazard for this key and `gradient` together. The fingerprint test
  gained the RateCapped key list.
- `inhomogen` is NOT deprecated: `Fixed` is the default and the
  existing randomized differential pins (Fixed rows bitwise vs the
  legacy oracle) keep passing unmodified - the proof the frozen path
  is intact.

## [0.6.37] - F1.2: gradient `RateCapped` - thickness-relative slope with caps

The second gradient profile mode: `f(z) = f_start + rate * (z /
ref_thickness)`, clamped to `[f_min, f_max]` - the slope is
thickness-relative, so a thick film saturates into a flat
pure-material tail where the cap binds.

### Added

- **`GradientMode::RateCapped { f_start, rate, ref_thickness, f_min,
  f_max }`** with the full validation set: `rate` finite (sign free),
  `ref_thickness > 0`, `f_start` and both caps in `[0, 1]`,
  `f_min <= f_max`, and `rate == 0` refused as the degenerate span
  (the same rationale as FixedSpan's `f_start == f_end`). Serde rides
  the one-key tagged form, pinned by a test.
- **The Python door's mode discrimination (A9.2):** the gradient dict
  needs exactly one slope spelling - `f_end` (FixedSpan) or `rate`
  (RateCapped, with optional `f_start`/`ref_thickness`/`f_min`/
  `f_max` defaults 0/100/0/1); supplying both is refused (`gradient:
  'rate' and 'f_end' are the same slope - give one.` - the plan's
  ASCII form; its 'delta' spelling maps to the FixedSpan endpoint key
  at this door).
- **The saturation physics, bitwise:** tail rows where the raw profile
  value exceeds the cap are EMA at exactly `f_max` - bitwise pure
  `material_b`; negative rates saturate at `f_min` the same way.
- **The N2 merge twin:** the saturated film through the nk-keyed
  merge - the tail rows collapse, the span stays one contiguous run of
  the carrier material, and the simulated merit is unchanged to float
  precision.

### Notes

- **The cap's binding is bitwise where the raw value strictly exceeds
  it, not where the arithmetic lands one ulp short:** `0.1 + 0.3*3.0`
  is `0.9999999999999999` in IEEE, so the plan's knee example
  ("saturates at 1.0 from z = 300 nm") is idealized - at z = 300 the
  profile value is the raw value. The twins use `rate = 0.25` (binary
  exact) where the saturation is unambiguous. The BEHAVIOUR (flat
  pure-material tail where the cap binds) is unchanged.
- **The merge is a no-op to ~1e-15 relative, not bitwise:** folding two
  rows into one makes `exp(i*d1)*exp(i*d2)` differ from
  `exp(i*(d1+d2))` at the last ulp. The plan's "bitwise equal spectra
  before and after the merge" is pinned as float-precision equality.
- RateCapped is NOT scale-free (F1.6's table): its fraction depends on
  absolute `z` - the scalable-span plumbing and the profile refresh
  are F1.6/F1.7's territory, untouched here.

## [0.6.36] - F1.1: the gradient data model, `FixedSpan` expansion, the homogenize path

A mixture gradient - a film whose composition interpolates between two
materials through an EMA kernel - becomes a first-class layer property,
safe from the first commit: the default is `None`, and the plan's own
defect ledger (B1's lesson) is why the bookkeeping half shipped first.

### Added

- **`GradientSpec`** (`structure/gradient.rs`): `material_a` (host,
  `f = 0`), `material_b` (inclusion, `f = 1`), `ema` (a whole `MixRule`
  - all six in-tree kernels work on day one, parameters live in the
  variants, A6), `mode` (`FixedSpan { f_start, f_end }` - `RateCapped`
  is F1.2), `shape` (`Linear` - enum-ready per D9.2), `sublayers`
  override (clamped `[2, 256]`, advisory outside). The docstring
  obligation (D9.4) lives on the module: `material_a` is ALWAYS the
  host and `f` is the volume fraction of B in A - Maxwell-Garnett and
  Mori-Tanaka are host/inclusion-asymmetric, so an A<->B swap is not
  the same physics. Kernels are called inclusion-first, mirroring
  `MaterialSpec`'s argument order.
- **`Layer.gradient: Option<GradientSpec>`** (default `None`), with the
  self-contained validation in `property_issues` (the one rule surface,
  N5): `gradient` + `inhomogen` refused (two profile engines, message
  names both and points at the alternative); `f_*` outside `[0, 1]`
  refused; `material_a == material_b` refused; `f_start == f_end`
  refused (N2 - a degenerate span merges away silently otherwise);
  sublayers outside `[2, 256]` advisory (clamped at expansion).
- **The expansion branch** beside the legacy graded branch, with the
  resolved sublayer rule `n = clamp(ceil(t / max_step), 3, 64)`,
  `max_step = min(20 nm, lambda_min / (10 * n_max_re))` (D3's "exact
  formula fixed at implementation"), the last sublayer absorbing the
  float remainder so `sum(d) == thickness` exactly (N7), INCLUSIVE
  endpoint sampling (first sublayer IS `f_start`, last IS `f_end` -
  the legacy graded convention, and the one that makes the endpoint
  rows bitwise thickness-independent), row-order inversion as the
  physical flip, and roughness/type on the first emitted sublayer
  only. Endpoints resolve through the provider under their own
  materials' groups (scaling, and the nk error channels when errors
  are on); the delta channels are inert on this branch (mixtures have
  no inh_delta).
- **The design path (A5):** `ArrayFilm.gradient`
  (`Option<GradientJson>`) with `nk_b` evaluated Python-side in
  `_film_dicts` (the design path has no materials library; the
  inclusion spectrum rides the film dict, registered under a
  `<film>~b` provider key); `material_a` defaults to the film's own
  name - its registered nk IS the host spectrum; absent `nk_b` or a
  shadowed `material_b` is refused, never a half-mixture fallback.
- **The homogenize path (R4/D4.3):** a non-background gradient film
  homogenizes to ONE row that is bitwise the direct EMA call at
  `f_mid` over both endpoint spectra, announced by a warning naming
  the mixture and `f_mid`. Never the `inhomogen = false` flag flip -
  a gradient film has no base nk. Background gradient films expand
  WITH the full profile, pinned, silent (A4's background predicate
  gains the `gradient.is_some()` disjunct).

### Notes

- Profile sampling is inclusive (`f_sublayer`), not the midpoint
  reading D3's prose suggests: with midpoint sampling no row sits at an
  endpoint and the gate's bitwise endpoint-row claim is unimplementable;
  the F1.6 table's "normalized depth between two endpoints" and the
  legacy factors' `i / (sub - 1)` convention both point at inclusive.
- `MixRule` gained `Serialize`/`Deserialize`/`PartialEq` (the serde
  representation is the schema-visible surface, pinned by a test);
  `Layer.gradient` is deliberately NOT yet in the serialized state map
  (F1.4 moves it in with the schema bump - no user-constructible
  gradient exists before F1.5's Python surface, so the temporary
  omission cannot lose user data).
- The structure path's provider-less row-count predictions (cold
  `Structure`/`Architect` fallbacks) count a gradient carrier as 1:
  documented at the sites; every production caller takes the
  provider-backed exact path.
- Both fingerprints hold (gradient absent, `None`, untouched paths);
  743 + 5 new gradient twins in `validation/regression/synthesis/
test_gradient_pipeline.py` and the Rust oracle twins in `expansion.rs`
/ `gradient.rs` / `structure.rs`.

## [0.6.35] - F0.3: `ThinLayerPolicy` - clamp up to the minimum instead of removing

"Remove a layer that got too thin" and "set that layer to the minimum
thickness you can actually deposit" are two different operations with two
different purposes; only the first existed.

### Added

- **`PipelineConfig.thin_layer_policy`** (`'remove'` | `'clamp_up_final'` |
  `'clamp_up_always'`; default `'remove'`, which reproduces 0.6.34 bit for
  bit and is what the fingerprints hold):
  - **`ClampUpFinal`** — the search runs exactly as today, elimination and
    all; only the run's final clamp pass sets a surviving sub-minimum film
    to `clamp_min_nm` instead of removing it. A film that survived the
    whole run and happens to sit at 0.8 nm comes out depositable. This is
    the variant the documentation leads with, and the one needle runs
    want.
  - **`ClampUpAlways`** — the floor sets films to `clamp_min_nm` during
    the run too, AND the LM lower bound moves with it (`lb =
    clamp_min_nm`). The coupling is the actual content of the feature:
    clamping up without raising the bound builds an oscillator — LM drives
    a film to 0.5 nm, the clamp puts it back to 2.0, the next solve drives
    it to 0.5 again — and the stagnation detector terminates the run with
    `STAGNATION_OSCILLATION`, a true report of a false condition. With the
    bound, a design whose optimum wants a 0.5 nm film converges to the
    floor in ONE optimization and the merit history is monotone (twin on
    the real machinery: a bare-substrate reflectance target pulls below
    the floor from a 10 nm start; five sweeps, monotone, film parked at
    exactly the floor). Refused at `NeedlePipeline::new` when
    `needles_per_cycle > 0` — with the floor as a hard bound nothing can
    ever be eliminated, so a rejected needle seed parks at the floor
    permanently and the layer count only ever grows; the refusal names
    both settings and points at `ClampUpFinal`.

### Notes

- On a span, clamping up is not a row operation: until F1.6 lands, an
  under-thickness graded span is removed whole with the F0.2 report under
  **every** policy (the enum's doc comment states the deferral; the span
  deferral twin inverts deliberately at F1.6).
- The default is bit-exact: both fingerprints hold with the policy absent
  and with it explicitly `Remove`; the default-equals-`Remove` identity is
  asserted.
- `clamp_min_nm` keeps its name (open decision 14): it does two jobs —
  elimination threshold under `Remove`, manufacturing floor under the
  clamp-up policies — and `ThinLayerPolicy`'s doc comment states both.

## [0.6.34] — F0.2: the floor, the cap and the layer budget become span quantities

The one release in the series licensed to change what a run produces, with
the licence enumerated before the code — seven items, a twin each, and an
eighth difference is a defect. Measured on the shipped 0.6.32 build, not
derived: the floor deleted a pinned 10 nm graded interlayer entirely, the
ceiling never fired on a graded film, the layer budget counted solver rows,
nobody was told about any of it, and the floor also deleted sub-floor
interface slices — losing the interface physics and a nanometre off the
film, both live bugs at 0.6.32.

### Changed

- **`clamp_all` judges spans, not rows** (`DesignStack::clamp_all` now
  returns a `ClampReport`). The comparison reads `D`, the span's
  slice-inclusive total: expansion carves the interface slice out of the
  carrier, so the slice IS part of the authored thickness — excluded from
  scaling, never from measuring (B1).
  - `D < clamp_min_nm` → the whole span is removed in one operation and
    the report names it (`"material (D nm)"`) — licence item 1.
  - `D > clamp_max_nm` on a multi-row span → **refused, not rescaled** —
    licence item 2. The pipeline pre-checks at `NeedlePipeline::new`, so
    a run never hits the refusal mid-flight except across a merge that
    grew a span.
  - A slice row is never a floor or a cap candidate on its own, in every
    branch — licence item 7 (B1, measured: 150.0 → 149.0 nm before, 150.0
    held after, control at 3.0 nm unchanged in both directions).
  - A one-row span is row and span in one: removal and capping behave
    exactly as before — a no-span run is bit-identical, full stop.
- **The layer budget counts spans** (licence item 3): `max_film_layers`
  is a manufacturability limit — one graded film is ONE physical layer, so
  a 57-row graded film no longer terminates the run on the pre-flight of
  cycle 1. `layer_count` in every phase result and `final_layer_count`
  read the physical-layer count (licence item 4).
- **`clamp_all`'s discarded return value becomes a report** (licence
  item 6, B4): `PipelinePhaseResult.clamp_report` and
  `PipelineResult.final_clamp_report` are `Option<ClampReport>` — `None`,
  and therefore absent from the result dicts, when nothing was removed
  and nothing was capped, so a no-clamp run's serialized output is
  byte-identical. The evaluator's clamp (after every thickness
  optimization, dozens of times per cycle) accumulates into the context;
  the pipeline drains it into the phase at record time.
- **`inflate_design` and `round_to_qwot` exclude non-singleton-bulk
  spans** from their selection lists (N3, licence item 5) — a pinned
  profile no longer comes out of an inflate pass with a distorted
  thickness distribution.
- **`remove_thin_layers` takes the span exemption** alongside the existing
  `optimize` filter: a row inside a multi-row span is not a removal
  candidate even when the flag says free.
- **`expand`'s cap refusal** moves nothing: `max_total_thickness_nm`
  stays a row sum (a sum is a sum either way).

### Fixed

- **The floor deletes a graded layer (A2) and the ceiling ignores it** —
  both span quantities now.
- **The floor deletes an interface slice (B1)** — losing the Looyenga
  mixed row AND shortening the carrier by the width the slice was carved
  from. Live at 0.6.32 on any stack with `interface = true` and a sub-floor
  slice; Python-reachable directly through `DesignStack.clamp_all`.

### Note

- The needle-run fingerprint pin's digest strips the `clamp_report` key
  (the one licensed difference), and with it stripped the recorded 0.6.32
  digest HOLDS — the pin proves F0.2's trajectory is bit-identical on its
  problem and that only licensed reporting appeared (three removals that
  0.6.32 performed silently are now named in the phase dicts).
- The Python `DesignStack.clamp_all` keeps its historical
  `(n_removed, n_capped)` tuple contract; a span-level refusal now raises
  `ValueError`.

## [0.6.33] — F0.1: span provenance on `DesignStack` (the bookkeeping half)

### Added

- **`DesignStack` now carries its spans** (`DesignStack::spans`, exposed as
  `spans()` alongside `films()`). A span is the contiguous run of solver
  rows one design layer expands to; `expand` computed this all along and
  every caller dropped it at the door. Nothing observable of any existing
  run reads it yet — this release is pure bookkeeping, and the gate was
  absolute: both fingerprints unmoved, full stop.

- **`Span::is_singleton_bulk()`** — the one predicate later items consult
  (N1): a span is singleton-bulk iff it holds exactly one non-slice row,
  so a plain film (1 row) and an interface-carrying plain film (2 rows,
  one of them the slice) both qualify while graded spans do not. Say
  "singleton-bulk", never "singleton". `Span` also now carries
  `bulk_start` (B9) — the bulk range `expand` computed and dropped — with
  the arithmetic predicate kept as the assertion that the two agree.

- **`DesignStack::span_of_row`**, **`needle_host_refusal`** (R3 — the
  needle-host admissibility rule now also lives at the public
  `insert_needle_seed` door, not only in the scan), and the named
  **`assert_spans_partition`** (R2), called at both constructors and every
  count-changing mutator.

- **The needle-run fingerprint pin**
  (`validation/regression/synthesis/test_needle_pin.py`): there was no
  in-tree test pinning a bit-exact needle run — the anchor lived in a
  scratchpad harness outside the repo (R1). The pin (fixed design, fixed
  targets, fixed seeds, hashed spectra) was committed before the first
  line of `DesignStack` changed, so "fingerprint unmoved" names a test
  that exists in the tree.

### Changed

- **The five row-count mutators maintain the partition**:
  `insert_needle_seed` (the host span splits locally; the seed becomes
  its own singleton; the bottom keeps the host's logical identity),
  `merge_adjacent` (films and spans rebuild in one pass — a merged row
  joins the FIRST span touched, matching "first layer's properties win";
  a span whose rows all merged away is dropped, N2), `remove_film`
  (spans shrink or drop with their rows), `clamp_all` (rows are still
  judged one by one — the span-quantity floor/cap is F0.2's licensed
  change), and `remove_thin_layers` (which maintains the partition
  transitively through `remove_film`).

- **`build_scan_sites` / `run_needle_pass` take the span partition**;
  scan sites are generated only inside singleton-bulk spans. An unbooked
  row list (`&[]`) keeps the flag-only rule the tree always ran — no
  existing stack has a needle-flagged multi-row span, so the scan is
  unchanged (pinned by the N1 twin).

- **`expand`'s emission loop body is extracted verbatim into
  `emit_entry`** (open decision 15, resolved: the extraction lands in
  F0.1 under the unconditional bit-exactness gate, and F1.7's
  `refresh_profiles` will call the same code construction does).
  `from_design` asserts the property the extraction's callers will rely
  on: the design path builds `inv = false` throughout, so emission never
  reaches backwards into previous spans.

### Note

`clamp_all`'s row-wise deletion — including the sub-floor interface-slice
row (B1, measured: 150.0 → 149.0 nm) — is deliberately still present at
this release; the fix is F0.2 licence item 7, and the clamp bookkeeping
test records the state honestly.

## [0.6.32] — README badges, and the MSRV that backs one of them

### Added

- **Six badges on the README**: crates.io version, PyPI version, Rust 1.88+,
  Python 3.12+, LGPL-3.0-or-later, and CI status. Each was fetched and read
  before being committed rather than pasted from a template — both registries
  return 200 and point their `repository`/`Homepage` at this repo, and
  `COPYING.LESSER` exists for the license link to land on.

- **`rust-version = "1.88"`**, declared in `[workspace.package]` and inherited
  by both crates. It was not declared anywhere, so the Rust badge would have
  been a hardcoded number with nothing keeping it honest.

  **Measured, not guessed.** 1.85.0 fails the engine crate with `` `let`
  expressions in this position are unstable`` (let chains, stabilized in
  1.88); 1.88.0 builds the whole workspace. With this declared, an old
  toolchain now says

      error: rustc 1.85.0 is not supported by the following package:
        navette@0.6.32 requires rustc 1.88

  instead of emitting a wall of syntax errors.

- **`Programming Language :: Python :: 3.12` / `3.13` classifiers.** shields.io's
  `pypi/pyversions` badge reads classifiers, *not* `requires-python`, and the
  published metadata carried only a bare `Python :: 3` — so the dynamic badge
  rendered "python | 3" for a package that needs 3.12. Caught by fetching the
  badge and reading its rendered text. The README uses a static `3.12+` badge
  until a release ships the corrected classifiers; it can be swapped back to
  the dynamic one afterwards.

### Note on what the version badges show

Both registries publish **0.5.0** while this tree is at 0.6.32, so the
crates.io and PyPI badges read 0.5.0. That is accurate, not stale markup —
the badges will follow the next release without an edit.

## [0.6.31] — CI had been red for 17 runs and the local battery could not see it

### Fixed

- **`cargo clippy -D warnings` fails on CI**, and has since 0.6.13
  (`perf(smatrix): R5.3`) — **17 consecutive red runs**, every one of them the
  same single finding:

      error: using `chunks_exact` with a constant chunk size
        --> rust/navette/src/smatrix/solver.rs:275

  `Solver::from_flat`'s index cache now uses `as_chunks::<2>()`, which hands
  the closure a `&[f64; 2]` and drops the bounds checks on `c[0]`/`c[1]`. The
  length is verified to be an exact multiple of 2 immediately above, so the
  remainder is empty by construction. Both bit-exactness fingerprints are
  unchanged.

### Why the local battery missed it

`chunks_exact_to_as_chunks` did not exist in the toolchain this work was being
verified against. Local was clippy **0.1.97 / rustc 1.97.1**; CI installs
`dtolnay/rust-toolchain@stable`, which was **1.98.0**. The local run was
genuinely clean — it simply could not see the lint, and nothing in the
verification battery compared the two versions. The whole 0.6.13-0.6.30 range
was verified against a stale linter.

The `rust (test + warnings)` job aborts at the `clippy` step, so everything
after it — the feature-gated builds and, since 0.6.30, the fmt gate — had
**never actually run on CI**. Checked explicitly this release with the 1.98
toolchain: clippy clean on the workspace and on all three feature
combinations, and `cargo fmt --all --check` reports zero diffs, so the gate
landed in 0.6.30 goes green rather than red on its first real execution.

### Added

- **`tools/check_toolchain.py`** — compares the active clippy/rustc against
  the newest toolchain `rustup` has installed, and against the CI workflow's
  declared channel. It prints what CI will use and fails when the local
  linter is older, so "clippy is clean" cannot again mean "clean on a version
  nobody runs". Documented in the README check list.

## [0.6.30] — rustfmt adopted, and the fmt gate goes blocking (R2.3b)

The last open item in the remediation plan's master table. No behaviour change:
both bit-exactness fingerprints are byte-identical across the reformat.

### The decision this was blocked on

R2.3b had been deferred since 0.6.12 pending a style choice: adopt rustfmt
defaults, or write a `rustfmt.toml` encoding the house style. Measuring it
settled the question — **there is no single house style to encode.** Classified
by each file's own indent unit:

| indent | files |
| --- | --- |
| 4 spaces (rustfmt default) | 69 |
| 2 spaces | 22 |

The 2-space cluster is `structure/`, parts of `synthesis/`, `config.rs`,
`sellmeier.rs` and `color/tables.rs`. Adopting the defaults is therefore the
*smaller* change as well as the conventional one; `tab_spaces = 2` would have
re-indented the 69-file majority to match the 22-file minority. **Chosen:
rustfmt defaults, no `rustfmt.toml`.**

### Corrected

- **"981 files differ from rustfmt defaults"**, stated in R2.3b and in the
  `ci.yml` header, was a count of diff *hunks*, not files. Measured before the
  reformat: **1110 hunks across 78 of the 98 `.rs` files**.

### Changed

- **One reformat commit, `cargo fmt --all` and nothing else** — 85 files. It
  is listed in the new **`.git-blame-ignore-revs`**, which GitHub honours
  automatically; activate it in a clone with

      git config blame.ignoreRevsFile .git-blame-ignore-revs

  Verified by blame that lines around recent work attribute to their authoring
  commit, not to the reformat.
- **`cargo fmt --all --check` is now blocking in CI.** With it,
  `continue-on-error` disappears from `ci.yml` entirely — the workflow has no
  advisory steps left.

### What the reformat actually touched

Most hunks in the already-4-space majority are wrapping past 100 columns and
`use` ordering. No `#[rustfmt::skip]` was needed anywhere: the numeric tables
in `color/golden.rs` are one row per line and inside 100 columns, so they came
through intact. Four places where a closing brace shared a line with the
following doc comment (`}    /// ...`) got split — a readability fix found for
free.

Checked explicitly that no doc-comment *prose* moved: every comment line
present before the reformat is present after, modulo leading indentation, the
four brace splits being the only additions.

## [0.6.29] — The Névot-Croce validity limit, written down where it is read

Documentation only; no behaviour change. 0.6.28 established that type-5
roughness produces energy outside its validity range and deliberately did not
gate it. This writes the limit down on every surface a caller might read, with
a usable threshold instead of a hand-wave.

### Why it fails

The reflection factor `exp(-2*kz1*kz2*sigma^2)` decays and is always well
behaved. The transmission factor `exp(+((kz1-kz2)*sigma)^2/2)` **grows**, and
it is the one that fails. Its exponent is driven by the index *contrast*, not
by `kz`: at normal incidence it is `(2*pi*dn*sigma/lambda)^2/2`.

That is the whole story of why a good model gives bad numbers here. Névot-Croce
is an X-ray reflectometry result, where `dn ~ 1e-5` and the transmission factor
is 1.000000... At optical contrast `dn` is five orders of magnitude larger, so
the term that is inert at its origin becomes the dominant error.

Usable budget — sigma injecting at most a fraction `eps` of spurious energy at
one interface:

    sigma_max = sqrt(ln(1+eps)) * lambda / (2*pi*dn)     eps=0.01 -> ~0.0159*lambda/dn

which is **4.7 nm** at dn = 1.35, lambda = 400 nm, and 29 nm at dn = 0.3.
Typical high-contrast visible coatings are already at the edge at single-digit
sigma.

Also worth stating: the physically correct direction is `R + T` slightly
*below* 1, the deficit being diffuse scatter this model does not track.
`R + T > 1` has no mechanism behind it and is the unambiguous signature of
having left the valid regime.

### Changed

- **`optics_core::nevot_croce_factors`** — the validity condition read
  `|kz*sigma| << 1`, which names the wrong quantity. Replaced with the contrast
  form, the closed-form sigma budget, and a measured table.
- **`structure::enums::RoughnessType` (Rust)** had a single line of docs while
  its Python twin had forty. The essentials now live on both.
- **`navette.structure.types.RoughnessType`** gains the threshold formula and a
  dn/lambda table, so "it degrades outside that regime" becomes a number a
  caller can apply.
- **`ScatterMatrix`'s `roughness_types` parameter** said NEVOT_CROCE "conserves
  specular energy (R+T=1)" with no qualifier. It does, to first order in
  sigma^2 — which is the half that matters when it stops being true.

### Fixed

- **`w_function` was documented as "Névot-Croce style interface attenuation"**
  in `_smatrix.pyi`. It is not: type 5 is not handled by `w_function` at all
  (it falls through to the catch-all and returns 1.0). `w_function` is the
  form factor for the analytic graded profiles, types 1-4. Corrected.

### Corrected

- **The figures `R + T = 1.068` and `1025` published in the 0.6.28 entry and in
  R3.5 do not reproduce on the stack those texts name.** Re-measured on
  air / 2.35 (120 nm) / 1.46 (200 nm) / 1.52, type 5 on all three interfaces,
  500-700 nm, `max(R + T)`:

  | sigma [nm] | s, normal | s, 0-89 deg | p, 0-89 deg |
  | --- | --- | --- | --- |
  | 5 | 1.00022 | 1.00022 | 1.01678 |
  | 10 | 1.00336 | 1.00336 | 1.06907 |
  | 20 | 1.04709 | 1.04709 | 1.31084 |
  | 100 | 49.4404 | 637.983 | 2547.9 |

  The conclusion of R3.5 is unchanged and if anything understated — p-pol
  degrades first and fastest, so a normal-incidence sanity check flatters the
  model. Both earlier texts are corrected in place with a note; the 0.6.28
  commit itself is left as published.

## [0.6.28] — A layer's numbers are judged where they are written (R3.5)

Roughness and graded-film parameters were accepted in silence and misbehaved
somewhere else, later. Measured before this release, on the current code:

| Written | What actually happened |
| --- | --- |
| `roughness = -20` | Solved byte-identically to `+20`. Every roughness form factor squares sigma, so a sign slip could never surface. |
| `roughness = NaN` | Every output NaN, with nothing pointing at the cause. |
| `inh_delta = -0.2` | Refinement factor goes negative, `as u32` saturates to 0, sub-layer count collapses to **1**. The grading was silently dropped and the film solved as homogeneous. |
| `inh_delta = 2.5` | Rows of `n = -0.59, k = -0.0125` went into the solver on a 2.35 + 0.05i film. Negative k is optical **gain** — the layer amplifies. |
| `rough_type = 6 / -1 / 99` | Already fail-closed via `try_from_i32` (no change). |

`Structure::validate` and `ScatterMatrix` also disagreed: the former refused a
negative roughness as an Error, the latter accepted it without a word. That is
resolved by keeping the refusal and moving it earlier, not by softening it.

### Added

- **`Layer::property_issues(&self, label)`** — one place that says what a
  layer's numbers may be, returning findings rather than a verdict so each
  caller decides. Finiteness and non-negativity for thickness, roughness and
  interface thickness; the half-open window `[0, 2)` for `inh_delta`; advisory
  warnings for an overhanging interface and for a graded layer with zero
  grading.
- **Every door that builds a layer now uses it**: the PyO3 constructor, the
  four numeric setters, `set_inhomogen` (which is what makes an inert
  `inh_delta` load-bearing), `set_properties`, `from_state`,
  `Structure::validate`, and the synthesis assembler — which builds its films
  from flag dicts and never touches the Python `Layer`, so gating only the
  constructor would have left the design path open. That was the exact mistake
  0.6.26 made and 0.6.27 fixed; it is not repeated here.
- `validation/smoke/test_layer_gate.py` (27 tests) and two engine tests,
  including one pinning every message as ASCII — these cross into Python and
  land on a cp1252 console, where an em dash raises `UnicodeEncodeError`.

### Changed — behaviour, read this

**Constructing an invalid `Layer` now raises `ValueError` instead of deferring
to `validate()`.** `Layer(-5.0, "TiO2")`, `Layer(100.0, "TiO2",
roughness=-20.0)` and friends previously returned an object and failed (or
didn't) later. Three regression tests deliberately built such a layer to prove
`validate()` caught it at solve time; each is re-pointed to assert the refusal
where it now happens, and the solve gate keeps a test of its own using an error
the layer gate cannot see — an unresolvable material, which is a property of
the structure and its provider, not of the layer.

`set_properties` applies to a copy and swaps it in only once the whole batch
passes, so a rejected batch leaves the layer exactly as it was rather than
half-written.

### Refused, not corrected

A negative sigma is not silently turned into `|sigma|`, and `inh_delta = 2.5`
is not clamped to 1.99. In both cases a correction would be indistinguishable
from the bug, and clamping only buys an index that grazes zero instead of
crossing it. Both messages say what the value would have done, not just what
range was expected.

### Deliberately out of scope

The flat-array surface — `ScatterMatrix(roughness_values=...)` and the native
`Solver` — stays permissive, as R3.1 promised. The gate is at layer
construction; raw arrays remain the escape hatch. What that leaves standing,
stated plainly rather than buried: Névot-Croce (type 5) at sigma = 20 nm still
produces **R + T > 1** on a 4-layer stack — energy created, no warning. (The
figures 1.068 and 1025 first written here do not reproduce on the stack named;
the 0.6.29 entry carries re-measured values and supersedes them.) A validity
band needs the wavelength and angle
grid, which a `Layer` does not have, so it cannot live at this gate.
`test_the_flat_array_surface_is_still_permissive` pins this as a decision so it
stays one.

## [0.6.27] — The layer-0 gate moves to where the stacks are actually built (R3.4)

0.6.26 taught `ScatterMatrix` to drop absorption in the incident medium and say
so. It turned out that covered one of five doors, and not the expensive one.
The synthesis and design surface never builds a `ScatterMatrix`: it hands its
own `ambient` straight to the native assembly. Measured before this release, a
`stack_from_layers` call with `ambient=(1.0+0.05j, "air")` produced **zero
warnings** and a stack still carrying `1+0.05j` on layer 0 — which the
optimizer then fitted a design against.

### Fixed

- **`run_design` was discarding its assembly warnings.** `driver.rs` did
  `let (stack, _warnings) = assemble_stack(...)`, and had done since it was
  written. The graded-film homogenization warning never reached anyone on the
  full-run path, making the most expensive path the quietest one. `run_design`
  now returns `(report, stack, warnings)` and they re-emit as Python warnings.
  This is a pre-existing bug in its own right, independent of the ambient work.

### Added

- **The rule now lives once, in the engine**, as
  `optics_core::sanitize_incident_index` — `None` on the common path, one scan,
  no allocation. It sits next to `forward_branch` because it is the half that
  disarms it.
- **Applied at `DesignStack::from_design`**, the single production constructor:
  `stack_from_layers`, `run_needle`, `design_from_config` and the PyO3
  `DesignStack.from_design` all funnel through it. The plain
  `navette._smatrix.DesignStack(...)` constructor reaches `with_films` instead
  and is gated at the PyO3 boundary with the same engine call.

### Why there and not nearer the solver

`from_design` runs once per stack assembly, and `DesignStack::ambient` is
private with no mutator — needle insertion, merging, clamping and thickness
steps cannot put the absorption back. One check covers an entire design run.
The alternative placement, `solver_arrays()`, is called per merit evaluation
(thousands of times per run) and would re-establish something that cannot have
changed. There is no per-evaluation cost.

### Unchanged

- The native `Solver` and `core_engine` stay permissive. R3.1 promised that
  escape hatch and its test still pins it.
- Both engine fingerprints are unmoved (`30d96909…3c6c`, `27db9666…fd85`).

### Corrected

0.6.26's notes claimed the existing branch test proved the sign flip was "dead
code on the supported path". It did not: that test never leaves `r0 < 1`, so it
never enters total internal reflection, which the supported path reaches
routinely. Measured instead, and kept as
`a_sanitized_ambient_never_reaches_the_flip`: ambients 1.0 / 1.52 / 2.35 / 4.0,
every angle from 0 to 89.9 degrees in 0.1 steps, against transparent through
metal-like layers — **28800 cases, evanescent regime included, zero flips**.
The same sweep with `k` on the ambient flips **10680 of 14256**, deciding on an
`Im(cos)` of order 1e-20. The flip is not dead code; it is a tripwire, and the
sanitizer is what disarms it.

### Tests

- `validation/smoke/test_ambient_gate.py` — 11 cases: all five doors, layer 0
  only (an absorbing substrate and film must survive), a transparent ambient
  stays silent, the correction survives every stack mutation, the corrected
  stack is bit-for-bit its transparent twin, and the engine stays permissive.
- `test_both_doors_explain_it_the_same_way` lifts the explanation out of the
  Python warning and requires it verbatim in the Rust one, so the two
  implementations cannot drift apart.
- Four new Rust tests: three on the sanitizer (including `-0.0` left alone — it
  is not absorption and does not trip the flip) and one on the gate inside
  `from_design`.

## [0.6.26] — An absorbing ambient is corrected and announced, not refused (R3.4)

0.6.21 refused to build a `ScatterMatrix` whose incident medium absorbs.
0.6.26 drops the absorption, says so, and solves the stack. **This reverses a
refusal: code that caught `ValueError` from that case will now get a
`UserWarning` and a result.**

The physics behind the refusal has not changed and is still why this item
exists — there is no reflectance to return for an absorbing ambient. What
changed is the remedy: refusing is the wrong trade for a stack whose ambient
absorption is incidental, such as a material table carrying a tiny `k` on air.

### Changed

- **`Im(n[0])` is set to zero, with a `UserWarning`.** The stack is then solved
  with a transparent ambient of index `Re(n[0])`. The warning names how many
  wavelengths were affected, the first index and its value, the largest
  `|Im(n)|`, what was dropped, why the absorbing-ambient problem has no
  reflectance, and what you are getting instead.
- The caller's array is never mutated. A 2-D `layer_indices` that is already
  `complex128` reaches the constructor without being copied, so the correction
  works on a copy — an in-place fix would have silently rewritten the caller's
  own data.
- Only layer 0 is touched. An absorbing substrate or interior layer is ordinary
  physics and is left exactly as it was.

### Why this is a correction, not an approximation

With `Im(n0) = 0` the transverse wavevector is real again, every layer gets the
standard branch, `R = |r|^2` is a true energy ratio and `R + T + A = 1` holds.
The engine solves the transparent-ambient stack *exactly*: the new test builds
the same stack twice, once with `n0 = 1.0 + 0.3j` and once with `n0 = 1.0`, and
requires bit-for-bit equality on all four channels rather than a tolerance. On
a lossless stack `R + T - 1` is now under 1e-12. Before this item it was 1.0096
at `k = 0.1`, 1.44 at `k = 1`, and `Rs = 417` at 10 degrees.

What is genuinely lost — and the warning says so — is the attenuation along the
path *through* the ambient before the light reaches the stack. That factor is
geometry-dependent and a semi-infinite ambient does not define it. Put the
absorbing medium on the substrate side if you need it carried.

There is still no tolerance band: `k = 1e-14` is corrected and warned about like
any other value, because the branch flip it used to cause had no threshold
either.

### Unchanged

- The native `Solver` stays permissive; nothing in the Rust engine moved. The
  branch rule, its doc comment and
  `absorbing_layers_under_a_real_ambient_never_reach_the_flip` all still hold.
  (That last claim was overstated here; corrected in 0.6.27.)
- The engine fingerprint is unmoved at `30d96909…3c6c`.

### Tests

- The two incident-medium rows moved from the reject table to the accept table.
- `test_absorbing_incident_medium_warns_and_explains_itself` — exactly one
  warning, plain ASCII, naming the drop, the transparent ambient, the energy
  argument, the substrate route and the offending value.
- `test_absorbing_ambient_is_solved_as_its_transparent_twin` — bit-for-bit
  against the equivalent transparent stack, plus the energy residual.
- `test_only_the_incident_row_is_touched_in_a_2d_index_array` — a per-wavelength
  grid is corrected on row 0 only, the caller's array is untouched, and the
  interior/substrate absorption survives.
- `garbage_in.py` gained a `warns` verdict distinct from `silent-clean`. The
  distinction is the point: a correction nobody is told about is the failure
  mode this whole item started from. No other row drifted.

## [0.6.25] — The plan closes, and one stub stops lying

Every code item in `docs/remediation_plan.md` is now done or explicitly closed;
the only open row left is R2.3b (rustfmt), deferred by maintainer decision at
0.6.12. Reconciling the master table against the sections it summarizes turned
up two real defects in R1.2, which are fixed here.

### Fixed

- **`solver_energy_conservation` accepts 1-D again.** R1.2 (0.5.2) widened the
  native binding to 2-D by changing its parameter type from
  `PyReadonlyArray1` to `PyReadonlyArray2` — which silently *dropped* the 1-D
  form it had always taken. Nothing in the repo called it that way so nothing
  broke here, but it is an exposed symbol, and PyO3 reports the mismatch as
  `TypeError: argument 'rs': 'ndarray' object is not an instance of 'ndarray'`,
  which tells a caller nothing at all. The arrays are now taken `Dyn` and rank
  1 and 2 are handled alike, with the result returned in the caller's own
  shape. Rank 3+ and shape disagreement raise `ValueError` naming what was
  expected and what arrived.
- **The `_smatrix.pyi` stub for it described a different function.** It said
  "``A = 1 - R - T`` per polarization, as a ``(2, n)`` array". It is
  ``max(|1-Rs-Ts|, |1-Rp-Tp|)`` — one residual per grid point, maxed *over* the
  two polarizations, in the input's shape. Both halves were wrong, in the
  direction that would send a reader indexing into a polarization axis that
  does not exist. (The `ScatterMatrix.energy_conservation` docstring was right
  all along; this was the stub only — the same failure mode as 0.6.21's
  index-layout fix.)

### Added

- Four rows in `validation/smoke/test_energy_conservation.py` pinning the
  native contract: 1-D in → 1-D out; the caller's own shape returned for
  `(3,7)`, `(5,2)`, `(1,4)` and `(6,)`; `ValueError` on mismatched shapes and
  on rank 3; and
  `test_native_is_the_max_over_polarizations_not_a_row_per_pol`, which is
  exactly the claim the stub used to get wrong.

### Changed

- **`docs/remediation_plan.md`: the master table now matches the document.**
  Thirteen rows (R1.1, R1.2, R2.1–R2.5, R3.1–R3.3, R4.1–R4.3) were still
  unstruck although every one of their sections had been marked DONE for
  several releases, so the summary a reader trusts said eleven P0/P1 items were
  outstanding when none were. Struck, each carrying the version that closed it.
  R1.1, R1.2 and R4.4's section headers gained the DONE markers they never got.
  28 of 29 rows are now struck; the one that is not is R2.3b.

## [0.6.24] — Sixteen anonymous slices get their names back (R6.1)

`needle_gradient` took thirty-four parameters, sixteen of them adjacent
`Option<&[f64]>` slices. That is not a signature, it is a memory test: nothing
stopped `targets_rb, weights_rb` from being passed where `targets_tb,
weights_tb` belong, and the result would have been a gradient that is wrong and
entirely plausible. §4.1 called it the codebase's #1 smell.

Nothing about the output changes. This is a pure refactor and it is held to
that: bit-identical on 234 all-channel calls.

### Changed

- **`NeedleDemands<'a>`** — the sixteen slices are now one `Default`-able
  struct with named fields, so a caller writes only what it has and a
  transposed pair is a compile error instead of a plausible number. The seven
  `(target, weight)` pairs become a table with two accessors keyed by a
  `Demand` enum: 16 `load_pair` calls and 14 hand-written accessor closures
  down to 4 and 2.
- **`PointOut`** is fourteen named fields no longer — one
  `[[Option<Vec<f64>>; 2]; N_CH]` indexed by a `Ch` enum. That is what lets the
  ladders become tables: the six-way multiblock ladder is now one loop over
  `mb_table`, and the thirteen seven-line `emit!` blocks are one loop over
  `emit_table`. The `emit!` macro — which existed only because it needed a
  field *identifier* — is gone. Emission order is preserved; callers index the
  returned maps by insertion order as well as by key.
- Free function 34 params → 19, `Solver::needle_gradient` 28 → 13, body
  644 → 543 lines, `if want_*` blocks 31 → 12. Adding an eighth demand type
  used to mean editing six ladders; it now means a variant and a table row.

### Unchanged (deliberately)

- **The Python signature.** All twenty-seven keyword arguments stay exactly as
  they were; the PyO3 wrappers assemble a `NeedleDemands` immediately before
  the call. A refactor whose whole gate is "no numeric drift" is the wrong
  place to also break a public API.
- **The coherent per-point ladder stays a direct-call ladder**, not a
  function-pointer table. It is the innermost loop and R6.3 measured what
  happens there (+6.5 % min, +13 % p10 for a bit-exact change). The demand
  lookup is table-driven; the dispatch is not, and the code says why.

### Gate

- **Bit-identity, measured rather than assumed.** `bitbase.py` covers the
  needle path but only sets `targets_r`/`weights_r`/`grads_r`, so a
  `Demand::Tb` ↔ `Demand::Rb` transposition would have slipped past it. A new
  harness drives all fourteen target/weight slices and both colour buckets
  with *different* random vectors, requests every channel plus the dispersion
  ladder, randomizes incoherent flags so the `Pmb` cascade runs, and sweeps
  s / p / sp × two dispersion channels — 234 calls, key order hashed with the
  values. Builds from before and after the refactor both return
  `27db9666…fd85`; `bitbase.py` is unchanged at `30d96909…3c6c`.
- `fd_step1.py`, `fd_rchannel.py`, `color_merit_check.py`,
  `color_grad_python.py` and the other six review harnesses all exit 0.
- **Perf inside noise, measured on the needle path** (`bench_backside_speed`
  does not run it). Six alternating A/B rounds, 250 reps, three request
  shapes: worst figure +0.88 %, and post-refactor wins the median on two of
  the three shapes. The abort threshold was 2 %.

## [0.6.23] — The solver duplication stays, and is now a test (R6.3)

R6.3 proposed folding `solve_point` and `solve_point_intensity` into one
const-generic body, gated on "abort if `bench_backside_speed` regresses > 2 %
or bit-identity breaks". It was built, it was bit-exact, and it was 3–13 %
slower. So it was aborted — and what §4.2 actually complains about ("deliberate,
but **undocumented** as a maintenance contract") shipped instead.

### Added

- **`intensity_path_matches_full_path_bitwise`** — the contract, enforced. 200
  randomized stacks (3–8 layers, absorbing layers, all six roughness types,
  random incoherent flags, all three coherence modes, s / p / both) drive both
  functions and compare all four intensity channels on `to_bits()`: 7200
  comparisons, no tolerance. A tolerance would pass exactly the reordering this
  is for. Negative-tested by reassociating `2.0 * PI * d_inc / lam` into
  `2.0 * PI / lam * d_inc` in the lean copy alone — the test fails on trial 27
  with a 1-ULP `ts` difference.
- **`intensity_path_leaves_the_complex_fields_unset`** — pins that the lean
  path returns NaN complex amplitudes and a zero coherency channel, so a future
  fill-in cannot make an unread field look meaningful.
- A maintenance contract on `solve_point` saying the copy is deliberate,
  carrying the measurement, and naming the test that holds the two in step.

### Measured (and rejected)

Two `.pyd` builds alternated in one session, eight rounds, 1500 reps each, on
the `ComplexAmps` path:

| | old | const-generic | delta |
|---|---|---|---|
| min | 0.0711 ms | 0.0757 ms | **+6.5 %** |
| p10 | 0.1136 ms | 0.1284 ms | **+13.0 %** |
| median | 0.1431 ms | 0.1477 ms | **+3.2 %** |

Old won on min in 7 of 8 rounds, p10 in 7 of 8, median in 6 of 8. `#[inline]`,
`#[inline(always)]` and no attribute were each built and measured; none closed
the gap. The bodies fold to textually identical code, and the fingerprint was
unchanged — the cost is codegen, not arithmetic.

R6.3's second half (single-pol delegating to dual) was not attempted: the dual
solver has no "one pol off" mode. It walks the interfaces once doing *both*
polarizations, so delegating means computing and discarding the unwanted one,
measured at 1.06–1.10× of total call time — an understatement of the in-solver
cost, since much of that call is fixed overhead.

### Fixed

- **A doc comment that documented an absent function, mid-sentence.**
  `coherent_block.rs` still carried `csqrt_fast`'s five-paragraph rationale
  after the function moved to `optics_core.rs` in the `crates/` → `rust/`
  reshuffle, and the move dropped four lines and glued the remainder onto
  `solve_pol_specialized`'s doc — so the sentence "The earlier naive form was
  chosen to mirror the reference's" simply ended, and `solve_pol_specialized`
  had no summary line of its own. The missing four lines were recovered from
  the pre-port tree; the rationale now sits on `csqrt_fast`, updated to point at
  `forward_branch` (its claim that "both solvers normalize the sign afterwards"
  went stale in 0.6.21), and `solve_pol_specialized` has its summary back.

## [0.6.22] — The color catalog gets its names back (R6.6)

`rust/navette/src/color/` held sixteen files called `func_01.rs` through
`func_16.rs` — the port-task numbering, frozen into the filenames. The only
decoder ring was a doc-comment list in `mod.rs`, which meant every
`use crate::color::func_11::...` in the tree was a lookup.

### Changed

- **16 modules renamed**, 207 references rewritten:
  `func_01`→`xyy`, `func_02`→`lch`, `func_03`→`luv`,
  `func_04`→`oklab_xyz`, `func_05`→`oklab_srgb`, `func_06`→`uvw1964`,
  `func_07`→`ucs1960`, `func_08`→`bradford`, `func_09`→`delta_e_76`,
  `func_10`→`delta_e_94`, `func_11`→`delta_e_cmc`, `func_12`→`din99`,
  `func_13`→`spectral_srgb`, `func_14`→`photometry`, `func_15`→`shapes`,
  `func_16`→`delta_e_2000`. Done with `git mv`, so `git log --follow` still
  walks each file's history.
- **The decoder ring is kept, pointing the other way.** `color/mod.rs` and
  `docs/plans/color/README.md` both carry an old→new migration table, so a
  `func_NN` reference in a commit message, a port record, or §11 of the code
  review still resolves. The dated documents (`docs/plans/color/func_NN_*.md`,
  the review's catalog) keep their original names on purpose.
- Six comments that named *ranges* of modules (`func_01–func_05 core`,
  `func_09..12 & 16: Delta-E metrics`, and four more) were rewritten into prose
  that does not depend on the numbering.
- The stale `// src/func_NN.rs` path header on line 2 of each file is now
  `// src/color/<name>.rs`.

### Unchanged

- No public API moves — the Python bindings call functions, not module paths.
  `check_exposure.py` keys on function names and passed untouched.
- The single-digit `func_0`–`func_5` references in `smatrix/` are a *different*
  numbering (the numba reference implementation's modules, imported by name by
  the parity tests) and were deliberately left alone. `color/mod.rs` now says
  so, so the distinction does not have to be re-derived.
- Bit-exactness fingerprint
  `30d9690992cfa8e9ebfd8a6b03b5d62a7cd511952d3bc50dd3fff9c18ade3c6c`,
  698 pytest, 458 + 22 + 15 cargo, 15 doc tests, clippy clean, ten review
  harnesses green.

## [0.6.21] — An absorbing incident medium is refused, not approximated (R3.4)

`ScatterMatrix(layer_indices=[1.52 + 0.062j, ...])` used to construct and solve
without a word. It returned `Rs = 86` and a negative `Rp`. The remediation plan
offered three ways out — refuse, warn, or renormalize `R` against the
incident-medium Poynting flux — and the third one turned out not to exist.

Against an absorbing ambient the total Poynting flux on the incident side
carries an interference term between the incident and reflected waves that does
not separate into "in" and "out", so there is no reference flux to divide by.
At oblique incidence it is worse than a missing normalization: `kx = k0 n0
sin(theta)` goes complex, the incident wave is inhomogeneous, and a real angle
of incidence does not say *which* inhomogeneous wave it is — the planes of
constant phase and constant amplitude come apart and the angle names only the
first. The input is under-determined. So it is refused.

### Changed

- **`ScatterMatrix` refuses `Im(n[0]) != 0`** — exact test, no tolerance band,
  no opt-out flag. The message names the offending wavelength index and value,
  says what breaks, and points at the two ways out: put the absorbing medium on
  the substrate side (fully supported and unchanged), or drive the native
  `navette._smatrix.Solver`, which R3.1 already established as the documented
  permissive path. Absorbing interior layers and absorbing *substrates* are
  untouched — `t_fwd`/`t_back` normalize by `Re(y)` ratios, so `T` stays an
  energy ratio on the exit side.
- **The `cos(theta)` branch rule now has one home**, `optics_core::forward_branch`,
  replacing five byte-identical inline copies in `coherent_block.rs` (four) and
  `needle_operator::cos_from_nsin`. The rule itself is unchanged, bit for bit;
  what is new is a doc comment recording why it is *not* extended to cover a
  complex incident index, and an ignored `n` parameter that exists to hold that
  explanation where the obvious "improvement" would be written.
- Three further branch sites were deliberately left alone: `cos_inc` in
  `core_engine.rs` and in `spacer_tau` have their flip neutralized by the
  `max(0, Im beta)` clamp that immediately follows, and `solver.rs:2084` is the
  eigenmode `n_eff` path — a different question with a different right answer.

### Fixed

- `src/navette/_smatrix.pyi` (new in 0.6.20) documented `Solver.indices` as
  wav-major. It is layer-major — `layer * n_wavs + wav`, what a C-contiguous
  `(n_layers, n_wavs)` array ravels to; the Rust parameter is literally named
  `indices_layer_major`, and `Solver::from_raw` transposes it into the wav-major
  cache. The `n_stack_cache` buffers elsewhere in that file really are
  wav-major, which is how the file-header generalization went wrong. Both
  corrected.

### Tests

- `optics_core.rs` — four tests. `forward_branch` is pinned against the exact
  inline expression the five sites carried, compared on `to_bits()` across a
  signed-zero sweep of all four quadrants (`-0.0` is not `< 0.0`, and the rule
  depends on that). One pins that `n` is ignored *and* exhibits a case where an
  `Im(n*cos)` rule would genuinely disagree, so the assertion is not vacuous.
  One covers the propagating, evanescent and exactly-critical real ambient.
  The last checks the doc comment's claim that an absorbing layer under a real
  ambient never reaches the flip at all.
- `validation/smoke/test_input_validation.py` — two reject rows (`k = 0.05` and
  `k = 1e-14`; there is no tolerance because 1e-14 was already enough to flip
  the branch), two accept rows (absorbing and strongly absorbing substrates),
  and three tests: the message carries its own justification and is ASCII, only
  row 0 of a 2-D index grid is judged and the offending wavelength index is
  named, and the native `Solver` escape hatch is still open.
- `validation/review/garbage_in.py` — four cases recording the new verdicts.

### Measured

- Normal incidence degrades *smoothly*, which is why it read as a
  normalization bug: `R + T` = 1.000053 at `k = 1e-3`, 1.0096 at `k = 0.1`,
  1.44 at `k = 1`.
- Oblique incidence is erratic. Deterministic and batch-independent, but `Rs`
  at 10 degrees alternates between 0.0024 and 417 as `k` moves over
  1e-16…1e-2, with no monotonicity — because `r0 = nsin * (1/n0)` should be
  exactly `sin(theta)` and instead carries a rounding residue of order 1e-31
  whose sign decides the branch, inverting `Re(cos)` from +0.9848 to -0.9848.
- Two attempts to repair the branch rather than refuse (`Im(n*cos) >= 0`, then
  the same with a 1e-12 relative tie-break) both produced *consistent* nonsense
  — `Rs = 417.6` at every `k` — and both moved the bit-exactness fingerprint.
  Reverted; the core is byte-identical to 0.6.20
  (`30d9690992cfa8e9ebfd8a6b03b5d62a7cd511952d3bc50dd3fff9c18ade3c6c`).

## [0.6.20] — Stubs for the other two engines, and a guard so they stay true (R6.5)

`_materials` was the only stubbed submodule, so IDE completion and any type
checker stopped at materials while the whole S-matrix and weaving surface —
the part anyone actually drives — had nothing.

### Added

- **`src/navette/_smatrix.pyi`** — 53 exported names: the 18 `NREQ_*` request
  bits and their seven builder functions, the eight kernels (`w_function`, the
  three Redheffer products, `solve_coherent_block_fields`, `core_engine`,
  `needle_engine`, `field_profile`, the eigenmode trio), `Solver` with all six
  methods including the 30-argument `needle_gradient`, the merit pair
  (`SimCurves`, `MeritSpec`) and the full synthesis pipeline (`LayerSpec`,
  `DesignStack`, `LmConfig`, `PipelineConfig`, `NeedleCycleConfig`,
  `SmatrixContext`, `NeedlePipeline`, `run_design`, `assemble_design`).
- **`src/navette/_spectralweave.pyi`** — `SpectralDataFrame`,
  `OpticalCollection`, `OpticalWeaver`, `TargetWeaver`, `calculate_merit`,
  including the dunder protocol methods (`__getitem__`, `__contains__`,
  `__len__`) that the registration list does not mention but callers use.
- **`tools/check_pyi_sync.py`**, wired into `ci.yml` beside the exposure and
  CIE guards. Two passes, because they rot differently:
  * **Surface** — parses each `#[pymodule]` body for its `wrap_pyfunction!` /
    `add_class::<>` / `m.add("…")` entries (resolving `#[pyclass(name = …)]`,
    so `PyDesignStack` is checked as `DesignStack`) and diffs against the
    stub's top-level names. Needs no build.
  * **Signatures** — when the extension imports, compares parameter *names*
    from PyO3's `__text_signature__` against the stub's, for every function,
    constructor and method. A renamed keyword argument passes the surface
    check and breaks every caller that used it. 64 of 64 `_smatrix` entries
    compare; nothing is skipped silently.

  Modules without a stub (`_color`, `_interpolate`, `_structure`) are listed
  explicitly in `UNSTUBBED` rather than inferred from which files exist, so a
  *new* unstubbed submodule fails the check instead of joining the gap in
  silence — and adding one of those three stubs fails until the name is
  removed from the list.

### Fixed

- `_materials.pyi` pointed at `crates/navette-py`, a path removed in the Rust
  consolidation, and had no SPDX header. Both corrected.

### Known

- No behaviour change; the fingerprint is unmoved at `30d96909…3c6c`.
- Return types in the new stubs are as precise as the boundary allows:
  `core_engine`, `Solver.solve`, `needle_gradient`, `run_design` and
  `NeedlePipeline.run` return `dict[str, Any]` because their key set is chosen
  at runtime by the request bitmask. The sync guard checks names and
  parameters, not return types — a `TypedDict` per request shape would be the
  next step, and is not one the review asked for.

## [0.6.19] — The doc-line batch, and the one line in it that was a bug (R6.4, part 2)

R6.4's remaining three items. Six of the seven doc lines were exactly that —
lines. The seventh, `extrap='error'`, turned out to describe behaviour nobody
would want, so it is a code change.

### Fixed

- **`extrap="error"` now errors.** It returned an array of NaN. The review
  reported "NaN plus an internal flag"; there is no flag — the NaN was the
  whole signal, and a caller who did not check `isnan` got a plausible-looking
  result with holes in it. `UniInterpolator::evaluate` is now
  `-> Result<Array2<f64>, String>`, refusing once up front with the offending
  index, its value, and the knot range; the binding raises `ValueError`.

  Made fallible rather than added as a parallel `evaluate_checked` because
  there are exactly two callers — the binding and `WeaverMaterialProvider` —
  so the honest signature is cheap, and it leaves no silent-NaN path for Rust
  callers either. The check is one pass before evaluating, not per point, so
  the inner kernels stay branch-free; the other two extrapolation modes cannot
  fail and pay nothing. Out-of-range NaN still exists *inside* the per-point
  kernel, which has nowhere to report to, and can no longer reach a caller.

### Documented

- **`"linear"` extrapolation is two different rules, and neither is scipy's.**
  Found while writing the test for the above. Hermite methods (`pchip`,
  `makima`) leave the knot range along the **endpoint derivative** they already
  computed — C¹ with the curve; secant methods (`linear`, `sprague`,
  `floater_hormann`) leave along the **end chord**. On `y = x³` sampled at
  0..3 that is 52 versus 46 at `x = 4`. scipy's `PchipInterpolator` continues
  the end *cubic* and matches neither. All three are defensible; assuming they
  agree is the trap. Pinned by `linear_extrapolation_means_two_different_lines`.
  Also recorded that derivatives are analytic only for `pchip` — the rest use a
  central difference with `h = 1e-6 · span`.
- **Bradford is row-vector**: `adapted = xyz @ M`, the transpose of the
  textbook column form. Writing the textbook `M @ xyz` in numpy does not raise;
  it returns a plausible, wrong colour.
- **The KK grid stops at 80 eV and clamps past it** — no extrapolation, no
  warning, so a UBF or Cody–Lorentz model queried in the EUV returns the 80 eV
  value. Also: near-resonance ε₁ accuracy is ~1 % (measured 0.07–1.1 % against
  an analytic Lorentz KK pair) and refining the quadrature does not close it,
  because ε₁ is deliberately the FFT-KK result on this grid rather than the
  closed-form Jellison–Modine expression.
- **The colour coverage rule is the opposite of the pointwise one.** A
  `SpectralTarget` whose grid does not line up with the sim grid is
  *interpolated* onto it, so every point still scores; a colour demand whose
  tables do not line up is *narrowed*, so points outside the overlap silently
  score nothing. Both are right for their own mathematics — interpolating an
  illuminant or CMF table would quietly change the colorimetry — but code
  ported from the pointwise side on the assumption that the grids get sorted
  out gets a narrower integral and no warning. Empty overlap still errors.
- **Weaver concurrency: per-frame locking is not per-key transactionality.**
  Two threads writing the *same* key interleave, last-writer-wins per fragment,
  so the reassembled curve can be thread A on the first frames and thread B on
  the rest — a curve nobody wrote. Distinct keys are fully independent, which is
  the intended one-key-per-step usage. No internal key lock exists to reach for.
- **`seed=None` is deliberately asymmetric**: errors on takes the thread RNG
  (the caller asked for randomness and did not pin it, so the run does not
  pretend to be reproducible); errors off takes a seeded(0) RNG nothing draws
  from.
- **`SpectralTarget` now carries the label vocabulary** — `R`/`T`/`A`,
  `RB`/`TB`/`AB` crossed with `s`/`p`/`u`, plus the differential-phase labels
  `PDts`/`PDtp`. Those two encode their own polarization, force raw-radian
  phase normalization, and are **transmission only**: `reference_phase`
  supports `passes = 2` for a reflection round trip but no reflection
  differential label exists, so there is no `PDrs`/`PDrp` to reach for.

### Changed

- **Naming drift resolved in favour of Navette**, in the four places it was
  real: `docs/materials-architecture.md`, `docs/plans/color/README.md` and
  `docs/plans/smatrix/REVIEW.md` are retitled and carry a mapping note, and the
  `src/navette/config/*.py` headers no longer say "Loom". Every remaining
  mention of Loom names the **numba reference implementation** —
  `loom_colorengine.py`, `loom_matrix.py`, `loom_unispline.py`, the `refs/`
  oracles — which really are called that and which the parity tests import by
  name. Renaming those would have broken imports to fix a cosmetic complaint.
  Also fixed a pointer to `LOOM_RUST_ARCHITECTURE.md`, a filename that has
  never existed in this repository.
- **The README has measured numbers instead of adjectives.** A table of the
  §5.2 release-build results with the bench script beside each row, and the two
  caveats kept visible: small grids are dispatch-bound and still behind the
  numba kernel, and `bench_refold` never exercises the needle *insertion* path.

### Added

- Five `UniInterpolator` tests in a module that had none of its own: in-range
  queries never fail in any mode, `error` mode errors where it used to return
  NaN, the other two modes still extend and clamp, an empty query is not an
  error even in `error` mode, and the two-lines extrapolation finding above.

### Known

- Bit-identity held: the differential fingerprint is unchanged at
  `30d96909…3c6c`, all ten review harnesses exit 0, both parity oracles pass.
- `extrap="error"` is a behaviour change for anyone who was relying on the NaN
  return. No caller in this repository was.

## [0.6.18] — Hygiene: the archive, the headers, and two docs that lied (R6.4, part 1)

No behaviour change; the solve fingerprint is unmoved at `30d96909…3c6c`. Four
of R6.4's seven items — the ones that are about the repository rather than
about the physics.

### Removed

- **`attic/`** — 26 retired files, import-isolated and never referenced by live
  code, but responsible for roughly 300 of the code review's 1 528 dead-code
  findings and most of its scaffolding-comment hits. Deleted; git history keeps
  them. `docs/plans/` and every `validation/**/refs/` stay — those are the
  parity oracles, not archive.

### Added

- **SPDX headers** on 213 sources that lacked them (`// SPDX-License-Identifier:
  LGPL-3.0-or-later`, `#` for Python), placed after a shebang and PEP 263 coding
  line where present so nothing moves ahead of a module docstring. 8 files
  already carried one. `COPYING` and `COPYING.LESSER` were already in the tree.

### Documented

- **`validation/README.md` rewritten.** Its tree predated `review/`,
  `parity/synthesis`, `parity/color`, `regression/**` and nine of the twelve
  smoke suites, and still pointed at `crates/navette-materials`, a path that has
  not existed since the Rust consolidation. Now carries the real tree, the
  per-directory test counts (287 smoke / 327 regression / 56 parity / 22
  goldens), the distinction between what `pytest validation` collects and the
  `review/` harnesses it does not, and the numba-thread-pool timing
  contamination that the two whole-engine parity scripts call `cooldown()` for.
  The "the rewrite is behind the kernel it replaced" note is now history with a
  date on it: that was true at 0.64–0.77× and R5.3 fixed it.
- **`docs/plans/bug_fix_plan.md` closed.** BUG-A, BUG-B and BUG-C had unticked
  checkboxes and, directly underneath, their own notes recording the fix — the
  tracker said open while the tree said fixed, which is precisely the audit-trail
  rot the review filed. Ticked against named passing tests
  (`test_asymmetric_interface_position_and_pair`,
  `test_inverted_roughness_follows_plane`,
  `test_partial_inversion_run_edge_clean`, all in
  `validation/regression/structure/test_inversion.py`), with a header saying the
  page is a record rather than a work item and that `expander.py` has been
  `rust/navette/src/structure/expansion.rs` since the port.

## [0.6.17] — The small physics nits, and which of them were real (R6.2)

Four items the review filed as minor. Two were: a missing clamp and a stale
docstring. One was the opposite of what it looked like — the `+ 0.0` the review
called "a no-op … harmless" is load-bearing, and removing it would have swung
Delta_R by 2π at normal incidence. The fourth needed a decision and got one.

### Fixed

- **`DOP_R` is clamped to ≤ 1**, as `DOP_T` already was. The asymmetry was the
  review's finding; the excess on a physical stack is pure round-off
  (`1.0000000000000004`, measured), but a degree of polarization above 1 is not
  a number a caller can do anything with. **This changes output bits** — see
  Changed.
- **`needle_slopes` docstring**: τ̂ is `iβ′(1+r₁₂²)/(1−r₁₂²)`, not `2iβ′/(1−r₁₂²)`.
  The code, the module header and the star-product test always had the right
  form; only the one-line summary was stale.

### Added

- **Sellmeier domain guard.** `n² = 1 + Σ Bᵢλ²/(λ²−Cᵢ)` is only a function
  between its poles: on a pole it is infinite, and inside the first resonance
  n² < 0 and `sqrt` returns NaN — which then travels through every layer of the
  stack and arrives as a NaN spectrum with nothing left pointing at the
  material that produced it. `sellmeier_domain_check` now refuses that grid with
  a message naming the wavelength and where the fit's resonances actually are.
  Both user-facing entry points run it (the PyO3 binding and the spec
  evaluator); the kernels themselves are untouched, so the check costs one scan
  for a non-finite value and the arithmetic that works out *which* pole only
  runs once something is already wrong.
  - The domain edge is not the pole. For BK7 the first resonance is at 77.5 nm,
    but n² stays positive down to ~70.7 nm (the other terms hold it up), so
    50 nm evaluates to a finite, physically meaningless n = 0.47. The guard
    catches non-finite values; it does not pretend to know where a fit stops
    being *trustworthy*. That is what the coefficient set's paper is for.

### Changed

- **Solve fingerprint moved**, for the first time in this series:
  `a99e8383…f154` → `30d96909…3c6c`. The delta is the `DOP_R` clamp and nothing
  else — verified by recomputing the unclamped ratio from the `S0_R`…`S3_R`
  channels (which stay unclamped) and asserting the emitted channel equals
  `min(raw, 1)` at every one of the fingerprint's 170 points. 20 of those
  exceeded 1: 18 by one ulp, and 2 by a lot (6.1 and 51.1) — those two come
  from the fingerprint harness randomizing the *incident* medium to an
  absorbing index, where `R = |r|²` is not an energy ratio at all and `Rs`
  reaches 86 and `Rp` goes negative. With a real ambient the same stack's
  worst `DOP_R` is 0.99997. Noted as an input-validation gap, not fixed here.

### Documented

- The `+ 0.0` on `s2r`/`s3r`/`s2t`/`s3t` is a negative-zero flush, not dead
  weight: `-2.0 * 0.0` is `-0.0`, and `atan2(-0.0, negative)` answers −π where
  `atan2(+0.0, negative)` answers +π. An isotropic stack at normal incidence
  has `cross_r` exactly real, so this is the ordinary case, not a corner. The
  numba reference carries the same flush at the same four places, with a
  comment; the Rust side now has one too, plus a test pinning Delta_R = +π.

## [0.6.16] — The two argmin backends, and what they are worth (R4.4c-argmin)

The last row left open under R4.4: `OptimizerBackend::ArgminGaussNewton` and
`ArgminTrustRegion`, behind a new `opt-argmin` cargo feature. Both are
reference points rather than candidates, and the interesting part of this
release is the measurement, not the plumbing.

### Added

- **`opt-argmin`** (off by default): `LmConfig(optimizer="argmin_gauss_newton")`
  and `"argmin_trust_region"`. Both unbounded, so both run on the same
  `IntervalMap` reparametrization `minpack_lm` uses — the map and the chain
  rule on the analytic Jacobian are shared, not re-derived. The feature reuses
  `nalgebra` through `argmin-math`'s `nalgebra_v0_34` backend, so enabling it
  alongside `opt-minpack-lm` pulls one linear-algebra crate, not two.
- **`OptimizerBackend::runs_to_max_iterations()`** — true for the argmin trust
  region alone. argmin's `TrustRegion::terminate` returns `NotTerminated`
  unconditionally: it has no convergence test, so `max_iterations` is its only
  exit and `MaxIterations` is its normal outcome, not a failure. A caller
  reading the termination reason needs that declared rather than inferred.
- **Part E of `validation/review/lm_check.py`**, skipped unless the wheel was
  built with the feature: the two backends on the same three refold starts the
  rest of the harness uses, against scipy.
- Ten cargo tests across both backends (interior optimum, box containment,
  cost convention, error propagation, analytic-Jacobian precedence, and the
  two failure modes below), plus name round-tripping and feature-hint coverage
  now driven from `ALL_BACKENDS` rather than a hand-kept list.

### Known — measured, on `validation/review/lm_check.py`'s refold starts

| backend | two films | two films, far | three films |
|---|---|---|---|
| `trf` | 8 evals | 21 evals | 32 evals |
| `argmin_trust_region` | 218 evals | 224 evals | 393 evals |
| `argmin_gauss_newton` | refuses | refuses | 201 evals, cost 2690.8 vs 2670.0 |

- **`argmin_trust_region` is correct and expensive.** It reaches the same
  optimum as scipy on all three starts (agreement 1e-13 to 1e-9) and pays
  7–27× `trf`'s evaluation count to do it, because it cannot stop early.
- **`argmin_gauss_newton` is not usable on this problem class.** The undamped
  step is `(JᵀJ)⁻¹Jᵀr`, so a singular `JᵀJ` makes it undefined — and a film
  driven toward zero thickness, or a bound the interior map has flattened, is
  enough to produce one. It refuses two of the three starts; on the third it
  runs to the iteration cap without converging. The refusal now carries a
  message that says which of those happened and names a damped backend to use
  instead, rather than surfacing argmin's bare "Non-invertible matrix".
- Neither backend changes anything on a default build: `available_optimizers()`
  still returns `["builtin", "trf"]`, and naming one the wheel lacks is still a
  `ValueError` with the rebuild command.

## [0.6.15] — The derive pass is no longer the serial half (R5.2)

R5.3's scaling bench made the gate for this item concrete: after the
construction fix, a rigorous twelve-channel request at 20 000 points spent
~1.1 ms of its 2.4 ms in one serial loop — the per-point derive, where every
requested channel's algebra and transcendentals (`sqrt`, `atan`, `atan2`,
`arg`) are evaluated. The solve itself was already on the rayon pool; this pass
was not.

### Changed

- **The derive pass splits across the rayon pool.** The point range and *every
  live destination buffer* are halved together, down to a 512-point leaf, so
  each output index is still written exactly once, by one task, with the
  expression it had when the pass was serial. There is no reduction and no
  shared accumulator: bit-identity is structural, not a tolerance.
- The 45 channel buffers and their halving are generated from one list, so a
  new channel cannot be added to the pass and forgotten by the split.

Measured with `bench_core_engine_scaling.py`, 6 layers, 1 angle, median of 25
(0.6.14 → 0.6.15):

| points | photometry | | rigorous | |
|---:|---:|---|---:|---|
| 500 | 0.183 → 0.177 ms | 1.03× | 0.231 → 0.231 ms | — |
| 2 000 | 0.364 → 0.344 ms | 1.06× | 0.440 → 0.401 ms | 1.10× |
| 5 000 | 0.592 → 0.540 ms | 1.10× | 0.769 → 0.666 ms | 1.15× |
| 20 000 | 1.539 → 1.499 ms | 1.03× | 2.387 → 1.684 ms | 1.42× |
| 60 000 | 4.229 → 3.700 ms | 1.14× | 6.977 → 4.022 ms | 1.73× |

The gain tracks how many channels the request asks for, which is what the pass
costs: a four-channel photometry request has almost nothing to divide, a
twelve-channel ellipsometric one has 1.7× at 60 000 points. Against the numba
kernel the rigorous request now runs at ~1.1–1.8× at 20 000–60 000 points,
where 0.6.13 had it at ~0.4–0.5×.

### Added

- Four solver tests: the split pass is bit-identical to a single serial call on
  a grid that is not a multiple of the leaf; the split reaches every index
  (buffers start as NaN, so a range nobody visits is caught directly rather
  than by comparing two equally-gapped runs); a grid below the leaf still
  derives; and absent channels — 41 of the 45 on a photometry request — survive
  the halving at every level.

### Known

- **The leaf size is a scheduling knob, not a numerical one.** Swept at 256 /
  512 / 2048 / 8192: the grid sizes that matter land within ~3% of each other,
  the one real effect being that a leaf above the grid means no split at all
  (2 000 points cost 0.428 ms at 2048 and 0.386 ms at 512). 512 is the smallest
  leaf that still buys something measurable.
- Small grids remain dispatch-bound and behind the numba kernel; see 0.6.14.
- Bit-identity held: the differential fingerprint is unchanged, all ten review
  harnesses exit 0, both parity oracles pass.

## [0.6.14] — `core_engine` was paying for its own construction (R5.3)

R2.4a timed the Rust engine at 0.5–0.65× the numba kernel it replaced and the
plan put the blame on the emit path — a `Vec<OpticalState>` walked once per
requested channel on the way out. The profile says otherwise. At 20 000 points
the emit is 0.01–0.25 ms of a 4.3–4.8 ms call; **`Solver::new` is 2.2–2.5 ms of
it**, and 1.77–2.24 ms of *that* is building the index cache as a
`Vec<Vec<Complex64>>` — roughly 40 000 heap allocations per call, on a path
whose callers construct a fresh solver for every evaluation.

### Changed

- **The index cache is one flat wav-major `Vec<Complex64>`** (stride
  `n_layers`) instead of a `Vec` of per-wavelength `Vec`s, and so is its
  reciprocal. Same numbers, same order, two allocations instead of 40 000.
- **`Solver::from_wav_major_flat`** — a second constructor that takes the
  interleaved `[re, im]` buffer the Python boundary already holds. `core_engine`
  was transposing that into layer-major so `Solver::new` could transpose it
  straight back; both passes are gone. `new` keeps its layer-major signature and
  its behaviour.
- **`flat_cache` is built on first use** (`OnceLock`), not in the constructor.
  Only `needle_gradient` reads it; every ordinary solve was paying for it.
- **The emit takes its `Solution` by value.** `PyArray::from_vec` moves, so
  borrowing meant cloning every channel — one extra full copy of the result.
- **Both parity oracles and the new scaling bench pause between engines.** See
  Known: the 0.5–0.65× figure was partly a measurement artefact.

### Added

- `validation/benches/smatrix/bench_core_engine_scaling.py` — the ratio across
  the grid (500 → 60 000 points) at two request breadths, which is what the
  single-size speed line in the parity scripts could not show.
- Four solver tests (the two constructors agree bit-for-bit on a dispersive
  stack; a mis-sized flat cache is refused; the flat constructor validates what
  `new` validates; the serial and parallel reciprocal paths agree across the
  threshold) and two Python tests that a wrong-length index cache raises instead
  of panicking through the binding.

Measured with `bench_core_engine_scaling.py`, 6 layers, 1 angle, median of 25:

| points | photometry before → after | | rigorous before → after | |
|---:|---:|---|---:|---|
| 500 | 0.254 → 0.183 ms | 1.39× | 0.324 → 0.230 ms | 1.41× |
| 2 000 | 0.598 → 0.364 ms | 1.64× | 0.718 → 0.440 ms | 1.63× |
| 5 000 | 1.233 → 0.592 ms | 2.08× | 1.474 → 0.769 ms | 1.92× |
| 20 000 | 4.149 → 1.539 ms | 2.70× | 5.023 → 2.387 ms | 2.10× |
| 60 000 | 12.63 → 4.229 ms | 2.99× | 15.74 → 6.977 ms | 2.26× |

Against the numba kernel the same bench now reports roughly 0.3–0.5× at 500
points, 0.6–0.8× at 5 000, and 0.8–1.2× at 20 000–60 000 — the rewrite passes
the kernel it replaced on large grids and is still behind on small ones, where
~150 µs of rayon dispatch is most of the call.

### Known

- **numba's default threading layer contaminated the original ratio.** Its
  workers spin-wait after a call returns, so timing the two engines back to back
  charges whichever runs second. The same Rust call at 20 000 points measures
  1.52 ms on its own, 2.00 ms immediately after a block of numba calls, and
  1.58 ms after a one-second pause. Forcing numba onto the non-spinning
  `workqueue` layer also removes the effect, but makes numba itself ~70% slower
  here, so both benches take the pause instead and leave each engine in its
  native configuration.
- **Small grids are still dispatch-bound.** Below ~5 000 points the parallel map
  costs more to start than the work it distributes, but serializing them is
  worse: a `SOLVE_PAR_THRESHOLD` experiment made 2 000 points 5× slower
  (0.35 → 1.75 ms), because serial is ~875 ns/point against ~80 ns/point
  effective under rayon. Not taken.
- **What is left is the derive loop.** On a rigorous 20 000-point request it is
  ~1.1 ms of a 2.4 ms call, and it is serial. That is plan item R5.2, whose gate
  ("only if the scaling bench shows the serial fraction matters") this bench now
  satisfies.
- Bit-identity held throughout: the differential fingerprint is unchanged, and
  all ten review harnesses and both parity oracles pass.

## [0.6.13] — What the batch unweave was actually spending (R5.1)

`unweave_collection` ran at 0.22–0.48× the Python reference and the plan blamed
per-key hash-map rebuilds, `String` clones and a whole-collection clone. None of
those were there any more. What it spends is a `memcpy`, and the reason the
reference does not is that the reference **does not own its data**.

### Changed

- **The batch unweave copies only the span the plan reaches**, not the whole
  input curve, and shares that one copy across a key's fragments. Frames rarely
  tile the entire curve handed to them; the part no fragment reads was being
  copied per key.
- **One collection-wide lock per key instead of one per fragment.** The key→
  frames map was taking its write lock once per (key, frame) pair — 16 384 times
  for a 512-key call, all of it on one lock. `OpticalCollection::
  map_frames_to_key` takes it once per key; the singular form now delegates to
  it.
- **`SpectralDataFrame::set_data` hashes the key once**, not twice
  (`contains_key` then `insert` → `insert` and read what it returns).
- **Above 8192 fragments the per-key loop runs on the rayon pool.** Gated on the
  fragment count rather than the byte volume, because the bytes are the part
  that does *not* parallelise — see Known.

Measured on the bench grid, Rust-side `unweave_collection`, median of 7:

| grid × frames × keys | before | after | |
|---|---:|---:|---|
| 1 000 × 32 × 128 | 0.296 ms | 0.234 ms | 1.27× |
| 10 000 × 64 × 256 | 2.585 ms | 1.346 ms | 1.92× |
| 10 000 × 256 × 64 | 1.024 ms | 0.446 ms | 2.29× |
| 100 000 × 512 × 32 | 2.957 ms | 2.057 ms | 1.44× |
| 100 000 × 32 × 512 | 41.29 ms | 35.92 ms | 1.15× |
| 100 000 × 16 × 32 | 1.848 ms | 1.8 ms | — |
| 500 000 × 32 × 16 | 3.711 ms | 3.7 ms | — |

### Fixed

- **A curve that is not one value per wavelength no longer panics in a batch.**
  The distribution plan addresses the source by position, so a short curve was
  indexed straight out of its buffer — a Rust panic surfacing as
  `pyo3_runtime.PanicException`, which does not even derive from `Exception`, so
  an `except Exception:` around the call did not catch it. The Python wrapper
  checked this for a single `unweave` and never for `unweave_batch`. It is a
  `ValueError` naming the key now, on both paths.
- **A rejected batch writes nothing.** Everything that depends only on the plan
  and the grid is checked before the first key is distributed, so a grid that
  misses part of a frame no longer rejects the call after an arbitrary prefix of
  it has already been applied — which, once the loop can run on a pool, would
  have meant a *scheduling-dependent* prefix.

### Added

- **`batch_equals_per_key`** in the spectral bench's correctness section, on both
  engines: a 40-frame × 256-key batch must write exactly what 256 single
  `unweave` calls write. The batch path resolves one plan, shares one
  materialised source per key and may run pooled; nothing in the per-op timings
  would have noticed the two drifting.
- 13 tests for `opticalweaver`, which had none of its own (it was covered only
  from Python), and four for the wrapper's batch guards.

### Known

- **`unweave_collection` is still ~0.28× the Python reference at 512-key scale,
  and that gap is not an optimization target.** The Python reference stores
  *views into the caller's numpy buffer*; the Rust engine owns its fragments.
  Mutating the caller's array after the call changes the stored curve on the
  Python engine and does not on the Rust one. Owning it means copying it: 512
  keys × 100 000 points is 400 MB in and 400 MB out, and 36 ms of that is
  22 GB/s of DRAM traffic — this machine's ceiling, not a code path. Measured
  1 / 4 / 32 rayon threads: 39.7 / 35.1 / 36.5 ms. The plan's stated goal for
  this item, "≥ 1× vs numba at 512-key scale", is reachable only by adopting the
  aliasing, which would make a fragment change under a caller who edited their
  own array afterwards. Recorded in `docs/remediation_plan.md`; not done.

## [0.6.12] — The method the docs have been naming (R4.6)

`thick_opt.rs` has said since the rewrite that it replaces
`scipy.optimize.least_squares(method="trf")`. What it implements is a
Levenberg-Marquardt that keeps `lb ≤ x ≤ ub` by vetoing and clamping the step
it already solved — bounds applied after the fact, to a step computed as
though they were not there. `synthesis::trf` is the method that sentence was
naming: trust-region reflective (Branch-Coleman-Li), where the bounds enter
the subproblem.

Hand-rolled, no new dependency, no cargo feature — so `optimizer="trf"` works
on a standard wheel. The built-in stays the default.

### Added

- **`LmConfig(optimizer="trf")`** — trust-region reflective. The trust region
  is reshaped every iteration by the Coleman-Li scaling `D = diag(√v)`, `v`
  being the distance to the bound the anti-gradient points at, and each step
  is the best of three candidates: the trust-region step cut back to the first
  bound it hits, that step **reflected** off the bound, and the constrained
  Cauchy step.
- **`synthesis::trf::trust_region_reflective`** — same signature and same
  `LmResult` as the built-in LM, so the two are interchangeable behind
  `run_optimizer`. The subproblem is Moré's, solved from the augmented QR
  R4.4b already built rather than from an SVD; the Newton recurrence on the
  secular equation is identical term for term.
- **`lm_check.py` part D** — the scipy comparison this repository could not
  make before. Parts A–C compare two *different* algorithms and can only ask
  for the same optimum; D is the same algorithm on both sides, so it asks for
  the same answer: costs to 1e-9 relative, thicknesses to 1e-4 nm, including
  on a merit whose optimum sits **on** the clamp. This also closes R4.4d's
  deferred "parametrize over every enabled backend" row.
- **`OptimizerBackend::feature()` is public** — the Python binding was
  carrying its own copy of the backend→cargo-feature mapping for its rebuild
  hint, which is the kind of pair that drifts.

### Known

- **TRF never returns a thickness exactly on a bound.** Its iterates must stay
  strictly interior — that is what keeps the Coleman-Li scaling
  differentiable — so where the built-in returns `50.0` it returns
  `49.999999999986834`, at an identical merit. The removal sweep compares
  against `clamp_min`, so a film driven to the bound is still removed, but a
  caller testing `x == ub` will be disappointed. This is the reason TRF is not
  made the default here.
- **`lambda_init`, `lambda_up`, `lambda_down`, `damping` and
  `gtol_scale_invariant` do nothing on this backend.** They are
  Levenberg-Marquardt settings; the trust-region radius plays their role and
  is not user-settable, exactly as in scipy.
- **The plan's algorithm sketch for this item described a different method**
  (a frozen MINPACK column-norm `D`, and "2-D subspace minimization" — which
  is scipy's sparse `tr_solver`, not its reflections). Both corrected in
  `docs/remediation_plan.md`; following either would have cost the scipy
  oracle that makes the item checkable.

## [0.6.11] — The two whole-engine parity oracles are back (R2.4a)

`parity/smatrix/test_core_engine_photometry_only.py` and
`test_core_engine_rigorous_ellipsometry.py` were the only tests that compared
the *entire* engine — inputs in, thirteen channels out — against the numba
kernel it replaced. Since R2.4 (0.5.8) they had been skipping, because they
loaded a `navette_matrix` extension that does not exist in this repository.
They now run against `navette._smatrix.core_engine`, and the parity layer has
no NEEDS PORT rows left: 13/13 channels, agreement ~1e-14 to ~1e-16.

### Changed

- **Both scripts ported to the request-driven `core_engine`.** The legacy
  kernels took `calc_s, calc_p` and a `debug_flag` and always returned a
  fixed-width tuple; the successor takes a request mask and emits only what
  was asked for. Every element of both tuples is still compared.
- **The coherence mode is `FRONT_BLOCK`, and the mapping is asserted.** The
  legacy ellipsometry kernel took reflection S₂/S₃ from the first coherent
  block and transmission S₂/S₃ from a Mueller cross-term product across
  blocks — exactly `core_engine.rs`'s `track_cross_channel == false` path.
  Half the cases flag an incoherent layer so the modes actually differ, and
  the script requires `COHERENCY_MATRIX` to *disagree* with the reference:
  without that, "mode A is the legacy treatment" would be an untested claim
  that happens to hold because nothing separates the modes.
- **A skipped channel is asserted absent, not compared against zeros.** The
  legacy kernel zero-filled the polarization it was told to skip. Comparing
  the missing channel against those zeros was tried deliberately: it passes,
  and would keep passing if the mask were ignored entirely. The port asserts
  the key is missing, and separately that `Rs` from an s-only request is
  bit-identical to `Rs` from the full one.
- **The 13th value moved rather than disappeared.** `conservation_err` is no
  longer an engine channel; it is `solver_energy_conservation` (R1.2), and
  the port reconstructs and compares it from the four intensities.

### Known

- **The rewrite is slower than the kernel it replaced** on this workload —
  0.5–0.65× across 500 to 60 000 grid points — and roughly half the cost
  looks to be outside the physics. Recorded as **R5.3 (P2)** rather than
  fixed here: a performance change does not belong in a parity port.
- **`Delta` is compared on the circle as a guard, not a finding.** The last
  ellipsometry case is built to straddle ±π, and the script asserts it really
  does — but both engines pick the same side everywhere, so a plain
  difference would pass this file today. Removing the circle handling was run
  as a deliberate break and was *not* caught; it is recorded in the plan that
  way rather than dressed up.

## [0.6.10] — Which solver runs is now a choice, and a named one (R4.4c, work item B)

The thickness optimizer's residual system was always solver-agnostic — a
closure over `x` and a `JacobianSource` over `x` — but only one solver was
ever reachable, and the plan's claim of scipy parity had no second
implementation to be parity *with* (review §3.6, §18.2). `synthesis::optimizer`
is now the single place that decides who gets handed that pair, and the
argmin ecosystem's reference LM is a `--features` flag away.

**Nothing changes for a standard build.** The default backend is the built-in
LM, the optional crates are default-off, and the whole cargo and Python suites
pass unchanged.

### Added

- **`LmConfig(optimizer=...)`** — `"builtin"` (default) or `"minpack_lm"`. The
  latter is the `levenberg-marquardt` crate (rust-cv, MINPACK `lmdif`-derived,
  MIT), available when the wheel was built `--features opt-minpack-lm`.
- **`navette._smatrix.available_optimizers()`** — what *this* wheel can run.
  Backends are cargo features, so the answer is a property of the build, and
  a caller should be able to ask rather than discover it by failing.
- **`synthesis::optimizer`** — `OptimizerBackend`, `OptimizerResult`,
  `run_optimizer`, and `IntervalMap`. A new solver is an arm here, not a
  second copy of the call site; R4.6's TRF backend takes the same seam.
- **`optimize_thicknesses_report`'s dict gained `"backend"`** — which solver
  produced the answer, beside which Jacobian path it took.

### Changed

- **`LmConfig` gained `backend`.** The plan's separate `OptimizerConfig`
  "mapping 1:1 onto `LmConfig`" would have been `LmConfig` plus one field, and
  two structs that must agree field for field are a synchronization bug
  waiting to be written. Every other knob already carries over — ftol, xtol,
  gtol and max_evals are MINPACK's own names.

### Notes

- **Bounds are not a shared contract, and this is the reason the built-in
  stays the default.** Navette's LM vetoes and clamps the step it just solved,
  so an optimum may sit *exactly* on a bound — a film driven to zero thickness
  is a real answer the synthesis loop then removes. Every reference
  implementation in the ecosystem is unbounded, so those backends run on
  `x = mid + half·tanh(u)` and their optima are *strictly inside* the box,
  with the gradient vanishing as a bound is approached. A contract
  difference, not a bug.
- **A missing backend is refused, never substituted.** Naming one the wheel
  was not built with raises `ValueError` carrying the rebuild command. A
  silent fall back to the built-in would make every "compared against the
  reference LM" claim a comparison with ourselves.
- **The plan's ndarray/nalgebra friction does not arise.** argmin-math 0.5
  ships a `nalgebra_0_34` backend and `levenberg-marquardt` 0.15 is built on
  nalgebra 0.34, so the two share one linear-algebra crate and navette's
  ndarray 0.17 pin is untouched — no adapter, no per-iteration conversion, no
  version bump. (The plan's alternative, argmin-math's default `Vec<f64>`
  backend, would not have worked: it has no `ArgminInv`, which GaussNewton
  requires.)
- **`stepbound` is not exposed.** It is a MINPACK-only knob; putting it on the
  shared config would add a setting that does nothing on the default backend.
  The crate's default stands until a real workload argues otherwise.
- Eight deliberate breaks, all caught: the interval map made linear, the chain
  rule dropped, the chain rule applied along rows instead of columns, the
  result left in `u`-space, MINPACK's ½‖r‖² not converted to Navette's ‖r‖²,
  an unavailable backend falling back to the built-in, differenced Jacobians
  counted as analytic, and the `atanh` clamp removed.

## [0.6.9] — The thickness optimizer stops guessing its own Jacobian (R4.5, increment ii)

`build_jacobian` was central differences: **2n full residual evaluations per
LM iteration**, each one a complete solver sweep over the angle × wavelength
grid, and each column carrying the difference noise that matters most exactly
where it hurts — small steps near an optimum (review §3.6). The Jacobian is
now assembled analytically:

```text
    J[i,k] = Σ_terms  ∂r_i/∂(curve value) · ∂(curve value)/∂d_k
```

with the left factor from 0.6.8's `MeritSpec::curve_sensitivity` and the right
from the same solver sweep that produces the curves. Analytic is the default;
`jacobian="fd"` restores the differences.

**No pinned optimum moved** — the whole cargo and Python suites pass unchanged,
so the §8 golden protocol had nothing to triage — and the S-matrix engine's
bit-exactness fingerprint is unchanged.

### Changed

- **`LmConfig.jacobian`** (`"analytic"` default, `"fd"`) chooses the path. The
  analytic one applies when every residual row is covered by the sensitivity
  chain; a spec with a phase target or a color demand declines *per run* and
  is differenced, so the setting is a preference, not a promise.
- **The thickness derivative is the needle operator with the needle material
  set to the host's own index.** Growing film *j* by δ is inserting a slab of
  *n_j* inside layer *j*: r₁₂ vanishes, ρ̂ with it, and τ̂ = iβ_j is the bare
  propagation slope. The existing dual-number composition `U ⊗ N ⊗ L` then
  differentiates the whole block through the same Redheffer star product the
  forward solver uses — every multiple-reflection path included, no
  hand-expanded algebra, and no way for the derivative to drift from the value
  it differentiates. Intensities follow the solver's own `finalize`: R = |r_f|²
  and T = |t_back|²·f_back, with the boundary-admittance factor constant under
  a film thickness and riding outside the derivative.

### Added

- **`SmatrixContext::simulate_with_deposits`** — the simulate that also reads
  off ∂(Rs, Rp, Ts, Tp)/∂(thickness) per grid point. The curves come back
  bit-identical to `simulate`'s, which is a test, not a hope.
- **`synthesis::jacobian`** — `CurveDeposits` and `assemble_jacobian`, the seam
  where the two halves meet. J is materialized m×n with per-destination writes
  and a fixed accumulation order per row: no cross-point reduction, nothing
  whose order depends on scheduling (§13).
- **`SmatrixContext.optimize_thicknesses_report(stack)`** → `(merit, report)`,
  with `iterations`, `evals`, `cost`, `termination`, `gain_ratio` and
  **`analytic_jacobians`** — how many of the run's Jacobians came from the
  analytic chain. Which path ran is now something a caller can *ask*, rather
  than infer from a timing; the bench and the fallback tests assert on it.
- **`levenberg_marquardt_with`** and the `JacobianSource` trait. A source may
  supply a Jacobian, decline (`Ok(None)` → difference this iteration), or fail
  (`Err` → abort). Declining and failing are deliberately different: a source
  that believed it had a Jacobian and was wrong should not hide behind a
  slower run.
- **Cargo tests** — the per-element analytic-vs-FD cross-check the plan names
  as required (deposits against `simulate` differenced in thickness space, and
  the assembled J against `residuals` differenced the same way: worst relative
  deviation ~1e-9), the coverage test at the seam the optimizer uses, the
  absorption two-channel row, and five driver-level tests of the hook.
- **`validation/smoke/test_analytic_jacobian.py`** — the end-to-end half, CI's
  subset: the two modes reach the same optimum, the report says which ran, a
  phase demand falls back and says so.
- **`bench_refold.py`** grew a Jacobian section. At ten films, iteration for
  iteration: **1.7× wall-clock, 16× fewer residual evaluations**.

### Notes

- **The speedup is not 2n, and the plan's estimate needs correcting.** Building
  the deposits costs one `StackFields` decomposition per point *and*
  polarization, so an analytic iteration is roughly three solver sweeps rather
  than one — against 2n + 1 for the differenced one. That is 1.7× at n = 10,
  not 20×. What does scale as promised is the residual evaluation count (16×
  fewer at n = 10, and it keeps growing with the stack) and, more importantly,
  the noise: the columns are exact.
- **`validation/review/lm_check.py`'s A2 stationarity check was testing the
  iteration cap.** At the default 200 iterations, both solvers are still
  crawling on the far-start cases — a run that stopped on `MaxIterations` has
  not claimed a stationary point, so asserting one of it says nothing about
  the optimizer. The harness now gives those cases room to converge and checks
  the termination criterion explicitly. (The analytic path is what exposed
  this: it takes a different path through the same flat valley and was a few
  parts in 10⁵ behind at iteration 200. With room, both converge to the same
  cost and both are fixed points.)
- Eight deliberate breaks, all caught: the T deposit losing its flux factor,
  differentiating the forward amplitude instead of the backward one, the film
  index not offset past the ambient, the factor of two in d|z|²/dd, the needle
  taking the ambient's index instead of the host's, a transposed J, a row
  keeping only its last term, and the coverage gate removed.

## [0.6.8] — Every residual row can say what it depends on (R4.5, increment i)

The merit half of the analytic Jacobian. `MeritSpec` can now report, for each
residual row it produces, which simulated curve values that row reads and with
what derivative — the factor `J[i,k] = Σ_terms ∂r_i/∂(curve value) ·
∂(curve value)/∂θ_k` needs on the left. The right-hand factor, the per-point
deposits from the solver sweep, is increment ii; nothing in the optimizer uses
this yet, so no synthesis result moves.

Splitting R4.5 this way is deliberate: the merit half is verifiable on its own
against a finite difference of `residuals()` **in curve space**, with no solver
in the loop, so when the two halves are joined a disagreement has only one
place left to be.

### Added

- **`MeritSpec::curve_sensitivity(&sim) -> MeritSensitivity`** — `rows[i]` is
  the list of `CurveTerm { curve, angle_row, wavelength, d_residual }` for
  residual `i`, in `residuals()` order, index for index. Covered: pointwise
  and integral-mean targets, over intensity and absorption curves, for every
  constraint kind and every real transform, on aligned and interpolated target
  grids.
- **`MeritSensitivity::uncovered` / `is_complete()`** — rows whose dependence
  this pass does not carry: phase targets (`arg()` of a complex row, plus a
  differential reference that depends on the stack's total thickness directly)
  and color demands (a spectrum-wide integral through the CIE chain, whose
  derivative `build_needle_targets` already emits in a different shape). They
  occupy their rows and are *listed*, never returned as empty-and-covered — a
  caller that finds one among its active demands must fall back to a
  finite-difference Jacobian, and `is_complete()` is how it asks.
- **Ten cargo tests**, led by `sensitivity_is_a_finite_difference_of_residuals`
  — five kinds × two transforms, central-differenced through `residuals()`
  itself. That is the anti-drift guard: this pass walks the target grids a
  second time rather than threading a sink through `residuals_into`, which is
  the bit-exactness-critical path, and the cross-check is what keeps the two
  walks in step.

### Fixed

- **`n_residuals()` over-counted integral targets.** It reported one component
  per target point; `residuals()` pushes exactly **one** row per integral
  frame, because the constraint kind applies to the mean, not to each point.
  A spec with one nine-point integral frame advertised 9 components and
  produced 1. The count is now what the vector actually holds, and its
  docstring says the remaining caveat out loud: a frame whose grid misses the
  simulated grid is skipped by `residuals_into` and cannot be counted from the
  spec alone.

### Notes

- The kinks are the constraint, not an approximation: `a`/`b` rows are exactly
  flat on their satisfied side and `r` rows exactly flat inside the band, so
  an inactive constraint comes back with *no terms at all* rather than terms
  that happen to be zero. `Log` is likewise flat below its 1e-12 clamp — the
  residual genuinely stops responding there, and so does a finite difference.
- Bit-exactness fingerprint unchanged; `bench_refold` unchanged (one LM
  optimize 2.1 ms — this increment is not yet on any hot path).

## [0.6.7] — The LM stops squaring its own condition number (R4.4b)

The bounded Levenberg-Marquardt behind every thickness optimization solved the
**normal equations** — `(JᵀJ + λ·diag(JᵀJ))δ = −Jᵀr` — which squares the
condition number of `J`. Thin-film stacks with correlated layers are exactly
where `JᵀJ` goes singular, and the λ-floor was what kept bailing the solve out
(review §3.6). The damping semantics are unchanged; how the step is obtained
is not.

Work item A of R4.4. Backend selection and the optional argmin-ecosystem
solvers (R4.4c) are a separate item; the built-in stays the only backend.

### Changed

- **The damped step is solved by QR.** The step is now the least-squares
  solution of the augmented system `[J; √λ·D] δ ≈ [−r; 0]`, whose condition
  number is the square root of the normal-equation matrix's. It is the *same*
  step — `RᵀR = JᵀJ` exactly — obtained from the square roots. Following
  MINPACK's `qrsolv` structure, `J` is factored once per iteration (m·n²) and
  each λ trial then costs a 2n×n QR (n³), so the per-iteration work does not
  grow; building `JᵀJ` is gone entirely, since only its diagonal was ever
  needed. Unpivoted Householder: the damped system is full rank for every
  λ > 0 with a floored `D`, so pivoting would add rank diagnostics, not
  solvability. The normal-equation solve survives as the fallback when the
  factorization degenerates, and as the oracle the QR is cross-checked
  against.
- **Damping is Nielsen's gain ratio**, not a fixed ×5 / ÷3 ladder: ρ = actual
  over predicted reduction; on acceptance λ ← λ·max(⅓, 1−(2ρ−1)³) and ν ← 2,
  on rejection λ ← λ·ν and ν ← 2ν. The old ladder is kept as
  `LmDamping::Fixed` (`damping="fixed"` from Python) so the two can be
  compared on the same problem.
- **The predicted reduction is computed for the step actually taken.** The
  bound veto and clamp rewrite δ *after* it is solved. Scoring it with the
  full LM step's prediction overstates what the model promised whenever a
  bound is active, so ρ comes out too small — over-damping, and premature
  ftol exits right at the boundary. This is the one place a naive MINPACK port
  goes wrong, and it is pinned by a test that fails with ρ = 0.049 instead of
  1.0 on an exactly-linear clamped problem.
- **ftol needs both reductions** (MINPACK semantics): actual *and* predicted
  relative reduction below `ftol`. A step clipped by a bound can deliver very
  little while the model still sees plenty of room — small progress is not the
  same fact as no progress left.
- **gtol gained the scale-invariant form**, `cos∠(J·e_j, r) ≤ gtol`, alongside
  the existing ‖Jᵀr‖∞ test (`gtol_scale_invariant`, default on;
  `gtol_scale_invariant=False` from Python restores the old behaviour alone).
  ‖Jᵀr‖∞ answers a different question after a parameter is rescaled — nm
  versus µm — and the cosine answers the same one. They are different notions
  of stationarity and can exit at different points on flat valleys.

### Added

- **`LmResult.gain_ratio`** — ρ of the last accepted step, against the
  prediction for the step actually taken. NaN when nothing was accepted. A
  diagnostic, not a control: it tells a synthesis stalling at its bounds apart
  from one that is finished.
- **`validation/review/lm_check.py`** — the scipy parity §18.2 says the plan
  docs have been claiming and nobody had reproduced. Runs
  `SmatrixContext.optimize_thicknesses` and
  `scipy.optimize.least_squares(method="trf")` over the *same* residual system
  and bounds, plus the cargo tests' pinned optima recomputed by scipy.
  `validation/smoke/test_lm_parity.py` is the subset CI runs.
- **Fourteen cargo tests** in `thick_opt.rs`, including the two the plan names
  as the guards that make the port trustworthy: the QR-vs-normal-equations
  step cross-check (rel ≤ 1e-8 where both are valid, and a Vandermonde case
  where the QR is measurably better), and the clipped-step prediction.

### Notes

- Two corrections to R4.4d, both found by running it:
  - **Cost parity from a far start is not a well-posed criterion.**
    Reflectance against thickness is oscillatory, so two local solvers started
    far from an optimum legitimately land in different basins. The harness
    compares them where the basin is unambiguous, and asserts the basin-free
    properties (improvement, stationarity) on the far starts.
  - **The QR's advantage is conditional, and the harness says so.** At the λ
    the solver starts from, Marquardt damping regularizes `JᵀJ` enough that
    the two formulations are indistinguishable. The QR matters where λ has
    decayed towards nothing near a good optimum — precisely where the last
    digits are decided.
- No behaviour changed in the S-matrix engine: the bit-exactness fingerprint
  is unchanged. `bench_refold` reports one LM optimize at 2.0 ms against the
  2.1 ms baseline.

## [0.6.6] — The clippy gate is blocking (R2.3a)

CI ran `cargo clippy` from 0.6.0 onward, but with `continue-on-error: true`
against 267 findings. An advisory gate over a tree that cannot go green is a
gate nobody reads. The tree is now clean and the step fails the build.

### Changed

- **`cargo clippy --workspace --all-targets -- -D warnings` is a blocking CI
  step.** 267 findings → 0. Nothing in the diff changes behaviour: the
  bit-exactness fingerprint (12 seeded random stacks × every `Request` bit,
  needle gradients, an eigenmode landscape, hashed over the raw f64 bytes)
  is byte-identical before and after, and was re-checked after each pass.
- **108 `excessive_precision` findings came from one generated file.**
  `rust/navette/src/color/matrices.rs` is emitted by
  `validation/benches/color/gen/gen_matrices.py`, which formatted with
  `f"{x:.17e}"` — enough digits to name every double, and ~5 more than needed
  for most of them. The generator now uses `repr(float(x))`, which is the
  shortest string that round-trips. All 108 literals were verified bit-identical
  through `struct.pack("<d", ...)` against the previous tree. Fixing only the
  output would have reintroduced the warnings on the next regeneration.
- **Remaining lint decisions are crate-level `allow`s with written rationale**,
  in both `lib.rs` files, rather than scattered per-site suppressions:
  - `neg_cmp_op_on_partial_ord` — all 15 sites are `!(x > 0.0)`-shaped
    validators. The negation is what rejects NaN; clippy's suggested
    `x <= 0.0` would wave NaN through every one of them. Load-bearing.
  - `too_many_arguments` — the wide kernel signatures are what R6.1 / R6.3
    exist to restructure. Silenced so the gate can go blocking now.
  - `needless_range_loop` — flat angle-major index arithmetic
    (`k = a * num_wavs + w`), where the variable indexes several arrays at
    different strides.
  - `type_complexity`, `should_implement_trait` — plan/cache tuples, and
    inherent `from_str` parsers returning the crate's own `String` error.

### Notes

- `cargo fmt --check` stays advisory (R2.3b): making it blocking requires a
  tree-wide reformat that rewrites `git blame` for the whole crate, and the
  plan asks for a style decision before that lands.

## [0.6.5] — The one shipped example did not run (R4.3)

`examples/spectralweave_example.py` failed on a fresh clone, at line 22 of 24.
Running it was how three separate defects came to light.

### Fixed

- **The example is rewritten against the documented API** and now runs end to
  end. It used the raw native `OpticalWeaver` with a key tuple; the documented
  path is `SimulationWeaver` + `OpticalFragment`. It also handed `unweave` a
  `linspace` over the same *range* as the stored frame, which is not the same
  thing as the same *grid* — it shared exactly one wavelength out of a hundred.
  The example now weaves three fragments on three grids, takes the continuous
  curve, unweaves an edited target back onto it, asserts the round trip, and
  shows the rejection path deliberately.
- **`unweave`'s error explains the contract instead of counting.** The failure
  surfaced as `Length mismatch` raised four levels down in
  `SpectralDataFrame::set_data`, which knows nothing about grids or targets.
  `unweave` and `unweave_collection` now check coverage where the context
  exists and say how many of the frame's wavelengths the supplied grid carried,
  which frame (by span), that matching is exact to 1e-12 rather than
  interpolated, and where to get a grid that works. `set_data`'s own message
  now names the frame and both counts.
- **`OpticalFragment` is hashable, so `unweave_batch` can be called at all.**
  Its signature is `dict[OpticalFragment, np.ndarray]`, but `frozen=True` on a
  dataclass holding two numpy arrays synthesises a `__hash__` over those
  arrays, which raises `unhashable type: numpy.ndarray` — and a `__eq__` that
  raises `truth value of an array is ambiguous`. `eq=False` gives identity
  semantics for both, which is the honest answer for a container of mutable
  buffers. Found by writing the example's last line.

### Added

- **`SimulationWeaver.unweave` rejects a non-`OpticalFragment` template** with
  a `TypeError` naming what it got and how to build one. The backend's key
  tuple is an implementation detail; passing one used to surface as
  `AttributeError: 'tuple' object has no attribute '_rust_key'`. A curve whose
  value count does not match its wavelength count is now caught here too,
  rather than becoming a distribution-plan error about something else.
- **`validation/smoke/test_examples.py`** — every `examples/*.py` runs in a
  subprocess (repo root as cwd, `PYTHONIOENCODING=utf-8`, `MPLBACKEND=Agg`,
  120 s timeout) and must exit 0 *and* print something. A subprocess keeps
  pytest's process state and any module-level side effects out of the suite.
  An empty `examples/` fails the guard rather than passing vacuously.
- **`validation/smoke/test_spectralweave_wrapper.py`** — 10 tests over the
  round trip, the coverage message's contents, the guards, and
  `OpticalFragment` as a dict key. Includes the positive control that matters:
  a grid covering *one* of two frames exactly is still accepted, so the new
  check cannot start over-rejecting partial updates.

### Evidence

Teeth proven by dropping two probe files into `examples/`: one exiting 3
(caught by the returncode assert) and one that runs silently (caught by the
output assert). 656 passed + 2 skipped; cargo 376; zero `cargo check`
warnings; all nine review harnesses exit 0.

### Known gaps

- `examples/` holds exactly one file. The runner is written to cover whatever
  lands there, but there is no second example to prove the parametrization
  against — the probe files above stood in for one.

## [0.6.4] — Color demands were optimized against zero in Python (R4.2)

The documented Python needle flow is `build_needle_targets` →
`needle_gradient`. Every pointwise demand came through it correctly. A
**color** demand — Lab/DE2000, White, Yellow, dominant wavelength — came
through as a gradient of exactly zero, with no warning and no error.

A color demand integrates the whole spectrum into one residual, so it has no
per-point target to fold into the `r`/`t` (target, weight) pairs. The native
fold emits its analytic `dF/dcurve` per solver point instead
(`NeedleTargets.grad_r`/`grad_t`), and `needle_pass.rs` adds those into the
same accumulator as the pointwise terms. But the **Python** dict from
`build_needle_targets` dropped both arrays, and `needle_gradient` had no
argument that could have accepted them — so the folded `r` pair was all
zeros, the engine faithfully returned zero, and a caller assembling their own
needle cycle in Python optimized a color demand against nothing.

Not a crash and not a wrong shape: a silent wrong-optimizer. `run_design` was
never affected — it folds and deposits entirely inside Rust.

### Added

- **`build_needle_targets` returns `"grads_r"` / `"grads_t"`** — flat `na*nw`
  f64 arrays, the same angle-major layout as the target/weight pairs. They are
  *not* a (target, weight) pair: each entry is the chain-rule factor
  `g = dF/dcurve` for one solver point, with the demand's weight, its current
  residual and the U-curve half already folded in. All-zero for a spec with no
  color demands, so passing them unconditionally is safe.
- **`needle_gradient` accepts `grads_r` / `grads_t`** (Python keyword, PyO3
  method, and the free `needle_engine` entry), deposited into `P` and `P_T`
  respectively via the existing `p_coherent_grad_r_from_fields` /
  `p_coherent_grad_t_from_fields` kernels — the same kernels, with the same
  zero-skip, that `needle_pass.rs:889-896` calls. The native internal path is
  untouched; its cargo tests still pin the deposit semantics.
- **A non-zero bucket handed to a channel that is not being computed is an
  error**, naming the missing bit. Dropping it silently would be this same bug
  wearing a new hat. An all-zero array is not a demand and passes, because the
  fold hands both arrays through unfiltered. Non-finite entries are rejected by
  index, matching the R3.1/R3.2 convention.
- `validation/review/color_grad_python.py` (new, 24 checks, exits 1 on drift)
  and `validation/smoke/test_color_needle_python.py` (new, 13 tests — the part
  CI runs).

### Changed

- `navette.synthesis`'s module docstring — the documented flow — now shows the
  color branch with a worked two-channel call, and says outright that omitting
  the two kwargs loses the color contribution.
- The Part C comment in `validation/review/color_merit_check.py` no longer
  describes the Python dict as omitting the buckets; it now says why that
  harness deliberately keeps hand-assembling the chain rule (it is the oracle
  the new harness is checked against, so it must not consume the engine's
  own answer).

### Evidence

`grads_r` vs a central-difference `dF/dR` on the sim row: max relative
deviation **7.1e-10** across all 31 points. The deposit against the pointwise
kernel driven to the same scalar (`g·P_ref/(2R)`, the two differ only in one
factor): **2.3e-15 … 3.1e-15** — 1e-12 is a real bound here, not a rounded
1e-6. End to end, `dF/d(thickness)` through the Python path against a
thickness FD of the merit, for all three layers: **1.6e-8 … 1.7e-7**; against
`color_merit_check.py`'s hand-assembled chain rule: **8.3e-16 … 9.4e-16**.
Superposition `P(pointwise + color) = P(pointwise) + P(color)` is bit-exact.

Teeth proven by disabling the R deposit and rebuilding: 2 of the 13 smoke
tests and 10 of the 24 harness checks fail; restored and re-verified.

### Known gaps

- Only the front R and T channels carry color buckets — that is the native
  fold's own v1 scope (`CurveId::Ru`/`Tu` take the ÷2 U-curve half; the back
  siblings and absorption fold nothing). The Python path now mirrors exactly
  what the native path computes, no more.
- The deposit rides `P`/`P_T`, so a color demand requires those bits. That is
  a real constraint, now stated in an error rather than discovered by a
  gradient that quietly reads zero.

## [0.6.3] — Every dependency floor was fiction (R4.1)

`requires-python = ">=3.12"`, and **not one** of the declared floors has a
cp312 wheel:

| declared | first cp312 wheel | now |
|---|---|---|
| `numpy>=1.22.0` | 1.26.0 | `numpy>=2.0` |
| `scipy>=1.8.0` | 1.11.2 | `scipy>=1.13.0` |
| `pyyaml>=6.0` | 6.0.1 | `pyyaml>=6.0.1` |
| `numba>=0.56.0` (extra) | 0.60.0 | `numba>=0.61.0` |

A floor is a promise that the package works with at least that version. These
could not be installed at all on the Python the project declares.

### Fixed

- **`numpy>=2.0`** — the floor R4.1 named. Two independent floors meet in this
  package and they are not the same thing: the extension is `abi3-py312` (one
  wheel for 3.12 and every later 3.x), while the `numpy` Rust crate (0.28)
  compiles against the numpy **2** C-API and makes no abi3-style promise of
  its own. A wheel built that way cannot load against numpy 1.x whatever the
  Python ABI says. Both floors are now documented next to each other in
  `pyproject.toml`.
- **`scipy>=1.13.0`** — the first scipy that both supports numpy 2 and ships
  cp312 wheels. Beyond R4.1's stated scope; found while making the gate
  blocking, and the same one-line defect.
- **`pyyaml>=6.0.1`**, **`numba>=0.61.0`** (extra) — same.

### Changed

- **`ci.yml`'s floor job is blocking, and no longer checks numpy by name.**
  It reads *every* `>=` floor out of `pyproject.toml`, pins them all at once
  with `--only-binary=:all:`, runs `pip check`, and then runs the full suite
  against that combination on the **lowest** supported Python. Checking numpy
  alone would have left three broken floors behind a green check — which is
  exactly how they survived. Renamed `numpy-floor` -> `dependency-floors`.
- The job now runs the whole suite rather than importing. A floor that imports
  but breaks a kernel is still a lie, and pinning the numba extra means the
  parity layer compares against its reference there instead of skipping.

### Verified

Locally, before pushing, on a clean Python 3.12.13 venv: the four floors
install together as wheels, `pip check` is clean, the PEP 517 build produces a
**release** extension, and `pytest validation` gives **631 passed, 2 skipped**
at `numpy 2.0.0 / scipy 1.13.0 / pyyaml 6.0.1 / numba 0.61.0` — identical to
the development environment.

## [0.6.2] — Eigenmode search is boxed, and says when it found nothing (R3.3)

`char_func` is `|1/r(n_eff)|^2`, so a pole drives it to zero — but so does
letting `|n_eff|` run away, and the search was unbounded. Seeded on a stack
with no s-polarized mode, it reached `n_eff = -2.0e8 + 8.4e6j` and reported a
characteristic value of `1.3e-217`: twelve orders *better* than the real
surface-plasmon pole it never found, for a trial index two hundred million
times the largest index in the stack. Nothing said anything.

### Fixed

- **The minimizer is confined to a physical box**, `3 x max |n|` over the
  stack. Every guided mode satisfies `min|n| <= |n_eff| <= max|n|`, so 3x is
  generous even for leaky and substrate-side modes. Candidates are *projected*
  into the box rather than rejected, so Nelder-Mead stays well-defined and
  slides along the boundary instead of reflecting off an invisible wall.
- **`refine_mode` refuses to call something a mode when the characteristic
  value says otherwise**: new `max_residual` kwarg, default `1e-6`, raising
  `ValueError` with the settled `n_eff`, the value, and what to do instead.
  Pass `max_residual=None` for the old unconditional behaviour. The margin is
  not tight — on the pinned SPP a real pole reaches `1e-17`, eleven orders
  below the threshold, while the boxed non-mode settles around `1e-3`.
- `char_func` now reports a non-finite `r` as `1e30` rather than propagating
  NaN into `find_minima`, whose comparisons are undefined on NaN.

### Changed

- **The box is in `char_func_xy`, not `char_func`.** The plan put it in the
  shared `char_func`, which the landscape scanner also calls — that would
  paint `1e30` across any part of a user-requested scan range lying outside
  the box, corrupting a diagnostic the caller explicitly asked for at the
  range they asked for it. The minimizer picks its own points and needs walls;
  the scanner is already bounded by its caller. A test pins the distinction.
- `find_eigenmodes` is untouched and does **not** go through the residual
  contract: its seeds come from a bounded landscape scan, which is what makes
  them seeds. On this stack it correctly returns `[]` for s-polarization — the
  documented workflow was never the broken one.

### Added

- `validation/smoke/test_eigenmodes.py` — 19 tests. Nothing pinned the
  eigenmode path before, which is why the runaway survived a full review
  cycle. Covers the SPP from three seeds, the field profile peaking at the
  metal/glass interface (the physics check a characteristic value cannot
  give you), the runaway from five seeds including two already outside the
  box, the flat s-polarized landscape, the threshold behaviour, and the
  landscape *not* being walled off. Removing the box fails exactly seven of
  them, including `DID NOT RAISE` on the original runaway seed.

### Known gaps

- The plan pinned the SPP at `n_eff = 1.0459458 + 0.0015949j, val < 1e-10`.
  That is not reproducible from the stack as described (the wavelength is not
  stated), and the nearby feature this code converges to — `1.0471183 + 0j`,
  val `1.8e-6` — is a shallow resonance, not a pole: `Im(n_eff) -> 0` on a
  lossy metal, and it sits below the substrate index. The genuine pole on this
  geometry is the glass-side plasmon at `1.7139417 + 0.0226069j`, val `1e-17`,
  and that is what is pinned. See the plan's R3.3 notes.

## [0.6.1] — Needle depth is range-checked in release builds (R3.2)

`needle_slopes4_ddz` guarded its host-layer invariant with `debug_assert!`,
which is compiled out of every release build — the one everybody runs. And
`locate_depth_in` does not reject an out-of-range depth either: it falls
through to the last layer of the range and returns a depth past that layer's
thickness. So `z = 1000` on a 400 nm stack, or `z = -10`, returned a gradient
for a needle that is not where the caller put it, with no error.

### Fixed

- **`z` outside the eligible span now raises `ValueError`**, naming the index,
  the value and the span. The two paths have different spans *within the same
  call* — the coherent kernels are confined to `[start_idx, end_idx]`, the
  multiblock cascade walks every non-ambient layer — so each is checked only
  when the request actually reaches it. A multiblock request is not held to
  the block's narrower bound; a call using both must satisfy both.
- **A non-finite `needle_n_per_wav` entry now raises**, naming the index. It
  previously made the whole gradient NaN with nothing to point at. (Beyond
  R3.2's stated scope, but the same check in the same place.)

Both endpoints stay legal and a `1e-9` tolerance absorbs round-off: a caller
writing `np.linspace(0, sum(d), k)` lands on the bottom endpoint with whatever
error the sum accumulated, and rejecting that would make the obvious way to
write the call fail intermittently.

### Changed

- The check lives in `solver::needle_gradient`, not in the two PyO3 entry
  points. The plan put it binding-side; one implementation in the core covers
  both bindings *and* Rust callers, and it returns `Result<_, String>` which
  the bindings already map to `PyValueError` — so the user-visible error is
  identical either way.
- The kernel keeps its `debug_assert!`s as Rust-caller invariants, with
  messages that now say which contract was skipped. A release `assert!` in an
  O(1) hot kernel would trade a silent wrong answer for a panic, which is
  worse on a library path.
- `validation/review/garbage_in.py`: the three needle rows flip from
  `silent-clean` / `NaN-in-output` to `raises`. The fourth (`n = 0.3`) stays
  permissive — an index below 1 is a real metallic index, not garbage.

### Added

- 13 more rows in `validation/smoke/test_input_validation.py`, including that
  `end_idx` narrows the span, that a multiblock request keeps the wide one,
  and that both endpoints plus the tolerance band are accepted.

## [0.6.0] — `ScatterMatrix` rejects malformed input (R3.1)

**Behaviour change.** The constructor now validates its inputs and raises
`ValueError`. Inputs that used to return plausible numbers for a stack nobody
asked for are errors; the message names the offending index and value.

### Fixed — these used to be silent

| input | old behaviour |
|---|---|
| negative thickness | same numbers as **deleting the layer** |
| NaN / inf thickness | same numbers as **deleting the layer** |
| NaN or inf refractive index | every output NaN, nothing naming the layer |
| refractive index with `\|n\|` past `sqrt(DBL_MAX)` | `n**2` overflowed inside the solve; NaN out |
| NaN wavelength | NaN output |
| wavelength <= 0 | `k = 2*pi/lambda` divides by zero |
| duplicated wavelength | NaN `GD`/`GDD`/`TOD`/`FOD` (the kernels divide by the grid spacing) |
| descending wavelength grid | worked, but silently a second convention |
| angle > 90 deg (e.g. 120 deg) | **aliased onto its mirror** (60 deg) — only `sin(theta)` reaches the engine |
| negative or NaN angle | aliased / NaN |

Descending grids are **rejected, not sorted**: silently reordering would
desynchronize the grid from the caller's own wavelength-indexed arrays, and
nothing would say so. Sorting is the caller's call.

### Unchanged — deliberately still accepted

Grazing incidence (90 deg, `R = 1`), a single-point wavelength grid, zero
thicknesses, a 1e9 nm layer (no upper cap — verified safe), and metallic
indices with `n < 1`. A validation layer that rejected any of these would
have broken the library to fix a bug; `test_input_validation.py` holds that
line with 9 positive controls.

The native `Solver` underneath stays permissive. It is also the optimizer's
inner loop and the Rust test surface, where a re-check per call is pure
overhead — so the checks live at the Python user surface, run once per
construction, and have no opt-out flag.

### Added

- `validation/smoke/test_input_validation.py` — 30 tests: 17 rejected inputs
  (each asserting the index and value appear in the message), 9 accepted
  ones, and four positive controls, including that a rejected construction
  does not mutate the caller's arrays and that the native `Solver` is never
  constructed for a stack that is about to be rejected.

### Changed

- `validation/review/garbage_in.py` recorded the old silent behaviour as
  prose. It now carries an expected verdict per case and **exits 1** when one
  drifts (it previously computed an `OK` flag it never used and always
  exited 0). The four `needle_gradient` rows are still permissive and are
  marked `[R3.2]` rather than quietly listed.
- `validation/parity/smatrix/test_physics_mirror.py` — the Nelder-Mead
  thickness search is unbounded and its minimum sits on the `d = 0` boundary,
  so the simplex reached for negative thickness. That was silently accepted
  before, meaning the optimizer was steered by the merit of a *different*
  stack; it now gets a sloped barrier at zero.

### Known gaps

- `needle_gradient` is unguarded: a NaN needle index still gives NaN, and
  `z` outside the stack or negative is accepted silently. That is R3.2.
- `roughness_values` and `incoherent_flags` are not validated (negative sigma
  is still accepted). Not in R3.1's scope; no silent-wrong-answer path is
  known for them.

## [0.5.9] — Request-bit and schema sync guards (R2.5)

Three constant tables are written once per language with nothing tying them
together: the 49 `REQ_*` solver bits, the 18 `NREQ_*` needle bits, and two
schema versions. Writing the guard found one of them already broken.

### Fixed

- **`PROGRAM_SCHEMA_VERSION` was decorative.** `config.rs`'s `gate_document`
  matched the envelope version against the literal `Some(1)` while quoting the
  constant in its rejection message, so bumping the constant would have made
  the gate reject exactly the version it claimed to read, and accept the one
  it claimed was stale. It compares against the constant now. Demonstrated:
  with the literal in place, setting the constant to 2 left the gate accepting
  1 and the new sync test green; with the fix, the same edit fails the test
  (`python PROGRAM_SCHEMA_VERSION=1, rust accepts 2`).

### Added

- `validation/smoke/test_request_bits.py` — 65 tests, guarding each table from
  both ends:
  - **Source-level**, parsing `pub const REQ_*` / `NREQ_*` out of
    `core_engine.rs` and `needle_engine.rs` and comparing name-for-name and
    bit-for-bit against `Request` / `NeedleRequest`. Catches a constant added
    on one side only — which the plan's semantic probe alone cannot see, since
    a bit Python never mirrored is a bit Python can never request. Skips when
    the Rust tree is absent (installed wheel).
  - **Behavioural**, driving all 49 single bits through the engine and
    asserting the returned channel keys exactly, plus the six convenience
    bundles and an all-bits-at-once call. Catches a renumbering both source
    tables agree on. Proven by swapping `RS`/`RP`: 3 failures.
  - **Density/uniqueness** on both flags, and completeness of the bit ->
    output-key map (two bits claiming one key is otherwise undetectable).
  - **Schema**, probing the two native gates at runtime rather than parsing
    them — `check_schema_version` and `load_document` are thin shims over
    Rust, so the version Rust accepts is observable from Python, wheel or
    checkout. Untagged states and untagged documents must be refused.

### Known gaps

- The `REQ_*` constants are still not exported by the extension, so the
  source-level half is the only exact check and it needs the Rust tree. The
  behavioural half covers the wheel case.

## [0.5.8] — Parity tests were passing without comparing anything (R2.4)

`validation/conftest.py` had ignored the whole `parity/` tree, so nobody had
run those files under pytest. Collecting them was supposed to be a one-line
change. It uncovered two defects that made five of the seven smatrix parity
scripts worthless as tests.

### Fixed

- **The numba reference never loaded, and the scripts called that a pass.**
  Five parity scripts imported `loom_matrix` from `validation/parity/loom/`,
  a directory that does not exist; the reference actually lives in
  `validation/parity/smatrix/refs/`. The `except ImportError` branch then
  scored every comparison `"PASS (rust-only)"` — a green verdict for a run
  that compared the Rust result against nothing. The import now points at
  `refs/`, and a genuinely missing reference **skips** instead of passing.
- **Failure exited 0.** The same scripts printed `OUTPUT_STATUS FAIL` and
  returned success, so no caller — CI, a shell loop, a human — could gate on
  them. Each now ends in a shared `report()` verdict, a `test_parity()`
  assertion for pytest, and `sys.exit(1)` for direct invocation. Verified by
  injecting a 1e-3 error into `test_w_function.py`: the script exits 1 and
  pytest reports FAILED.
- **`test_needle_t_a_phi.py` collected zero tests** — it had a `main()` and
  no test function. Given a pytest wrapper that asserts `main() == 0`.
- **Two `core_engine_*` scripts crashed collection** with a module-level
  `sys.exit(1)`, which pytest reports as INTERNALERROR. They now skip with
  the real reason (see Known gaps).

### Added

- `validation/parity/_parity.py` — shared preamble for the parity layer:
  `console_utf8()` (cp1252 consoles choked on `μ`), `require_reference()`
  (skip under pytest, exit 0 when run directly), and `report()`.

### Changed

- `validation/conftest.py` collects `parity/` and ignores only what is not a
  test: `benches`, `refs/*` (the reference implementations themselves),
  `gen_*.py` fixture generators, and `golden_mirror.py`.
- Suite is now **504 passed, 2 skipped** (was 450). The plan estimated the
  directory held 11 pytest-style tests; it holds 17 files and 54 tests.

### Known gaps

- **R2.4a**: `test_core_engine_photometry_only.py` and
  `test_core_engine_rigorous_ellipsometry.py` load a standalone
  `navette_matrix` extension from `validation/parity/target/release/`. That
  crate does not exist anywhere in this repository, so the plan's remedy
  ("document its build step") is not possible — the two files need porting
  onto `navette._smatrix.core_engine` and its request mask. They are the only
  whole-engine parity oracles; every other parity file covers one kernel.
- The `edge_nc` case in `test_solve_coherent_block_fields.py` is
  self-referential: `refs/loom_matrix.py` was edited in lockstep with the
  Rust change it is supposed to check independently.

## [0.5.7] — CI fixes found by running it (R2.3)

The first live run of `ci.yml` was green on all three blocking jobs (rust,
python/windows, python/ubuntu) and surfaced two defects in the workflow
itself. Inspection would not have caught either.

### Fixed

- **The numpy-floor job failed for the wrong reason.** With no cp312 wheel,
  pip fell back to building numpy 1.22.0 from sdist and died with
  `BackendUnavailable: Cannot import 'setuptools.build_meta'` — a build-env
  error that buries the actual finding. It now installs with
  `--only-binary=:all:` and, on failure, prints the finding as a
  `::error::`: the declared floor has no wheel for the declared
  `requires-python`. (Still advisory; R4.1 fixes the floor itself.)
- **`actions/checkout@v4` and `actions/setup-python@v5` are Node 20**, which
  the runners now force onto Node 24 with a deprecation warning. Bumped to
  `@v5` / `@v6` (both `node24`) in `ci.yml` **and** `release.yml`.

### Known gaps

- `release.yml` still pins `actions/upload-artifact@v4` and
  `download-artifact@v4`. Left alone deliberately: that workflow only runs on
  a tag, so a major bump there cannot be verified before it matters, and it
  is the workflow that publishes to PyPI.

## [0.5.6] — push/PR CI gate (R2.3, §7)

### Added

- **`.github/workflows/ci.yml`** — runs on every push and pull request.
  Until now the only workflow was `release.yml`, and it only ran on a tag:
  a broken example, unsynced request bits and a tree full of live compiler
  warnings all survived the 0.5.0 release with nothing to catch them.

  Blocking: `cargo test --workspace`; `cargo check --workspace --all-targets`
  under `RUSTFLAGS=-D warnings` (possible as of 0.5.5); `pytest validation` on
  Windows **and** Linux; `tools/check_exposure.py` — the README has claimed
  this is "enforced in CI" for some time and it now is; `tools/check_cie_sync.py`;
  and an assertion that the PEP 517 install really is a release build, which
  exercises R2.1's `build_profile()` probe end to end.

  Advisory (`continue-on-error`), each with the reason recorded in the file:
  `cargo clippy` (248 findings), `cargo fmt --check` (981 files differ from
  rustfmt defaults), and the numpy-floor job (see below). Making any of them
  blocking means landing a large mechanical diff first; doing that inside the
  commit that introduces the gate would bury the gate. A CI that cannot go
  green teaches people to ignore CI.

### Fixed

- **`release.yml`'s `check-versions` could not catch the version skew that
  actually happened.** It compared the tag against `pyproject.toml`,
  `Cargo.toml` and `__about__.py` — but not `Cargo.lock` (which lagged at
  0.5.1 while `Cargo.toml` said 0.5.2, see 0.5.3) nor the internal `navette`
  path-dependency requirement in `Cargo.toml`. All six are now checked, each
  mismatch reported individually as a `::error::` rather than one opaque
  "version mismatch". Verified both ways against the working tree.

- `tools/check_exposure.py`: allowlisted `max_disp_order`. It had counted as
  exposed only because `navette-py/src/smatrix.rs` carried an **unused `use`**
  of it; deleting that dead import in 0.5.5 surfaced it. A bare `use`
  satisfies this lint, so an unused-import cleanup can legitimately turn it
  red — worth knowing before assuming a red run means a lost binding.

### Known gaps

- The **declared numpy floor is fiction**: `pyproject.toml` says
  `numpy>=1.22.0` while `requires-python` is `>=3.12`, and numpy only gained
  cp312 wheels at 1.26. The CI job that would catch this is present but
  advisory; **R4.1 must raise the floor and delete its `continue-on-error`**.

## [0.5.5] — Zero compiler warnings (R2.3 prerequisite, §7)

Prerequisite for the push/PR CI gate: `-D warnings` cannot be enabled while
the tree emits any. All 25 are gone; the build is warning-free.

The §7 inventory listed 13. `cargo check --workspace --all-targets` finds 25 —
the extra 12 are four `unused_mut`, one never-read enum field, one
non-snake-case binding, two private-interface warnings, and four pyo3
deprecations. As with R1.1's site list, the inventory was a lead, not a census.

### Fixed

- **Two dead computations, not just dead names.** `solver::find_minima` copied
  the whole landscape into `land_vec` and bound `land` to it; neither was ever
  read (the median is taken from a second copy). The copy is gone — one fewer
  full-grid allocation per eigenmode scan. `IntGap::EdgeMean`'s first field
  (the op mean) was never read either: the edge form uses `d_bar`, re-derived
  from `opm` at the consumption site. The variant now carries only the
  bandwidth.
- Unused imports (`NeedleTargets`, `cplx`, `PI`, `max_disp_order`,
  `rayon::prelude`, `ArraySeed`, two glob imports), unused bindings (`land`,
  `m` ×3, `wavelengths`, `num_angles`, `ok`), and four redundant `mut`.
  `NeedleTargets` is used only by `cycle.rs`'s tests, so it moved into the
  test module rather than being deleted.
- `PyFilmInput` was private while appearing in the `pub(crate)` signatures of
  `assemble_design` and `run_design`; it is now `pub(crate)` to match.
- pyo3 0.28 deprecation: `#[pyclass]` types deriving `Clone` get an automatic
  `FromPyObject` that is becoming opt-in. The five `config_type!` classes and
  `PyLayerSpec` now declare `from_py_object` explicitly, which **preserves
  today's behavior** rather than silently losing the conversion at the next
  pyo3 bump.

Behavior-neutral by construction, and verified: eigenmode scan / coarse-minima
/ landscape checksums are **bit-identical** across the change (three
resolutions × three median factors, built before and after and diffed), and
450 pytest + 376 cargo tests pass.

## [0.5.4] — Build-profile probe + bench guard, UTF-8 benches (R2.1, R2.2, §5.1.1, §17)

### Added

- **`navette.build_profile()`** — reports the Cargo profile the installed
  extension was compiled with, `"release"` or `"debug"`. A debug build imports
  and computes correctly but runs several times slower, so nothing detects it
  by eye; this makes it a one-line check.
- **`validation/benches/_bench_common.py`** — shared bench preamble.
  `require_release()` exits rather than time a debug build; `setup_bench()`
  reconfigures stdout/stderr to UTF-8 and puts the repo's `src/` on
  `sys.path`; `bench_provenance()` returns the profile/version/platform block
  now stamped into bench results JSON.
- All **seven** benches (`smatrix/bench_backside_speed.py`,
  `structure/bench_grid_assert.py`, `synthesis/bench_refold.py`,
  `spectralweave/navette_{spectral,target}_bench.py`,
  `interpolate/1dinterpol_test_bench.py`, `color/bench_validate_color.py`)
  call the preamble before importing `navette`.
- `validation/smoke/test_build_profile.py` (12 tests) — asserts the gate
  refuses `"debug"` and an unbuilt extension, that the package-level accessor
  matches the native symbol, and that **no bench can be added without
  `require_release()`**. It deliberately does *not* assert the installed
  profile: testing against a debug build is fine, *timing* one is not.

### Fixed

- **Benches died with `UnicodeEncodeError` on a cp1252 console** (§17) — the
  interpolate, target and color benches print `->` arrows and emoji and needed
  `PYTHONIOENCODING=utf-8` to run at all. `setup_bench()` fixes this at the
  source; verified by running every bench with the variable unset.
- Benches that did `sys.path.insert(0, "src")` only worked when launched from
  the repo root; they now resolve `src/` from `__file__`.

### Changed

- **Every `maturin develop` build instruction now says `--release`** — the
  README, `navette/__init__.py`, and the sixteen `ImportError` build hints in
  the wrapper modules. The README carries an explicit warning; following any
  of the old hints produced exactly the debug build of §5.1.1.
- Committed bench results (`validation/benches/smatrix/results/*.json`)
  predate the gate and are of **unknown profile** — `bench_backside_speed.py`'s
  `A_legacy` case measures ~0.06 ms on a release build against the 0.50 ms
  recorded there, so those files are almost certainly debug-build numbers and
  must not be compared against new runs. Noted in `validation/README.md`.
- `tools/check_exposure.py`: allowlisted `nevot_croce_factors` (internal, via
  every solve/field/needle path). The lint had been **failing since 0.5.3** —
  it is not wired into any workflow yet, which is R2.3.

## [0.5.3] — Névot-Croce: remaining two sites + single source of truth (R1.1, §3.2)

### Fixed

- **Two further Névot-Croce sites still applied the reflection factor to
  transmission.** R1.1 (0.5.1) fixed the two sites the remediation plan listed;
  the type-5 branch actually exists in **four** places. Still broken were:

  | Site | Reached from |
  |------|--------------|
  | `solver::field_prof` | `ScatterMatrix.field_profile()` |
  | `needle_operator::interface_matrix` | `needle_gradient()` / needle synthesis |

  The needle case mattered most: synthesis merit is evaluated through
  `coherent_block` (fixed in 0.5.1) while its needle sensitivity ran through
  `needle_operator` (not fixed), so the P-function stopped being a derivative
  of the merit it optimizes. Measured on a 21-layer σ = 4 nm stack, the needle
  gradient moved by **20.9 %**; the field profile by 0.08 %.

### Changed

- **`optics_core::nevot_croce_factors()` is now the single source of truth**
  for the type-5 factors. All four interface builders call it instead of
  open-coding the exponentials, which is what allowed a partial fix to land.
  Bit-identical output on every path (verified by checksum: roughness type 0
  is unchanged everywhere, and `coherent_block` type 5 is unchanged).
- `RoughnessType` and `ScatterMatrix.energy_conservation()` docstrings now
  state the observable consequence of Névot-Croce's perturbative unitarity:
  **`T` can exceed 1 and the residual `A = 1 - R - T` can go negative** outside
  `|kz·σ| ≪ 1`; neither is clamped, and `energy_conservation()` reports a
  magnitude, so an energy *gain* is indistinguishable from absorption.
  Measured: n 1 → 4.28, σ = 20 nm, 550 nm gives T = 1.0768, A = −0.2347.
- `Cargo.lock` was one version behind (it recorded 0.5.1 while `Cargo.toml`
  said 0.5.2). `release.yml`'s `check-versions` greps `pyproject.toml`,
  `Cargo.toml` and `__about__.py` but never the lock, so this passed CI
  silently while dirtying the tree on any `cargo build`.

### Added

- `validation/smoke/test_rtype5_cross_path.py` (15 tests) — makes a partial
  fix impossible to repeat, three independent ways:
  - **index-matched invisibility**: a type-5 interface between two identical
    media must be bit-identical to no interface at all (r = 0 and the
    transmission factor collapses to exactly 1), checked on `compute()`,
    `field_profile()` and `needle_gradient()` at zero tolerance;
  - **needle ↔ solver agreement**: `2·P_T` must equal a central difference of
    `ΣT²` taken through the ordinary solver, to 1e-6, on a rough stack;
  - **source guard**: every `rtype == 5` branch must call the shared helper,
    and the reflection exponential may be written out in `optics_core.rs` only.

  Verified to have teeth: reintroducing the old factor fails 11 of the 15,
  while the pre-existing 423-test suite passes on that same broken build.

## [0.5.2] — energy_conservation() accepts 2-D input (R1.2, §15)

### Fixed

- **`ScatterMatrix.energy_conservation()` no longer raises `TypeError`.** The
  wrapper feeds 2-D `[n_angles, n_wavs]` arrays to the native
  `solver_energy_conservation`, which previously accepted only 1-D input, so
  the advertised API failed on every call. The binding now accepts 2-D, checks
  shape agreement, and reshapes the elementwise residual back to 2-D; the
  single-angle path still squeezes to 1-D.

### Added

- `validation/smoke/test_energy_conservation.py`: lossless (~0 residual),
  absorbing (residual == absorptance, monotone in thickness), 2-D shape, and
  1-D squeeze coverage.

### Changed

- `test_package_version` validates semver format instead of pinning a literal,
  so it survives per-item version bumps (release workflow still enforces sync).

## [0.5.1] — Névot-Croce transmission factor (R1.1, §3.2)

### Fixed

- **Energy conservation in the Névot-Croce roughness path (roughness type 5).**
  The type-5 branch previously applied the reflection Debye-Waller factor
  `f = exp(-2·kz1·kz2·σ²)` to the transmission amplitudes as well, silently
  destroying transmitted energy at rough interfaces. Transmission now takes the
  canonical wavevector-difference factor `exp(+((kz1-kz2)·σ)²/2)` (Névot & Croce
  1980; de Boer, Phys. Rev. B 44, 498 (1991); Stearns, J. Appl. Phys. 65, 491
  (1989)), so `R + T → 1` for lossless interfaces.

  Measured single interface, n = 1.0 → 1.5001, λ = 550 nm, normal incidence
  (no absorption — energy must be conserved):

  | σ (nm) | old R+T (buggy) | new R+T (fixed) |
  |--------|-----------------|-----------------|
  | 0      | 1.000000        | 1.000000        |
  | 3      | 0.992977        | 1.000001        |
  | 6      | 0.972202        | 1.000016        |
  | 10     | 0.924678        | 1.000125        |

  **This changes numerical results for every stack using roughness type 5**
  (XRR is the main consumer). Fixed at both `coherent_block.rs` call sites; the
  in-repo loom parity mirror was patched in lockstep. **Incomplete** — the
  remediation plan's site inventory said "two verbatim sites" but there are
  four; `solver::field_prof` and `needle_operator::interface_matrix` were
  missed and are fixed in 0.5.3.

### Added

- `validation/smoke/test_rtype5_energy.py`: energy-conservation sweep
  (rtype × σ × contrast × angle), a regression guard against re-damping
  transmission, and an independent closed-form Névot-Croce literature anchor
  checking engine R and T individually to machine precision.
