> **SUPERSEDED (phasing only) — see `docs/implementation_plan.md`.**
> The §7 phasing and its 0.4.32-0.4.35 version reservations are dead: those
> numbers shipped as other work. The executable plan renumbers them to 0.6.39+
> and sequences them after the gradient work, because a design film that
> expands to a multi-row span cannot be routed by row index. Everything else
> in this document — the segment model (§3.1), the five rejected alternatives
> (§3.2), the detailed design (§4), the eight v1 limitations (§6) — remains
> authoritative.

# Multi-environment optimization — proposal & plan

STATUS: proposal (requested 2026-09-07, **rewritten same day**: the first
draft scoped environments to half-space pairs only — insufficient. The
design below scopes them to full surrounding stacks, which subsumes the
halves-only case as degenerate). Not implemented. No version bump.
Reserves **0.4.32+** on acceptance; the deferred opacity spec (§D11 of
`color_targets_optionB_plan.md`) moves to next-free on recall.

## 0. TL;DR

Optimize **shared designated segment(s)** evaluated inside **K different
surrounding stacks** in one joint run. Each environment sandwiches the
same optimizable stack(s) between its own fixed front/back layer
sequences (different materials, counts, thicknesses per env — e.g. the
same coating behind different cover glasses, on different substrates,
in different housings). New machinery vs today: segmented assembly
(fixed + design segments concatenated per env), parameter identity by
**film name** across envs (positions differ — the fold/buckets are
positional today), needle insertion restricted to design segments via
the existing per-film flags, and K solves per eval with env-tagged
demands. Merit/LM/fold kernels untouched — only assembly, routing, and
the driver loop change.

## 1. Goal & use cases

Shared design variables, differing contexts, one joint optimum:

- **U1 — same coating, different covers:** one AR stack deployed bare
  *and* laminated under different cover glasses (different front
  sequences, not just different ambient n).
- **U2 — same coating, different substrates/systems:** front coating
  identical, back surroundings differ (different substrate, different
  detector stack behind it).
- **U3 — two coatable segments, different housings:** front-side AR +
  back-side mirror optimized jointly, each environment adding its own
  fixed window/absorber layers around them.
- **U4 — immersion with different windows:** dry vs oil vs water paths
  where each path includes its own window stack (again: sequences, not
  scalars).

Non-goal: per-environment design variants (that is K independent
designs). The design segments are **one object** shared by all envs —
insertions and thickness steps propagate everywhere by construction.

## 2. What already exists (no new functionality needed)

### 2.1 Per-demand angle / illuminant / band / curve — DONE
`angle=`, `illuminant=` (+ own white), `wavelength_range=[lo, hi]`,
`curve=` incl. backside (`RB*`/`TB*` ride `SimCurves::back` from the
same solve). All five "multi" axes except the surroundings are solved
patterns; the env-tagged demand model from §4 slots straight onto them.

### 2.2 Per-film optimize/needle flags — DONE
`LayerSpec/ArrayFilm` carry `optimize` + `needle` booleans; the scan
sites skip non-admissible films and LM freezes non-optimized ones. The
surroundings will ride exactly this mechanism (forced `false`, §4.2).

### 2.3 Name-keyed contrast map — DONE
Needle seeds key on host names (`film{i}` / custom). Design-segment
film names are global across envs, so the contrast map needs no
env dimension — surroundings are simply absent from it.

### 2.4 Single monolithic stack per run — the boundary
`run_design(ambient, substrate, films, …)` assembles one flat film list
into one `DesignStack`, solves once, folds once. Scan sites
(`ScanSite.film_idx`) and `NeedleTargets` buckets are **positional**
(0-based within that one stack). There is no segment concept, no
second assembly, and Python cannot inject either into the native
LM/needle loop — so the addition is native, at assembly + driver level.

**Gap (exactly two):** (a) no segmented assembly — one flat list per
run; (b) positional-only parameter identity — nothing maps "the same
film" across two different-length assemblies.

## 3. Proposal

### 3.1 Chosen: named design segments + per-env segment lists

- Global **`design` segments**, one definition each:
  `design: {"coat": [layers…], "back-mirror": [layers…]}`. Film names
  inside design segments are **globally unique and stable** — the cross-
  env parameter identity (§2.4b).
- Each environment defines an ordered **segment list**; a segment is
  either inline fixed layers or a `{design: id}` reference:
  ```json
  "environments": [
    {"name": "bare", "stack": [{"layers": []}, {"design": "coat"},
                               {"layers": [substrate_half]}]},
    {"name": "laminated", "stack": [{"layers": [cover_glass, adhesive]},
                                    {"design": "coat"},
                                    {"layers": [substrate_half]}]}
  ]
  ```
  Full layer list per env = concatenation. Halves stay the outer layers
  (thickness-0 convention) — so the halves-only env of the first draft
  is the degenerate case (empty fixed segments), subsumed, not separate.
- **Shared parameter vector** = design films with `optimize=true`;
  fixed-segment films are assembled with `optimize=false,
  needle=false` **forced** (`true` there refuses, named — silent
  freezing would lie about the optimum).
- Per eval: K assemblies + K solves → each demand evaluates against
  its tagged env's curves (`environment=` on all three demand kinds,
  `None` = first env) → residuals concatenate env-major,
  insertion-minor (deterministic).
- Needle insertion locus = **(design segment, intra-segment position)**,
  translated to absolute indices per env. The design object is single,
  so every insertion propagates to all envs automatically. Fold
  deposits route **by design-film name** into the shared buckets and
  sum across envs; surroundings never deposit (inadmissible films).
  LM steps the shared thickness vector against the joint residuals —
  unchanged machinery.

### 3.2 Rejected alternatives

- **(a) Sequential per-env runs** — oscillates, never joint. Rejected.
- **(b) User-space joint merit** — LM + needle live behind the native
  boundary. Rejected.
- **(c) Halves-only environments** (first draft of this file) —
  cannot express different cover/substrate *sequences*; a different
  ambient scalar is not a different cover glass. Subsumed as the
  degenerate case, not the model. Rejected as the scope.
- **(d) Full-layer-list per env with design-by-flag** (each env repeats
  all layers, design films marked `design: true`) — repeats the shared
  object K times, invites accidental divergence (different initial d
  per env = ill-posed shared vector), and needs consistency validation
  that the segment model gets for free (one definition). Rejected —
  single definition is the invariant.
- **(e) Per-env grids** — K solver contexts for an unstated use case.
  Shared grid stays a §6 limitation.

## 4. Detailed design (minimal additions only)

### 4.1 Schema (additive-optional → old calls parse unchanged)

- `DesignSegmentJson { layers: [LayerJson…] }` under `"design": {}`
  (absent = today's flat `films` behavior, bit-identical path).
- `EnvSegmentJson = {layers: […]} | {design: "<id>"}`;
  `EnvironmentJson { name, stack: [EnvSegmentJson] }`.
- `TargetSet.environments` + per-demand `environment: Option<String>`
  on all three kinds (same as first draft — unchanged).
- Fixed-segment films are **auto-named** per env+position
  (`{env}.fixed[{seg}][{i}]`); user names there are refused on
  collision with any design name, ignored otherwise (they address
  nothing — no buckets, no contrast keys). Design-name duplicates
  across segments refuse at compile, named.

### 4.2 Compile & validation (all native, all named)

- Every env must reference **each design segment exactly once** (v1;
  omission/double-ref refuse naming env + segment — §7.4 records the
  relaxation).
- `optimize/needle: true` inside a fixed segment refuses
  (`"<env>.fixed[…]: optimize/needle must be false — surroundings are
  not design variables"`).
- Unknown `environment=` on a demand refuses naming the demand +
  known env names (first-draft rule, kept).
- Store `env_idx` on compiled demands; store the per-env
  absolute-index map `design_slot → (env, abs_idx)` for fold routing.

### 4.3 Eval (assembly + driver loop only)

- Assembler concatenates per env: fixed layers verbatim (flags forced
  false) + the live design objects (current thicknesses — one object,
  read K times). K `DesignStack`s → K solves per eval.
- New thin `residuals_multi(&[SimCurves])` routing each demand to
  `sims[d.env_idx]`; single-env delegates with `&[sim]` (bitwise
  today's path). Missing-curve penalty per `(env, key)` group. No new
  `NREQ_*`, no merit-formula change.

### 4.4 Needle + LM (routing only, no new kernels)

- Scan sites built per env over the full stack, filtered to design
  films; candidate locus translated back to (segment, intra-position)
  so one insertion edits the shared object.
- Fold per env with existing arms (incl. color `grad_r/grad_t` +
  mirror kernels), deposits mapped via `design_slot` into shared
  buckets and summed; assert cross-env bucket alignment (same design
  object ⇒ same layout — the assert is the proof).
- Contrast map, insertion physics, thickness-LM, defrag: untouched.

### 4.5 Python surface (thin, per binding principles)

- `run_needle(layers | design={…}, …, environments=[…])`; old
  flat-films + `ambient/substrate` kwargs keep working (singleton
  default env, bitwise path).
- Demand `environment=` tags pass through `_dump`; native
  `__post_init__` refuses unknown names (no Python pre-checks).
- Program schema: `design:` + `environments:` sections; file-first
  refs, names resolved in the live context like materials/groups.

### 4.6 Integration & speed (unified path — assess & verify)

DECISION: one unified path, not a fork. Dual maintenance is rejected
by house rule, and it is unnecessary: every added cost sits outside
the hot path or behind a single predictable branch.

- **Compile time** (run start): name registry, exact-once validation,
  `env_idx` tagging, the `design_slot → (env, abs_idx)` routing table.
  Validated once, never consulted mid-eval.
- **Assembly time** (run start + after each needle insertion — NOT per
  eval): K concatenations of fixed segments + the live design object.
  Fixed surroundings assemble **once** (they never change); per eval
  only the design films are re-read. Concatenation is O(layers)
  pointer copies against an O(angles × wavelengths × layers) solve:
  noise.
- **Untouched, by contract:** the solver, `SimCurves` layout, merit
  formulas, fold arms, `SimCurves::back` routing. None learns that
  environments exist. Added struct fields (`env_idx`, segment
  metadata) live in cold structs the hot loops never read.
- **The one branch:** per eval the driver checks `K == 1` vs `K > 1`
  exactly once, *outside* all loops. K=1 calls the existing functions
  with existing signatures (`residuals(&sim)`, current fold) — same
  call sequence op-for-op. K>1 calls `residuals_multi(&[sims])` +
  per-env folds + bucket sum (one slice index per demand, O(params)
  adds — dwarfed by the K solves, which are inherent, not overhead).
- **Why a runtime branch, not const-generics:** K arrives as runtime
  data (JSON), so monomorphization would need a match-dispatch to
  pick the instantiation — the same branch plus code bloat. One
  predictable branch per eval retires at ~0 cycles.

VERIFY (gates, not claims — wired into S1/S2, §7):

- (a) **Bitwise twin:** K=1-segmented binary vs pre-seg binary, same
  seeds — proves the untouched call sequence. Any diff fails S2.
- (b) **Bench gate:** K=1-segmented vs pre-seg on the existing eval
  benchmarks (`validation/benches/…`), threshold **<1%** — proves the
  branch costs nothing. Any excess fails S2 and means the phase is
  over-scoped (new code leaked into the hot path).
- (c) **Scaling check:** K=2-identical-surroundings merit ≡ K=1
  exactly (proves routing adds nothing but the second solve); K=2
  wall-time ≈ 2× single-solve ± concatenation noise (proves no
  hidden per-env fixed cost).

## 5. What is possible day one (coverage matrix)

| Demand kind | × K surroundings | × front/back | × N angles | × M illuminants | × band |
|---|---|---|---|---|---|
| Spectral | **new** | already ✓ | already ✓ | n/a | pointwise |
| Angular | **new** | already ✓ | n/a | n/a | pointwise |
| Color | **new** | refused (unchanged) | already ✓ | already ✓ | already ✓ |

U1–U3 for all kinds (color×back still refused); U1+backside composed
(e.g. air-front + laminated-back spectral) falls out with zero extra
code; the first draft's halves-only cases work as empty-fixed-segment
envs.

## 6. Required limitations (non-negotiable v1)

1. **Shared grid.** One wavelengths + angles for all envs (§3.2e).
2. **One design definition.** Same layer count/materials/initial-d per
   segment across envs — there is only one object (§3.2d).
3. **Surroundings fully fixed.** No env-local free variables: a film
   is either a shared design variable or fixed everywhere. (Per-name
   sharing *of surroundings* is the named follow-up, not a phase.)
4. **Insertions propagate.** Needle cannot add a layer in one env only
   — that would fork the shared object. (This is the point, not a bug;
   state it so nobody files it.)
5. **Each env refs each design exactly once** (§4.2; relax in §7.4).
6. **Cost ×K full solves** per eval and per candidate. Documented
   follow-up: S-matrix embedding — surroundings are fixed, so their
   S-matrices precompute once and the per-eval solve shrinks to the
   design segments via the in-tree Redheffer star product (benchmarks
   exist: `bench` redheffer suites). v1 concatenates; embedding is
   pure speedup, no semantic change.
7. **Residual order env-major** (positional consumers re-check
   bookkeeping; single-env order unchanged).
8. **Color × back refused; opacity composes multiplicatively**
   (paired driver × K envs) — both unchanged/parked as in draft one.

## 7. Phased implementation (each: gate → twin → commit to `dev`)

- **S1 — schema + compile (0.4.32).** Segment model, name registry,
  exact-once rule, flag refusals, unknown-env refusals, `None`→env 0.
  Twins: old flat calls assemble bitwise-identical stacks; each
  refusal names env/segment/demand; K=1 segmented ≡ flat (same
  `DesignStack`, `assert_eq` on films). Bench baseline recorded here
  (pre-seg eval timings — the S2 <1% gate in §4.6 measures against it).
- **S2 — driver + joint merit (0.4.33).** K assemblies/solves,
  `residuals_multi`, per-(env,key) penalties. Twins: K=1 segmented run
  BITWISE vs pre-seg binary (same seeds) — §4.6 gate (a); K=2 identical surroundings ≡
  K=1 merit exactly (§4.6 gate (c) scaling part); bench gate §4.6 (b)
  <1% vs the S1 baseline; bitwise run is §4.6 gate (a); U1 hand oracle
  (single film behind two different
  cover sequences vs numpy transfer-matrix) 1e-12.
- **S3 — needle + LM joint (0.4.34).** Locus translation, name-routed
  fold sum, alignment assert. Twins: K=1 fold BITWISE vs pre-seg;
  K=2-identical fold ≡ K=1; FD with differing surroundings (joint grad
  = sum of per-env analytics) 1e-12; 2-env U1 needle run improves
  both envs monotonically (anti-§3.2a proof); insertion lands in the
  design segment in every env assembly (positions differ, names match).
- **S4 — Python surface + program (0.4.35).** `design=`/`environments=`,
  demand tags, program sections. Twins: surface ≡ JSON HEX; refusals
  at construction (`ValueError`, named); docs + matrix row.
- **S5 — docs + release.** Kinds-doc environment section, worked U1
  example (same AR bare vs laminated), ×K cost note, exposure re-audit.
  Tag per policy.

Estimated total: ~350 Rust lines + ~100 Python lines + tests — assembly
+ routing, no physics. Any kernel/fold-arm/merit-formula diff means the
  phase is over-scoped.

## 8. Open questions (resolve in S1, all one-liners)

1. **Default-env name** (`"default"` explicit vs `None`-only)?
   Recommend named — explicit beats magic (kept from draft one).
2. **Env list home** (call-only vs also `TargetSet.environments`,
   JSON-wins-on-clash)? Recommend both (kept from draft one).
3. **Per-env `missing_penalty` scale?** Recommend one floor for all
   (kept from draft one).
4. **Exact-once relaxation:** allow envs to omit a design segment
   (it just doesn't contribute there)? Recommend refuse in v1 —
   omission today is almost certainly a typo'd env, and silent
   non-contribution would lie about jointness. Relax with a loud
   `contributes: false` marker if a real case appears.
5. **Fixed-film user names:** auto-only, or allow non-colliding custom
   names (for callbacks/debugging)? Recommend auto-only v1 — names
   that address nothing invite `contrast` entries that silently never
   split.
