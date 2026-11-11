# Synthesis: open items, inconsistencies and decisions

Status 2026-11-11, branch `needle-material-pool`, after W6 (first cut) and
the integral phase fold. Sources: `synthesis_targets_plan.md` (W1-W8, §4-§7),
the W6 work and a read of the code paths it touched. Each item says where
it lives, what is wrong or missing, and what would close it.

Severity: **bug** = a wrong answer today; **gap** = a refusal or a missing
feature, said out loud; **debt** = code, docs or pins out of step.

---

## 1. Bugs and inconsistencies found

### 1.1 `p_function_multiblock` misaligns its output under a host mask (bug)

`needle_operator.rs`, `locate_hosts_multiblock` + `p_function_multiblock`.
The locator `continue`s past a masked z without pushing an entry, so
`locs` is shorter than `z_grid`. `p_multiblock_point` returns one value
per `loc`, and `p_function_multiblock` adds `contrib[zi]` into `p_out[zi]`:
every z after the first masked one gets the P of a later site, and the
tail stays 0.

- Failure: `z_grid = [z_a (masked), z_b, z_c]` returns `[P(z_b), P(z_c), 0]`.
- Reach: the public Rust function only. It is not bound to Python, and
  synthesis does not call it (W6 uses its own site-to-host map in
  `needle_pass.rs`). Its tests pass no mask.
- Fix: carry the z index in `locs` (or push a sentinel) and scatter by it,
  and add a masked case to the multiblock FD test.

### 1.2 Linux needle-pin digest predates W1 (debt)

`validation/regression/synthesis/test_needle_pin.py` holds a Linux digest
from before the W1 fold; Windows was re-recorded at W1. W4, W5, W7, the
integral fold and W6 did not move the Windows digest. The Linux value
must be re-read from CI (`ubuntu-24.04`, release) and recorded with
provenance, or the Linux job fails on the first CI run of this branch.

### 1.3 Two notions of "total thickness" (debt, W6)

`DesignStack::total_thickness_nm` sums every film;
`coating_thickness_nm` (new) leaves out incoherent films. The thickness
budget, the pipeline's `total_thickness_nm` and `final_total_thickness_nm`
and the inflate report use the coating sum (equal bit for bit on coherent
stacks). The Python `DesignStack.total_thickness()` (and the stack's
`repr`) still call `total_thickness_nm`, so they count a plate. Decide whether it should be renamed, or the
report keys say "coating". Decided (D7): keep both, name both; not yet
done.

### 1.4 The PD reference thickness counts an incoherent plate (latent)

`simulate_inner` takes `D = Σ thicknesses[1..nl-1]` for differential
phase. On a stack with a plate this would include it. It cannot be
reached today: every phase demand is refused on an incoherent stack.
It matters if phase demands inside the first block are ever allowed
(§3, W6 follow-up 3).

### 1.5 Plan text out of step (debt)

- `synthesis_targets_plan.md` §1.6 still says construction refuses
  `coherent = false`; W6 lifted that for fixed films. Updated in this
  change, together with the W6 status note.
- §1.3b's row "any channel through an incoherent layer: blocked by W6"
  and §3's step 9: updated to "first cut, intensity only".
- §6 follow-ups (target-kinds doc, needle docs with a generated
  supported-target matrix) have not been done for any workstream.

### 1.6 A doc comment sits on the wrong item (debt)

`needle_operator.rs`, around `PmbQuantity`: the doc comment and the
`#[allow(clippy::too_many_arguments)]` meant for `p_multiblock_point`
sit above `enum PmbQuantity`. So `PmbQuantity`'s docs open with
`p_multiblock_point`'s description, and `p_multiblock_point` has no
docs. It is harmless at runtime. Fix: move both lines down to the
function.

### 1.7 Inflate and clamp touched pinned films (fixed, D10)

`inflate.rs` picked every single-layer span as a candidate whatever its
`optimize` flag, and the clamp sweep removed, clamped up or capped a
pinned film like any other (only cleanup's removal checked the flag).
D10: a pinned film is untouched. Inflate and `round_to_qwot` skip it,
and `clamp_all_policy` passes a pinned span (every bulk row
`optimize = false`) through unchanged, below the floor or above the
ceiling, and `merge_adjacent` no longer fuses a pinned film with a free
one of the same index (interface slices, whose flag is derived, still
merge). That also covers incoherent plates (always pinned) and fixed
surround films in a joint run; the design's pinned copies follow their
primary through `sync`/`expand` as before. A pinned film that hosts
needles (`needle = true`) can still be split by an insertion; its
halves stay pinned, so a sliver left by an insertion is kept.

---

## 2. Deferred and not-done items, by workstream

### W3, capability check
- **Item 2 not done.** `merit()` still adds the 1e6 missing-curve
  penalty silently; D2 decided to make it a hard error.

### W5, transforms in the needle
- `complex` transform: refused at the first fold (no fold form exists).
- Ratio, extinction and contrast targets: not started; would reuse the
  tangent-pair construction.

### W6, incoherent synthesis (first cut, this change)
Done: fixed incoherent films anywhere in the film list; forward through
the engine's block cascade (front rows bit-identical to `solve_point`);
needle P through the cascade adjoint for every intensity and colour
bucket (FD oracle, 0° and 30°); LM by finite differences; budget,
clamp, merge and inflate leave the plate alone. Not done:

1. **Analytic Jacobian through the cascade.** LM differences on any
   stack with a plate (n_par extra solves per iteration). The adjoint
   weights exist (`CascadeTrack.g`); deposits per block times them would
   close it.
2. **Multi-environment runs: done.** `simulate_all` no longer refuses
   (its check only ever looked at environment 0's assembly, so a plate
   in another environment's fixed segment already got through). Each
   assembly solves on its own; the joint needle sums each environment's
   slope through its own cascade. Checked: the joint-slope oracle
   (analytic sum vs a finite difference of the joint merit) with the
   coat in front of a plate, behind it, and once on each face (a copy in
   each block); joint runs with interior and boundary needles in those
   layouts plus an inverted copy on the back face; a plate as a pinned
   film of the design segment itself, routed into every environment and
   left at its thickness; a Python joint run with the coat on both faces
   of a plate. Phase demands in such runs: D4 step 1, done (§3.1).
3. **Phase demands** on a stack with a spacer are dropped from that
   environment with a warning (D4 step 1), including phases that live
   entirely in the first block (front `r`) or the last one. Keeping
   those single-block phases is D4's open remainder (needs 1.4).
4. **An uncoated face cannot grow a coating.** The plate is never a host,
   and a boundary site with no coherent film above it gets P = 0, so a
   face with no film on it stays bare. A face that should be coated
   needs at least one seed film. Documented here; not yet in user docs.
5. **No explicit exit medium.** The plan proposed
   `DesignStack.exit: Option<LayerSpec>` with the substrate as a finite
   incoherent layer. Implemented instead as "the plate is a film, the
   exit medium is the substrate". Same physics, no new field; the
   Python doors already pass `coherent` per film. Record it as the
   adopted model or schedule the field (decision D5).
6. **Edge materials** (`outer`/`inner`) next to a plate: the edge rule
   names the first/last film, which may now be the plate. Not tested.
7. **Roughness at a plate face**: carried to the solver unchanged, as the
   engine does; no synthesis test covers a rough plate interface.
8. **Back-incidence from an absorbing plate**: `check_back_exit` looks at
   the exit medium only. A strongly absorbing plate makes back demands
   ill-conditioned (tiny T_b) but not undefined; no guard.
9. **Test coverage**: the needle oracle runs through a plate at 0° and
   30°, not 60°. No test covers two plates, a plate as the first or last
   film, a plate with only back films, Log or weighted demands through a
   plate. The joint needle path is covered since W6 item 2.

### W7, material identity
- `pool_share_tol` (near-equal indices sharing a sweep): not done.
- Shared-sweep and merge **reports** (which entries shared a sweep,
  which films fused): not done.
- Both are off the to-do list until a user asks (D9).
- The joint (F2) scan sweeps each entry separately (its boundary loci
  differ per entry); only the flat scan shares sweeps.

### W8, new target types
Done: W8a retardance (RetR/RetT, DeltaR/DeltaT). Open, in plan order:
1. Rs/Rp ratio, extinction ratio, diattenuation (tangent pair or grad
   bucket).
2. Psi_R/Psi_T: needs a complex gradient operator
   (`p_coherent_grad_c_from_fields`).
3. DOP / Stokes for a given input state: needs an input-polarization
   field on the target.
4. GD…FOD: port the engine's dispersion forward into `simulate_inner`,
   add curve ids, fold onto the existing Q-row sensitivities; LM by FD
   first.
5. PDrs/PDrp: deferred (reference plane of the round trip, D3).

### Further gaps (plan §5)
- Angle-cone and bandwidth-averaged targets (W4 point weights are the
  building block).
- Colour on back incidence (compile refusal); colour with a spectral
  emphasis profile.
- Weight profiles per environment.
- A randomized standing property test P = dM/dδ (the oracle covers fixed
  cases only).
- General-fit interplay: group `n_factor` as a fit parameter must keep
  W7 labels stable.

---

## 3. Decisions required

Decided 2026-11-11 unless marked open.

| # | Question | Decision |
|---|---|---|
| D1 | First fold of a run: static (`from_spec` folds with no simulation) or an initial simulation? | **Keep static, but no silent fallback.** Refold is on by default and refolds every cycle, the first included, so the static fold only serves `refold_per_cycle = false`. The weak spot is `cycle.rs`'s `if let Ok(sim) … && let Ok(f) …`: a failed refold falls back to the previous or static fold without a word. Make that failure an error (or at least a run-report warning). |
| D2 | `merit()`'s 1e6 missing-curve penalty after W3. | **Hard error** (or penalty plus a warning naming the demand). Since W3 a missing curve can only come from a spec/evaluator mismatch, i.e. a bug; it is added silently today (`merit.rs` `merit_multi`). |
| D3 | PDrs/PDrp round-trip reference plane. | **Deferred** until a user case. |
| D4 | W6: phase demands on a stack that loses phase (an incoherent plate). | **Warn, never silently scramble**; let a phase demand live on a phase-keeping environment in parallel with a phase-losing one. Step 1 done (§3.1); single-block phase on a plate stack still open. |
| D5 | W6 model: plate as an incoherent film vs `DesignStack.exit`. | **Keep the current model.** The exit medium is whichever half-space row ends an environment (`layer_type` 2); in multi-environment runs fixed segments may carry further media behind a substrate, so nothing inside a stack has to be marked "exit". Possible later hard rename: `DesignStack.substrate` → `exit`. |
| D6 | W6 follow-ups: Jacobian vs multi-environment. | **Multi-environment support first, then the analytic Jacobian through the cascade.** |
| D7 | `total_thickness_nm` vs `coating_thickness_nm` (1.3). | **Keep both, name both.** Python gets `DesignStack.coating_thickness()` beside `total_thickness()`; the run-result keys are renamed to say coating (`final_coating_thickness_nm`, per-phase `coating_thickness_nm`), hard rename. |
| D8 | Fix 1.1 now, or refuse `host_mask`? | Fix (small); not yet done. |
| D9 | W7 `pool_share_tol` and the sweep/merge reports. | **Off the to-do list until asked** (no code is removed; neither affects results). |
| D10 | Pinned films (`optimize = false`) under inflate and the clamp sweep. | **Untouched.** Inflate and `round_to_qwot` skip them; the clamp sweep neither removes, clamps up, caps nor refuses them. Done (see 1.7). |

### 3.1 D4 proposal: phase demands and phase-losing environments

The problem: an incoherent plate (or any spacer) scrambles the phase of
everything that crosses it. A phase demand on such a stack is either
undefined (transmission through the plate) or defined only for one
coherent block (front reflection off the first block, back reflection
off the last). Today every phase demand on such a stack is refused.

Proposal, in two steps:

1. **Per-environment routing, with a warning.** A phase demand names the
   environments it applies to (default: all). At compile time each
   (demand, environment) pair is classified:
   - *phase kept*: the environment is fully coherent, or the demand
     lives entirely in one coherent block (front `r` of the first block,
     back `r` of the last) - the demand is evaluated there;
   - *phase lost*: the light path crosses a spacer - the pair is
     **dropped from the merit**, and the run reports a warning naming
     the demand and the environment ("phase of RetT is not defined
     through the incoherent plate in 'with_plate'; dropped from this
     environment").
   A demand that ends up dropped from **every** environment is a hard
   error, not a warning: nothing would optimize it.
2. **The intended use** is then a two-environment run: the same design
   segment once on a coherent surround (e.g. the coating on a
   semi-infinite substrate, where phase is defined) carrying the phase
   demands, and once on the real plate carrying the intensity demands.
   W6 item 2 (multi-environment with a plate) is in place, so this
   layout runs today for intensity demands; step 1 is what lets the
   phase demands join it.

**Step 1 as built (2026-11-11).** `run_environments` drops every phase
demand of an environment whose assembly holds an incoherent film
(`MeritSpec::drop_phase_demands`) and returns one warning naming those
environments; the Python door raises it as a `UserWarning`. Each
environment then builds complex rows only for its own demands
(`SmatrixContext::simulate_env`, also used by the per-cycle refold), so
the plate environment solves through its cascade while the coherent one
keeps its phase rows. Two simplifications against the proposal above:
- Every phase demand in a plate environment is dropped, single-block
  ones (front `r`) included; keeping those is the open remainder.
- The hard error is run-level: the run refuses when every phase demand
  was dropped. A compiled spec holds one demand per environment and no
  longer knows which ones came from the same user target, so "this
  target was dropped everywhere" cannot be told apart from "this target
  only ever named the plate environment" while another phase target
  still binds elsewhere; that case warns.
The flat path (`run_design`, one environment) keeps its own refusal.

Single-block phase on a plate stack (front `r` of the first block)
needs 1.4 fixed first (the differential-phase reference thickness must
count only that block).

Plan open question 4 (coated back side) is answered by W6: both faces
can be coated, provided each starts with at least one film.

---

## 4. Order of work

1. D10, pinned films untouched (done).
2. W6 item 2: multi-environment runs with an incoherent film (D6) (done).
3. D4: per-environment phase routing with warnings (§3.1) (step 1 done).
4. W6 item 1: analytic Jacobian through the cascade (D6).
5. Small items whenever convenient: D1 (refold failure is an error),
   D2 (missing curve is an error), D7 (coating-thickness names), 1.1 /
   D8, 1.6.

## 5. Housekeeping before merge

- Re-read the Linux needle-pin digest from CI (1.2).
- CHANGELOG entries for W7, the integral fold and W6 are in place;
  release ritual (7 version sites, fingerprints) not yet run.
- User docs: the supported-target matrix, the "seed a film on each face"
  rule for plates (W6 item 4), and the refused combination (phase
  through a plate, and any phase demand in a joint run with a plate).
- Python door fix with W6 item 2: `run_needle(environments=…)` added
  the run's `ambient`/`substrate` row to every environment, so an
  environment with its own `layer_type` 0/2 row was refused ("at most
  one substrate"). It now adds the default only where the environment
  has none, as its comment always said.
