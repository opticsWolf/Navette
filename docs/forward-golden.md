# Forward-model golden baseline

A regression baseline for every spectral and phase quantity Navette
computes, built to guard performance work on the solver (reuse across
optimisation loops) against any numeric drift. It lives in
`validation/regression/forward/`:

| File | Role |
|---|---|
| `_cases.py` | materials, grids, the 34 cases, evaluator, pack/unpack, tolerances |
| `golden/forward_golden.npz` | the pinned arrays (798 stacked members, ~3 MB) |
| `golden/manifest.json` | format version, Navette version + git commit, grids, case definitions, request masks, block layouts, file SHA-256 |
| `test_forward_golden.py` | recomputes every case, compares, reports the worst deviation per quantity |
| `test_forward_physics.py` | independent references (Fresnel, Airy, energy, identities, phase continuity, door agreement) and the known defects as strict xfails |
| `gen_forward_golden.py` | the only writer: dry run by default, `--write` to rewrite |

Baseline written at Navette 0.8.1, commit `b24b06c`; regenerated at
`22cd32e` (defect 2) and `4c2d969` (defect 3, section 4).

## 1. Survey: what can be computed

### 1.1 Engine observables (`ScatterMatrix.compute` / native `Solver.solve`)

All 51 request bits of `core_engine.rs` / `navette.smatrix.Request`. Arrays
are `[n_angles, n_wavs]`. Every one is pinned in every case where the door
allows it (the cross channel is refused under `FRONT_BLOCK` on a stack with
an interior incoherent layer, C2).

| Family | Keys (bit) | Convention, as verified | Independent check |
|---|---|---|---|
| Reflectance, transmittance | `Rs` `Rp` `Ts` `Tp` (0-3) | T is the Poynting flux ratio, `Re(q_exit)/Re(q_0)` (s), `Re(n*_exit cos)/Re(n*_0 cos)` (p) | Fresnel, Airy |
| Unpolarised | `R_avg` `T_avg` (4-5) | `(s + p) / 2` | identity |
| Absorptance | `A_s` `A_p` `A_avg` (6-8) | `1 - R - T`, unclamped (Névot-Croce can make it negative) | identity, sign |
| Ellipsometry | `Psi_R` `Psi_T` (9-10), `Delta_R` `Delta_T` (11-12) | radians; `tan Psi = sqrt(Rp/Rs)`; `Delta_R = arg(-rp/rs)`, `Delta_T = arg(tp/ts)` | amplitudes |
| Polarisation | `DOP_R` `DOP_T` (13-14), `Diattenuation_R/T` (15-16), `Retardance_R/T` (35-36) | 45-degree-input Stokes; `D = (Rp-Rs)/(Rp+Rs)`; `Retardance_R = Delta_R - pi`, `Retardance_T = Delta_T` | identities |
| Stokes | `S0..S3_R` (17-20), `S0..S3_T` (21-24) | `S0 = Rs+Rp`, `S1 = Rp-Rs`, `S2_R = -2 Re(cross_R)`, `S3_R = -2 Im(cross_R)`, `S2_T = 2 Re(cross_T)`, `S3_T = 2 Im(cross_T)` | pure-state identity |
| Absolute phase | `phi_rs` `phi_rp` `phi_ts` `phi_tp` (25-28) | `arg` of the amplitude, principal value | Airy, continuity |
| Complex amplitudes | `rs_c` `rp_c` `ts_c` `tp_c` (29-32) | exp(+i beta); `rp_c` = minus textbook r_p; `tp_c` = textbook t_p x cos(theta_exit)/cos(theta_0) | Fresnel, Airy |
| Coherency | `cross_R` `cross_T` (33-34) | `rp rs*`; `tp ts*` scaled by the real `sqrt(f_s f_p)`, so `abs(cross_T)^2 = Ts Tp` (0 without flux) | amplitudes, echo sum |
| Dispersion | `GD/GDD/TOD/FOD_{R,T}_{s,p}` (37-40) | fs, fs^2, ...; `GD = +d phi / d omega` | finite difference |
| Back incidence | `phi_rbs/rbp/tbs/tbp` (41-44), `rbs_c ... tbp_c` (45-48) | phase = `arg` of the back amplitude | identity |
| Differential phase | `PDts` `PDtp` (49-50) | `arg(t * exp(-i 2 pi n0 D cos(theta0) / lambda))`, D = interior thickness | reference rotation |

`energy_conservation()` is `max(|1-Rs-Ts|, |1-Rp-Tp|)` over pinned
intensities, so it is covered implicitly.

### 1.2 Synthesis forward model (`SmatrixContext.simulate`)

A second solver: the dual coherent-block solve that every optimisation loop
calls. It produces `Rs Rp Ts Tp` and the complex forward `t` for phase and
differential-phase demands. `SimCurves` has no getters, so the baseline
reads them back through `MeritSpec.residuals` with neutral targets
(target 0, tolerance 1), one context per angle: `sim_Rs/Rp/Ts/Tp`,
`sim_phi_ts/tp`, `sim_PDts/tp`. Pinned for every single-block case.

### 1.3 Layer features

| Feature | Parameters | Where it acts | Cases |
|---|---|---|---|
| Absorbing film / substrate | complex `n` | engine | `film_tio2_abs`, `film_ag`, `mdm_ag`, `bare_si`, `bare_ag` |
| Total / frustrated internal reflection | dense ambient | engine | `bare_tir`, `film_ftir` |
| Interface roughness | `rough_type` 1 LINEAR, 2 STEP, 3 EXPONENTIAL, 4 GAUSSIAN, 5 NEVOT_CROCE; `roughness` sigma | engine (per interface) | `rough_*`, `rough_mixed_ml` |
| Roughness cap | `roughness_cap` (default 2, 0 off) | engine | `rough_cap_binds`, `rough_cap_off` |
| Inhomogeneity (graded index) | `inhomogen`, `inh_delta`; `inh_mode` Fixed / RateCapped | structure expansion into sub-layers | `inh_fixed`, `inh_ratecapped`, `inh_multilayer` |
| Mixture gradient (EMA) | `gradient`: Bruggeman / Looyenga / MaxwellGarnett; FixedSpan / RateCapped; `sublayers` | structure expansion | `gradient_*` |
| Interface slices | `interface`, `interface_thickness` | structure expansion | `interface_slices` |
| Group scaling | thickness factor/summand, n/k factor, inh/roughness/interface summands | structure expansion | `group_scaled`, `combined_*` |
| Incoherent (thick) layer | `coherent=False` | engine | `incoh_*`, `combined_*` |
| Coherence mode | 0 FRONT_BLOCK, 1 COHERENCY_MATRIX, 2 FULLY_COHERENT | engine | `incoh_substrate_fb/cm`, `fully_coherent_thick`, `combined_cm/fb` |

The expanded solver rows (`indices`, `thicknesses`, flags, roughness) are
pinned too, so a drift in the structure layer is reported as such and is
not confused with one in the solver.

### 1.4 Not covered (out of scope for a spectral/phase baseline)

Guided-mode tools (`eigenmode_landscape`, `refine_mode`,
`find_eigenmodes`, `field_profile`); needle gradients and the needle/LM
pipeline (pinned by `regression/synthesis/test_needle_pin.py` and the bench
fingerprint); merit values; error sampling (`get_error_solver_inputs`,
which is random). The raw `core_engine` pyfunction runs the same kernels as
`Solver` and is compared against numba in `parity/smatrix/`.

## 2. Cases, grids, tolerances

* **Wavelength ranges**, solved separately: UV 250-400, VIS 400-800,
  NIR 800-2000 nm, 8 points each (`np.linspace`).
* **Angles**: 0, 15, 30, 45, 60, 75, 89 degrees, passed as `angles_deg=`.
* **Materials**: closed-form dispersion (Cauchy dielectrics, a TiO2 with a
  UV absorption tail, a Drude silver, a dispersive Si), so every channel is
  wavelength dependent.
* **34 cases**: bare substrates, single films, quarter-wave, broadband AR,
  metal-dielectric and 21-layer edge stacks, each feature alone, and two
  stacks with every feature at once (`combined_cm`, `combined_fb`).

| Quantity | Tolerance |
|---|---|
| intensities, amplitudes, Stokes, coherency, solver rows | 1e-12 absolute |
| phases, Delta, retardance, PD (wrapped difference), Psi, DOP, diattenuation | 1e-9 absolute |
| GD, GDD, TOD, FOD | 1e-8 relative to the channel's largest value |

The looser groups are ratios that become ill-conditioned when the signal
vanishes (the phase of a near-zero amplitude, the diattenuation of two
near-zero reflectances) or finite differences that amplify rounding.
At the time of writing the baseline reproduces **bit for bit**;
`NAVETTE_GOLDEN_EXACT=1` asserts that. `NAVETTE_GOLDEN_REPORT=1` prints the
worst deviation per quantity across all cases. A failing case lists the
worst deviation per quantity, with the offending ones marked.

## 3. Regenerating

```
python validation/regression/forward/gen_forward_golden.py          # dry run, writes nothing
python validation/regression/forward/gen_forward_golden.py --write  # rewrite baseline + manifest
```

Rewrite only for an intended change to the forward model, and say why in
the commit message. A performance change must pass against the baseline
unchanged. The test refuses a manifest whose case definitions, grids or file
hash do not match the code.

## 4. Known defects found while building it

The baseline pins today's numbers, defects included: its job is to catch
drift. Each defect is a strict xfail in `test_forward_physics.py`, so a fix
flips it. Drop the mark and regenerate in the same change.

1. **Multi-angle synthesis reads the wrong angle (fixed).** `SmatrixContext`
   (and the `run_design` driver) gave `SimCurves` sin(theta) as its angle
   axis, but `MeritSpec` keys are in degrees and `SimCurves::angle_row`
   compares the two directly. With more than one angle, every key at a
   non-zero angle scored against the last angle's curves. Keys and the
   axis are now both `Radians` (`docs/plans/angle_units_plan.md`); the
   solver still runs on sin(theta).
2. **Synthesis T into an absorbing substrate (fixed).**
   `SmatrixContext.simulate`'s `Ts`/`Tp` disagreed with the Poynting
   transmittance the engine returns: T = 84 on bare Ag at 400 nm, 1e-5 off
   on Si. The evaluator took the backward channel, |t_back|²·Re(y_0)/Re(y_exit),
   which matches the forward one only when both half-spaces are
   transparent; into an absorbing substrate it divides by a near-zero
   Re(y_exit). It now takes the forward |t_fwd|²·Re(y_exit)/Re(y_0), like the
   engine and the needle gradient, and so does its thickness gradient. On
   transparent substrates every `sim_T*` moved by up to 1.1e-14, so the
   baseline was regenerated.
3. **Transmitted Stokes S2/S3 lacked the exit-medium factor (fixed).**
   `S2_T`/`S3_T` came from `tp ts*` while `S0_T`/`S1_T` carry the
   admittance ratio, so `DOP_T` of a coherent stack was below 1 whenever
   the exit index differed from the ambient (0.658 for a film on glass at
   normal incidence, 0.19-0.28 on bare Si). Each block's transmitted cross
   term is now scaled by the real `sqrt(f_s f_p)` at its source, forward
   and (in `COHERENCY_MATRIX`) backward, so `abs(cross_T)^2 = Ts Tp` for
   one block and the phase is kept; into an absorbing exit `f_s != f_p`,
   so no single factor on the fields would do. Without flux (past TIR)
   the term is exactly +0: `Retardance_T` reads 0 there instead of the
   evanescent field's phase, and `DOP_T` the 0 convention for an undefined
   state. `Psi_T` and `Delta_T` did not move beyond rounding; the
   `COHERENCY_MATRIX` reflection echoes moved by 1e-16. The baseline was
   regenerated.
