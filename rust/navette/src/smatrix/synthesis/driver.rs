// SPDX-License-Identifier: LGPL-3.0-or-later
//! synthesis::driver — array-level design assembly + end-to-end runs.
//!
//! Rust-first port of the Python `stack_from_layers`/`run_needle`
//! orchestration: evaluated nk arrays + flag structs in, expanded
//! `DesignStack` or full `PipelineResult` out. Config-file assembly
//! (`DesignRequest`) lives in `design_config`; both converge on
//! `DesignStack::from_design` here through [`assemble_stack`].

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use num_complex::Complex64;

use super::config::PipelineConfig;
use super::cycle::{
    MaterialPool, NeedleCycleConfig, NeedleRules, label_pool, pool_entry, pool_from_layers,
};
use super::evaluator::SmatrixContext;
use super::merit::{MeritSpec, Radians};
use super::pipeline::SpectralInputs;
use super::pipeline::{NeedlePipeline, PipelineResult};
use super::structure::{ClampReport, DesignStack, LayerSpec};
use super::thick_opt::LmConfig;
use crate::structure::{Group, Layer, LayerType};

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------

/// One film: evaluated nk on the grid + authoring flags.
#[derive(Clone, Debug)]
pub struct ArrayFilm {
    pub name: String,
    pub nk: Vec<Complex64>,
    pub d_nm: f64,
    pub coherent: bool,
    pub roughness: f64,
    pub rough_type: i32,
    pub inhomogen: bool,
    pub inh_delta: f64,
    pub interface: bool,
    pub interface_thickness: f64,
    pub optimize: bool,
    pub needle: bool,
    /// Mixture-gradient profile (F1.1, A5). The JSON carries the
    /// inclusion spectrum (`nk_b`) Python-side; the host is expected to
    /// name the film itself, whose nk is already registered.
    pub gradient: Option<crate::structure::gradient::GradientJson>,
}

/// One needle pool material: a label and its evaluated nk on the grid.
#[derive(Clone, Debug)]
pub struct ArrayMaterial {
    pub name: String,
    pub nk: Vec<Complex64>,
}

// ---------------------------------------------------------------------------
// Assembly + runs
// ---------------------------------------------------------------------------

fn fixed_half(name: &str, nk: Vec<Complex64>) -> LayerSpec {
    LayerSpec {
        material: Arc::from(name),
        nk: Arc::from(nk),
        d_nm: 0.0,
        coherent: true,
        rough_type: 0,
        rough_val: 0.0,
        optimize: false,
        needle: false,
        label: None,
    }
}

/// Expand evaluated arrays into a `DesignStack` (single `from_design`).
/// Returns `(stack, warnings)` — warnings re-emit upstream.
#[allow(clippy::too_many_arguments)]
pub fn assemble_stack(
    ambient_name: &str,
    ambient_nk: Vec<Complex64>,
    substrate_name: &str,
    substrate_nk: Vec<Complex64>,
    films: &[ArrayFilm],
    groups: &HashMap<String, Group>,
    wavelengths: &[f64],
) -> Result<(DesignStack, Vec<String>), String> {
    use crate::structure::RoughnessType;
    let mut seen = HashSet::new();
    for f in films {
        if !seen.insert(f.name.as_str()) {
            return Err(format!(
                "duplicate film name {:?} (film names key the nk table)",
                f.name
            ));
        }
        if f.nk.len() != wavelengths.len() {
            return Err(format!(
                "film {:?}: nk length {} != {} wavelengths",
                f.name,
                f.nk.len(),
                wavelengths.len()
            ));
        }
    }
    let mut design = Vec::with_capacity(films.len());
    let mut nk_map: HashMap<Arc<str>, Vec<Complex64>> = HashMap::new();
    // The design surface builds its films from flag dicts and never goes
    // through the Python `Layer`, so the layer-construction gate has to be
    // applied here too or this door is simply open. Same rule, same messages.
    let mut film_warnings: Vec<String> = Vec::new();
    for f in films {
        let mut layer = Layer::film(f.d_nm, &f.name);
        layer.layer_type = LayerType::Film;
        layer.coherent = f.coherent;
        layer.roughness = f.roughness;
        layer.rough_type = RoughnessType::try_from_i32(f.rough_type)
            .map_err(|e| format!("film {:?}: {e}", f.name))?;
        layer.inhomogen = f.inhomogen;
        layer.inh_delta = f.inh_delta;
        layer.interface = f.interface;
        layer.interface_thickness = f.interface_thickness;
        layer.optimize = f.optimize;
        layer.needle = f.needle;
        if let Some(grad) = &f.gradient {
            // A5: nk_b must ride the film dict (never a fallback to the
            // film's own nk - a half-mixture is worse than a refusal),
            // on the grid, and the host is expected to name the film
            // itself (whose nk is already registered).
            if grad.nk_b.len() != wavelengths.len() {
                return Err(format!(
                    "film {:?}: gradient nk_b length {} != {} wavelengths",
                    f.name,
                    grad.nk_b.len(),
                    wavelengths.len()
                ));
            }
            if grad.material_a != f.name {
                return Err(format!(
                    "film {:?}: gradient host '{}' is not the film's own name; the design path \
                     resolves the host spectrum through the film's registered nk",
                    f.name, grad.material_a
                ));
            }
            layer.gradient = Some(grad.to_spec());
        }
        let issues = layer.property_issues(&format!("film {:?}", f.name));
        crate::structure::validation::ValidationIssue::gate(&issues, "assemble_stack")?;
        film_warnings.extend(issues.iter().map(|i| i.message.clone()));
        nk_map.insert(Arc::from(f.name.as_str()), f.nk.clone());
        if let Some(grad) = &f.gradient {
            nk_map.insert(Arc::from(grad.material_b.as_str()), grad.nk_b.clone());
        }
        design.push(layer);
    }
    // A gradient's material_b entry must still be ITS spectrum when the
    // loop ends: a later film named like an earlier gradient's inclusion
    // would otherwise shadow the endpoint silently.
    for f in films {
        if let Some(grad) = &f.gradient
            && nk_map.get(grad.material_b.as_str()) != Some(&grad.nk_b)
        {
            return Err(format!(
                "film {:?}: gradient material_b {:?} is shadowed by another film of the same name",
                f.name, grad.material_b
            ));
        }
    }
    // Background is implied, not declared (mirrors the other drivers);
    // A4: a gradient carrier is a profiled film for this rule too.
    let background: HashSet<String> = design
        .iter()
        .filter(|l| (l.inhomogen || l.gradient.is_some()) && !l.optimize && !l.needle)
        .map(|l| l.material.clone())
        .collect();
    DesignStack::from_design(
        fixed_half(ambient_name, ambient_nk),
        fixed_half(substrate_name, substrate_nk),
        &design,
        &nk_map,
        groups,
        wavelengths,
        &background,
    )
    .map(|(mut stack, warnings)| {
        // W7: film names here are free, not material labels.
        stack.clear_labels();
        film_warnings.extend(warnings);
        (stack, film_warnings)
    })
}

/// End-to-end design run: assemble, fold demands, execute the macro-loop.
/// `callback` fires per macro-cycle.
///
/// `materials` is the needle pool; `None` takes one entry per distinct
/// index among the assembled layers ([`pool_from_layers`]). `rules` narrow
/// which needles go into which materials and which materials may form the
/// edge layers ([`NeedleRules`]); their labels are pool labels.
///
/// Returns `(report, stack, warnings)`. The warnings are the assembly's --
/// a dropped ambient absorption (R3.4), a homogenized graded film -- and they
/// must be surfaced upstream. They were dropped on the floor here until
/// 0.6.27, which made this the one design path that corrected inputs in
/// silence. The list also carries the run setup's feasibility drops
/// (phase through an incoherent film, back incidence from an absorbing
/// exit; [`drop_infeasible`]), made once before the pipeline starts.
#[allow(clippy::too_many_arguments)]
pub fn run_design(
    ambient_name: &str,
    ambient_nk: Vec<Complex64>,
    substrate_name: &str,
    substrate_nk: Vec<Complex64>,
    films: &[ArrayFilm],
    groups: &HashMap<String, Group>,
    materials: Option<&[ArrayMaterial]>,
    rules: &NeedleRules,
    wavelengths: &[f64],
    angles: &[Radians],
    spec: &MeritSpec,
    cfg: PipelineConfig,
    needle_cfg: NeedleCycleConfig,
    lm: LmConfig,
    mut callback: impl FnMut(usize, &super::pipeline::PipelinePhaseResult) -> Result<(), String>,
) -> Result<(PipelineResult, DesignStack, Vec<String>), String> {
    SmatrixContext::check_spec(spec)?;
    let (mut stack, mut warnings) = assemble_stack(
        ambient_name,
        ambient_nk,
        substrate_name,
        substrate_nk,
        films,
        groups,
        wavelengths,
    )?;
    let pool = match materials {
        None => pool_from_layers(&stack, |_| true)?,
        Some(ms) => {
            for m in ms {
                if m.nk.len() != wavelengths.len() {
                    return Err(format!(
                        "needle material '{}': nk length {} != {} wavelengths",
                        m.name,
                        m.nk.len(),
                        wavelengths.len()
                    ));
                }
            }
            label_pool(
                ms.iter()
                    .map(|m| pool_entry(&m.name, Arc::from(m.nk.clone()))),
            )?
        }
    };
    let pool = rules.apply(&mut stack, pool)?;
    let (spec, dropped) = drop_infeasible(
        "run_design",
        spec,
        std::slice::from_ref(&stack),
        &[],
        wavelengths,
    )?;
    let spec = &spec;
    warnings.extend(dropped);
    let spectral = SpectralInputs::from_spec(spec, angles, wavelengths)?;
    let sin_theta = spectral.sin_theta.clone();
    let cfg = cfg.validated()?;
    let mut ctx = SmatrixContext {
        wavls: wavelengths.to_vec(),
        sin_theta,
        spec: spec.clone(),
        clamp_min_nm: cfg.clamp_min_nm,
        clamp_max_nm: cfg.clamp_max_nm,
        thin_layer_policy: cfg.thin_layer_policy,
        lm,
        clamp_accumulator: ClampReport::default(),
        envs: None,
    };
    let mut pipe = NeedlePipeline::new(stack, spectral, cfg, needle_cfg, pool)?;
    let report = pipe.run(&mut ctx, |cycle, phase, _det| callback(cycle, phase))?;
    Ok((report, pipe.stack, warnings))
}

/// F2.2: an end-to-end JOINT run over K environments.
///
/// `envs` is the compile of F2.1 — the shared design assembled once per
/// environment, plus the `design_slot → (env, span)` routing. The pipeline
/// carries environment 0's stack, which IS the shared design object; every
/// other environment is re-expressed from it at each eval
/// (`CompiledEnvironments::expand`), so a thickness step propagates to all
/// K by construction rather than by a synchronization step that could be
/// forgotten.
///
/// K = 1 runs the flat path: the context's `envs` is left `None`, so the
/// call sequence is the pre-F2.2 one op for op (§4.6's single branch, taken
/// here, once, outside every loop).
///
/// **Needles run jointly (F2.3).** The sweep scans every environment's own
/// assembly, sums P into shared buckets by design parameter, and inserts
/// once — into the shared design object and into every environment's
/// template at the same time, so the assemblies split together rather than
/// drifting apart.
///
/// **Removal and inflate are still refused while K > 1.** Both change the
/// design's span layout the other way — a slot disappears, or a new one
/// appears between two others — and the compile has no inverse of
/// `insert_seed` yet. Running them here would desynchronize the
/// environments or silently optimize environment 0 alone; refusing says so.
///
/// The same rule reaches the thin-layer policy, which is where it actually
/// bites: `remove` eliminates a sub-floor film on the pipeline's FINAL
/// sweep, and an elimination is a removal whoever ordered it. Joint runs
/// take `clamp_up_final` — the documented default for needle runs anyway —
/// and the context pins the floor as an LM bound for the duration
/// ([`SmatrixContext::optimize_thicknesses`]), so the span layout moves
/// only through `insert_seed`.
///
/// **Demands an environment cannot answer are dropped, once.** A phase
/// demand in an environment whose assembly holds an incoherent film (D4),
/// or a back demand where its exit medium absorbs, is dropped with one
/// warning per kind in the returned list; the same run keeps those
/// demands in the environments that can answer them. Only a run left with
/// no demand at all refuses ([`drop_infeasible`]).
#[allow(clippy::too_many_arguments)]
pub fn run_environments(
    envs: super::environments::CompiledEnvironments,
    pool: MaterialPool,
    wavelengths: &[f64],
    angles: &[Radians],
    spec: &MeritSpec,
    cfg: PipelineConfig,
    needle_cfg: NeedleCycleConfig,
    lm: LmConfig,
    mut callback: impl FnMut(usize, &super::pipeline::PipelinePhaseResult) -> Result<(), String>,
) -> Result<(PipelineResult, DesignStack, Vec<String>), String> {
    SmatrixContext::check_spec(spec)?;
    let k = envs.n_envs();
    // A segment's edge lists live in the compile, not on the flat stack,
    // so boundary needles that could reach a segment edge take the joint
    // path even in one environment.
    let joint = envs.is_joint() || (needle_cfg.boundary_needles && envs.has_segment_edges());
    if spec.n_envs() != k {
        return Err(format!(
            "run_environments: the merit spec has {} environment(s) but the \
             design compiles to {} ({}) - a demand would score against the \
             wrong surroundings, or against none",
            spec.n_envs(),
            k,
            envs.names().join(", ")
        ));
    }
    // M7 (review PB): the count matching proves nothing about WHICH
    // environment each demand scores against - a demand's `env_idx` is a
    // POSITION in this list, resolved when the spec was compiled and
    // unrecoverable afterwards. A spec that remembers its roster gets the
    // binding checked; one that does not keeps the count check above,
    // which is every pre-F2.4 caller.
    let named = spec.env_names();
    if !named.is_empty() && named != envs.names() {
        return Err(format!(
            "run_environments: the merit spec was compiled against \
             environments ({}) but the design compiles to ({}) - a demand's \
             environment is a POSITION in that list, so this scores every \
             demand against the wrong surroundings. Pass the same roster in \
             the same order, or hand the TargetCollection to the run door \
             and let it build the spec from the request",
            named.join(", "),
            envs.names().join(", ")
        ));
    }
    if joint {
        let mut blocked: Vec<&str> = Vec::new();
        if cfg.enable_cleanup {
            blocked.push("enable_cleanup");
        }
        if cfg.enable_inflate {
            blocked.push("enable_inflate");
        }
        if cfg.thin_layer_policy == crate::smatrix::synthesis::config::ThinLayerPolicy::Remove {
            blocked.push("thin_layer_policy 'remove'");
        }
        if !blocked.is_empty() {
            return Err(format!(
                "run_environments: {} environment(s){} with {} - removing a \
                 design parameter, or inserting one outside the needle path, \
                 changes the shared design's span layout in a direction the \
                 compile cannot yet follow (F2.3 gave it `insert_seed`, not \
                 its inverse). Joint thickness optimization and joint needles \
                 run today: use thin_layer_policy 'clamp_up_final' and leave \
                 cleanup and inflate off, or run one environment.",
                k,
                if envs.has_copies() {
                    ", a design segment used more than once,"
                } else {
                    ""
                },
                blocked.join(" and ")
            ));
        }
    }

    // Once, here: what each environment's complete stack cannot answer
    // is dropped with a warning; the same demand still binds the
    // environments that can.
    let (spec, warnings) = drop_infeasible(
        "run_environments",
        spec,
        envs.stacks(),
        envs.names(),
        wavelengths,
    )?;
    let spec = &spec;

    let stack = envs.stacks()[0].clone();
    let spectral = SpectralInputs::from_spec(spec, angles, wavelengths)?;
    let sin_theta = spectral.sin_theta.clone();
    let cfg = cfg.validated()?;
    let mut ctx = SmatrixContext {
        wavls: wavelengths.to_vec(),
        sin_theta,
        spec: spec.clone(),
        clamp_min_nm: cfg.clamp_min_nm,
        clamp_max_nm: cfg.clamp_max_nm,
        thin_layer_policy: cfg.thin_layer_policy,
        lm,
        clamp_accumulator: ClampReport::default(),
        // The branch. One environment whose design films are each used
        // once keeps the flat path exactly.
        envs: if joint { Some(Arc::new(envs)) } else { None },
    };
    let mut pipe = NeedlePipeline::new(stack, spectral, cfg, needle_cfg, pool)?;
    let report = pipe.run(&mut ctx, |cycle, phase, _det| callback(cycle, phase))?;
    // The final clamp moves primaries only; the returned stack's copies
    // follow them.
    if let Some(envs) = &ctx.envs {
        envs.sync(&mut pipe.stack)?;
    }
    Ok((report, pipe.stack, warnings))
}

/// `" in environment(s) 'a', 'b'"` for the environments `mask` marks, or
/// `""` when there are no names to give (the flat door's one stack).
fn env_list(names: &[String], mask: &[bool]) -> String {
    let marked: Vec<String> = names
        .iter()
        .zip(mask)
        .filter(|(_, m)| **m)
        .map(|(n, _)| format!("'{n}'"))
        .collect();
    if marked.is_empty() {
        String::new()
    } else {
        format!(" in environment(s) {}", marked.join(", "))
    }
}

/// The run doors' one-time feasibility pass: drop the demands the
/// assembled stacks cannot answer, and say so once.
///
/// Targets belong to the whole design; whether a demand can be computed
/// depends on the complete stack it scores against (`stacks[env]`), so
/// this is checked here, at run setup, and never inside the pipeline:
///
/// - **Phase across an incoherent film (D4).** Only intensities cross an
///   incoherent layer; a phase demand on such a stack is dropped.
/// - **Back incidence from an absorbing exit medium (W2).** A back demand
///   whose wavelengths reach a point where the exit medium has
///   k > [`BACK_EXIT_MAX_K`] is dropped: incidence from inside a
///   semi-infinite absorber has no defined incident flux.
///
/// Each kind gives ONE warning per run, however many demands it drops.
/// The run refuses only when nothing the caller asked for is left. The
/// evaluator keeps its per-solve refusals as a safety net for a
/// `SmatrixContext` built by hand; after this pass the run doors never
/// reach them, so the loops carry no check and emit nothing.
fn drop_infeasible(
    door: &str,
    spec: &MeritSpec,
    stacks: &[DesignStack],
    names: &[String],
    wavelengths: &[f64],
) -> Result<(MeritSpec, Vec<String>), String> {
    use super::evaluator::BACK_EXIT_MAX_K;
    let mut spec = spec.clone();
    let had_demands = spec.has_demands();
    let mut warnings = Vec::new();

    let lose_phase: Vec<bool> = stacks.iter().map(|st| st.has_incoherent()).collect();
    let dropped = spec.drop_phase_demands(&lose_phase);
    if !dropped.is_empty() {
        warnings.push(format!(
            "{} phase demand(s) dropped{}: the stack holds an incoherent film, \
             and phase is not defined across an incoherent layer (only \
             intensity, absorption and colour demands reach through one)",
            dropped.len(),
            env_list(names, &lose_phase)
        ));
    }

    // Per environment, the grid points where the exit medium absorbs.
    let absorbs: Vec<Vec<bool>> = stacks
        .iter()
        .map(|st| {
            st.substrate()
                .nk
                .iter()
                .map(|n| n.im.abs() > BACK_EXIT_MAX_K)
                .collect()
        })
        .collect();
    let mut hit = vec![false; stacks.len()];
    let dropped = spec.drop_targets(|t, curve| {
        let e = t.env_idx as usize;
        if !curve.is_back() || t.wavelengths.is_empty() {
            return false;
        }
        let Some(abs) = absorbs.get(e) else {
            return false;
        };
        let lo = t.wavelengths.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = t
            .wavelengths
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let bad = wavelengths
            .iter()
            .zip(abs)
            .any(|(&lam, &a)| a && lam >= lo && lam <= hi);
        if bad {
            hit[e] = true;
        }
        bad
    });
    if !dropped.is_empty() {
        let exits: Vec<String> = stacks
            .iter()
            .zip(&hit)
            .filter(|(_, h)| **h)
            .map(|(st, _)| format!("'{}'", st.substrate().material))
            .collect();
        warnings.push(format!(
            "{} back-incidence demand(s) dropped{}: the exit medium ({}) absorbs \
             (k > {BACK_EXIT_MAX_K:e}) over their wavelengths, and incidence from \
             inside a semi-infinite absorber is undefined (RB/TB/AB demands need \
             a transparent exit medium)",
            dropped.len(),
            env_list(names, &hit),
            exits.join(", ")
        ));
    }

    if had_demands && !spec.has_demands() {
        return Err(format!(
            "{door}: no demand is left to synthesize - {}",
            warnings.join("; ")
        ));
    }
    Ok((spec, warnings))
}

// ---------------------------------------------------------------------------
// Tests (standalone: no Python)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn films() -> Vec<ArrayFilm> {
        vec![
            ArrayFilm {
                name: "L".to_string(),
                nk: vec![Complex64::new(1.45, 0.0); 2],
                d_nm: 100.0,
                coherent: true,
                roughness: 0.0,
                rough_type: 0,
                inhomogen: false,
                inh_delta: 0.1,
                interface: false,
                interface_thickness: 0.0,
                optimize: true,
                needle: true,
                gradient: None,
            },
            ArrayFilm {
                name: "H".to_string(),
                nk: vec![Complex64::new(2.1, 0.0); 2],
                d_nm: 60.0,
                coherent: true,
                roughness: 0.0,
                rough_type: 0,
                inhomogen: false,
                inh_delta: 0.1,
                interface: false,
                interface_thickness: 0.0,
                optimize: true,
                needle: true,
                gradient: None,
            },
        ]
    }

    #[test]
    fn assemble_flat_stack() {
        let wl = vec![500.0, 600.0];
        let (stack, warns) = assemble_stack(
            "air",
            vec![Complex64::new(1.0, 0.0); 2],
            "sub",
            vec![Complex64::new(1.52, 0.0); 2],
            &films(),
            &HashMap::new(),
            &wl,
        )
        .unwrap();
        assert!(warns.is_empty());
        assert_eq!(stack.films().len(), 2);
    }

    #[test]
    fn assemble_refuses() {
        let wl = vec![500.0, 600.0];
        let air = vec![Complex64::new(1.0, 0.0); 2];
        let sub = vec![Complex64::new(1.52, 0.0); 2];
        let mut dup = films();
        dup[1].name = "L".to_string();
        assert!(
            assemble_stack(
                "air",
                air.clone(),
                "sub",
                sub.clone(),
                &dup,
                &HashMap::new(),
                &wl
            )
            .is_err()
        );
        let mut bad = films();
        bad[0].nk = vec![Complex64::new(1.45, 0.0); 3];
        assert!(assemble_stack("air", air, "sub", sub, &bad, &HashMap::new(), &wl).is_err());
    }

    /// One demand on `curve` at 0°, over both grid points.
    fn demand(spec: &mut MeritSpec, curve: super::super::merit::CurveId, phase: bool) {
        use super::super::merit::{ConstraintKind, MeritKey, MeritTarget, SimTransform};
        let key = spec.add_key(MeritKey {
            angle: Radians::ZERO,
            curve,
        }) as u32;
        spec.add_target(MeritTarget {
            key_idx: key,
            env_idx: 0,
            wavelengths: Arc::from([500.0, 600.0].as_slice()),
            kind: ConstraintKind::Exact,
            transform: if phase {
                SimTransform::Phase
            } else {
                SimTransform::Linear
            },
            norm_factor: 1.0,
            normalized_targets: Arc::from([0.0, 0.0].as_slice()),
            tolerances: Arc::from([0.1, 0.1].as_slice()),
            band: Arc::from([].as_slice()),
            phase,
            differential_passes: None,
            retardance: false,
            ratio: None,
            dispersion: None,
            point_weights: Vec::new().into(),
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
    }

    /// The setup-time pass drops what the complete stack cannot answer -
    /// phase through an incoherent film, back incidence from an absorbing
    /// exit - with one warning per kind, keeps the rest untouched, and
    /// refuses only when nothing is left.
    #[test]
    fn drop_infeasible_warns_once_per_kind() {
        use super::super::merit::CurveId;
        let wl = vec![500.0, 600.0];
        let mut fs = films();
        fs[1].coherent = false;
        fs[1].d_nm = 1.0e6;
        fs[1].optimize = false;
        fs[1].needle = false;
        let (stack, _) = assemble_stack(
            "air",
            vec![Complex64::new(1.0, 0.0); 2],
            "metal",
            vec![Complex64::new(0.2, 3.0); 2],
            &fs,
            &HashMap::new(),
            &wl,
        )
        .unwrap();
        let stacks = std::slice::from_ref(&stack);

        let mut spec = MeritSpec::new();
        demand(&mut spec, CurveId::Rs, false);
        demand(&mut spec, CurveId::Ts, true);
        demand(&mut spec, CurveId::Tp, true);
        demand(&mut spec, CurveId::TBs, false);
        let (kept, warnings) = drop_infeasible("t", &spec, stacks, &[], &wl).unwrap();
        assert_eq!(kept.targets().len(), 1);
        assert_eq!(warnings.len(), 2, "{warnings:?}");
        assert!(
            warnings[0].starts_with("2 phase demand(s) dropped:"),
            "{warnings:?}"
        );
        assert!(warnings[1].contains("'metal'"), "{warnings:?}");

        // A spec the stack can answer passes through unchanged and silent.
        let mut front = MeritSpec::new();
        demand(&mut front, CurveId::Rs, false);
        let (kept, warnings) = drop_infeasible("t", &front, stacks, &[], &wl).unwrap();
        assert!(warnings.is_empty());
        assert_eq!(kept.targets().len(), 1);

        // Nothing left: refuse, carrying the reasons.
        let mut back = MeritSpec::new();
        demand(&mut back, CurveId::RBp, false);
        let err = drop_infeasible("t", &back, stacks, &[], &wl).unwrap_err();
        assert!(
            err.contains("no demand is left") && err.contains("'metal'"),
            "{err}"
        );
    }
}
