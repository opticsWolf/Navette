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
use super::cycle::ContrastMap;
use super::cycle::NeedleCycleConfig;
use super::evaluator::SmatrixContext;
use super::merit::MeritSpec;
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

/// One contrast seed: host name → fresh seed carrier.
#[derive(Clone, Debug)]
pub struct ArraySeed {
    pub host: String,
    pub seed_name: String,
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
    .map(|(stack, warnings)| {
        film_warnings.extend(warnings);
        (stack, film_warnings)
    })
}

/// End-to-end design run: assemble, fold demands, execute the macro-loop.
/// `angles_deg` in degrees; `callback` fires per macro-cycle.
///
/// Returns `(report, stack, warnings)`. The warnings are the assembly's --
/// a dropped ambient absorption (R3.4), a homogenized graded film -- and they
/// must be surfaced upstream. They were dropped on the floor here until
/// 0.6.27, which made this the one design path that corrected inputs in
/// silence.
#[allow(clippy::too_many_arguments)]
pub fn run_design(
    ambient_name: &str,
    ambient_nk: Vec<Complex64>,
    substrate_name: &str,
    substrate_nk: Vec<Complex64>,
    films: &[ArrayFilm],
    groups: &HashMap<String, Group>,
    seeds: &[ArraySeed],
    wavelengths: &[f64],
    angles_deg: &[f64],
    spec: &MeritSpec,
    cfg: PipelineConfig,
    needle_cfg: NeedleCycleConfig,
    lm: LmConfig,
    mut callback: impl FnMut(usize, &super::pipeline::PipelinePhaseResult) -> Result<(), String>,
) -> Result<(PipelineResult, DesignStack, Vec<String>), String> {
    let (stack, warnings) = assemble_stack(
        ambient_name,
        ambient_nk,
        substrate_name,
        substrate_nk,
        films,
        groups,
        wavelengths,
    )?;
    let mut cmap = ContrastMap::new();
    for s in seeds {
        if s.nk.len() != wavelengths.len() {
            return Err(format!(
                "contrast '{}': nk length {} != {} wavelengths",
                s.host,
                s.nk.len(),
                wavelengths.len()
            ));
        }
        cmap.insert(
            Arc::from(s.host.as_str()),
            LayerSpec {
                material: Arc::from(s.seed_name.as_str()),
                nk: Arc::from(s.nk.clone()),
                d_nm: 0.0,
                coherent: true,
                rough_type: 0,
                rough_val: 0.0,
                optimize: true,
                needle: true,
            },
        );
    }
    let spectral = SpectralInputs::from_spec(spec, angles_deg, wavelengths)?;
    let cfg = cfg.validated()?;
    let sin_theta: Vec<f64> = angles_deg.iter().map(|a| a.to_radians().sin()).collect();
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
    let mut pipe = NeedlePipeline::new(stack, spectral, cfg, needle_cfg, cmap)?;
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
#[allow(clippy::too_many_arguments)]
pub fn run_environments(
    envs: super::environments::CompiledEnvironments,
    contrast: ContrastMap,
    wavelengths: &[f64],
    angles_deg: &[f64],
    spec: &MeritSpec,
    cfg: PipelineConfig,
    needle_cfg: NeedleCycleConfig,
    lm: LmConfig,
    mut callback: impl FnMut(usize, &super::pipeline::PipelinePhaseResult) -> Result<(), String>,
) -> Result<(PipelineResult, DesignStack, Vec<String>), String> {
    let k = envs.n_envs();
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
    if k > 1 {
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
                "run_environments: {} environments with {} - removing a \
                 design parameter, or inserting one outside the needle path, \
                 changes the shared design's span layout in a direction the \
                 compile cannot yet follow (F2.3 gave it `insert_seed`, not \
                 its inverse). Joint thickness optimization and joint needles \
                 run today: use thin_layer_policy 'clamp_up_final' and leave \
                 cleanup and inflate off, or run one environment.",
                k,
                blocked.join(" and ")
            ));
        }
    }

    let stack = envs.stacks()[0].clone();
    let spectral = SpectralInputs::from_spec(spec, angles_deg, wavelengths)?;
    let cfg = cfg.validated()?;
    let sin_theta: Vec<f64> = angles_deg.iter().map(|a| a.to_radians().sin()).collect();
    let mut ctx = SmatrixContext {
        wavls: wavelengths.to_vec(),
        sin_theta,
        spec: spec.clone(),
        clamp_min_nm: cfg.clamp_min_nm,
        clamp_max_nm: cfg.clamp_max_nm,
        thin_layer_policy: cfg.thin_layer_policy,
        lm,
        clamp_accumulator: ClampReport::default(),
        // The branch. One environment keeps the flat path exactly.
        envs: if k == 1 { None } else { Some(Arc::new(envs)) },
    };
    let mut pipe = NeedlePipeline::new(stack, spectral, cfg, needle_cfg, contrast)?;
    let report = pipe.run(&mut ctx, |cycle, phase, _det| callback(cycle, phase))?;
    Ok((report, pipe.stack, Vec::new()))
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
}
