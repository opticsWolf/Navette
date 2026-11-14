// SPDX-License-Identifier: LGPL-3.0-or-later
//! The general fit engine: any scalar of a model against measured spectra.
//!
//! Needle synthesis refines film thicknesses against design targets
//! (`smatrix::synthesis::thickness_refine`). This module fits *any* scalar
//! of a model (layer thickness, grading, roughness, interface width,
//! gradient profile, group factors and summands, material dispersion
//! parameters, EMA fractions) against measured curves, through the general
//! solver (`expand` → `solve_arrays`), so every coherence mode, inversion,
//! group and grading rule the solver knows applies unchanged.
//!
//! - [`params`]: addresses ([`ParamAddr`]), [`Parameter`] with bounds,
//!   scale and transform.
//! - [`model`]: [`FitModel`], the `Send + Sync` mirror of an architect.
//! - [`data`]: [`SpectrumData`], [`NkData`] and the residual row convention.
//! - [`forward`]: [`StackForward`], one evaluation through the solver;
//!   [`DispersionForward`], material specs on measured nk grids. Both feed
//!   one residual vector, so a joint fit is one problem.
//! - [`fit()`]: one least-squares stage over the free parameters, and the
//!   [`FitReport`] (χ², AIC/BIC, covariance, correlations, parameters at
//!   bounds).
//! - [`dispersion`]: the doors [`fit_dispersion`] (one spec against nk
//!   data) and [`table_to_model`] (rank model candidates for a table).
//!
//! - [`strategy`], [`sets`], [`explore`]: staged strategies, parameter
//!   sets, generated strategies, the tournament and the benchmark.
//! - [`request`]: the fit as data, a program document's `fit` section,
//!   and the doors the Python package calls.
//! - [`Watcher`]: progress per stage and per lower χ², and a way to stop.

pub mod convert;
pub mod data;
pub mod dispersion;
pub mod doors;
pub mod explore;
pub mod forward;
pub mod global;
mod jacobian;
pub mod model;
pub mod params;
pub mod request;
pub mod sets;
mod stats;
pub mod strategy;

use std::collections::HashSet;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

pub use convert::{
    AngleUnit, DeltaRange, DeltaSign, EllForm, InstrumentConvention, RawDataset, RawQuantity, Side,
    SigmaSource, ValueScale, WavelengthUnit, from_solver, to_solver,
};
pub use data::{Instrument, NkColumn, NkData, NkQuantity, SpectrumData, channel_request};
pub use dispersion::{CandidateFit, fit_dispersion, table_to_model};
pub use explore::{
    BenchmarkEntry, BenchmarkReport, DEFAULT_ORDERING_CAP, Generated, RunScore, StrategyGen,
    TournamentEntry, Truth, benchmark, draw_truth, generate, tournament,
};
pub use forward::{DispersionForward, StackForward};
pub use global::GlobalMethod;
pub use model::{DEFAULT_ENV, FREEFORM_VALUES, FitModel, Segment};
pub use params::{InstrumentField, ParamAddr, ParamKind, Parameter, Tie, Transform};
pub use request::{FitRequest, FitSection, fit_program, request_from_json};
pub use sets::{ParamSet, Selector, select};
pub use strategy::{
    Accept, Repeat, Sampler, Stage, StageKind, StageLsq, StageReport, Strategy, StrategyReport,
    run_strategy, run_strategy_watched,
};

use crate::lsq::backends::{LsqBackend, solve_least_squares};
use crate::lsq::lm::{JacobianMode, JacobianSource, LmConfig, LmTermination, build_jacobian};
use model::Target;
use params::{default_bounds, default_scale, material_bounds};

/// What to fit: a model, the data it must reproduce, and which of its
/// scalars may move.
#[derive(Debug, Clone)]
pub struct FitProblem {
    pub model: FitModel,
    pub data: Vec<SpectrumData>,
    /// Measured optical constants of the model's materials; their rows
    /// follow the spectra's.
    pub nk_data: Vec<NkData>,
    /// Free and fixed parameters. A fixed one is written once and held.
    pub parameters: Vec<Parameter>,
    /// Scalars that follow a driver; never parameters themselves.
    pub ties: Vec<Tie>,
    /// Smoothness rows over free-form materials' values; they follow the
    /// priors.
    pub smoothness: Vec<Smoothness>,
    /// `solve_arrays`' mode: 0 front block, 1 coherency matrix, 2 fully
    /// coherent.
    pub coherence_mode: i32,
    /// Solver settings; `backend` picks the solver.
    pub lsq: LmConfig,
}

impl FitProblem {
    pub fn new(model: FitModel, data: Vec<SpectrumData>, parameters: Vec<Parameter>) -> Self {
        Self {
            model,
            data,
            nk_data: Vec::new(),
            parameters,
            ties: Vec::new(),
            smoothness: Vec::new(),
            coherence_mode: 0,
            lsq: LmConfig::default(),
        }
    }

    /// A problem against nk data only (a dispersion fit).
    pub fn dispersion(model: FitModel, nk_data: Vec<NkData>, parameters: Vec<Parameter>) -> Self {
        Self {
            nk_data,
            ..Self::new(model, Vec::new(), parameters)
        }
    }

    /// The model's curves at the start values, one vector per dataset at
    /// its kept points (the dataset values are not read): the spectra, then
    /// the nk datasets with their columns in order.
    pub fn simulate(&self) -> Result<Vec<Vec<f64>>, String> {
        let c = Compiled::new(self)?;
        c.simulate(&c.base)
    }
}

/// A smoothness prior on a free-form material (§3.5.4): one residual row
/// `(v[i-1] - 2 v[i] + v[i+1]) / sigma` per inner knot of the list
/// `material[material].param` (a [`model::FREEFORM_VALUES`] list), fixed
/// and free knots alike.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Smoothness {
    pub material: String,
    pub param: String,
    pub sigma: f64,
}

impl Smoothness {
    pub fn new(material: impl Into<String>, param: impl Into<String>, sigma: f64) -> Self {
        Self {
            material: material.into(),
            param: param.into(),
            sigma,
        }
    }
}

/// The values of `material[material].param` in `model`.
fn knot_values(model: &FitModel, s: &Smoothness) -> Result<Vec<f64>, String> {
    let what = format!("smoothness material[{}].{}", s.material, s.param);
    let Some(crate::structure::Entry::Spec(spec)) = model.materials().get(&s.material) else {
        return Err(format!("{what}: no material spec {:?}", s.material));
    };
    if !model::FREEFORM_VALUES.contains(&(spec.model.as_str(), s.param.as_str())) {
        return Err(format!(
            "{what}: {} {:?} is not a free-form list; those are {:?}",
            spec.model,
            s.param,
            model::FREEFORM_VALUES
        ));
    }
    spec.params
        .get(&s.param)
        .and_then(|v| v.as_array())
        .and_then(|a| a.iter().map(|x| x.as_f64()).collect::<Option<Vec<f64>>>())
        .ok_or_else(|| format!("{what}: not a list of numbers"))
}

/// One parameter after the fit.
#[derive(Debug, Clone, PartialEq)]
pub struct FittedParameter {
    pub addr: String,
    pub start: f64,
    pub value: f64,
    /// `sqrt(C_ii)`; `None` for a fixed parameter or without a covariance.
    pub std_error: Option<f64>,
    pub lb: f64,
    pub ub: f64,
    pub fixed: bool,
    /// The solution sits on `lb` or `ub`: the covariance there describes a
    /// constrained point and understates the uncertainty.
    pub at_bound: bool,
}

/// One dataset at the solution.
#[derive(Debug, Clone, PartialEq)]
pub struct DatasetFit {
    pub id: String,
    /// How the data was read from the instrument (`SpectrumData::source`).
    pub source: Option<String>,
    pub n_points: usize,
    /// This dataset's share of χ².
    pub chi2: f64,
    /// RMS of `sim - value` (wrapped for angular channels), in the
    /// channel's own unit.
    pub rms: f64,
}

/// Result of [`fit`].
#[derive(Debug, Clone)]
pub struct FitReport {
    /// In the problem's order, fixed ones included.
    pub parameters: Vec<FittedParameter>,
    pub datasets: Vec<DatasetFit>,
    /// Σ r², r = sqrt(weight)·(sim − value)/σ.
    pub chi2: f64,
    /// χ² / (m − n_free); NaN when m <= n_free.
    pub chi2_reduced: f64,
    /// `m·ln(χ²/m) + 2·n_free`: the Gaussian AIC with the noise scale
    /// estimated by the fit, so a ranking of models on the same data does
    /// not depend on how σ was set. Only differences mean anything.
    pub aic: f64,
    /// `m·ln(χ²/m) + n_free·ln m`, likewise.
    pub bic: f64,
    pub n_residuals: usize,
    pub n_free: usize,
    /// `n_free × n_free`, row-major, free parameters in problem order, in
    /// the parameters' own units: `(JᵀJ)⁻¹ · χ²_red`.
    pub covariance: Option<Vec<f64>>,
    pub correlation: Option<Vec<f64>>,
    /// Why there is no covariance, when there is none.
    pub covariance_note: Option<String>,
    pub termination: LmTermination,
    pub iterations: usize,
    pub evals: usize,
    pub backend: LsqBackend,
    /// The model with the fitted values written in.
    pub model: FitModel,
}

impl FitReport {
    pub fn parameter(&self, addr: &str) -> Option<&FittedParameter> {
        self.parameters.iter().find(|p| p.addr == addr)
    }
}

/// One free parameter, compiled.
#[derive(Debug, Clone)]
struct Free {
    target: Target,
    transform: Transform,
    /// `t(x0)`.
    t0: f64,
    scale: f64,
    /// `(p0, sigma)`.
    prior: Option<(f64, f64)>,
}

/// A tie, resolved: `(target, factor, offset)` per follower.
#[derive(Debug, Clone)]
struct CompiledTie {
    driver: Target,
    targets: Vec<(Target, f64, f64)>,
}

fn apply_ties(model: &mut FitModel, ties: &[CompiledTie]) -> Result<(), String> {
    for t in ties {
        let x = model.read(&t.driver)?;
        for (target, factor, offset) in &t.targets {
            model.write(target, factor * x + offset)?;
        }
    }
    Ok(())
}

impl Free {
    fn x(&self, u: f64) -> f64 {
        self.transform.inverse(self.t0 + self.scale * u)
    }

    fn u(&self, x: f64) -> f64 {
        (self.transform.forward(x) - self.t0) / self.scale
    }
}

/// One dataset's targets at its kept points.
#[derive(Debug, Clone)]
struct Rows {
    values: Vec<f64>,
    /// sqrt(weight) / σ per point.
    inv_sigma: Vec<f64>,
    angular: bool,
}

/// A problem ready to evaluate.
struct Compiled {
    /// The model with every start value written in.
    base: FitModel,
    stack: Option<StackForward>,
    dispersion: Option<DispersionForward>,
    free: Vec<Free>,
    ties: Vec<CompiledTie>,
    /// Per parameter (problem order): start, lb, ub, index into `free`.
    resolved: Vec<(f64, f64, f64, Option<usize>)>,
    /// Spectra, then nk datasets.
    /// Dataset id and source, data then nk data.
    ids: Vec<(String, Option<String>)>,
    rows: Vec<Rows>,
    /// Smoothness priors and their row counts.
    smooth: Vec<(Smoothness, usize)>,
}

impl Compiled {
    fn new(problem: &FitProblem) -> Result<Self, String> {
        let mut ids = HashSet::new();
        if let Some(d) = problem.data.iter().find(|d| !ids.insert(d.id.as_str())) {
            return Err(format!("fit: two datasets are named {:?}", d.id));
        }
        let mut base = problem.model.clone();
        base.seed_instruments(&problem.data);
        let mut targets: Vec<Target> = Vec::new();
        let mut free = Vec::new();
        let mut resolved = Vec::new();
        let mut live = HashSet::new();
        let mut moving = HashSet::new();
        let mut note = |t: &Target| match t {
            Target::Material { material, .. } => {
                live.insert(material.clone());
            }
            Target::Data { dataset, field } if field.moves_geometry() => {
                moving.insert(dataset.clone());
            }
            _ => {}
        };
        for p in &problem.parameters {
            let addr = &p.addr;
            let target = base.resolve(addr)?;
            if targets.contains(&target) {
                return Err(format!("{addr}: the same scalar is named twice"));
            }
            note(&target);
            let x0 = match p.value {
                Some(v) => v,
                None => base.read(&target)?,
            };
            base.write(&target, x0)
                .map_err(|e| format!("{addr}: {e}"))?;
            let (dlb, dub) = match base.material_slot(&target) {
                Some((model, leaf)) => material_bounds(&model, &leaf),
                None => default_bounds(addr, x0),
            };
            let mut lb = p.lb.unwrap_or(dlb);
            let ub = p.ub.unwrap_or(dub);
            if p.transform == Transform::Log {
                if p.lb.is_none() && lb <= 0.0 {
                    lb = f64::MIN_POSITIVE;
                }
                if !(lb > 0.0 && x0 > 0.0) {
                    return Err(format!(
                        "{addr}: a log parameter needs a positive start and lower bound (start {x0}, lb {lb})"
                    ));
                }
            }
            if !(lb <= x0 && x0 <= ub) {
                return Err(format!("{addr}: start {x0} is outside [{lb}, {ub}]"));
            }
            let index = if p.fixed {
                None
            } else {
                if !(lb < ub) {
                    return Err(format!("{addr}: bounds [{lb}, {ub}] leave nothing to fit"));
                }
                let scale = p.scale.unwrap_or_else(|| match addr {
                    ParamAddr::Data { field, .. } if p.transform == Transform::Linear => {
                        field.scale()
                    }
                    _ => default_scale(p.transform, x0),
                });
                if !(scale > 0.0 && scale.is_finite()) {
                    return Err(format!(
                        "{addr}: scale {scale} is not a positive finite number"
                    ));
                }
                if let Some((_, s)) = p.prior
                    && !(s > 0.0 && s.is_finite())
                {
                    return Err(format!(
                        "{addr}: prior sigma {s} is not a positive finite number"
                    ));
                }
                free.push(Free {
                    target: target.clone(),
                    transform: p.transform,
                    t0: p.transform.forward(x0),
                    scale,
                    prior: p.prior,
                });
                Some(free.len() - 1)
            };
            targets.push(target);
            resolved.push((x0, lb, ub, index));
        }
        let mut followers: Vec<Target> = Vec::new();
        let mut ties = Vec::new();
        for tie in &problem.ties {
            let driver = base.resolve(&tie.driver).map_err(|e| format!("tie: {e}"))?;
            let mut compiled = Vec::new();
            for (addr, factor, offset) in &tie.targets {
                let target = base.resolve(addr).map_err(|e| format!("tie: {e}"))?;
                if targets.contains(&target) {
                    return Err(format!(
                        "tie: {addr} is a parameter; a tied scalar follows its driver"
                    ));
                }
                if followers.contains(&target) || target == driver {
                    return Err(format!("tie: {addr} follows more than one driver"));
                }
                if !(factor.is_finite() && offset.is_finite()) {
                    return Err(format!("tie: {addr}: factor and offset must be finite"));
                }
                note(&target);
                followers.push(target.clone());
                compiled.push((target, *factor, *offset));
            }
            ties.push(CompiledTie {
                driver,
                targets: compiled,
            });
        }
        if ties.iter().any(|t| followers.contains(&t.driver)) {
            return Err(
                "tie: a driver is itself tied; tie the followers to the first driver".into(),
            );
        }
        apply_ties(&mut base, &ties)?;
        base.check()?;
        let mut smooth = Vec::new();
        for s in &problem.smoothness {
            if !(s.sigma.is_finite() && s.sigma > 0.0) {
                return Err(format!(
                    "smoothness material[{}].{}: sigma {} is not a positive finite number",
                    s.material, s.param, s.sigma
                ));
            }
            let n = knot_values(&base, s)?.len();
            if n < 3 {
                return Err(format!(
                    "smoothness material[{}].{}: needs at least 3 knots, has {n}",
                    s.material, s.param
                ));
            }
            smooth.push((s.clone(), n - 2));
        }
        if problem.data.is_empty() && problem.nk_data.is_empty() {
            return Err("fit: no data".into());
        }
        let stack = if problem.data.is_empty() {
            None
        } else {
            Some(StackForward::new(
                &base,
                &problem.data,
                problem.coherence_mode,
                &live,
                &moving,
            )?)
        };
        let dispersion = if problem.nk_data.is_empty() {
            None
        } else {
            Some(DispersionForward::new(&base, &problem.nk_data)?)
        };
        let mut rows: Vec<Rows> = problem
            .data
            .iter()
            .map(|d| {
                let kept = d.kept();
                let sw = d.weight.sqrt();
                Rows {
                    values: kept.iter().map(|&i| d.values[i]).collect(),
                    inv_sigma: kept.iter().map(|&i| sw / d.sigma_at(i)).collect(),
                    angular: data::is_angular(&d.channel),
                }
            })
            .collect();
        rows.extend(problem.nk_data.iter().map(|d| {
            let kept = d.kept();
            let sw = d.weight.sqrt();
            Rows {
                values: d
                    .columns
                    .iter()
                    .flat_map(|c| kept.iter().map(|&i| c.values[i]))
                    .collect(),
                inv_sigma: d
                    .columns
                    .iter()
                    .flat_map(|c| kept.iter().map(move |&i| sw / c.sigma_at(i)))
                    .collect(),
                angular: false,
            }
        }));
        let ids = problem
            .data
            .iter()
            .map(|d| (d.id.clone(), d.source.clone()))
            .chain(problem.nk_data.iter().map(|d| (d.id.clone(), None)))
            .collect();
        Ok(Self {
            base,
            stack,
            dispersion,
            free,
            ties,
            resolved,
            ids,
            rows,
            smooth,
        })
    }

    /// χ² at free values `x` (free order), `None` where the model cannot
    /// be evaluated.
    #[cfg_attr(not(feature = "opt-basin"), allow(dead_code))]
    fn chi2_at(&self, x: &[f64]) -> Option<f64> {
        let u: Vec<f64> = self.free.iter().zip(x).map(|(f, &x)| f.u(x)).collect();
        let mut r = Vec::new();
        self.residuals(&u, &mut r).ok()?;
        Some(r.iter().map(|v| v * v).sum())
    }

    /// Prior rows `(x - p0) / sigma` of the free parameters that have one.
    fn prior_rows(&self, u: &[f64], out: &mut Vec<f64>) {
        for (f, &uj) in self.free.iter().zip(u) {
            if let Some((p0, s)) = f.prior {
                out.push((f.x(uj) - p0) / s);
            }
        }
    }

    fn n_priors(&self) -> usize {
        self.free.iter().filter(|f| f.prior.is_some()).count()
    }

    /// Smoothness rows of `model`, prior by prior.
    fn smooth_rows(&self, model: &FitModel, out: &mut Vec<f64>) -> Result<(), String> {
        for (s, _) in &self.smooth {
            let v = knot_values(model, s)?;
            out.extend(v.windows(3).map(|w| (w[0] - 2.0 * w[1] + w[2]) / s.sigma));
        }
        Ok(())
    }

    fn n_smooth(&self) -> usize {
        self.smooth.iter().map(|(_, n)| n).sum()
    }

    /// `target` is a knot of a smoothed list.
    fn smooths(&self, target: &Target) -> bool {
        matches!(target, Target::Material { material, key, .. }
            if self.smooth.iter().any(|(s, _)| s.material == *material && s.param == *key))
    }

    /// Spectra, then nk datasets.
    fn simulate(&self, model: &FitModel) -> Result<Vec<Vec<f64>>, String> {
        let mut sims = match &self.stack {
            Some(s) => s.simulate(model)?,
            None => Vec::new(),
        };
        if let Some(d) = &self.dispersion {
            sims.extend(d.simulate(model)?);
        }
        Ok(sims)
    }

    fn model_at(&self, u: &[f64]) -> Result<FitModel, String> {
        let mut model = self.base.clone();
        for (f, &uj) in self.free.iter().zip(u) {
            model.write(&f.target, f.x(uj))?;
        }
        apply_ties(&mut model, &self.ties)?;
        Ok(model)
    }

    /// Raw `sim - value` per dataset (wrapped where angular).
    fn deviations(&self, model: &FitModel) -> Result<Vec<Vec<f64>>, String> {
        let sims = self.simulate(model)?;
        Ok(sims
            .iter()
            .zip(&self.rows)
            .map(|(sim, rows)| {
                sim.iter()
                    .zip(&rows.values)
                    .map(|(s, v)| {
                        let d = s - v;
                        if rows.angular { data::wrap_angle(d) } else { d }
                    })
                    .collect()
            })
            .collect())
    }

    fn residuals(&self, u: &[f64], out: &mut Vec<f64>) -> Result<(), String> {
        out.clear();
        let model = self.model_at(u)?;
        let dev = self.deviations(&model)?;
        for (d, rows) in dev.iter().zip(&self.rows) {
            out.extend(d.iter().zip(&rows.inv_sigma).map(|(d, w)| d * w));
        }
        self.prior_rows(u, out);
        self.smooth_rows(&model, out)?;
        if let Some(r) = out.iter().find(|r| !r.is_finite()) {
            return Err(format!("fit: non-finite residual {r}"));
        }
        Ok(())
    }
}

/// The [`DatasetFit`] id of the prior rows, when any parameter has one.
pub const PRIOR_ID: &str = "prior";

/// The [`DatasetFit`] id of the smoothness rows, when there are any.
pub const SMOOTH_ID: &str = "smoothness";

/// Every row of a trial point the forward model refuses: large enough
/// that no step lands there, small enough that χ² stays finite.
const REFUSED_RESIDUAL: f64 = 1e100;

/// What a running fit tells its [`Watcher`].
#[derive(Debug, Clone, PartialEq)]
pub enum Progress {
    /// A stage starts; `stage` counts from 0.
    StageStart {
        stage: usize,
        name: String,
        free: Vec<String>,
    },
    /// A least-squares run of the stage reached a lower χ² than any of
    /// its points before: about once per accepted solver step. `evals`
    /// counts this run's evaluations so far.
    Improved {
        stage: usize,
        evals: usize,
        chi2: f64,
    },
    /// A `global` stage's search found a lower χ² than any point before
    /// it (the start included). `evals` counts the search's evaluations.
    Searched {
        stage: usize,
        evals: usize,
        chi2: f64,
    },
    /// A stage ended; `chi2` is its kept result on its own data.
    StageEnd {
        stage: usize,
        name: String,
        chi2: f64,
        accepted: bool,
    },
}

/// A progress callback and the stop it can request.
///
/// The callback may run on several threads at once (a multi-start stage
/// fits its starts in parallel). Returning `Err` stops the whole fit: every
/// solver run still going is starved (its trial points are refused without
/// being evaluated) and the fit returns that error.
pub struct Watcher<'a> {
    notify: Box<dyn Fn(&Progress) -> Result<(), String> + Sync + 'a>,
    stopped: Mutex<Option<String>>,
}

impl<'a> Watcher<'a> {
    pub fn new(notify: impl Fn(&Progress) -> Result<(), String> + Sync + 'a) -> Self {
        Self {
            notify: Box::new(notify),
            stopped: Mutex::new(None),
        }
    }

    /// Pass `p` on, unless the fit is already stopping.
    pub fn tell(&self, p: &Progress) {
        if self.is_stopped() {
            return;
        }
        if let Err(e) = (self.notify)(p) {
            let mut s = self.stopped.lock().unwrap_or_else(|e| e.into_inner());
            s.get_or_insert(e);
        }
    }

    pub fn is_stopped(&self) -> bool {
        self.stopped
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .is_some()
    }

    /// `Err` with the callback's message when it asked to stop.
    pub fn check(&self) -> Result<(), String> {
        match &*self.stopped.lock().unwrap_or_else(|e| e.into_inner()) {
            Some(e) => Err(format!("fit stopped by the progress callback: {e}")),
            None => Ok(()),
        }
    }
}

impl std::fmt::Debug for Watcher<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Watcher")
            .field("stopped", &self.stopped)
            .finish_non_exhaustive()
    }
}

/// Fit the free parameters of `problem` in one least-squares stage.
pub fn fit(problem: &FitProblem) -> Result<FitReport, String> {
    run(problem, true, None)
}

/// [`fit`], telling `watch` (as stage `stage`) each lower χ².
pub fn fit_watched(
    problem: &FitProblem,
    watch: Option<(&Watcher<'_>, usize)>,
) -> Result<FitReport, String> {
    run(problem, true, watch)
}

/// The report at the start values, without moving: χ², per-dataset
/// shares, and the covariance of the free parameters there. What a
/// strategy reports at its end, over everything any stage freed.
pub fn evaluate(problem: &FitProblem) -> Result<FitReport, String> {
    run(problem, false, None)
}

/// χ² of a problem as a function of its free parameters' values (in
/// problem order), compiled once: what a global search samples.
#[cfg_attr(not(feature = "opt-basin"), allow(dead_code))]
pub(crate) struct Chi2Fn(Compiled);

#[cfg_attr(not(feature = "opt-basin"), allow(dead_code))]
impl Chi2Fn {
    pub(crate) fn new(problem: &FitProblem) -> Result<Self, String> {
        Compiled::new(problem).map(Self)
    }

    /// `None` where the model cannot be evaluated.
    pub(crate) fn chi2(&self, x: &[f64]) -> Option<f64> {
        self.0.chi2_at(x)
    }
}

fn run(
    problem: &FitProblem,
    solve: bool,
    watch: Option<(&Watcher<'_>, usize)>,
) -> Result<FitReport, String> {
    let c = Compiled::new(problem)?;
    let n = c.free.len();
    // The start is evaluated strictly: a model that cannot be evaluated
    // there is the caller's error. A trial point the forward model refuses
    // later (a Sellmeier pole moved into the range, say) costs a large
    // residual instead, so the solver rejects that step rather than the
    // whole fit. The solution and its Jacobian are strict again.
    let strict = |u: &[f64], out: &mut Vec<f64>| c.residuals(u, out);
    let mut r0 = Vec::new();
    strict(&vec![0.0; n], &mut r0)?;
    let m0 = r0.len();
    let evals_seen = AtomicUsize::new(0);
    let best_seen = Mutex::new(f64::INFINITY);
    let trial = |u: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
        let refuse = |out: &mut Vec<f64>| {
            out.clear();
            out.resize(m0, REFUSED_RESIDUAL);
        };
        if watch.is_some_and(|(w, _)| w.is_stopped()) {
            refuse(out);
            return Ok(());
        }
        if c.residuals(u, out).is_err() {
            refuse(out);
            return Ok(());
        }
        if let Some((w, stage)) = watch {
            let evals = evals_seen.fetch_add(1, Ordering::Relaxed) + 1;
            let chi2: f64 = out.iter().map(|x| x * x).sum();
            let lower = {
                let mut best = best_seen.lock().unwrap_or_else(|e| e.into_inner());
                let lower = chi2 < *best;
                if lower {
                    *best = chi2;
                }
                lower
            };
            if lower {
                w.tell(&Progress::Improved { stage, evals, chi2 });
            }
        }
        Ok(())
    };
    let residuals = strict;
    let source = match problem.lsq.jacobian {
        JacobianMode::Analytic if n > 0 => {
            Some(jacobian::FitJacobian::new(&c, watch.map(|(w, _)| w))?)
        }
        _ => None,
    };

    let (u, termination, iterations, evals, backend) = if n == 0 || !solve {
        (
            vec![0.0; n],
            LmTermination::NoFreeParameters,
            0,
            0,
            problem.lsq.backend,
        )
    } else {
        let mut lb = Vec::with_capacity(n);
        let mut ub = Vec::with_capacity(n);
        for (_, plb, pub_, index) in &c.resolved {
            if let Some(j) = index {
                lb.push(c.free[*j].u(*plb));
                ub.push(c.free[*j].u(*pub_));
            }
        }
        let r = solve_least_squares(
            &trial,
            source.as_ref(),
            &vec![0.0; n],
            &lb,
            &ub,
            &problem.lsq,
        )?;
        if let Some((w, _)) = watch {
            w.check()?;
        }
        (r.x, r.termination, r.iterations, r.evals, r.backend)
    };

    // Instruments are seeded per dataset so `data[..]` resolves; the
    // report keeps only those that depart from the identity.
    let mut model = c.model_at(&u)?;
    model.prune_instruments(&problem.data);
    let mut r = Vec::new();
    residuals(&u, &mut r)?;
    let m = r.len();
    let chi2: f64 = r.iter().map(|x| x * x).sum();
    let chi2_reduced = if m > n {
        chi2 / (m - n) as f64
    } else {
        f64::NAN
    };
    let (mf, nf) = (m as f64, n as f64);
    let aic = mf * (chi2 / mf).ln() + 2.0 * nf;
    let bic = mf * (chi2 / mf).ln() + nf * mf.ln();

    let dev = c.deviations(&model)?;
    let mut offset = 0;
    let mut datasets: Vec<DatasetFit> = c
        .ids
        .iter()
        .zip(&dev)
        .map(|((id, source), dv)| {
            let k = dv.len();
            let share = r[offset..offset + k].iter().map(|x| x * x).sum();
            offset += k;
            DatasetFit {
                id: id.clone(),
                source: source.clone(),
                n_points: k,
                chi2: share,
                rms: (dv.iter().map(|x| x * x).sum::<f64>() / k as f64).sqrt(),
            }
        })
        .collect();
    let k = c.n_priors();
    if k > 0 {
        // RMS in units of each prior's sigma: the rows themselves.
        let rows = &r[offset..offset + k];
        datasets.push(DatasetFit {
            id: PRIOR_ID.to_string(),
            source: None,
            n_points: k,
            chi2: rows.iter().map(|x| x * x).sum(),
            rms: (rows.iter().map(|x| x * x).sum::<f64>() / k as f64).sqrt(),
        });
    }
    let ks = c.n_smooth();
    if ks > 0 {
        let rows = &r[offset + k..offset + k + ks];
        datasets.push(DatasetFit {
            id: SMOOTH_ID.to_string(),
            source: None,
            n_points: ks,
            chi2: rows.iter().map(|x| x * x).sum(),
            rms: (rows.iter().map(|x| x * x).sum::<f64>() / ks as f64).sqrt(),
        });
    }

    // Covariance in u, then x = t⁻¹(t0 + scale·u): dx/du = scale·t⁻¹'.
    let (covariance, covariance_note) = if n == 0 {
        (None, None)
    } else {
        let mut jac = Vec::new();
        // The watcher may have stopped the fit: the source then declines,
        // and the covariance falls back to central differences.
        let filled = match &source {
            Some(s) => s.fill(&u, &mut jac),
            None => Ok(None),
        };
        let filled = match filled {
            Ok(Some(_)) => Ok(()),
            Ok(None) => build_jacobian(&residuals, &u, &mut jac).map(|_| ()),
            Err(e) => Err(e),
        };
        match filled.and_then(|_| stats::inverse_normal_matrix(&jac, m, n)) {
            Ok(cu) if chi2_reduced.is_finite() => {
                let dxdu: Vec<f64> = c
                    .free
                    .iter()
                    .zip(&u)
                    .map(|(f, &uj)| f.scale * f.transform.slope(f.x(uj)))
                    .collect();
                let mut cx = cu;
                for i in 0..n {
                    for k in 0..n {
                        cx[i * n + k] *= dxdu[i] * dxdu[k] * chi2_reduced;
                    }
                }
                (Some(cx), None)
            }
            Ok(_) => (None, Some(format!("{m} residuals for {n} free parameters"))),
            Err(e) => (None, Some(e)),
        }
    };
    let correlation = covariance.as_ref().map(|cv| stats::correlation(cv, n));

    let parameters = problem
        .parameters
        .iter()
        .zip(&c.resolved)
        .map(|(p, &(x0, lb, ub, index))| {
            let (value, std_error, at_bound) = match index {
                Some(j) => {
                    let f = &c.free[j];
                    let (ulb, uub) = (f.u(lb), f.u(ub));
                    let near =
                        |b: f64| b.is_finite() && (u[j] - b).abs() <= 1e-10 * b.abs().max(1.0);
                    (
                        f.x(u[j]),
                        covariance.as_ref().map(|cv| cv[j * n + j].sqrt()),
                        near(ulb) || near(uub),
                    )
                }
                None => (x0, None, false),
            };
            FittedParameter {
                addr: p.addr.to_string(),
                start: x0,
                value,
                std_error,
                lb,
                ub,
                fixed: index.is_none(),
                at_bound,
            }
        })
        .collect();

    Ok(FitReport {
        parameters,
        datasets,
        chi2,
        chi2_reduced,
        aic,
        bic,
        n_residuals: m,
        n_free: n,
        covariance,
        correlation,
        covariance_note,
        termination,
        iterations,
        evals,
        backend,
        model,
    })
}
