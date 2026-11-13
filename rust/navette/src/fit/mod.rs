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
//! Phases 1-2 of `docs/plans/general_fit_plan.md`: one stage, FD Jacobian,
//! one sample. Stages, ties, priors and several samples are phase 3.

pub mod convert;
pub mod data;
pub mod dispersion;
pub mod forward;
pub mod model;
pub mod params;
pub mod sets;
mod stats;
pub mod strategy;

use std::collections::HashSet;

pub use convert::{
    AngleUnit, DeltaRange, DeltaSign, EllForm, InstrumentConvention, RawDataset, RawQuantity, Side,
    SigmaSource, ValueScale, WavelengthUnit, from_solver, to_solver,
};
pub use data::{NkColumn, NkData, NkQuantity, SpectrumData, channel_request};
pub use dispersion::{CandidateFit, fit_dispersion, table_to_model};
pub use forward::{DispersionForward, StackForward};
pub use model::{DEFAULT_ENV, FitModel, Segment};
pub use params::{ParamAddr, ParamKind, Parameter, Tie, Transform};
pub use sets::{ParamSet, Selector, select};
pub use strategy::{
    Accept, Repeat, Sampler, Stage, StageKind, StageLsq, StageReport, Strategy, StrategyReport,
    run_strategy,
};

use crate::lsq::backends::{LsqBackend, solve_least_squares};
use crate::lsq::lm::{LmConfig, LmTermination, NoJacobian, build_jacobian};
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
}

impl Compiled {
    fn new(problem: &FitProblem) -> Result<Self, String> {
        let mut base = problem.model.clone();
        let mut targets: Vec<Target> = Vec::new();
        let mut free = Vec::new();
        let mut resolved = Vec::new();
        let mut live = HashSet::new();
        for p in &problem.parameters {
            let addr = &p.addr;
            let target = base.resolve(addr)?;
            if targets.contains(&target) {
                return Err(format!("{addr}: the same scalar is named twice"));
            }
            if let Target::Material { material, .. } = &target {
                live.insert(material.clone());
            }
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
                let scale = p.scale.unwrap_or_else(|| default_scale(p.transform, x0));
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
                if let Target::Material { material, .. } = &target {
                    live.insert(material.clone());
                }
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
        })
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
        let dev = self.deviations(&self.model_at(u)?)?;
        for (d, rows) in dev.iter().zip(&self.rows) {
            out.extend(d.iter().zip(&rows.inv_sigma).map(|(d, w)| d * w));
        }
        self.prior_rows(u, out);
        if let Some(r) = out.iter().find(|r| !r.is_finite()) {
            return Err(format!("fit: non-finite residual {r}"));
        }
        Ok(())
    }
}

/// The [`DatasetFit`] id of the prior rows, when any parameter has one.
pub const PRIOR_ID: &str = "prior";

/// Every row of a trial point the forward model refuses: large enough
/// that no step lands there, small enough that χ² stays finite.
const REFUSED_RESIDUAL: f64 = 1e100;

/// Fit the free parameters of `problem` in one least-squares stage.
pub fn fit(problem: &FitProblem) -> Result<FitReport, String> {
    run(problem, true)
}

/// The report at the start values, without moving: χ², per-dataset
/// shares, and the covariance of the free parameters there. What a
/// strategy reports at its end, over everything any stage freed.
pub fn evaluate(problem: &FitProblem) -> Result<FitReport, String> {
    run(problem, false)
}

fn run(problem: &FitProblem, solve: bool) -> Result<FitReport, String> {
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
    let trial = |u: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
        if c.residuals(u, out).is_err() {
            out.clear();
            out.resize(m0, REFUSED_RESIDUAL);
        }
        Ok(())
    };
    let residuals = strict;

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
            None::<&NoJacobian>,
            &vec![0.0; n],
            &lb,
            &ub,
            &problem.lsq,
        )?;
        (r.x, r.termination, r.iterations, r.evals, r.backend)
    };

    let model = c.model_at(&u)?;
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

    // Covariance in u, then x = t⁻¹(t0 + scale·u): dx/du = scale·t⁻¹'.
    let (covariance, covariance_note) = if n == 0 {
        (None, None)
    } else {
        let mut jac = Vec::new();
        match build_jacobian(&residuals, &u, &mut jac)
            .and_then(|_| stats::inverse_normal_matrix(&jac, m, n))
        {
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
