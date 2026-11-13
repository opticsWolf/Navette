// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::forward::stack — a fit model through the general solver.
//!
//! One evaluation: the model (parameters already written) → materials on
//! the solve grid → [`expand`] → [`solve_arrays`] → the channel each
//! dataset reads. Groups, grading, gradients, interfaces, roughness,
//! inversion and every coherence mode come with `expand` and the solver;
//! nothing here re-implements them.
//!
//! **One solve per environment per evaluation.** An environment's grid is
//! the sorted union of its datasets' (windowed) wavelengths and its angles
//! the union of their angles, so each dataset's points are grid points and
//! no curve is ever interpolated; the request is their bits OR-ed.
//!
//! Materials no parameter can reach are evaluated once per grid, here, and
//! served as arrays; the rest are evaluated per evaluation, once each (a
//! small per-evaluation cache: `expand` asks per row).

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};

use num_complex::Complex64;

use crate::smatrix::angle::Radians;
use crate::smatrix::solver::solve_arrays;
use crate::structure::expansion::{ExpandOptions, expand};
use crate::structure::providers::{Entry, MaterialProvider};

use crate::fit::data::{SpectrumData, channel_request};
use crate::fit::model::FitModel;

/// Where one dataset's points sit in its environment's solution.
#[derive(Debug, Clone)]
struct Pick {
    /// Index into `StackForward::solves`.
    solve: usize,
    channel: String,
    angle: usize,
    /// Grid index of each kept point.
    wav: Vec<usize>,
}

/// One environment's solve.
#[derive(Debug, Clone)]
struct EnvSolve {
    env: String,
    grid: Vec<f64>,
    angles: Vec<Radians>,
    requested: u64,
    /// Materials no parameter reaches, evaluated on `grid`.
    frozen: HashMap<String, Vec<Complex64>>,
}

/// The compiled forward model of one problem.
#[derive(Debug, Clone)]
pub struct StackForward {
    solves: Vec<EnvSolve>,
    coherence_mode: i32,
    picks: Vec<Pick>,
}

impl StackForward {
    /// Compile `data` against `model`. `live` names the materials a
    /// parameter can change; every other material is evaluated now.
    pub fn new(
        model: &FitModel,
        data: &[SpectrumData],
        coherence_mode: i32,
        live: &HashSet<String>,
    ) -> Result<Self, String> {
        if data.is_empty() {
            return Err("fit: no data".into());
        }
        let mut solves: Vec<EnvSolve> = Vec::new();
        for d in data {
            d.check()?;
            model
                .entries(&d.env)
                .map_err(|e| format!("dataset {:?}: {e}", d.id))?;
            let s = match solves.iter().position(|s| s.env == d.env) {
                Some(s) => s,
                None => {
                    solves.push(EnvSolve {
                        env: d.env.clone(),
                        grid: Vec::new(),
                        angles: Vec::new(),
                        requested: 0,
                        frozen: HashMap::new(),
                    });
                    solves.len() - 1
                }
            };
            let s = &mut solves[s];
            s.requested |= channel_request(&d.channel)?;
            s.grid
                .extend(d.kept().into_iter().map(|i| d.wavelengths[i]));
            if !s
                .angles
                .iter()
                .any(|a| a.get().to_bits() == d.angle.get().to_bits())
            {
                s.angles.push(d.angle);
            }
        }
        for s in &mut solves {
            s.grid.sort_by(f64::total_cmp);
            s.grid.dedup_by(|a, b| a.to_bits() == b.to_bits());
            // A material that cannot be evaluated on this grid is left
            // live: it fails only if a row actually asks for it.
            for (name, entry) in model.materials() {
                if live.contains(name) {
                    continue;
                }
                if let Ok(nk) = evaluate(entry, &s.grid) {
                    s.frozen.insert(name.clone(), nk);
                }
            }
        }

        let picks = data
            .iter()
            .map(|d| {
                let solve = solves.iter().position(|s| s.env == d.env).unwrap();
                let s = &solves[solve];
                Pick {
                    solve,
                    channel: d.channel.clone(),
                    angle: s
                        .angles
                        .iter()
                        .position(|a| a.get().to_bits() == d.angle.get().to_bits())
                        .unwrap(),
                    wav: d
                        .kept()
                        .into_iter()
                        .map(|i| {
                            s.grid
                                .binary_search_by(|g| g.total_cmp(&d.wavelengths[i]))
                                .unwrap()
                        })
                        .collect(),
                }
            })
            .collect();

        Ok(Self {
            solves,
            coherence_mode,
            picks,
        })
    }

    /// The solve grid of environment `env` (nm), if a dataset uses it.
    pub fn grid(&self, env: &str) -> Option<&[f64]> {
        self.solves
            .iter()
            .find(|s| s.env == env)
            .map(|s| s.grid.as_slice())
    }

    /// Simulated values, one vector per dataset, at its kept points.
    pub fn simulate(&self, model: &FitModel) -> Result<Vec<Vec<f64>>, String> {
        let sols = self
            .solves
            .iter()
            .map(|s| {
                let provider = FitProvider {
                    frozen: &s.frozen,
                    live: model.materials(),
                    grid: &s.grid,
                    cache: RefCell::new(HashMap::new()),
                };
                let (sa, _spans) = expand(
                    model.entries(&s.env)?,
                    &provider,
                    &s.grid,
                    model.groups(),
                    ExpandOptions::deterministic(),
                )?;
                let (sol, _warnings) = solve_arrays(
                    &sa.indices,
                    &sa.thicknesses,
                    &sa.incoherent,
                    &sa.rough_types,
                    &sa.rough_vals,
                    &s.grid,
                    &s.angles,
                    s.requested,
                    self.coherence_mode,
                )?;
                Ok(sol)
            })
            .collect::<Result<Vec<_>, String>>()?;
        self.picks
            .iter()
            .map(|p| {
                let n_w = self.solves[p.solve].grid.len();
                let map = sols[p.solve]
                    .f64maps
                    .iter()
                    .find(|(k, _)| *k == p.channel)
                    .map(|(_, v)| v)
                    .ok_or_else(|| format!("fit: the solver did not return {}", p.channel))?;
                Ok(p.wav.iter().map(|&w| map[p.angle * n_w + w]).collect())
            })
            .collect()
    }
}

pub(super) fn evaluate(entry: &Entry, grid: &[f64]) -> Result<Vec<Complex64>, String> {
    match entry {
        Entry::Spec(spec) => spec.evaluate(grid),
        Entry::Array(nk) if nk.len() == grid.len() => Ok(nk.clone()),
        Entry::Array(nk) => Err(format!(
            "an nk array of {} points cannot serve a {}-point fit grid; give it as a Table spec",
            nk.len(),
            grid.len()
        )),
    }
}

/// Per-evaluation provider: frozen arrays, else the model's own entries.
struct FitProvider<'a> {
    frozen: &'a HashMap<String, Vec<Complex64>>,
    live: &'a HashMap<String, Entry>,
    grid: &'a [f64],
    cache: RefCell<HashMap<String, Vec<Complex64>>>,
}

impl MaterialProvider for FitProvider<'_> {
    fn nk(&self, name: &str, _wavelengths: &[f64]) -> Result<Vec<Complex64>, String> {
        if let Some(nk) = self.frozen.get(name) {
            return Ok(nk.clone());
        }
        if let Some(nk) = self.cache.borrow().get(name) {
            return Ok(nk.clone());
        }
        let entry = self
            .live
            .get(name)
            .ok_or_else(|| format!("fit: unknown material '{name}'"))?;
        let nk = evaluate(entry, self.grid).map_err(|e| format!("material {name}: {e}"))?;
        self.cache.borrow_mut().insert(name.to_string(), nk.clone());
        Ok(nk)
    }

    fn contains(&self, name: &str) -> bool {
        self.frozen.contains_key(name) || self.live.contains_key(name)
    }

    fn grid(&self) -> Option<&[f64]> {
        Some(self.grid)
    }

    fn names(&self) -> Vec<String> {
        self.live.keys().cloned().collect()
    }
}
