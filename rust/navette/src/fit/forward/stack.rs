// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::forward::stack — a fit model through the general solver.
//!
//! One evaluation: the model (parameters already written) → materials on
//! the solve grid → [`expand`] → [`solve_arrays`] → the channel each
//! dataset reads. Groups, grading, gradients, interfaces, roughness,
//! inversion and every coherence mode come with `expand` and the solver;
//! nothing here re-implements them.
//!
//! **One solve per evaluation.** The grid is the sorted union of every
//! dataset's (windowed) wavelengths and the angles are the union of their
//! angles, so each dataset's points are grid points and no curve is ever
//! interpolated; the request is every dataset's bit OR-ed.
//!
//! Materials no parameter can reach are evaluated once, here, and served
//! as arrays; the rest are evaluated per evaluation, once each (a small
//! per-evaluation cache: `expand` asks per row).

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};

use num_complex::Complex64;

use crate::smatrix::angle::Radians;
use crate::smatrix::solver::solve_arrays;
use crate::structure::expansion::{ExpandOptions, expand};
use crate::structure::providers::{Entry, MaterialProvider};

use crate::fit::data::{SpectrumData, channel_request};
use crate::fit::model::FitModel;

/// Where one dataset's points sit in a solution.
#[derive(Debug, Clone)]
struct Pick {
    channel: String,
    angle: usize,
    /// Grid index of each kept point.
    wav: Vec<usize>,
}

/// The compiled forward model of one problem.
#[derive(Debug, Clone)]
pub struct StackForward {
    grid: Vec<f64>,
    angles: Vec<Radians>,
    requested: u64,
    coherence_mode: i32,
    picks: Vec<Pick>,
    /// Materials no parameter reaches, evaluated on `grid`.
    frozen: HashMap<String, Vec<Complex64>>,
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
        let mut grid: Vec<f64> = Vec::new();
        let mut angles: Vec<Radians> = Vec::new();
        let mut requested = 0u64;
        for d in data {
            d.check()?;
            requested |= channel_request(&d.channel)?;
            grid.extend(d.kept().into_iter().map(|i| d.wavelengths[i]));
            if !angles
                .iter()
                .any(|a| a.get().to_bits() == d.angle.get().to_bits())
            {
                angles.push(d.angle);
            }
        }
        grid.sort_by(f64::total_cmp);
        grid.dedup_by(|a, b| a.to_bits() == b.to_bits());

        let picks = data
            .iter()
            .map(|d| Pick {
                channel: d.channel.clone(),
                angle: angles
                    .iter()
                    .position(|a| a.get().to_bits() == d.angle.get().to_bits())
                    .unwrap(),
                wav: d
                    .kept()
                    .into_iter()
                    .map(|i| {
                        grid.binary_search_by(|g| g.total_cmp(&d.wavelengths[i]))
                            .unwrap()
                    })
                    .collect(),
            })
            .collect();

        // A material that cannot be evaluated on this grid is left live: it
        // fails only if a row actually asks for it.
        let mut frozen = HashMap::new();
        for (name, entry) in model.materials() {
            if live.contains(name) {
                continue;
            }
            if let Ok(nk) = evaluate(entry, &grid) {
                frozen.insert(name.clone(), nk);
            }
        }

        Ok(Self {
            grid,
            angles,
            requested,
            coherence_mode,
            picks,
            frozen,
        })
    }

    /// The solve grid (nm).
    pub fn grid(&self) -> &[f64] {
        &self.grid
    }

    /// Simulated values, one vector per dataset, at its kept points.
    pub fn simulate(&self, model: &FitModel) -> Result<Vec<Vec<f64>>, String> {
        let provider = FitProvider {
            frozen: &self.frozen,
            live: model.materials(),
            grid: &self.grid,
            cache: RefCell::new(HashMap::new()),
        };
        let (sa, _spans) = expand(
            model.entries(),
            &provider,
            &self.grid,
            model.groups(),
            ExpandOptions::deterministic(),
        )?;
        let (sol, _warnings) = solve_arrays(
            &sa.indices,
            &sa.thicknesses,
            &sa.incoherent,
            &sa.rough_types,
            &sa.rough_vals,
            &self.grid,
            &self.angles,
            self.requested,
            self.coherence_mode,
        )?;
        let n_w = self.grid.len();
        self.picks
            .iter()
            .map(|p| {
                let map = sol
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
