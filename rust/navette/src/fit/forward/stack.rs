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
//! **A dataset whose instrument moves the geometry** (an angle offset, a
//! wavelength shift or a slit, fixed or fitted) is solved on its own, on
//! the grid and angle its instrument gives at this evaluation: its points
//! shifted, each spread over the slit's nodes, the result averaged back.
//! Every dataset's reading then passes through `scale · value + offset`.
//!
//! Materials no parameter can reach are evaluated once per grid, here, and
//! served as arrays; the rest are evaluated per evaluation, once each (a
//! small per-evaluation cache: `expand` asks per row).

use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};

use num_complex::Complex64;

use crate::smatrix::angle::Radians;
use crate::smatrix::solver::solve_arrays;
use crate::structure::expansion::{ExpandOptions, SolverArrays, expand};
use crate::structure::providers::{Entry, MaterialProvider};

use crate::fit::data::{Instrument, SpectrumData, channel_request, is_intensity};
use crate::fit::model::FitModel;

/// Where one dataset's points sit in its environment's solution.
#[derive(Debug, Clone)]
pub(crate) struct Pick {
    /// Index into `StackForward::solves`.
    pub(crate) solve: usize,
    /// The dataset's id (its instrument in the model).
    pub(crate) id: String,
    pub(crate) channel: String,
    pub(crate) angle: usize,
    /// Grid index of each kept point (a dataset solved on its own: the
    /// point's index, its slit nodes following).
    pub(crate) wav: Vec<usize>,
}

/// A dataset solved on its own grid.
#[derive(Debug, Clone)]
pub(crate) struct Own {
    /// The dataset's id.
    id: String,
    /// Slit nodes per point: 1, or 7 when the slit is (or may become)
    /// wider than zero.
    nodes: usize,
}

/// One environment's solve.
#[derive(Debug, Clone)]
pub(crate) struct EnvSolve {
    pub(crate) env: String,
    /// The solve grid; for a dataset solved on its own, its kept
    /// wavelengths as read (the instrument moves them per evaluation).
    pub(crate) grid: Vec<f64>,
    /// The solve angles; for a dataset solved on its own, its angle as read.
    pub(crate) angles: Vec<Radians>,
    requested: u64,
    /// Materials no parameter reaches, evaluated on `grid`.
    frozen: HashMap<String, Vec<Complex64>>,
    /// `Some` for a dataset solved on its own.
    pub(crate) own: Option<Own>,
}

/// Where and how one solve runs at one model.
struct Geometry<'a> {
    grid: Cow<'a, [f64]>,
    angles: Cow<'a, [Radians]>,
    /// Slit weight of each grid point (a dataset solved on its own).
    weights: Option<Vec<f64>>,
    nodes: usize,
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
    /// `moving` names the datasets whose instrument geometry a parameter
    /// can change: those, and any whose instrument already departs from
    /// the curve's geometry, are solved on their own.
    pub fn new(
        model: &FitModel,
        data: &[SpectrumData],
        coherence_mode: i32,
        live: &HashSet<String>,
        moving: &HashSet<String>,
    ) -> Result<Self, String> {
        if data.is_empty() {
            return Err("fit: no data".into());
        }
        let mut solves: Vec<EnvSolve> = Vec::new();
        let mut picks = Vec::with_capacity(data.len());
        for d in data {
            d.check()?;
            model
                .entries(&d.env)
                .map_err(|e| format!("dataset {:?}: {e}", d.id))?;
            let instrument = model.instrument(&d.id).copied().unwrap_or(d.instrument);
            instrument.check(&d.id, &d.channel)?;
            let request = channel_request(&d.channel)?;
            let kept = d.kept();
            if moving.contains(&d.id) || instrument.moves_geometry() {
                let slit = instrument.bandwidth != 0.0
                    || (moving.contains(&d.id) && is_intensity(&d.channel));
                solves.push(EnvSolve {
                    env: d.env.clone(),
                    grid: kept.iter().map(|&i| d.wavelengths[i]).collect(),
                    angles: vec![d.angle],
                    requested: request,
                    frozen: HashMap::new(),
                    own: Some(Own {
                        id: d.id.clone(),
                        nodes: if slit { 7 } else { 1 },
                    }),
                });
                picks.push(Pick {
                    solve: solves.len() - 1,
                    id: d.id.clone(),
                    channel: d.channel.clone(),
                    angle: 0,
                    wav: (0..kept.len()).collect(),
                });
                continue;
            }
            let s = match solves
                .iter()
                .position(|s| s.env == d.env && s.own.is_none())
            {
                Some(s) => s,
                None => {
                    solves.push(EnvSolve {
                        env: d.env.clone(),
                        grid: Vec::new(),
                        angles: Vec::new(),
                        requested: 0,
                        frozen: HashMap::new(),
                        own: None,
                    });
                    solves.len() - 1
                }
            };
            let es = &mut solves[s];
            es.requested |= request;
            es.grid.extend(kept.iter().map(|&i| d.wavelengths[i]));
            if !es
                .angles
                .iter()
                .any(|a| a.get().to_bits() == d.angle.get().to_bits())
            {
                es.angles.push(d.angle);
            }
            picks.push(Pick {
                solve: s,
                id: d.id.clone(),
                channel: d.channel.clone(),
                angle: 0,
                wav: Vec::new(),
            });
        }
        for s in solves.iter_mut().filter(|s| s.own.is_none()) {
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
        for (p, d) in picks.iter_mut().zip(data) {
            let s = &solves[p.solve];
            if s.own.is_some() {
                continue;
            }
            p.angle = s
                .angles
                .iter()
                .position(|a| a.get().to_bits() == d.angle.get().to_bits())
                .unwrap();
            p.wav = d
                .kept()
                .into_iter()
                .map(|i| {
                    s.grid
                        .binary_search_by(|g| g.total_cmp(&d.wavelengths[i]))
                        .unwrap()
                })
                .collect();
        }

        Ok(Self {
            solves,
            coherence_mode,
            picks,
        })
    }

    /// The solve grid of environment `env` (nm), if a dataset uses it
    /// (the shared one; a dataset solved on its own is not listed).
    pub fn grid(&self, env: &str) -> Option<&[f64]> {
        self.solves
            .iter()
            .find(|s| s.env == env && s.own.is_none())
            .map(|s| s.grid.as_slice())
    }

    /// Simulated values, one vector per dataset, at its kept points.
    pub fn simulate(&self, model: &FitModel) -> Result<Vec<Vec<f64>>, String> {
        let mut out: Vec<Vec<f64>> = vec![Vec::new(); self.picks.len()];
        for s in 0..self.solves.len() {
            for (d, values) in self.simulate_solve(model, s)? {
                out[d] = values;
            }
        }
        Ok(out)
    }

    /// The environments' solves, in order.
    pub(crate) fn solves(&self) -> &[EnvSolve] {
        &self.solves
    }

    /// Where each dataset reads its environment's solution.
    pub(crate) fn picks(&self) -> &[Pick] {
        &self.picks
    }

    /// Solve `s`'s grid and angles at `model`.
    fn geometry(&self, model: &FitModel, s: usize) -> Result<Geometry<'_>, String> {
        let es = &self.solves[s];
        let Some(own) = &es.own else {
            return Ok(Geometry {
                grid: Cow::Borrowed(&es.grid),
                angles: Cow::Borrowed(&es.angles),
                weights: None,
                nodes: 1,
            });
        };
        let i = instrument(model, &own.id);
        let (grid, weights) = i.sample(&es.grid, own.nodes);
        if let Some(w) = grid.iter().find(|w| !(**w > 0.0)) {
            return Err(format!(
                "dataset {:?}: the instrument puts a point at {w} nm",
                own.id
            ));
        }
        let angle = Radians::new(es.angles[0].get() + i.angle_offset)
            .map_err(|e| format!("dataset {:?}: angle offset: {e}", own.id))?;
        Ok(Geometry {
            grid: Cow::Owned(grid),
            angles: Cow::Owned(vec![angle]),
            weights: Some(weights),
            nodes: own.nodes,
        })
    }

    /// Environment `s`'s chain expanded at `model`, with the material names
    /// the expansion asked for.
    pub(crate) fn expand_solve(
        &self,
        model: &FitModel,
        s: usize,
    ) -> Result<(SolverArrays, HashSet<String>), String> {
        let es = &self.solves[s];
        let g = self.geometry(model, s)?;
        let provider = FitProvider {
            frozen: &es.frozen,
            live: model.materials(),
            grid: &g.grid,
            cache: RefCell::new(HashMap::new()),
            asked: RefCell::new(HashSet::new()),
        };
        let (sa, _spans) = expand(
            model.entries(&es.env)?,
            &provider,
            &g.grid,
            model.groups(),
            ExpandOptions::deterministic(),
        )?;
        Ok((sa, provider.asked.into_inner()))
    }

    /// One environment's datasets at `model`: `(dataset index, values)`
    /// for each dataset this solve serves, in dataset order.
    pub(crate) fn simulate_solve(
        &self,
        model: &FitModel,
        s: usize,
    ) -> Result<Vec<(usize, Vec<f64>)>, String> {
        let (sa, _) = self.expand_solve(model, s)?;
        self.read_solve(model, &sa, s)
    }

    /// [`Self::simulate_solve`] from arrays expanded at `model`.
    pub(crate) fn read_solve(
        &self,
        model: &FitModel,
        sa: &SolverArrays,
        s: usize,
    ) -> Result<Vec<(usize, Vec<f64>)>, String> {
        let es = &self.solves[s];
        let g = self.geometry(model, s)?;
        let (sol, _warnings) = solve_arrays(
            &sa.indices,
            &sa.thicknesses,
            &sa.incoherent,
            &sa.rough_types,
            &sa.rough_vals,
            &g.grid,
            &g.angles,
            es.requested,
            self.coherence_mode,
        )?;
        let n_w = g.grid.len();
        self.picks
            .iter()
            .enumerate()
            .filter(|(_, p)| p.solve == s)
            .map(|(d, p)| {
                let map = sol
                    .f64maps
                    .iter()
                    .find(|(k, _)| *k == p.channel)
                    .map(|(_, v)| v)
                    .ok_or_else(|| format!("fit: the solver did not return {}", p.channel))?;
                let i = instrument(model, &p.id);
                // The identity reading leaves the solver's bits alone.
                let read = |v: f64| {
                    if i.scale == 1.0 && i.offset == 0.0 {
                        v
                    } else {
                        i.scale * v + i.offset
                    }
                };
                let values = match &g.weights {
                    None => p
                        .wav
                        .iter()
                        .map(|&w| read(map[p.angle * n_w + w]))
                        .collect(),
                    Some(weights) => p
                        .wav
                        .iter()
                        .map(|&w| {
                            let at = w * g.nodes..(w + 1) * g.nodes;
                            read(at.map(|k| weights[k] * map[k]).sum())
                        })
                        .collect(),
                };
                Ok((d, values))
            })
            .collect()
    }
}

/// Dataset `id`'s instrument at `model`; the identity before any seeding.
fn instrument(model: &FitModel, id: &str) -> Instrument {
    model.instrument(id).copied().unwrap_or_default()
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
    /// Every name `expand` asked for: what this environment depends on.
    asked: RefCell<HashSet<String>>,
}

impl MaterialProvider for FitProvider<'_> {
    fn nk(&self, name: &str, _wavelengths: &[f64]) -> Result<Vec<Complex64>, String> {
        if !self.asked.borrow().contains(name) {
            self.asked.borrow_mut().insert(name.to_string());
        }
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
