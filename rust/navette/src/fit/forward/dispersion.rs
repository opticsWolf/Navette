// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::forward::dispersion — material specs on measured nk grids.
//!
//! One evaluation: each [`NkData`]'s material (parameters already written)
//! evaluated on that dataset's kept wavelengths, then the columns it
//! measures. No solve, no stack: microseconds per evaluation, so FD
//! Jacobians are cheap. It feeds the same residual vector as
//! [`StackForward`](super::StackForward), which is what makes a joint
//! fit (a film's spectrum and a literature nk table) one problem.

use num_complex::Complex64;

use crate::fit::data::{NkData, NkQuantity};
use crate::fit::model::FitModel;

use super::stack::evaluate;

/// One dataset's material, kept wavelengths and measured quantities.
#[derive(Debug, Clone)]
pub(crate) struct NkPick {
    pub(crate) material: String,
    pub(crate) wavelengths: Vec<f64>,
    pub(crate) quantities: Vec<NkQuantity>,
}

/// The compiled nk forward model of one problem.
#[derive(Debug, Clone)]
pub struct DispersionForward {
    picks: Vec<NkPick>,
}

impl DispersionForward {
    /// Compile `data` against `model`: every dataset must name one of the
    /// model's materials.
    pub fn new(model: &FitModel, data: &[NkData]) -> Result<Self, String> {
        if data.is_empty() {
            return Err("fit: no nk data".into());
        }
        let picks = data
            .iter()
            .map(|d| {
                d.check()?;
                if !model.materials().contains_key(&d.material) {
                    return Err(format!(
                        "nk dataset {:?}: the model has no material {:?}",
                        d.id, d.material
                    ));
                }
                Ok(NkPick {
                    material: d.material.clone(),
                    wavelengths: d.kept().into_iter().map(|i| d.wavelengths[i]).collect(),
                    quantities: d.columns.iter().map(|c| c.quantity).collect(),
                })
            })
            .collect::<Result<_, String>>()?;
        Ok(Self { picks })
    }

    /// Simulated values, one vector per dataset: its columns in order, each
    /// at its kept points.
    pub fn simulate(&self, model: &FitModel) -> Result<Vec<Vec<f64>>, String> {
        (0..self.picks.len())
            .map(|i| self.simulate_one(model, i))
            .collect()
    }

    pub(crate) fn picks(&self) -> &[NkPick] {
        &self.picks
    }

    /// Dataset `i`'s material on its wavelengths.
    pub(crate) fn nk_of(&self, model: &FitModel, i: usize) -> Result<Vec<Complex64>, String> {
        let p = &self.picks[i];
        let entry = &model.materials()[&p.material];
        evaluate(entry, &p.wavelengths).map_err(|e| format!("material {}: {e}", p.material))
    }

    /// Dataset `i` alone: its columns in order, each at its kept points.
    pub(crate) fn simulate_one(&self, model: &FitModel, i: usize) -> Result<Vec<f64>, String> {
        let nk = self.nk_of(model, i)?;
        Ok(self.picks[i]
            .quantities
            .iter()
            .flat_map(|q| nk.iter().map(move |&z| q.of(z)))
            .collect())
    }
}
