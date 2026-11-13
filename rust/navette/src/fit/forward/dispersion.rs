// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::forward::dispersion — material specs on measured nk grids.
//!
//! One evaluation: each [`NkData`]'s material (parameters already written)
//! evaluated on that dataset's kept wavelengths, then the columns it
//! measures. No solve, no stack: microseconds per evaluation, so FD
//! Jacobians are cheap. It feeds the same residual vector as
//! [`StackForward`](super::StackForward), which is what makes a joint
//! fit (a film's spectrum and a literature nk table) one problem.

use crate::fit::data::{NkData, NkQuantity};
use crate::fit::model::FitModel;

use super::stack::evaluate;

/// One dataset's material, kept wavelengths and measured quantities.
#[derive(Debug, Clone)]
struct NkPick {
    material: String,
    wavelengths: Vec<f64>,
    quantities: Vec<NkQuantity>,
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
        self.picks
            .iter()
            .map(|p| {
                let entry = &model.materials()[&p.material];
                let nk = evaluate(entry, &p.wavelengths)
                    .map_err(|e| format!("material {}: {e}", p.material))?;
                Ok(p.quantities
                    .iter()
                    .flat_map(|q| nk.iter().map(move |&z| q.of(z)))
                    .collect())
            })
            .collect()
    }
}
