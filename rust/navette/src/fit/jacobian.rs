// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::jacobian — the fit's own Jacobian (plan §3.7, phase 6).
//!
//! Central differences over the whole residual vector cost two full
//! evaluations per column: every environment expanded and solved, every
//! nk table evaluated, whatever the column can actually move. This source
//! builds the same matrix from less:
//!
//! - **Reach.** A column recomputes only the blocks its parameter can
//!   move. A layer reaches the environments whose chain holds it; a group
//!   or a material, the environments whose expansion asks for that name;
//!   a material also its own nk tables; a tie driver what its followers
//!   reach. Every other row of the column is zero, which is exactly what
//!   central differences give there (the same inputs, the same bits).
//!   Blocks it does reach are differenced with the same step and the same
//!   arithmetic as `lsq::lm`'s central differences, so those entries are
//!   bit-identical to the old Jacobian.
//! - **Thickness columns, analytically.** When a column moves only the
//!   slot thicknesses of an environment's expansion (a film's thickness,
//!   a group's thickness factor or summand, graded and sliced layers
//!   included: the expansion is linear in them), and the environment is
//!   coherent and read through R, T or A, the column is the chain
//!   `Σ_slot ∂X/∂d_slot · ∂d_slot/∂u`. The first factor is the needle
//!   operator's thickness slope (`smatrix::needle_operator`), read off one
//!   field decomposition per point for every slot at once; the second is
//!   the expansion differenced (linear, so exact to rounding). The whole
//!   thickness block costs about one solve instead of `2n`.
//! - **Dispersion columns, analytically.** A top-level number of a
//!   Konstant, Cauchy or Sellmeier material differentiates in closed form
//!   on its nk tables.
//!
//! Anything the source cannot express it differences, block by block; a
//! block that fails to evaluate at a probe hands the whole Jacobian back
//! to the driver's central differences, so refused trial points behave as
//! before. `LmConfig::jacobian = Fd` skips this source entirely.

use std::collections::HashSet;

use num_complex::Complex64;
use rayon::prelude::*;

use crate::lsq::lm::JacobianSource;
use crate::smatrix::needle_operator::{
    block_flux_factors, build_stack_fields_range_dual, thickness_slopes_ddz,
};
use crate::structure::expansion::SolverArrays;
use crate::structure::{Entry, MaterialSpec};

use super::model::{FitModel, Target};
use super::{Compiled, Watcher, data};

/// The derivative of a channel the analytic chain serves, from the four
/// intensity slopes `[dRs, dRp, dTs, dTp]`.
fn channel_slope(channel: &str, d: &[f64; 4]) -> Option<f64> {
    Some(match channel {
        "Rs" => d[0],
        "Rp" => d[1],
        "Ts" => d[2],
        "Tp" => d[3],
        "R_avg" => 0.5 * (d[0] + d[1]),
        "T_avg" => 0.5 * (d[2] + d[3]),
        "A_s" => -(d[0] + d[2]),
        "A_p" => -(d[1] + d[3]),
        "A_avg" => -0.5 * (d[0] + d[1] + d[2] + d[3]),
        _ => return None,
    })
}

/// `[dRs, dRp, dTs, dTp]` per point (`angle · n_w + w`, the solver's
/// order) and per slot of `slots`, at the expansion `sa`. One field
/// decomposition per point and polarization serves every slot.
fn intensity_slopes(
    sa: &SolverArrays,
    grid: &[f64],
    sin_theta: &[f64],
    slots: &[usize],
) -> Vec<[f64; 4]> {
    let nl = sa.n_rows();
    let nw = grid.len();
    let ns_slots = slots.len();
    let per_point: Vec<Vec<[f64; 4]>> = (0..sin_theta.len() * nw)
        .into_par_iter()
        .map(|k| {
            let (a, w) = (k / nw, k % nw);
            let lam = grid[w];
            let ns: Vec<Complex64> = (0..nl).map(|l| sa.indices[l * nw + w]).collect();
            let nsin_fi = ns[0] * Complex64::new(sin_theta[a], 0.0);
            let (fs, fp) = build_stack_fields_range_dual(
                0,
                nl - 1,
                &ns,
                &sa.thicknesses,
                &sa.rough_vals,
                &sa.rough_types,
                lam,
                nsin_fi,
            );
            let mut out = vec![[0.0; 4]; ns_slots];
            for (pol, fields) in [(0usize, &fs), (1usize, &fp)] {
                let m = fields.s_left[fields.end];
                let f_fwd = block_flux_factors(fields, pol as i32)[2];
                let (r_conj, t_conj) = (m.0.conj(), m.2.conj());
                for (o, &j) in out.iter_mut().zip(slots) {
                    let [dr, dt] = thickness_slopes_ddz(fields, j, lam);
                    o[pol] = 2.0 * (r_conj * dr).re;
                    o[2 + pol] = 2.0 * f_fwd * (t_conj * dt).re;
                }
            }
            out
        })
        .collect();
    per_point.into_iter().flatten().collect()
}

/// `∂(n + ik)/∂p` on `wavelengths` for a top-level number `key` of
/// `spec`, at `nk` (the spec's own values there), in closed form; `None`
/// for a model or a key this does not cover.
fn nk_slope(
    spec: &MaterialSpec,
    key: &str,
    wavelengths: &[f64],
    nk: &[Complex64],
) -> Option<Vec<Complex64>> {
    let num = |k: &str| spec.params.get(k).and_then(|v| v.as_f64());
    let um2 = |w: f64| {
        let um = w * 1.0e-3;
        um * um
    };
    let real = |f: &dyn Fn(f64, Complex64) -> f64| -> Vec<Complex64> {
        wavelengths
            .iter()
            .zip(nk)
            .map(|(&w, &z)| Complex64::new(f(um2(w), z), 0.0))
            .collect()
    };
    match (spec.model.as_str(), key) {
        ("Konstant", "n") => Some(vec![Complex64::new(1.0, 0.0); wavelengths.len()]),
        ("Konstant", "k") => Some(vec![Complex64::new(0.0, 1.0); wavelengths.len()]),
        ("Cauchy", "A") => Some(real(&|_, _| 1.0)),
        ("Cauchy", "B") => Some(real(&|l2, _| 1.0 / l2)),
        ("Cauchy", "C") => Some(real(&|l2, _| 1.0 / (l2 * l2))),
        ("Drude", "omega_p" | "gamma" | "epsilon_inf") => {
            // eps = eps_inf − wp²/(E² + iγE), E guarded as the model does;
            // ∂(n + ik) = ∂eps/(2(n + ik)).
            let wp = num("omega_p")?;
            let gamma = num("gamma")?;
            Some(
                wavelengths
                    .iter()
                    .zip(nk)
                    .map(|(&w, &z)| {
                        let e = crate::materials::units::energy_ev(w) + 1.0e-12;
                        let d = Complex64::new(e * e, gamma * e);
                        let deps = match key {
                            "omega_p" => -2.0 * wp / d,
                            "gamma" => wp * wp * Complex64::new(0.0, e) / (d * d),
                            _ => Complex64::new(1.0, 0.0),
                        };
                        deps / (2.0 * z)
                    })
                    .collect(),
            )
        }
        ("Sellmeier", k) => {
            // n² = 1 + Σ B·l2/(l2 − C): ∂n = ∂(n²)/(2n).
            let i = match k.get(1..) {
                Some("1") => 1,
                Some("2") => 2,
                Some("3") => 3,
                _ => return None,
            };
            let b = num(&format!("B{i}"))?;
            let c = num(&format!("C{i}"))?;
            match k.get(..1).unwrap_or("") {
                "B" => Some(real(&|l2, z| l2 / (l2 - c) / (2.0 * z.re))),
                // The third term counts only when B3 is nonzero.
                "C" if i == 3 && b == 0.0 => Some(real(&|_, _| 0.0)),
                "C" => Some(real(&|l2, z| b * l2 / ((l2 - c) * (l2 - c)) / (2.0 * z.re))),
                _ => None,
            }
        }
        _ => None,
    }
}

/// How one column moves one environment.
enum SolveColumn {
    /// The residual rows of the environment's datasets, differenced.
    Rows(Vec<(usize, f64)>),
    /// Only slot thicknesses move: `∂d_slot/∂u` per moving slot.
    Thickness(Vec<(usize, f64)>),
}

/// What a column computes, before the analytic sweeps are folded in.
struct Column {
    /// Rows already final (differenced blocks, closed-form nk, priors).
    rows: Vec<(usize, f64)>,
    /// `(solve, ∂d_slot/∂u)` for the environments read analytically.
    thickness: Vec<(usize, Vec<(usize, f64)>)>,
}

/// The fit's [`JacobianSource`]: see the module docs.
pub(super) struct FitJacobian<'a, 'w> {
    c: &'a Compiled,
    watch: Option<&'a Watcher<'w>>,
    /// Per free parameter: the solves and the nk datasets it moves.
    reach_solves: Vec<Vec<usize>>,
    reach_nk: Vec<Vec<usize>>,
    /// A free parameter drives tie followers.
    drives: Vec<bool>,
    /// Per solve: every dataset reads a channel the chain serves.
    chain_ok: Vec<bool>,
    /// First residual row of each dataset (spectra, then nk), then of the
    /// priors; and the total.
    offsets: Vec<usize>,
    prior_offset: usize,
    smooth_offset: usize,
    /// Per free parameter: it (or a follower) is a smoothed knot.
    smooths: Vec<bool>,
    m: usize,
}

impl<'a, 'w> FitJacobian<'a, 'w> {
    pub(super) fn new(c: &'a Compiled, watch: Option<&'a Watcher<'w>>) -> Result<Self, String> {
        let mut offsets = Vec::with_capacity(c.rows.len());
        let mut m = 0;
        for r in &c.rows {
            offsets.push(m);
            m += r.values.len();
        }
        let prior_offset = m;
        m += c.n_priors();
        let smooth_offset = m;
        m += c.n_smooth();

        // What each environment depends on: its rows' materials, gradient
        // ends, and whatever its expansion asked the provider for.
        let solves = c.stack.as_ref().map(|s| s.solves()).unwrap_or(&[]);
        let mut names: Vec<HashSet<String>> = Vec::with_capacity(solves.len());
        let mut ranges = Vec::with_capacity(solves.len());
        for (s, es) in solves.iter().enumerate() {
            let mut set = HashSet::new();
            for (layer, _) in c.base.entries(&es.env)? {
                set.insert(layer.material.clone());
                if let Some(g) = &layer.gradient {
                    set.insert(g.material_a.clone());
                    set.insert(g.material_b.clone());
                }
            }
            let (_, asked) = c.stack.as_ref().unwrap().expand_solve(&c.base, s)?;
            set.extend(asked);
            names.push(set);
            ranges.push(
                c.base
                    .env_range(&es.env)
                    .ok_or("fit: unknown environment")?,
            );
        }
        let nk_materials: Vec<&str> = c
            .dispersion
            .as_ref()
            .map(|d| d.picks().iter().map(|p| p.material.as_str()).collect())
            .unwrap_or_default();
        let reach_of = |t: &Target, solves_out: &mut Vec<usize>, nk_out: &mut Vec<usize>| {
            for s in 0..names.len() {
                let hit = match t {
                    Target::Layer { entries, .. } => entries.iter().any(|e| ranges[s].contains(e)),
                    Target::Group { group, .. } => names[s].contains(group),
                    Target::Material { material, .. } => names[s].contains(material),
                    Target::Data { dataset, .. } => c
                        .stack
                        .as_ref()
                        .unwrap()
                        .picks()
                        .iter()
                        .any(|p| p.solve == s && p.id == *dataset),
                };
                if hit && !solves_out.contains(&s) {
                    solves_out.push(s);
                }
            }
            if let Target::Material { material, .. } = t {
                for (i, m) in nk_materials.iter().enumerate() {
                    if m == material && !nk_out.contains(&i) {
                        nk_out.push(i);
                    }
                }
            }
        };
        let mut reach_solves = Vec::with_capacity(c.free.len());
        let mut reach_nk = Vec::with_capacity(c.free.len());
        let mut drives = Vec::with_capacity(c.free.len());
        let smooths = c
            .free
            .iter()
            .map(|f| {
                c.smooths(&f.target)
                    || c.ties
                        .iter()
                        .filter(|t| t.driver == f.target)
                        .any(|t| t.targets.iter().any(|(x, _, _)| c.smooths(x)))
            })
            .collect();
        for f in &c.free {
            let (mut rs, mut rn) = (Vec::new(), Vec::new());
            reach_of(&f.target, &mut rs, &mut rn);
            let mut driver = false;
            for tie in c.ties.iter().filter(|t| t.driver == f.target) {
                driver = true;
                for (t, _, _) in &tie.targets {
                    reach_of(t, &mut rs, &mut rn);
                }
            }
            rs.sort_unstable();
            rn.sort_unstable();
            reach_solves.push(rs);
            reach_nk.push(rn);
            drives.push(driver);
        }
        // A dataset solved on its own reads a moving grid: differenced.
        let chain_ok = (0..solves.len())
            .map(|s| {
                let picks = c.stack.as_ref().unwrap().picks();
                solves[s].own.is_none()
                    && picks
                        .iter()
                        .filter(|p| p.solve == s)
                        .all(|p| channel_slope(&p.channel, &[0.0; 4]).is_some())
            })
            .collect();
        Ok(Self {
            c,
            watch,
            reach_solves,
            reach_nk,
            drives,
            chain_ok,
            offsets,
            prior_offset,
            smooth_offset,
            smooths,
            m,
        })
    }

    /// Residual rows of dataset `d` from simulated values, exactly as
    /// `Compiled::residuals` forms them.
    fn rows_of<'s>(&'s self, d: usize, sim: &'s [f64]) -> impl Iterator<Item = f64> + 's {
        let r = &self.c.rows[d];
        sim.iter()
            .zip(&r.values)
            .zip(&r.inv_sigma)
            .map(move |((s, v), w)| {
                let dv = s - v;
                (if r.angular { data::wrap_angle(dv) } else { dv }) * w
            })
    }

    /// `(rp − rm)·inv2h` per row of dataset `d`, into `out`.
    fn difference(
        &self,
        d: usize,
        sp: &[f64],
        sm: &[f64],
        inv2h: f64,
        out: &mut Vec<(usize, f64)>,
    ) {
        let o = self.offsets[d];
        for (i, (rp, rm)) in self.rows_of(d, sp).zip(self.rows_of(d, sm)).enumerate() {
            out.push((o + i, (rp - rm) * inv2h));
        }
    }

    /// Column `j` at `u`, the analytic thickness parts left open. `None`
    /// when a probe fails: the driver then differences everything.
    fn column(
        &self,
        u: &[f64],
        j: usize,
        base: &FitModel,
        base_sa: &[Option<SolverArrays>],
    ) -> Option<Column> {
        let c = self.c;
        let f = &c.free[j];
        let h = f64::EPSILON.cbrt() * u[j].abs().max(1.0);
        let inv2h = 1.0 / (2.0 * h);
        let (mut up, mut um) = (u.to_vec(), u.to_vec());
        up[j] += h;
        um[j] -= h;
        let mut col = Column {
            rows: Vec::new(),
            thickness: Vec::new(),
        };
        // Closed-form nk blocks first: a column they cover entirely needs
        // no probe models.
        let n_spectra = c.stack.as_ref().map_or(0, |s| s.picks().len());
        let mut nk_differenced = Vec::new();
        for &i in &self.reach_nk[j] {
            match self.nk_closed_form(base, j, i, u) {
                Some(rows) => {
                    let o = self.offsets[n_spectra + i];
                    col.rows.extend(rows.into_iter().map(|(r, v)| (o + r, v)));
                }
                None => nk_differenced.push(i),
            }
        }
        let probes = if self.reach_solves[j].is_empty() && nk_differenced.is_empty() {
            None
        } else {
            Some((c.model_at(&up).ok()?, c.model_at(&um).ok()?))
        };
        if let (Some(stack), Some((mp, mm))) = (&c.stack, &probes) {
            for &s in &self.reach_solves[j] {
                let (sap, _) = stack.expand_solve(mp, s).ok()?;
                let (sam, _) = stack.expand_solve(mm, s).ok()?;
                // An instrument moves what the arrays cannot show.
                let how = match &base_sa[s] {
                    Some(b) if !matches!(f.target, Target::Data { .. }) => {
                        thickness_only(b, &sap, &sam, inv2h)
                    }
                    _ => None,
                };
                let how = match how {
                    Some(dt) => SolveColumn::Thickness(dt),
                    None => {
                        let (p, m) = rayon::join(
                            || stack.read_solve(mp, &sap, s),
                            || stack.read_solve(mm, &sam, s),
                        );
                        let (p, m) = (p.ok()?, m.ok()?);
                        let mut rows = Vec::new();
                        for ((d, sp), (_, sm)) in p.iter().zip(&m) {
                            self.difference(*d, sp, sm, inv2h, &mut rows);
                        }
                        SolveColumn::Rows(rows)
                    }
                };
                match how {
                    SolveColumn::Rows(rows) => col.rows.extend(rows),
                    SolveColumn::Thickness(dt) => col.thickness.push((s, dt)),
                }
            }
        }
        if let (Some(disp), Some((mp, mm))) = (&c.dispersion, &probes) {
            for &i in &nk_differenced {
                let sp = disp.simulate_one(mp, i).ok()?;
                let sm = disp.simulate_one(mm, i).ok()?;
                self.difference(n_spectra + i, &sp, &sm, inv2h, &mut col.rows);
            }
        }
        if let Some((p0, s)) = f.prior {
            let k = c.free[..j].iter().filter(|g| g.prior.is_some()).count();
            let rp = (f.x(up[j]) - p0) / s;
            let rm = (f.x(um[j]) - p0) / s;
            col.rows.push((self.prior_offset + k, (rp - rm) * inv2h));
        }
        if self.smooths[j] {
            let owned;
            let (mp, mm) = match &probes {
                Some((mp, mm)) => (mp, mm),
                None => {
                    owned = (c.model_at(&up).ok()?, c.model_at(&um).ok()?);
                    (&owned.0, &owned.1)
                }
            };
            let (mut rp, mut rm) = (Vec::new(), Vec::new());
            c.smooth_rows(mp, &mut rp).ok()?;
            c.smooth_rows(mm, &mut rm).ok()?;
            for (i, (p, m)) in rp.iter().zip(&rm).enumerate() {
                if p != m {
                    col.rows.push((self.smooth_offset + i, (p - m) * inv2h));
                }
            }
        }
        Some(col)
    }

    /// Column `j`'s rows on nk dataset `i` in closed form, row index within
    /// the dataset; `None` when the model or the parameter is not covered.
    fn nk_closed_form(
        &self,
        base: &FitModel,
        j: usize,
        i: usize,
        u: &[f64],
    ) -> Option<Vec<(usize, f64)>> {
        if self.drives[j] {
            return None;
        }
        let c = self.c;
        let f = &c.free[j];
        let Target::Material {
            material,
            key,
            pointer,
        } = &f.target
        else {
            return None;
        };
        if !pointer.is_empty() {
            return None;
        }
        let Entry::Spec(spec) = base.materials().get(material)? else {
            return None;
        };
        let disp = c.dispersion.as_ref()?;
        let pick = &disp.picks()[i];
        let nk = disp.nk_of(base, i).ok()?;
        let dnk = nk_slope(spec, key, &pick.wavelengths, &nk)?;
        let x = f.x(u[j]);
        let dxdu = f.scale * f.transform.slope(x);
        let n_spectra = c.stack.as_ref().map_or(0, |s| s.picks().len());
        let inv_sigma = &c.rows[n_spectra + i].inv_sigma;
        let npts = nk.len();
        let mut out = Vec::with_capacity(pick.quantities.len() * npts);
        for (qi, q) in pick.quantities.iter().enumerate() {
            for (p, (&z, &dz)) in nk.iter().zip(&dnk).enumerate() {
                let r = qi * npts + p;
                out.push((r, q.slope(z, dz) * dxdu * inv_sigma[r]));
            }
        }
        out.iter().all(|(_, v)| v.is_finite()).then_some(out)
    }
}

/// `∂d_slot/∂u` when the two probes differ from `base` only in slot
/// thicknesses, and only inside the stack (row 0 and the last are the
/// half-spaces); `None` otherwise.
fn thickness_only(
    base: &SolverArrays,
    p: &SolverArrays,
    m: &SolverArrays,
    inv2h: f64,
) -> Option<Vec<(usize, f64)>> {
    let same = |a: &SolverArrays| {
        a.thicknesses.len() == base.thicknesses.len()
            && a.incoherent == base.incoherent
            && a.rough_types == base.rough_types
            && a.rough_vals == base.rough_vals
            && a.indices == base.indices
    };
    if !(same(p) && same(m)) {
        return None;
    }
    let nl = base.thicknesses.len();
    let mut dt = Vec::new();
    for (slot, (tp, tm)) in p.thicknesses.iter().zip(&m.thicknesses).enumerate() {
        let d = (tp - tm) * inv2h;
        if d != 0.0 {
            if slot == 0 || slot + 1 == nl {
                return None;
            }
            dt.push((slot, d));
        }
    }
    Some(dt)
}

impl JacobianSource for FitJacobian<'_, '_> {
    fn fill(&self, u: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        if self.watch.is_some_and(|w| w.is_stopped()) {
            return Ok(None);
        }
        let c = self.c;
        let n = c.free.len();
        let Ok(base) = c.model_at(u) else {
            return Ok(None);
        };
        // The expansion at `u` of each environment the chain may read.
        let solves = c.stack.as_ref().map(|s| s.solves()).unwrap_or(&[]);
        let mut base_sa: Vec<Option<SolverArrays>> = Vec::with_capacity(solves.len());
        for s in 0..solves.len() {
            let wanted = self.chain_ok[s] && self.reach_solves.iter().any(|r| r.contains(&s));
            let sa = if wanted {
                match c.stack.as_ref().unwrap().expand_solve(&base, s) {
                    Ok((sa, _)) => {
                        let nl = sa.n_rows();
                        let interior = if nl > 2 {
                            &sa.incoherent[1..nl - 1]
                        } else {
                            &[][..]
                        };
                        (nl >= 2 && !interior.contains(&true)).then_some(sa)
                    }
                    Err(_) => return Ok(None),
                }
            } else {
                None
            };
            base_sa.push(sa);
        }
        let columns: Vec<Option<Column>> = (0..n)
            .into_par_iter()
            .map(|j| self.column(u, j, &base, &base_sa))
            .collect();
        let Some(columns) = columns.into_iter().collect::<Option<Vec<Column>>>() else {
            return Ok(None);
        };

        jac.clear();
        jac.resize(self.m * n, 0.0);
        for (j, col) in columns.iter().enumerate() {
            for &(r, v) in &col.rows {
                jac[r * n + j] = v;
            }
        }

        // The analytic environments: one sweep each over the slots any
        // column moves, then every column's chain.
        let stack = c.stack.as_ref();
        for (s, sa) in base_sa.iter().enumerate() {
            let Some(sa) = sa else { continue };
            let mut slots: Vec<usize> = columns
                .iter()
                .flat_map(|col| col.thickness.iter().filter(|(t, _)| *t == s))
                .flat_map(|(_, dt)| dt.iter().map(|(slot, _)| *slot))
                .collect();
            if slots.is_empty() {
                continue;
            }
            slots.sort_unstable();
            slots.dedup();
            let es = &solves[s];
            let sin_theta: Vec<f64> = es.angles.iter().map(|a| a.get().sin()).collect();
            let slopes = intensity_slopes(sa, &es.grid, &sin_theta, &slots);
            let ns = slots.len();
            let nw = es.grid.len();
            let picks = stack.unwrap().picks();
            for (j, col) in columns.iter().enumerate() {
                for (_, dt) in col.thickness.iter().filter(|(t, _)| *t == s) {
                    let at: Vec<(usize, f64)> = dt
                        .iter()
                        .map(|(slot, d)| (slots.binary_search(slot).unwrap(), *d))
                        .collect();
                    for (d, p) in picks.iter().enumerate().filter(|(_, p)| p.solve == s) {
                        let inv_sigma = &c.rows[d].inv_sigma;
                        let scale = base.instrument(&p.id).map_or(1.0, |i| i.scale);
                        for (q, &w) in p.wav.iter().enumerate() {
                            let k = p.angle * nw + w;
                            let mut v = 0.0;
                            for &(i, dd) in &at {
                                v += channel_slope(&p.channel, &slopes[k * ns + i]).unwrap() * dd;
                            }
                            jac[(self.offsets[d] + q) * n + j] = scale * v * inv_sigma[q];
                        }
                    }
                }
            }
        }
        if jac.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        Ok(Some(self.m))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, HashMap};

    use serde_json::json;

    use super::FitJacobian;
    use crate::fit::model::Segment;
    use crate::fit::{
        Compiled, FitModel, FitProblem, NkColumn, NkData, NkQuantity, Parameter, SpectrumData, Tie,
    };
    use crate::lsq::lm::{JacobianSource, build_jacobian};
    use crate::smatrix::angle::Radians;
    use crate::structure::enums::{LayerType, RoughnessType};
    use crate::structure::{Entry, Group, Layer, MaterialSpec, Structure};

    fn spec(model: &str, params: serde_json::Value) -> Entry {
        let params: BTreeMap<String, serde_json::Value> = serde_json::from_value(params).unwrap();
        Entry::Spec(MaterialSpec::new(model, params))
    }

    fn half_space(material: &str, role: LayerType) -> Layer {
        Layer {
            layer_type: role,
            ..Layer::film(0.0, material)
        }
    }

    fn wavelengths() -> Vec<f64> {
        (0..31).map(|i| 400.0 + 12.0 * i as f64).collect()
    }

    /// air | TiO2 80 | SiO2 120 | TiO2 60 | glass, a group per material,
    /// plus a bare-glass environment.
    fn model() -> FitModel {
        model_with(vec![
            Layer::film(80.0, "TiO2"),
            Layer::film(120.0, "SiO2"),
            Layer::film(60.0, "TiO2"),
        ])
    }

    /// [`model`] with rough interfaces: Névot-Croce under the SiO2 film,
    /// Gaussian under the last TiO2.
    fn rough_model() -> FitModel {
        model_with(vec![
            Layer::film(80.0, "TiO2"),
            Layer {
                rough_type: RoughnessType::NevotCroce,
                roughness: 3.0,
                ..Layer::film(120.0, "SiO2")
            },
            Layer {
                rough_type: RoughnessType::Gaussian,
                roughness: 2.0,
                ..Layer::film(60.0, "TiO2")
            },
        ])
    }

    fn model_with(films: Vec<Layer>) -> FitModel {
        let mut layers = vec![half_space("air", LayerType::Ambient)];
        layers.extend(films);
        layers.push(half_space("glass", LayerType::Substrate));
        let groups = HashMap::from([
            ("TiO2".to_string(), Group::new("TiO2")),
            ("SiO2".to_string(), Group::new("SiO2")),
        ]);
        let materials = HashMap::from([
            ("air".to_string(), spec("Konstant", json!({"n": 1.0}))),
            (
                "glass".to_string(),
                spec("Cauchy", json!({"A": 1.50, "B": 0.004, "C": 0.0})),
            ),
            (
                "TiO2".to_string(),
                spec("Cauchy", json!({"A": 2.2, "B": 0.02, "C": 0.001})),
            ),
            (
                "SiO2".to_string(),
                spec(
                    "Sellmeier",
                    json!({"B1": 0.696, "C1": 0.0047, "B2": 0.408, "C2": 0.0135, "B3": 0.0, "C3": 97.9}),
                ),
            ),
            (
                "metal".to_string(),
                spec("Konstant", json!({"n": 0.2, "k": 3.0})),
            ),
            (
                "free".to_string(),
                spec(
                    "Drude",
                    json!({"omega_p": 9.0, "gamma": 0.07, "epsilon_inf": 1.2}),
                ),
            ),
            (
                "knotty".to_string(),
                spec(
                    "Spline",
                    json!({"wavelengths": [380.0, 500.0, 650.0, 800.0],
                           "n": [2.4, 2.25, 2.2, 2.15], "k": [0.05, 0.01, 0.002, 0.001]}),
                ),
            ),
            (
                "kkfilm".to_string(),
                spec(
                    "BSplineEps",
                    json!({"energies": [1.0, 2.5, 3.5, 4.5, 6.0, 8.0],
                           "eps2": [0.0, 0.3, 1.5, 3.0, 1.0, 0.0], "epsilon_inf": 2.5}),
                ),
            ),
        ]);
        let mut m =
            FitModel::from_structure(&Structure::new(layers, groups), "main", materials).unwrap();
        m.add_environment(
            "bare",
            vec![
                Segment::Fixed(half_space("air", LayerType::Ambient)),
                Segment::Fixed(half_space("glass", LayerType::Substrate)),
            ],
        )
        .unwrap();
        m
    }

    fn curves(spec: &[(&str, &str, f64)]) -> Vec<SpectrumData> {
        spec.iter()
            .enumerate()
            .map(|(i, &(env, ch, deg))| {
                let wl = wavelengths();
                // Values off the model so residuals are not all zero.
                let v = wl.iter().map(|w| 0.1 + 1e-4 * w).collect();
                SpectrumData::new(
                    format!("d{i}"),
                    ch,
                    Radians::from_degrees(deg).unwrap(),
                    wl,
                    v,
                    0.01,
                )
                .on(env)
            })
            .collect()
    }

    fn problem(data: &[(&str, &str, f64)], addrs: &[&str]) -> FitProblem {
        FitProblem::new(
            model(),
            curves(data),
            addrs.iter().map(|a| Parameter::parse(a).unwrap()).collect(),
        )
    }

    /// The source's Jacobian and central differences at a point off the
    /// start: `(analytic, fd, m, n)`.
    fn both(p: &FitProblem) -> (Vec<f64>, Vec<f64>, usize, usize) {
        let c = Compiled::new(p).unwrap();
        let n = c.free.len();
        let u: Vec<f64> = (0..n).map(|j| 0.03 * (j as f64 + 1.0)).collect();
        let src = FitJacobian::new(&c, None).unwrap();
        let mut a = Vec::new();
        let m = src.fill(&u, &mut a).unwrap().expect("the source declined");
        let mut f = Vec::new();
        let r = |u: &[f64], out: &mut Vec<f64>| c.residuals(u, out);
        build_jacobian(&r, &u, &mut f).unwrap();
        assert_eq!(a.len(), f.len());
        assert_eq!(m * n, f.len());
        (a, f, m, n)
    }

    /// Every column within `rel` of its largest finite-difference entry.
    fn close(p: &FitProblem, rel: f64) {
        let (a, f, m, n) = both(p);
        for j in 0..n {
            let scale = (0..m).map(|i| f[i * n + j].abs()).fold(0.0, f64::max);
            assert!(scale > 0.0, "column {j} is zero");
            for i in 0..m {
                let (x, y) = (a[i * n + j], f[i * n + j]);
                assert!(
                    (x - y).abs() <= rel * scale,
                    "column {j} row {i}: {x} vs {y} (column scale {scale})"
                );
            }
        }
    }

    #[test]
    fn differenced_columns_are_the_old_bits() {
        // Index columns on two environments: every column is differenced,
        // reached blocks with the old arithmetic, the rest exactly zero.
        let p = problem(
            &[
                ("default", "R_avg", 8.0),
                ("bare", "Rs", 45.0),
                ("default", "Psi_R", 60.0),
            ],
            &[
                "material[TiO2].A",
                "material[glass].A",
                "group[SiO2].n_factor",
            ],
        );
        let (a, f, _, _) = both(&p);
        assert!(a.iter().zip(&f).all(|(x, y)| x.to_bits() == y.to_bits()));
    }

    #[test]
    fn thickness_columns_match_differences_in_every_coherence_mode() {
        for mode in [0, 1, 2] {
            let mut p = problem(
                &[
                    ("default", "R_avg", 8.0),
                    ("default", "Tp", 45.0),
                    ("default", "A_avg", 0.0),
                    ("bare", "Rp", 60.0),
                ],
                &[
                    "layer[main:1].thickness",
                    "layer[main:2].thickness",
                    "group[TiO2].thick_factor",
                    "group[SiO2].thick_summand",
                    "material[TiO2].B",
                ],
            );
            p.coherence_mode = mode;
            // A summand started at 0 takes a 1e-3 scale: central
            // differences then step 6e-9 nm and read rounding noise.
            p.parameters[3] = p.parameters[3].clone().start(2.0);
            close(&p, 1e-5);
            // The chain served the thickness columns (differenced ones
            // would be the old bits) and differenced the index column.
            let (a, f, m, n) = both(&p);
            let same = |j: usize| (0..m).all(|i| a[i * n + j].to_bits() == f[i * n + j].to_bits());
            assert!(
                (0..4).all(|j| !same(j)),
                "mode {mode}: a thickness column was differenced"
            );
            assert!(same(4), "mode {mode}");
        }
    }

    #[test]
    fn rough_interfaces_keep_the_thickness_chain_right() {
        let mut p = problem(
            &[("default", "R_avg", 8.0), ("default", "Ts", 50.0)],
            &[
                "layer[main:1].thickness",
                "layer[main:2].thickness",
                "layer[main:3].thickness",
            ],
        );
        p.model = rough_model();
        close(&p, 1e-5);
        let (a, f, m, n) = both(&p);
        assert!((0..n).all(|j| (0..m).any(|i| a[i * n + j].to_bits() != f[i * n + j].to_bits())));
    }

    #[test]
    fn a_tied_thickness_and_a_mixed_read_stay_right() {
        // Psi has no analytic chain: that environment is differenced.
        let mut p = problem(
            &[("default", "Rs", 20.0), ("default", "Psi_R", 65.0)],
            &["layer[main:1].thickness"],
        );
        p.ties = vec![
            Tie::parse("layer[main:1].thickness")
                .unwrap()
                .follow("layer[main:3].thickness", 0.75, 0.0)
                .unwrap(),
        ];
        close(&p, 1e-5);
        let mut p = problem(&[("default", "T_avg", 0.0)], &["layer[main:1].thickness"]);
        p.ties = vec![
            Tie::parse("layer[main:1].thickness")
                .unwrap()
                .equal("layer[main:3].thickness")
                .unwrap(),
        ];
        close(&p, 1e-5);
    }

    #[test]
    fn instrument_columns_and_the_datasets_they_move_stay_right() {
        // d0 is read through a gain, d1 is tilted and blurred (solved on
        // its own), d2 is plain and keeps the thickness chain.
        let mut p = problem(
            &[
                ("default", "R_avg", 8.0),
                ("default", "Tp", 45.0),
                ("default", "Rs", 20.0),
            ],
            &[
                "layer[main:1].thickness",
                "data[d0].scale",
                "data[d0].offset",
                "data[d1].angle_offset",
                "data[d1].bandwidth",
                "data[d1].wavelength_shift",
            ],
        );
        p.data[0].instrument.scale = 0.9;
        p.data[0].instrument.offset = 0.01;
        p.data[1].instrument.bandwidth = 4.0;
        close(&p, 1e-5);
    }

    #[test]
    fn free_form_knots_and_their_smoothness_rows_stay_right() {
        let mut p = problem(
            &[("default", "R_avg", 8.0), ("default", "Tp", 50.0)],
            &[
                "material[knotty].n[1]",
                "material[knotty].k[0]",
                "material[kkfilm].eps2[2]",
                "layer[main:1].thickness",
            ],
        );
        p.model = model_with(vec![
            Layer::film(80.0, "knotty"),
            Layer::film(120.0, "SiO2"),
            Layer::film(60.0, "kkfilm"),
        ]);
        p.smoothness = vec![
            crate::fit::Smoothness::new("knotty", "n", 0.01),
            crate::fit::Smoothness::new("knotty", "k", 0.01),
            crate::fit::Smoothness::new("kkfilm", "eps2", 0.1),
        ];
        close(&p, 1e-5);
        // Two rows for each list of four, four for the six coefficients.
        let c = Compiled::new(&p).unwrap();
        assert_eq!(c.n_smooth(), 2 + 2 + 4);
    }

    #[test]
    fn nk_columns_in_closed_form_match_differences() {
        let wl = wavelengths();
        let col = |q| NkColumn::new(q, vec![1.0; wl.len()], 1e-3);
        let mut p = problem(
            &[("default", "R_avg", 8.0)],
            &[
                "material[TiO2].A",
                "material[TiO2].B",
                "material[TiO2].C",
                "material[SiO2].B1",
                "material[SiO2].C2",
                "material[metal].n",
                "material[metal].k",
                "layer[main:2].thickness",
                "material[free].omega_p",
                "material[free].gamma",
                "material[free].epsilon_inf",
            ],
        );
        p.nk_data = vec![
            NkData::new(
                "t",
                "TiO2",
                wl.clone(),
                vec![col(NkQuantity::N), col(NkQuantity::K)],
            ),
            NkData::new("s", "SiO2", wl.clone(), vec![col(NkQuantity::N)]),
            NkData::new(
                "m",
                "metal",
                wl.clone(),
                vec![
                    col(NkQuantity::Eps1),
                    col(NkQuantity::Eps2),
                    col(NkQuantity::K),
                ],
            ),
            NkData::new(
                "f",
                "free",
                wl.clone(),
                vec![col(NkQuantity::N), col(NkQuantity::Eps2)],
            ),
        ];
        // eps_inf leaves eps2 alone; at its default scale the central
        // difference reads eps ~ -30 rounded through sqrt and back.
        p.parameters[10] = p.parameters[10].clone().scale(0.5);
        close(&p, 1e-5);
        // The closed forms served them: differenced, they would be the
        // old bits.
        let (a, f, m, n) = both(&p);
        let same = |j: usize| (0..m).all(|i| a[i * n + j].to_bits() == f[i * n + j].to_bits());
        for j in [1, 2, 3, 8, 9] {
            assert!(!same(j), "column {j} was differenced");
        }
    }

    #[test]
    fn priors_and_log_parameters() {
        let mut p = problem(
            &[("default", "R_avg", 8.0)],
            &["layer[main:2].thickness", "group[TiO2].thick_factor"],
        );
        p.parameters[0] = p.parameters[0].clone().prior(118.0, 3.0);
        p.parameters[1] = p.parameters[1].clone().log();
        close(&p, 1e-5);
    }
}
