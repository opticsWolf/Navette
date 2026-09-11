//! smatrix::solver — configured thin-film solver over `core_engine`.
//!
//! Rust-first port of the Python `ScatterMatrix` driver: input validation,
//! the parallel per-point solve, per-point derives (Stokes, Psi/Delta,
//! phases), the cross-wavelength dispersion post-pass, and keyed output.
//! The PyO3 `core_engine` binding thins onto this; Rust consumers solve
//! with no Python.

use std::f64::consts::PI;
use std::sync::OnceLock;

use num_complex::{Complex64, ComplexFloat};
use rayon::prelude::*;

use super::core_engine::*;
use super::needle_engine::*;
use super::optimizer::{char_func, char_func_xy, n_eff_bound};
use super::optics_core::{nevot_croce_factors, redheffer_product_complex_field_inner, w_function_inner};
use super::needle_operator::*;
use super::optics_core::C_NM_PER_FS;

// ---------------------------------------------------------------------------
// Solver
// ---------------------------------------------------------------------------

/// Above this many index-cache entries (`n_layers * n_wavs`) the reciprocal
/// cache is built on the rayon pool. Below it the hand-off costs more than the
/// divisions; the map is elementwise either way, so the two agree bit for bit.
const RECIP_PAR_THRESHOLD: usize = 4096;

/// Configured solver: validated inputs + precomputed index caches.
/// `indices_layer_major` is `(n_layers, n_wavs)` row-major complex.
pub struct Solver {
  wavls: Vec<f64>,
  sin_theta: Vec<f64>,
  /// Layer indices per wavelength, wav-major and **flat**: wavelength `w`'s
  /// layers are `n_cache[w * n_layers..][..n_layers]`.
  ///
  /// Flat rather than `Vec<Vec<Complex64>>`, which is what it was. The nested
  /// form is one heap allocation per wavelength -- 40 000 of them for a
  /// 20 000-point grid counting the inverse cache -- rebuilt on every
  /// `core_engine` call, and that alone was about 45% of what the call spent
  /// end to end (R5.3).
  n_cache: Vec<Complex64>,
  /// Reciprocals of [`Solver::n_cache`], same layout.
  inv_n_cache: Vec<Complex64>,
  /// Wav-major interleaved `[re, im]` view of `n_cache`, built on first use.
  ///
  /// Only `needle_gradient` reads it. Building it eagerly charged every
  /// `solve` for a layout nothing on that path looks at.
  flat_cache: OnceLock<Vec<f64>>,
  thicknesses: Vec<f64>,
  incoherent_flags: Vec<i32>,
  rough_types: Vec<i32>,
  rough_vals: Vec<f64>,
  coherence_mode: i32,
  n_layers: usize,
}

/// Request masks for the `ScatterMatrix` convenience views (mirror the
/// Python view methods exactly; `pol` is 's', 'p' or 'u').
pub fn rt_request(pol: &str) -> Result<u64, String> {
  let mut req = REQ_R_AVG | REQ_T_AVG;
  match pol {
    "s" => req |= REQ_RS | REQ_TS,
    "p" => req |= REQ_RP | REQ_TP,
    "u" => req |= REQ_RS | REQ_TS | REQ_RP | REQ_TP,
    _ => return Err("pol must be 's', 'p', or 'u'".to_string()),
  }
  Ok(req)
}

pub fn ellipsometry_request(transmission: bool) -> u64 {
  let mut req = REQ_PSI_R | REQ_DELTA_R | REQ_DOP_R | REQ_RS | REQ_RP | REQ_R_AVG;
  if transmission {
    req |= REQ_PSI_T | REQ_DELTA_T | REQ_DOP_T | REQ_TS | REQ_TP | REQ_T_AVG;
  }
  req
}

pub fn absorption_request() -> u64 {
  REQ_A_S | REQ_A_P | REQ_A_AVG
}

pub fn amplitudes_request() -> u64 {
  REQ_RS_C | REQ_RP_C | REQ_TS_C | REQ_TP_C
}

pub fn stokes_request(reflection: bool, transmission: bool) -> Result<u64, String> {
  let mut req = 0;
  if reflection {
    req |= REQ_S0_R | REQ_S1_R | REQ_S2_R | REQ_S3_R;
  }
  if transmission {
    req |= REQ_S0_T | REQ_S1_T | REQ_S2_T | REQ_S3_T;
  }
  if req == 0 {
    return Err("select reflection and/or transmission".to_string());
  }
  Ok(req)
}

pub fn dispersion_request(
  reflection: bool,
  transmission: bool,
  s_pol: bool,
  p_pol: bool,
) -> Result<u64, String> {
  let mut req = 0;
  if reflection && s_pol {
    req |= REQ_DISP_R_S;
  }
  if reflection && p_pol {
    req |= REQ_DISP_R_P;
  }
  if transmission && s_pol {
    req |= REQ_DISP_T_S;
  }
  if transmission && p_pol {
    req |= REQ_DISP_T_P;
  }
  if req == 0 {
    return Err("select at least one channel for dispersion".to_string());
  }
  Ok(req)
}

/// `max(|1-Rs-Ts|, |1-Rp-Tp|)` per grid point (energy-conservation residual).
/// Slices must share length; empty input is refused.
pub fn energy_conservation(
  rs: &[f64],
  rp: &[f64],
  ts: &[f64],
  tp: &[f64],
) -> Result<Vec<f64>, String> {
  if rs.len() != rp.len() || rs.len() != ts.len() || rs.len() != tp.len() {
    return Err("energy_conservation: R/T slices must share length".to_string());
  }
  if rs.is_empty() {
    return Err("energy_conservation: empty input".to_string());
  }
  Ok(rs.iter().zip(rp).zip(ts).zip(tp)
    .map(|(((s, p), t_s), t_p)| ((1.0 - s - t_s).abs()).max((1.0 - p - t_p).abs()))
    .collect())
}

impl Solver {
  #[allow(clippy::too_many_arguments)]
  /// Shape and range checks shared by both constructors. `n_indices` is the
  /// index cache's element count, whichever layout it arrived in.
  #[allow(clippy::too_many_arguments)]
  fn validate(
    wavelengths: &[f64],
    sin_theta: &[f64],
    n_layers: usize,
    n_indices: usize,
    thicknesses: &[f64],
    incoherent_flags: &[i32],
    rough_types: &[i32],
    rough_vals: &[f64],
    coherence_mode: i32,
  ) -> Result<(), String> {
    if wavelengths.is_empty() {
      return Err("wavelengths must be non-empty".to_string());
    }
    if sin_theta.is_empty() {
      return Err("angles must be non-empty".to_string());
    }
    if n_layers < 2 {
      return Err("need at least 2 layers (ambient + substrate)".to_string());
    }
    if n_indices != n_layers * wavelengths.len() {
      return Err(format!(
        "indices length {} must equal n_layers ({}) × n_wavs ({})",
        n_indices,
        n_layers,
        wavelengths.len()
      ));
    }
    for (name, v) in [
      ("thicknesses", thicknesses.len()),
      ("incoherent_flags", incoherent_flags.len()),
      ("roughness_types", rough_types.len()),
      ("roughness_values", rough_vals.len()),
    ] {
      if v != n_layers {
        return Err(format!("{name} length {v} must equal n_layers {n_layers}"));
      }
    }
    if !(MODE_A..=MODE_C).contains(&coherence_mode) {
      return Err(
        "coherence_mode must be 0 (front_block), 1 (coherency_matrix), or 2 (fully_coherent)."
          .to_string(),
      );
    }
    Ok(())
  }

  pub fn new(
    wavelengths: &[f64],
    sin_theta: &[f64],
    indices_layer_major: &[Complex64],
    n_layers: usize,
    thicknesses: &[f64],
    incoherent_flags: &[i32],
    rough_types: &[i32],
    rough_vals: &[f64],
    coherence_mode: i32,
  ) -> Result<Self, String> {
    Self::validate(
      wavelengths,
      sin_theta,
      n_layers,
      indices_layer_major.len(),
      thicknesses,
      incoherent_flags,
      rough_types,
      rough_vals,
      coherence_mode,
    )?;
    let n_wavs = wavelengths.len();
    let mut n_cache = Vec::with_capacity(n_layers * n_wavs);
    for w in 0..n_wavs {
      for l in 0..n_layers {
        n_cache.push(indices_layer_major[l * n_wavs + w]);
      }
    }
    Ok(Self::assemble(
      wavelengths,
      sin_theta,
      n_cache,
      n_layers,
      thicknesses,
      incoherent_flags,
      rough_types,
      rough_vals,
      coherence_mode,
    ))
  }

  /// Same solver from the wav-major interleaved `[re, im]` cache the Python
  /// boundary already holds: `n_stack_flat[(w * n_layers + l) * 2]` and `+ 1`.
  ///
  /// `new` wants `(n_layers, n_wavs)` complex, so reaching it from that layout
  /// meant transposing into layer-major and letting `new` transpose straight
  /// back -- two passes over the whole cache to arrive where the caller
  /// started (R5.3).
  #[allow(clippy::too_many_arguments)]
  pub fn from_wav_major_flat(
    wavelengths: &[f64],
    sin_theta: &[f64],
    n_stack_flat: &[f64],
    n_layers: usize,
    thicknesses: &[f64],
    incoherent_flags: &[i32],
    rough_types: &[i32],
    rough_vals: &[f64],
    coherence_mode: i32,
  ) -> Result<Self, String> {
    let n_wavs = wavelengths.len();
    if n_stack_flat.len() != n_layers * n_wavs * 2 {
      return Err(format!(
        "index cache length {} must equal n_layers ({n_layers}) x n_wavs ({n_wavs}) x 2",
        n_stack_flat.len()
      ));
    }
    let n_cache: Vec<Complex64> = n_stack_flat
      .chunks_exact(2)
      .map(|c| Complex64::new(c[0], c[1]))
      .collect();
    Self::validate(
      wavelengths,
      sin_theta,
      n_layers,
      n_cache.len(),
      thicknesses,
      incoherent_flags,
      rough_types,
      rough_vals,
      coherence_mode,
    )?;
    Ok(Self::assemble(
      wavelengths,
      sin_theta,
      n_cache,
      n_layers,
      thicknesses,
      incoherent_flags,
      rough_types,
      rough_vals,
      coherence_mode,
    ))
  }

  /// The shared tail of both constructors: take an already wav-major flat
  /// index cache, derive the reciprocals, and own the rest.
  #[allow(clippy::too_many_arguments)]
  fn assemble(
    wavelengths: &[f64],
    sin_theta: &[f64],
    n_cache: Vec<Complex64>,
    n_layers: usize,
    thicknesses: &[f64],
    incoherent_flags: &[i32],
    rough_types: &[i32],
    rough_vals: &[f64],
    coherence_mode: i32,
  ) -> Self {
    // Elementwise and pure, so serial and parallel are bit-identical; above the
    // threshold the 120 000 complex divisions a 20 000-point grid needs are a
    // third of what building the solver costs.
    let inv_n_cache: Vec<Complex64> = if n_cache.len() >= RECIP_PAR_THRESHOLD {
      n_cache.par_iter().map(|n| n.recip()).collect()
    } else {
      n_cache.iter().map(|n| n.recip()).collect()
    };
    Self {
      wavls: wavelengths.to_vec(),
      sin_theta: sin_theta.to_vec(),
      n_cache,
      inv_n_cache,
      flat_cache: OnceLock::new(),
      thicknesses: thicknesses.to_vec(),
      incoherent_flags: incoherent_flags.to_vec(),
      rough_types: rough_types.to_vec(),
      rough_vals: rough_vals.to_vec(),
      coherence_mode,
      n_layers,
    }
  }

  /// Layer indices for wavelength `w`.
  #[inline]
  fn layer_n(&self, w: usize) -> &[Complex64] {
    &self.n_cache[w * self.n_layers..(w + 1) * self.n_layers]
  }

  /// Reciprocal layer indices for wavelength `w`.
  #[inline]
  fn layer_inv_n(&self, w: usize) -> &[Complex64] {
    &self.inv_n_cache[w * self.n_layers..(w + 1) * self.n_layers]
  }

  /// Wav-major interleaved `[re, im]` cache, built on first use.
  fn flat_cache(&self) -> &[f64] {
    self.flat_cache.get_or_init(|| {
      let mut v = Vec::with_capacity(self.n_cache.len() * 2);
      for nv in &self.n_cache {
        v.push(nv.re);
        v.push(nv.im);
      }
      v
    })
  }

  /// Raw-input constructor mirroring the Python driver: angles in degrees
  /// unless `radians` is set; `indices` is either per-layer (len ==
  /// `n_layers`, broadcast over wavelengths) or full `(n_layers, n_wavs)`
  /// row-major; per-layer option slices default when `None`.
  #[allow(clippy::too_many_arguments)]
  pub fn from_raw(
    wavelengths: &[f64],
    angles: &[f64],
    radians: bool,
    indices: &[Complex64],
    n_layers: usize,
    thicknesses: Option<&[f64]>,
    incoherent_flags: Option<&[i32]>,
    roughness_types: Option<&[i32]>,
    roughness_values: Option<&[f64]>,
    coherence_mode: i32,
  ) -> Result<Self, String> {
    if wavelengths.is_empty() {
      return Err("`wavelengths` must be non-empty.".to_string());
    }
    if angles.is_empty() {
      return Err("`angles` must be non-empty.".to_string());
    }
    let nw = wavelengths.len();
    let full = if indices.len() == n_layers {
      let mut v = Vec::with_capacity(n_layers * nw);
      for &n in indices.iter().take(n_layers) {
        for _ in 0..nw {
          v.push(n);
        }
      }
      v
    } else if indices.len() == n_layers * nw {
      indices.to_vec()
    } else {
      return Err("`indices` must be per-layer or (n_layers, n_wavs).".to_string());
    };
    let theta: Vec<f64> = angles
      .iter()
      .map(|a| if radians { *a } else { a.to_radians() })
      .collect();
    let sin_theta: Vec<f64> = theta.iter().map(|t| t.sin()).collect();
    let zero_f = vec![0.0; n_layers];
    let zero_i = vec![0; n_layers];
    Self::new(
      wavelengths,
      &sin_theta,
      &full,
      n_layers,
      thicknesses.unwrap_or(&zero_f),
      incoherent_flags.unwrap_or(&zero_i),
      roughness_types.unwrap_or(&zero_i),
      roughness_values.unwrap_or(&zero_f),
      coherence_mode,
    )
  }

  pub fn n_angles(&self) -> usize {
    self.sin_theta.len()
  }

  pub fn n_wavs(&self) -> usize {
    self.wavls.len()
  }

  /// Needle-operator gradients over the solver's own grid and stack.
  /// Per-point target/weight inputs are scalars (broadcast) or full
  /// angle-major vectors; `None` means target 0 / weight 1.
  ///
  /// `grads_r`/`grads_t` are the Option-B color buckets (`dF/dcurve` per
  /// point, from `build_needle_targets`), not a target/weight pair; they
  /// default to 0 and are deposited into `P` / `P_T`.
  #[allow(clippy::too_many_arguments)]
  pub fn needle_gradient(
    &self,
    needle_n_per_wav: &[Complex64],
    z_grid: &[f64],
    requested: u64,
    incoherent_flags: Option<&[i32]>,
    targets_r: Option<&[f64]>,
    weights_r: Option<&[f64]>,
    targets_t: Option<&[f64]>,
    weights_t: Option<&[f64]>,
    targets_a: Option<&[f64]>,
    weights_a: Option<&[f64]>,
    targets_phi: Option<&[f64]>,
    weights_phi: Option<&[f64]>,
    targets_tb: Option<&[f64]>,
    weights_tb: Option<&[f64]>,
    targets_rb: Option<&[f64]>,
    weights_rb: Option<&[f64]>,
    targets_ab: Option<&[f64]>,
    weights_ab: Option<&[f64]>,
    grads_r: Option<&[f64]>,
    grads_t: Option<&[f64]>,
    start_idx: usize,
    end_idx: Option<usize>,
    channel: usize,
    calc_s: bool,
    calc_p: bool,
    host_mask: Option<&[bool]>,
    gain_shift_phi: f64,
  ) -> Result<NeedleSolution, String> {
    needle_gradient(
      &self.wavls,
      &self.sin_theta,
      self.n_layers,
      self.flat_cache(),
      &self.thicknesses,
      &self.rough_types,
      &self.rough_vals,
      needle_n_per_wav,
      z_grid,
      requested,
      incoherent_flags.or(Some(&self.incoherent_flags)),
      targets_r,
      weights_r,
      targets_t,
      weights_t,
      targets_a,
      weights_a,
      targets_phi,
      weights_phi,
      targets_tb,
      weights_tb,
      targets_rb,
      weights_rb,
      targets_ab,
      weights_ab,
      grads_r,
      grads_t,
      start_idx,
      end_idx,
      channel,
      calc_s,
      calc_p,
      host_mask,
      gain_shift_phi,
    )
  }

  /// Solve every (angle, wavelength) point for `requested` and derive the
  /// keyed output maps (flat `[n_angles × n_wavs]` angle-major).
  pub fn solve(&self, requested: u64) -> Result<Solution, String> {
    if requested == 0 {
      return Err("empty request mask: select at least one observable".to_string());
    }
    let num_wavs = self.wavls.len();
    let num_angles = self.sin_theta.len();
    let total_points = num_wavs * num_angles;
    let idx_n = self.n_layers - 1;

    let Plan { need_s, need_p, need_cross, level } = resolve_plan(requested);
    let want_phi_rs = requested & (REQ_PHI_RS | REQ_DISP_R_S) != 0;
    let want_phi_rp = requested & (REQ_PHI_RP | REQ_DISP_R_P) != 0;
    let want_phi_ts = requested & (REQ_PHI_TS | REQ_DISP_T_S) != 0;
    let want_phi_tp = requested & (REQ_PHI_TP | REQ_DISP_T_P) != 0;

    let states: Vec<OpticalState> = (0..total_points)
      .into_par_iter()
      .map(|k| {
        let a = k / num_wavs;
        let w = k % num_wavs;
        match level {
          Level::Intensities => solve_point_intensity(
            idx_n,
            self.wavls[w],
            self.sin_theta[a],
            self.layer_n(w),
            self.layer_inv_n(w),
            &self.thicknesses,
            &self.incoherent_flags,
            &self.rough_types,
            &self.rough_vals,
            self.coherence_mode,
            need_s,
            need_p,
          ),
          Level::ComplexAmps | Level::Cross => solve_point(
            idx_n,
            self.wavls[w],
            self.sin_theta[a],
            self.layer_n(w),
            self.layer_inv_n(w),
            &self.thicknesses,
            &self.incoherent_flags,
            &self.rough_types,
            &self.rough_vals,
            self.coherence_mode,
            need_s,
            need_p,
            need_cross,
          ),
        }
      })
      .collect();

    macro_rules! f64buf {
      ($name:ident, $cond:expr) => {
        let mut $name: Option<Vec<f64>> =
          if $cond { Some(vec![0.0; total_points]) } else { None };
      };
    }
    macro_rules! cbuf {
      ($name:ident, $bit:expr) => {
        let mut $name: Option<Vec<Complex64>> = if requested & $bit != 0 {
          Some(vec![Complex64::new(0.0, 0.0); total_points])
        } else {
          None
        };
      };
    }
    f64buf!(b_rs, requested & REQ_RS != 0);
    f64buf!(b_rp, requested & REQ_RP != 0);
    f64buf!(b_ts, requested & REQ_TS != 0);
    f64buf!(b_tp, requested & REQ_TP != 0);
    f64buf!(b_ravg, requested & REQ_R_AVG != 0);
    f64buf!(b_tavg, requested & REQ_T_AVG != 0);
    f64buf!(b_as, requested & REQ_A_S != 0);
    f64buf!(b_ap, requested & REQ_A_P != 0);
    f64buf!(b_aavg, requested & REQ_A_AVG != 0);
    f64buf!(b_psi_r, requested & REQ_PSI_R != 0);
    f64buf!(b_psi_t, requested & REQ_PSI_T != 0);
    f64buf!(b_delta_r, requested & REQ_DELTA_R != 0);
    f64buf!(b_delta_t, requested & REQ_DELTA_T != 0);
    f64buf!(b_dop_r, requested & REQ_DOP_R != 0);
    f64buf!(b_dop_t, requested & REQ_DOP_T != 0);
    f64buf!(b_diatt_r, requested & REQ_DIATT_R != 0);
    f64buf!(b_diatt_t, requested & REQ_DIATT_T != 0);
    f64buf!(b_s0r, requested & REQ_S0_R != 0);
    f64buf!(b_s1r, requested & REQ_S1_R != 0);
    f64buf!(b_s2r, requested & REQ_S2_R != 0);
    f64buf!(b_s3r, requested & REQ_S3_R != 0);
    f64buf!(b_s0t, requested & REQ_S0_T != 0);
    f64buf!(b_s1t, requested & REQ_S1_T != 0);
    f64buf!(b_s2t, requested & REQ_S2_T != 0);
    f64buf!(b_s3t, requested & REQ_S3_T != 0);
    f64buf!(b_retard_r, requested & REQ_RETARD_R != 0);
    f64buf!(b_retard_t, requested & REQ_RETARD_T != 0);
    f64buf!(b_phi_rs, want_phi_rs);
    f64buf!(b_phi_rp, want_phi_rp);
    f64buf!(b_phi_ts, want_phi_ts);
    f64buf!(b_phi_tp, want_phi_tp);
    f64buf!(b_phi_rbs, requested & REQ_PHI_RBS != 0);
    f64buf!(b_phi_rbp, requested & REQ_PHI_RBP != 0);
    f64buf!(b_phi_tbs, requested & REQ_PHI_TBS != 0);
    f64buf!(b_phi_tbp, requested & REQ_PHI_TBP != 0);

    cbuf!(b_rs_c, REQ_RS_C);
    cbuf!(b_rp_c, REQ_RP_C);
    cbuf!(b_ts_c, REQ_TS_C);
    cbuf!(b_tp_c, REQ_TP_C);
    cbuf!(b_rbs_c, REQ_RBS_C);
    cbuf!(b_rbp_c, REQ_RBP_C);
    cbuf!(b_tbs_c, REQ_TBS_C);
    cbuf!(b_tbp_c, REQ_TBP_C);
    cbuf!(b_cross_r, REQ_CROSS_R);
    cbuf!(b_cross_t, REQ_CROSS_T);

    // The derive pass -- per-point algebra plus up to twelve
    // transcendentals per point -- split across the rayon pool by halving the
    // point range and every live buffer together (R5.2). Each destination
    // index is written once, by one task, with the same expression the serial
    // pass used, so the result is bit-identical by construction.
    {
      let sinks = Sinks {
        b_rs: b_rs.as_deref_mut(),
        b_rp: b_rp.as_deref_mut(),
        b_ts: b_ts.as_deref_mut(),
        b_tp: b_tp.as_deref_mut(),
        b_ravg: b_ravg.as_deref_mut(),
        b_tavg: b_tavg.as_deref_mut(),
        b_as: b_as.as_deref_mut(),
        b_ap: b_ap.as_deref_mut(),
        b_aavg: b_aavg.as_deref_mut(),
        b_psi_r: b_psi_r.as_deref_mut(),
        b_psi_t: b_psi_t.as_deref_mut(),
        b_delta_r: b_delta_r.as_deref_mut(),
        b_delta_t: b_delta_t.as_deref_mut(),
        b_dop_r: b_dop_r.as_deref_mut(),
        b_dop_t: b_dop_t.as_deref_mut(),
        b_diatt_r: b_diatt_r.as_deref_mut(),
        b_diatt_t: b_diatt_t.as_deref_mut(),
        b_s0r: b_s0r.as_deref_mut(),
        b_s1r: b_s1r.as_deref_mut(),
        b_s2r: b_s2r.as_deref_mut(),
        b_s3r: b_s3r.as_deref_mut(),
        b_s0t: b_s0t.as_deref_mut(),
        b_s1t: b_s1t.as_deref_mut(),
        b_s2t: b_s2t.as_deref_mut(),
        b_s3t: b_s3t.as_deref_mut(),
        b_retard_r: b_retard_r.as_deref_mut(),
        b_retard_t: b_retard_t.as_deref_mut(),
        b_phi_rs: b_phi_rs.as_deref_mut(),
        b_phi_rp: b_phi_rp.as_deref_mut(),
        b_phi_ts: b_phi_ts.as_deref_mut(),
        b_phi_tp: b_phi_tp.as_deref_mut(),
        b_phi_rbs: b_phi_rbs.as_deref_mut(),
        b_phi_rbp: b_phi_rbp.as_deref_mut(),
        b_phi_tbs: b_phi_tbs.as_deref_mut(),
        b_phi_tbp: b_phi_tbp.as_deref_mut(),
        b_rs_c: b_rs_c.as_deref_mut(),
        b_rp_c: b_rp_c.as_deref_mut(),
        b_ts_c: b_ts_c.as_deref_mut(),
        b_tp_c: b_tp_c.as_deref_mut(),
        b_rbs_c: b_rbs_c.as_deref_mut(),
        b_rbp_c: b_rbp_c.as_deref_mut(),
        b_tbs_c: b_tbs_c.as_deref_mut(),
        b_tbp_c: b_tbp_c.as_deref_mut(),
        b_cross_r: b_cross_r.as_deref_mut(),
        b_cross_t: b_cross_t.as_deref_mut(),
      };
      derive_par(&states, sinks);
    }

    let omega: Vec<f64> =
      self.wavls.iter().map(|&l| 2.0 * PI * C_NM_PER_FS / l).collect();
    let disp = |phi: &Option<Vec<f64>>| -> Option<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
      phi.as_ref().map(|p| dispersion_channel(p, &omega, num_angles, num_wavs))
    };
    let disp_r_s = if requested & REQ_DISP_R_S != 0 { disp(&b_phi_rs) } else { None };
    let disp_r_p = if requested & REQ_DISP_R_P != 0 { disp(&b_phi_rp) } else { None };
    let disp_t_s = if requested & REQ_DISP_T_S != 0 { disp(&b_phi_ts) } else { None };
    let disp_t_p = if requested & REQ_DISP_T_P != 0 { disp(&b_phi_tp) } else { None };

    let mut f64maps: Vec<(String, Vec<f64>)> = Vec::new();
    macro_rules! keep_f64 {
      ($name:expr, $buf:expr) => {
        if let Some(b) = $buf {
          f64maps.push(($name.to_string(), b));
        }
      };
    }
    keep_f64!("Rs", b_rs);
    keep_f64!("Rp", b_rp);
    keep_f64!("Ts", b_ts);
    keep_f64!("Tp", b_tp);
    keep_f64!("R_avg", b_ravg);
    keep_f64!("T_avg", b_tavg);
    keep_f64!("A_s", b_as);
    keep_f64!("A_p", b_ap);
    keep_f64!("A_avg", b_aavg);
    keep_f64!("Psi_R", b_psi_r);
    keep_f64!("Psi_T", b_psi_t);
    keep_f64!("Delta_R", b_delta_r);
    keep_f64!("Delta_T", b_delta_t);
    keep_f64!("DOP_R", b_dop_r);
    keep_f64!("DOP_T", b_dop_t);
    keep_f64!("Diattenuation_R", b_diatt_r);
    keep_f64!("Diattenuation_T", b_diatt_t);
    keep_f64!("S0_R", b_s0r);
    keep_f64!("S1_R", b_s1r);
    keep_f64!("S2_R", b_s2r);
    keep_f64!("S3_R", b_s3r);
    keep_f64!("S0_T", b_s0t);
    keep_f64!("S1_T", b_s1t);
    keep_f64!("S2_T", b_s2t);
    keep_f64!("S3_T", b_s3t);
    keep_f64!("Retardance_R", b_retard_r);
    keep_f64!("Retardance_T", b_retard_t);
    // Phases surface only when explicitly requested (not merely
    // dispersion-needed) — mirrors the binding contract.
    if requested & REQ_PHI_RS != 0 { keep_f64!("phi_rs", b_phi_rs); }
    if requested & REQ_PHI_RP != 0 { keep_f64!("phi_rp", b_phi_rp); }
    if requested & REQ_PHI_TS != 0 { keep_f64!("phi_ts", b_phi_ts); }
    if requested & REQ_PHI_TP != 0 { keep_f64!("phi_tp", b_phi_tp); }
    if requested & REQ_PHI_RBS != 0 { keep_f64!("phi_rbs", b_phi_rbs); }
    if requested & REQ_PHI_RBP != 0 { keep_f64!("phi_rbp", b_phi_rbp); }
    if requested & REQ_PHI_TBS != 0 { keep_f64!("phi_tbs", b_phi_tbs); }
    if requested & REQ_PHI_TBP != 0 { keep_f64!("phi_tbp", b_phi_tbp); }

    let mut c64maps: Vec<(String, Vec<Complex64>)> = Vec::new();
    macro_rules! keep_c {
      ($name:expr, $buf:expr) => {
        if let Some(b) = $buf {
          c64maps.push(($name.to_string(), b));
        }
      };
    }
    keep_c!("rs_c", b_rs_c);
    keep_c!("rp_c", b_rp_c);
    keep_c!("ts_c", b_ts_c);
    keep_c!("tp_c", b_tp_c);
    keep_c!("rbs_c", b_rbs_c);
    keep_c!("rbp_c", b_rbp_c);
    keep_c!("tbs_c", b_tbs_c);
    keep_c!("tbp_c", b_tbp_c);
    keep_c!("cross_R", b_cross_r);
    keep_c!("cross_T", b_cross_t);

    let mut dispmaps: Vec<(String, [Vec<f64>; 4])> = Vec::new();
    // NOTE: dispersion key names are assembled by the caller (binding)
    // from the 4-tuple order [GD, GDD, TOD, FOD]; see `DISP_SUFFIXES`.
    if let Some(t) = disp_r_s { dispmaps.push(("R_s".to_string(), [t.0, t.1, t.2, t.3])); }
    if let Some(t) = disp_r_p { dispmaps.push(("R_p".to_string(), [t.0, t.1, t.2, t.3])); }
    if let Some(t) = disp_t_s { dispmaps.push(("T_s".to_string(), [t.0, t.1, t.2, t.3])); }
    if let Some(t) = disp_t_p { dispmaps.push(("T_p".to_string(), [t.0, t.1, t.2, t.3])); }

    Ok(Solution { n_angles: num_angles, n_wavs: num_wavs, f64maps, c64maps, dispmaps })
  }
}

/// Leaf size of the derive pass's parallel split: a range of at most this many
/// grid points is derived on one thread.
///
/// The split halves the point range and **every live destination buffer**
/// together, so each output index is still written exactly once, by one task,
/// with the same expression it had when the pass was serial. Bit-identity is
/// therefore structural and holds at any leaf size (§13; §0.1 rule 5) -- this
/// constant trades scheduling overhead against parallelism and nothing else.
/// Swept at 256 / 512 / 2048 / 8192 on a rigorous twelve-channel request: the
/// grid sizes that matter land within ~3% of each other, except that a leaf
/// above the grid means no split at all -- 2 000 points cost 0.428 ms at 2048
/// and 0.386 ms at 256/512. 512 is the smallest leaf that still buys something
/// measurable; below ~500 points the whole pass is tens of microseconds and
/// there is nothing left to divide.
const DERIVE_PAR_LEAF: usize = 512;

/// The derive pass's destinations: one optional slice per channel, so a
/// channel nobody requested costs a `None` rather than a buffer.
///
/// Generated from a single list because the halving is the part that must not
/// miss one, and a hand-written `split_at` over 45 fields is exactly the kind
/// of thing that silently does.
macro_rules! derive_sinks {
  (f64: $($f:ident),* $(,)? ; c64: $($c:ident),* $(,)?) => {
    struct Sinks<'a> {
      $( $f: Option<&'a mut [f64]>, )*
      $( $c: Option<&'a mut [Complex64]>, )*
    }

    impl<'a> Sinks<'a> {
      /// Every destination absent -- the starting point for a request that
      /// asks for only a few channels.
      #[cfg(test)]
      fn none() -> Sinks<'a> {
        Sinks { $( $f: None, )* $( $c: None, )* }
      }

      /// Two sink sets over disjoint index ranges: `..mid` and `mid..`.
      fn split_at(self, mid: usize) -> (Sinks<'a>, Sinks<'a>) {
        $(
          let $f = match self.$f {
            Some(b) => { let (lo, hi) = b.split_at_mut(mid); (Some(lo), Some(hi)) }
            None => (None, None),
          };
        )*
        $(
          let $c = match self.$c {
            Some(b) => { let (lo, hi) = b.split_at_mut(mid); (Some(lo), Some(hi)) }
            None => (None, None),
          };
        )*
        (
          Sinks { $( $f: $f.0, )* $( $c: $c.0, )* },
          Sinks { $( $f: $f.1, )* $( $c: $c.1, )* },
        )
      }
    }
  };
}

derive_sinks!(
  f64:
  b_rs, b_rp, b_ts, b_tp,
  b_ravg, b_tavg, b_as, b_ap,
  b_aavg, b_psi_r, b_psi_t, b_delta_r,
  b_delta_t, b_dop_r, b_dop_t, b_diatt_r,
  b_diatt_t, b_s0r, b_s1r, b_s2r,
  b_s3r, b_s0t, b_s1t, b_s2t,
  b_s3t, b_retard_r, b_retard_t, b_phi_rs,
  b_phi_rp, b_phi_ts, b_phi_tp, b_phi_rbs,
  b_phi_rbp, b_phi_tbs, b_phi_tbp,
  ;
  c64:
  b_rs_c, b_rp_c, b_ts_c, b_tp_c,
  b_rbs_c, b_rbp_c, b_tbs_c, b_tbp_c,
  b_cross_r, b_cross_t,
);

/// Derive `states` into `sinks`, halving onto the rayon pool down to
/// `DERIVE_PAR_LEAF`. `states` and `sinks` index the same points.
fn derive_par(states: &[OpticalState], sinks: Sinks<'_>) {
  if states.len() <= DERIVE_PAR_LEAF {
    derive_range(states, sinks);
    return;
  }
  let mid = states.len() / 2;
  let (lo_states, hi_states) = states.split_at(mid);
  let (lo, hi) = sinks.split_at(mid);
  rayon::join(|| derive_par(lo_states, lo), || derive_par(hi_states, hi));
}

/// One thread's share of the derive pass: the per-point algebra and
/// transcendentals, written to whichever destinations exist.
fn derive_range(states: &[OpticalState], mut sinks: Sinks<'_>) {
  macro_rules! put {
    ($buf:ident, $k:expr, $val:expr) => {
      if let Some(b) = sinks.$buf.as_mut() {
        b[$k] = $val;
      }
    };
  }

  for (k, s) in states.iter().enumerate() {
    let rs = s.rs;
    let rp = s.rp;
    let ts = s.ts;
    let tp = s.tp;

    put!(b_rs, k, rs);
    put!(b_rp, k, rp);
    put!(b_ts, k, ts);
    put!(b_tp, k, tp);
    put!(b_ravg, k, 0.5 * (rs + rp));
    put!(b_tavg, k, 0.5 * (ts + tp));
    put!(b_as, k, 1.0 - rs - ts);
    put!(b_ap, k, 1.0 - rp - tp);
    put!(b_aavg, k, 1.0 - 0.5 * (rs + rp) - 0.5 * (ts + tp));

    let s0r = rp + rs;
    let s1r = rp - rs;
    // The `+ 0.0` is not dead weight: `-2.0 * 0.0` is `-0.0`, and Delta is
    // `atan2(s3, s2)`, which answers -pi for a negative zero where it answers
    // +pi for a positive one. An isotropic stack at normal incidence has
    // `cross_r` exactly real, so this is the ordinary case, not a corner --
    // dropping the term flips Delta_R by 2*pi there. IEEE 754 addition maps
    // -0.0 + 0.0 -> +0.0, which is the normalization. The numba reference does
    // the same at the same four places; parity depends on it.
    let s2r = -2.0 * s.cross_r.re + 0.0;
    let s3r = -2.0 * s.cross_r.im + 0.0;
    put!(b_s0r, k, s0r);
    put!(b_s1r, k, s1r);
    put!(b_s2r, k, s2r);
    put!(b_s3r, k, s3r);
    let s0t = tp + ts;
    let s1t = tp - ts;
    // Same negative-zero flush as the reflected pair above.
    let s2t = 2.0 * s.cross_t.re + 0.0;
    let s3t = 2.0 * s.cross_t.im + 0.0;
    put!(b_s0t, k, s0t);
    put!(b_s1t, k, s1t);
    put!(b_s2t, k, s2t);
    put!(b_s3t, k, s3t);

    put!(b_diatt_r, k, s1r / (s0r + 1e-20));
    put!(b_diatt_t, k, s1t / (s0t + 1e-20));

    // Both DOPs are clamped to 1. For transmission the excess is physical
    // bookkeeping -- s2t/s3t are the single-pass cross term while s0t collects
    // the multi-bounce intensity, so their ratio can drift above 1. For
    // reflection the algebra is exact (s1r^2 + s2r^2 + s3r^2 = s0r^2 for a
    // single coherent block) and the excess is pure round-off, measured at
    // 1.0000000000000004. Either way a degree of polarization above 1 is not a
    // number anyone can use, and the asymmetry was a review finding (R6.2).
    put!(b_dop_r, k, ((s1r * s1r + s2r * s2r + s3r * s3r).sqrt() / (s0r + 1e-20)).min(1.0));
    put!(b_dop_t, k, ((s1t * s1t + s2t * s2t + s3t * s3t).sqrt() / (s0t + 1e-20)).min(1.0));

    put!(b_psi_r, k, if rs < RS_FLOOR { PI / 2.0 } else { (rp / rs).sqrt().atan() });
    put!(b_delta_r, k, if rs < RS_FLOOR { 0.0 } else { s3r.atan2(s2r) });
    put!(b_psi_t, k, if ts < TS_FLOOR { PI / 2.0 } else { (tp / ts).sqrt().atan() });
    put!(b_delta_t, k, if ts < TS_FLOOR { 0.0 } else { s3t.atan2(s2t) });

    put!(b_retard_r, k, s.cross_r.arg());
    put!(b_retard_t, k, s.cross_t.arg());

    put!(b_phi_rs, k, s.rs_c.arg());
    put!(b_phi_rp, k, s.rp_c.arg());
    put!(b_phi_ts, k, s.ts_c.arg());
    put!(b_phi_tp, k, s.tp_c.arg());
    put!(b_phi_rbs, k, s.rbs_c.arg());
    put!(b_phi_rbp, k, s.rbp_c.arg());
    put!(b_phi_tbs, k, s.tbs_c.arg());
    put!(b_phi_tbp, k, s.tbp_c.arg());

    put!(b_rs_c, k, s.rs_c);
    put!(b_rp_c, k, s.rp_c);
    put!(b_ts_c, k, s.ts_c);
    put!(b_tp_c, k, s.tp_c);
    put!(b_rbs_c, k, s.rbs_c);
    put!(b_rbp_c, k, s.rbp_c);
    put!(b_tbs_c, k, s.tbs_c);
    put!(b_tbp_c, k, s.tbp_c);
    put!(b_cross_r, k, s.cross_r);
    put!(b_cross_t, k, s.cross_t);
  }
}

/// Solved output: flat angle-major buffers plus grid shape.
pub struct Solution {
  pub n_angles: usize,
  pub n_wavs: usize,
  pub f64maps: Vec<(String, Vec<f64>)>,
  pub c64maps: Vec<(String, Vec<Complex64>)>,
  /// `(channel, [GD, GDD, TOD, FOD])`; key suffixes per `DISP_SUFFIXES`.
  pub dispmaps: Vec<(String, [Vec<f64>; 4])>,
}

/// Key suffixes for the dispersion 4-tuple, in order.
pub const DISP_SUFFIXES: [(&str, usize); 4] =
  [("GD", 0), ("GDD", 1), ("TOD", 2), ("FOD", 3)];

/// Solve pre-expanded arrays (native `solve_structure` core): grid check,
/// half-space convention warning, one solve. `indices` is `(n_rows,
/// n_wavs)` row-major complex; returns `(solution, warnings)`.
#[allow(clippy::too_many_arguments)]
pub fn solve_arrays(
  indices: &[Complex64],
  thicknesses: &[f64],
  incoherent: &[bool],
  rough_types: &[i32],
  rough_vals: &[f64],
  wavelengths: &[f64],
  angles: &[f64],
  radians: bool,
  requested: u64,
  coherence_mode: i32,
) -> Result<(Solution, Vec<String>), String> {
  let n_rows = thicknesses.len();
  if indices.len() != n_rows * wavelengths.len() {
    return Err(format!(
      "solve_arrays: provider grid mismatch: solver arrays vs {n_rows} layers x {} wavelengths.",
      wavelengths.len()
    ));
  }
  let mut warnings = Vec::new();
  if n_rows >= 2 && (thicknesses[0] != 0.0 || thicknesses[n_rows - 1] != 0.0) {
    warnings.push(
      "solve_arrays: first/last thickness is not 0 (ambient/substrate convention); \
       the engine treats row 0/last as half-spaces."
        .to_string(),
    );
  }
  let inc: Vec<i32> = incoherent.iter().map(|b| i32::from(*b)).collect();
  let solver = Solver::from_raw(
    wavelengths,
    angles,
    radians,
    indices,
    n_rows,
    Some(thicknesses),
    Some(&inc),
    Some(rough_types),
    Some(rough_vals),
    coherence_mode,
  )?;
  Ok((solver.solve(requested)?, warnings))
}

// ---------------------------------------------------------------------------
// Needle gradients (moved verbatim from the PyO3 binding)
// ---------------------------------------------------------------------------

/// Needle-gradient output: flat `[n_points × n_depths]` buffers plus shape.
pub struct NeedleSolution {
  pub n_points: usize,
  pub n_depths: usize,
  pub maps: Vec<(String, Vec<f64>)>,
}

#[allow(clippy::too_many_arguments)]
pub fn needle_gradient(
    wavls: &[f64],
    sin_theta: &[f64],
    n_layers: usize,
    n_stack_cache: &[f64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    needle_n_per_wav: &[Complex64],
    z_grid: &[f64],
    requested: u64,
    incoherent_flags: Option<&[i32]>,
    targets_r: Option<&[f64]>,
    weights_r: Option<&[f64]>,
    targets_t: Option<&[f64]>,
    weights_t: Option<&[f64]>,
    targets_a: Option<&[f64]>,
    weights_a: Option<&[f64]>,
    targets_phi: Option<&[f64]>,
    weights_phi: Option<&[f64]>,
    targets_tb: Option<&[f64]>,
    weights_tb: Option<&[f64]>,
    targets_rb: Option<&[f64]>,
    weights_rb: Option<&[f64]>,
    targets_ab: Option<&[f64]>,
    weights_ab: Option<&[f64]>,
    grads_r: Option<&[f64]>,
    grads_t: Option<&[f64]>,
    start_idx: usize,
    end_idx: Option<usize>,
    channel: usize,
    calc_s: bool,
    calc_p: bool,
    host_mask: Option<&[bool]>,
    gain_shift_phi: f64,
) -> Result<NeedleSolution, String>
{
    if requested == 0 {
        return Err(String::from("empty request mask"));
    }

    let num_wavs = wavls.len();
    let num_angles = sin_theta.len();
    let total_points = num_wavs * num_angles;
    let nz = z_grid.len();

    let nl = n_layers;
    if !(0..nl).contains(&start_idx) {
        return Err(String::from("start_idx out of range"));
    }
    let idx_end = end_idx.unwrap_or(nl - 1);
    if idx_end < start_idx + 2 || idx_end >= nl {
        return Err(String::from(
            "end_idx must leave at least one host layer inside [start_idx, end_idx]",
        ));
    }
    if num_wavs == 0 || num_angles == 0 || nz == 0 {
        return Err(String::from("empty grid"));
    }
    if needle_n_per_wav.len() != num_wavs {
        return Err(String::from(
            "needle_n_per_wav must have one complex index per wavelength",
        ));
    }
    if n_stack_cache.len() != num_wavs * nl * 2 {
        return Err(String::from("n_stack_cache layout mismatch"));
    }
    let want_p = requested & NREQ_P != 0;
    let want_pmb = requested & NREQ_P_MB != 0;
    let want_pmb_t = requested & NREQ_P_MB_T != 0;
    let want_pmb_a = requested & NREQ_P_MB_A != 0;
    let want_ptb = requested & NREQ_P_TB != 0;
    let want_prb = requested & NREQ_P_RB != 0;
    let want_pab = requested & NREQ_P_AB != 0;
    let want_pmb_tb = requested & NREQ_P_MB_TB != 0;
    let want_pmb_rb = requested & NREQ_P_MB_RB != 0;
    let want_pmb_ab = requested & NREQ_P_MB_AB != 0;
    let want_pt = requested & NREQ_P_T != 0;
    let want_pa = requested & NREQ_P_A != 0;
    let want_pphi = requested & NREQ_P_PHI != 0;
    let want_disp = max_disp_order(requested).is_some();
    if !calc_s && !calc_p {
        return Err(String::from("no polarization branch enabled"));
    }
    if channel > 3 {
        return Err(String::from("channel must be 0..=3"));
    }

    // ---- R3.2: needle depth and index range -------------------------------
    // `z` outside the eligible span used to be a `debug_assert!` only. Release
    // builds located a host anyway -- `locate_depth_in` falls through to the
    // last layer of the range and hands back a xi past its thickness -- and
    // returned a plausible-looking gradient for a needle that is not where the
    // caller put it.
    //
    // The two paths have different spans within the same call: the coherent
    // kernels are confined to the block `[start_idx, idx_end]`, while the
    // multiblock cascade walks every non-ambient layer. Each span is checked
    // only when the request actually reaches that path, so a multiblock-only
    // call is not held to the block's narrower bound -- and a call that uses
    // both must satisfy both.
    let span_of = |lo: usize, hi: usize| -> f64 {
        thicknesses.iter().take(hi.min(thicknesses.len())).skip(lo).sum()
    };
    let mut spans: Vec<(&str, f64)> = Vec::new();
    if want_p || want_pt || want_pa || want_pphi || want_ptb || want_prb
        || want_pab || want_disp
    {
        spans.push(("coherent block", span_of(start_idx + 1, idx_end)));
    }
    if want_pmb || want_pmb_t || want_pmb_a || want_pmb_tb || want_pmb_rb
        || want_pmb_ab
    {
        spans.push(("multiblock cascade", span_of(1, nl - 1)));
    }
    const Z_TOL: f64 = 1e-9;
    for (i, &z) in z_grid.iter().enumerate() {
        if !z.is_finite() {
            return Err(format!("z_grid[{i}] is not finite ({z})"));
        }
        for (what, span) in &spans {
            if z < -Z_TOL || z > span + Z_TOL {
                return Err(format!(
                    "z_grid[{i}] = {z} lies outside the needle-eligible depth span [0, {span}] of the {what}; depth is absolute, measured from the top of the first eligible layer"
                ));
            }
        }
    }
    for (i, n) in needle_n_per_wav.iter().enumerate() {
        if !n.re.is_finite() || !n.im.is_finite() {
            return Err(format!(
                "needle_n_per_wav[{i}] is not finite ({n}); a NaN needle index makes the whole gradient NaN with nothing naming the point"
            ));
        }
    }

    // Optional per-point merit inputs (default: target 0, weight 1).
    // Scalars broadcast; full vectors are angle-major.
    let load_pair = |a: &Option<&[f64]>, name: &str| -> Result<Option<Vec<f64>>, String> {
        match a {
            Some(arr) => {
                if arr.len() == 1 {
                    Ok(Some(vec![arr[0]; total_points]))
                } else if arr.len() == total_points {
                    Ok(Some(arr.to_vec()))
                } else {
                    Err(format!(
                        "{name} must be a scalar or have num_angles*num_wavs entries (angle-major)",
                    ))
                }
            }
            None => Ok(None),
        }
    };
    let tgt = load_pair(&targets_r, "targets_r")?;
    let wgt = load_pair(&weights_r, "weights_r")?;
    let target_of = |k: usize| tgt.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_of = |k: usize| wgt.as_ref().map(|t| t[k]).unwrap_or(1.0);
    let tgt_t = load_pair(&targets_t, "targets_t")?;
    let wgt_t = load_pair(&weights_t, "weights_t")?;
    let tgt_a = load_pair(&targets_a, "targets_a")?;
    let wgt_a = load_pair(&weights_a, "weights_a")?;
    let tgt_phi = load_pair(&targets_phi, "targets_phi")?;
    let wgt_phi = load_pair(&weights_phi, "weights_phi")?;
    let target_t_of = |k: usize| tgt_t.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_t_of = |k: usize| wgt_t.as_ref().map(|t| t[k]).unwrap_or(1.0);
    let target_a_of = |k: usize| tgt_a.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_a_of = |k: usize| wgt_a.as_ref().map(|t| t[k]).unwrap_or(1.0);
    let target_phi_of = |k: usize| tgt_phi.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_phi_of = |k: usize| wgt_phi.as_ref().map(|t| t[k]).unwrap_or(1.0);
    let tgt_tb = load_pair(&targets_tb, "targets_tb")?;
    let wgt_tb = load_pair(&weights_tb, "weights_tb")?;
    let tgt_rb = load_pair(&targets_rb, "targets_rb")?;
    let wgt_rb = load_pair(&weights_rb, "weights_rb")?;
    let tgt_ab = load_pair(&targets_ab, "targets_ab")?;
    let wgt_ab = load_pair(&weights_ab, "weights_ab")?;
    let target_tb_of = |k: usize| tgt_tb.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_tb_of = |k: usize| wgt_tb.as_ref().map(|t| t[k]).unwrap_or(1.0);
    let target_rb_of = |k: usize| tgt_rb.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_rb_of = |k: usize| wgt_rb.as_ref().map(|t| t[k]).unwrap_or(1.0);
    let target_ab_of = |k: usize| tgt_ab.as_ref().map(|t| t[k]).unwrap_or(0.0);
    let weight_ab_of = |k: usize| wgt_ab.as_ref().map(|t| t[k]).unwrap_or(1.0);

    // ---- R4.2: Option-B color gradient buckets ----------------------------
    // `grads_r`/`grads_t` are NOT a (target, weight) pair: each entry is the
    // chain-rule factor g = dF/dcurve for one solver point, already carrying
    // the demand's weight, its current residual and the U-curve half. They
    // ride the R and T channels, exactly as `needle_pass::needle_pass_scan`
    // adds them into its single accumulator, so a caller that folds a color
    // demand in Python gets the same number the native pipeline computes.
    //
    // Default 0.0, not 1.0: absent color demands must deposit nothing.
    let grd_r = load_pair(&grads_r, "grads_r")?;
    let grd_t = load_pair(&grads_t, "grads_t")?;
    let grad_r_of = |k: usize| grd_r.as_ref().map(|g| g[k]).unwrap_or(0.0);
    let grad_t_of = |k: usize| grd_t.as_ref().map(|g| g[k]).unwrap_or(0.0);
    // A color bucket can only be deposited onto a channel that is being
    // computed. Dropping it silently would reproduce the very bug this
    // exists to fix -- a color demand optimized against, with no error and
    // a gradient of zero -- so say so instead. An all-zero array is not a
    // demand and passes: callers hand the fold's arrays through unfiltered.
    let any_nonzero = |v: &Option<Vec<f64>>| {
        v.as_ref().is_some_and(|a| a.iter().any(|&x| x != 0.0))
    };
    if any_nonzero(&grd_r) && !want_p {
        return Err(String::from(
            "grads_r carries non-zero color deposits but NREQ_P was not requested; the R-channel gradient they belong to is not being computed"
        ));
    }
    if any_nonzero(&grd_t) && !want_pt {
        return Err(String::from(
            "grads_t carries non-zero color deposits but NREQ_P_T was not requested; the T-channel gradient they belong to is not being computed"
        ));
    }
    for (name, v) in [("grads_r", &grd_r), ("grads_t", &grd_t)] {
        if let Some(a) = v
            && let Some(i) = a.iter().position(|x| !x.is_finite()) {
            return Err(format!("{name}[{i}] is not finite ({})", a[i]));
        }
    }

    // Incoherent flags only needed for the multiblock path.
    let want_any_pmb =
        want_pmb || want_pmb_t || want_pmb_a || want_pmb_tb || want_pmb_rb || want_pmb_ab;
    let inc = match (&incoherent_flags, want_any_pmb) {
        (_, false) => None,
        (None, true) => {
            return Err(String::from(
                "NREQ_P_MB requires incoherent_flags",
            ))
        }
        (Some(a), true) => {
            let v = a;
            if v.len() != nl {
                return Err(String::from(
                    "incoherent_flags must have n_layers entries",
                ));
            }
            Some(v.to_vec())
        }
    };
    let mask = match &host_mask {
        Some(a) => {
            let v = a;
            if v.len() != nl {
                return Err(String::from(
                    "host_mask must have n_layers entries",
                ));
            }
            Some(v.to_vec())
        }
        None => None,
    };

    // Host maps are geometry-only: compute once, share across all points.
    let mb_locs = match &inc {
        Some(flags) => Some(locate_hosts_multiblock(thicknesses, flags, z_grid, mask.as_deref())
            ?),
        None => None,
    };
    let coh_locs: Vec<(usize, f64)> =
        if want_p || want_pt || want_pa || want_pphi || want_ptb || want_prb || want_pab || want_disp {
        z_grid
            .iter()
            .map(|&z| locate_depth_in(thicknesses, start_idx, idx_end, z))
            .collect()
    } else {
        Vec::new()
    };

    struct PointOut {
        p: [Option<Vec<f64>>; 2],
        pmb: [Option<Vec<f64>>; 2],
        q: [Option<Vec<f64>>; 2], // Q rows (order 0), flattened nz
        pt: [Option<Vec<f64>>; 2],
        pa: [Option<Vec<f64>>; 2],
        pphi: [Option<Vec<f64>>; 2],
        pmb_t: [Option<Vec<f64>>; 2],
        pmb_a: [Option<Vec<f64>>; 2],
        ptb: [Option<Vec<f64>>; 2],
        prb: [Option<Vec<f64>>; 2],
        pab: [Option<Vec<f64>>; 2],
        pmb_tb: [Option<Vec<f64>>; 2],
        pmb_rb: [Option<Vec<f64>>; 2],
        pmb_ab: [Option<Vec<f64>>; 2],
    }
    impl PointOut {
        fn empty() -> Self {
            PointOut {
                p: [None, None], pmb: [None, None], q: [None, None],
                pt: [None, None], pa: [None, None], pphi: [None, None],
                pmb_t: [None, None], pmb_a: [None, None],
                ptb: [None, None], prb: [None, None], pab: [None, None],
                pmb_tb: [None, None], pmb_rb: [None, None], pmb_ab: [None, None],
            }
        }
    }

    let pol_on = [calc_s, calc_p];

    // ── Phase A: everything expressible per point, in parallel ──
    let outs: Vec<PointOut> =
        (0..total_points)
            .into_par_iter()
            .map(|k| {
                let a = k / num_wavs;
                let w = k % num_wavs;
                let lam = wavls[w];
                let sin_t = sin_theta[a];
                let base = w * nl * 2;
                let ns: Vec<Complex64> = (0..nl)
                    .map(|l| Complex64::new(n_stack_cache[base + l * 2], n_stack_cache[base + l * 2 + 1]))
                    .collect();
                let nsin_fi = ns[0] * Complex64::new(sin_t, 0.0);
                let np_c = needle_n_per_wav[w];
                let tgt_k = target_of(k);
                let wgt_k = weight_of(k);

                let mut o = PointOut::empty();

                // Coherent observables share ONE fields build per polarization.
                if want_p || want_pt || want_pa || want_pphi || want_ptb || want_prb || want_pab || want_disp {
                    for (pi, &on) in pol_on.iter().enumerate() {
                        if !on {
                            continue;
                        }
                        let pol = pi as i32;
                        let fields = build_stack_fields_range(
                            start_idx, idx_end, &ns, thicknesses, rough_vals, rough_types,
                            lam, nsin_fi, pol,
                        );
                        if want_p {
                            let mut v = p_coherent_from_fields(
                                &fields, nsin_fi, lam, pol, np_c, tgt_k, wgt_k,
                                thicknesses, start_idx, idx_end, z_grid,
                            );
                            // Color deposit onto the same channel. Zero-skip
                            // is the native pattern (needle_pass.rs:889): a
                            // color-free call pays one float compare a point.
                            let g = grad_r_of(k);
                            if g != 0.0 {
                                for (zi, c) in p_coherent_grad_r_from_fields(
                                    &fields, nsin_fi, lam, pol, np_c, g,
                                    thicknesses, start_idx, idx_end, z_grid,
                                ).into_iter().enumerate() {
                                    v[zi] += c;
                                }
                            }
                            o.p[pi] = Some(v);
                        }
                        if want_pt {
                            let mut v = p_coherent_t_from_fields(
                                &fields, nsin_fi, lam, pol, np_c,
                                target_t_of(k), weight_t_of(k),
                                thicknesses, start_idx, idx_end, z_grid,
                            );
                            let g = grad_t_of(k);
                            if g != 0.0 {
                                for (zi, c) in p_coherent_grad_t_from_fields(
                                    &fields, nsin_fi, lam, pol, np_c, g,
                                    thicknesses, start_idx, idx_end, z_grid,
                                ).into_iter().enumerate() {
                                    v[zi] += c;
                                }
                            }
                            o.pt[pi] = Some(v);
                        }
                        if want_pa {
                            o.pa[pi] = Some(p_coherent_a_from_fields(
                                &fields, nsin_fi, lam, pol, np_c,
                                target_a_of(k), weight_a_of(k),
                                thicknesses, start_idx, idx_end, z_grid,
                            ));
                        }
                        if want_pphi {
                            o.pphi[pi] = Some(p_coherent_phi_from_fields(
                                &fields, nsin_fi, lam, pol, np_c, channel,
                                target_phi_of(k), weight_phi_of(k),
                                thicknesses, start_idx, idx_end, z_grid,
                            ));
                        }
                        if want_ptb {
                            o.ptb[pi] = Some(p_coherent_tb_from_fields(
                                &fields, nsin_fi, lam, pol, np_c,
                                target_tb_of(k), weight_tb_of(k),
                                thicknesses, start_idx, idx_end, z_grid,
                            ));
                        }
                        if want_prb {
                            o.prb[pi] = Some(p_coherent_rb_from_fields(
                                &fields, nsin_fi, lam, pol, np_c,
                                target_rb_of(k), weight_rb_of(k),
                                thicknesses, start_idx, idx_end, z_grid,
                            ));
                        }
                        if want_pab {
                            o.pab[pi] = Some(p_coherent_ab_from_fields(
                                &fields, nsin_fi, lam, pol, np_c,
                                target_ab_of(k), weight_ab_of(k),
                                thicknesses, start_idx, idx_end, z_grid,
                            ));
                        }
                        if want_disp {
                            let m = fields.s_left[idx_end];
                            let amp = [m.0, m.1, m.2, m.3][channel];
                            let r2 = amp.norm_sqr();
                            let mut qv = vec![0.0_f64; nz];
                            if r2 > 1e-20 {
                                for (zi, &(j, xi)) in coh_locs.iter().enumerate() {
                                    // Per-channel slope (channel-0 here would
                                    // mix r-motion into t-phase).
                                    let da = needle_slopes4_ddz(
                                        &fields, nsin_fi, j, xi, np_c, pol, lam)[channel];
                                    qv[zi] = (amp.conj() * da).im / r2;
                                }
                            }
                            o.q[pi] = Some(qv);
                        }
                    }
                }

                if let (Some(flags), Some(locs)) = (&inc, &mb_locs) {
                    for (pi, &on) in pol_on.iter().enumerate() {
                        if !on {
                            continue;
                        }
                        if want_pmb {
                            o.pmb[pi] = Some(p_multiblock_point(
                                lam, sin_t, &ns, thicknesses, flags, rough_vals, rough_types,
                                np_c, PmbQuantity::R, tgt_k, wgt_k, locs, pi as i32,
                            ));
                        }
                        if want_pmb_t {
                            o.pmb_t[pi] = Some(p_multiblock_point(
                                lam, sin_t, &ns, thicknesses, flags, rough_vals, rough_types,
                                np_c, PmbQuantity::T,
                                target_t_of(k), weight_t_of(k), locs, pi as i32,
                            ));
                        }
                        if want_pmb_a {
                            o.pmb_a[pi] = Some(p_multiblock_point(
                                lam, sin_t, &ns, thicknesses, flags, rough_vals, rough_types,
                                np_c, PmbQuantity::A,
                                target_a_of(k), weight_a_of(k), locs, pi as i32,
                            ));
                        }
                        if want_pmb_tb {
                            o.pmb_tb[pi] = Some(p_multiblock_point(
                                lam, sin_t, &ns, thicknesses, flags, rough_vals, rough_types,
                                np_c, PmbQuantity::TB,
                                target_tb_of(k), weight_tb_of(k), locs, pi as i32,
                            ));
                        }
                        if want_pmb_rb {
                            o.pmb_rb[pi] = Some(p_multiblock_point(
                                lam, sin_t, &ns, thicknesses, flags, rough_vals, rough_types,
                                np_c, PmbQuantity::RB,
                                target_rb_of(k), weight_rb_of(k), locs, pi as i32,
                            ));
                        }
                        if want_pmb_ab {
                            o.pmb_ab[pi] = Some(p_multiblock_point(
                                lam, sin_t, &ns, thicknesses, flags, rough_vals, rough_types,
                                np_c, PmbQuantity::AB,
                                target_ab_of(k), weight_ab_of(k), locs, pi as i32,
                            ));
                        }
                    }
                }

                o
            })
            .collect::<Vec<_>>();

    // ── Phase B: spectral differentiation chain (crosses wavelengths) ──
    let max_order = max_disp_order(requested);
    // chains[pol][order][k*nz+zi]
    let disp_chain: Vec<Option<Vec<Vec<Vec<f64>>>>> = match max_order {
        None => vec![None, None],
        Some(mo) => {
            let omega: Vec<f64> =
                wavls.iter().map(|&l| 2.0 * std::f64::consts::PI * C_NM_PER_FS / l).collect();
            pol_on
                .iter()
                .enumerate()
                .map(|(pi, &on)| {
                    if !on || !want_disp {
                        return None;
                    }
                    if outs.iter().any(|o| o.q[pi].is_none()) {
                        return None;
                    }
                    let q0: Vec<Vec<f64>> =
                        outs.iter().map(|o| o.q[pi].clone().unwrap()).collect();
                    let mut chain = vec![q0.clone()];
                    for _ in 0..mo {
                        let prev = chain.last().unwrap();
                        chain.push(spectral_gradient_step(prev, &omega, num_wavs, num_angles, nz));
                    }
                    Some(chain)
                })
                .collect()
        }
    };
    let _ = channel;

    // ── Assemble dict ──
    let mut maps: Vec<(String, Vec<f64>)> = Vec::new();

    macro_rules! emit {
        ($name:expr, $field:ident, $pi:expr) => {{
            let name: String = $name;
            let mut flat: Vec<f64> = Vec::with_capacity(total_points * nz);
            for o in &outs {
                match &o.$field[$pi] {
                    Some(v) => flat.extend_from_slice(v),
                    None => {
                        return Err(String::from(
                            "internal error: missing output buffer",
                        ))
                    }
                }
            }
            if gain_shift_phi != 0.0 && name.starts_with("P_PHI") {
                for v in flat.iter_mut() {
                    *v -= gain_shift_phi;
                }
            }
            maps.push((name, flat));
        }};
    }

    let pol_suffix = |pi: usize| if pi == 0 { "s" } else { "p" };
    if want_p {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_{}", pol_suffix(pi)), p, pi);
            }
        }
    }
    if want_pt {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_T_{}", pol_suffix(pi)), pt, pi);
            }
        }
    }
    if want_pa {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_A_{}", pol_suffix(pi)), pa, pi);
            }
        }
    }
    if want_pphi {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_PHI_{}", pol_suffix(pi)), pphi, pi);
            }
        }
    }
    if want_pmb {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("Pmb_{}", pol_suffix(pi)), pmb, pi);
            }
        }
    }
    if want_pmb_t {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("Pmb_T_{}", pol_suffix(pi)), pmb_t, pi);
            }
        }
    }
    if want_pmb_a {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("Pmb_A_{}", pol_suffix(pi)), pmb_a, pi);
            }
        }
    }
    if want_ptb {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_TB_{}", pol_suffix(pi)), ptb, pi);
            }
        }
    }
    if want_prb {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_RB_{}", pol_suffix(pi)), prb, pi);
            }
        }
    }
    if want_pab {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("P_AB_{}", pol_suffix(pi)), pab, pi);
            }
        }
    }
    if want_pmb_tb {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("Pmb_TB_{}", pol_suffix(pi)), pmb_tb, pi);
            }
        }
    }
    if want_pmb_rb {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("Pmb_RB_{}", pol_suffix(pi)), pmb_rb, pi);
            }
        }
    }
    if want_pmb_ab {
        for (pi, &on) in pol_on.iter().enumerate() {
            if on {
                emit!(format!("Pmb_AB_{}", pol_suffix(pi)), pmb_ab, pi);
            }
        }
    }
    const DISP_KEYS: [&str; 5] = ["dphi", "dgd", "dgdd", "dtod", "dfod"];
    if let Some(mo) = max_order {
        for pi in 0..2 {
            if !pol_on[pi] {
                continue;
            }
            if let Some(chain) = &disp_chain[pi] {
                for order in 0..=mo {
                    let key = format!("{}_{}", DISP_KEYS[order], pol_suffix(pi));
                    let mut flat: Vec<f64> = Vec::with_capacity(total_points * nz);
                    for row in &chain[order] {
                        flat.extend_from_slice(row);
                    }
                    maps.push((key, flat));
                }
            }
        }
    }

    Ok(NeedleSolution { n_points: total_points, n_depths: nz, maps })
}

// ---------------------------------------------------------------------------
// Eigenmode tools (moved verbatim from the PyO3 binding)
// ---------------------------------------------------------------------------

impl Solver {
  /// `(lam, per-layer complex indices)` for one wavelength: explicit
  /// index wins, else nearest grid point, else the single grid point.
  pub fn index_column(
    &self,
    wavelength: Option<f64>,
    wav_index: Option<usize>,
  ) -> Result<(f64, Vec<Complex64>), String> {
    let w = match wav_index {
      Some(i) => {
        if i >= self.wavls.len() {
          return Err(format!("wav_index {i} out of range."));
        }
        i
      }
      None => match wavelength {
        Some(lam) => {
          let mut best = 0;
          for (i, w) in self.wavls.iter().enumerate() {
            if (w - lam).abs() < (self.wavls[best] - lam).abs() {
              best = i;
            }
          }
          best
        }
        None => {
          if self.wavls.len() != 1 {
            return Err(format!(
              "specify wavelength or wav_index (grid has {} wavelengths).",
              self.wavls.len()
            ));
          }
          0
        }
      },
    };
    Ok((self.wavls[w], self.layer_n(w).to_vec()))
  }

  /// Scan `|1/r(n_eff)|^2` over a complex effective-index box.
  /// Returns `(real_vals, imag_vals, flat imag-major values)`.
  pub fn landscape(
    &self,
    real_range: (f64, f64),
    imag_range: (f64, f64),
    points_real: usize,
    points_imag: usize,
    pol: i32,
    wavelength: Option<f64>,
    wav_index: Option<usize>,
  ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), String> {
    let (lam, col) = self.index_column(wavelength, wav_index)?;
    Ok(scan_box(
      &col,
      &self.thicknesses,
      &self.rough_types,
      &self.rough_vals,
      lam,
      pol,
      real_range.0,
      real_range.1,
      imag_range.0,
      imag_range.1,
      points_real,
      points_imag,
    ))
  }

  /// Coarse minima of a landscape flat (imag-major `n_imag × n_real`).
  pub fn local_minima(
    flat: &[f64],
    n_real: usize,
    n_imag: usize,
    real_vals: &[f64],
    imag_vals: &[f64],
    median_factor: f64,
  ) -> Vec<(f64, f64)> {
    find_minima(flat, n_real, n_imag, real_vals, imag_vals, median_factor)
  }

  /// Nelder-Mead refine of one complex eigenmode guess.
  /// Returns `(n_eff, characteristic_value)`.
  pub fn refine_mode(
    &self,
    guess: Complex64,
    pol: i32,
    wavelength: Option<f64>,
    wav_index: Option<usize>,
    step: f64,
    tol: f64,
    max_iter: usize,
  ) -> Result<(Complex64, f64), String> {
    let (lam, col) = self.index_column(wavelength, wav_index)?;
    let (re, im, val) = nelder_refine(
      &col,
      &self.thicknesses,
      &self.rough_types,
      &self.rough_vals,
      lam,
      pol,
      (guess.re, guess.im),
      step,
      tol,
      max_iter,
    );
    Ok((Complex64::new(re, im), val))
  }

  /// Scan, locate coarse minima, optionally refine each.
  pub fn find_eigenmodes(
    &self,
    real_range: (f64, f64),
    imag_range: (f64, f64),
    points_real: usize,
    points_imag: usize,
    median_factor: f64,
    refine: bool,
    pol: i32,
    wavelength: Option<f64>,
    wav_index: Option<usize>,
  ) -> Result<Vec<Complex64>, String> {
    let (real_vals, imag_vals, flat) =
      self.landscape(real_range, imag_range, points_real, points_imag, pol, wavelength, wav_index)?;
    let seeds = Self::local_minima(&flat, points_real, points_imag, &real_vals, &imag_vals, median_factor);
    if !refine {
      return Ok(seeds.into_iter().map(|(re, im)| Complex64::new(re, im)).collect());
    }
    let mut out = Vec::with_capacity(seeds.len());
    for (re, im) in seeds {
      let (n_eff, _) =
        self.refine_mode(Complex64::new(re, im), pol, wavelength, wav_index, 1e-3, 1e-9, 200)?;
      out.push(n_eff);
    }
    Ok(out)
  }

  /// `|E(z)|` profile for one eigenmode: `(z, E, start, end, layer_n)`.
  pub fn field_profile(
    &self,
    n_eff: Complex64,
    pol: i32,
    wavelength: Option<f64>,
    wav_index: Option<usize>,
    points_per_layer: usize,
  ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<Complex64>), String> {
    let (lam, col) = self.index_column(wavelength, wav_index)?;
    field_prof(
      &col,
      &self.thicknesses,
      &self.rough_types,
      &self.rough_vals,
      lam,
      n_eff,
      pol,
      points_per_layer,
    )
  }

} // impl Solver (eigen drivers)

/// Landscape scan over the complex effective-index box.
/// Returns `(real_vals, imag_vals, flat imag-major values)`.
pub fn scan_box(
    n_stack: &[Complex64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    lam: f64,
    pol: i32,
    real_min: f64,
    real_max: f64,
    imag_min: f64,
    imag_max: f64,
    points_real: usize,
    points_imag: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (n_slice, d_slice, rt_slice, rv_slice) = (n_stack, thicknesses, rough_types, rough_vals);

    // Reciprocals computed ONCE per wavelength and shared read-only across all
    // grid points (and all rayon threads). Previously this Vec was allocated
    // inside every char_func call — thousands of heap allocations per scan.
    let inv_n: Vec<Complex64> = n_slice.iter().map(|n| n.recip()).collect();

    let real_vals: Vec<f64> = (0..points_real)
        .map(|i| real_min + (i as f64) * (real_max - real_min) / ((points_real - 1) as f64))
        .collect();
    let imag_vals: Vec<f64> = (0..points_imag)
        .map(|i| imag_min + (i as f64) * (imag_max - imag_min) / ((points_imag - 1) as f64))
        .collect();

    let landscape: Vec<f64> = (0..points_imag * points_real)
            .into_par_iter()
            .map(|idx| {
                let i = idx / points_real;
                let j = idx % points_real;
                let nr = real_vals[j];
                let ni = imag_vals[i];
                let n_eff = Complex64::new(nr, ni);
                char_func(n_slice, &inv_n, d_slice, rt_slice, rv_slice, lam, n_eff, pol)
            })
            .collect();
    (real_vals, imag_vals, landscape)
}

/// Coarse local minima below `median_factor * median(values)`.
/// `flat` is imag-major with shape `(n_imag, n_real)`; the first/last
/// real columns are skipped (matches the reference sweep).
pub fn find_minima(
    flat: &[f64],
    n_real: usize,
    n_imag: usize,
    real_vals: &[f64],
    imag_vals: &[f64],
    median_factor: f64,
) -> Vec<(f64, f64)> {
    let at = |i: usize, j: usize| flat[i * n_real + j];

    // True median of the landscape (the previous code averaged, which the
    // `median_factor` name and the Python reference (`np.median`) do not).
    // Sentinel 1e30 cells sort to the top and so don't perturb the median,
    // whereas they badly skewed the mean.
    let mut sorted: Vec<f64> = flat.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let len = sorted.len();
    let median = if len == 0 {
        0.0
    } else if len % 2 == 1 {
        sorted[len / 2]
    } else {
        0.5 * (sorted[len / 2 - 1] + sorted[len / 2])
    };
    let threshold = median * median_factor;

    let mut candidates = Vec::new();
    if n_real < 2 {
        return candidates;
    }
    for i in 0..n_imag {
        // Skip the first/last real columns, matching the reference
        // (`for j in range(1, len(Nr) - 1)`); all imag rows are scanned so
        // lossless modes on the Im=0 edge are still detected.
        for j in 1..n_real - 1 {
            let val = at(i, j);
            if val >= threshold {
                continue;
            }
            let i0 = i.saturating_sub(1);
            let i1 = (i + 1).min(n_imag - 1);
            let j0 = j.saturating_sub(1);
            let j1 = (j + 1).min(n_real - 1);
            let mut is_min = true;
            'neighbors: for ii in i0..=i1 {
                for jj in j0..=j1 {
                    if ii == i && jj == j {
                        continue;
                    }
                    if at(ii, jj) <= val {
                        is_min = false;
                        break 'neighbors;
                    }
                }
            }
            if is_min {
                candidates.push((real_vals[j], imag_vals[i]));
            }
        }
    }
    candidates
}

/// Nelder-Mead refine of one complex eigenmode guess.
/// Returns `(re, im, characteristic_value)`.
pub fn nelder_refine(
    n_stack: &[Complex64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    lam: f64,
    pol: i32,
    x0: (f64, f64),
    step: f64,
    tol: f64,
    max_iter: usize,
) -> (f64, f64, f64) {
    let (n_slice, d_slice, rt_slice, rv_slice) = (n_stack, thicknesses, rough_types, rough_vals);

    // Reciprocals computed once and reused across every simplex evaluation.
    let inv_n: Vec<Complex64> = n_slice.iter().map(|n| n.recip()).collect();

    // R3.3: project every candidate into the physical box instead of letting
    // the simplex leave it. `char_func_xy` already reports 1e30 outside, but a
    // rejection alone leaves Nelder-Mead reflecting against an infinite wall
    // it cannot see past; projection keeps the search well-defined and simply
    // slides along the boundary. Inert for any seed with a real mode nearby.
    let bound = n_eff_bound(n_slice);
    let clamp = |p: [f64; 2]| -> [f64; 2] {
        [
            if p[0].is_finite() { p[0].clamp(-bound, bound) } else { 0.0 },
            if p[1].is_finite() { p[1].clamp(-bound, bound) } else { 0.0 },
        ]
    };

    let mut simplex = [clamp([x0.0, x0.1]),
        clamp([x0.0 + step, x0.1]),
        clamp([x0.0, x0.1 + step * 0.1])];
    let mut values: Vec<f64> = simplex
        .iter()
        .map(|x| char_func_xy(x, n_slice, &inv_n, d_slice, rt_slice, rv_slice, lam, pol))
        .collect();

    let alpha = 1.0;
    let gamma = 2.0;
    let rho = 0.5;
    let sigma = 0.5;
    let mut iter = 0;

    loop {
        let mut indices: Vec<usize> = (0..3).collect();
        indices.sort_by(|&i, &j| values[i].partial_cmp(&values[j]).unwrap());
        let (best, good, worst) = (indices[0], indices[1], indices[2]);

        let centroid = [
            (simplex[best][0] + simplex[good][0]) / 2.0,
            (simplex[best][1] + simplex[good][1]) / 2.0,
        ];
        let reflected = clamp([
            centroid[0] + alpha * (centroid[0] - simplex[worst][0]),
            centroid[1] + alpha * (centroid[1] - simplex[worst][1]),
        ]);
        let f_ref = char_func_xy(&reflected, n_slice, &inv_n, d_slice, rt_slice, rv_slice, lam, pol);

        if f_ref < values[best] {
            let expanded = clamp([
                centroid[0] + gamma * (reflected[0] - centroid[0]),
                centroid[1] + gamma * (reflected[1] - centroid[1]),
            ]);
            let f_exp = char_func_xy(&expanded, n_slice, &inv_n, d_slice, rt_slice, rv_slice, lam, pol);
            if f_exp < f_ref {
                simplex[worst] = expanded;
                values[worst] = f_exp;
            } else {
                simplex[worst] = reflected;
                values[worst] = f_ref;
            }
        } else if f_ref < values[good] {
            simplex[worst] = reflected;
            values[worst] = f_ref;
        } else {
            let contracted = clamp([
                centroid[0] + rho * (simplex[worst][0] - centroid[0]),
                centroid[1] + rho * (simplex[worst][1] - centroid[1]),
            ]);
            let f_con = char_func_xy(&contracted, n_slice, &inv_n, d_slice, rt_slice, rv_slice, lam, pol);
            if f_con < values[worst] {
                simplex[worst] = contracted;
                values[worst] = f_con;
            } else {
                for i in 0..3 {
                    if i != best {
                        simplex[i][0] = simplex[best][0] + sigma * (simplex[i][0] - simplex[best][0]);
                        simplex[i][1] = simplex[best][1] + sigma * (simplex[i][1] - simplex[best][1]);
                        values[i] = char_func_xy(&simplex[i], n_slice, &inv_n, d_slice, rt_slice, rv_slice, lam, pol);
                    }
                }
            }
        }

        iter += 1;
        let size = ((simplex[0][0] - simplex[1][0]).powi(2) + (simplex[0][1] - simplex[1][1]).powi(2)).sqrt()
                + ((simplex[1][0] - simplex[2][0]).powi(2) + (simplex[1][1] - simplex[2][1]).powi(2)).sqrt()
                + ((simplex[2][0] - simplex[0][0]).powi(2) + (simplex[2][1] - simplex[0][1]).powi(2)).sqrt();
        if size < tol || iter >= max_iter {
            break;
        }
    }

    let best_idx = (0..3).min_by(|&i, &j| values[i].partial_cmp(&values[j]).unwrap()).unwrap();
    (simplex[best_idx][0], simplex[best_idx][1], values[best_idx])
}

/// Per-layer data for the field-profile sweep.
struct LayerData {
  n: Complex64,
  cos: Complex64,
  thickness: f64,
}

/// `|E(z)|` profile through the stack for one eigenmode.
/// Returns `(z, E, layer_start, layer_end, layer_index)`.
pub fn field_prof(
    n_stack: &[Complex64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    lam: f64,
    n_eff: Complex64,
    pol: i32,
    points_per_layer: usize,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<Complex64>), String>
{
    let (n_slice, d_slice, rt_slice, rv_slice) = (n_stack, thicknesses, rough_types, rough_vals);

    let n_layers = n_slice.len();
    if n_layers < 2 {
        return Err("field profile needs at least 2 layers".to_string());
    }

    let two_pi_lam = 2.0 * PI / lam;

    // Precompute layer data: n, cosθ, thickness
    let mut layers: Vec<LayerData> = Vec::with_capacity(n_layers);
    for i in 0..n_layers {
        let n = n_slice[i];
        let r0 = n_eff * n.recip();
        let v = Complex64::new(1.0, 0.0) - r0 * r0;
        let mut cos = v.sqrt();
        if cos.im < 0.0 {
            cos = -cos;
        }
        layers.push(LayerData {
            n,
            cos,
            thickness: d_slice[i],
        });
    }

    // Helper for Fresnel + roughness at an interface (i -> i+1)
    let interface_props = |i: usize| -> (Complex64, Complex64, Complex64, Complex64) {
        let n_curr = layers[i].n;
        let cos_curr = layers[i].cos;
        let y_curr = if pol == 0 {
            n_curr * cos_curr
        } else {
            let c = if cos_curr.norm() < 1e-12 { Complex64::new(1e-12, 0.0) } else { cos_curr };
            n_curr / c
        };
        let n_next = layers[i+1].n;
        let cos_next = layers[i+1].cos;
        let y_next = if pol == 0 {
            n_next * cos_next
        } else {
            let c = if cos_next.norm() < 1e-12 { Complex64::new(1e-12, 0.0) } else { cos_next };
            n_next / c
        };

        let den = y_curr + y_next;
        let den_safe = if den.norm() < 1e-100 { Complex64::new(1e-100, 1e-100) } else { den };
        let inv_den = den_safe.recip();
        let r12 = (y_curr - y_next) * inv_den;
        let t12 = y_curr * 2.0 * inv_den;
        let t21 = y_next * 2.0 * inv_den;
        let r21 = -r12;

        let sigma = rv_slice[i+1];
        let rtype = rt_slice[i+1];
        if rtype != 0 && sigma > 0.0 {
            let kz1 = two_pi_lam * n_curr * cos_curr;
            let kz2 = two_pi_lam * n_next * cos_next;
            if rtype == 5 {
                // Névot-Croce — shared factors, see `optics_core`. R1.1.
                let (f, ga) = nevot_croce_factors(kz1, kz2, sigma);
                (r12 * f, r21 * f, t12 * ga, t21 * ga)
            } else {
                let al = w_function_inner(2.0 * kz1 * sigma, rtype);
                let be = w_function_inner(2.0 * kz2 * sigma, rtype);
                let ga = w_function_inner((kz1 - kz2) * sigma, rtype);
                (r12 * al, r21 * be, t12 * ga, t21 * ga)
            }
        } else {
            (r12, r21, t12, t21)
        }
    };

    // Propagation phase through a layer (i)
    let prop_phase = |i: usize| -> Complex64 {
        let d = layers[i].thickness;
        if d <= 1e-12 {
            return Complex64::new(1.0, 0.0);
        }
        let mut beta = two_pi_lam * d * layers[i].n * layers[i].cos;
        if beta.im < 0.0 {
            beta = Complex64::new(beta.re, -beta.im);
        }
        (Complex64::new(0.0, 1.0) * beta).exp()
    };

    // ---------- Build left and right S‑matrices ----------
    // S_left[i] = S‑matrix from ambient up to the left side of layer i (i from 1 to n_layers-1)
    let mut s_left: Vec<(Complex64, Complex64, Complex64, Complex64)> = Vec::with_capacity(n_layers);
    s_left.push((Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0))); // identity before ambient

    for i in 0..n_layers-1 {
        let mut sg = *s_left.last().unwrap();
        if i > 0 && layers[i].thickness > 1e-12 {
            let phi = prop_phase(i);
            sg = redheffer_product_complex_field_inner(
                sg.0, sg.1, sg.2, sg.3,
                Complex64::new(0.0, 0.0), phi, phi, Complex64::new(0.0, 0.0),
            );
        }
        let iface = interface_props(i);
        sg = redheffer_product_complex_field_inner(sg.0, sg.1, sg.2, sg.3, iface.0, iface.1, iface.2, iface.3);
        s_left.push(sg);
    }

    // S_right[i] = S‑matrix from substrate up to the right side of layer i (i from n_layers-2 down to 0)
    let mut s_right: Vec<Option<(Complex64, Complex64, Complex64, Complex64)>> = vec![None; n_layers];
    s_right[n_layers-1] = Some((Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)));

    for i in (0..n_layers-1).rev() {
        let mut sg = s_right[i+1].unwrap();
        if i+1 < n_layers-1 && layers[i+1].thickness > 1e-12 {
            let phi = prop_phase(i+1);
            sg = redheffer_product_complex_field_inner(
                Complex64::new(0.0, 0.0), phi, phi, Complex64::new(0.0, 0.0),
                sg.0, sg.1, sg.2, sg.3,
            );
        }
        let iface = interface_props(i);
        sg = redheffer_product_complex_field_inner(iface.0, iface.1, iface.2, iface.3, sg.0, sg.1, sg.2, sg.3);
        s_right[i] = Some(sg);
    }

    // ---------- Compute field inside each layer ----------
    let mut z_pos = Vec::new();
    let mut e_mag = Vec::new();
    let mut layer_start = Vec::new();
    let mut layer_end = Vec::new();
    let mut layer_n = Vec::new();

    let mut z_cursor = 0.0;

    for i in 1..n_layers-1 {
        let d = layers[i].thickness;
        if d <= 1e-12 {
            continue;
        }
        let sl = &s_left[i];
        let sr = s_right[i].as_ref().unwrap();
        let denom = Complex64::new(1.0, 0.0) - sl.3 * sr.0;
        let denom_safe = if denom.norm() < 1e-100 {
            Complex64::new(1e-100, 1e-100)
        } else {
            denom
        };
        let inv_denom = denom_safe.recip();
        let e_plus = sl.2 * inv_denom;
        let e_minus = sr.0 * e_plus;
        let mut beta = two_pi_lam * d * layers[i].n * layers[i].cos;
        if beta.im < 0.0 {
            beta = Complex64::new(beta.re, -beta.im);
        }

        let step = d / (points_per_layer as f64);
        for k in 0..=points_per_layer {
            let zz = k as f64 * step;
            let xi = zz / d;
            let e_z = e_plus * (Complex64::new(0.0, 1.0) * beta * xi).exp()
                    + e_minus * (-Complex64::new(0.0, 1.0) * beta * xi).exp();
            z_pos.push(z_cursor + zz);
            e_mag.push(e_z.norm());
        }
        layer_start.push(z_cursor);
        layer_end.push(z_cursor + d);
        layer_n.push(layers[i].n);
        z_cursor += d;
    }

    // Normalise E‑field to max = 1
    let max_e = e_mag.iter().copied().fold(0.0, f64::max);
    if max_e > 0.0 {
        for val in &mut e_mag {
            *val /= max_e;
        }
    }

    Ok((z_pos, e_mag, layer_start, layer_end, layer_n))
}

// ---------------------------------------------------------------------------
// Tests (standalone: no Python)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
  use super::super::core_engine::{REQ_R_AVG, REQ_RS, REQ_T_AVG};
  use super::*;

  /// (n_layers, n_wavs) row-major complex indices.
  fn ambient_substrate(n_amb: f64, n_sub: f64, nw: usize) -> Vec<Complex64> {
    let mut out = Vec::with_capacity(2 * nw);
    for _ in 0..nw {
      out.push(Complex64::new(n_amb, 0.0));
    }
    for _ in 0..nw {
      out.push(Complex64::new(n_sub, 0.0));
    }
    out
  }

  /// The same dispersive two-layer stack in both index layouts: layer-major
  /// complex for `new`, wav-major interleaved `[re, im]` for
  /// `from_wav_major_flat`. The indices vary with both `l` and `w` on purpose
  /// -- a constant stack would not notice a transposed cache.
  fn both_layouts(nw: usize) -> (Vec<Complex64>, Vec<f64>) {
    let n = |l: usize, w: usize| {
      Complex64::new(1.0 + l as f64 * 0.5 + w as f64 * 0.01, l as f64 * 0.001 * w as f64)
    };
    let mut layer_major = Vec::with_capacity(2 * nw);
    for l in 0..2 {
      for w in 0..nw {
        layer_major.push(n(l, w));
      }
    }
    let mut flat = Vec::with_capacity(4 * nw);
    for w in 0..nw {
      for l in 0..2 {
        flat.push(n(l, w).re);
        flat.push(n(l, w).im);
      }
    }
    (layer_major, flat)
  }

  #[test]
  fn bare_interface_fresnel_hex() {
    // Air (n=1) → glass (n=1.5): Rs(0°) = ((1-1.5)/(1+1.5))² = 0.04.
    let wl = vec![500.0, 600.0];
    let s = Solver::new(
      &wl,
      &[0.0],
      &ambient_substrate(1.0, 1.5, 2),
      2,
      &[0.0, 0.0],
      &[0, 0],
      &[0, 0],
      &[0.0, 0.0],
      2,
    )
    .unwrap();
    let sol = s.solve(REQ_RS | REQ_R_AVG).unwrap();
    let rs = sol.f64maps.iter().find(|(k, _)| k == "Rs").unwrap().1.clone();
    for v in &rs {
      assert!((v - 0.04).abs() < 1e-15, "Rs must be 0.04, got {v}");
    }
    let ra = sol.f64maps.iter().find(|(k, _)| k == "R_avg").unwrap().1.clone();
    for v in &ra {
      assert!((v - 0.04).abs() < 1e-15);
    }
  }

  #[test]
  fn view_masks_match_python() {
    use super::super::core_engine::{
      REQ_A_AVG, REQ_A_P, REQ_A_S, REQ_DELTA_R, REQ_DOP_R, REQ_PSI_R, REQ_RP, REQ_RS,
      REQ_R_AVG, REQ_S0_R, REQ_S1_R, REQ_S2_R, REQ_S3_R, REQ_TS,
    };
    assert_eq!(super::rt_request("s").unwrap(), REQ_R_AVG | REQ_T_AVG | REQ_RS | REQ_TS);
    assert_eq!(
      super::rt_request("u").unwrap(),
      REQ_R_AVG | REQ_T_AVG | REQ_RS | REQ_TS | REQ_RP | super::super::core_engine::REQ_TP
    );
    assert!(super::rt_request("x").is_err());
    assert_eq!(
      super::ellipsometry_request(false),
      REQ_PSI_R | REQ_DELTA_R | REQ_DOP_R | REQ_RS | REQ_RP | REQ_R_AVG
    );
    assert_eq!(super::absorption_request(), REQ_A_S | REQ_A_P | REQ_A_AVG);
    assert_eq!(
      super::stokes_request(true, false).unwrap(),
      REQ_S0_R | REQ_S1_R | REQ_S2_R | REQ_S3_R
    );
    assert!(super::stokes_request(false, false).is_err());
    assert!(super::dispersion_request(false, false, false, false).is_err());
  }

  #[test]
  fn energy_conservation_values() {
    let e = super::energy_conservation(&[0.04], &[0.04], &[0.96], &[0.96]).unwrap();
    assert!((e[0] - 0.0).abs() < 1e-15);
    let e = super::energy_conservation(&[0.5], &[0.4], &[0.3], &[0.4]).unwrap();
    assert!((e[0] - 0.2).abs() < 1e-15);
    assert!(super::energy_conservation(&[0.5], &[0.4], &[0.3], &[]).is_err());
  }

  #[test]
  fn from_raw_broadcast_matches_explicit() {
    let wl = vec![500.0, 600.0];
    let per_layer = vec![Complex64::new(1.0, 0.0), Complex64::new(1.5, 0.0)];
    let a = super::Solver::from_raw(&wl, &[0.0], false, &per_layer, 2, None, None, None, None, 2)
      .unwrap();
    let full = vec![
      Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0),
      Complex64::new(1.5, 0.0), Complex64::new(1.5, 0.0),
    ];
    let b = super::Solver::new(&wl, &[0.0], &full, 2, &[0.0, 0.0], &[0, 0], &[0, 0], &[0.0, 0.0], 2)
      .unwrap();
    use super::super::core_engine::REQ_RS;
    let ra = a.solve(REQ_RS).unwrap();
    let rb = b.solve(REQ_RS).unwrap();
    let fa = ra.f64maps.iter().find(|(k, _)| k == "Rs").unwrap().1.clone();
    let fb = rb.f64maps.iter().find(|(k, _)| k == "Rs").unwrap().1.clone();
    assert_eq!(fa.len(), fb.len());
    for (x, y) in fa.iter().zip(fb.iter()) {
      assert_eq!(x.to_bits(), y.to_bits());
    }
  }

  #[test]
  fn needle_gradient_standalone() {
    use super::super::needle_engine::NREQ_P;
    let wl = vec![500.0, 600.0];
    let idx = vec![
      Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0),
      Complex64::new(2.0, 0.0), Complex64::new(2.0, 0.0),
      Complex64::new(1.5, 0.0), Complex64::new(1.5, 0.0),
    ];
    let s = super::Solver::new(
      &wl, &[0.0], &idx, 3, &[0.0, 100.0, 0.0], &[0, 0, 0], &[0, 0, 0], &[0.0, 0.0, 0.0], 2,
    )
    .unwrap();
    let sol = s
      .needle_gradient(
        &[Complex64::new(2.1, 0.0), Complex64::new(2.1, 0.0)],
        &[10.0, 50.0, 90.0],
        NREQ_P,
        None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None,
        0, Some(2), 0, true, true, None, 0.0,
      )
      .unwrap();
    assert_eq!(sol.n_depths, 3);
    assert_eq!(sol.n_points, 2);
    let ps = sol.maps.iter().find(|(k, _)| k == "P_s").expect("P_s").1.clone();
    assert_eq!(ps.len(), 6);
    assert!(ps.iter().all(|v| v.is_finite()));
  }

  #[test]
  fn eigen_tools_standalone() {
    // 3-layer waveguide-ish stack at one wavelength.
    let col = vec![
      Complex64::new(1.0, 0.0),
      Complex64::new(2.0, 0.0),
      Complex64::new(1.5, 0.0),
    ];
    let th = vec![0.0, 500.0, 0.0];
    let rt = vec![0, 0, 0];
    let rv = vec![0.0, 0.0, 0.0];
    let (re, im, flat) = super::scan_box(&col, &th, &rt, &rv, 600.0, 0, 1.5, 2.0, 0.0, 0.05, 4, 3);
    assert_eq!(re.len(), 4);
    assert_eq!(im.len(), 3);
    assert_eq!(flat.len(), 12);
    assert!(flat.iter().all(|v| v.is_finite()));
    // Synthetic valley: minimum at (re[1], im[1]).
    let synth = vec![5.0, 5.0, 5.0, 5.0, 1.0, 5.0, 5.0, 5.0, 5.0];
    let rr = vec![1.0, 2.0, 3.0];
    let ii = vec![0.0, 0.1, 0.2];
    let mins = super::find_minima(&synth, 3, 3, &rr, &ii, 0.5);
    assert_eq!(mins, vec![(2.0, 0.1)]);
    let (r, i, v) = super::nelder_refine(&col, &th, &rt, &rv, 600.0, 0, (1.7, 0.01), 1e-3, 1e-9, 50);
    assert!(r.is_finite() && i.is_finite() && v.is_finite());
    let prof = super::field_prof(&col, &th, &rt, &rv, 600.0, Complex64::new(1.7, 0.01), 0, 4);
    assert!(prof.is_ok());
    let (z, e, _, _, _) = prof.unwrap();
    assert_eq!(z.len(), e.len());
    assert!(!z.is_empty());
  }

  #[test]
  fn solve_arrays_contract() {
    use super::super::core_engine::REQ_RS;
    let wl = vec![500.0, 600.0];
    let idx = vec![
      Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0),
      Complex64::new(1.5, 0.0), Complex64::new(1.5, 0.0),
    ];
    // Grid mismatch refused.
    assert!(super::solve_arrays(&idx, &[0.0, 0.0], &[false, false], &[0, 0], &[0.0, 0.0],
      &[500.0], &[0.0], false, REQ_RS, 2).is_err());
    // Clean half-spaces: no warnings, Rs = 0.04.
    let (sol, warns) = super::solve_arrays(&idx, &[0.0, 0.0], &[false, false], &[0, 0], &[0.0, 0.0],
      &wl, &[0.0], false, REQ_RS, 2).unwrap();
    assert!(warns.is_empty());
    let rs = sol.f64maps.iter().find(|(k, _)| k == "Rs").unwrap().1.clone();
    assert!((rs[0] - 0.04).abs() < 1e-15);
    // Nonzero half-spaces: warned, not refused.
    let (_, warns) = super::solve_arrays(&idx, &[5.0, 0.0], &[false, false], &[0, 0], &[0.0, 0.0],
      &wl, &[0.0], false, REQ_RS, 2).unwrap();
    assert_eq!(warns.len(), 1);
  }

  #[test]
  fn the_flat_constructor_solves_what_the_layer_major_one_solves() {
    let nw = 64;
    let wl: Vec<f64> = (0..nw).map(|i| 400.0 + i as f64).collect();
    let (layer_major, flat) = both_layouts(nw);
    let a = Solver::new(&wl, &[0.2], &layer_major, 2, &[0.0, 0.0], &[0, 0], &[0, 0],
      &[0.0, 0.0], 2).unwrap();
    let b = Solver::from_wav_major_flat(&wl, &[0.2], &flat, 2, &[0.0, 0.0], &[0, 0], &[0, 0],
      &[0.0, 0.0], 2).unwrap();
    let va = a.solve(REQ_RS).unwrap().f64maps.iter().find(|(k, _)| k == "Rs").unwrap().1.clone();
    let vb = b.solve(REQ_RS).unwrap().f64maps.iter().find(|(k, _)| k == "Rs").unwrap().1.clone();
    // Bit-identical, not merely close: the two paths must differ only in how
    // they walked the caller's buffer (§0.1 rule 5).
    assert_eq!(va, vb);
    // And the curve must actually disperse, or a transposed cache would pass.
    assert!(va[0] != va[nw - 1]);
  }

  #[test]
  fn a_flat_cache_of_the_wrong_length_is_rejected() {
    let wl = vec![500.0, 600.0];
    let (_, flat) = both_layouts(2);
    let ok = Solver::from_wav_major_flat(&wl, &[0.0], &flat, 2, &[0.0, 0.0], &[0, 0], &[0, 0],
      &[0.0, 0.0], 2);
    assert!(ok.is_ok());
    let err = Solver::from_wav_major_flat(&wl, &[0.0], &flat[..flat.len() - 2], 2, &[0.0, 0.0],
      &[0, 0], &[0, 0], &[0.0, 0.0], 2)
    .err()
    .expect("a short cache must be refused");
    assert!(err.contains("index cache length"), "{err}");
  }

  #[test]
  fn the_flat_constructor_refuses_what_new_refuses() {
    let wl = vec![500.0];
    let (_, flat) = both_layouts(1);
    let call = |sin: &[f64], thick: &[f64], mode: i32| {
      Solver::from_wav_major_flat(&wl, sin, &flat, 2, thick, &[0, 0], &[0, 0], &[0.0, 0.0], mode)
    };
    assert!(call(&[0.0], &[0.0, 0.0], 2).is_ok());
    assert!(call(&[], &[0.0, 0.0], 2).is_err(), "no angles");
    assert!(call(&[0.0], &[0.0, 0.0], 5).is_err(), "unknown coherence mode");
    assert!(call(&[0.0], &[0.0], 2).is_err(), "one thickness for two layers");
  }

  #[test]
  fn both_reciprocal_paths_give_the_same_cache() {
    // The parallel branch is a plain map over independent elements -- no
    // reduction, so it owes the serial one bit-identity. Cross the threshold
    // in both directions rather than assume it (R5.3).
    for nw in [8usize, RECIP_PAR_THRESHOLD] {
      let wl: Vec<f64> = (0..nw).map(|i| 400.0 + i as f64 * 0.1).collect();
      let (_, flat) = both_layouts(nw);
      let s = Solver::from_wav_major_flat(&wl, &[0.0], &flat, 2, &[0.0, 0.0], &[0, 0], &[0, 0],
        &[0.0, 0.0], 2)
      .unwrap();
      assert_eq!(s.n_cache.len(), 2 * nw);
      for w in 0..nw {
        for (n, inv) in s.layer_n(w).iter().zip(s.layer_inv_n(w)) {
          assert_eq!(*inv, n.recip(), "wavelength {w}");
        }
      }
    }
  }

  /// A deterministic spread of states: every field varies with the index, so a
  /// split that drops, duplicates or misaligns a range shows up as a wrong
  /// value rather than as a plausible-looking constant.
  fn states_for(n: usize) -> Vec<OpticalState> {
    (0..n)
      .map(|i| {
        let x = i as f64;
        let c = |a: f64, b: f64| Complex64::new(a + x * 1e-4, b - x * 7e-5);
        OpticalState {
          // Kept inside (0, 1) and away from the Rs/Ts floors so the
          // Psi/Delta branches take their transcendental arm.
          rs: 0.05 + 0.4 * ((x * 0.013).sin() * 0.5 + 0.5),
          rp: 0.07 + 0.4 * ((x * 0.017).cos() * 0.5 + 0.5),
          ts: 0.11 + 0.4 * ((x * 0.019).sin() * 0.5 + 0.5),
          tp: 0.13 + 0.4 * ((x * 0.023).cos() * 0.5 + 0.5),
          rs_c: c(0.2, 0.3),
          rp_c: c(-0.4, 0.1),
          ts_c: c(0.5, -0.2),
          tp_c: c(0.6, 0.7),
          rbs_c: c(-0.1, -0.8),
          rbp_c: c(0.9, 0.4),
          tbs_c: c(0.3, -0.6),
          tbp_c: c(-0.7, 0.2),
          cross_r: c(0.15, -0.25),
          cross_t: c(-0.35, 0.45),
        }
      })
      .collect()
  }

  /// Every f64 channel then every complex channel, in the order the sink
  /// struct declares them, all destinations live.
  type AllChannels = (Vec<Vec<f64>>, Vec<Vec<Complex64>>);

  /// Derive `states` into a full set of buffers, either through the halving
  /// (`split = true`) or in one serial call.
  fn derive_all(states: &[OpticalState], split: bool) -> AllChannels {
    let n = states.len();
      let mut b_rs = vec![f64::NAN; n];
      let mut b_rp = vec![f64::NAN; n];
      let mut b_ts = vec![f64::NAN; n];
      let mut b_tp = vec![f64::NAN; n];
      let mut b_ravg = vec![f64::NAN; n];
      let mut b_tavg = vec![f64::NAN; n];
      let mut b_as = vec![f64::NAN; n];
      let mut b_ap = vec![f64::NAN; n];
      let mut b_aavg = vec![f64::NAN; n];
      let mut b_psi_r = vec![f64::NAN; n];
      let mut b_psi_t = vec![f64::NAN; n];
      let mut b_delta_r = vec![f64::NAN; n];
      let mut b_delta_t = vec![f64::NAN; n];
      let mut b_dop_r = vec![f64::NAN; n];
      let mut b_dop_t = vec![f64::NAN; n];
      let mut b_diatt_r = vec![f64::NAN; n];
      let mut b_diatt_t = vec![f64::NAN; n];
      let mut b_s0r = vec![f64::NAN; n];
      let mut b_s1r = vec![f64::NAN; n];
      let mut b_s2r = vec![f64::NAN; n];
      let mut b_s3r = vec![f64::NAN; n];
      let mut b_s0t = vec![f64::NAN; n];
      let mut b_s1t = vec![f64::NAN; n];
      let mut b_s2t = vec![f64::NAN; n];
      let mut b_s3t = vec![f64::NAN; n];
      let mut b_retard_r = vec![f64::NAN; n];
      let mut b_retard_t = vec![f64::NAN; n];
      let mut b_phi_rs = vec![f64::NAN; n];
      let mut b_phi_rp = vec![f64::NAN; n];
      let mut b_phi_ts = vec![f64::NAN; n];
      let mut b_phi_tp = vec![f64::NAN; n];
      let mut b_phi_rbs = vec![f64::NAN; n];
      let mut b_phi_rbp = vec![f64::NAN; n];
      let mut b_phi_tbs = vec![f64::NAN; n];
      let mut b_phi_tbp = vec![f64::NAN; n];
      let mut b_rs_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_rp_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_ts_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_tp_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_rbs_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_rbp_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_tbs_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_tbp_c = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_cross_r = vec![Complex64::new(f64::NAN, f64::NAN); n];
      let mut b_cross_t = vec![Complex64::new(f64::NAN, f64::NAN); n];
    {
      let sinks = Sinks {
          b_rs: Some(&mut b_rs),
          b_rp: Some(&mut b_rp),
          b_ts: Some(&mut b_ts),
          b_tp: Some(&mut b_tp),
          b_ravg: Some(&mut b_ravg),
          b_tavg: Some(&mut b_tavg),
          b_as: Some(&mut b_as),
          b_ap: Some(&mut b_ap),
          b_aavg: Some(&mut b_aavg),
          b_psi_r: Some(&mut b_psi_r),
          b_psi_t: Some(&mut b_psi_t),
          b_delta_r: Some(&mut b_delta_r),
          b_delta_t: Some(&mut b_delta_t),
          b_dop_r: Some(&mut b_dop_r),
          b_dop_t: Some(&mut b_dop_t),
          b_diatt_r: Some(&mut b_diatt_r),
          b_diatt_t: Some(&mut b_diatt_t),
          b_s0r: Some(&mut b_s0r),
          b_s1r: Some(&mut b_s1r),
          b_s2r: Some(&mut b_s2r),
          b_s3r: Some(&mut b_s3r),
          b_s0t: Some(&mut b_s0t),
          b_s1t: Some(&mut b_s1t),
          b_s2t: Some(&mut b_s2t),
          b_s3t: Some(&mut b_s3t),
          b_retard_r: Some(&mut b_retard_r),
          b_retard_t: Some(&mut b_retard_t),
          b_phi_rs: Some(&mut b_phi_rs),
          b_phi_rp: Some(&mut b_phi_rp),
          b_phi_ts: Some(&mut b_phi_ts),
          b_phi_tp: Some(&mut b_phi_tp),
          b_phi_rbs: Some(&mut b_phi_rbs),
          b_phi_rbp: Some(&mut b_phi_rbp),
          b_phi_tbs: Some(&mut b_phi_tbs),
          b_phi_tbp: Some(&mut b_phi_tbp),
          b_rs_c: Some(&mut b_rs_c),
          b_rp_c: Some(&mut b_rp_c),
          b_ts_c: Some(&mut b_ts_c),
          b_tp_c: Some(&mut b_tp_c),
          b_rbs_c: Some(&mut b_rbs_c),
          b_rbp_c: Some(&mut b_rbp_c),
          b_tbs_c: Some(&mut b_tbs_c),
          b_tbp_c: Some(&mut b_tbp_c),
          b_cross_r: Some(&mut b_cross_r),
          b_cross_t: Some(&mut b_cross_t),
      };
      if split {
        derive_par(states, sinks);
      } else {
        derive_range(states, sinks);
      }
    }
    (
      vec![
        b_rs,
        b_rp,
        b_ts,
        b_tp,
        b_ravg,
        b_tavg,
        b_as,
        b_ap,
        b_aavg,
        b_psi_r,
        b_psi_t,
        b_delta_r,
        b_delta_t,
        b_dop_r,
        b_dop_t,
        b_diatt_r,
        b_diatt_t,
        b_s0r,
        b_s1r,
        b_s2r,
        b_s3r,
        b_s0t,
        b_s1t,
        b_s2t,
        b_s3t,
        b_retard_r,
        b_retard_t,
        b_phi_rs,
        b_phi_rp,
        b_phi_ts,
        b_phi_tp,
        b_phi_rbs,
        b_phi_rbp,
        b_phi_tbs,
        b_phi_tbp,
      ],
      vec![
        b_rs_c,
        b_rp_c,
        b_ts_c,
        b_tp_c,
        b_rbs_c,
        b_rbp_c,
        b_tbs_c,
        b_tbp_c,
        b_cross_r,
        b_cross_t,
      ],
    )
  }

  #[test]
  fn dop_r_is_clamped_like_dop_t() {
    // A degree of polarization above 1 is not a number anyone can use. The
    // reflected one was left unclamped while the transmitted one was not
    // (review §3.3); both are clamped since R6.2. Driven here with a cross
    // term far larger than the intensities, which no physical stack produces
    // -- the point is that the guard is in the expression, not that this
    // state is reachable.
    let mut st = states_for(1)[0];
    st.rs = 0.25;
    st.rp = 0.25;
    st.cross_r = Complex64::new(10.0, 3.0);
    st.cross_t = Complex64::new(10.0, 3.0);
    st.ts = 0.25;
    st.tp = 0.25;
    let states = [st];
    let mut dop_r = [f64::NAN];
    let mut dop_t = [f64::NAN];
    {
      let mut sinks = Sinks::none();
      sinks.b_dop_r = Some(&mut dop_r);
      sinks.b_dop_t = Some(&mut dop_t);
      derive_range(&states, sinks);
    }
    assert_eq!(dop_r[0], 1.0, "DOP_R must be clamped");
    assert_eq!(dop_t[0], 1.0, "DOP_T must still be clamped");
  }

  #[test]
  fn a_real_cross_term_puts_delta_r_at_plus_pi_not_minus_pi() {
    // The `+ 0.0` on s2r/s3r is a negative-zero flush, not the no-op the
    // review took it for: `-2.0 * 0.0` is `-0.0`, and `atan2(-0.0, negative)`
    // is -pi where `atan2(+0.0, negative)` is +pi. An isotropic stack at
    // normal incidence has `cross_r` exactly real, so this is the ordinary
    // case; dropping the term would swing Delta_R by 2*pi there and break
    // parity with the numba reference, which flushes at the same four places.
    let mut st = states_for(1)[0];
    st.rs = 0.2;
    st.rp = 0.2;
    st.cross_r = Complex64::new(0.2, 0.0);
    let states = [st];
    let mut delta_r = [f64::NAN];
    {
      let mut sinks = Sinks::none();
      sinks.b_delta_r = Some(&mut delta_r);
      derive_range(&states, sinks);
    }
    assert_eq!(delta_r[0], PI, "Delta_R = {}, expected +pi", delta_r[0]);

    // The transmitted pair carries the same flush with the opposite sign
    // convention, so a real negative cross_t is its +pi case.
    let mut st = states_for(1)[0];
    st.ts = 0.2;
    st.tp = 0.2;
    st.cross_t = Complex64::new(-0.2, 0.0);
    let states = [st];
    let mut delta_t = [f64::NAN];
    {
      let mut sinks = Sinks::none();
      sinks.b_delta_t = Some(&mut delta_t);
      derive_range(&states, sinks);
    }
    assert_eq!(delta_t[0], PI, "Delta_T = {}, expected +pi", delta_t[0]);
  }

  #[test]
  fn the_split_derive_is_bit_identical_to_the_serial_one() {
    // Not a round multiple of the leaf, so the halving ends on odd sizes and
    // an off-by-one in the split would land somewhere.
    let states = states_for(4 * DERIVE_PAR_LEAF + 37);
    let (par_f, par_c) = derive_all(&states, true);
    let (ser_f, ser_c) = derive_all(&states, false);
    // Bit-for-bit: there is no reduction here, so anything less is a bug in
    // the split, not a floating-point tolerance question (§0.1 rule 5).
    assert_eq!(par_f, ser_f);
    assert_eq!(par_c, ser_c);
  }

  #[test]
  fn the_split_reaches_every_index() {
    // The buffers start as NaN. A range the halving never visits stays NaN,
    // and NaN != NaN would not be caught by comparing two equally-gapped
    // runs -- so check coverage directly.
    let states = states_for(3 * DERIVE_PAR_LEAF + 5);
    let (f, c) = derive_all(&states, true);
    for (ch, buf) in f.iter().enumerate() {
      assert_eq!(buf.len(), states.len());
      for (i, v) in buf.iter().enumerate() {
        assert!(!v.is_nan(), "f64 channel {ch} index {i} was never written");
      }
    }
    for (ch, buf) in c.iter().enumerate() {
      for (i, v) in buf.iter().enumerate() {
        assert!(!v.re.is_nan() && !v.im.is_nan(),
          "complex channel {ch} index {i} was never written");
      }
    }
  }

  #[test]
  fn a_grid_below_the_leaf_never_splits_and_still_derives() {
    let states = states_for(DERIVE_PAR_LEAF - 1);
    let (par_f, par_c) = derive_all(&states, true);
    let (ser_f, ser_c) = derive_all(&states, false);
    assert_eq!(par_f, ser_f);
    assert_eq!(par_c, ser_c);
    assert!(!par_f[0].iter().any(|v| v.is_nan()));
  }

  #[test]
  fn a_channel_nobody_asked_for_is_left_alone() {
    // `None` sinks are the common case -- a four-channel photometry request
    // leaves 41 of the 45 unasked -- so the halving must tolerate them at
    // every level, not just the top.
    let states = states_for(2 * DERIVE_PAR_LEAF + 3);
    let n = states.len();
    let mut rs = vec![f64::NAN; n];
    let mut cross = vec![Complex64::new(f64::NAN, f64::NAN); n];
    {
      let mut sinks = Sinks::none();
      sinks.b_rs = Some(&mut rs);
      sinks.b_cross_r = Some(&mut cross);
      derive_par(&states, sinks);
    }
    for (i, s) in states.iter().enumerate() {
      assert_eq!(rs[i], s.rs);
      assert_eq!(cross[i], s.cross_r);
    }
  }

  #[test]
  fn validation_refuses() {
    let wl = vec![500.0];
    let idx = vec![Complex64::new(1.0, 0.0)];
    assert!(Solver::new(&[], &[0.0], &idx, 1, &[0.0], &[0], &[0], &[0.0], 2).is_err());
    assert!(Solver::new(&wl, &[], &idx, 1, &[0.0], &[0], &[0], &[0.0], 2).is_err());
    assert!(Solver::new(&wl, &[0.0], &idx, 1, &[0.0], &[0], &[0], &[0.0], 2).is_err());
    assert!(Solver::new(&wl, &[0.0], &idx, 2, &[0.0], &[0], &[0], &[0.0], 5).is_err());
    let s = Solver::new(&wl, &[0.0], &[Complex64::new(1.0, 0.0), Complex64::new(1.5, 0.0)], 2,
      &[0.0, 0.0], &[0, 0], &[0, 0], &[0.0, 0.0], 2).unwrap();
    assert!(s.solve(0).is_err());
  }
}
