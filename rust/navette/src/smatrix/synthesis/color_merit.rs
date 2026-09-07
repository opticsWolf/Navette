//! synthesis::color_merit — colorimetric demand kernel (Option B plan D2).
//!
//! P1 serves `Lab | XyY` × `DeltaE2000 | DeltaE76 | Channels`. The P2
//! variants (`LCh | Oklab | Y`) exist in the enums so the schema is stable,
//! but evaluate through them is refused here until 0.4.27.
//!
//! All fns are `pub(crate)`: the whole demand lifecycle (compile → merit →
//! needle fold) lives in-crate; the PyO3 surface binds those arms, never
//! this kernel directly (exposure lint stays green by construction).

use serde::{Deserialize, Serialize};

use crate::color::common::{xyz_to_lab, REF_WHITE_D65};
use crate::color::common::xyz_to_srgb;
use crate::color::func_02::{lab_to_lch, lch_to_lab};
use crate::color::func_03::xyz_to_luv;
use crate::color::func_12::lab_to_din99;
use crate::color::func_04::xyz_to_oklab;
use crate::color::func_08::adapt;
use crate::color::func_01::xyz_to_xyy;
use crate::color::func_09::delta_e_76_single;
use crate::color::func_16::delta_e_2000_single;

/// P1: Lab | XyY. P2 adds LCh | Oklab | Y (scalar) — refused in `new`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ColorQuantity {
  Lab,
  XyY,
  LCh,
  Oklab,
  Y,
  /// Dominant wavelength + purity (P3): [wavelength_nm, purity] pair ref,
  /// Channels-style residual off tol[0] (nm) / tol[1] (purity).
  DomWl,
  /// Display-referred sRGB triple (P3): D65-adapted, gamma-encoded,
  /// UNCLIPPED (linear extension past [0,1] keeps gradients honest;
  /// display clipping is presentation, not optimization).
  Srgb,
  /// CIELUV under the demand white (P3).
  Luv,
  /// Raw tristimulus XYZ (P3): linear sensor-space matching.
  Xyz,
  /// DIN99 coordinates (P3, graphics ke = kch = 1): euclideanised Lab.
  Din99,
  /// CIE whiteness [W, Tw] pair (P3): W on the 0-100 scale
  /// (Y is 0-1 here, hence the x100).
  White,
  /// ASTM E313 yellowness index, scalar (P3).
  Yellow,
}

/// E313 YI coefficients for D65/10 deg (ASTM E313): the schema default.
/// Other illuminant/observer geometries need their own table values —
/// pass them explicitly (cx, cz), never silently reuse these.
pub(crate) const E313_CX_D65_10: f64 = 1.3013;
pub(crate) const E313_CZ_D65_10: f64 = 1.1498;

/// P1: all three live; XyY×ΔE and Oklab×ΔE are refused (see matrix).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ColorDistance {
  DeltaE2000,
  DeltaE76,
  Channels,
}

/// Triple reference — or scalar for `Y` (luminance-only, P2).
/// Untagged (`[62, 18, -34]` vs `12.5`); validated against the quantity
/// in `new`. Malformed shapes never reach this type: the compile arm
/// catches them via `ReferenceJson::Other` with a named refusal.
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum ColorReference {
  Triple([f64; 3]),
  /// Dominant-wavelength pair: [wavelength_nm, purity] (purity < 0 flags
  /// the complementary branch — explicit in the value, never a hidden mode).
  Pair([f64; 2]),
  Scalar(f64),
}

/// One compiled color demand. Tables ride NATIVE grids (resampled per eval
/// into scratch — no memo, keeps future holders `Send+Sync`). `white` is
/// the demand illuminant's own white, integrated once here at construction
/// (never adapt foreign-white numbers).
#[derive(Clone, Debug)]
pub struct ColorDemand {
  /// Index into `MeritSpec::keys` (set at compile; groups the demand for
  /// missing-penalty + residual ordering, exactly like pointwise targets).
  pub key_idx: u32,
  pub cmf: Vec<[f64; 3]>,
  pub cmf_wl: Vec<f64>,
  pub illuminant: Vec<f64>,
  pub illum_wl: Vec<f64>,
  pub white: [f64; 3],
  /// DomWl locus polyline ((x, y), wavelength), built from the demand CMF
  /// in `new` (monochromatic chromaticities — illuminant-independent).
  /// Empty for other quantities.
  pub locus: Vec<([f64; 2], f64)>,
  /// E313 coefficients for `Yellow` (defaults: D65/10 deg table).
  pub yi_cx: f64,
  pub yi_cz: f64,
  pub quantity: ColorQuantity,
  pub reference: ColorReference,
  pub distance: ColorDistance,
  /// Per-channel tolerances (Channels mode). Unit default (documented).
  pub tol: [f64; 3],
  pub weight: f64,
}

fn check_grid(name: &str, wl: &[f64], cols: usize) -> Result<(), String> {
  if wl.is_empty() {
    return Err(format!("color: {name} wavelength grid is empty."));
  }
  if wl.windows(2).any(|w| !(w[0] < w[1])) {
    return Err(format!("color: {name} wavelengths not strictly increasing."));
  }
  if wl.iter().any(|v| !v.is_finite()) {
    return Err(format!("color: {name} has non-finite wavelengths."));
  }
  let _ = cols;
  Ok(())
}

impl ColorDemand {
  pub(crate) fn new(
    key_idx: u32,
    cmf: Vec<[f64; 3]>,
    cmf_wl: Vec<f64>,
    illuminant: Vec<f64>,
    illum_wl: Vec<f64>,
    quantity: ColorQuantity,
    reference: ColorReference,
    distance: ColorDistance,
    weight: f64,
    yi_cx: f64,
    yi_cz: f64,
  ) -> Result<Self, String> {
    // Reference shape gates on quantity (all three shapes, both directions).
    match (&quantity, &reference) {
      (ColorQuantity::DomWl | ColorQuantity::White, ColorReference::Pair(_)) => {}
      (ColorQuantity::Y | ColorQuantity::Yellow, ColorReference::Scalar(_)) => {}
      (_, ColorReference::Triple(_))
        if !matches!(quantity, ColorQuantity::DomWl | ColorQuantity::White | ColorQuantity::Y | ColorQuantity::Yellow) => {}
      (ColorQuantity::DomWl, _) => {
        return Err(
          "color: quantity 'DomWl' needs a [2] [wavelength_nm, purity] reference.".to_string(),
        )
      }
      (ColorQuantity::Y, _) => {
        return Err("color: quantity 'Y' needs a scalar reference.".to_string())
      }
      (ColorQuantity::White, _) => {
        return Err("color: quantity 'White' needs a [2] [whiteness, tint] reference.".to_string())
      }
      (ColorQuantity::Yellow, _) => {
        return Err("color: quantity 'Yellow' needs a scalar reference.".to_string())
      }
      (_, ColorReference::Pair(_)) => {
        return Err("color: [2] pair reference needs quantity 'DomWl'|'White'.".to_string())
      }
      (_, _) => return Err("color: scalar reference needs quantity 'Y'.".to_string()),
    }
    // Compat matrix (§D2): ΔE is Lab-space (XyY takes Channels);
    // Oklab-ΔE would double-count (equal-tol Channels is mathematically
    // identical to unweighted Euclidean); Y is scalar (Channels only).
    match (quantity, distance) {
      (ColorQuantity::XyY, ColorDistance::DeltaE2000) | (ColorQuantity::XyY, ColorDistance::DeltaE76) => {
        return Err(format!(
          "color: quantity 'XyY' with distance '{distance:?}' refused (DeltaE is Lab-space; use Channels)."
        ))
      }
      (ColorQuantity::Oklab, ColorDistance::DeltaE2000) | (ColorQuantity::Oklab, ColorDistance::DeltaE76) => {
        return Err(format!(
          "color: quantity 'Oklab' with distance '{distance:?}' refused (use equal-tol Channels)."
        ))
      }
      (ColorQuantity::Y, ColorDistance::DeltaE2000) | (ColorQuantity::Y, ColorDistance::DeltaE76) => {
        return Err(format!(
          "color: quantity 'Y' with distance '{distance:?}' refused (scalar demand takes Channels)."
        ))
      }
      (ColorQuantity::Srgb, ColorDistance::DeltaE2000)
      | (ColorQuantity::Srgb, ColorDistance::DeltaE76)
      | (ColorQuantity::Luv, ColorDistance::DeltaE2000)
      | (ColorQuantity::Luv, ColorDistance::DeltaE76)
      | (ColorQuantity::Xyz, ColorDistance::DeltaE2000)
      | (ColorQuantity::Xyz, ColorDistance::DeltaE76) => {
        return Err(format!(
          "color: quantity '{quantity:?}' with distance '{distance:?}' refused (display/coordinate quantities take Channels)."
        ))
      }
      (ColorQuantity::Din99, ColorDistance::DeltaE2000)
      | (ColorQuantity::Din99, ColorDistance::DeltaE76)
      | (ColorQuantity::White, ColorDistance::DeltaE2000)
      | (ColorQuantity::White, ColorDistance::DeltaE76)
      | (ColorQuantity::Yellow, ColorDistance::DeltaE2000)
      | (ColorQuantity::Yellow, ColorDistance::DeltaE76) => {
        return Err(format!(
          "color: quantity '{quantity:?}' with distance '{distance:?}' refused (DIN99 is already euclidean; whiteness/yellowness are scalar-or-pair indices — use Channels)."
        ))
      }
      (ColorQuantity::DomWl, ColorDistance::DeltaE2000) | (ColorQuantity::DomWl, ColorDistance::DeltaE76) => {
        return Err(format!(
          "color: quantity 'DomWl' with distance '{distance:?}' refused (pair demand takes Channels)."
        ))
      }
      _ => {}
    }
    check_grid("CMF", &cmf_wl, cmf.len())?;
    check_grid("illuminant", &illum_wl, illuminant.len())?;
    if cmf.len() != cmf_wl.len() || illuminant.len() != illum_wl.len() {
      return Err("color: table values length != grid length.".to_string());
    }
    if cmf.iter().any(|c| c.iter().any(|v| !v.is_finite()))
      || illuminant.iter().any(|v| !v.is_finite())
    {
      return Err("color: non-finite CMF/illuminant table value.".to_string());
    }
    if !weight.is_finite() || weight < 0.0 {
      return Err(format!("color: weight must be finite and >= 0, got {weight}."));
    }
    if let ColorReference::Triple(t) = &reference {
      if t.iter().any(|v| !v.is_finite()) {
        return Err("color: non-finite reference triple.".to_string());
      }
    }
    if !yi_cx.is_finite() || !yi_cz.is_finite() {
      return Err("color: non-finite E313 coefficients.".to_string());
    }
    // DomWl locus: monochromatic chromaticities of the demand CMF
    // (illuminant-independent); degenerate tables refused here.
    let locus: Vec<([f64; 2], f64)> = if quantity == ColorQuantity::DomWl {
      let mut loc = Vec::with_capacity(cmf.len());
      for (i, c) in cmf.iter().enumerate() {
        let s = c[0] + c[1] + c[2];
        if s > 0.0 {
          loc.push(([c[0] / s, c[1] / s], cmf_wl[i]));
        }
      }
      if loc.len() < 2 {
        return Err("color: degenerate CMF locus (need >= 2 chromatic points).".to_string());
      }
      loc
    } else {
      Vec::new()
    };
    // Own-white rule: integrate the illuminant as a perfect diffuser on
    // its NATIVE grid (no resample involved — both tables native here).
    let ones = vec![1.0; illuminant.len()];
    let white = xyz_of_spectrum(&ones, &illum_wl, &cmf, &cmf_wl, &illuminant, &illum_wl)?;
    if white.iter().any(|v| !v.is_finite()) {
      return Err("color: degenerate illuminant (non-finite white point).".to_string());
    }
    Ok(ColorDemand {
      key_idx,
      locus,
      yi_cx,
      yi_cz,
      cmf,
      cmf_wl,
      illuminant,
      illum_wl,
      white,
      quantity,
      reference,
      distance,
      tol: [1.0, 1.0, 1.0],
      weight,
    })
  }
}

/// Linear interpolation of a native table at `x`; `None` outside coverage.
fn resample(tbl_wl: &[f64], tbl: &[f64], x: f64) -> Option<f64> {
  if x < tbl_wl[0] || x > tbl_wl[tbl_wl.len() - 1] {
    return None;
  }
  let mut lo = 0usize;
  let mut hi = tbl_wl.len() - 1;
  while hi - lo > 1 {
    let mid = (lo + hi) / 2;
    if tbl_wl[mid] <= x {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  let t = (x - tbl_wl[lo]) / (tbl_wl[hi] - tbl_wl[lo]);
  Some(tbl[lo] + t * (tbl[hi] - tbl[lo]))
}

/// Integration workspace: resampled tables on the covered sim points +
/// the forward-difference Δλ weights. On a uniform fully-covered grid the
/// weights are a constant interval and the sums below reproduce
/// `func_13` summation op-for-op (R1 pins this).
struct XyzWorkspace {
  /// Covered sim indices (ascending) parallel to dw/e_res/cmf_res.
  idx: Vec<usize>,
  xyz: [f64; 3],
  k: f64,
  dw: Vec<f64>,
  e_res: Vec<f64>,
  cmf_res: Vec<[f64; 3]>,
}

fn xyz_workspace(
  sim_row: &[f64],
  sim_wl: &[f64],
  cmf: &[[f64; 3]],
  cmf_wl: &[f64],
  illum: &[f64],
  illum_wl: &[f64],
) -> Result<XyzWorkspace, String> {
  if sim_row.len() != sim_wl.len() || sim_row.is_empty() {
    return Err("color: spectrum length != grid length (or empty).".to_string());
  }
  if sim_row.iter().any(|v| !v.is_finite()) {
    return Err("color: non-finite spectrum value.".to_string());
  }
  let cmf_c: Vec<Vec<f64>> = (0..3).map(|c| cmf.iter().map(|t| t[c]).collect()).collect();
  let mut idx = Vec::new();
  let mut e_res = Vec::new();
  let mut cmf_res = Vec::new();
  for (i, &w) in sim_wl.iter().enumerate() {
    let e = resample(illum_wl, illum, w);
    let x = resample(cmf_wl, &cmf_c[0], w);
    let y = resample(cmf_wl, &cmf_c[1], w);
    let z = resample(cmf_wl, &cmf_c[2], w);
    if let (Some(e), Some(x), Some(y), Some(z)) = (e, x, y, z) {
      idx.push(i);
      e_res.push(e);
      cmf_res.push([x, y, z]);
    }
  }
  if idx.is_empty() {
    return Err("color: no overlap between sim grid and CMF/illuminant tables.".to_string());
  }
  if idx.len() < 2 {
    return Err("color: table overlap is a single point (need >= 2).".to_string());
  }
  // Forward-difference Δλ on the covered subset (uniform grid → constant).
  let wl: Vec<f64> = idx.iter().map(|&i| sim_wl[i]).collect();
  let m = wl.len();
  let mut dw = vec![0.0; m];
  for i in 0..m - 1 {
    dw[i] = wl[i + 1] - wl[i];
  }
  dw[m - 1] = wl[m - 1] - wl[m - 2];
  if dw.iter().any(|&d| !(d > 0.0)) {
    return Err("color: sim grid not strictly increasing over the overlap.".to_string());
  }
  let denom: f64 = (0..m).map(|i| e_res[i] * cmf_res[i][1] * dw[i]).sum();
  if !(denom > 1e-300) {
    return Err("color: degenerate illuminant overlap (k normalization failed).".to_string());
  }
  let k = 1.0 / denom;
  let mut xyz = [0.0; 3];
  for i in 0..m {
    let r = sim_row[idx[i]];
    let w = r * e_res[i] * k * dw[i];
    xyz[0] += w * cmf_res[i][0];
    xyz[1] += w * cmf_res[i][1];
    xyz[2] += w * cmf_res[i][2];
  }
  Ok(XyzWorkspace { idx, xyz, k, dw, e_res, cmf_res })
}

/// `XYZ = Σ R·E·cmf·k·Δλ`, `k = 1/ΣE·ȳ·Δλ` (perfect diffuser ⇒ Y ≡ 1).
pub(crate) fn xyz_of_spectrum(
  sim_row: &[f64],
  sim_wl: &[f64],
  cmf: &[[f64; 3]],
  cmf_wl: &[f64],
  illum: &[f64],
  illum_wl: &[f64],
) -> Result<[f64; 3], String> {
  xyz_workspace(sim_row, sim_wl, cmf, cmf_wl, illum, illum_wl).map(|w| w.xyz)
}

/// XYZ → quantity triple. Lab white = the demand illuminant's own white.
pub(crate) fn color_of_xyz(
  xyz: &[f64; 3],
  white: &[f64; 3],
  quantity: ColorQuantity,
) -> Result<[f64; 3], String> {
  if xyz.iter().any(|v| !v.is_finite()) || white.iter().any(|v| !v.is_finite()) {
    return Err("color: non-finite XYZ/white in quantity map.".to_string());
  }
  match quantity {
    ColorQuantity::Lab => {
      let mut out = [[0.0; 3]];
      xyz_to_lab(&[*xyz], white, &mut out);
      Ok(out[0])
    }
    ColorQuantity::XyY => {
      let mut out = [[0.0; 3]];
      xyz_to_xyy(&[*xyz], &mut out);
      Ok(out[0])
    }
    ColorQuantity::LCh => {
      let mut lab = [[0.0; 3]];
      xyz_to_lab(&[*xyz], white, &mut lab);
      let mut out = [[0.0; 3]];
      lab_to_lch(&lab, &mut out);
      Ok(out[0])
    }
    ColorQuantity::Oklab => {
      // Oklab is D65-defined: Bradford-adapt non-D65 XYZ first via the
      // same `adapt` the bindings (and func_13) use — no clipping, so
      // the FD gradient keeps the smooth map.
      let mut adapted = [[0.0; 3]];
      adapt(&[*xyz], white, &REF_WHITE_D65, false, &mut adapted);
      let mut out = [[0.0; 3]];
      xyz_to_oklab(&adapted, &mut out);
      Ok(out[0])
    }
    ColorQuantity::Y => Err("color: quantity 'Y' is scalar (no triple form).".to_string()),
    ColorQuantity::DomWl => {
      Err("color: quantity 'DomWl' is a [2] pair (no triple form).".to_string())
    }
    ColorQuantity::Srgb => {
      // sRGB is D65-defined (same adapt as Oklab); unclipped gamma
      // extension past [0,1] — exact in gamut, smooth outside it.
      let mut adapted = [[0.0; 3]];
      adapt(&[*xyz], white, &REF_WHITE_D65, false, &mut adapted);
      let mut out = [[0.0; 3]];
      xyz_to_srgb(&adapted, false, &mut out);
      Ok(out[0])
    }
    ColorQuantity::Luv => {
      let mut out = [[0.0; 3]];
      xyz_to_luv(&[*xyz], white, &mut out);
      Ok(out[0])
    }
    ColorQuantity::Xyz => Ok(*xyz),
    ColorQuantity::Din99 => {
      let mut lab = [[0.0; 3]];
      xyz_to_lab(&[*xyz], white, &mut lab);
      let mut out = [[0.0; 3]];
      lab_to_din99(&lab, 1.0, 1.0, &mut out);
      Ok(out[0])
    }
    ColorQuantity::White | ColorQuantity::Yellow => {
      Err("color: quantity 'White'|'Yellow' is a pair/scalar index (no triple form).".to_string())
    }
  }
}

/// CIE whiteness [W, Tw] (W on 0-100: Y here is 0-1). Pure closed form on
/// XYZ + the demand white — illuminant-agnostic by construction.
pub(crate) fn cie_whiteness(xyz: &[f64; 3], white: &[f64; 3]) -> [f64; 2] {
  let s = xyz[0] + xyz[1] + xyz[2];
  let ws = white[0] + white[1] + white[2];
  let (x, y) = (xyz[0] / s, xyz[1] / s);
  let (xn, yn) = (white[0] / ws, white[1] / ws);
  let w = 100.0 * xyz[1] + 800.0 * (xn - x) + 1700.0 * (yn - y);
  let tw = 1000.0 * (xn - x) - 650.0 * (yn - y);
  [w, tw]
}

/// ASTM E313 yellowness index with explicit coefficients.
pub(crate) fn e313_yellowness(xyz: &[f64; 3], cx: f64, cz: f64) -> f64 {
  100.0 * (cx * xyz[0] - cz * xyz[2]) / xyz[1]
}

/// Wrap a hue difference in degrees to [-180, 180] BEFORE scaling
/// (179 vs -179 is 2 deg, not 358).
pub(crate) fn wrap_deg(d: f64) -> f64 {
  d - 360.0 * (d / 360.0).round()
}

/// Scalar objective `F(xyz)`: `w·ΔE²`, resp. `w·Σ((c−c_t)/tol)²`.
/// Dominant wavelength (nm) + purity from XYZ against the demand locus.
///
/// Forward ray (white -> sample) vs locus-only segments: a hit past the
/// sample is spectral (purity = 1/t). A miss means a purple-direction ray;
/// the backward extension then hits the locus (complementary branch,
/// purity NEGATIVE — the branch rule is explicit in the value, never a
/// hidden mode). Achromatic samples (|d| < 1e-9 in xy) return (0, 0):
/// hue carries no information at white, purity does the work (documented
/// kink — any rule kinks there, lambda is undefined at white).
fn dom_wl_purity(demand: &ColorDemand, xyz: &[f64; 3]) -> Result<(f64, f64), String> {
  let s = xyz[0] + xyz[1] + xyz[2];
  if !(s > 0.0) || xyz.iter().any(|v| !v.is_finite()) {
    return Err("color: non-finite/non-positive XYZ in DomWl map.".to_string());
  }
  let ws = demand.white[0] + demand.white[1] + demand.white[2];
  let w = [demand.white[0] / ws, demand.white[1] / ws];
  let px = xyz[0] / s;
  let py = xyz[1] / s;
  let d = [px - w[0], py - w[1]];
  if d[0] * d[0] + d[1] * d[1] < 1e-18 {
    return Ok((0.0, 0.0));
  }
  // Locus pairs carry their own wavelengths (zero-sum CMF rows were
  // skipped at construction — no index alignment needed).
  let locus = &demand.locus;
  let mut forward: Option<(f64, f64)> = None;
  let mut backward: Option<(f64, f64)> = None;
  for wseg in locus.windows(2) {
    let (a, la) = (wseg[0].0, wseg[0].1);
    let (b, lb) = (wseg[1].0, wseg[1].1);
    let e = [b[0] - a[0], b[1] - a[1]];
    let den = d[0] * e[1] - d[1] * e[0];
    if den.abs() < 1e-300 {
      continue;
    }
    let aw = [a[0] - w[0], a[1] - w[1]];
    let t = (aw[0] * e[1] - aw[1] * e[0]) / den;
    let sg = (aw[0] * d[1] - aw[1] * d[0]) / den;
    if !(0.0..=1.0).contains(&sg) {
      continue;
    }
    let lam = la + sg * (lb - la);
    // Forward: smallest t past the sample (1-eps keeps monochrome edge).
    if t > 1.0 - 1e-9 && forward.map_or(true, |(bt, _)| t < bt) {
      forward = Some((t, lam));
    }
    // Backward ray: u = -t form, smallest positive u.
    let u = -t;
    if u > 1e-9 && backward.map_or(true, |(bu, _)| u < bu) {
      backward = Some((u, lam));
    }
  }
  if let Some((t, lam)) = forward {
    return Ok((lam, 1.0 / t));
  }
  if let Some((u, lam)) = backward {
    return Ok((lam, -1.0 / u));
  }
  Err("color: DomWl ray misses the locus both ways (degenerate tables?).".to_string())
}

fn objective_of_xyz(demand: &ColorDemand, xyz: &[f64; 3]) -> Result<f64, String> {
  // Scalar Y: single residual off tol[0] (no triple form involved).
  if demand.quantity == ColorQuantity::Y {
    let ColorReference::Scalar(t) = &demand.reference else {
      return Err("color: quantity 'Y' needs a scalar reference.".to_string());
    };
    let r = (xyz[1] - t) / demand.tol[0];
    return Ok(demand.weight * r * r);
  }
  if demand.quantity == ColorQuantity::DomWl {
    let ColorReference::Pair(t) = &demand.reference else {
      return Err(
        "color: quantity 'DomWl' needs a [2] [wavelength_nm, purity] reference.".to_string(),
      );
    };
    let (dl, dp) = dom_wl_purity(demand, xyz)?;
    let rl = (dl - t[0]) / demand.tol[0];
    let rp = (dp - t[1]) / demand.tol[1];
    return Ok(demand.weight * (rl * rl + rp * rp));
  }
  if demand.quantity == ColorQuantity::White {
    let ColorReference::Pair(t) = &demand.reference else {
      return Err("color: quantity 'White' needs a [2] [whiteness, tint] reference.".to_string());
    };
    if !(xyz[1] > 0.0) || xyz.iter().any(|v| !v.is_finite()) {
      return Err("color: non-finite/non-positive XYZ in whiteness map.".to_string());
    }
    let [w, tw] = cie_whiteness(xyz, &demand.white);
    let rw = (w - t[0]) / demand.tol[0];
    let rt = (tw - t[1]) / demand.tol[1];
    return Ok(demand.weight * (rw * rw + rt * rt));
  }
  if demand.quantity == ColorQuantity::Yellow {
    let ColorReference::Scalar(t) = &demand.reference else {
      return Err("color: quantity 'Yellow' needs a scalar reference.".to_string());
    };
    if !(xyz[1] > 1e-12) || xyz.iter().any(|v| !v.is_finite()) {
      return Err("color: non-finite/near-zero-Y XYZ in yellowness map.".to_string());
    }
    let r = (e313_yellowness(xyz, demand.yi_cx, demand.yi_cz) - t) / demand.tol[0];
    return Ok(demand.weight * r * r);
  }
  let c = color_of_xyz(xyz, &demand.white, demand.quantity)?;
  match demand.distance {
    ColorDistance::DeltaE2000 | ColorDistance::DeltaE76 => {
      let ColorReference::Triple(t) = &demand.reference else {
        return Err("color: scalar reference needs quantity 'Y'.".to_string());
      };
      // ΔE lives in Lab: the LCh ref converts (exact map), and so does
      // the op point (recomputed from XYZ — no trig roundtrip).
      let t_lab;
      let t_ref: &[f64; 3] = if demand.quantity == ColorQuantity::LCh {
        let mut lab = [[0.0; 3]];
        lch_to_lab(&[*t], &mut lab);
        t_lab = lab[0];
        &t_lab
      } else {
        t
      };
      let c_lab_hold;
      let c_lab: &[f64; 3] = if demand.quantity == ColorQuantity::LCh {
        let mut lab = [[0.0; 3]];
        xyz_to_lab(&[*xyz], &demand.white, &mut lab);
        c_lab_hold = lab[0];
        &c_lab_hold
      } else {
        &c
      };
      let d = match demand.distance {
        ColorDistance::DeltaE2000 => delta_e_2000_single(c_lab, t_ref, 1.0, 1.0, 1.0),
        _ => delta_e_76_single(c_lab, t_ref),
      };
      Ok(demand.weight * d * d)
    }
    ColorDistance::Channels => {
      let ColorReference::Triple(t) = &demand.reference else {
        return Err("color: scalar reference needs quantity 'Y'.".to_string());
      };
      let mut s = 0.0;
      for i in 0..3 {
        let mut diff = c[i] - t[i];
        // Hue wraps BEFORE scaling (LCh channel 2, degrees).
        if demand.quantity == ColorQuantity::LCh && i == 2 {
          diff = wrap_deg(diff);
        }
        let r = diff / demand.tol[i];
        s += r * r;
      }
      Ok(demand.weight * s)
    }
  }
}

/// `(residual = √F, grad = ∂F/∂R(λ))`: analytic `dXYZ/dR` (exact, linear)
/// chained with a central 3-pt FD of the XYZ→objective map
/// (`h = 1e-6·(1+|XYZ|)`). Cost: 1 XYZ + 6 tiny map evals.
/// `(residual = √F, covered)` with covered sim indices: the needle fold
/// deposits each gradient value at its own sim point (solver-mapped).
/// `eval_color` is the dense projection (covered order = ascending).
pub(crate) fn eval_color_covered(
  demand: &ColorDemand,
  sim_row: &[f64],
  sim_wl: &[f64],
) -> Result<(f64, Vec<(usize, f64)>), String> {
  let ws = xyz_workspace(sim_row, sim_wl, &demand.cmf, &demand.cmf_wl, &demand.illuminant, &demand.illum_wl)?;
  let f0 = objective_of_xyz(demand, &ws.xyz)?;
  if !f0.is_finite() {
    return Err("color: non-finite objective at op point.".to_string());
  }
  // Central FD of F over XYZ.
  let mut df = [0.0; 3];
  for j in 0..3 {
    let h = 1e-6 * (1.0 + ws.xyz[j].abs());
    let mut xp = ws.xyz;
    xp[j] += h;
    let mut xm = ws.xyz;
    xm[j] -= h;
    let fp = objective_of_xyz(demand, &xp)?;
    let fm = objective_of_xyz(demand, &xm)?;
    if !(fp.is_finite() && fm.is_finite()) {
      return Err("color: non-finite objective under FD step (singular map point).".to_string());
    }
    df[j] = (fp - fm) / (2.0 * h);
  }
  // Chain with analytic dXYZ/dR(λᵢ) = E·cmf·k·Δλ.
  let mut grad = Vec::with_capacity(ws.dw.len());
  for (j, &sim_idx) in ws.idx.iter().enumerate() {
    let s = ws.e_res[j] * ws.k * ws.dw[j];
    let g = s * (df[0] * ws.cmf_res[j][0] + df[1] * ws.cmf_res[j][1] + df[2] * ws.cmf_res[j][2]);
    grad.push((sim_idx, g));
  }
  Ok((f0.sqrt(), grad))
}

pub(crate) fn eval_color(
  demand: &ColorDemand,
  sim_row: &[f64],
  sim_wl: &[f64],
) -> Result<(f64, Vec<f64>), String> {
  eval_color_covered(demand, sim_row, sim_wl)
    .map(|(r, v)| (r, v.into_iter().map(|(_, g)| g).collect()))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::color::tables::default_tables;

  /// Tiny native toy tables (uniform 10 nm): independent of the embed.
  fn toy() -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>, Vec<f64>) {
    let wl: Vec<f64> = (0..8).map(|i| 500.0 + 10.0 * i as f64).collect();
    let cmf: Vec<[f64; 3]> = (0..8)
      .map(|i| [0.10 + 0.01 * i as f64, 0.20 + 0.01 * i as f64, 0.05 + 0.005 * i as f64])
      .collect();
    let illum = vec![1.0; 8];
    (wl.clone(), cmf, wl, illum)
  }

  fn toy_demand() -> ColorDemand {
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    ColorDemand::new(
      0,
      cmf,
      cmf_wl,
      illuminant,
      illum_wl,
      ColorQuantity::Lab,
      ColorReference::Triple([60.0, 10.0, -20.0]),
      ColorDistance::DeltaE2000,
      2.0, E313_CX_D65_10, E313_CZ_D65_10,
    )
    .unwrap()
  }

  #[test]
  fn white_identity_y_exactly_one() {
    // k-normalization self-check: perfect diffuser ⇒ Y ≡ 1 (same sums,
    // same order, numerator and denominator — bitwise, not approximate).
    let d = toy_demand();
    let ones = vec![1.0; 8];
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    let xyz = xyz_of_spectrum(&ones, &cmf_wl, &cmf, &cmf_wl, &illuminant, &illum_wl).unwrap();
    // 1-ulp, not bitwise: the per-term k multiplication rounds once.
    assert!((xyz[1] - 1.0).abs() < 1e-15);
    assert!((d.white[1] - 1.0).abs() < 1e-15);
    assert!(xyz[0].is_finite() && xyz[2].is_finite());
  }

  #[test]
  fn ramp_matches_independent_trapezoid() {
    // Independent reimplementation (trapezoid, reverse summation): same
    // math, different formula AND different rounding path — agreement
    // bounds both formula and summation-order error.
    let (wl, cmf, _, illum) = toy();
    let ramp: Vec<f64> = (0..8).map(|i| 0.2 + 0.05 * i as f64).collect();
    let xyz = xyz_of_spectrum(&ramp, &wl, &cmf, &wl, &illum, &wl).unwrap();
    let h: f64 = 10.0;
    let denom: f64 = (0..8).rev().map(|i| illum[i] * cmf[i][1] * h).sum();
    let k = 1.0 / denom;
    for c in 0..3 {
      // Trapezoid differs (endpoint halves) — bound, don't equal.
      let trap_sum: f64 = (0..8)
        .map(|i| {
          let w: f64 = if i == 0 || i == 7 { 0.5 } else { 1.0 };
          w * ramp[i] * illum[i] * cmf[i][c]
        })
        .sum();
      let trap = h * trap_sum;
      // Coarse 8-pt grid: endpoint-halving legitimately differs ~15%
      // (the rectangular reverse twin below is the strict one).
      let rel = ((xyz[c] / k) - trap).abs() / trap.abs().max(1e-300);
      assert!(rel < 0.20, "c={c} rel={rel}");
    }
    // Rectangular twin (same formula, reverse order): tight.
    for c in 0..3 {
      let rect: f64 = (0..8).rev().map(|i| ramp[i] * illum[i] * cmf[i][c] * h).sum();
      let rel = ((xyz[c] / k) - rect).abs() / rect.abs().max(1e-300);
      assert!(rel < 1e-12, "c={c} rel={rel}");
    }
  }

  #[test]
  fn quantity_wrappers_are_bitwise() {
    let d = toy_demand();
    let xyz = [0.3, 0.4, 0.2];
    let mut lab = [[0.0; 3]];
    xyz_to_lab(&[xyz], &d.white, &mut lab);
    assert_eq!(color_of_xyz(&xyz, &d.white, ColorQuantity::Lab).unwrap(), lab[0]);
    let mut xyy = [[0.0; 3]];
    xyz_to_xyy(&[xyz], &mut xyy);
    assert_eq!(color_of_xyz(&xyz, &d.white, ColorQuantity::XyY).unwrap(), xyy[0]);
  }

  #[test]
  fn de_wrappers_are_bitwise() {
    let d = toy_demand();
    let xyz = [0.3, 0.4, 0.2];
    let lab = color_of_xyz(&xyz, &d.white, ColorQuantity::Lab).unwrap();
    let t = [60.0, 10.0, -20.0];
    assert_eq!(delta_e_2000_single(&lab, &t, 1.0, 1.0, 1.0), {
      let f = objective_of_xyz(&d, &xyz).unwrap();
      (f / 2.0).sqrt()
    });
  }

  #[test]
  fn gradient_matches_brute_force_bump() {
    // Analytic-Jacobian+3pt-FD g(λ) vs per-λ bumped full evals.
    let d = toy_demand();
    let (wl, _, _, _) = toy();
    let row = vec![0.5, 0.55, 0.6, 0.65, 0.6, 0.55, 0.5, 0.45];
    let (resid, grad) = eval_color(&d, &row, &wl).unwrap();
    assert!(resid.is_finite() && resid > 0.0);
    assert_eq!(grad.len(), 8);
    let delta = 1e-5;
    for i in 0..8 {
      let mut rp = row.clone();
      rp[i] += delta;
      let mut rm = row.clone();
      rm[i] -= delta;
      let fp = objective_of_xyz(&d, &xyz_of_spectrum(&rp, &wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap()).unwrap();
      let fm = objective_of_xyz(&d, &xyz_of_spectrum(&rm, &wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap()).unwrap();
      let brute = (fp - fm) / (2.0 * delta);
      let rel = (grad[i] - brute).abs() / brute.abs().max(1e-300);
      assert!(rel < 1e-6, "λ{i} rel={rel}");
    }
  }

  #[test]
  fn refusals_name_the_culprit() {
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    let mk = |q, r, dist| {
      ColorDemand::new(
        0, cmf.clone(), cmf_wl.clone(), illuminant.clone(), illum_wl.clone(), q, r, dist, 1.0, E313_CX_D65_10, E313_CZ_D65_10,
      )
    };
    assert!(mk(ColorQuantity::XyY, ColorReference::Triple([0.3, 0.3, 0.5]), ColorDistance::DeltaE2000)
      .unwrap_err()
      .contains("XyY"));
    // P2 quantities construct (compat matrix below gates the pairs).
    for q in [ColorQuantity::LCh, ColorQuantity::Oklab] {
      assert!(mk(q, ColorReference::Triple([0.0; 3]), ColorDistance::Channels).is_ok());
    }
    assert!(mk(ColorQuantity::Y, ColorReference::Scalar(0.5), ColorDistance::Channels).is_ok());
    assert!(mk(ColorQuantity::Oklab, ColorReference::Triple([0.0; 3]), ColorDistance::DeltaE2000)
      .unwrap_err()
      .contains("Oklab"));
    assert!(mk(ColorQuantity::Y, ColorReference::Scalar(0.5), ColorDistance::DeltaE76)
      .unwrap_err()
      .contains("'Y'"));
    // Shape gate fires before the P2 gate (malformed regardless of phase).
    assert!(mk(ColorQuantity::Y, ColorReference::Triple([0.0; 3]), ColorDistance::Channels)
      .unwrap_err()
      .contains("scalar"));
    assert!(mk(ColorQuantity::Lab, ColorReference::Scalar(1.0), ColorDistance::Channels)
      .unwrap_err()
      .contains("scalar reference"));
    // Empty overlap / single point / non-finite.
    let d = toy_demand();
    let far: Vec<f64> = (0..8).map(|i| 1000.0 + 10.0 * i as f64).collect();
    assert!(xyz_of_spectrum(&vec![0.5; 8], &far, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl)
      .unwrap_err()
      .contains("no overlap"));
    let mut bad = vec![0.5; 8];
    bad[3] = f64::NAN;
    assert!(eval_color(&d, &bad, &toy().0).unwrap_err().contains("non-finite"));
  }

  #[test]
  fn embedded_defaults_recover_d65_white() {
    // End-to-end embed sanity: D65 white from embedded tables ≈ the
    // textbook D65 chromaticity (Y ≡ 1 by construction).
    let t = default_tables();
    assert_eq!(t.cmf_wl.len(), 471);
    assert_eq!(t.cmf_wl[0], 360.0);
    assert_eq!(t.cmf_wl[470], 830.0);
    let ones = vec![1.0; t.illum_wl.len()];
    let white =
      xyz_of_spectrum(&ones, &t.illum_wl, &t.cmf_xyz, &t.cmf_wl, &t.illum, &t.illum_wl).unwrap();
    assert!((white[1] - 1.0).abs() < 1e-12);
    assert!((white[0] - 0.9505).abs() < 1e-3, "X={}", white[0]);
    assert!((white[2] - 1.0890).abs() < 1e-3, "Z={}", white[2]);
  }

  fn p2_demand(q: ColorQuantity, r: ColorReference) -> ColorDemand {
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    ColorDemand::new(0, cmf, cmf_wl, illuminant, illum_wl, q, r, ColorDistance::Channels, 1.0, E313_CX_D65_10, E313_CZ_D65_10)
      .unwrap()
  }

  #[test]
  fn lch_wrappers_are_bitwise() {
    let d = p2_demand(ColorQuantity::LCh, ColorReference::Triple([60.0, 20.0, 100.0]));
    let xyz = [0.3, 0.4, 0.2];
    let mut lab = [[0.0; 3]];
    xyz_to_lab(&[xyz], &d.white, &mut lab);
    let mut lch = [[0.0; 3]];
    crate::color::func_02::lab_to_lch(&lab, &mut lch);
    assert_eq!(color_of_xyz(&xyz, &d.white, ColorQuantity::LCh).unwrap(), lch[0]);
    // ΔE path converts the LCh ref to Lab (roundtrip identity pins it).
    let mut back = [[0.0; 3]];
    crate::color::func_02::lch_to_lab(&[lch[0]], &mut back);
    for i in 0..3 {
      assert!((back[0][i] - lab[0][i]).abs() < 1e-12);
    }
  }

  #[test]
  fn oklab_wrapper_is_bitwise_and_white_is_unit() {
    use crate::color::common::REF_WHITE_D65;
    let d = p2_demand(ColorQuantity::Oklab, ColorReference::Triple([0.5, 0.01, -0.02]));
    let xyz = [0.3, 0.4, 0.2];
    let mut adapted = [[0.0; 3]];
    adapt(&[xyz], &d.white, &REF_WHITE_D65, false, &mut adapted);
    let mut ok = [[0.0; 3]];
    crate::color::func_04::xyz_to_oklab(&adapted, &mut ok);
    assert_eq!(color_of_xyz(&xyz, &d.white, ColorQuantity::Oklab).unwrap(), ok[0]);
    // Full Bradford path (clearly non-D65 toy white): maps onto D65 to
    // rounding. (Near-D65 whites take adapt's copy short-circuit instead —
    // the common D65-demand regime, exact by construction.)
    let (twl, tcmf, tewl, till) = toy();
    let ones = vec![1.0; 8];
    let tw = xyz_of_spectrum(&ones, &tewl, &tcmf, &twl, &till, &tewl).unwrap();
    let mut w65 = [[0.0; 3]];
    adapt(&[tw], &tw, &REF_WHITE_D65, false, &mut w65);
    for i in 0..3 {
      assert!((w65[0][i] - REF_WHITE_D65[i]).abs() < 1e-12, "{w65:?}");
    }
    // D65 white is Oklab near-neutral — up to a PRE-EXISTING in-tree
    // constant/matrix mismatch (REF_WHITE_D65 is the 2-dp-chromaticity
    // derivation [0.95045593, 1, 1.08905775]; the Oklab matrix is
    // native-white, hence b ≈ -1.2e-4). Systematic, ~1e-4 in b, far below
    // demand relevance; kernel and oracle agree bitwise regardless
    // (proven via the bound `_color` twins in R3).
    let mut wok = [[0.0; 3]];
    crate::color::func_04::xyz_to_oklab(&[REF_WHITE_D65], &mut wok);
    assert!((wok[0][0] - 1.0).abs() < 3e-6, "{wok:?}");
    assert!(wok[0][1].abs() < 3e-4 && wok[0][2].abs() < 3e-4, "{wok:?}");
  }

  #[test]
  fn lch_de_measures_in_lab() {
    // The ΔE arm converts BOTH sides to Lab (regression: c once leaked
    // through in LCh — the R3 twin caught it).
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    let wl = cmf_wl.clone();
    let d = ColorDemand::new(
      0, cmf, cmf_wl, illuminant, illum_wl,
      ColorQuantity::LCh,
      ColorReference::Triple([60.0, 20.0, 100.0]),
      ColorDistance::DeltaE76,
      1.0, E313_CX_D65_10, E313_CZ_D65_10,
    )
    .unwrap();
    let row = vec![0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25];
    let (r, _) = eval_color(&d, &row, &wl).unwrap();
    let xyz = xyz_of_spectrum(&row, &wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap();
    let mut lab = [[0.0; 3]];
    xyz_to_lab(&[xyz], &d.white, &mut lab);
    let mut ref_lab = [[0.0; 3]];
    lch_to_lab(&[[60.0, 20.0, 100.0]], &mut ref_lab);
    assert_eq!(r, delta_e_76_single(&lab[0], &ref_lab[0]));
  }

  #[test]
  fn y_scalar_is_exact_and_hue_wraps() {
    // wrap unit: 179 vs -179 differ by 2 deg, not 358.
    assert_eq!(wrap_deg(358.0), -2.0);
    assert_eq!(wrap_deg(-358.0), 2.0);
    // Y residual == hand formula, bitwise.
    let d = p2_demand(ColorQuantity::Y, ColorReference::Scalar(0.5));
    let (wl, _, _, _) = toy();
    let row = vec![0.6; 8];
    let (r, g) = eval_color(&d, &row, &wl).unwrap();
    let xyz = xyz_of_spectrum(&row, &wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap();
    let e = ((xyz[1] - 0.5) / 1.0).abs();
    assert_eq!(r, e);
    assert_eq!(g.len(), 8);
    // Hue wrap through the objective: h=179 vs ref h=-179 -> 2 deg gap.
    let dh = p2_demand(ColorQuantity::LCh, ColorReference::Triple([60.0, 20.0, -179.0]));
    let xyz_h = [0.3, 0.4, 0.2];
    let f = objective_of_xyz(&dh, &xyz_h).unwrap();
    let c = color_of_xyz(&xyz_h, &dh.white, ColorQuantity::LCh).unwrap();
    let gap = (wrap_deg(c[2] - -179.0) / 1.0).powi(2);
    let rest = ((c[0] - 60.0).powi(2) + (c[1] - 20.0).powi(2)).max(0.0);
    assert!((f - (rest + gap)).abs() < 1e-12, "{f}");
  }

  #[test]
  fn achromatic_gradients_stay_finite() {
    // C=0 singularity (LCh) + cbrt neutrals (Oklab): uniform spectra land
    // on the achromatic axis — gradients must stay finite (FD-twin class).
    let (wl, _, _, _) = toy();
    let row = vec![0.5; 8];
    for (q, r) in [
      (ColorQuantity::LCh, ColorReference::Triple([50.0, 5.0, 10.0])),
      (ColorQuantity::Oklab, ColorReference::Triple([0.5, 0.0, 0.0])),
      (ColorQuantity::Y, ColorReference::Scalar(0.4)),
    ] {
      let (cmf_wl, cmf, illum_wl, illuminant) = toy();
      let d = ColorDemand::new(0, cmf, cmf_wl, illuminant, illum_wl, q, r, ColorDistance::Channels, 1.0, E313_CX_D65_10, E313_CZ_D65_10)
        .unwrap();
      let (resid, grad) = eval_color(&d, &row, &wl).unwrap();
      assert!(resid.is_finite());
      assert!(grad.iter().all(|v| v.is_finite()), "{q:?}");
    }
  }

  fn domwl_demand() -> ColorDemand {
    let t = default_tables();
    ColorDemand::new(
      0,
      t.cmf_xyz.clone(),
      t.cmf_wl.clone(),
      t.illum.clone(),
      t.illum_wl.clone(),
      ColorQuantity::DomWl,
      ColorReference::Pair([550.0, 0.9]),
      ColorDistance::Channels,
      1.0, E313_CX_D65_10, E313_CZ_D65_10,
    )
    .unwrap()
  }

  #[test]
  fn domwl_spectral_complementary_and_white() {
    let d = domwl_demand();
    assert!(d.locus.len() >= 400);
    let wl = &d.illum_wl;
    // Narrow lobe at 550 nm: spectral, high purity.
    let gauss: Vec<f64> = wl.iter().map(|w| (-((w - 550.0) / 15.0).powi(2)).exp()).collect();
    let xyz = xyz_of_spectrum(&gauss, wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap();
    let (lam, p) = dom_wl_purity(&d, &xyz).unwrap();
    assert!((lam - 550.0).abs() < 3.0, "{lam}");
    assert!(p > 0.8, "{p}");
    // Twin magenta lobes (450 + 650): purple direction -> complementary.
    let mag: Vec<f64> = wl
      .iter()
      .map(|w| {
        (-((w - 450.0) / 12.0).powi(2)).exp() + (-((w - 650.0) / 12.0).powi(2)).exp()
      })
      .collect();
    let xyz_m = xyz_of_spectrum(&mag, wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap();
    let (lam_c, p_c) = dom_wl_purity(&d, &xyz_m).unwrap();
    assert!(p_c < 0.0, "{p_c}");
    assert!((490.0..570.0).contains(&lam_c), "{lam_c}");
    // Achromatic: hue carries nothing (rule returns (0, 0)).
    let ones = vec![1.0; wl.len()];
    let xyz_w = xyz_of_spectrum(&ones, wl, &d.cmf, &d.cmf_wl, &d.illuminant, &d.illum_wl).unwrap();
    assert_eq!(dom_wl_purity(&d, &xyz_w).unwrap(), (0.0, 0.0));
    // Gradient finite on the smooth (spectral, mid-locus) point.
    let (resid, grad) = eval_color(&d, &gauss, wl).unwrap();
    assert!(resid.is_finite());
    assert!(grad.iter().all(|v| v.is_finite()));
  }

  #[test]
  fn domwl_shape_and_pair_gates() {
    let t = default_tables();
    let mk = |q, r: ColorReference| {
      ColorDemand::new(
        0,
        t.cmf_xyz.clone(),
        t.cmf_wl.clone(),
        t.illum.clone(),
        t.illum_wl.clone(),
        q,
        r,
        ColorDistance::Channels,
        1.0, E313_CX_D65_10, E313_CZ_D65_10,
      )
    };
    assert!(mk(ColorQuantity::DomWl, ColorReference::Pair([550.0, 0.9])).is_ok());
    assert!(mk(ColorQuantity::DomWl, ColorReference::Triple([0.0; 3])).unwrap_err().contains("DomWl"));
    assert!(mk(ColorQuantity::Lab, ColorReference::Pair([550.0, 0.9])).unwrap_err().contains("pair reference"));
    let de = ColorDemand::new(
      0,
      t.cmf_xyz.clone(),
      t.cmf_wl.clone(),
      t.illum.clone(),
      t.illum_wl.clone(),
      ColorQuantity::DomWl,
      ColorReference::Pair([550.0, 0.9]),
      ColorDistance::DeltaE2000,
      1.0, E313_CX_D65_10, E313_CZ_D65_10,
    );
    assert!(de.unwrap_err().contains("pair demand takes Channels"));
  }

  #[test]
  fn coordinate_quantities_are_bitwise_and_smooth() {
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    let ones = vec![1.0; 8];
    let white =
      xyz_of_spectrum(&ones, &illum_wl, &cmf, &cmf_wl, &illuminant, &illum_wl).unwrap();
    let xyz = [0.3, 0.4, 0.2];
    // sRGB: same adapt + unclipped map (bitwise); negatives extend
    // linearly (finite, monotone — no clip flat-spots, no NaN).
    let mut adapted = [[0.0; 3]];
    adapt(&[xyz], &white, &crate::color::common::REF_WHITE_D65, false, &mut adapted);
    let mut rgb = [[0.0; 3]];
    xyz_to_srgb(&adapted, false, &mut rgb);
    assert_eq!(color_of_xyz(&xyz, &white, ColorQuantity::Srgb).unwrap(), rgb[0]);
    let far = color_of_xyz(&[-0.5, 0.1, 2.0], &white, ColorQuantity::Srgb).unwrap();
    assert!(far.iter().all(|v| v.is_finite()));
    assert!(far[0] < 0.0 && far[2] > 1.0); // honest extension, not clipped
    // Luv under the demand white (bitwise); XYZ identity.
    let mut luv = [[0.0; 3]];
    xyz_to_luv(&[xyz], &white, &mut luv);
    assert_eq!(color_of_xyz(&xyz, &white, ColorQuantity::Luv).unwrap(), luv[0]);
    assert_eq!(color_of_xyz(&xyz, &white, ColorQuantity::Xyz).unwrap(), xyz);
    // FD gradients finite on a smooth point (affine toy is fine here).
    let (wl, _, _, _) = toy();
    let row = vec![0.5, 0.55, 0.6, 0.65, 0.6, 0.55, 0.5, 0.45];
    for (q, r) in [
      (ColorQuantity::Srgb, ColorReference::Triple([0.5, 0.4, 0.6])),
      (ColorQuantity::Luv, ColorReference::Triple([50.0, 10.0, -20.0])),
      (ColorQuantity::Xyz, ColorReference::Triple([0.3, 0.4, 0.2])),
    ] {
      let (cmf_wl2, cmf2, illum_wl2, illuminant2) = toy();
      let d = ColorDemand::new(
        0, cmf2, cmf_wl2, illuminant2, illum_wl2, q, r, ColorDistance::Channels, 1.0, E313_CX_D65_10, E313_CZ_D65_10,
      )
      .unwrap();
      let (resid, grad) = eval_color(&d, &row, &wl).unwrap();
      assert!(resid.is_finite(), "{q:?}");
      assert!(grad.iter().all(|v| v.is_finite()), "{q:?}");
    }
  }

  fn idx_demand(q: ColorQuantity, r: ColorReference) -> ColorDemand {
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    ColorDemand::new(
      0, cmf, cmf_wl, illuminant, illum_wl, q, r, ColorDistance::Channels, 1.0,
      E313_CX_D65_10, E313_CZ_D65_10,
    )
    .unwrap()
  }

  #[test]
  fn din99_matches_delta_path_and_white_yellow_hand() {
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    let wl = cmf_wl.clone();
    let ones = vec![1.0; 8];
    let white =
      xyz_of_spectrum(&ones, &illum_wl, &cmf, &cmf_wl, &illuminant, &illum_wl).unwrap();
    let xyz = [0.3, 0.4, 0.2];
    // Din99 wrapper bitwise vs direct calls.
    let mut lab = [[0.0; 3]];
    xyz_to_lab(&[xyz], &white, &mut lab);
    let mut d99 = [[0.0; 3]];
    crate::color::func_12::lab_to_din99(&lab, 1.0, 1.0, &mut d99);
    assert_eq!(color_of_xyz(&xyz, &white, ColorQuantity::Din99).unwrap(), d99[0]);
    // Whiteness/yellowness closed forms, bitwise vs hand evaluation.
    let [w, tw] = cie_whiteness(&xyz, &white);
    let s = xyz[0] + xyz[1] + xyz[2];
    let ws = white[0] + white[1] + white[2];
    assert_eq!(w, 100.0 * xyz[1] + 800.0 * (white[0] / ws - xyz[0] / s) + 1700.0 * (white[1] / ws - xyz[1] / s));
    assert_eq!(tw, 1000.0 * (white[0] / ws - xyz[0] / s) - 650.0 * (white[1] / ws - xyz[1] / s));
    assert_eq!(e313_yellowness(&xyz, 1.3013, 1.1498), 100.0 * (1.3013 * xyz[0] - 1.1498 * xyz[2]) / xyz[1]);
    // Perfect diffuser hits W = 100, Tw = 0 (bitwise — same fractions).
    let [w1, t1] = cie_whiteness(&white, &white);
    // 1-ulp via white[1] (per-term k rounding, same class as the Y identity).
    assert!((w1 - 100.0).abs() < 1e-12, "{w1}");
    assert_eq!(t1, 0.0);
    // End-to-end residuals finite on a smooth row (all three).
    let row = vec![0.5, 0.55, 0.6, 0.65, 0.6, 0.55, 0.5, 0.45];
    for d in [
      idx_demand(ColorQuantity::Din99, ColorReference::Triple([50.0, 5.0, 5.0])),
      idx_demand(ColorQuantity::White, ColorReference::Pair([90.0, 2.0])),
      idx_demand(ColorQuantity::Yellow, ColorReference::Scalar(5.0)),
    ] {
      let (resid, grad) = eval_color(&d, &row, &wl).unwrap();
      assert!(resid.is_finite());
      assert!(grad.iter().all(|v| v.is_finite()));
    }
  }

  #[test]
  fn index_gates_name_all_sides() {
    let bad_q = |q, r: ColorReference, d: ColorDistance| {
      let (cmf_wl, cmf, illum_wl, illuminant) = toy();
      ColorDemand::new(0, cmf, cmf_wl, illuminant, illum_wl, q, r, d, 1.0, 1.3013, 1.1498)
    };
    assert!(bad_q(ColorQuantity::White, ColorReference::Triple([0.0; 3]), ColorDistance::Channels)
      .unwrap_err()
      .contains("White"));
    assert!(bad_q(ColorQuantity::Yellow, ColorReference::Triple([0.0; 3]), ColorDistance::Channels)
      .unwrap_err()
      .contains("scalar"));
    assert!(bad_q(ColorQuantity::Din99, ColorReference::Triple([50.0; 3]), ColorDistance::DeltaE76)
      .unwrap_err()
      .contains("Channels"));
    // Explicit coefficients ride through (F2 geometry would pass its own).
    let (cmf_wl, cmf, illum_wl, illuminant) = toy();
    let d = ColorDemand::new(
      0, cmf, cmf_wl, illuminant, illum_wl,
      ColorQuantity::Yellow, ColorReference::Scalar(5.0), ColorDistance::Channels, 1.0,
      1.2769, 1.0592,
    )
    .unwrap();
    assert_eq!((d.yi_cx, d.yi_cz), (1.2769, 1.0592));
  }
}
