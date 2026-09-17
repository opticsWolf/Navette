// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::merit — MeritSpec: flat optimization targets + residual kernel.
//!
//! Bridge between navette_spectralweave's TargetWeaver and the synthesis
//! loop. The residual math below is lifted VERBATIM from
//! `calculate_merit` (navette_spectralweave/src/targetweaver.rs):
//!   * Exact/Above/Below/Range/CenterBand activation (`kind`)
//!   * Linear/Log/Phase/Complex sim-side transforms (`transform`)
//!   * bit-exact aligned-grid fast path + two-pointer monotone interpolation
//!   * overlap-skip and missing-key penalty semantics
//!
//! Residual space is the square root of merit space: `merit()` sums `r²`,
//! so the Range/CenterBand arms return `0`/`(ad-bw)/tol` and `d/bw` /
//! `sqrt(((ad-bw)/tol)² + 1)` respectively — the `+1` under the root is the
//! band-edge continuity level from `calculate_merit`.
//! Normalization itself is NOT re-implemented: it happened at ingestion in
//! TargetWeaver::register_metadata, and the Python converter copies the
//! finished (normalized_targets, norm_factor, floored tolerances) over.
//!
//! Differences from calculate_merit (by design):
//!   * No OpticalWeaver: simulation curves arrive as SimCurves — plain
//!     [n_angles, n_wav] row-major slices indexed by (pol, channel).
//!   * Entries are pre-flattened into (key-group, target) lists; the
//!     missing-curve penalty is applied ONCE PER KEY GROUP, matching the
//!     per-key `continue` in calculate_merit.
//!   * Angle rows resolved with argmin(|angles − key_angle|), identical to
//!     `_collect_target_angles` + row selection in needle_synthesis.py.

use std::sync::Arc;

use num_complex::Complex64;

use super::color_merit::{ColorDemand, eval_color};

// ---------------------------------------------------------------------------
// Vocabulary types
// ---------------------------------------------------------------------------

/// Polarization vocabulary shared with `_RESULT_KEY_MAP`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pol {
    S,
    P,
    U,
}

/// Spectral channel.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Channel {
    R,
    T,
}

/// Identifies one simulated curve (mirrors `_RESULT_KEY_MAP`), one derived
/// absorptance demand, or one back-incidence demand.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CurveId {
    Rs,
    Rp,
    Ru,
    Ts,
    Tp,
    Tu,
    /// Absorptance demand (s): derived as A = 1 − Rs − Ts, never stored.
    As,
    /// Absorptance demand (p): derived as A = 1 − Rp − Tp, never stored.
    Ap,
    /// Absorptance demand (unpolarized): A = 1 − Ru − Tu, never stored.
    Au,
    /// Back-reflectance (s): back-incidence experiment, `back` rows.
    RBs,
    /// Back-reflectance (p).
    RBp,
    /// Back-reflectance (unpolarized).
    RBu,
    /// Back-transmittance (s).
    TBs,
    /// Back-transmittance (p).
    TBp,
    /// Back-transmittance (unpolarized).
    TBu,
    /// Back-absorptance (s): A = 1 − RBs − TBs, never stored.
    ABs,
    /// Back-absorptance (p).
    ABp,
    /// Back-absorptance (unpolarized).
    ABu,
}

impl CurveId {
    pub const ALL: [CurveId; 18] = [
        CurveId::Rs,
        CurveId::Rp,
        CurveId::Ru,
        CurveId::Ts,
        CurveId::Tp,
        CurveId::Tu,
        CurveId::As,
        CurveId::Ap,
        CurveId::Au,
        CurveId::RBs,
        CurveId::RBp,
        CurveId::RBu,
        CurveId::TBs,
        CurveId::TBp,
        CurveId::TBu,
        CurveId::ABs,
        CurveId::ABp,
        CurveId::ABu,
    ];

    pub fn new(pol: Pol, channel: Channel) -> Self {
        match (pol, channel) {
            (Pol::S, Channel::R) => CurveId::Rs,
            (Pol::P, Channel::R) => CurveId::Rp,
            (Pol::U, Channel::R) => CurveId::Ru,
            (Pol::S, Channel::T) => CurveId::Ts,
            (Pol::P, Channel::T) => CurveId::Tp,
            (Pol::U, Channel::T) => CurveId::Tu,
        }
    }

    /// True for the derived absorptance demands (`As`/`Ap`/`Au`/`ABs`/`ABp`/`ABu`).
    pub fn is_absorption(self) -> bool {
        matches!(
            self,
            CurveId::As | CurveId::Ap | CurveId::Au | CurveId::ABs | CurveId::ABp | CurveId::ABu
        )
    }

    /// True for back-incidence demands (`RB*`/`TB*`/`AB*`).
    pub fn is_back(self) -> bool {
        matches!(
            self,
            CurveId::RBs
                | CurveId::RBp
                | CurveId::RBu
                | CurveId::TBs
                | CurveId::TBp
                | CurveId::TBu
                | CurveId::ABs
                | CurveId::ABp
                | CurveId::ABu
        )
    }

    /// Companion intensity curves an absorption demand derives from.
    /// `None` for plain simulated curves.
    pub fn absorption_companions(self) -> Option<(CurveId, CurveId)> {
        match self {
            CurveId::As => Some((CurveId::Rs, CurveId::Ts)),
            CurveId::Ap => Some((CurveId::Rp, CurveId::Tp)),
            CurveId::Au => Some((CurveId::Ru, CurveId::Tu)),
            CurveId::ABs => Some((CurveId::RBs, CurveId::TBs)),
            CurveId::ABp => Some((CurveId::RBp, CurveId::TBp)),
            CurveId::ABu => Some((CurveId::RBu, CurveId::TBu)),
            _ => None,
        }
    }

    /// Index into `SimCurves::curves`
    /// ([Rs, Rp, Ru, Ts, Tp, Tu, As, Ap, Au]; absorption slots stay `None`).
    pub fn index(self) -> usize {
        match self {
            CurveId::Rs => 0,
            CurveId::Rp => 1,
            CurveId::Ru => 2,
            CurveId::Ts => 3,
            CurveId::Tp => 4,
            CurveId::Tu => 5,
            CurveId::As => 6,
            CurveId::Ap => 7,
            CurveId::Au => 8,
            _ => panic!("back-incidence curves live in SimCurves::back (see back_index)"),
        }
    }

    /// Index into `SimCurves::back` ([RBs, RBp, RBu, TBs, TBp, TBu]).
    /// `None` for front demands (see `index`) and absorption demands.
    pub fn back_index(self) -> Option<usize> {
        match self {
            CurveId::RBs => Some(0),
            CurveId::RBp => Some(1),
            CurveId::RBu => Some(2),
            CurveId::TBs => Some(3),
            CurveId::TBp => Some(4),
            CurveId::TBu => Some(5),
            _ => None,
        }
    }

    /// S-matrix element a phase demand on this curve tracks: front R → 0,
    /// front T → 2, back R → 3, back T → 1. `None` for absorption demands
    /// (phase of absorption is meaningless) and unpolarized keys
    /// (phase of averaged intensities is ill-defined).
    pub fn phase_channel(self) -> Option<usize> {
        match self {
            CurveId::Rs | CurveId::Rp => Some(0),
            CurveId::Ts | CurveId::Tp => Some(2),
            CurveId::RBs | CurveId::RBp => Some(3),
            CurveId::TBs | CurveId::TBp => Some(1),
            _ => None,
        }
    }

    /// Parse curve codes: front `Rs/Rp/Ru/Ts/Tp/Tu`, absorption
    /// `As/Ap/Au`, back `RBs/RBp/RBu/TBs/TBp/TBu` and `ABs/ABp/ABu`.
    /// None for anything else.
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "Rs" => Some(CurveId::Rs),
            "Rp" => Some(CurveId::Rp),
            "Ru" => Some(CurveId::Ru),
            "Ts" => Some(CurveId::Ts),
            "Tp" => Some(CurveId::Tp),
            "Tu" => Some(CurveId::Tu),
            "As" => Some(CurveId::As),
            "Ap" => Some(CurveId::Ap),
            "Au" => Some(CurveId::Au),
            "RBs" => Some(CurveId::RBs),
            "RBp" => Some(CurveId::RBp),
            "RBu" => Some(CurveId::RBu),
            "TBs" => Some(CurveId::TBs),
            "TBp" => Some(CurveId::TBp),
            "TBu" => Some(CurveId::TBu),
            "ABs" => Some(CurveId::ABs),
            "ABp" => Some(CurveId::ABp),
            "ABu" => Some(CurveId::ABu),
            _ => None,
        }
    }
}

/// Constraint activation — mirrors spectralweave `TargetKind`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ConstraintKind {
    /// Always active.
    Exact,
    /// Active only while sim is BELOW target (driving it up).
    Above,
    /// Active only while sim is ABOVE target (driving it down).
    Below,
    /// Hard box of half-width `band`: zero inside, exceedance outside
    /// (paired `Above`/`Below` at centre∓band; bare band falls back to tol).
    Range,
    /// Soft box of half-width `band`: reduced `(d/band)` inside,
    /// exceedance plus continuity level outside (bare band falls back to exact).
    /// NOTE for the needle fold: the `+1` level is gradient-free and dropped
    /// there by design (see `docs/spectralweave-target-kinds.md`).
    CenterBand,
}

impl ConstraintKind {
    /// Parse `"e"/"a"/"b"/"r"/"c"` (spectralweave `TargetKind` codes);
    /// None for anything else.
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "e" => Some(ConstraintKind::Exact),
            "a" => Some(ConstraintKind::Above),
            "b" => Some(ConstraintKind::Below),
            "r" => Some(ConstraintKind::Range),
            "c" => Some(ConstraintKind::CenterBand),
            _ => None,
        }
    }
}

/// Sim-side transform — mirrors spectralweave `ResolvedNormMode`.
///
/// Linear-mode folding note: `(nf·(sim − t)/tol)² ≡ ((sim − t)/(tol/nf))²`,
/// so purely-linear specs may equivalently store Linear + original targets
/// and eff_tol. We keep the spectralweave representation verbatim instead so
/// the parity test compares like-for-like.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SimTransform {
    Linear,
    Log,
    Phase,
    Complex,
}

impl SimTransform {
    /// Parse `"linear"/"log"/"phase"/"complex"`; None for anything else.
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "linear" => Some(SimTransform::Linear),
            "log" => Some(SimTransform::Log),
            "phase" => Some(SimTransform::Phase),
            "complex" => Some(SimTransform::Complex),
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation curves
// ---------------------------------------------------------------------------

/// Simulated R/T curves on the fixed solver grid.
///
/// Each curve is row-major [n_angles, n_wavs]; `angles` / `wavelengths`
/// describe the axes. Arc slices make this cheap to share across rayon
/// workers during Jacobian assembly.
///
/// The indices are **per wavelength** (PD1): length `wavelengths.len()`
/// or 1 (a non-dispersive medium broadcasts, so `Default`'s `[1.0]` keeps
/// `total_d = 0` ≡ absolute phase bit-for-bit). Any other length is a
/// refusal (`reference_length_issue`), naming both numbers.
#[derive(Clone, Debug)]
pub struct SimCurves {
    pub angles: Arc<[f64]>,
    pub wavelengths: Arc<[f64]>,
    pub total_d: f64,
    /// Real incidence index per wavelength (differential-phase reference).
    /// Length `wavelengths.len()` or 1 (broadcast). Default `[1.0]`.
    pub n_front_re: Arc<[f64]>,
    /// Real exit index per wavelength; same shape rule. Default `[1.0]`.
    pub n_back_re: Arc<[f64]>,
    /// Order: [Rs, Rp, Ru, Ts, Tp, Tu, As, Ap, Au]. Absorption slots stay
    /// `None`: absorptance is derived from the companion R/T curves.
    pub curves: [Option<Arc<[f64]>>; 9],
    /// Back-incidence intensity rows: [RBs, RBp, RBu, TBs, TBp, TBu].
    /// Back-absorptance derives from back companions; no slots of its own.
    pub back: [Option<Arc<[f64]>>; 6],
    /// Complex amplitudes for phase demands, front R/T curves only
    /// (same indexing as `curves`; row-major [n_angles, n_wavs]).
    pub cplx: [Option<Arc<[Complex64]>>; 6],
    /// Complex back amplitudes for back-phase demands: [RBs, RBp, TBs, TBp]
    /// (no unpolarized complex rows — back-phase needs s/p keys).
    pub cplx_back: [Option<Arc<[Complex64]>>; 4],
}

impl Default for SimCurves {
    fn default() -> Self {
        Self {
            angles: Arc::from(Vec::new()),
            wavelengths: Arc::from(Vec::new()),
            total_d: 0.0,
            n_front_re: Arc::from([1.0_f64]),
            n_back_re: Arc::from([1.0_f64]),
            curves: Default::default(),
            back: Default::default(),
            cplx: Default::default(),
            cplx_back: Default::default(),
        }
    }
}

impl MeritSpec {
    /// Key-group / target-frame counts (tests + converters).
    pub fn key_count(&self) -> usize {
        self.keys.len()
    }
    pub fn target_count(&self) -> usize {
        self.targets.len()
    }
}

impl SimCurves {
    /// Grid points covered by one row (`n_angles × n_wavs`).
    pub fn grid_points(&self) -> usize {
        self.angles.len() * self.wavelengths.len()
    }

    /// Refusal text for a differential-phase reference row whose length is
    /// neither 1 (broadcast) nor `wavelengths.len()`, naming both numbers.
    /// `None` = both rows well-shaped. Cheap (two compares); wired into the
    /// consumers so a hand-built `SimCurves` fails closed instead of
    /// mis-sampling (ground rule 7 - ASCII).
    pub fn reference_length_issue(&self) -> Option<String> {
        let n_wavs = self.wavelengths.len();
        for (side, row) in [
            ("n_front_re", &self.n_front_re),
            ("n_back_re", &self.n_back_re),
        ] {
            if row.len() != 1 && row.len() != n_wavs {
                return Some(format!(
                    "{side} has length {}, but the grid has {} wavelengths - \
                     supply length 1 (constant medium) or {}",
                    row.len(),
                    n_wavs,
                    n_wavs
                ));
            }
        }
        None
    }

    /// Store one float row (absorptance derives from companions —
    /// nothing to store; lengths are fail-closed).
    pub fn set_curve(&mut self, id: CurveId, values: Arc<[f64]>) -> Result<(), String> {
        if id.is_absorption() {
            return Err(format!(
                "absorption id '{id:?}': absorptance derives from companions, nothing to store"
            ));
        }
        if values.len() != self.grid_points() {
            return Err(format!(
                "curve '{id:?}': {} entries != {} grid points.",
                values.len(),
                self.grid_points()
            ));
        }
        match id.back_index() {
            Some(i) => self.back[i] = Some(values),
            None => self.curves[id.index()] = Some(values),
        }
        Ok(())
    }

    /// Store one complex-amplitude row for phase demands (front s/p R/T
    /// or back s/p keys; lengths are fail-closed).
    pub fn set_complex(&mut self, id: CurveId, values: Arc<[Complex64]>) -> Result<(), String> {
        if values.len() != self.grid_points() {
            return Err(format!(
                "complex curve '{id:?}': {} entries != {} grid points.",
                values.len(),
                self.grid_points()
            ));
        }
        if id.is_back() {
            let i = match id {
                CurveId::RBs => 0,
                CurveId::RBp => 1,
                CurveId::TBs => 2,
                CurveId::TBp => 3,
                _ => {
                    return Err(format!(
                        "complex row for '{id:?}': back-phase needs s/p keys"
                    ));
                }
            };
            self.cplx_back[i] = Some(values);
        } else {
            match id {
                CurveId::Rs | CurveId::Rp | CurveId::Ts | CurveId::Tp => {
                    self.cplx[id.index()] = Some(values)
                }
                _ => {
                    return Err(format!(
                        "complex row for '{id:?}': phase needs s/p R/T keys"
                    ));
                }
            }
        }
        Ok(())
    }

    pub fn curve(&self, id: CurveId) -> Option<&Arc<[f64]>> {
        // Back-incidence rows live in `back` (see back_curve); absorption
        // slots stay None (derived from companions).
        match id.back_index() {
            Some(_) => None,
            None => self.curves.get(id.index()).and_then(|c| c.as_ref()),
        }
    }

    /// Back-incidence intensity row, if supplied.
    pub fn back_curve(&self, id: CurveId) -> Option<&Arc<[f64]>> {
        id.back_index()
            .and_then(|i| self.back.get(i))
            .and_then(|c| c.as_ref())
    }

    /// Complex-amplitude row for phase demands (front R/T s/p curves).
    pub fn complex_curve(&self, id: CurveId) -> Option<&Arc<[Complex64]>> {
        if id.is_absorption() || id.is_back() {
            return None;
        }
        match id {
            CurveId::Ru | CurveId::Tu => None, // unpolarized phase is ill-defined
            _ => self.cplx.get(id.index()).and_then(|c| c.as_ref()),
        }
    }

    /// Complex back-amplitude row for back-phase demands (s/p only).
    pub fn complex_back_curve(&self, id: CurveId) -> Option<&Arc<[Complex64]>> {
        let i = match id {
            CurveId::RBs => 0,
            CurveId::RBp => 1,
            CurveId::TBs => 2,
            CurveId::TBp => 3,
            _ => return None,
        };
        self.cplx_back.get(i).and_then(|c| c.as_ref())
    }

    /// Row index for an angle value: argmin(|angles − a|), first minimum
    /// wins — identical to numpy argmin semantics used by
    /// `_simulate_to_weaver`.
    pub fn angle_row(&self, a: f64) -> usize {
        let mut best = 0usize;
        let mut best_d = f64::INFINITY;
        for (i, &ang) in self.angles.iter().enumerate() {
            let d = (ang - a).abs();
            if d < best_d {
                best_d = d;
                best = i;
            }
        }
        best
    }
}

// ---------------------------------------------------------------------------
// MeritSpec
// ---------------------------------------------------------------------------

/// One residual from a scaled diff (signed forms), shared by the pointwise
/// path and the integral mean path (which calls it once with mean
/// diff/tol/band). Extracted verbatim — formulas bit-identical.
fn kind_residual(kind: ConstraintKind, scaled_diff: f64, tol: f64, bw: f64) -> f64 {
    match kind {
        ConstraintKind::Exact => scaled_diff / tol,
        ConstraintKind::Above if scaled_diff < 0.0 => scaled_diff / tol,
        ConstraintKind::Below if scaled_diff > 0.0 => scaled_diff / tol,
        ConstraintKind::Range => {
            // Bare `r` without a band falls back to the tolerance
            // as half-width (paired a/b at centre∓tol).
            let bw_eff = if bw <= 0.0 { tol } else { bw };
            let ad = scaled_diff.abs();
            if ad <= bw_eff {
                0.0
            } else {
                (ad - bw_eff) / tol
            }
        }
        ConstraintKind::CenterBand => {
            if bw <= 0.0 {
                scaled_diff / tol
            } else {
                let ad = scaled_diff.abs();
                if ad <= bw {
                    scaled_diff / bw
                } else {
                    (((ad - bw) / tol).powi(2) + 1.0).sqrt()
                }
            }
        }
        _ => 0.0,
    }
}

/// One unique (angle, curve) demand — mirrors an `OpticalKey`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MeritKey {
    pub angle: f64,
    pub curve: CurveId,
}

/// One target frame flattened onto its own wavelength grid — mirrors a
/// `(frame uid, OpticalKey) → TargetEntry` pair in spectralweave.
#[derive(Clone, Debug)]
pub struct MeritTarget {
    /// Index into `MeritSpec::keys`.
    pub key_idx: u32,
    /// F2.1: index into the compiled environment list - which of the K
    /// solves this demand reads. `0` for every single-environment run,
    /// which is every run until F2.2 wires the driver loop.
    pub env_idx: u32,
    /// This entry's target wavelength grid (may be any sub/superset of the
    /// solver grid; two-pointer interpolation handles the rest).
    pub wavelengths: Arc<[f64]>,
    pub kind: ConstraintKind,
    pub transform: SimTransform,
    pub norm_factor: f64,
    /// Post-normalization targets from TargetWeaver.
    pub normalized_targets: Arc<[f64]>,
    /// Floored tolerances copied verbatim from the TargetEntry.
    pub tolerances: Arc<[f64]>,
    /// Scaled band half-widths for `Range`/`CenterBand`, copied verbatim
    /// from the TargetEntry (empty means all-zero = unused).
    pub band: Arc<[f64]>,
    /// Phase demand: sample `arg()` of the complex row for the key curve's
    /// element (see `CurveId::phase_channel`) instead of intensities.
    /// Absorption keys with `phase` are rejected at registration.
    pub phase: bool,
    /// Differential phase (`PDts`/`PDtp`): subtract `passes × reference_phase`
    /// (equivalent incidence-medium layer of `SimCurves::total_d`) from the
    /// sampled `arg()` before wrapping. `None` = absolute phase. Requires
    /// `phase` (rejected otherwise); `passes` is 1 for transmitted, 2 for a
    /// reflection round trip; must be finite and non-negative.
    pub differential_passes: Option<f64>,
    /// User weight (default 1): multiplies this frame's merit sum.
    /// Applied at the residual level (`r × √(weight/count)`) so values,
    /// LM Jacobians and the needle fold stay consistent by construction.
    pub weight: f64,
    /// Count-normalization divisor (default None = off): the frame's sum
    /// is divided by this (target-level point count — the converter
    /// resolves it, since angular targets span many single-point frames).
    /// `Some(n)` needs `n > 0` (rejected otherwise).
    pub count_norm: Option<f64>,
    /// Integral target: constrain the MEAN of the scaled diffs (single
    /// residual `R = mean(d)/mean(tol)`), not each point. Kinds apply once
    /// to the mean. Rejected with `count_norm` (the mean already is one).
    pub integral: bool,
}

/// Flat, immutable, Send+Sync target description.
#[derive(Clone, Debug)]
pub struct MeritSpec {
    keys: Vec<MeritKey>,
    targets: Vec<MeritTarget>,
    color: Vec<ColorDemand>,
    /// F2.1: how many environments this spec was compiled against - the
    /// K of "K assemblies, K solves". One unless a segmented request
    /// said otherwise, and `Default` says one rather than zero because
    /// every spec built by hand (the tests, the Rust consumers) is a
    /// single-environment spec.
    n_envs: usize,
    /// M7 (review PB): the environment NAMES this spec was compiled
    /// against, in roster order — when the request named them.
    ///
    /// After `compile_merit_spec` resolves every `environment=` tag to
    /// an index the roster is discarded, so the only thing left binding
    /// a spec to a design was `n_envs`: a COUNT. Two rosters with the
    /// same names in a different order therefore matched, and every
    /// demand scored against the wrong surroundings with no refusal
    /// anywhere. Keeping the names lets the run door check the binding
    /// it actually depends on.
    ///
    /// Empty means "not recorded", not "no environments": a spec built
    /// by hand and every spec whose `TargetSet` named nothing leaves it
    /// empty and keeps the count check alone. Recording names only when
    /// the caller wrote some is what makes this additive — a K=1 request
    /// whose demands carry no tag must keep resolving to environment 0
    /// whatever that environment is called.
    env_names: Vec<String>,
}

impl Default for MeritSpec {
    fn default() -> Self {
        Self {
            keys: Vec::new(),
            targets: Vec::new(),
            color: Vec::new(),
            n_envs: 1,
            env_names: Vec::new(),
        }
    }
}

impl MeritSpec {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of environments this spec routes across.
    pub fn n_envs(&self) -> usize {
        self.n_envs
    }

    /// The roster this spec was compiled against, or empty if the
    /// request named no environments (M7).
    ///
    /// Empty is not an error and not a claim of zero environments: it
    /// means the binding can only be checked by count, which is what
    /// every caller before F2.4 got.
    pub fn env_names(&self) -> &[String] {
        &self.env_names
    }

    /// Record the roster, and with it the environment count (M7).
    ///
    /// Refuses an empty roster, an empty name and a duplicate: a name
    /// that repeats would make a position comparison pass on a roster
    /// whose own `resolve_env` is ambiguous, which is the quiet
    /// wrong-surroundings result this exists to prevent.
    /// `compile_merit_spec` already refuses both upstream — refused here
    /// too because this is the door a Rust consumer reaches.
    pub fn set_env_roster(&mut self, names: &[String]) -> Result<(), String> {
        if names.is_empty() {
            return Err("MeritSpec: the environment roster must be non-empty \
                        (leave it unset to keep the count check)"
                .to_string());
        }
        for (i, name) in names.iter().enumerate() {
            if name.is_empty() {
                return Err(format!("MeritSpec: roster[{i}] is empty"));
            }
            if names[..i].contains(name) {
                return Err(format!("MeritSpec: duplicate roster name {name:?}"));
            }
        }
        self.set_n_envs(names.len())?;
        self.env_names = names.to_vec();
        Ok(())
    }

    /// Declare the environment count, before any demand is added.
    ///
    /// Refuses zero (there is always at least the default environment)
    /// and refuses to shrink below a demand already registered - an
    /// `env_idx` pointing past the end of the solve list is the one way
    /// this bookkeeping could become a panic in F2.2's driver, so it is
    /// made unrepresentable here instead.
    pub fn set_n_envs(&mut self, n: usize) -> Result<(), String> {
        if n == 0 {
            return Err("MeritSpec: n_envs must be >= 1".to_string());
        }
        // M7: a recorded roster IS the count. Letting the two disagree
        // would leave the run door checking a name list of one length
        // against a design of another, and reporting whichever mismatch
        // it happened to test first.
        if !self.env_names.is_empty() && self.env_names.len() != n {
            return Err(format!(
                "MeritSpec: n_envs={n} contradicts the recorded roster of {} \
                 ({}) - set the roster instead",
                self.env_names.len(),
                self.env_names.join(", ")
            ));
        }
        let used = self
            .targets
            .iter()
            .map(|t| t.env_idx)
            .chain(self.color.iter().map(|d| d.env_idx))
            .max();
        if let Some(hi) = used
            && hi as usize >= n
        {
            return Err(format!(
                "MeritSpec: n_envs={n} would orphan a demand already \
                 registered against environment {hi}"
            ));
        }
        self.n_envs = n;
        Ok(())
    }

    /// Register a key group; returns its index.
    pub fn add_key(&mut self, key: MeritKey) -> usize {
        self.keys.push(key);
        self.keys.len() - 1
    }

    /// Append one target frame to key group `key_idx`.
    ///
    /// Wavelengths/targets/tolerances must have equal length; `band` must
    /// either be empty (all-zero = unused) or match that length.
    pub fn add_target(&mut self, target: MeritTarget) -> Result<(), String> {
        if target.env_idx as usize >= self.n_envs {
            return Err(format!(
                "MeritTarget env_idx={} but the spec has {} environment(s)",
                target.env_idx, self.n_envs
            ));
        }
        let n = target.wavelengths.len();
        if target.normalized_targets.len() != n || target.tolerances.len() != n {
            return Err(format!(
                "MeritTarget length mismatch: wl={} targets={} tols={}",
                n,
                target.normalized_targets.len(),
                target.tolerances.len()
            ));
        }
        if !target.band.is_empty() && target.band.len() != n {
            return Err(format!(
                "MeritTarget band mismatch: wl={} band={}",
                n,
                target.band.len()
            ));
        }
        if target.key_idx as usize >= self.keys.len() {
            return Err(format!(
                "key_idx {} out of range ({} keys registered)",
                target.key_idx,
                self.keys.len()
            ));
        }
        if target.phase
            && self.keys[target.key_idx as usize]
                .curve
                .phase_channel()
                .is_none()
        {
            return Err(format!(
                "phase demand on {:?}: absorption/unpolarized curves have no phase",
                self.keys[target.key_idx as usize].curve
            ));
        }
        // Mirror invariant (see spectralweave `register_metadata`): the phase
        // arm scales nothing, so phase demands must carry raw values with
        // norm_factor == 1 (converters: divide the resolved triple by nf).
        if target.transform == SimTransform::Phase && (target.norm_factor - 1.0).abs() > 1e-12 {
            return Err(format!(
                "phase transform needs norm_factor == 1 (got {}); pass raw values",
                target.norm_factor
            ));
        }
        // Differential phase is phase-only (subtracting a propagation
        // reference from an intensity is meaningless) with a sane pass count.
        if let Some(passes) = target.differential_passes {
            if !target.phase {
                return Err("differential_passes without phase: PD demands are phase-only".into());
            }
            if !(passes >= 0.0) || !passes.is_finite() {
                return Err(format!(
                    "differential passes must be finite and >= 0 (got {passes})"
                ));
            }
        }
        // Weight/count trust boundary (bindings + converter pass user
        // values straight through; garbage here means NaN merits).
        if !target.weight.is_finite() || target.weight < 0.0 {
            return Err(format!(
                "weight must be finite and >= 0 (got {})",
                target.weight
            ));
        }
        if let Some(n) = target.count_norm
            && (!n.is_finite() || n <= 0.0)
        {
            return Err(format!("count_norm must be finite and > 0 (got {n})"));
        }
        // Integral targets already are means — a count divisor would
        // double-dilute silently.
        if target.integral && target.count_norm.is_some() {
            return Err("integral targets reject count_norm (the mean already is one)".into());
        }
        self.targets.push(target);
        Ok(())
    }

    /// True if any demand samples complex amplitudes (phase or
    /// differential-phase). The evaluator uses this to skip complex-row
    /// assembly for intensity-only specs — virtually-free values are always
    /// built, requested ones only on demand.
    pub fn uses_phase(&self) -> bool {
        self.targets.iter().any(|t| t.phase)
    }

    /// True if any demand subtracts the equivalent-medium reference.
    /// Gates the (trivial but non-zero) stack-metadata computation.
    pub fn uses_differential(&self) -> bool {
        self.targets.iter().any(|t| t.differential_passes.is_some())
    }

    /// Which reference sides the spec's differential demands actually
    /// read, as `(front, back)`. `n_back_re` is consumed only under
    /// `key.curve.is_back()` (merit.rs residuals, needle_pass gain
    /// shifts), and no differential label maps to a back curve yet (PD
    /// plan §7 decision 6: `differential()` maps only `PDts`/`PDtp`, both
    /// front) — so the back side can only be wrong once a label maps to
    /// it, and the scalar guard must not claim otherwise (review PD G1).
    pub fn demanded_reference_sides(&self) -> (bool, bool) {
        let mut front = false;
        let mut back = false;
        for t in &self.targets {
            if t.differential_passes.is_some() {
                if self.keys[t.key_idx as usize].curve.is_back() {
                    back = true;
                } else {
                    front = true;
                }
            }
        }
        (front, back)
    }

    pub fn keys(&self) -> &[MeritKey] {
        &self.keys
    }

    pub fn targets(&self) -> &[MeritTarget] {
        &self.targets
    }

    /// Append one compiled color demand (1 residual). `key_idx` range is
    /// checked like `add_target` — the demand groups with its key for
    /// missing-penalty + residual ordering.
    pub fn add_color_demand(&mut self, demand: ColorDemand) -> Result<(), String> {
        if demand.env_idx as usize >= self.n_envs {
            return Err(format!(
                "ColorDemand env_idx={} but the spec has {} environment(s)",
                demand.env_idx, self.n_envs
            ));
        }
        if demand.key_idx as usize >= self.keys.len() {
            return Err(format!(
                "key_idx {} out of range ({} keys registered)",
                demand.key_idx,
                self.keys.len()
            ));
        }
        self.color.push(demand);
        Ok(())
    }

    pub fn color_demands(&self) -> &[ColorDemand] {
        &self.color
    }

    /// Number of residual components `residuals()` produces: one per target
    /// point for pointwise frames, ONE per integral frame (the kind applies
    /// to the mean, not to each point), plus one per color demand. Inactive
    /// constraints still occupy their rows — they contribute zeros, so the
    /// length does not depend on the values.
    ///
    /// It does depend on the *grids*: a frame whose wavelengths do not
    /// overlap the simulated grid is skipped entirely by `residuals_into`
    /// and contributes nothing, which cannot be known from the spec alone.
    /// This is therefore the count for a simulation that covers every frame.
    pub fn n_residuals(&self) -> usize {
        self.targets
            .iter()
            .map(|t| if t.integral { 1 } else { t.wavelengths.len() })
            .sum::<usize>()
            + self.color.len()
    }

    /// Scalar merit function: Σ residual² + missing-key penalties.
    ///
    /// Semantics identical to `calculate_merit(sim, tw, missing_penalty)`:
    /// a missing curve costs `missing_penalty` ONCE per key group, and
    /// target grids that do not overlap the simulated grid are skipped
    /// silently (zero contribution).
    ///
    /// # Panics
    ///
    /// Panics on a malformed reference-row length (`n_front_re`/`n_back_re`
    /// neither 1 nor the grid length) — same contract as `residuals`,
    /// unreachable from Python (the FFI constructor refuses), fail-closed
    /// for the hand-built Rust caller (review PD G4).
    pub fn merit(&self, sim: &SimCurves, missing_penalty: f64) -> f64 {
        self.merit_multi(std::slice::from_ref(sim), missing_penalty)
    }

    /// F2.2: the joint merit over K environments — Σ_env Σ_key residual²,
    /// with the missing-curve penalty charged per `(env, key)` pair.
    ///
    /// `sims[e]` is environment `e`'s solve; a demand is scored against the
    /// curves of the environment it was tagged with at compile time, and a
    /// `(env, key)` pair no demand belongs to is skipped without charge —
    /// it has no curve to miss.
    ///
    /// At K = 1 this IS the pre-F2.2 body: one environment, the same key
    /// loop, the same accumulation order (multi-env plan §4.6's single
    /// branch, taken here by the slice length).
    ///
    /// # Panics
    ///
    /// On a malformed reference-row length (as `merit`), and on
    /// `sims.len() != n_envs()`. The second is the same class of
    /// construction contract: scoring a two-environment spec against one
    /// solve would silently drop every demand tagged with the other
    /// environment, which is exactly the quiet wrong-surroundings result
    /// the F2.1 compile refusals exist to prevent.
    pub fn merit_multi(&self, sims: &[SimCurves], missing_penalty: f64) -> f64 {
        self.check_sims(sims, "merit_multi");
        let mut total = 0.0;
        let mut buf: Vec<f64> = Vec::new();
        for (e, sim) in sims.iter().enumerate() {
            for k in 0..self.keys.len() {
                buf.clear();
                match self.residuals_into(sim, k, e as u32, &mut buf) {
                    Ok(()) => total += buf.iter().map(|r| r * r).sum::<f64>(),
                    Err(_) => total += missing_penalty,
                }
            }
        }
        total
    }

    /// Fixed-length residual vector (zeros where constraints are inactive).
    ///
    /// Returns `Err(missing CurveId)` if a demanded curve was not supplied.
    /// Component order: keys in registration order, targets per key in
    /// insertion order, points along each target grid — deterministic,
    /// which the thickness optimizer relies on.
    ///
    /// # Panics
    ///
    /// Panics on a malformed reference-row length (`n_front_re`/`n_back_re`
    /// neither 1 nor the grid length). The FFI constructor refuses that at
    /// build time, so this is unreachable from Python; a Rust caller
    /// hand-building `SimCurves` — the caller this exists for — should
    /// expect the panic rather than a silent mis-sample (review PD G4).
    pub fn residuals(&self, sim: &SimCurves, out: &mut Vec<f64>) -> Result<(), CurveId> {
        self.residuals_multi(std::slice::from_ref(sim), out)
    }

    /// F2.2: the joint residual vector over K environments.
    ///
    /// Component order is **environment-major, insertion-minor**: every
    /// residual of environment 0 in `residuals()` order, then every
    /// residual of environment 1, and so on. `n_residuals()` is already the
    /// total across environments (each target carries its own `env_idx`),
    /// so it is the reservation and the length here too.
    ///
    /// A missing curve is a hard error, as in the single-environment door —
    /// the penalty substitution is `merit_multi`'s, not the LM's.
    ///
    /// # Panics
    ///
    /// As `merit_multi`: malformed reference rows, or `sims.len()` not
    /// equal to `n_envs()`.
    pub fn residuals_multi(&self, sims: &[SimCurves], out: &mut Vec<f64>) -> Result<(), CurveId> {
        self.check_sims(sims, "residuals_multi");
        out.clear();
        out.reserve(self.n_residuals());
        for (e, sim) in sims.iter().enumerate() {
            for k in 0..self.keys.len() {
                self.residuals_into(sim, k, e as u32, out)?;
            }
        }
        Ok(())
    }

    // -- internals -----------------------------------------------------------

    /// Shared precondition of both multi doors: one solve per registered
    /// environment, each with well-formed reference rows.
    fn check_sims(&self, sims: &[SimCurves], who: &str) {
        if sims.len() != self.n_envs {
            panic!(
                "MeritSpec::{who}: {} simulation(s) for {} environment(s) - \
                 every environment must be solved, or demands tagged with \
                 the missing one would score against nothing",
                sims.len(),
                self.n_envs
            );
        }
        for sim in sims {
            // Reference-row shape is a construction contract (PD1): a bad
            // length would mis-sample silently, so fail closed here too.
            // Unreachable from the engine fill and the FFI doors; a
            // hand-built struct is the only way in.
            if let Some(msg) = sim.reference_length_issue() {
                panic!("{msg}");
            }
        }
    }

    /// Residuals for every target belonging to one `(environment, key)`
    /// group.
    ///
    /// Inner loop is a verbatim lift of the calculate_merit body:
    /// overlap skip → aligned fast path / two-pointer interpolation →
    /// mode transform → kind activation. Absorption demands (`As`/`Ap`)
    /// derive A = 1 − R − T from the companion curves on the shared grid;
    /// a missing companion fails the whole key group like a missing curve.
    fn residuals_into(
        &self,
        sim: &SimCurves,
        key_idx: usize,
        env: u32,
        out: &mut Vec<f64>,
    ) -> Result<(), CurveId> {
        let key = &self.keys[key_idx];
        let ang_row = sim.angle_row(key.angle);
        let sim_wl: &[f64] = &sim.wavelengths;
        let n_wav = sim_wl.len();
        // One intensity row, either side (front `curves` or `back`).
        let irow = |id: CurveId| -> Result<&[f64], CurveId> {
            let arc = if id.is_back() {
                sim.back_curve(id)
            } else {
                sim.curve(id)
            };
            arc.map(|c| &c[ang_row * n_wav..(ang_row + 1) * n_wav])
                .ok_or(key.curve)
        };
        // F2.2: `(env, key)` is the group, not `key`. The environment
        // predicate is the ONLY change to this body - no new formula, no
        // new activation kind - and at K = 1 every target carries
        // `env_idx == 0`, so it never rejects anything.
        for t in self
            .targets
            .iter()
            .filter(|t| t.key_idx as usize == key_idx && t.env_idx == env)
        {
            // Resolve this target's simulated input BEFORE pushing anything,
            // so missing rows leave `out` untouched. Phase demands sample
            // arg() of the complex row for the key's element.
            enum TargetInput<'a> {
                Intensity(&'a [f64]),
                Absorption(&'a [f64], &'a [f64]),
                Phase(&'a [Complex64]),
            }
            let input: TargetInput = if t.phase {
                let crow = if key.curve.is_back() {
                    sim.complex_back_curve(key.curve)
                } else {
                    sim.complex_curve(key.curve)
                };
                match crow {
                    Some(c) => TargetInput::Phase(&c[ang_row * n_wav..(ang_row + 1) * n_wav]),
                    None => return Err(key.curve),
                }
            } else if let Some((rc, tc)) = key.curve.absorption_companions() {
                TargetInput::Absorption(irow(rc)?, irow(tc)?)
            } else {
                TargetInput::Intensity(irow(key.curve)?)
            };
            let t_wl: &[f64] = &t.wavelengths;
            if t_wl.is_empty() {
                continue;
            }
            // Skip frames whose grid does not overlap the simulated curve.
            if sim_wl.last().is_none_or(|&l| l < t_wl[0])
                || sim_wl.first().zip(t_wl.last()).is_none_or(|(&f, &l)| f > l)
            {
                continue;
            }

            // Fast path: when the target grid coincides bit-for-bit with a
            // contiguous block of the simulated grid, read simulated values
            // directly — no interpolation, no per-point division.
            let offset = sim_wl.partition_point(|&x| x < t_wl[0]);
            let aligned = offset + t_wl.len() <= n_wav
                && t_wl
                    .iter()
                    .zip(&sim_wl[offset..offset + t_wl.len()])
                    .all(|(&a, &b)| a.to_bits() == b.to_bits());

            // Two-pointer state advances monotonically across the sorted
            // target grid (O(n + m), never reset inside this entry). The
            // sampler is shared by both companion rows (identical grids),
            // so absorption stays consistent point-for-point.
            let mut sim_idx = 0usize;
            let sample = |row: &[f64], i: usize, sim_idx: &mut usize| -> f64 {
                if aligned {
                    return row[offset + i];
                }
                let target_w = t_wl[i];
                while *sim_idx + 1 < n_wav && sim_wl[*sim_idx + 1] < target_w {
                    *sim_idx += 1;
                }
                if *sim_idx + 1 < n_wav && sim_wl[*sim_idx] <= target_w {
                    let w0 = sim_wl[*sim_idx];
                    let w1 = sim_wl[*sim_idx + 1];
                    let v0 = row[*sim_idx];
                    let v1 = row[*sim_idx + 1];
                    if (w1 - w0).abs() < 1e-14 {
                        v0
                    } else {
                        v0 + (target_w - w0) * (v1 - v0) / (w1 - w0)
                    }
                } else if *sim_idx < n_wav {
                    row[*sim_idx]
                } else {
                    row[n_wav - 1]
                }
            };
            // Complex twin for phase demands (shared two-pointer state is
            // safe: identical grids take identical paths).
            let sample_c = |row: &[Complex64], i: usize, sim_idx: &mut usize| -> Complex64 {
                if aligned {
                    return row[offset + i];
                }
                let target_w = t_wl[i];
                while *sim_idx + 1 < n_wav && sim_wl[*sim_idx + 1] < target_w {
                    *sim_idx += 1;
                }
                if *sim_idx + 1 < n_wav && sim_wl[*sim_idx] <= target_w {
                    let w0 = sim_wl[*sim_idx];
                    let w1 = sim_wl[*sim_idx + 1];
                    let v0 = row[*sim_idx];
                    let v1 = row[*sim_idx + 1];
                    if (w1 - w0).abs() < 1e-14 {
                        v0
                    } else {
                        v0 + (v1 - v0) * ((target_w - w0) / (w1 - w0))
                    }
                } else if *sim_idx < n_wav {
                    row[*sim_idx]
                } else {
                    row[n_wav - 1]
                }
            };
            // Per-wavelength reference index, same bracket as the complex
            // row (PD1): the shared two-pointer advances once per point and
            // the second walk finds its bracket already reached. The aligned
            // path reads directly, so a broadcast (length-1) row and a
            // constant full row both yield the scalar bit-for-bit - the
            // non-dispersive case is arithmetic-identical to pre-PD1.
            let sample_n = |row: &[f64], i: usize, sim_idx: &mut usize| -> f64 {
                if row.len() == 1 {
                    return row[0]; // constant medium: broadcast, no state
                }
                if aligned {
                    return row[offset + i];
                }
                let target_w = t_wl[i];
                while *sim_idx + 1 < n_wav && sim_wl[*sim_idx + 1] < target_w {
                    *sim_idx += 1;
                }
                if *sim_idx + 1 < n_wav && sim_wl[*sim_idx] <= target_w {
                    let w0 = sim_wl[*sim_idx];
                    let w1 = sim_wl[*sim_idx + 1];
                    let v0 = row[*sim_idx];
                    let v1 = row[*sim_idx + 1];
                    if (w1 - w0).abs() < 1e-14 {
                        v0
                    } else {
                        v0 + (target_w - w0) * (v1 - v0) / (w1 - w0)
                    }
                } else if *sim_idx < n_wav {
                    row[*sim_idx]
                } else {
                    row[n_wav - 1]
                }
            };
            // Differential-phase reference for this key (front/back medium
            // + total thickness from the sim metadata; None = absolute).
            // `key.angle` is degrees (converter convention); the reference
            // helper converts internally. The medium lookup stays INSIDE
            // the `if let` below so non-differential keys pay nothing.
            let diff_passes = t.differential_passes;
            // Weight + count normalization at the residual level (once per
            // frame): merit scales by weight/count exactly, and LM
            // Jacobians differentiate the scaled residuals consistently.
            // Defaults (1.0/None) are the identity — legacy paths bit-safe.
            let rscale = (t.weight / t.count_norm.unwrap_or(1.0)).sqrt();
            let mut acc_d = 0.0;
            let mut acc_tol = 0.0;
            let mut acc_bw = 0.0;
            for i in 0..t_wl.len() {
                let sim_raw = match &input {
                    TargetInput::Intensity(row) => sample(row, i, &mut sim_idx),
                    TargetInput::Absorption(r, tt) => {
                        1.0 - sample(r, i, &mut sim_idx) - sample(tt, i, &mut sim_idx)
                    }
                    TargetInput::Phase(crow) => {
                        let mut a = sample_c(crow, i, &mut sim_idx).arg();
                        if let Some(passes) = diff_passes {
                            // Per-λ reference index (PD1): sample the row at
                            // the TARGET wavelength, same bracket/interpolation
                            // class as the complex row above. The front/back
                            // pick stays a branch - non-differential keys pay
                            // nothing (the `if let` guards it).
                            let n_row: &[f64] = if key.curve.is_back() {
                                &sim.n_back_re
                            } else {
                                &sim.n_front_re
                            };
                            a -= crate::smatrix::optics_core::reference_phase(
                                t_wl[i],
                                sample_n(n_row, i, &mut sim_idx),
                                key.angle,
                                sim.total_d,
                                passes,
                            );
                        }
                        a
                    }
                };

                let target_scaled = t.normalized_targets[i];

                let scaled_diff = match t.transform {
                    SimTransform::Phase => {
                        // Wrap the residual to [-pi, pi] without trig — cheaper
                        // than the equivalent sin/cos/atan2 formulation.
                        let diff = sim_raw - target_scaled;
                        diff - std::f64::consts::TAU * (diff / std::f64::consts::TAU).round()
                    }
                    SimTransform::Log => sim_raw.max(1e-12).log10() * t.norm_factor - target_scaled,
                    SimTransform::Linear | SimTransform::Complex => {
                        sim_raw * t.norm_factor - target_scaled
                    }
                };

                let tol = t.tolerances[i];
                let bw = t.band.get(i).copied().unwrap_or(0.0);

                if t.integral {
                    // Mean branch: accumulate raw ingredients; the kind
                    // applies ONCE to the mean after the loop.
                    acc_d += scaled_diff;
                    acc_tol += tol;
                    acc_bw += bw;
                } else {
                    out.push(kind_residual(t.kind, scaled_diff, tol, bw) * rscale);
                }
            }
            if t.integral {
                // Single mean residual R = mean(d)/mean(tol); kinds
                // constrain the MEAN (integral-`a` = lower bound on the
                // average). Weight multiplies (count rejected at intake).
                let n = t_wl.len() as f64;
                let tol_bar = (acc_tol / n).max(1e-300);
                out.push(kind_residual(t.kind, acc_d / n, tol_bar, acc_bw / n) * rscale);
            }
        }
        // Color demands of this key group, in demand insertion order
        // (deterministic; documented). One residual each (see `n_residuals`).
        for d in self
            .color
            .iter()
            .filter(|d| d.key_idx as usize == key_idx && d.env_idx == env)
        {
            out.push(self.color_residual_into(sim, d)?);
        }
        Ok(())
    }

    /// One color residual (√F, weight inside F — applied once, at the
    /// objective level, unlike the pointwise `rscale`).
    ///
    /// Failure policy (documented divergence from pointwise grid-miss
    /// skips): a missing curve fails the whole key group (`Err(key.curve)`
    /// → missing-penalty in `merit()`, hard error in `residuals()`/LM —
    /// parity with pointwise missing-curve handling). An empty table/sim
    /// overlap likewise errors instead of skipping: a color demand that
    /// sees nothing is a spec bug, never a silent zero.
    fn color_residual_into(&self, sim: &SimCurves, d: &ColorDemand) -> Result<f64, CurveId> {
        let key = &self.keys[d.key_idx as usize];
        // Compile refused back curves; the kernel never sees them.
        if key.curve.is_back() {
            return Err(key.curve);
        }
        let ang_row = sim.angle_row(key.angle);
        let n_wav = sim.wavelengths.len();
        let row = sim.curve(key.curve).ok_or(key.curve)?;
        let sim_row = row
            .get(ang_row * n_wav..(ang_row + 1) * n_wav)
            .ok_or(key.curve)?;
        eval_color(d, sim_row, &sim.wavelengths)
            .map(|(r, _)| r)
            .map_err(|_| key.curve)
    }
}

// ---------------------------------------------------------------------------
// Row-wise curve sensitivity (R4.5, work item i)
// ---------------------------------------------------------------------------

/// How one residual row depends on one simulated curve value.
///
/// `curve`/`angle_row`/`wavelength` address a single entry of `SimCurves`
/// (`curve[angle_row * n_wav + wavelength]`); `d_residual` is ∂r/∂(that
/// value). A row with an empty term list does not move when the curves move
/// — an inactive constraint, a target grid that misses the simulation, or a
/// kind sitting in its dead zone.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CurveTerm {
    pub curve: CurveId,
    pub angle_row: usize,
    pub wavelength: usize,
    pub d_residual: f64,
}

/// ∂r/∂(simulated curve) for every residual row, in `residuals()` order.
///
/// The merit half of the analytic Jacobian (R4.5): with ∂(curve)/∂θ from the
/// solver sweep, J[i,k] = Σ_terms d_residual · ∂(curve)/∂θ_k. Splitting it
/// here means this half is verifiable on its own, against a finite difference
/// of `residuals()` in curve space — no solver in the loop.
#[derive(Clone, Debug, Default)]
pub struct MeritSensitivity {
    /// `rows[i]` are the terms of residual `i`. Length == `n_residuals()`.
    pub rows: Vec<Vec<CurveTerm>>,
    /// Rows whose dependence this pass cannot express, in the same indexing.
    /// A caller that finds one among its active demands must fall back to a
    /// finite-difference Jacobian rather than treat the row as constant —
    /// see `is_complete`.
    pub uncovered: Vec<usize>,
}

impl MeritSensitivity {
    /// True when every row's dependence is expressed. `false` means some row
    /// is a phase or color demand, whose chain this pass does not carry.
    pub fn is_complete(&self) -> bool {
        self.uncovered.is_empty()
    }
}

/// Where one target point reads the simulated grid, and with what weights.
enum SampleRef {
    /// Reads `sim[idx]` outright (aligned grids, or a clamped edge).
    Direct(usize),
    /// Reads between `idx` and `idx + 1` with fraction `f`:
    /// value = (1−f)·sim[idx] + f·sim[idx+1].
    Interp(usize, f64),
}

impl SampleRef {
    fn read(&self, row: &[f64]) -> f64 {
        match *self {
            SampleRef::Direct(i) => row[i],
            SampleRef::Interp(i, f) => row[i] + f * (row[i + 1] - row[i]),
        }
    }
}

/// ∂(kind_residual)/∂(scaled_diff) — `kind_residual` differentiated arm by arm.
///
/// The kinks are the constraint, not an approximation: `a`/`b` are exactly
/// zero on their satisfied side and `r` is exactly zero inside its band. At a
/// kink this returns the derivative of the arm `kind_residual` itself takes
/// there, so the two agree on which side of the boundary the point is on.
fn d_kind_residual(kind: ConstraintKind, scaled_diff: f64, tol: f64, bw: f64) -> f64 {
    match kind {
        ConstraintKind::Exact => 1.0 / tol,
        ConstraintKind::Above if scaled_diff < 0.0 => 1.0 / tol,
        ConstraintKind::Below if scaled_diff > 0.0 => 1.0 / tol,
        ConstraintKind::Above | ConstraintKind::Below => 0.0,
        ConstraintKind::Range => {
            let bw_eff = if bw <= 0.0 { tol } else { bw };
            let ad = scaled_diff.abs();
            if ad <= bw_eff {
                0.0
            } else {
                scaled_diff.signum() / tol
            }
        }
        ConstraintKind::CenterBand => {
            if bw <= 0.0 {
                1.0 / tol
            } else {
                let ad = scaled_diff.abs();
                if ad <= bw {
                    1.0 / bw
                } else {
                    // d/dx sqrt(((|x|−bw)/tol)² + 1)
                    let u = (ad - bw) / tol;
                    scaled_diff.signum() * u / (tol * (u * u + 1.0).sqrt())
                }
            }
        }
    }
}

/// ∂(scaled_diff)/∂(sim_raw) for one transform.
///
/// `Phase` wraps by subtracting a multiple of 2π, which is piecewise the
/// identity. `Log` is flat below the 1e-12 clamp — the residual genuinely
/// stops responding there, and so does a finite difference.
fn d_transform(transform: SimTransform, sim_raw: f64, norm_factor: f64) -> f64 {
    match transform {
        SimTransform::Phase => 1.0,
        SimTransform::Log => {
            if sim_raw > 1e-12 {
                norm_factor / (sim_raw * std::f64::consts::LN_10)
            } else {
                0.0
            }
        }
        SimTransform::Linear | SimTransform::Complex => norm_factor,
    }
}

/// Does the target grid coincide bit-for-bit with a contiguous block of the
/// simulated grid, and where does that block start? The same expression
/// `residuals_into` uses — integer and bit comparisons only, no arithmetic
/// that could round differently.
fn alignment(t_wl: &[f64], sim_wl: &[f64]) -> (bool, usize) {
    let offset = sim_wl.partition_point(|&x| x < t_wl[0]);
    let aligned = offset + t_wl.len() <= sim_wl.len()
        && t_wl
            .iter()
            .zip(&sim_wl[offset..offset + t_wl.len()])
            .all(|(&a, &b)| a.to_bits() == b.to_bits());
    (aligned, offset)
}

/// Which simulated points target point `i` reads. Mirrors the two-pointer
/// sampler in `residuals_into`, including its edge clamps; `sim_idx` advances
/// monotonically across the sorted target grid exactly as it does there.
fn resolve_sample(
    t_wl: &[f64],
    sim_wl: &[f64],
    i: usize,
    sim_idx: &mut usize,
    aligned: bool,
    offset: usize,
) -> SampleRef {
    let n_wav = sim_wl.len();
    if aligned {
        return SampleRef::Direct(offset + i);
    }
    let target_w = t_wl[i];
    while *sim_idx + 1 < n_wav && sim_wl[*sim_idx + 1] < target_w {
        *sim_idx += 1;
    }
    if *sim_idx + 1 < n_wav && sim_wl[*sim_idx] <= target_w {
        let w0 = sim_wl[*sim_idx];
        let w1 = sim_wl[*sim_idx + 1];
        if (w1 - w0).abs() < 1e-14 {
            SampleRef::Direct(*sim_idx)
        } else {
            SampleRef::Interp(*sim_idx, (target_w - w0) / (w1 - w0))
        }
    } else if *sim_idx < n_wav {
        SampleRef::Direct(*sim_idx)
    } else {
        SampleRef::Direct(n_wav - 1)
    }
}

impl MeritSpec {
    /// ∂r/∂(simulated curve) for every residual row (R4.5).
    ///
    /// Mirrors `residuals()` exactly in row order and count, so the two index
    /// together. Covered: pointwise and integral targets over intensity and
    /// absorption curves, every constraint kind and every real transform.
    /// **Not** covered, and listed in `uncovered`: phase targets (`arg()` of a
    /// complex row, plus a differential reference that depends on the stack's
    /// total thickness directly rather than through any curve) and color
    /// demands (a spectrum-wide integral through the CIE chain —
    /// `build_needle_targets` already emits *that* derivative, in the
    /// per-solver-point shape the needle pass wants rather than per row).
    ///
    /// This deliberately walks the target grids a second time instead of
    /// being folded into `residuals_into`. That function is the
    /// bit-exactness-critical path — its arithmetic is pinned point for point
    /// against the Python original — and threading a sink through it to save
    /// one traversal would put every residual in the engine at risk to buy a
    /// Jacobian. The two are kept in step by
    /// `sensitivity_is_a_finite_difference_of_residuals`, which differences
    /// `residuals()` itself: any drift between them fails it.
    pub fn curve_sensitivity(&self, sim: &SimCurves) -> Result<MeritSensitivity, CurveId> {
        // No reference-row-length guard here, unlike `merit`/`residuals`,
        // and that is deliberate, not an oversight (review PD G4): the
        // reference is additive and independent of the curve, so it drops
        // out of d(residual)/d(curve) — `curve_sensitivity_into` never
        // reads `n_front_re`/`n_back_re`, so a malformed length cannot
        // hurt this path.
        let mut s = MeritSensitivity {
            rows: Vec::with_capacity(self.n_residuals()),
            uncovered: Vec::new(),
        };
        for k in 0..self.keys.len() {
            self.curve_sensitivity_into(sim, k, &mut s)?;
        }
        // The invariant is row-for-row agreement with `residuals()` itself,
        // not with `n_residuals()`: a frame that misses the simulated grid
        // is skipped by both and counted by neither. Debug builds check it
        // against the real thing, so any future divergence in either pass
        // trips everywhere, not only in the one cross-check test.
        debug_assert!(
            {
                let mut v = Vec::new();
                self.residuals(sim, &mut v).map(|()| v.len()) == Ok(s.rows.len())
            },
            "curve_sensitivity produced {} rows; residuals() disagrees",
            s.rows.len()
        );
        Ok(s)
    }

    fn curve_sensitivity_into(
        &self,
        sim: &SimCurves,
        key_idx: usize,
        out: &mut MeritSensitivity,
    ) -> Result<(), CurveId> {
        let key = &self.keys[key_idx];
        let ang_row = sim.angle_row(key.angle);
        let sim_wl: &[f64] = &sim.wavelengths;
        let n_wav = sim_wl.len();
        let irow = |id: CurveId| -> Result<&[f64], CurveId> {
            let arc = if id.is_back() {
                sim.back_curve(id)
            } else {
                sim.curve(id)
            };
            arc.map(|c| &c[ang_row * n_wav..(ang_row + 1) * n_wav])
                .ok_or(key.curve)
        };

        for t in self
            .targets
            .iter()
            .filter(|t| t.key_idx as usize == key_idx)
        {
            let t_wl: &[f64] = &t.wavelengths;
            let n_rows = if t.integral { 1 } else { t_wl.len() };

            // sim_raw = const + Σ sgn·curve. Absorption is A = 1 − R − T, so
            // the constant is 1 and both companions enter with −1. Resolve
            // the rows BEFORE anything else, so a missing curve errors here
            // exactly where `residuals_into` errors.
            let (konst, channels): (f64, Vec<(CurveId, &[f64], f64)>) = if t.phase {
                // Phase rows end up uncovered, but the *failure* has to
                // match: `residuals_into` errors on a missing complex row
                // before it pushes anything, so a spec that cannot be
                // evaluated at all must not come back as a tidy list of
                // empty rows.
                let have = if key.curve.is_back() {
                    sim.complex_back_curve(key.curve).is_some()
                } else {
                    sim.complex_curve(key.curve).is_some()
                };
                if !have {
                    return Err(key.curve);
                }
                (0.0, Vec::new())
            } else if let Some((rc, tc)) = key.curve.absorption_companions() {
                (1.0, vec![(rc, irow(rc)?, -1.0), (tc, irow(tc)?, -1.0)])
            } else {
                (0.0, vec![(key.curve, irow(key.curve)?, 1.0)])
            };

            if t_wl.is_empty() {
                continue;
            }
            // Grid miss: `residuals_into` pushes nothing for this target, so
            // neither does this.
            if sim_wl.last().is_none_or(|&l| l < t_wl[0])
                || sim_wl.first().zip(t_wl.last()).is_none_or(|(&f, &l)| f > l)
            {
                continue;
            }
            // A phase target occupies its rows like any other; they are just
            // unexplained. Reported as uncovered rather than left empty, so a
            // phase demand is never quietly read as covered-and-constant.
            if t.phase {
                for _ in 0..n_rows {
                    out.uncovered.push(out.rows.len());
                    out.rows.push(Vec::new());
                }
                continue;
            }

            let (aligned, offset) = alignment(t_wl, sim_wl);
            let rscale = (t.weight / t.count_norm.unwrap_or(1.0)).sqrt();

            // Re-derive the values the residual pass saw; every derivative
            // factor below is evaluated at those same points.
            let mut sim_idx = 0usize;
            let mut refs: Vec<SampleRef> = Vec::with_capacity(t_wl.len());
            let mut raws: Vec<f64> = Vec::with_capacity(t_wl.len());
            let mut diffs: Vec<f64> = Vec::with_capacity(t_wl.len());
            for i in 0..t_wl.len() {
                let r = resolve_sample(t_wl, sim_wl, i, &mut sim_idx, aligned, offset);
                let raw = konst
                    + channels
                        .iter()
                        .map(|(_, row, sgn)| sgn * r.read(row))
                        .sum::<f64>();
                let target_scaled = t.normalized_targets[i];
                let diff = match t.transform {
                    SimTransform::Phase => {
                        let d = raw - target_scaled;
                        d - std::f64::consts::TAU * (d / std::f64::consts::TAU).round()
                    }
                    SimTransform::Log => raw.max(1e-12).log10() * t.norm_factor - target_scaled,
                    SimTransform::Linear | SimTransform::Complex => {
                        raw * t.norm_factor - target_scaled
                    }
                };
                refs.push(r);
                raws.push(raw);
                diffs.push(diff);
            }

            let n = t_wl.len() as f64;
            if t.integral {
                // One row over the mean: ∂R/∂value = rscale · kind'(mean)
                //                        · (1/n) · transform'(raw_i) · channel
                //                        · sample weight.
                let mean_d = diffs.iter().sum::<f64>() / n;
                let mean_tol = (t.tolerances.iter().take(t_wl.len()).sum::<f64>() / n).max(1e-300);
                let mean_bw = t.band.iter().take(t_wl.len()).sum::<f64>() / n;
                let dk = d_kind_residual(t.kind, mean_d, mean_tol, mean_bw);
                let mut terms: Vec<CurveTerm> = Vec::new();
                for i in 0..t_wl.len() {
                    let dt = d_transform(t.transform, raws[i], t.norm_factor);
                    push_terms(
                        &mut terms,
                        &channels,
                        ang_row,
                        &refs[i],
                        rscale * dk * dt / n,
                    );
                }
                out.rows.push(terms);
            } else {
                for i in 0..t_wl.len() {
                    let tol = t.tolerances[i];
                    let bw = t.band.get(i).copied().unwrap_or(0.0);
                    let dk = d_kind_residual(t.kind, diffs[i], tol, bw);
                    let dt = d_transform(t.transform, raws[i], t.norm_factor);
                    let mut terms: Vec<CurveTerm> = Vec::new();
                    push_terms(&mut terms, &channels, ang_row, &refs[i], rscale * dk * dt);
                    out.rows.push(terms);
                }
            }
        }

        // One row per color demand, all unexplained here: a color residual
        // integrates the whole spectrum through the CIE chain, and
        // `build_needle_targets` already emits that derivative as
        // `grad_r`/`grad_t`, per solver point rather than per row.
        for _ in self.color.iter().filter(|d| d.key_idx as usize == key_idx) {
            out.uncovered.push(out.rows.len());
            out.rows.push(Vec::new());
        }
        Ok(())
    }
}

/// Append ∂r/∂value for every (channel × simulated point) this sample reads.
/// Terms that are structurally zero are dropped rather than stored: a caller
/// multiplying them by ∂(curve)/∂θ would only be adding zeros.
fn push_terms(
    terms: &mut Vec<CurveTerm>,
    channels: &[(CurveId, &[f64], f64)],
    angle_row: usize,
    r: &SampleRef,
    scale: f64,
) {
    if scale == 0.0 {
        return;
    }
    for &(curve, _, sgn) in channels {
        let mut emit = |wavelength: usize, w: f64| {
            if w != 0.0 {
                terms.push(CurveTerm {
                    curve,
                    angle_row,
                    wavelength,
                    d_residual: w,
                });
            }
        };
        match *r {
            SampleRef::Direct(i) => emit(i, scale * sgn),
            SampleRef::Interp(i, f) => {
                emit(i, scale * sgn * (1.0 - f));
                emit(i + 1, scale * sgn * f);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests — pin the lifted kernel against hand-computed values
// ---------------------------------------------------------------------------

/// Reference-rotation factors for differential-phase demands:
/// `exp(-i·ref)` with `ref = passes·2π·n_inc·total_d·cosθ/λ` per
/// wavelength. `arg(a·factor)` is the differential phase `Δφ`.
///
/// `n_inc` is per-λ (PD1/PD2): length 1 broadcasts (the scalar door),
/// otherwise exactly `wavelengths.len()`; anything else refuses, naming
/// both numbers. The per-λ arithmetic is identical to the pre-PD2 scalar
/// kernel when the row is constant — same `f64` operand, same order.
pub fn reference_rotation(
    wavelengths: &[f64],
    angle_deg: f64,
    n_inc: &[f64],
    total_d: f64,
    passes: f64,
) -> Result<Vec<num_complex::Complex64>, String> {
    use std::f64::consts::PI;
    if n_inc.is_empty() {
        return Err("reference_rotation: n_inc is empty - supply length 1 \
             (constant medium) or one value per wavelength"
            .into());
    }
    if n_inc.len() != 1 && n_inc.len() != wavelengths.len() {
        return Err(format!(
            "reference_rotation: n_inc has length {}, but the grid has {} \
             wavelengths - supply length 1 (constant medium) or {}",
            n_inc.len(),
            wavelengths.len(),
            wavelengths.len()
        ));
    }
    let cos_t = angle_deg.to_radians().cos();
    Ok(wavelengths
        .iter()
        .enumerate()
        .map(|(i, w)| {
            let n_at = if n_inc.len() == 1 { n_inc[0] } else { n_inc[i] };
            let r = passes * 2.0 * PI * n_at * total_d * cos_t / w;
            num_complex::Complex64::new(0.0, -r).exp()
        })
        .collect())
}

/// Apply per-wavelength rotation factors to flat rows in place.
/// `rows.len()` must be a multiple of `rot.len()` (last axis wavelength).
pub fn rotate_rows(
    rows: &mut [num_complex::Complex64],
    rot: &[num_complex::Complex64],
) -> Result<(), String> {
    if rot.is_empty() || !rows.len().is_multiple_of(rot.len()) {
        return Err(format!(
            "rotate_rows: {} entries not a multiple of {} wavelengths.",
            rows.len(),
            rot.len()
        ));
    }
    for (i, v) in rows.iter_mut().enumerate() {
        *v *= rot[i % rot.len()];
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const NW: usize = 5;

    /// Sim grid 400..800 step 100 nm, single angle row, R_s supplied.
    fn sim_one_angle(vals: &[f64; NW]) -> SimCurves {
        SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: [
                Some(Arc::from(vals.to_vec())),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            ..Default::default()
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn entry(
        key_idx: u32,
        wl: Vec<f64>,
        targets: Vec<f64>,
        tols: Vec<f64>,
        kind: ConstraintKind,
        transform: SimTransform,
        norm_factor: f64,
    ) -> MeritTarget {
        entry_banded(
            key_idx,
            wl,
            targets,
            tols,
            vec![],
            kind,
            transform,
            norm_factor,
        )
    }

    fn entry_banded(
        key_idx: u32,
        wl: Vec<f64>,
        targets: Vec<f64>,
        tols: Vec<f64>,
        band: Vec<f64>,
        kind: ConstraintKind,
        transform: SimTransform,
        norm_factor: f64,
    ) -> MeritTarget {
        MeritTarget {
            key_idx,
            env_idx: 0,
            wavelengths: wl.into(),
            kind,
            transform,
            norm_factor,
            normalized_targets: targets.into(),
            tolerances: tols.into(),
            band: band.into(),
            phase: false,
            differential_passes: None,
            integral: false,
            weight: 1.0,
            count_norm: None,
        }
    }

    fn entry_phase(
        key_idx: u32,
        wl: Vec<f64>,
        targets: Vec<f64>,
        tols: Vec<f64>,
        kind: ConstraintKind,
        transform: SimTransform,
        norm_factor: f64,
    ) -> MeritTarget {
        MeritTarget {
            key_idx,
            env_idx: 0,
            wavelengths: wl.into(),
            kind,
            transform,
            norm_factor,
            normalized_targets: targets.into(),
            tolerances: tols.into(),
            band: Vec::new().into(),
            phase: true,
            differential_passes: None,
            integral: false,
            weight: 1.0,
            count_norm: None,
        }
    }

    #[test]
    fn reference_rows_broadcast_and_refuse_by_length() {
        // PD1: the reference rows are per-λ (length nw) or length-1
        // (broadcast); anything else is a refusal naming both numbers.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0, 500.0],
            vec![0.0, 0.0],
            vec![0.05, 0.05],
            ConstraintKind::Exact,
            SimTransform::Phase,
            1.0,
        ))
        .unwrap();
        // Length-1 broadcast and Default both accept (Default is [1.0]).
        let mut sim = sim_one_angle(&[0.5, 1.0, 0.9, 0.8, 0.7]);
        sim.n_front_re = Arc::from([1.33_f64]);
        sim.n_back_re = Arc::from([1.52_f64]);
        assert!(sim.reference_length_issue().is_none());
        let _ = spec.merit(&sim, 1e6); // must not panic
        // Length 2 against a 5-wavelength grid refuses, naming both.
        let mut sim = sim_one_angle(&[0.5, 1.0, 0.9, 0.8, 0.7]);
        sim.n_front_re = Arc::from([1.3_f64, 1.4_f64]);
        assert_eq!(
            sim.reference_length_issue().as_deref(),
            Some(
                "n_front_re has length 2, but the grid has 5 wavelengths - \
                 supply length 1 (constant medium) or 5"
            )
        );
        let merit_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = spec.merit(&sim, 1e6);
        }));
        assert!(
            merit_result.is_err(),
            "bad reference length must fail closed"
        );
    }

    #[test]
    fn reference_rows_default_and_full_length_valid() {
        // Default ([1.0]) and a full nw-length row (the engine fill shape)
        // are both valid; total_d = 0 keeps differential ≡ absolute.
        let sim = SimCurves::default();
        assert!(sim.reference_length_issue().is_none());
        assert_eq!(sim.n_front_re.len(), 1);
        assert_eq!(sim.n_front_re[0], 1.0);
        assert_eq!(sim.n_back_re[0], 1.0);
        let mut sim = sim_one_angle(&[0.5, 1.0, 0.9, 0.8, 0.7]);
        sim.n_front_re = Arc::from(vec![1.0_f64; NW]);
        assert!(sim.reference_length_issue().is_none());
    }

    #[test]
    fn linear_fold_exact_zero() {
        // targets [0.5, 1.0]: avg = 0.75 → nf = 4/3 (register_metadata math)
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        let nf = 1.0 / 0.75;
        spec.add_target(entry(
            k as u32,
            vec![400.0, 500.0],
            vec![0.5 * nf, 1.0 * nf],
            vec![0.01, 0.01],
            ConstraintKind::Exact,
            SimTransform::Linear,
            nf,
        ))
        .unwrap();

        let sim = sim_one_angle(&[0.5, 1.0, 0.9, 0.8, 0.7]);
        assert_eq!(spec.merit(&sim, 1e6), 0.0);
    }

    #[test]
    fn exact_residual_hand_computed() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0],
            vec![0.5], // already-normalized target
            vec![0.1], // tol
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let sim = sim_one_angle(&[0.55, 0., 0., 0., 0.]);
        // r = (0.55 − 0.5)/0.1 = 0.5 → merit 0.25
        assert!((spec.merit(&sim, 1e6) - 0.25).abs() < 1e-14);
    }

    #[test]
    fn above_below_activation() {
        // Above: active only while sim < target; Below: mirror image.
        let mk = |kind| {
            let mut s = MeritSpec::new();
            let k = s.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            s.add_target(entry(
                k as u32,
                vec![400.0],
                vec![0.5],
                vec![0.1],
                kind,
                SimTransform::Linear,
                1.0,
            ))
            .unwrap();
            s
        };
        let above = mk(ConstraintKind::Above);
        let below = mk(ConstraintKind::Below);

        let sim_hi = sim_one_angle(&[0.6, 0., 0., 0., 0.]); // sim > target
        let sim_lo = sim_one_angle(&[0.3, 0., 0., 0., 0.]); // sim < target

        assert_eq!(above.merit(&sim_hi, 0.0), 0.0); // satisfied → inactive
        assert!((above.merit(&sim_lo, 0.0) - 4.0).abs() < 1e-14); // ((0.3−0.5)/0.1)²
        assert_eq!(below.merit(&sim_lo, 0.0), 0.0); // satisfied → inactive
        assert!((below.merit(&sim_hi, 0.0) - 1.0).abs() < 1e-14); // ((0.6−0.5)/0.1)²
    }

    #[test]
    fn residual_vector_zeros_inactive_but_fixed_length() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0, 500.0],
            vec![0.5, 0.5],
            vec![0.1, 0.1],
            ConstraintKind::Above,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let sim = sim_one_angle(&[0.6, 0.3, 0., 0., 0.]);
        let mut out = Vec::new();
        spec.residuals(&sim, &mut out).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], 0.0); // Above satisfied at 400
        assert!((out[1] - (-2.0)).abs() < 1e-14); // active at 500
        assert_eq!(spec.n_residuals(), 2);
    }

    #[test]
    fn log_transform_matches_register_metadata() {
        // Log mode: nf = 1/avg(|log10 t|); scaled diff = log10(sim)·nf − tgt·nf
        let targets = [0.01_f64, 1.0];
        let log_sum: f64 = targets.iter().map(|v| v.max(1e-12).log10().abs()).sum();
        let avg = log_sum / 2.0;
        let nf = 1.0 / avg.max(1e-12);

        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0],
            vec![targets[0].max(1e-12).log10() * nf],
            vec![0.05],
            ConstraintKind::Exact,
            SimTransform::Log,
            nf,
        ))
        .unwrap();
        // sim == target exactly → zero residual
        let sim = sim_one_angle(&[0.01, 0., 0., 0., 0.]);
        assert!(spec.merit(&sim, 0.0).abs() < 1e-14);
        // sim off by ×10 in raw space → log10 diff = 1 · nf / tol
        let sim10 = sim_one_angle(&[0.1, 0., 0., 0., 0.]);
        let expect = (1.0 * nf / 0.05).powi(2);
        assert!((spec.merit(&sim10, 0.0) - expect).abs() < 1e-9 * expect);
    }

    #[test]
    fn phase_wrap_no_trig() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0],
            vec![0.0],
            vec![1.0],
            ConstraintKind::Exact,
            SimTransform::Phase,
            1.0,
        ))
        .unwrap();
        // sim = π + 0.1 → wrapped diff −(π − 0.1)
        let sim = sim_one_angle(&[std::f64::consts::PI + 0.1, 0., 0., 0., 0.]);
        let expect = (std::f64::consts::PI - 0.1).powi(2);
        assert!((spec.merit(&sim, 0.0) - expect).abs() < 1e-12);
    }

    #[test]
    fn misaligned_interpolation_two_pointer() {
        // Coarse target grid over a fine sim ramp — linear interp hits exactly.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k as u32,
            vec![450.0, 550.0],
            vec![0.5, 1.5],
            vec![1.0, 1.0],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let sim = sim_one_angle(&[0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!(spec.merit(&sim, 0.0).abs() < 1e-14);
    }

    #[test]
    fn extrapolation_clamps_and_overlap_skips() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        // Below sim range → overlap holds (sim covers [400,800]) → clamps to first val 0.3.
        spec.add_target(entry(
            k as u32,
            vec![350.0],
            vec![0.3],
            vec![1.0],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        // Above sim range entirely (900 > 800) → skipped by overlap rule,
        // but still occupies its residual slot.
        spec.add_target(entry(
            k as u32,
            vec![900.0],
            vec![0.0],
            vec![1.0],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();

        let sim = sim_one_angle(&[0.3, 0., 0., 0., 0.]);
        assert_eq!(spec.merit(&sim, 0.0), 0.0);
        assert_eq!(spec.n_residuals(), 2);
    }

    #[test]
    fn missing_curve_penalty_once_per_key() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Tp,
        }); // not supplied
        spec.add_target(entry(
            k as u32,
            vec![400.0],
            vec![0.0],
            vec![1.0],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        spec.add_target(entry(
            k as u32,
            vec![500.0],
            vec![0.0],
            vec![1.0],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let k2 = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }); // supplied
        spec.add_target(entry(
            k2 as u32,
            vec![400.0],
            vec![0.0],
            vec![1.0],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();

        let sim = sim_one_angle(&[0.0, 0., 0., 0., 0.]);
        // penalty once for the Tp group + zero from Rs
        assert_eq!(spec.merit(&sim, 123.0), 123.0);
        assert!(matches!(
            spec.residuals(&sim, &mut Vec::new()),
            Err(CurveId::Tp)
        ));
    }

    #[test]
    fn angle_row_argmin_semantics() {
        let mut sim = SimCurves {
            angles: vec![0.0, 30.0, 60.0].into(),
            wavelengths: vec![500.0].into(),
            curves: [None, None, None, None, None, None, None, None, None],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            ..Default::default()
        };
        sim.curves[CurveId::Ru.index()] = Some(Arc::from(vec![10.0, 20.0, 30.0]));

        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 25.0,
            curve: CurveId::Ru,
        }); // → row 30°
        spec.add_target(entry(
            k as u32,
            vec![500.0],
            vec![20.0],
            vec![0.5],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();

        let mut out = Vec::new();
        spec.residuals(&sim, &mut out).unwrap();
        assert_eq!(out.len(), 1);
        assert!(out[0].abs() < 1e-14); // picked the 30° row value 20
    }

    #[test]
    fn aligned_fast_path_agrees_with_interp_path() {
        let vals = [0.31, 0.52, 0.66, 0.71, 0.90];

        // Aligned: target wl exactly on the sim grid (bit-equal).
        let mut spec_a = MeritSpec::new();
        let ka = spec_a.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec_a
            .add_target(entry(
                ka as u32,
                vec![500.0],
                vec![0.42],
                vec![0.07],
                ConstraintKind::Exact,
                SimTransform::Linear,
                1.0,
            ))
            .unwrap();

        let sim = sim_one_angle(&vals);
        let mut out = Vec::new();
        spec_a.residuals(&sim, &mut out).unwrap();
        let r_aligned = out[0];

        // Interp: shift the sim grid by tiny per-point amounts so no value is
        // bit-equal, then request the same physical wavelength.
        let shifted: SimCurves = SimCurves {
            angles: sim.angles.clone(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64 + 1e-11 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: sim.curves.clone(),
            back: sim.back.clone(),
            cplx: sim.cplx.clone(),
            cplx_back: sim.cplx_back.clone(),
            ..Default::default()
        };
        let mut spec_u = MeritSpec::new();
        let ku = spec_u.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec_u
            .add_target(entry(
                ku as u32,
                vec![500.0 + 1e-11],
                vec![0.42],
                vec![0.07],
                ConstraintKind::Exact,
                SimTransform::Linear,
                1.0,
            ))
            .unwrap();

        let mut out2 = Vec::new();
        spec_u.residuals(&shifted, &mut out2).unwrap();
        let r_interp = out2[0];

        assert!((r_aligned - r_interp).abs() < 1e-9);
    }

    #[test]
    fn constraint_kind_from_str() {
        assert_eq!(ConstraintKind::from_str("e"), Some(ConstraintKind::Exact));
        assert_eq!(ConstraintKind::from_str("a"), Some(ConstraintKind::Above));
        assert_eq!(ConstraintKind::from_str("b"), Some(ConstraintKind::Below));
        assert_eq!(ConstraintKind::from_str("r"), Some(ConstraintKind::Range));
        assert_eq!(
            ConstraintKind::from_str("c"),
            Some(ConstraintKind::CenterBand)
        );
        assert_eq!(ConstraintKind::from_str("x"), None);
    }

    #[test]
    fn range_box_dead_band() {
        // target 0.5 (nf=1), tol 0.1, band 0.05.
        let mk = |band: Vec<f64>| {
            let mut s = MeritSpec::new();
            let k = s.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            s.add_target(entry_banded(
                k as u32,
                vec![400.0],
                vec![0.5],
                vec![0.1],
                band,
                ConstraintKind::Range,
                SimTransform::Linear,
                1.0,
            ))
            .unwrap();
            s
        };
        let spec = mk(vec![0.05]);
        assert_eq!(spec.merit(&sim_one_angle(&[0.5, 0., 0., 0., 0.]), 0.0), 0.0);
        assert_eq!(
            spec.merit(&sim_one_angle(&[0.53, 0., 0., 0., 0.]), 0.0),
            0.0
        ); // d=0.03 in band
        // d=0.1 → ((0.1−0.05)/0.1)² = 0.25
        assert!((spec.merit(&sim_one_angle(&[0.6, 0., 0., 0., 0.]), 0.0) - 0.25).abs() < 1e-14);
        // Bare band falls back to tol as half-width: d=0.05 inside → 0.
        let bare = mk(vec![]);
        assert_eq!(
            bare.merit(&sim_one_angle(&[0.55, 0., 0., 0., 0.]), 0.0),
            0.0
        );
        // d=0.2 → ((0.2−0.1)/0.1)² = 1.
        assert!((bare.merit(&sim_one_angle(&[0.7, 0., 0., 0., 0.]), 0.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn centerband_soft_box() {
        // target 0.5 (nf=1), tol 0.1, band 0.05.
        let mk = |band: Vec<f64>| {
            let mut s = MeritSpec::new();
            let k = s.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            s.add_target(entry_banded(
                k as u32,
                vec![400.0],
                vec![0.5],
                vec![0.1],
                band,
                ConstraintKind::CenterBand,
                SimTransform::Linear,
                1.0,
            ))
            .unwrap();
            s
        };
        let spec = mk(vec![0.05]);
        assert_eq!(spec.merit(&sim_one_angle(&[0.5, 0., 0., 0., 0.]), 0.0), 0.0);
        // Inside: (0.02/0.05)² = 0.16.
        assert!((spec.merit(&sim_one_angle(&[0.52, 0., 0., 0., 0.]), 0.0) - 0.16).abs() < 1e-14);
        // At the edge: 1.0 from both sides (continuity).
        assert!((spec.merit(&sim_one_angle(&[0.55, 0., 0., 0., 0.]), 0.0) - 1.0).abs() < 1e-14);
        // Outside: ((0.05)/0.1)² + 1 = 1.25.
        assert!((spec.merit(&sim_one_angle(&[0.6, 0., 0., 0., 0.]), 0.0) - 1.25).abs() < 1e-14);
        // Bare band degrades to exact: (0.1/0.1)² = 1.
        let bare = mk(vec![]);
        assert!((bare.merit(&sim_one_angle(&[0.6, 0., 0., 0., 0.]), 0.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn merit_two_keys_no_double_count() {
        // Regression: merit() reused its scratch buffer across keys,
        // over-counting every key after the first.
        let mut spec = MeritSpec::new();
        let k0 = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k0 as u32,
            vec![400.0],
            vec![0.5],
            vec![0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let k1 = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rp,
        });
        spec.add_target(entry(
            k1 as u32,
            vec![400.0],
            vec![0.5],
            vec![0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let mut sim = sim_one_angle(&[0.6, 0., 0., 0., 0.]);
        sim.curves[CurveId::Rp.index()] = Some(Arc::from(vec![0.6, 0., 0., 0., 0.]));
        // Each key: ((0.6−0.5)/0.1)² = 1 → total 2 (was 3 with stale buffer).
        assert!((spec.merit(&sim, 1e6) - 2.0).abs() < 1e-14);
    }

    #[test]
    fn absorption_derived_from_companions() {
        // R row 0.6, T row 0.3 → A = 0.1; demand A = 0.1, tol 0.05 → 0.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::As,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0],
            vec![0.1],
            vec![0.05],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let mut sim = SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: [None, None, None, None, None, None, None, None, None],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            ..Default::default()
        };
        sim.curves[CurveId::Rs.index()] = Some(Arc::from(vec![0.6, 0., 0., 0., 0.]));
        sim.curves[CurveId::Ts.index()] = Some(Arc::from(vec![0.3, 0., 0., 0., 0.]));
        assert!(spec.merit(&sim, 1e6) < 1e-28); // A = 1−0.6−0.3 ≈ 0.1
        let mut out = Vec::new();
        spec.residuals(&sim, &mut out).unwrap();
        assert_eq!(out.len(), 1);
        assert!(out[0].abs() < 1e-14);
        // Missing companion → penalty once, residuals Err on the demand.
        sim.curves[CurveId::Ts.index()] = None;
        assert_eq!(spec.merit(&sim, 123.0), 123.0);
        assert!(matches!(
            spec.residuals(&sim, &mut Vec::new()),
            Err(CurveId::As)
        ));
        // Unpolarized demand derives from the Ru/Tu companions the same way.
        let mut spec_u = MeritSpec::new();
        let ku = spec_u.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Au,
        });
        spec_u
            .add_target(entry(
                ku as u32,
                vec![400.0],
                vec![0.2],
                vec![0.1],
                ConstraintKind::Exact,
                SimTransform::Linear,
                1.0,
            ))
            .unwrap();
        sim.curves[CurveId::Ru.index()] = Some(Arc::from(vec![0.5, 0., 0., 0., 0.]));
        sim.curves[CurveId::Tu.index()] = Some(Arc::from(vec![0.3, 0., 0., 0., 0.]));
        assert!(spec_u.merit(&sim, 1e6) < 1e-28); // A = 1−0.5−0.3 ≈ 0.2
        sim.curves[CurveId::Ru.index()] = None;
        assert!(matches!(
            spec_u.residuals(&sim, &mut Vec::new()),
            Err(CurveId::Au)
        ));
    }

    #[test]
    fn phase_demand_samples_argument() {
        // Complex row 0.5·e^{i·0.3}; demand phase 0.3 (nf=1) → zero.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry_phase(
            k as u32,
            vec![400.0],
            vec![0.3],
            vec![0.05],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let mut sim = SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: [None, None, None, None, None, None, None, None, None],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            ..Default::default()
        };
        sim.cplx[0] = Some(Arc::from(vec![
            Complex64::from_polar(0.5, 0.3),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]));
        assert!(spec.merit(&sim, 1e6) < 1e-28);
        // Missing complex row → penalty once, Err on the demand.
        sim.cplx[0] = None;
        assert_eq!(spec.merit(&sim, 123.0), 123.0);
        assert!(matches!(
            spec.residuals(&sim, &mut Vec::new()),
            Err(CurveId::Rs)
        ));
    }

    #[test]
    fn phase_demand_wraps_in_phase_mode() {
        // Sim phase 0.3 + 2π − 0.01 vs target 0.3: wrapped diff −0.01.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec.add_target(entry_phase(
            k as u32,
            vec![400.0],
            vec![0.3],
            vec![0.05],
            ConstraintKind::Exact,
            SimTransform::Phase,
            1.0,
        ))
        .unwrap();
        let mut sim = SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: [None, None, None, None, None, None, None, None, None],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            ..Default::default()
        };
        sim.cplx[3] = Some(Arc::from(vec![
            Complex64::from_polar(0.7, 0.3 - 0.01 + std::f64::consts::TAU),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]));
        // (−0.01/0.05)² = 0.04
        assert!((spec.merit(&sim, 1e6) - 0.04).abs() < 1e-12);
    }

    #[test]
    fn phase_on_absorption_rejected() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::As,
        });
        let mut tgt = entry_phase(
            k as u32,
            vec![400.0],
            vec![0.0],
            vec![0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        tgt.key_idx = k as u32;
        assert!(spec.add_target(tgt).is_err());
        let k2 = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ru,
        });
        let tgt2 = entry_phase(
            k2 as u32,
            vec![400.0],
            vec![0.0],
            vec![0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        assert!(spec.add_target(tgt2).is_err());
    }

    /// Differential-phase entry: absolute-phase entry + `passes`.
    fn entry_pd(
        key_idx: u32,
        wl: Vec<f64>,
        targets: Vec<f64>,
        tols: Vec<f64>,
        kind: ConstraintKind,
        passes: f64,
    ) -> MeritTarget {
        let mut t = entry_phase(key_idx, wl, targets, tols, kind, SimTransform::Phase, 1.0);
        t.differential_passes = Some(passes);
        t
    }

    fn sim_pd() -> SimCurves {
        // One angle (0°), NW wavelengths from 400 nm; Ts complex row
        // 0.7·e^{i·0.3} at 400 nm; stack D = 100 nm of air (n = 1).
        let mut sim = SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: [None, None, None, None, None, None, None, None, None],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            total_d: 100.0,
            n_front_re: Arc::from([1.0_f64]),
            n_back_re: Arc::from([1.0_f64]),
        };
        sim.cplx[3] = Some(Arc::from(vec![
            Complex64::from_polar(0.7, 0.3),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]));
        sim
    }

    #[test]
    fn differential_phase_subtracts_reference() {
        // λ = 400, D = 100, n = 1, θ = 0: ref = 2π·100/400 = π/2 ≈ 1.5707963.
        // Δφ = 0.3 − π/2; demanding exactly that with tol 0.05 → zero.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        let delta = 0.3 - std::f64::consts::PI / 2.0;
        spec.add_target(entry_pd(
            k as u32,
            vec![400.0],
            vec![delta],
            vec![0.05],
            ConstraintKind::Exact,
            1.0,
        ))
        .unwrap();
        let sim = sim_pd();
        assert!(spec.merit(&sim, 1e6) < 1e-28);
        // Same demand as absolute (passes path off): residual is −π/2 →
        // (−π/2/0.05)² ≈ 986.96 — the reference is doing the work.
        let mut abs_spec = MeritSpec::new();
        let ka = abs_spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        abs_spec
            .add_target(entry_phase(
                ka as u32,
                vec![400.0],
                vec![delta],
                vec![0.05],
                ConstraintKind::Exact,
                SimTransform::Phase,
                1.0,
            ))
            .unwrap();
        let m_abs = abs_spec.merit(&sim, 1e6);
        let expect = ((0.3 - delta) / 0.05).powi(2);
        assert!(
            (m_abs - expect).abs() < 1e-9,
            "m_abs={m_abs} expect={expect}"
        );
    }

    #[test]
    fn differential_phase_zero_d_is_absolute() {
        // D = 0 kills the reference: differential ≡ absolute bit-for-bit.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec.add_target(entry_pd(
            k as u32,
            vec![400.0],
            vec![0.3],
            vec![0.05],
            ConstraintKind::Exact,
            1.0,
        ))
        .unwrap();
        let mut sim = sim_pd();
        sim.total_d = 0.0;
        assert!(spec.merit(&sim, 1e6) < 1e-28);
    }

    #[test]
    fn differential_phase_passes_scale() {
        // passes = 2 doubles the subtracted reference (round-trip).
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        let delta = 0.3 - std::f64::consts::PI; // 2 × π/2
        spec.add_target(entry_pd(
            k as u32,
            vec![400.0],
            vec![delta],
            vec![0.05],
            ConstraintKind::Exact,
            2.0,
        ))
        .unwrap();
        assert!(spec.merit(&sim_pd(), 1e6) < 1e-28);
    }

    #[test]
    fn differential_validation() {
        // Without phase → Err; negative/NaN passes → Err.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        let mut t = entry(
            k as u32,
            vec![400.0],
            vec![0.0],
            vec![0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        t.differential_passes = Some(1.0);
        assert!(spec.add_target(t).is_err());
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let mut t2 = entry_pd(
                k as u32,
                vec![400.0],
                vec![0.0],
                vec![0.1],
                ConstraintKind::Exact,
                bad,
            );
            t2.differential_passes = Some(bad);
            assert!(spec.add_target(t2).is_err(), "passes={bad}");
        }
    }

    #[test]
    fn weight_and_count_scale_merit() {
        // Base: two Exact points, nf = 1, tol 0.1, sim 0.1 off →
        // r = ±1/point → merit 2.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry(
            k as u32,
            vec![400.0, 500.0],
            vec![0.5, 0.5],
            vec![0.1, 0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let mut sim = sim_one_angle(&[0.6, 0.4, 0.0, 0.0, 0.0]);
        assert!((spec.merit(&sim, 1e6) - 2.0).abs() < 1e-12);
        // weight 2 → merit 4 (residuals scale by √2).
        let mut sw = MeritSpec::new();
        let kw = sw.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        let mut tw = entry(
            kw as u32,
            vec![400.0, 500.0],
            vec![0.5, 0.5],
            vec![0.1, 0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        tw.weight = 2.0;
        sw.add_target(tw).unwrap();
        assert!((sw.merit(&sim, 1e6) - 4.0).abs() < 1e-12);
        let mut out = Vec::new();
        sw.residuals(&sim, &mut out).unwrap();
        assert!((out[0].abs() - 2.0f64.sqrt()).abs() < 1e-12);
        // count 2 → merit 1 (mean, not sum).
        let mut sc = MeritSpec::new();
        let kc = sc.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        let mut tc = entry(
            kc as u32,
            vec![400.0, 500.0],
            vec![0.5, 0.5],
            vec![0.1, 0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        tc.count_norm = Some(2.0);
        sc.add_target(tc).unwrap();
        assert!((sc.merit(&sim, 1e6) - 1.0).abs() < 1e-12);
        // weight 3 + count 2 → 3.
        let mut sb = MeritSpec::new();
        let kb = sb.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        let mut tb = entry(
            kb as u32,
            vec![400.0, 500.0],
            vec![0.5, 0.5],
            vec![0.1, 0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        tb.weight = 3.0;
        tb.count_norm = Some(2.0);
        sb.add_target(tb).unwrap();
        assert!((sb.merit(&sim, 1e6) - 3.0).abs() < 1e-12);
        // Trust boundary: negative/NaN weight, non-positive count rejected.
        for (w, c) in [
            (-1.0, None),
            (f64::NAN, None),
            (1.0, Some(0.0)),
            (1.0, Some(-2.0)),
        ] {
            let mut sx = MeritSpec::new();
            let kx = sx.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            let mut tx = entry(
                kx as u32,
                vec![400.0],
                vec![0.5],
                vec![0.1],
                ConstraintKind::Exact,
                SimTransform::Linear,
                1.0,
            );
            tx.weight = w;
            tx.count_norm = c;
            assert!(sx.add_target(tx).is_err(), "w={w} c={c:?}");
        }
        let _ = &mut sim;
    }

    // --- F2.2: the multi-environment doors -----------------------------

    /// `Rs = 0.5 ± 0.1` on one point, at environment `env`.
    fn env_entry(key_idx: u32, env: u32) -> MeritTarget {
        let mut t = entry(
            key_idx,
            vec![400.0],
            vec![0.5],
            vec![0.1],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        );
        t.env_idx = env;
        t
    }

    /// The residual vector concatenates ENVIRONMENT-MAJOR: every residual
    /// of environment 0 in `residuals()` order, then environment 1's.
    /// Insertion order stays minor inside each environment.
    #[test]
    fn residuals_concatenate_environment_major() {
        let mut spec = MeritSpec::new();
        spec.set_n_envs(2).unwrap();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }) as u32;
        // Two demands in env 0 (so insertion order is observable), one in
        // env 1 - a deliberately ragged split: the concatenation must not
        // assume the environments contribute equally.
        spec.add_target(env_entry(k, 0)).unwrap();
        let mut second = env_entry(k, 0);
        second.tolerances = vec![0.2].into();
        spec.add_target(second).unwrap();
        spec.add_target(env_entry(k, 1)).unwrap();

        let a = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]); // diff +0.1
        let b = sim_one_angle(&[0.9, 0.0, 0.0, 0.0, 0.0]); // diff +0.4
        let mut out = Vec::new();
        spec.residuals_multi(&[a, b], &mut out).unwrap();
        assert_eq!(out.len(), 3, "{out:?}");
        assert!((out[0] - 1.0).abs() < 1e-12, "{out:?}"); // env0, tol 0.1
        assert!((out[1] - 0.5).abs() < 1e-12, "{out:?}"); // env0, tol 0.2
        assert!((out[2] - 4.0).abs() < 1e-12, "{out:?}"); // env1, tol 0.1
        assert_eq!(out.len(), spec.n_residuals());
    }

    /// The penalty group is the `(environment, key)` PAIR. A curve missing
    /// in one environment costs one penalty, not one per environment and
    /// not none.
    #[test]
    fn a_missing_curve_costs_one_penalty_per_environment_key_pair() {
        let mut spec = MeritSpec::new();
        spec.set_n_envs(2).unwrap();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rp, // sim_one_angle fills Rs only
        }) as u32;
        spec.add_target(env_entry(k, 0)).unwrap();
        spec.add_target(env_entry(k, 1)).unwrap();
        let a = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        let b = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(spec.merit_multi(&[a, b], 7.0), 14.0);

        // Only environment 1 demands the missing curve: one penalty.
        let mut half = MeritSpec::new();
        half.set_n_envs(2).unwrap();
        let ok = half.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }) as u32;
        let bad = half.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rp,
        }) as u32;
        half.add_target(env_entry(ok, 0)).unwrap();
        half.add_target(env_entry(bad, 1)).unwrap();
        let a = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        let b = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        // env0: Rs present, residual 1 → 1.0. env0/Rp: no demand → free.
        // env1/Rs: no demand → free. env1/Rp: demanded, missing → 7.0.
        // The residual is 1.0 only to rounding (0.6 − 0.5 is not exact in
        // binary); the PENALTY half of the sum is what this pins.
        assert!((half.merit_multi(&[a, b], 7.0) - 8.0).abs() < 1e-12);
    }

    /// An `(environment, key)` pair no demand belongs to is skipped
    /// outright. It has no curve to miss, so charging it a penalty would
    /// make the merit depend on which keys some OTHER environment happens
    /// to have registered.
    #[test]
    fn an_empty_environment_key_pair_is_free() {
        let mut spec = MeritSpec::new();
        spec.set_n_envs(3).unwrap();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }) as u32;
        spec.add_target(env_entry(k, 1)).unwrap();
        let s = || sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        // 1.0 to rounding: one residual, no penalty from the two empty
        // pairs. A charged pair would land at 8.0 or 15.0, not near 1.
        assert!((spec.merit_multi(&[s(), s(), s()], 7.0) - 1.0).abs() < 1e-12);
        let mut out = Vec::new();
        spec.residuals_multi(&[s(), s(), s()], &mut out).unwrap();
        assert_eq!(out.len(), 1);
    }

    /// `merit`/`residuals` are the one-element case, not a separate path.
    #[test]
    fn the_single_environment_doors_delegate() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }) as u32;
        spec.add_target(env_entry(k, 0)).unwrap();
        let sim = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(
            spec.merit(&sim, 1e6).to_bits(),
            spec.merit_multi(std::slice::from_ref(&sim), 1e6).to_bits()
        );
        let (mut a, mut b) = (Vec::new(), Vec::new());
        spec.residuals(&sim, &mut a).unwrap();
        spec.residuals_multi(std::slice::from_ref(&sim), &mut b)
            .unwrap();
        assert_eq!(a, b);
    }

    /// Scoring a two-environment spec against one solve would silently
    /// drop every demand tagged with the other environment. Fail closed.
    #[test]
    #[should_panic(expected = "1 simulation(s) for 2 environment(s)")]
    fn one_solve_for_two_environments_panics() {
        let mut spec = MeritSpec::new();
        spec.set_n_envs(2).unwrap();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }) as u32;
        spec.add_target(env_entry(k, 1)).unwrap();
        let sim = sim_one_angle(&[0.6, 0.0, 0.0, 0.0, 0.0]);
        spec.merit(&sim, 1e6);
    }

    fn entry_integral(key_idx: u32, kind: ConstraintKind) -> MeritTarget {
        let mut t = entry(
            key_idx,
            vec![400.0, 500.0, 600.0],
            vec![0.5, 0.5, 0.5],
            vec![0.1, 0.1, 0.1],
            kind,
            SimTransform::Linear,
            1.0,
        );
        t.integral = true;
        t
    }

    #[test]
    fn integral_mean_single_residual() {
        // Sim [0.6, 0.5, 0.4] vs 0.5: mean diff 0 → merit 0, ONE residual.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry_integral(k as u32, ConstraintKind::Exact))
            .unwrap();
        let sim = sim_one_angle(&[0.6, 0.5, 0.4, 0.0, 0.0]);
        assert!(spec.merit(&sim, 1e6) < 1e-28);
        let mut out = Vec::new();
        spec.residuals(&sim, &mut out).unwrap();
        assert_eq!(out.len(), 1);
        // Sim [0.7, 0.6, 0.5]: mean diff 0.1, tol 0.1 → R = 1 → merit 1.
        let sim2 = sim_one_angle(&[0.7, 0.6, 0.5, 0.0, 0.0]);
        assert!((spec.merit(&sim2, 1e6) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn integral_kinds_mask_the_mean() {
        // Above: mean above target → silent; mean below → active.
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_target(entry_integral(k as u32, ConstraintKind::Above))
            .unwrap();
        // Mean 0.6 ≥ 0.5 → 0 (even though point 600 dips to 0.4!).
        let sim_hi = sim_one_angle(&[0.7, 0.7, 0.4, 0.0, 0.0]);
        assert!(spec.merit(&sim_hi, 1e6) < 1e-28);
        // Mean 0.4 < 0.5 → ((0.4−0.5)/0.1)² = 1.
        let sim_lo = sim_one_angle(&[0.4, 0.4, 0.4, 0.0, 0.0]);
        assert!((spec.merit(&sim_lo, 1e6) - 1.0).abs() < 1e-12);
        // Range with band 0.05 (raw, nf = 1): mean inside → 0.
        let mut spec_r = MeritSpec::new();
        let kr = spec_r.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        let mut tr = entry_integral(kr as u32, ConstraintKind::Range);
        tr.band = vec![0.05, 0.05, 0.05].into();
        spec_r.add_target(tr).unwrap();
        let sim_in = sim_one_angle(&[0.54, 0.54, 0.54, 0.0, 0.0]);
        assert!(spec_r.merit(&sim_in, 1e6) < 1e-28);
        // Mean 0.6: exceedance (0.1−0.05)/0.1 = 0.5 → 0.25.
        let sim_out = sim_one_angle(&[0.6, 0.6, 0.6, 0.0, 0.0]);
        assert!((spec_r.merit(&sim_out, 1e6) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn integral_rejects_count_norm() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        let mut t = entry_integral(k as u32, ConstraintKind::Exact);
        t.count_norm = Some(3.0);
        assert!(spec.add_target(t).is_err());
    }

    #[test]
    fn back_intensity_and_absorption() {
        // RBs row 0.4, TBs row 0.5: R demand 0.4 → 0; ABs demand 0.1 → 0.
        let mut spec = MeritSpec::new();
        let kr = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::RBs,
        });
        spec.add_target(entry(
            kr as u32,
            vec![400.0],
            vec![0.4],
            vec![0.05],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let ka = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::ABs,
        });
        spec.add_target(entry(
            ka as u32,
            vec![400.0],
            vec![0.1],
            vec![0.05],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        let mut sim = SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            curves: [None, None, None, None, None, None, None, None, None],
            back: [None, None, None, None, None, None],
            cplx: [None, None, None, None, None, None],
            cplx_back: [None, None, None, None],
            ..Default::default()
        };
        sim.back[0] = Some(Arc::from(vec![0.4, 0., 0., 0., 0.]));
        sim.back[3] = Some(Arc::from(vec![0.5, 0., 0., 0., 0.]));
        assert!(spec.merit(&sim, 1e6) < 1e-28); // 0 + (1−0.4−0.5−0.1)²
    }

    #[test]
    fn band_length_mismatch_rejected() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        assert!(
            spec.add_target(entry_banded(
                k as u32,
                vec![400.0],
                vec![0.0],
                vec![1.0],
                vec![0.1, 0.2],
                ConstraintKind::Range,
                SimTransform::Linear,
                1.0
            ))
            .is_err()
        );
    }

    #[test]
    fn rotation_kernel_values() {
        // λ=1000, θ=0, n=1, d=250, passes=1: ref = 2π·250/1000 = π/2.
        let rot = super::reference_rotation(&[1000.0], 0.0, &[1.0], 250.0, 1.0).unwrap();
        assert!((rot[0].re - 0.0).abs() < 1e-15);
        assert!((rot[0].im + 1.0).abs() < 1e-15);
        let mut rows = vec![num_complex::Complex64::new(1.0, 0.0); 3];
        super::rotate_rows(&mut rows, &rot).unwrap();
        assert!((rows[0].im + 1.0).abs() < 1e-15);
        assert!(super::rotate_rows(&mut rows[..2], &rot).is_ok());
        let mut bad = vec![num_complex::Complex64::new(1.0, 0.0); 2];
        let rot2 =
            super::reference_rotation(&[1000.0, 500.0, 250.0], 0.0, &[1.0], 0.0, 1.0).unwrap();
        assert!(super::rotate_rows(&mut bad, &rot2).is_err());
    }

    /// PD2: the kernel broadcasts a length-1 index bitwise and refuses a
    /// length that is neither 1 nor the grid length, naming both numbers.
    #[test]
    fn reference_rotation_broadcasts_and_refuses_by_length() {
        use std::f64::consts::TAU;
        let wls = [400.0_f64, 500.0, 600.0];
        let n_col = [1.5_f64, 1.7, 1.9];
        // Broadcast (length 1) == the pre-PD2 scalar kernel, bitwise:
        // identical f64 operand, identical order (exp(-i*r) itself).
        let rot_b = super::reference_rotation(&wls, 0.0, &[1.7], 100.0, 1.0).unwrap();
        let rot_x = super::reference_rotation(&wls, 0.0, &[1.7, 1.7, 1.7], 100.0, 1.0).unwrap();
        assert_eq!(rot_b.to_vec(), rot_x);
        // Hand: arg(exp(-i*r)) = -r (wrapped to (-pi, pi]).
        let wrapped = |x: f64| x - TAU * (x / TAU + 0.5).floor();
        let expect = wrapped(-TAU * 1.7 * 100.0 / 400.0);
        assert!(
            (rot_b[0].arg() - expect).abs() < 1e-12,
            "arg={:?}",
            rot_b[0]
        );
        // Dispersive row: per-lambda values.
        let rot_d = super::reference_rotation(&wls, 0.0, &n_col, 200.0, 1.0).unwrap();
        let expect_args: Vec<f64> = n_col
            .iter()
            .zip(wls)
            .map(|(n, w)| wrapped(-TAU * n * 200.0 / w))
            .collect();
        for (r, e) in rot_d.iter().zip(expect_args.iter()) {
            assert!((r.arg() - e).abs() < 1e-12, "arg={r:?} expect={e}");
        }
        // Refusals: length 2 against a 3-wavelength grid (both named), empty.
        let err = super::reference_rotation(&wls, 0.0, &[1.5, 1.7], 200.0, 1.0).unwrap_err();
        assert!(
            err.contains("length 2") && err.contains("3 wavelengths"),
            "err={err}"
        );
        let err = super::reference_rotation(&wls, 0.0, &[], 200.0, 1.0).unwrap_err();
        assert!(err.contains("empty"), "err={err}");
    }

    /// PD2 / review PD G1: the scalar guard reports only the sides a
    /// demand can actually read. `differential()` maps only front curves
    /// (PDts/PDtp), so a front-only spec with a default scalar `n_back`
    /// is correct and must be silent; a scalar `n_front` on the same spec
    /// is still reported, alone.
    #[test]
    fn demanded_reference_sides_tracks_the_labels() {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec.add_target(entry_pd(
            k as u32,
            vec![400.0],
            vec![0.0],
            vec![0.1],
            ConstraintKind::Exact,
            1.0,
        ))
        .unwrap();
        assert_eq!(spec.demanded_reference_sides(), (true, false));
        // A back curve (hypothetical label, decision 6) flips the side.
        let mut spec_b = MeritSpec::new();
        let kb = spec_b.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::TBs,
        });
        spec_b
            .add_target(entry_pd(
                kb as u32,
                vec![400.0],
                vec![0.0],
                vec![0.1],
                ConstraintKind::Exact,
                1.0,
            ))
            .unwrap();
        assert_eq!(spec_b.demanded_reference_sides(), (false, true));
        // A non-differential spec demands nothing.
        let mut spec_p = MeritSpec::new();
        let kp = spec_p.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec_p
            .add_target(entry_phase(
                kp as u32,
                vec![400.0],
                vec![0.0],
                vec![0.1],
                ConstraintKind::Exact,
                SimTransform::Phase,
                1.0,
            ))
            .unwrap();
        assert_eq!(spec_p.demanded_reference_sides(), (false, false));
    }

    #[test]
    fn curve_setters_validate() {
        use super::{CurveId, SimCurves};
        use std::sync::Arc;
        let mut sim = SimCurves {
            angles: Arc::from([0.0]),
            wavelengths: Arc::from([500.0, 600.0]),
            total_d: 0.0,
            n_front_re: Arc::from([1.0_f64]),
            n_back_re: Arc::from([1.0_f64]),
            curves: Default::default(),
            back: Default::default(),
            cplx: Default::default(),
            cplx_back: Default::default(),
        };
        assert!(sim.set_curve(CurveId::Rs, Arc::from([0.1, 0.2])).is_ok());
        assert!(sim.set_curve(CurveId::Rs, Arc::from([0.1])).is_err());
        assert!(sim.set_curve(CurveId::As, Arc::from([0.1, 0.2])).is_err());
        assert!(
            sim.set_complex(
                CurveId::Rs,
                Arc::from([num_complex::Complex64::new(1.0, 0.0); 2])
            )
            .is_ok()
        );
        assert!(
            sim.set_complex(
                CurveId::Ru,
                Arc::from([num_complex::Complex64::new(1.0, 0.0); 2])
            )
            .is_err()
        );
    }

    #[test]
    fn length_mismatch_and_bad_key_rejected() {
        let mut spec = MeritSpec::new();
        let _k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        assert!(
            spec.add_target(entry(
                0,
                vec![400.0],
                vec![0.0, 0.0],
                vec![1.0],
                ConstraintKind::Exact,
                SimTransform::Linear,
                1.0
            ))
            .is_err()
        );
        assert!(
            spec.add_target(entry(
                7,
                vec![400.0],
                vec![0.0],
                vec![1.0],
                ConstraintKind::Exact,
                SimTransform::Linear,
                1.0
            ))
            .is_err()
        );
    }

    /// Toy native tables on the sim grid (node-exact resample).
    fn color_toy() -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>, Vec<f64>) {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + 100.0 * i as f64).collect();
        let cmf: Vec<[f64; 3]> = (0..NW)
            .map(|i| {
                [
                    0.10 + 0.01 * i as f64,
                    0.20 + 0.01 * i as f64,
                    0.05 + 0.005 * i as f64,
                ]
            })
            .collect();
        (wl.clone(), cmf, wl, vec![1.0; NW])
    }

    fn color_spec_xyy(target_y: f64) -> MeritSpec {
        use crate::smatrix::synthesis::color_merit::{
            ColorDemand, ColorDistance, ColorQuantity, ColorReference,
        };
        let (cmf_wl, cmf, illum_wl, illuminant) = color_toy();
        // Achromatic reference: white chromaticity + target luminance.
        let ones = vec![1.0; NW];
        let white = crate::smatrix::synthesis::color_merit::xyz_of_spectrum(
            &ones,
            &illum_wl,
            &cmf,
            &cmf_wl,
            &illuminant,
            &illum_wl,
            None,
        )
        .unwrap();
        let s = white[0] + white[1] + white[2];
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        });
        spec.add_color_demand(
            ColorDemand::new(
                k as u32,
                cmf,
                cmf_wl,
                illuminant,
                illum_wl,
                ColorQuantity::XyY,
                ColorReference::Triple([white[0] / s, white[1] / s, target_y]),
                ColorDistance::Channels,
                1.0,
                crate::smatrix::synthesis::color_merit::E313_CX_D65_10,
                crate::smatrix::synthesis::color_merit::E313_CZ_D65_10,
                None,
            )
            .unwrap(),
        )
        .unwrap();
        spec
    }

    #[test]
    fn color_missing_curve_takes_one_penalty() {
        // Shared key: color + spectral demands fail as ONE group (no double
        // penalty) — parity with pointwise missing-curve handling.
        let mut spec = color_spec_xyy(0.5);
        spec.add_target(entry(
            0,
            vec![400.0],
            vec![0.5],
            vec![0.01],
            ConstraintKind::Exact,
            SimTransform::Linear,
            1.0,
        ))
        .unwrap();
        assert_eq!(spec.n_residuals(), 1 + 1);
        let full = sim_one_angle(&[0.5; NW]);
        assert!(spec.merit(&full, 1e6).is_finite());
        let mut out = Vec::new();
        assert!(spec.residuals(&full, &mut out).is_ok());
        assert_eq!(out.len(), 2);
        let empty = SimCurves {
            angles: vec![0.0].into(),
            wavelengths: (0..NW)
                .map(|i| 400.0 + 100.0 * i as f64)
                .collect::<Vec<_>>()
                .into(),
            ..Default::default()
        };
        assert_eq!(spec.merit(&empty, 1e6), 1e6);
        assert_eq!(spec.residuals(&empty, &mut out).unwrap_err(), CurveId::Rs);
    }

    #[test]
    fn lm_drives_color_residual_to_zero() {
        // Thick-opt consumes color residuals with zero further changes: a
        // 1-param uniform-reflector problem (Y(R) = R by k-normalization)
        // converges from 0.2 to the 0.5 target.
        use crate::smatrix::synthesis::thick_opt::{LmConfig, levenberg_marquardt};
        let spec = color_spec_xyy(0.5);
        let r0 = {
            let mut out = Vec::new();
            spec.residuals(&sim_one_angle(&[0.2; NW]), &mut out)
                .unwrap();
            out.iter().map(|r| r * r).sum::<f64>()
        };
        let res = levenberg_marquardt(
            &|x: &[f64], out: &mut Vec<f64>| {
                let row = [x[0]; NW];
                spec.residuals(&sim_one_angle(&row), out)
                    .map_err(|id| format!("{id:?}"))
            },
            &[0.2],
            &[0.01],
            &[1.0],
            &LmConfig::default(),
        )
        .unwrap();
        assert!(res.cost < r0);
        assert!((res.x[0] - 0.5).abs() < 1e-6, "x={}", res.x[0]);
        // Floor is 1-ulp chromaticity noise (x/y of a uniform spectrum
        // reproduce white to rounding), not optimizer failure.
        assert!(res.cost < 1e-14, "cost={}", res.cost);
        assert!(res.x[0].is_finite());
    }

    // -------------------------------------------------------------------------
    // R4.5: row-wise curve sensitivity, against a finite difference of
    // `residuals()` itself. This is also the anti-drift guard between
    // `curve_sensitivity` and `residuals_into`, which walk the target grids
    // separately on purpose (see the doc comment on `curve_sensitivity`).
    // -------------------------------------------------------------------------

    mod sensitivity {
        use super::*;
        use std::sync::Arc;

        const NW: usize = 9;
        const NA: usize = 2;

        fn wl() -> Vec<f64> {
            (0..NW).map(|i| 500.0 + 25.0 * i as f64).collect()
        }

        /// Two angle rows of smooth, distinguishable curves.
        fn sim_of(rs: &[f64], ts: &[f64]) -> SimCurves {
            let mut s = SimCurves {
                angles: vec![0.0, 30.0].into(),
                wavelengths: wl().into(),
                total_d: 400.0,
                n_front_re: Arc::from([1.0_f64]),
                n_back_re: Arc::from([1.52_f64]),
                ..Default::default()
            };
            s.set_curve(CurveId::Rs, Arc::from(rs.to_vec())).unwrap();
            s.set_curve(CurveId::Ts, Arc::from(ts.to_vec())).unwrap();
            s
        }

        fn base_sim() -> SimCurves {
            let rs: Vec<f64> = (0..NA * NW)
                .map(|k| 0.08 + 0.05 * ((k as f64) * 0.7).sin())
                .collect();
            let ts: Vec<f64> = (0..NA * NW)
                .map(|k| 0.80 + 0.04 * ((k as f64) * 0.4).cos())
                .collect();
            sim_of(&rs, &ts)
        }

        fn target(
            key_idx: u32,
            grid: Vec<f64>,
            kind: ConstraintKind,
            transform: SimTransform,
            integral: bool,
            weight: f64,
        ) -> MeritTarget {
            let n = grid.len();
            MeritTarget {
                key_idx,
                env_idx: 0,
                wavelengths: grid.into(),
                kind,
                transform,
                norm_factor: 1.0,
                normalized_targets: vec![0.05; n].into(),
                tolerances: vec![0.02; n].into(),
                band: vec![0.01; n].into(),
                phase: false,
                differential_passes: None,
                integral,
                weight,
                count_norm: None,
            }
        }

        /// Rebuild the sim with one (curve, flat index) entry shifted by `h`.
        fn bumped(base: &SimCurves, curve: CurveId, flat: usize, h: f64) -> SimCurves {
            let mut rs = base.curve(CurveId::Rs).unwrap().to_vec();
            let mut ts = base.curve(CurveId::Ts).unwrap().to_vec();
            match curve {
                CurveId::Rs => rs[flat] += h,
                CurveId::Ts => ts[flat] += h,
                other => panic!("unexpected curve {other:?}"),
            }
            sim_of(&rs, &ts)
        }

        /// Dense ∂r/∂value by central difference, for one (curve, flat) entry.
        fn fd_column(spec: &MeritSpec, sim: &SimCurves, curve: CurveId, flat: usize) -> Vec<f64> {
            let h = 1e-7;
            let mut plus = Vec::new();
            let mut minus = Vec::new();
            spec.residuals(&bumped(sim, curve, flat, h), &mut plus)
                .unwrap();
            spec.residuals(&bumped(sim, curve, flat, -h), &mut minus)
                .unwrap();
            plus.iter()
                .zip(&minus)
                .map(|(p, m)| (p - m) / (2.0 * h))
                .collect()
        }

        /// The analytic sensitivity as a dense column for one (curve, flat).
        fn analytic_column(
            s: &MeritSensitivity,
            n_wav: usize,
            curve: CurveId,
            flat: usize,
        ) -> Vec<f64> {
            s.rows
                .iter()
                .map(|terms| {
                    terms
                        .iter()
                        .filter(|t| t.curve == curve && t.angle_row * n_wav + t.wavelength == flat)
                        .map(|t| t.d_residual)
                        .sum()
                })
                .collect()
        }

        fn compare(spec: &MeritSpec, sim: &SimCurves, tol: f64, label: &str) {
            let s = spec.curve_sensitivity(sim).unwrap();
            // Against the residual vector itself, not `n_residuals()`: the two
            // differ by design for a frame that misses the simulated grid, and
            // the contract this pass owes its caller is index-for-index with
            // `residuals()`.
            let mut r0 = Vec::new();
            spec.residuals(sim, &mut r0).unwrap();
            assert_eq!(s.rows.len(), r0.len(), "{label}: row count");
            assert!(s.is_complete(), "{label}: uncovered {:?}", s.uncovered);

            let mut checked = 0usize;
            let mut worst = 0.0f64;
            for curve in [CurveId::Rs, CurveId::Ts] {
                for flat in 0..NA * NW {
                    let fd = fd_column(spec, sim, curve, flat);
                    let an = analytic_column(&s, NW, curve, flat);
                    let scale = fd.iter().fold(1.0f64, |a, v| a.max(v.abs()));
                    for (i, (a, f)) in an.iter().zip(&fd).enumerate() {
                        let dev = (a - f).abs() / scale;
                        worst = worst.max(dev);
                        assert!(dev < tol, "{label}: row {i}, {curve:?}[{flat}]: {a} vs {f}");
                        if f.abs() > 1e-9 {
                            checked += 1;
                        }
                    }
                }
            }
            assert!(
                checked > 0,
                "{label}: every finite difference was zero — vacuous"
            );
            println!("  {label}: {checked} live entries, worst rel dev {worst:.3e}");
        }

        #[test]
        fn sensitivity_is_a_finite_difference_of_residuals() {
            // Every kind x every real transform, pointwise, on the aligned grid.
            for kind in [
                ConstraintKind::Exact,
                ConstraintKind::Above,
                ConstraintKind::Below,
                ConstraintKind::Range,
                ConstraintKind::CenterBand,
            ] {
                for transform in [SimTransform::Linear, SimTransform::Log] {
                    let mut spec = MeritSpec::new();
                    let k = spec.add_key(MeritKey {
                        angle: 30.0,
                        curve: CurveId::Rs,
                    });
                    let mut t = target(k as u32, wl(), kind, transform, false, 1.7);
                    // Aim at the middle of the R sweep, in whatever units the
                    // transform works in. A level the curve never reaches would leave
                    // `a`/`b` inactive at every point and the comparison vacuous —
                    // which is exactly what a flat 0.05 does to `b` under `Log`.
                    let level = if transform == SimTransform::Log {
                        0.08f64.log10()
                    } else {
                        0.08
                    };
                    t.normalized_targets = vec![level; NW].into();
                    spec.add_target(t).unwrap();
                    compare(&spec, &base_sim(), 1e-6, &format!("{kind:?}/{transform:?}"));
                }
            }
        }

        #[test]
        fn sensitivity_covers_the_integral_mean_row() {
            // One residual over the mean of nine points: every point contributes
            // 1/n, and the kind applies once, to the mean.
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            spec.add_target(target(
                k as u32,
                wl(),
                ConstraintKind::Exact,
                SimTransform::Linear,
                true,
                3.0,
            ))
            .unwrap();
            assert_eq!(spec.n_residuals(), 1);
            compare(&spec, &base_sim(), 1e-6, "integral/Exact");
        }

        #[test]
        fn sensitivity_covers_absorption_through_both_companions() {
            // A = 1 − R − T: the row must depend on BOTH curves, each with the
            // opposite sign of an R row. A one-channel bug would still pass a
            // single-curve check.
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 30.0,
                curve: CurveId::As,
            });
            spec.add_target(target(
                k as u32,
                wl(),
                ConstraintKind::Exact,
                SimTransform::Linear,
                false,
                1.0,
            ))
            .unwrap();
            let sim = base_sim();
            compare(&spec, &sim, 1e-6, "absorption");

            let s = spec.curve_sensitivity(&sim).unwrap();
            let curves: Vec<CurveId> = s.rows[0].iter().map(|t| t.curve).collect();
            assert!(
                curves.contains(&CurveId::Rs) && curves.contains(&CurveId::Ts),
                "absorption row reads {curves:?}"
            );
            assert!(
                s.rows[0].iter().all(|t| t.d_residual < 0.0),
                "A = 1 − R − T: raising either companion must lower the residual"
            );
        }

        #[test]
        fn sensitivity_covers_an_interpolated_target_grid() {
            // Off-grid target points read two simulated points with weights
            // (1−f, f). A pass that only handled the aligned fast path would look
            // perfect on every other test in this module.
            let grid: Vec<f64> = (0..7).map(|i| 512.0 + 31.0 * i as f64).collect();
            let sim = base_sim();
            // Under `Linear`/`Exact` the derivative is the same number wherever
            // the row is evaluated, so that pair pins the interpolation WEIGHTS
            // and nothing else. `Log` makes ∂r/∂value depend on the interpolated
            // value itself, which pins the point the derivative is taken AT —
            // reading the left neighbour instead of interpolating passes the
            // first and fails the second.
            for transform in [SimTransform::Linear, SimTransform::Log] {
                let mut spec = MeritSpec::new();
                let k = spec.add_key(MeritKey {
                    angle: 0.0,
                    curve: CurveId::Rs,
                });
                let mut t = target(
                    k as u32,
                    grid.clone(),
                    ConstraintKind::Exact,
                    transform,
                    false,
                    1.0,
                );
                if transform == SimTransform::Log {
                    t.normalized_targets = vec![0.08f64.log10(); grid.len()].into();
                }
                spec.add_target(t).unwrap();
                compare(&spec, &sim, 1e-6, &format!("interpolated/{transform:?}"));

                let s = spec.curve_sensitivity(&sim).unwrap();
                assert!(
                    s.rows.iter().any(|r| r.len() == 2),
                    "no row read two simulated points — the grid was not off-grid"
                );
            }
        }

        #[test]
        fn sensitivity_covers_several_keys_targets_and_angles_at_once() {
            // Row ORDER is the contract: keys in registration order, targets per
            // key in insertion order, points along the grid. A mis-ordered pass
            // would still be elementwise correct on any single-target spec.
            let mut spec = MeritSpec::new();
            let k0 = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            let k1 = spec.add_key(MeritKey {
                angle: 30.0,
                curve: CurveId::Ts,
            });
            spec.add_target(target(
                k0 as u32,
                wl(),
                ConstraintKind::Exact,
                SimTransform::Linear,
                false,
                1.0,
            ))
            .unwrap();
            spec.add_target(target(
                k0 as u32,
                wl()[2..6].to_vec(),
                ConstraintKind::Below,
                SimTransform::Linear,
                true,
                2.5,
            ))
            .unwrap();
            spec.add_target(target(
                k1 as u32,
                wl(),
                ConstraintKind::CenterBand,
                SimTransform::Linear,
                false,
                0.4,
            ))
            .unwrap();
            compare(&spec, &base_sim(), 1e-6, "multi-key");
        }

        #[test]
        fn a_target_grid_that_misses_the_simulation_contributes_no_rows() {
            // `residuals_into` pushes nothing at all for a non-overlapping frame.
            // If this pass pushed empty rows instead, every row after it would be
            // misaligned — silently, since the values would all still be finite.
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            spec.add_target(target(
                k as u32,
                vec![1200.0, 1300.0],
                ConstraintKind::Exact,
                SimTransform::Linear,
                false,
                1.0,
            ))
            .unwrap();
            spec.add_target(target(
                k as u32,
                wl(),
                ConstraintKind::Exact,
                SimTransform::Linear,
                false,
                1.0,
            ))
            .unwrap();
            let sim = base_sim();
            let mut r = Vec::new();
            spec.residuals(&sim, &mut r).unwrap();
            let s = spec.curve_sensitivity(&sim).unwrap();
            assert_eq!(s.rows.len(), r.len(), "row counts diverged on a grid miss");
            compare(&spec, &sim, 1e-6, "grid-miss");
        }

        #[test]
        fn a_phase_target_is_reported_uncovered_not_zero() {
            // The whole point of `uncovered`: a caller must fall back to a finite
            // difference, not conclude that the phase rows are constant.
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            let mut t = target(
                k as u32,
                wl(),
                ConstraintKind::Exact,
                SimTransform::Phase,
                false,
                1.0,
            );
            t.phase = true;
            spec.add_target(t).unwrap();

            let mut sim = base_sim();
            let cplx: Vec<Complex64> = (0..NA * NW)
                .map(|k| {
                    Complex64::new(
                        0.3 * ((k as f64) * 0.3).cos(),
                        0.2 * ((k as f64) * 0.5).sin(),
                    )
                })
                .collect();
            sim.set_complex(CurveId::Rs, Arc::from(cplx)).unwrap();
            let s = spec.curve_sensitivity(&sim).unwrap();
            assert!(!s.is_complete());
            assert_eq!(s.uncovered.len(), NW);
            assert_eq!(s.rows.len(), spec.n_residuals());
            assert!(s.rows.iter().all(|r| r.is_empty()));
        }

        #[test]
        fn a_missing_curve_errors_where_the_residual_pass_errors() {
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rp,
            });
            spec.add_target(target(
                k as u32,
                wl(),
                ConstraintKind::Exact,
                SimTransform::Linear,
                false,
                1.0,
            ))
            .unwrap();
            let sim = base_sim();
            let mut r = Vec::new();
            assert_eq!(spec.residuals(&sim, &mut r).unwrap_err(), CurveId::Rp);
            assert_eq!(spec.curve_sensitivity(&sim).unwrap_err(), CurveId::Rp);
        }

        #[test]
        fn the_weight_and_count_normalization_ride_through() {
            // rscale = sqrt(weight / count_norm) multiplies the residual, so it
            // multiplies its derivative too. Two specs differing only in weight
            // must differ in sensitivity by exactly that ratio.
            let build = |weight: f64, count: Option<f64>| {
                let mut spec = MeritSpec::new();
                let k = spec.add_key(MeritKey {
                    angle: 0.0,
                    curve: CurveId::Rs,
                });
                let mut t = target(
                    k as u32,
                    wl(),
                    ConstraintKind::Exact,
                    SimTransform::Linear,
                    false,
                    weight,
                );
                t.count_norm = count;
                spec.add_target(t).unwrap();
                spec
            };
            let sim = base_sim();
            let a = build(1.0, None).curve_sensitivity(&sim).unwrap();
            let b = build(9.0, Some(4.0)).curve_sensitivity(&sim).unwrap();
            let ratio = (9.0f64 / 4.0).sqrt();
            for (ra, rb) in a.rows.iter().zip(&b.rows) {
                for (ta, tb) in ra.iter().zip(rb) {
                    assert!(
                        (tb.d_residual - ta.d_residual * ratio).abs()
                            <= 1e-14 * tb.d_residual.abs().max(1.0),
                        "{} vs {}",
                        tb.d_residual,
                        ta.d_residual * ratio
                    );
                }
            }
            compare(&build(9.0, Some(4.0)), &sim, 1e-6, "weighted");
        }

        #[test]
        fn an_inactive_constraint_has_no_terms_at_all() {
            // A `b` (below) demand already satisfied everywhere is flat: the rows
            // must carry no terms, not terms that happen to be zero, so a caller
            // assembling J does no work for them.
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            let mut t = target(
                k as u32,
                wl(),
                ConstraintKind::Below,
                SimTransform::Linear,
                false,
                1.0,
            );
            t.normalized_targets = vec![10.0; NW].into(); // R is nowhere near 10
            spec.add_target(t).unwrap();

            let sim = base_sim();
            let s = spec.curve_sensitivity(&sim).unwrap();
            assert!(s.is_complete());
            assert!(s.rows.iter().all(|r| r.is_empty()), "{:?}", s.rows);

            let mut r = Vec::new();
            spec.residuals(&sim, &mut r).unwrap();
            assert!(
                r.iter().all(|v| *v == 0.0),
                "the constraint should be inactive"
            );
        }
    }
}
