//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::optimizer — which least-squares solver runs, and on what
//! problem.
//!
//! The residual system is already solver-agnostic: `MeritSpec::residuals` is
//! a closure over `x`, and `JacobianSource` is a closure over `x` that
//! returns `∂r/∂x`. This module is the one place that decides *who* is handed
//! that pair, so an alternative solver is a new arm here rather than a second
//! copy of the call site (R4.4c, review §3.6 / §18.2).
//!
//! # Bounds are not a shared contract
//!
//! Navette's built-in LM keeps `lb ≤ x ≤ ub` by vetoing and clamping the step
//! it just solved, so an optimum **may sit exactly on a bound** — a film
//! driven to zero thickness is a real answer the synthesis loop then removes.
//! Every reference implementation in the ecosystem is *unbounded*
//! ([`OptimizerBackend::MinpackLm`] included), so those backends run on an
//! interior reparametrization ([`IntervalMap`]) and their optima are
//! **strictly inside** the box. That is a contract difference, not a bug, and
//! it is why the built-in stays the default for bounded problems.
//!
//! # Feature gating
//!
//! Backends that need a dependency are declared unconditionally and
//! *implemented* behind a cargo feature, so a build without the feature can
//! still name the backend and be told how to get it. Selecting one that is
//! not compiled in is an error with a rebuild hint, never a silent fallback
//! to a different solver.

use super::thick_opt::{
    levenberg_marquardt_with, JacobianSource, LmConfig, LmResult, LmTermination,
};

// ---------------------------------------------------------------------------
// Backend selection
// ---------------------------------------------------------------------------

/// Which least-squares solver runs.
///
/// Variants exist whether or not their implementation was compiled in — a
/// name that cannot be spelled cannot be diagnosed, and "rebuild with
/// `--features opt-minpack-lm`" is the useful answer to selecting one that is
/// missing.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum OptimizerBackend {
    /// Navette's own bounded Levenberg-Marquardt (`thick_opt`): QR step
    /// solve, gain-ratio damping, veto+clamp bounds, analytic or differenced
    /// Jacobian. The default, and the only backend with bounds semantics of
    /// its own.
    #[default]
    BuiltinLm,
    /// The `levenberg-marquardt` crate (rust-cv, MINPACK `lmdif`-derived) —
    /// the ecosystem's reference LM, behind `opt-minpack-lm`. Unbounded: runs
    /// on the interior reparametrization.
    MinpackLm,
    /// Trust-region reflective (Branch-Coleman-Li), the reference method for
    /// *bounded* least squares and the one `thick_opt`'s docs have been
    /// naming since the rewrite (`synthesis::trf`, R4.6). Hand-rolled, so it
    /// is always available; bounds enter the subproblem rather than clipping
    /// its answer.
    Trf,
}

impl OptimizerBackend {
    /// The name this backend answers to on the Python surface and in configs.
    pub fn as_str(self) -> &'static str {
        match self {
            OptimizerBackend::BuiltinLm => "builtin",
            OptimizerBackend::MinpackLm => "minpack_lm",
            OptimizerBackend::Trf => "trf",
        }
    }

    /// Parse a backend name. `Err` carries the list of names, because the
    /// caller that got this wrong is usually a user typing a string.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "builtin" => Ok(OptimizerBackend::BuiltinLm),
            "minpack_lm" => Ok(OptimizerBackend::MinpackLm),
            "trf" => Ok(OptimizerBackend::Trf),
            other => Err(format!(
                "unknown optimizer backend {other:?} (expected one of: builtin, minpack_lm, trf)"
            )),
        }
    }

    /// Whether this build can actually run it.
    pub fn is_available(self) -> bool {
        match self {
            OptimizerBackend::BuiltinLm => true,
            OptimizerBackend::MinpackLm => cfg!(feature = "opt-minpack-lm"),
            OptimizerBackend::Trf => true,
        }
    }

    /// Whether the backend honours `lb`/`ub` natively, or has to be wrapped
    /// in [`IntervalMap`] and therefore ends strictly inside the box.
    pub fn bounds_are_native(self) -> bool {
        matches!(self, OptimizerBackend::BuiltinLm | OptimizerBackend::Trf)
    }

    /// The cargo feature that supplies this backend, if any. `None` for the
    /// two that are always compiled in.
    ///
    /// Public because the Python binding builds its own rebuild hint and a
    /// second copy of this mapping is exactly the kind of pair that drifts.
    pub fn feature(self) -> Option<&'static str> {
        match self {
            OptimizerBackend::BuiltinLm => None,
            OptimizerBackend::MinpackLm => Some("opt-minpack-lm"),
            OptimizerBackend::Trf => None,
        }
    }

    fn unavailable(self) -> String {
        match self.feature() {
            Some(f) => format!(
                "optimizer backend {:?} is not compiled in — rebuild with \
                 `--features {}` (maturin: `maturin develop --release \
                 --features {}`)",
                self.as_str(),
                f,
                f
            ),
            None => format!("optimizer backend {:?} is unavailable", self.as_str()),
        }
    }
}

/// Every backend this crate knows about, in declaration order.
///
/// The list is the whole enum, not the compiled-in subset — pair it with
/// [`OptimizerBackend::is_available`] to answer "what can this build run?",
/// which is the question a caller has before it picks one.
pub const ALL_BACKENDS: [OptimizerBackend; 3] =
    [OptimizerBackend::BuiltinLm, OptimizerBackend::MinpackLm, OptimizerBackend::Trf];

/// The names this build can actually run. `["builtin"]` on a standard build.
pub fn available_backends() -> Vec<&'static str> {
    ALL_BACKENDS.iter().filter(|b| b.is_available()).map(|b| b.as_str()).collect()
}

/// What every backend returns.
///
/// The shape is the built-in LM's, because that is what the synthesis loop
/// already consumes. Two fields are LM-specific and carry what a foreign
/// backend can honestly say: `gain_ratio` is `NaN` from a solver with no such
/// notion, and `analytic_jacobians` counts Jacobians served by the analytic
/// source on *any* backend, since that is a property of the problem, not the
/// solver.
#[derive(Clone, Debug)]
pub struct OptimizerResult {
    pub x: Vec<f64>,
    pub cost: f64,
    pub iterations: usize,
    pub evals: usize,
    pub termination: LmTermination,
    pub gain_ratio: f64,
    pub analytic_jacobians: usize,
    /// Which backend produced this. The same struct comes back from all of
    /// them, so the answer has to travel with it.
    pub backend: OptimizerBackend,
}

impl OptimizerResult {
    fn from_lm(r: LmResult, backend: OptimizerBackend) -> Self {
        OptimizerResult {
            x: r.x,
            cost: r.cost,
            iterations: r.iterations,
            evals: r.evals,
            termination: r.termination,
            gain_ratio: r.gain_ratio,
            analytic_jacobians: r.analytic_jacobians,
            backend,
        }
    }
}

// ---------------------------------------------------------------------------
// Interior reparametrization
// ---------------------------------------------------------------------------

/// `x = mid + half·tanh(u)`: an unbounded `u` for every backend that has no
/// bounds of its own.
///
/// The map is elementwise, so an analytic Jacobian passes through it by
/// scaling column *k* by `dx_k/du_k` and a differenced one needs nothing at
/// all — `u` is just where the differences are taken.
///
/// Two properties worth stating because they are the cost of using it:
///
/// * **The optimum is strictly interior.** `tanh` never reaches ±1, so a
///   parameter whose true optimum is the bound converges to it only in the
///   limit.
/// * **The gradient vanishes as the bound is approached.** `dx/du = half·
///   sech²(u)` → 0, so a parameter pressed against a bound stalls rather than
///   stopping. On a merit whose optimum is interior — which is the case the
///   plan wraps these backends for — neither matters.
#[derive(Clone, Debug)]
pub struct IntervalMap {
    mid: Vec<f64>,
    half: Vec<f64>,
}

/// `tanh(u)` is clamped to this before `atanh`, so an `x0` sitting exactly on
/// a bound maps to a large finite `u` instead of ±∞. The resulting `x` is
/// inside the bound by `half · 1e-14`, which is below the width of any
/// thickness the synthesis loop can represent.
const ATANH_CLAMP: f64 = 1.0 - 1e-14;

impl IntervalMap {
    /// Build the map for `[lb, ub]`. Bounds must be finite; a degenerate
    /// `lb == ub` pins that parameter (slope 0), which is the only sensible
    /// reading of "optimize a variable that cannot move".
    pub fn new(lb: &[f64], ub: &[f64]) -> Result<Self, String> {
        if lb.len() != ub.len() {
            return Err("IntervalMap: bound length mismatch".into());
        }
        let mut mid = Vec::with_capacity(lb.len());
        let mut half = Vec::with_capacity(lb.len());
        for (k, (&l, &u)) in lb.iter().zip(ub.iter()).enumerate() {
            if !l.is_finite() || !u.is_finite() {
                return Err(format!(
                    "IntervalMap: bound {k} is not finite ([{l}, {u}]) — an \
                     unbounded backend needs a finite box to be reparametrized \
                     into"
                ));
            }
            if l > u {
                return Err(format!("IntervalMap: lb[{k}] > ub[{k}]"));
            }
            mid.push(0.5 * (l + u));
            half.push(0.5 * (u - l));
        }
        Ok(IntervalMap { mid, half })
    }

    pub fn len(&self) -> usize {
        self.mid.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mid.is_empty()
    }

    /// `u → x`, always inside the box.
    pub fn to_bounded(&self, u: &[f64]) -> Vec<f64> {
        u.iter()
            .enumerate()
            .map(|(k, &uk)| self.mid[k] + self.half[k] * uk.tanh())
            .collect()
    }

    /// `x → u`. `x` outside the box saturates rather than producing NaN.
    pub fn to_unbounded(&self, x: &[f64]) -> Vec<f64> {
        x.iter()
            .enumerate()
            .map(|(k, &xk)| {
                if self.half[k] == 0.0 {
                    0.0
                } else {
                    let t = ((xk - self.mid[k]) / self.half[k]).clamp(-ATANH_CLAMP, ATANH_CLAMP);
                    t.atanh()
                }
            })
            .collect()
    }

    /// `dx_k/du_k` at `u`.
    pub fn slope(&self, k: usize, uk: f64) -> f64 {
        let c = uk.cosh();
        self.half[k] / (c * c)
    }
}

/// A [`JacobianSource`] in `x`, seen from `u`: same rows, each column scaled
/// by `dx_k/du_k`.
///
/// Compiled in only where something uses it: the unbounded backends, and the
/// tests that pin it against differences taken in `u` directly.
#[cfg(any(test, feature = "opt-minpack-lm"))]
struct MappedJacobian<'a, J: ?Sized> {
    inner: &'a J,
    map: &'a IntervalMap,
}

#[cfg(any(test, feature = "opt-minpack-lm"))]
impl<J: JacobianSource + ?Sized> JacobianSource for MappedJacobian<'_, J> {
    fn fill(&self, u: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        let x = self.map.to_bounded(u);
        let m = match self.inner.fill(&x, jac)? {
            Some(m) => m,
            None => return Ok(None),
        };
        let n = u.len();
        for (k, &uk) in u.iter().enumerate() {
            let s = self.map.slope(k, uk);
            for i in 0..m {
                jac[i * n + k] *= s;
            }
        }
        Ok(Some(m))
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Minimize ½‖r(x)‖² over `lb ≤ x ≤ ub` on the backend `cfg.backend` names.
///
/// `residuals` and `jacobian` are the same pair [`levenberg_marquardt_with`]
/// takes; a backend without bounds of its own receives them composed with
/// [`IntervalMap`] and never sees an infeasible `x`.
pub fn run_optimizer<F, J>(
    residuals: &F,
    jacobian: Option<&J>,
    x0: &[f64],
    lb: &[f64],
    ub: &[f64],
    cfg: &LmConfig,
) -> Result<OptimizerResult, String>
where
    F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
    J: JacobianSource + ?Sized,
{
    let backend = cfg.backend;
    if !backend.is_available() {
        return Err(backend.unavailable());
    }
    match backend {
        OptimizerBackend::BuiltinLm => {
            levenberg_marquardt_with(residuals, jacobian, x0, lb, ub, cfg)
                .map(|r| OptimizerResult::from_lm(r, backend))
        },
        OptimizerBackend::Trf => {
            crate::smatrix::synthesis::trf::trust_region_reflective(
                residuals, jacobian, x0, lb, ub, cfg,
            )
            .map(|r| OptimizerResult::from_lm(r, backend))
        },
        OptimizerBackend::MinpackLm => {
            #[cfg(feature = "opt-minpack-lm")]
            {
                minpack::run(residuals, jacobian, x0, lb, ub, cfg)
            }
            #[cfg(not(feature = "opt-minpack-lm"))]
            {
                let _ = (residuals, jacobian, x0, lb, ub);
                Err(backend.unavailable())
            }
        },
    }
}

// ---------------------------------------------------------------------------
// MINPACK-derived reference LM (feature `opt-minpack-lm`)
// ---------------------------------------------------------------------------

#[cfg(feature = "opt-minpack-lm")]
mod minpack {
    use super::super::thick_opt::build_jacobian;
    use super::*;
    use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt, TerminationReason};
    use nalgebra::storage::Owned;
    use nalgebra::{DMatrix, DVector, Dyn};
    use std::cell::{Cell, RefCell};

    /// The `levenberg-marquardt` crate's view of our residual pair.
    ///
    /// The crate has no error channel — `residuals()`/`jacobian()` return
    /// `Option` and a `None` becomes `TerminationReason::User` — so a failure
    /// is stashed here and re-raised after `minimize` returns. Losing the
    /// message would turn "the merit spec is missing a curve" into "the
    /// solver gave up", which are not the same fact.
    pub(super) struct Problem<'a> {
        residuals: &'a (dyn Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync + 'a),
        jacobian: Option<&'a dyn JacobianSource>,
        x: DVector<f64>,
        /// (residual calls, Jacobian builds, Jacobians served analytically)
        counts: Cell<(usize, usize, usize)>,
        err: RefCell<Option<String>>,
    }

    impl Problem<'_> {
        fn fail(&self, e: String) {
            let mut slot = self.err.borrow_mut();
            if slot.is_none() {
                *slot = Some(e);
            }
        }
    }

    impl LeastSquaresProblem<f64, Dyn, Dyn> for Problem<'_> {
        type ResidualStorage = Owned<f64, Dyn>;
        type JacobianStorage = Owned<f64, Dyn, Dyn>;
        type ParameterStorage = Owned<f64, Dyn>;

        fn set_params(&mut self, x: &DVector<f64>) {
            self.x.copy_from(x);
        }

        fn params(&self) -> DVector<f64> {
            self.x.clone()
        }

        fn residuals(&self) -> Option<DVector<f64>> {
            let mut out = Vec::new();
            if let Err(e) = (self.residuals)(self.x.as_slice(), &mut out) {
                self.fail(e);
                return None;
            }
            let (r, j, a) = self.counts.get();
            self.counts.set((r + 1, j, a));
            Some(DVector::from_vec(out))
        }

        fn jacobian(&self) -> Option<DMatrix<f64>> {
            let n = self.x.len();
            let mut jac: Vec<f64> = Vec::new();

            let (m, analytic) = match self.jacobian.map(|s| s.fill(self.x.as_slice(), &mut jac)) {
                Some(Err(e)) => {
                    self.fail(e);
                    return None;
                },
                Some(Ok(Some(m))) => (m, true),
                // No source, or a source that declined: difference it here.
                // The crate never differences on its own — its
                // `differentiate_numerically` is a checker, not a fallback.
                _ => match build_jacobian(&self.residuals, self.x.as_slice(), &mut jac) {
                    Ok(added) => {
                        let (r, j, a) = self.counts.get();
                        self.counts.set((r + added, j, a));
                        (jac.len() / n, false)
                    },
                    Err(e) => {
                        self.fail(e);
                        return None;
                    },
                },
            };

            let (r, j, a) = self.counts.get();
            self.counts.set((r, j + 1, a + usize::from(analytic)));
            // `jac` is row-major m×n; nalgebra is column-major.
            Some(DMatrix::from_fn(m, n, |i, k| jac[i * n + k]))
        }
    }

    fn map_termination(t: TerminationReason) -> Result<LmTermination, String> {
        match t {
            TerminationReason::Orthogonal => Ok(LmTermination::Gradient),
            TerminationReason::ResidualsZero => Ok(LmTermination::Cost),
            // MINPACK reports both tests in one variant; ftol is the stronger
            // statement, so it wins when both fired.
            TerminationReason::Converged { ftol: true, .. } => Ok(LmTermination::Cost),
            TerminationReason::Converged { xtol: true, .. } => Ok(LmTermination::Step),
            TerminationReason::Converged { .. } => Ok(LmTermination::Cost),
            TerminationReason::LostPatience => Ok(LmTermination::MaxIterations),
            TerminationReason::NoImprovementPossible(_) => Ok(LmTermination::Stalled),
            TerminationReason::User(w) => Err(format!("minpack_lm: problem reported failure ({w})")),
            TerminationReason::Numerical(w) => {
                Err(format!("minpack_lm: non-finite value in {w}"))
            },
            TerminationReason::WrongDimensions(w) => {
                Err(format!("minpack_lm: inconsistent problem shape ({w})"))
            },
            TerminationReason::NoParameters => Err("minpack_lm: empty parameter vector".into()),
            TerminationReason::NoResiduals => Err("minpack_lm: empty residual vector".into()),
        }
    }

    pub(super) fn run<F, J>(
        residuals: &F,
        jacobian: Option<&J>,
        x0: &[f64],
        lb: &[f64],
        ub: &[f64],
        cfg: &LmConfig,
    ) -> Result<OptimizerResult, String>
    where
        F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
        J: JacobianSource + ?Sized,
    {
        let n = x0.len();
        if n == 0 {
            return Err("minpack_lm: empty parameter vector".into());
        }
        if lb.len() != n || ub.len() != n {
            return Err("minpack_lm: bound length mismatch".into());
        }
        // The crate is unbounded, so the whole run happens in `u`. Every
        // piece of that is named and separately tested: the map, the chain
        // rule on the analytic Jacobian, and — for the differenced path —
        // nothing at all, since differences taken in `u` need no correction.
        let map = IntervalMap::new(lb, ub)?;
        let u0 = map.to_unbounded(x0);
        let mapped_r =
            |u: &[f64], out: &mut Vec<f64>| -> Result<(), String> { residuals(&map.to_bounded(u), out) };
        let mapped_j = jacobian.map(|j| MappedJacobian { inner: j, map: &map });

        let problem = Problem {
            residuals: &mapped_r,
            jacobian: mapped_j.as_ref().map(|j| j as &dyn JacobianSource),
            x: DVector::from_row_slice(&u0),
            counts: Cell::new((0, 0, 0)),
            err: RefCell::new(None),
        };
        let driver = LevenbergMarquardt::new()
            .with_ftol(cfg.ftol)
            .with_xtol(cfg.xtol)
            .with_gtol(cfg.gtol)
            .with_patience(cfg.max_evals);
        let (problem, report) = driver.minimize(problem);

        if let Some(e) = problem.err.borrow_mut().take() {
            return Err(e);
        }
        let termination = map_termination(report.termination)?;
        let (evals, iterations, analytic) = problem.counts.get();
        Ok(OptimizerResult {
            x: map.to_bounded(problem.x.as_slice()),
            // The crate reports ½‖r‖²; Navette's cost is ‖r‖².
            cost: 2.0 * report.objective_function,
            // The crate reports evaluations, not iterations. MINPACK builds
            // the Jacobian exactly once per outer iteration, so counting
            // builds *is* the iteration count.
            iterations,
            evals,
            termination,
            gain_ratio: f64::NAN,
            analytic_jacobians: analytic,
            backend: OptimizerBackend::MinpackLm,
        })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::thick_opt::{build_jacobian, NoJacobian};
    use super::*;

    /// r_i(x) = a_i·x_0 + b_i·x_1 − y_i — an exactly linear least-squares
    /// problem with a unique interior optimum at (1, −2).
    fn linear(x: &[f64], out: &mut Vec<f64>) -> Result<(), String> {
        out.clear();
        for i in 0..6 {
            let a = 1.0 + i as f64;
            let b = 2.0 - 0.3 * i as f64;
            let y = a * 1.0 + b * -2.0;
            out.push(a * x[0] + b * x[1] - y);
        }
        Ok(())
    }

    struct LinearJac;

    impl JacobianSource for LinearJac {
        fn fill(&self, _x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
            jac.clear();
            for i in 0..6 {
                jac.push(1.0 + i as f64);
                jac.push(2.0 - 0.3 * i as f64);
            }
            Ok(Some(6))
        }
    }

    // -- backend names ------------------------------------------------------

    #[test]
    fn every_backend_name_round_trips() {
        for b in [OptimizerBackend::BuiltinLm, OptimizerBackend::MinpackLm] {
            assert_eq!(OptimizerBackend::parse(b.as_str()), Ok(b));
        }
    }

    #[test]
    fn the_available_list_always_contains_the_default_and_only_real_entries() {
        let av = available_backends();
        assert!(av.contains(&"builtin"));
        for name in &av {
            assert!(OptimizerBackend::parse(name).unwrap().is_available(), "{name}");
        }
        assert_eq!(av.contains(&"minpack_lm"), cfg!(feature = "opt-minpack-lm"));
    }

    #[test]
    fn an_unknown_backend_name_lists_the_known_ones() {
        let e = OptimizerBackend::parse("scipy").unwrap_err();
        assert!(e.contains("builtin") && e.contains("minpack_lm"), "got: {e}");
    }

    #[test]
    fn a_missing_backend_says_how_to_get_it_rather_than_falling_back() {
        // The message is the whole point of declaring a variant the build
        // cannot run, so it is checked either way.
        let b = OptimizerBackend::MinpackLm;
        let msg = b.unavailable();
        assert!(msg.contains("opt-minpack-lm"), "got: {msg}");
        // And on a build without it, the dispatcher refuses rather than
        // quietly handing the problem to the built-in.
        if !b.is_available() {
            let cfg = LmConfig { backend: b, ..Default::default() };
            let err = run_optimizer(
                &linear,
                None::<&NoJacobian>,
                &[0.0, 0.0],
                &[-5.0, -5.0],
                &[5.0, 5.0],
                &cfg,
            )
            .unwrap_err();
            assert!(err.contains("opt-minpack-lm"), "got: {err}");
        }
    }

    // -- the interior reparametrization -------------------------------------

    #[test]
    fn the_interval_map_round_trips_interior_points() {
        let map = IntervalMap::new(&[0.0, -3.0], &[100.0, 7.0]).unwrap();
        let x = [12.5, 2.25];
        let back = map.to_bounded(&map.to_unbounded(&x));
        for (a, b) in back.iter().zip(x.iter()) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }

    #[test]
    fn a_point_on_the_bound_maps_to_a_finite_u_just_inside() {
        let map = IntervalMap::new(&[0.0], &[100.0]).unwrap();
        for &x in &[0.0, 100.0] {
            let u = map.to_unbounded(&[x]);
            assert!(u[0].is_finite(), "u = {} at x = {x}", u[0]);
            let back = map.to_bounded(&u)[0];
            assert!(back > 0.0 && back < 100.0, "x={x} came back as {back}");
            assert!((back - x).abs() < 1e-9, "x={x} came back as {back}");
        }
    }

    #[test]
    fn the_map_never_leaves_the_box_however_large_u_is() {
        let map = IntervalMap::new(&[2.0], &[5.0]).unwrap();
        for &u in &[-1e6, -40.0, 0.0, 40.0, 1e6] {
            let x = map.to_bounded(&[u])[0];
            assert!((2.0..=5.0).contains(&x), "u={u} gave x={x}");
        }
    }

    #[test]
    fn the_slope_is_the_derivative_of_the_map() {
        let map = IntervalMap::new(&[0.0, -1.0], &[10.0, 1.0]).unwrap();
        for &u in &[-2.0, -0.3, 0.0, 0.7, 3.0] {
            for k in 0..2 {
                let h = 1e-6;
                let mut up = vec![0.0; 2];
                let mut um = vec![0.0; 2];
                up[k] = u + h;
                um[k] = u - h;
                let fd = (map.to_bounded(&up)[k] - map.to_bounded(&um)[k]) / (2.0 * h);
                let an = map.slope(k, u);
                assert!((fd - an).abs() < 1e-6 * an.abs().max(1e-3), "k={k} u={u}: {fd} vs {an}");
            }
        }
    }

    #[test]
    fn a_degenerate_interval_pins_its_parameter() {
        let map = IntervalMap::new(&[4.0], &[4.0]).unwrap();
        assert_eq!(map.to_unbounded(&[4.0]), vec![0.0]);
        assert_eq!(map.to_bounded(&[123.0]), vec![4.0]);
        assert_eq!(map.slope(0, 123.0), 0.0);
    }

    #[test]
    fn an_infinite_bound_is_refused_with_the_reason() {
        let e = IntervalMap::new(&[0.0], &[f64::INFINITY]).unwrap_err();
        assert!(e.contains("finite"), "got: {e}");
    }

    #[test]
    fn the_mapped_jacobian_is_the_chain_rule() {
        let map = IntervalMap::new(&[-4.0, -4.0], &[4.0, 4.0]).unwrap();
        let mapped = MappedJacobian { inner: &LinearJac, map: &map };
        let u = [0.4, -0.9];
        let mut jac = Vec::new();
        let m = mapped.fill(&u, &mut jac).unwrap().unwrap();
        assert_eq!(m, 6);

        // Differenced directly in u-space — the same thing by construction,
        // which is exactly why the analytic path has to agree with it.
        let ru = |u: &[f64], out: &mut Vec<f64>| linear(&map.to_bounded(u), out);
        let mut fd = Vec::new();
        build_jacobian(&ru, &u, &mut fd).unwrap();
        for (a, b) in jac.iter().zip(fd.iter()) {
            assert!((a - b).abs() < 1e-6 * a.abs().max(1.0), "{a} vs {b}");
        }
    }

    // -- dispatch -----------------------------------------------------------

    #[test]
    fn the_builtin_backend_is_the_default_and_is_the_lm_itself() {
        let cfg = LmConfig::default();
        assert_eq!(cfg.backend, OptimizerBackend::BuiltinLm);
        let x0 = [0.0, 0.0];
        let lb = [-5.0, -5.0];
        let ub = [5.0, 5.0];
        let direct = levenberg_marquardt_with(&linear, Some(&LinearJac), &x0, &lb, &ub, &cfg).unwrap();
        let via = run_optimizer(&linear, Some(&LinearJac), &x0, &lb, &ub, &cfg).unwrap();
        assert_eq!(via.backend, OptimizerBackend::BuiltinLm);
        assert_eq!(via.x, direct.x);
        assert_eq!(via.cost, direct.cost);
        assert_eq!(via.evals, direct.evals);
        assert_eq!(via.iterations, direct.iterations);
    }

    #[cfg(feature = "opt-minpack-lm")]
    mod minpack_backend {
        use super::*;

        fn cfg() -> LmConfig {
            LmConfig { backend: OptimizerBackend::MinpackLm, ..Default::default() }
        }

        #[test]
        fn it_finds_the_same_interior_optimum_as_the_builtin() {
            let x0 = [0.0, 0.0];
            let lb = [-5.0, -5.0];
            let ub = [5.0, 5.0];
            let a = run_optimizer(&linear, Some(&LinearJac), &x0, &lb, &ub, &LmConfig::default())
                .unwrap();
            let b = run_optimizer(&linear, Some(&LinearJac), &x0, &lb, &ub, &cfg()).unwrap();
            assert_eq!(b.backend, OptimizerBackend::MinpackLm);
            for (p, q) in a.x.iter().zip(b.x.iter()) {
                assert!((p - q).abs() < 1e-7, "{a:?} vs {b:?}");
            }
            assert!((b.x[0] - 1.0).abs() < 1e-7 && (b.x[1] + 2.0).abs() < 1e-7, "{:?}", b.x);
            assert!(b.cost < 1e-18, "cost {}", b.cost);
        }

        #[test]
        fn it_reports_the_analytic_jacobians_it_was_given() {
            let x0 = [0.0, 0.0];
            let lb = [-5.0, -5.0];
            let ub = [5.0, 5.0];
            let an = run_optimizer(&linear, Some(&LinearJac), &x0, &lb, &ub, &cfg()).unwrap();
            let fd = run_optimizer(&linear, None::<&NoJacobian>, &x0, &lb, &ub, &cfg()).unwrap();
            assert!(an.analytic_jacobians > 0, "{an:?}");
            assert_eq!(fd.analytic_jacobians, 0);
            assert_eq!(an.iterations, an.analytic_jacobians);
            // Differencing costs 2n evaluations per Jacobian; the analytic
            // source costs none. Same problem, so the gap is the accounting.
            assert!(fd.evals > an.evals, "{} vs {}", fd.evals, an.evals);
        }

        #[test]
        fn the_optimum_stays_inside_the_box_when_the_true_one_is_outside() {
            // The unconstrained optimum is (1, -2); squeeze the box so both
            // components want to leave it. Interior parametrization means
            // "arbitrarily close to the bound", never past it.
            let lb = [1.5, -1.0];
            let ub = [4.0, 3.0];
            let r = run_optimizer(&linear, Some(&LinearJac), &[2.0, 0.0], &lb, &ub, &cfg()).unwrap();
            for (k, &v) in r.x.iter().enumerate() {
                assert!(v > lb[k] && v < ub[k], "x[{k}]={v} left [{}, {}]", lb[k], ub[k]);
            }
            assert!((r.x[0] - 1.5).abs() < 1e-3, "expected the bound, got {}", r.x[0]);
        }

        #[test]
        fn a_failing_residual_reports_its_own_message_not_the_solvers() {
            let boom = |_x: &[f64], _o: &mut Vec<f64>| -> Result<(), String> {
                Err("missing curve Rs".into())
            };
            let e = run_optimizer(&boom, None::<&NoJacobian>, &[0.0], &[-1.0], &[1.0], &cfg())
                .unwrap_err();
            assert!(e.contains("missing curve Rs"), "got: {e}");
        }

        #[test]
        fn the_cost_is_navettes_convention_not_minpacks() {
            // Start at the optimum of a problem whose residual is a known
            // constant vector, so ‖r‖² is known by hand: r = (3, 4) at any x
            // ⇒ cost 25, and MINPACK's own objective would be 12.5.
            let flat = |_x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
                out.clear();
                out.push(3.0);
                out.push(4.0);
                Ok(())
            };
            let r = run_optimizer(&flat, None::<&NoJacobian>, &[0.5], &[0.0], &[1.0], &cfg()).unwrap();
            assert!((r.cost - 25.0).abs() < 1e-9, "cost {}", r.cost);
        }
    }
}
