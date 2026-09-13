// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::config — PipelineConfig / TerminationReason / cycle knobs.
//!
//! Verbatim port of needle_pipeline.py section 1–2 (defaults included).

// ---------------------------------------------------------------------------
// TerminationReason
// ---------------------------------------------------------------------------

/// Why the pipeline stopped — mirrors Python `TerminationReason`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TerminationReason {
    LayerBudgetReached,
    ThicknessBudgetReached,
    StagnationPlateau,
    StagnationOscillation,
    StagnationDivergence,
    MaxIterationsReached,
    MeritTargetReached,
    UserAbort,
}

impl TerminationReason {
    pub fn name(self) -> &'static str {
        match self {
            TerminationReason::LayerBudgetReached => "LAYER_BUDGET_REACHED",
            TerminationReason::ThicknessBudgetReached => "THICKNESS_BUDGET_REACHED",
            TerminationReason::StagnationPlateau => "STAGNATION_PLATEAU",
            TerminationReason::StagnationOscillation => "STAGNATION_OSCILLATION",
            TerminationReason::StagnationDivergence => "STAGNATION_DIVERGENCE",
            TerminationReason::MaxIterationsReached => "MAX_ITERATIONS_REACHED",
            TerminationReason::MeritTargetReached => "MERIT_TARGET_REACHED",
            TerminationReason::UserAbort => "USER_ABORT",
        }
    }
}

// ---------------------------------------------------------------------------
// PipelineConfig
// ---------------------------------------------------------------------------

/// What the floor does to a sub-minimum layer (F0.3, U1).
///
/// "Remove a layer that got too thin" and "set that layer to the minimum
/// thickness you can actually deposit" are two different operations with
/// two different purposes. `Remove` is today's behaviour and what the
/// fingerprints hold.
///
/// `clamp_min_nm` genuinely does two jobs and this enum selects which one
/// it does: an *elimination threshold* under `Remove` (the optimizer
/// driving a film to zero is the optimizer saying it wants that layer
/// gone), and a *manufacturing floor* under the clamp-up policies. The
/// config key keeps its name; this comment is where both jobs are stated.
///
/// Span handling (F0.3's deferral, landed in three steps): a scalable
/// multi-row span is SCALED to `clamp_min_nm`, fractions preserved
/// (F1.6's scale operation, the deferral's stated future). A plain
/// singleton-bulk span clamps up with its interface slice untouched
/// (C2: the slice belongs to the authored thickness, so the bulk row
/// takes `clamp_min_nm - t_slice`). Every other multi-row span is
/// removed whole with the F0.2 report under EVERY policy - a profile
/// the user did not offer for scaling has no principled clamp-up.
/// Under `ClampUpAlways` the LM bound binds the BULK row (the
/// parameter), so an interface-carrying film settles one slice
/// thickness above the floor: `bulk >= clamp_min_nm`, total
/// `>= clamp_min_nm + t_slice` - conservative by construction, not a
/// bug (C2/V3).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ThinLayerPolicy {
    /// Sub-minimum layers are removed, during the run and at the end.
    /// The search is exactly today's, elimination and all.
    #[default]
    Remove,
    /// The search runs exactly as today; only the run's final clamp pass
    /// sets a surviving sub-minimum layer to `clamp_min_nm` instead of
    /// removing it. A film the optimizer genuinely wanted gone is already
    /// gone by then; a film that survived and happens to sit at 0.8 nm
    /// comes out depositable. This is the variant the documentation leads
    /// with.
    ClampUpFinal,
    /// The floor sets films to `clamp_min_nm` during the run too, and the
    /// LM lower bound moves with it (`lb = clamp_min_nm`) - without the
    /// bound, LM drives a film to 0.5 nm, the clamp puts it back to 2.0,
    /// and the pair oscillates until the stagnation detector terminates
    /// the run with a true report of a false condition. Refused at
    /// `NeedlePipeline::new` when `needles_per_cycle > 0`: with the floor
    /// as a hard bound nothing can ever be eliminated, so a bad needle
    /// seed parks at the floor permanently and layer count only ever
    /// grows. The mode for re-optimizing a fixed architecture.
    ClampUpAlways,
}

impl ThinLayerPolicy {
    /// The Python-facing spelling (see `PyPipelineConfig`).
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "remove" => Ok(Self::Remove),
            "clamp_up_final" => Ok(Self::ClampUpFinal),
            "clamp_up_always" => Ok(Self::ClampUpAlways),
            other => Err(format!(
                "thin_layer_policy must be 'remove', 'clamp_up_final' or \
                 'clamp_up_always', got {other:?}"
            )),
        }
    }
}

/// All loop-control parameters — mirrors Python `PipelineConfig`.
#[derive(Clone, Debug)]
pub struct PipelineConfig {
    // -- budgets --
    pub max_film_layers: usize,
    pub max_total_thickness_nm: f64,
    pub max_macro_cycles: usize,
    pub merit_target: f64,

    // -- clamping --
    pub clamp_min_nm: f64,
    pub clamp_max_nm: f64,
    /// F0.3 (U1): what the floor does to a sub-minimum layer. `Remove`
    /// reproduces today bit for bit; the clamp-up variants and their
    /// couplings are documented on the enum.
    pub thin_layer_policy: ThinLayerPolicy,

    // -- needle --
    pub needles_per_cycle: usize,

    // -- cleanup --
    pub enable_cleanup: bool,
    /// `None` → falls back to `clamp_min_nm` (`validated()` resolves it).
    pub cleanup_min_nm: Option<f64>,
    /// `None` → uncapped.
    pub cleanup_max_removals: Option<usize>,

    // -- inflate --
    pub enable_inflate: bool,
    pub inflate_addon_qwot: f64,
    pub inflate_reference_wl: f64,
    /// `None` → all layers.
    pub inflate_max_layers: Option<usize>,

    // -- stagnation --
    pub stagnation_window: usize,
    pub stagnation_gradient_tol: f64,
    pub stagnation_oscillation_ratio: f64,
    pub stagnation_divergence_count: usize,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        PipelineConfig {
            max_film_layers: 40,
            max_total_thickness_nm: 5000.0,
            max_macro_cycles: 50,
            merit_target: 0.0,
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            thin_layer_policy: ThinLayerPolicy::Remove,
            needles_per_cycle: 3,
            enable_cleanup: true,
            cleanup_min_nm: None, // resolved to clamp_min_nm by validated()
            cleanup_max_removals: None,
            enable_inflate: false,
            inflate_addon_qwot: 2.0,
            inflate_reference_wl: 550.0,
            inflate_max_layers: None,
            stagnation_window: 5,
            stagnation_gradient_tol: 1e-4,
            stagnation_oscillation_ratio: 0.75,
            stagnation_divergence_count: 3,
        }
    }
}

impl PipelineConfig {
    /// Apply the `__post_init__` semantics: resolve `cleanup_min_nm`
    /// fallback and validate invariants. Returns an owned, resolved copy.
    pub fn validated(mut self) -> Result<Self, String> {
        if self.cleanup_min_nm.is_none() {
            self.cleanup_min_nm = Some(self.clamp_min_nm);
        }
        if self.clamp_min_nm < 0.0 {
            return Err("clamp_min_nm must be non-negative.".into());
        }
        if self.clamp_max_nm <= self.clamp_min_nm {
            return Err("clamp_max_nm must be greater than clamp_min_nm.".into());
        }
        if self.stagnation_window < 2 {
            return Err("stagnation_window must be >= 2.".into());
        }
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f03_default_is_remove_and_parse_round_trips() {
        assert_eq!(ThinLayerPolicy::default(), ThinLayerPolicy::Remove);
        for (s, want) in [
            ("remove", ThinLayerPolicy::Remove),
            ("clamp_up_final", ThinLayerPolicy::ClampUpFinal),
            ("clamp_up_always", ThinLayerPolicy::ClampUpAlways),
        ] {
            assert_eq!(ThinLayerPolicy::parse(s).unwrap(), want);
        }
        assert!(ThinLayerPolicy::parse("Remove").is_err()); // case-sensitive
        assert!(ThinLayerPolicy::parse("scale").is_err());
    }

    #[test]
    fn defaults_match_python() {
        let c = PipelineConfig::default();
        assert_eq!(c.max_film_layers, 40);
        assert!((c.max_total_thickness_nm - 5000.0).abs() < 1e-12);
        assert_eq!(c.max_macro_cycles, 50);
        assert_eq!(c.merit_target, 0.0);
        assert!((c.clamp_min_nm - 2.0).abs() < 1e-12);
        assert!((c.clamp_max_nm - 1000.0).abs() < 1e-12);
        assert_eq!(c.needles_per_cycle, 3);
        assert!(c.enable_cleanup && !c.enable_inflate);
        assert!((c.inflate_addon_qwot - 2.0).abs() < 1e-12);
        assert!((c.inflate_reference_wl - 550.0).abs() < 1e-12);
        assert_eq!(c.stagnation_window, 5);
        assert!((c.stagnation_gradient_tol - 1e-4).abs() < 1e-15);
        assert!((c.stagnation_oscillation_ratio - 0.75).abs() < 1e-12);
        assert_eq!(c.stagnation_divergence_count, 3);
    }

    #[test]
    fn validated_resolves_cleanup_fallback_and_checks_invariants() {
        let c = PipelineConfig::default().validated().unwrap();
        assert!((c.cleanup_min_nm.unwrap() - c.clamp_min_nm).abs() < 1e-12);

        assert!(
            PipelineConfig {
                clamp_min_nm: -1.0,
                ..Default::default()
            }
            .validated()
            .is_err()
        );
        assert!(
            PipelineConfig {
                clamp_min_nm: 5.0,
                clamp_max_nm: 5.0,
                ..Default::default()
            }
            .validated()
            .is_err()
        );
        assert!(
            PipelineConfig {
                stagnation_window: 1,
                ..Default::default()
            }
            .validated()
            .is_err()
        );
    }
}
