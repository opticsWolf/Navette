# SPDX-License-Identifier: LGPL-3.0-or-later
# Type stubs for the compiled extension navette._smatrix
# (submodule of the aggregated navette._navette extension, built from
# rust/navette-py over the navette Rust crate).
#
# Wavelengths are nanometres, thicknesses nanometres, angles degrees unless a
# parameter says otherwise. The `n_stack_cache` buffers are wav-major flat:
# layer `j` at wavelength `i` lives at `i * n_layers + j`. `Solver.__init__` is
# the exception -- its `indices` is layer-major (`layer * n_wavs + wav`), which
# is what a C-contiguous `(n_layers, n_wavs)` array ravels to; the Rust side
# transposes it into the wav-major cache.
#
# Kept in sync with the PyO3 registration list by tools/check_pyi_sync.py.
from __future__ import annotations

from typing import Any, Callable, Mapping, Sequence

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]
IntArray = npt.NDArray[np.int32]
BoolArray = npt.NDArray[np.bool_]

# ---------------------------------------------------------------------------
# Request bits
# ---------------------------------------------------------------------------
# `requested` is a bitmask over these. The helper builders below assemble the
# usual combinations; OR them directly for anything else. The returned dict of
# `Solver.solve` / `core_engine` carries one entry per requested channel.

NREQ_P: int
NREQ_P_MB: int
NREQ_P_T: int
NREQ_P_A: int
NREQ_P_PHI: int
NREQ_P_MB_T: int
NREQ_P_MB_A: int
NREQ_P_TB: int
NREQ_P_RB: int
NREQ_P_AB: int
NREQ_P_MB_TB: int
NREQ_P_MB_RB: int
NREQ_P_MB_AB: int
NREQ_DPHI: int
NREQ_DGD: int
NREQ_DGDD: int
NREQ_DTOD: int
NREQ_DFOD: int

def solver_rt_request(pol: str) -> int: ...
def solver_ellipsometry_request(transmission: bool) -> int: ...
def solver_absorption_request() -> int: ...
def solver_amplitudes_request() -> int: ...
def solver_stokes_request(reflection: bool, transmission: bool) -> int: ...
def solver_dispersion_request(
    reflection: bool, transmission: bool, s_pol: bool, p_pol: bool
) -> int: ...
def solver_energy_conservation(
    rs: FloatArray, rp: FloatArray, ts: FloatArray, tp: FloatArray
) -> FloatArray:
    """`max(|1 - Rs - Ts|, |1 - Rp - Tp|)` per grid point.

    One residual per point, maxed *over* the two polarizations -- not one
    row per polarization. The four inputs must share a shape, 1-D `(n,)` or
    2-D `(n_angles, n_wavs)`, and the result comes back in that same shape.
    """
    ...

# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------

def w_function(q: complex, rough_type: int) -> complex:
    """Névot-Croce style interface attenuation for roughness model
    `rough_type`."""
    ...

def redheffer_product_complex_field(
    r_a_front: complex, t_a_back: complex, t_a_fwd: complex, r_a_back: complex,
    r_b_front: complex, t_b_back: complex, t_b_fwd: complex, r_b_back: complex,
) -> tuple[complex, complex, complex, complex]: ...
def redheffer_product_real(
    ra_rf: float, ra_tb: float, ra_tf: float, ra_rb: float,
    rb_rf: float, rb_tb: float, rb_tf: float, rb_rb: float,
) -> tuple[float, float, float, float]: ...
def redheffer_product_cross(
    a_cf: complex, a_db: complex, a_df: complex, a_cb: complex,
    b_cf: complex, b_db: complex, b_df: complex, b_cb: complex,
) -> tuple[complex, complex, complex, complex]: ...
def solve_coherent_block_fields(
    start_idx: int,
    end_idx: int,
    n_stack: ComplexArray,
    d_stack: FloatArray,
    rough_vals: FloatArray,
    rough_types: IntArray,
    lam: float,
    nsin_fi: complex,
    pol: int,
) -> tuple[complex, complex, complex, complex, float, float, float, float]:
    """`(r_front, t_back, t_fwd, r_back, R_front, T_back, T_fwd, R_back)`."""
    ...

def core_engine(
    wavls: FloatArray,
    sin_theta_arr: FloatArray,
    n_layers: int,
    n_stack_cache: FloatArray,
    thicknesses: FloatArray,
    incoherent_flags: IntArray,
    rough_types: IntArray,
    rough_vals: FloatArray,
    coherence_mode: int,
    requested: int,
) -> dict[str, Any]:
    """One-shot solve. `n_stack_cache` is the wav-major interleaved
    `[re, im]` index buffer. Prefer `Solver` when solving the same stack more
    than once — it keeps the index cache."""
    ...

def needle_engine(
    wavls: FloatArray,
    sin_theta_arr: FloatArray,
    n_layers: int,
    n_stack_cache: FloatArray,
    thicknesses: FloatArray,
    rough_types: IntArray,
    rough_vals: FloatArray,
    needle_n_per_wav: ComplexArray,
    z_grid: FloatArray,
    requested: int,
    incoherent_flags: IntArray | None = ...,
    targets_r: FloatArray | None = ...,
    weights_r: FloatArray | None = ...,
    targets_t: FloatArray | None = ...,
    weights_t: FloatArray | None = ...,
    targets_a: FloatArray | None = ...,
    weights_a: FloatArray | None = ...,
    targets_phi: FloatArray | None = ...,
    weights_phi: FloatArray | None = ...,
    targets_tb: FloatArray | None = ...,
    weights_tb: FloatArray | None = ...,
    targets_rb: FloatArray | None = ...,
    weights_rb: FloatArray | None = ...,
    targets_ab: FloatArray | None = ...,
    weights_ab: FloatArray | None = ...,
    grads_r: FloatArray | None = ...,
    grads_t: FloatArray | None = ...,
    start_idx: int = ...,
    end_idx: int | None = ...,
    channel: int = ...,
    calc_s: bool = ...,
    calc_p: bool = ...,
    host_mask: BoolArray | None = ...,
) -> dict[str, Any]: ...
def scan_landscape(
    n_stack: ComplexArray,
    thicknesses: FloatArray,
    rough_types: IntArray,
    rough_vals: FloatArray,
    lam: float,
    pol: int,
    real_min: float,
    real_max: float,
    imag_min: float,
    imag_max: float,
    points_real: int,
    points_imag: int,
) -> tuple[list[float], list[float], FloatArray]:
    """`(real_vals, imag_vals, landscape)` for an eigenmode search."""
    ...

def find_local_minima(
    landscape: FloatArray,
    real_vals: Sequence[float],
    imag_vals: Sequence[float],
    median_factor: float,
) -> list[tuple[float, float]]: ...
def nelder_mead(
    n_stack: ComplexArray,
    thicknesses: FloatArray,
    rough_types: IntArray,
    rough_vals: FloatArray,
    lam: float,
    pol: int,
    x0: tuple[float, float],
    step: float,
    tol: float,
    max_iter: int,
) -> tuple[float, float, float]:
    """`(n_eff_real, n_eff_imag, residual)`."""
    ...

def field_profile(
    n_stack: ComplexArray,
    thicknesses: FloatArray,
    rough_types: IntArray,
    rough_vals: FloatArray,
    lam: float,
    n_eff: complex,
    pol: int,
    points_per_layer: int,
) -> tuple[list[float], list[float], list[float], list[float], list[complex]]: ...

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

class Solver:
    """A stack, held so it can be solved repeatedly without rebuilding the
    index cache.

    `indices` is layer-major flat (`layer * n_wavs + wav`), i.e. exactly a
    C-contiguous `(n_layers, n_wavs)` array raveled. `angles` are degrees
    unless `angles_in_radians`.
    """

    def __init__(
        self,
        wavelengths: FloatArray,
        angles: FloatArray,
        indices: ComplexArray,
        n_layers: int,
        thicknesses: FloatArray | None = ...,
        incoherent_flags: IntArray | None = ...,
        roughness_types: IntArray | None = ...,
        roughness_values: FloatArray | None = ...,
        coherence_mode: int = ...,
        angles_in_radians: bool = ...,
    ) -> None: ...
    @property
    def n_angles(self) -> int: ...
    @property
    def n_wavs(self) -> int: ...
    def solve(self, requested: int) -> dict[str, Any]:
        """One entry per channel in the `requested` bitmask."""
        ...

    def landscape(
        self,
        real_range: tuple[float, float],
        imag_range: tuple[float, float],
        points_real: int,
        points_imag: int,
        pol: int,
        wavelength: float | None = ...,
        wav_index: int | None = ...,
    ) -> tuple[list[float], list[float], list[float]]: ...
    def refine_mode(
        self,
        guess: complex,
        pol: int,
        wavelength: float | None = ...,
        wav_index: int | None = ...,
        step: float = ...,
        tol: float = ...,
        max_iter: int = ...,
    ) -> tuple[complex, float]: ...
    def find_eigenmodes(
        self,
        real_range: tuple[float, float],
        imag_range: tuple[float, float],
        points: tuple[int, int] = ...,
        median_factor: float = ...,
        refine: bool = ...,
        pol: int = ...,
        wavelength: float | None = ...,
        wav_index: int | None = ...,
    ) -> list[complex]: ...
    def field_profile(
        self,
        n_eff: complex,
        pol: int,
        wavelength: float | None = ...,
        wav_index: int | None = ...,
        points_per_layer: int = ...,
    ) -> tuple[list[float], list[float], list[float], list[float], list[complex]]: ...
    def needle_gradient(
        self,
        needle_n_per_wav: ComplexArray,
        z_grid: FloatArray,
        requested: int,
        incoherent_flags: IntArray | None = ...,
        targets_r: FloatArray | None = ...,
        weights_r: FloatArray | None = ...,
        targets_t: FloatArray | None = ...,
        weights_t: FloatArray | None = ...,
        targets_a: FloatArray | None = ...,
        weights_a: FloatArray | None = ...,
        targets_phi: FloatArray | None = ...,
        weights_phi: FloatArray | None = ...,
        targets_tb: FloatArray | None = ...,
        weights_tb: FloatArray | None = ...,
        targets_rb: FloatArray | None = ...,
        weights_rb: FloatArray | None = ...,
        targets_ab: FloatArray | None = ...,
        weights_ab: FloatArray | None = ...,
        grads_r: FloatArray | None = ...,
        grads_t: FloatArray | None = ...,
        start_idx: int = ...,
        end_idx: int | None = ...,
        channel: int = ...,
        calc_s: bool = ...,
        calc_p: bool = ...,
        host_mask: BoolArray | None = ...,
        gain_shift_phi: float = ...,
    ) -> dict[str, Any]: ...

# ---------------------------------------------------------------------------
# Merit
# ---------------------------------------------------------------------------

class SimCurves:
    """Simulated curves on one (angle, wavelength) grid, keyed by curve id.

    `total_d` / `n_front` / `n_back` feed the differential-phase reference;
    leaving them at their defaults reproduces absolute phase bit-for-bit.
    """

    def __init__(
        self,
        angles: FloatArray,
        wavelengths: FloatArray,
        total_d: float = ...,
        n_front: float = ...,
        n_back: float = ...,
    ) -> None: ...
    def set_curve(self, curve_id: str, values: FloatArray) -> None: ...
    def set_complex(self, curve_id: str, values: ComplexArray) -> None: ...

class MeritSpec:
    """Compiled merit: keys (angle, curve) plus the targets hung off them."""

    def __init__(self) -> None: ...
    def add_key(self, angle: float, curve_id: str) -> int: ...
    def add_target(
        self,
        key_idx: int,
        wavelengths: FloatArray,
        normalized: FloatArray,
        tolerances: FloatArray,
        kind: str,
        transform: str,
        norm_factor: float,
        band: FloatArray | None = ...,
        phase: bool = ...,
        differential_passes: float | None = ...,
        weight: float = ...,
        count_norm: float | None = ...,
        integral: bool = ...,
    ) -> None: ...
    def merit(self, sim: SimCurves, missing_penalty: float) -> float: ...
    def residuals(self, sim: SimCurves) -> FloatArray: ...
    def n_residuals(self) -> int: ...

def compile_merit_spec(request_json: str) -> MeritSpec: ...
def reference_rotation(
    wavelengths: FloatArray,
    angle_deg: float,
    n_inc: float = ...,
    total_d: float = ...,
    passes: float = ...,
) -> ComplexArray:
    """`exp(-i * passes * k * D * cos)` per wavelength — the differential-phase
    reference."""
    ...

def rotate_rows(rows: Sequence[complex], rot: Sequence[complex]) -> list[complex]: ...
def build_needle_targets(
    spec: MeritSpec,
    angles: FloatArray,
    wavelengths: FloatArray,
    sim: SimCurves | None = ...,
) -> dict[str, Any]:
    """Per-channel `targets`/`weights` surfaces for `needle_gradient`.

    Colour demands deposit into internal `grad_r`/`grad_t` buckets that this
    dict does **not** carry; a hand-assembled Python needle cycle therefore
    sees zero colour gradient. Use the `NeedlePipeline` / `run_design` path
    when the spec has colour demands.
    """
    ...

# ---------------------------------------------------------------------------
# Synthesis pipeline
# ---------------------------------------------------------------------------

class LayerSpec:
    def __init__(
        self,
        material: str,
        nk: ComplexArray,
        thickness: float,
        coherent: bool = ...,
        rough_type: int = ...,
        rough_val: float = ...,
        optimize: bool = ...,
        needle: bool = ...,
    ) -> None: ...
    @property
    def material(self) -> str: ...
    @property
    def nk(self) -> ComplexArray: ...
    @property
    def thickness(self) -> float: ...
    @property
    def coherent(self) -> bool: ...
    @property
    def rough_type(self) -> int: ...
    @property
    def rough_val(self) -> float: ...
    @property
    def optimize(self) -> bool: ...
    @property
    def needle(self) -> bool: ...

class DesignStack:
    """Ambient / films / substrate, with the film edits synthesis performs."""

    def __init__(
        self, ambient: LayerSpec, substrate: LayerSpec, films: Sequence[LayerSpec]
    ) -> None: ...
    @staticmethod
    def from_design(
        ambient: LayerSpec,
        substrate: LayerSpec,
        films: Sequence[Any],
        nk: Mapping[str, ComplexArray],
        groups: Mapping[str, Any],
        wavelengths: FloatArray,
        background: Sequence[str] | None = ...,
    ) -> DesignStack: ...
    @staticmethod
    def design_from_config(
        request_json: str, wavelengths: FloatArray
    ) -> tuple[DesignStack, dict[str, Any]]: ...
    def film_count(self) -> int: ...
    def total_thickness(self) -> float: ...
    def films(self) -> list[dict[str, Any]]: ...
    def to_dict(self) -> dict[str, Any]: ...
    def set_thickness(self, film_idx: int, thickness: float) -> None: ...
    def insert_needle_seed(
        self, film_idx: int, depth_into_layer_nm: float, seed: LayerSpec
    ) -> None: ...
    def merge_adjacent(self) -> int: ...
    def remove_film(self, film_idx: int) -> LayerSpec: ...
    def clamp_all(self, min_nm: float, max_nm: float) -> tuple[int, int]: ...

class LmConfig:
    """Levenberg-Marquardt settings. `optimizer` selects the backend;
    `available_optimizers()` lists what this build carries."""

    def __init__(
        self,
        max_iterations: int = ...,
        max_evals: int = ...,
        ftol: float = ...,
        xtol: float = ...,
        gtol: float = ...,
        lambda_init: float = ...,
        lambda_up: float = ...,
        lambda_down: float = ...,
        damping: str = ...,
        gtol_scale_invariant: bool = ...,
        jacobian: str = ...,
        optimizer: str = ...,
    ) -> None: ...

class PipelineConfig:
    def __init__(
        self,
        max_film_layers: int = ...,
        max_total_thickness_nm: float = ...,
        max_macro_cycles: int = ...,
        merit_target: float = ...,
        clamp_min_nm: float = ...,
        clamp_max_nm: float = ...,
        needles_per_cycle: int = ...,
        enable_cleanup: bool = ...,
        cleanup_min_nm: float | None = ...,
        cleanup_max_removals: int | None = ...,
        enable_inflate: bool = ...,
        inflate_addon_qwot: float = ...,
        inflate_reference_wl: float = ...,
        inflate_max_layers: int | None = ...,
        stagnation_window: int = ...,
        stagnation_gradient_tol: float = ...,
        stagnation_oscillation_ratio: float = ...,
        stagnation_divergence_count: int = ...,
    ) -> None: ...

class NeedleCycleConfig:
    def __init__(
        self,
        max_needles: int = ...,
        convergence_threshold: float = ...,
        needle_seed_thickness_nm: float = ...,
        scan_step_nm: float = ...,
        refold_per_cycle: bool = ...,
    ) -> None: ...

class SmatrixContext:
    """A merit spec bound to a grid — simulate, score, or optimize a stack
    against it."""

    def __init__(
        self,
        spec: MeritSpec,
        angles_deg: FloatArray,
        wavelengths: FloatArray,
        clamp_min: float = ...,
        clamp_max: float = ...,
        lm: LmConfig | None = ...,
    ) -> None: ...
    def simulate(self, stack: DesignStack) -> SimCurves: ...
    def evaluate_merit(self, stack: DesignStack) -> float: ...
    def optimize_thicknesses(self, stack: DesignStack) -> float: ...
    def optimize_thicknesses_report(
        self, stack: DesignStack
    ) -> tuple[float, dict[str, Any] | None]: ...

class NeedlePipeline:
    """Needle synthesis over one contrast pair. `run` reports through
    `callback`, which may abort by raising."""

    def __init__(
        self,
        stack: DesignStack,
        spec: MeritSpec,
        angles_deg: FloatArray,
        wavelengths: FloatArray,
        contrast: Mapping[str, LayerSpec],
        pipeline_config: PipelineConfig | None = ...,
        needle_config: NeedleCycleConfig | None = ...,
        lm: LmConfig | None = ...,
    ) -> None: ...
    def run(self, callback: Callable[..., Any] | None = ...) -> dict[str, Any]: ...

def available_optimizers() -> list[str]:
    """Backends this build carries: always `builtin` and `trf`, plus
    `minpack_lm` / `argmin_*` when their cargo features are on."""
    ...

def assemble_design(
    ambient_nk: Sequence[complex],
    ambient_name: str,
    substrate_nk: Sequence[complex],
    substrate_name: str,
    films: Sequence[Any],
    groups: Mapping[str, Any],
    seeds: Sequence[tuple[str, str, Sequence[complex]]],
    wavelengths: FloatArray,
) -> tuple[DesignStack, dict[str, Any]]: ...
def run_design(
    ambient_nk: Sequence[complex],
    ambient_name: str,
    substrate_nk: Sequence[complex],
    substrate_name: str,
    films: Sequence[Any],
    groups: Mapping[str, Any],
    seeds: Sequence[tuple[str, str, Sequence[complex]]],
    wavelengths: FloatArray,
    angles_deg: FloatArray,
    spec: MeritSpec,
    pipeline_config: PipelineConfig | None = ...,
    needle_config: NeedleCycleConfig | None = ...,
    lm: LmConfig | None = ...,
    callback: Callable[..., Any] | None = ...,
) -> dict[str, Any]: ...
