# SPDX-License-Identifier: LGPL-3.0-or-later
# Type stubs for the compiled extension navette._spectralweave
# (submodule of the aggregated navette._navette extension, built from
# rust/navette-py over the navette Rust crate).
#
# Keys are the triple (angle_deg, polarization, spectral) — e.g.
# (0.0, "s", "R"). Wavelengths are in nanometres throughout.
#
# Kept in sync with the PyO3 registration list by tools/check_pyi_sync.py.
from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]

#: (angle_deg, polarization, spectral)
Key = tuple[float, str, str]

class SpectralDataFrame:
    """One wavelength grid plus every curve sampled on it.

    Frames are handed out by `OpticalCollection` / `OpticalWeaver`; there is
    no public constructor.
    """

    @property
    def uid(self) -> int: ...
    @property
    def wavelength(self) -> FloatArray: ...
    @property
    def wl_bounds(self) -> tuple[float, float]: ...
    def keys(self) -> list[Key]: ...
    def set_data(self, key: Key, value: FloatArray, wavelength: FloatArray) -> bool: ...
    def remove(self, key: Key) -> None: ...
    def __getitem__(self, key: Key) -> FloatArray: ...
    def __contains__(self, key: Key) -> bool: ...
    def __len__(self) -> int: ...

class OpticalCollection:
    """Frames grouped by grid fingerprint, so curves on identical grids share
    storage."""

    def __init__(self) -> None: ...
    @property
    def display_spectral(self) -> str: ...
    @display_spectral.setter
    def display_spectral(self, unit_str: str) -> None: ...
    @property
    def display_intensity(self) -> str: ...
    @display_intensity.setter
    def display_intensity(self, unit_str: str) -> None: ...
    @property
    def frame_count(self) -> int: ...
    @property
    def frames(self) -> list[SpectralDataFrame]: ...
    def frame(self, index: int) -> SpectralDataFrame: ...
    def frames_for_key(self, key: Key) -> list[SpectralDataFrame]: ...
    def keys(self) -> list[Key]: ...
    def get_converted(self, key: Key) -> tuple[list[FloatArray], list[FloatArray]]: ...
    def set_data(
        self,
        key: Key,
        value: FloatArray,
        wavelength: FloatArray,
        input_spectral: str | None = ...,
        input_intensity: str | None = ...,
    ) -> None: ...
    def __contains__(self, key: Key) -> bool: ...
    def __len__(self) -> int: ...

class OpticalWeaver:
    """Distributes long curves across frames and reassembles them on demand.

    Cheap to clone and safe to share across threads. Locking is **per frame**:
    concurrent writers to *distinct* keys are independent, but two threads
    writing the same key interleave fragment-wise (last writer wins per
    fragment). See the module docs in `opticalweaver.rs`.
    """

    def __init__(self, cache_size: int = ...) -> None: ...
    @property
    def display_spectral(self) -> str: ...
    @display_spectral.setter
    def display_spectral(self, unit_str: str) -> None: ...
    @property
    def display_intensity(self) -> str: ...
    @display_intensity.setter
    def display_intensity(self, unit_str: str) -> None: ...
    @property
    def frame_count(self) -> int: ...
    @property
    def frames(self) -> list[SpectralDataFrame]: ...
    @property
    def generation(self) -> int: ...
    def frame(self, index: int) -> SpectralDataFrame: ...
    def frames_for_key(self, key: Key) -> list[SpectralDataFrame]: ...
    def keys(self) -> list[Key]: ...
    def get_converted(self, key: Key) -> tuple[list[FloatArray], list[FloatArray]]: ...
    def set_data(
        self,
        key: Key,
        value: FloatArray,
        wavelength: FloatArray,
        input_spectral: str | None = ...,
        input_intensity: str | None = ...,
    ) -> None: ...
    def get_weaved(self, key: Key) -> tuple[FloatArray, FloatArray]:
        """`(wavelength, value)` for one key. Raises `ValueError` if the key
        has never been written."""
        ...

    def get_weaved_collections(self) -> list[tuple[FloatArray, dict[Key, Any]]]: ...
    def unweave(self, key: Key, full_wavelength: FloatArray, full_data: FloatArray) -> int:
        """Distribute one long curve across the frames it spans; returns the
        fragment count. The grid must contain every point of each frame it
        overlaps (exact match to 1e-12, not interpolation)."""
        ...

    def unweave_collection(
        self, common_wavelength: FloatArray, data_batch: Mapping[Key, FloatArray]
    ) -> int:
        """`unweave` for many keys on one shared grid. Rejected batches write
        nothing."""
        ...

    def invalidate_cache(self) -> None: ...
    def __contains__(self, key: Key) -> bool: ...
    def __len__(self) -> int: ...

class TargetWeaver:
    """Target curves in merit space, distributed over the same frame store."""

    def __init__(self, cache_size: int = ..., tolerance_floor: float = ...) -> None: ...
    def add_spectral_target(
        self,
        wavelengths: FloatArray,
        values: FloatArray,
        tolerances: FloatArray,
        angle: float,
        polarization: str,
        spectral: str,
        kind: str,
        norm_mode: str,
        band: FloatArray | None = ...,
        weight: float = ...,
        normalize_count: bool = ...,
        integral: bool = ...,
    ) -> None: ...
    def add_angular_target(
        self,
        wavelength: float,
        angles: FloatArray,
        values: FloatArray,
        tolerances: FloatArray,
        polarization: str,
        spectral: str,
        kind: str,
        norm_mode: str,
        band: FloatArray | None = ...,
        weight: float = ...,
        normalize_count: bool = ...,
        integral: bool = ...,
    ) -> None: ...
    def export_entries(self) -> list[dict[str, Any]]: ...

def calculate_merit(
    sim_weaver: OpticalWeaver, target_weaver: TargetWeaver, missing_penalty: float = ...
) -> float: ...
