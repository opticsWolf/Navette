// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fragment-weaving store: many overlapping spectral curves over shared grids.
//!
//! A [`SpectralDataFrame`] owns one wavelength grid plus the curves sampled on
//! it. An [`OpticalCollection`] groups frames by grid fingerprint so curves on
//! identical grids share storage, and an [`OpticalWeaver`] distributes long
//! input curves across frames (with an LRU-cached distribution plan) and
//! re-assembles continuous curves on demand. All state is `Arc` + interior
//! mutability, so weavers are cheap to clone and share across threads.

use std::borrow::Cow;
use std::hash::{Hash, Hasher};
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use ahash::AHashMap;
use lru::LruCache;
use rayon::prelude::*;
use smallvec::SmallVec;
use parking_lot::RwLock;


// ---------------------------------------------------------------------------
// Unit system
// ---------------------------------------------------------------------------
/// Unit tag for spectral (wavelength) and intensity axes.
/// Conversions are currently identity (grids are stored in nm, data raw);
/// the tag exists so callers can be explicit and future units can convert.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum Unit {
    NM,
    RAW,
}

/// Convert `value` from one unit to another, borrowing when no-op.
/// Returns a borrowed view when `from == to`, otherwise an owned copy.
#[inline]
pub fn convert_unit<'a>(value: &'a [f64], from: Unit, to: Unit) -> Cow<'a, [f64]> {
    if from == to {
        Cow::Borrowed(value)
    } else {
        Cow::Owned(value.to_vec())
    }
}

/// 128-bit fingerprint of a wavelength grid (xxh3 over the raw bits).
/// Used as the LRU/cache key so identical grids share frames and plans.
pub type WlSig = u128;

#[inline]
/// Fingerprint `wl` for frame/plan lookup.
pub fn wl_signature(wl: &[f64]) -> WlSig {
    xxhash_rust::xxh3::xxh3_128(bytemuck::cast_slice::<f64, u8>(wl))
}

#[inline]
/// Bit-exact grid equality (length + raw f64 bits, no tolerance).
/// Guards against xxh3 collisions before frames are shared.
pub fn wl_bits_eq(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && bytemuck::cast_slice::<f64, u8>(a) == bytemuck::cast_slice::<f64, u8>(b)
}

#[inline]
/// Canonical display name of a unit (`"NM"` / `"RAW"`).
pub fn unit_to_str(u: Unit) -> &'static str {
    match u {
        Unit::NM => "NM",
        Unit::RAW => "RAW",
    }
}

#[inline]
/// Parse a spectral-unit label (defaults to nm).
pub fn parse_spectral(s: Option<&str>) -> Unit {
    match s {
        Some("NM") | None => Unit::NM,
        _ => Unit::NM,
    }
}

#[inline]
/// Parse an intensity-unit label (defaults to raw).
pub fn parse_intensity(s: Option<&str>) -> Unit {
    match s {
        Some("RAW") | None => Unit::RAW,
        _ => Unit::RAW,
    }
}

// ---------------------------------------------------------------------------
// Optical key
// ---------------------------------------------------------------------------
/// Identity of one curve: base wavelength plus data-type and polarisation labels.
/// Equality and hashing use the exact f64 bits of `wavelength`.
#[derive(Debug, Clone)]
pub struct OpticalKey {
    pub wavelength: f64,
    pub data_type: Arc<str>,
    pub polarisation: Arc<str>,
}

impl PartialEq for OpticalKey {
    fn eq(&self, other: &Self) -> bool {
        self.wavelength.to_bits() == other.wavelength.to_bits()
            && self.data_type == other.data_type
            && self.polarisation == other.polarisation
    }
}
impl Eq for OpticalKey {}

impl Hash for OpticalKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.wavelength.to_bits().hash(state);
        self.data_type.hash(state);
        self.polarisation.hash(state);
    }
}

impl From<(f64, String, String)> for OpticalKey {
    fn from(t: (f64, String, String)) -> Self {
        OpticalKey {
            wavelength: t.0,
            data_type: t.1.into(),
            polarisation: t.2.into(),
        }
    }
}

impl OpticalKey {
    #[inline]
    /// `(wavelength, data_type, polarisation)` triple for the Python boundary.
    pub fn as_tuple(&self) -> (f64, String, String) {
        (
            self.wavelength,
            self.data_type.to_string(),
            self.polarisation.to_string(),
        )
    }
}

// ---------------------------------------------------------------------------
// SpectralData & SpectralDataFrame
// ---------------------------------------------------------------------------
/// A (possibly strided) view onto a shared data buffer: `buf[start..start+len]`.
/// Contiguous slices can alias the source `Arc` (zero-copy); strided gathers copy.
/// Derefs to `[f64]` so callers treat it as a plain slice.
#[derive(Clone)]
pub struct SpectralData {
    buf: Arc<[f64]>,
    start: usize,
    len: usize,
}

impl SpectralData {
    #[inline]
    /// Wrap a whole owned buffer as one contiguous fragment.
    pub fn from_arc(buf: Arc<[f64]>) -> Self {
        let len = buf.len();
        SpectralData { buf, start: 0, len }
    }
}

impl std::ops::Deref for SpectralData {
    type Target = [f64];
    #[inline]
    fn deref(&self) -> &[f64] {
        &self.buf[self.start..self.start + self.len]
    }
}

/// One wavelength grid plus every curve sampled on it.
/// `uid` is process-unique; `set_data` rejects grids that differ bit-for-bit.
pub struct SpectralDataFrame {
    pub uid: usize,
    data: RwLock<AHashMap<OpticalKey, SpectralData>>,
    wavelength: Arc<[f64]>,
    wl_min: f64,
    wl_max: f64,
}

impl SpectralDataFrame {
    /// Create a frame on a strictly-increasing grid. Errors on empty or
    /// non-monotonic grids.
    pub fn new(wavelength: &[f64]) -> Result<Self, String> {
        static UID_GEN: AtomicUsize = AtomicUsize::new(0);
        let uid = UID_GEN.fetch_add(1, Ordering::SeqCst);

        if wavelength.is_empty() {
            return Err("SpectralDataFrame: wavelength array must be non-empty.".to_string());
        }
        for i in 0..wavelength.len() - 1 {
            if !(wavelength[i] < wavelength[i + 1]) {
                return Err("SpectralDataFrame: wavelength array must be strictly monotonically increasing.".to_string());
            }
        }

        let wl_min = wavelength[0];
        let wl_max = wavelength[wavelength.len() - 1];

        Ok(SpectralDataFrame {
            uid,
            data: RwLock::new(AHashMap::new()),
            wavelength: Arc::from(wavelength),
            wl_min,
            wl_max,
        })
    }

    /// Store `value` under `key`. Returns whether the key is new.
    /// Errors when `value` misaligns with the grid (or the optional grid
    /// conflicts with it).
    pub fn set_data(
        &self,
        key: OpticalKey,
        value: SpectralData,
        wavelength: Option<&[f64]>,
    ) -> Result<bool, String> {
        if let Some(wl) = wavelength
            && !wl_bits_eq(&self.wavelength, wl) {
            return Err(format!(
                "SpectralDataFrame(uid={}): wavelength grid conflict.",
                self.uid
            ));
        }
        if value.len() != self.wavelength.len() {
            return Err(format!(
                "SpectralDataFrame(uid={}): {} value(s) for a {}-point grid.",
                self.uid,
                value.len(),
                self.wavelength.len()
            ));
        }

        // `insert` already tells us whether the key was there, so this is one
        // hash of `key` rather than the `contains_key` + `insert` pair's two.
        // A batch unweave does this once per (key, frame) -- tens of thousands
        // of times for a 512-key call.
        Ok(self.data.write().insert(key, value).is_none())
    }

    /// Cloned view of the curve under `key`, if present.
    pub fn get_data(&self, key: &OpticalKey) -> Option<SpectralData> {
        self.data.read().get(key).cloned()
    }

    /// All keys stored in this frame.
    pub fn keys(&self) -> Vec<OpticalKey> {
        self.data.read().keys().cloned().collect()
    }

    /// Number of curves stored in this frame.
    pub fn len(&self) -> usize {
        self.data.read().len()
    }

    /// Whether the frame holds no curves. A frame can legitimately be empty:
    /// it is created for a wavelength grid and populated afterwards.
    pub fn is_empty(&self) -> bool {
        self.data.read().is_empty()
    }

    /// Borrowed grid shared by all curves in this frame.
    pub fn wavelength(&self) -> &[f64] {
        &self.wavelength
    }

    /// `(min, max)` grid endpoints (cached at construction).
    pub fn wl_bounds(&self) -> (f64, f64) {
        (self.wl_min, self.wl_max)
    }

    /// Drop the curve under `key`. Errors when the key is absent.
    pub fn remove(&self, key: &OpticalKey) -> Result<(), String> {
        if self.data.write().remove(key).is_none() {
            return Err("Key not found".to_string());
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// OpticalCollection
// ---------------------------------------------------------------------------
/// A set of frames indexed by grid fingerprint plus a key→frames map.
/// Curves on identical grids share one frame (fingerprints + bit-equality
/// guard); `key_map` tracks which frames hold fragments of each key.
pub struct OpticalCollection {
    pub frames: RwLock<Vec<Arc<SpectralDataFrame>>>,
    wl_fingerprints: RwLock<AHashMap<WlSig, Arc<SpectralDataFrame>>>,
    key_map: RwLock<AHashMap<OpticalKey, SmallVec<[Arc<SpectralDataFrame>; 2]>>>,
    display_spectral: RwLock<Unit>,
    display_intensity: RwLock<Unit>,
}

impl Default for OpticalCollection {
    fn default() -> Self {
        Self::new()
    }
}

impl OpticalCollection {
    /// Empty collection with nm/raw display units.
    pub fn new() -> Self {
        OpticalCollection {
            frames: RwLock::new(Vec::new()),
            wl_fingerprints: RwLock::new(AHashMap::new()),
            key_map: RwLock::new(AHashMap::new()),
            display_spectral: RwLock::new(Unit::NM),
            display_intensity: RwLock::new(Unit::RAW),
        }
    }

    /// Display unit for wavelength axes returned to callers.
    pub fn set_display_spectral(&self, unit: Unit) {
        *self.display_spectral.write() = unit;
    }
    /// Current wavelength display unit.
    pub fn display_spectral(&self) -> Unit {
        *self.display_spectral.read()
    }
    /// Display unit for data values returned to callers.
    pub fn set_display_intensity(&self, unit: Unit) {
        *self.display_intensity.write() = unit;
    }
    /// Current data display unit.
    pub fn display_intensity(&self) -> Unit {
        *self.display_intensity.read()
    }
    /// Number of frames (distinct grids) held.
    pub fn frame_count(&self) -> usize {
        self.frames.read().len()
    }
    /// Number of distinct keys across all frames.
    pub fn len_keys(&self) -> usize {
        self.key_map.read().len()
    }
    /// All distinct keys across all frames.
    pub fn keys(&self) -> Vec<OpticalKey> {
        self.key_map.read().keys().cloned().collect()
    }
    /// True when any frame holds fragments of `key`.
    pub fn contains_key(&self, key: &OpticalKey) -> bool {
        self.key_map.read().contains_key(key)
    }
    /// Frame by insertion index, if in range.
    pub fn frame_at(&self, index: usize) -> Option<Arc<SpectralDataFrame>> {
        self.frames.read().get(index).cloned()
    }
    /// Cloned list of all frames (lock-free iteration for callers).
    pub fn frames_snapshot(&self) -> Vec<Arc<SpectralDataFrame>> {
        self.frames.read().clone()
    }
    /// Frames holding fragments of `key` (inline for the common ≤ 2 case).
    pub fn frames_for_key(&self, key: &OpticalKey) -> Option<SmallVec<[Arc<SpectralDataFrame>; 2]>> {
        self.key_map.read().get(key).cloned()
    }

    /// Per-fragment `(data, wavelength)` pairs under `key`, converted to the
    /// display units. Returns `None` when the key is unknown.
    pub fn get_converted(&self, key: &OpticalKey) -> Option<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        let frames = self.key_map.read().get(key)?.clone();
        let int_unit = *self.display_intensity.read();
        let spec_unit = *self.display_spectral.read();

        let mut data_list = Vec::with_capacity(frames.len());
        let mut wl_list = Vec::with_capacity(frames.len());
        for frm in frames {
            let data_raw = frm
                .get_data(key)
                .unwrap_or_else(|| SpectralData::from_arc(Arc::from(&[][..])));
            data_list.push(convert_unit(&data_raw, Unit::RAW, int_unit).into_owned());
            wl_list.push(convert_unit(frm.wavelength(), Unit::NM, spec_unit).into_owned());
        }
        Some((data_list, wl_list))
    }

    /// Record that `frame` holds fragments of `key`. Returns false when the
    /// mapping already existed.
    pub fn map_frame_to_key(&self, key: &OpticalKey, frame: &Arc<SpectralDataFrame>) -> bool {
        self.map_frames_to_key(key, std::slice::from_ref(frame)) == 1
    }

    /// Record that every frame in `frames` holds fragments of `key`, under a
    /// single lock. Returns how many were not already recorded.
    ///
    /// The singular form above takes the collection-wide write lock once per
    /// *fragment*; a batch unweave writes one fragment per (key, frame) pair,
    /// which is tens of thousands for a 512-key call and all of it contending
    /// on one lock. This form takes it once per key. Frames keep the order
    /// they arrive in -- the distribution plan's order -- which is what
    /// `get_converted` hands back, so it must not depend on how the batch was
    /// scheduled; each key is still processed by one thread, start to finish.
    pub fn map_frames_to_key(&self, key: &OpticalKey, frames: &[Arc<SpectralDataFrame>]) -> usize {
        let mut key_map = self.key_map.write();
        let mapped = key_map.entry(key.clone()).or_default();
        let mut added = 0;
        for frame in frames {
            if mapped.iter().any(|f| Arc::ptr_eq(f, frame)) {
                continue;
            }
            mapped.push(Arc::clone(frame));
            added += 1;
        }
        added
    }

    pub fn set_data(
        &self,
        key: OpticalKey,
        value: &[f64],
        wavelength: &[f64],
        input_spectral: Unit,
        input_intensity: Unit,
    ) -> Result<bool, String> {
        let base_wl = convert_unit(wavelength, input_spectral, Unit::NM);
        let base_data = convert_unit(value, input_intensity, Unit::RAW);
        let (target_frame, frame_created) = self.get_or_create_frame(&base_wl)?;
        let data_arc: Arc<[f64]> = match base_data {
            Cow::Borrowed(b) => Arc::from(b),
            Cow::Owned(o) => Arc::from(o),
        };
        let is_new = target_frame.set_data(key.clone(), SpectralData::from_arc(data_arc), None)?;
        if is_new {
            self.map_frame_to_key(&key, &target_frame);
        }
        Ok(frame_created)
    }

    /// Fetch the frame for `wl_arr` (bit-exact grid match) or create it.
    /// Returns the frame plus whether it is newly created.
    fn get_or_create_frame(&self, wl_arr: &[f64]) -> Result<(Arc<SpectralDataFrame>, bool), String> {
        let sig = wl_signature(wl_arr);
        if let Some(frm) = self.wl_fingerprints.read().get(&sig)
            && wl_bits_eq(frm.wavelength(), wl_arr) {
            return Ok((frm.clone(), false));
        }
        let mut fp = self.wl_fingerprints.write();
        if let Some(frm) = fp.get(&sig)
            && wl_bits_eq(frm.wavelength(), wl_arr) {
            return Ok((frm.clone(), false));
        }
        let new_frame = Arc::new(SpectralDataFrame::new(wl_arr)?);
        self.frames.write().push(new_frame.clone());
        fp.insert(sig, new_frame.clone());
        Ok((new_frame, true))
    }
}

// ---------------------------------------------------------------------------
// Distribution plan & OpticalWeaver
// ---------------------------------------------------------------------------
#[derive(Clone)]
enum SliceOrIndices {
    Slice(usize, usize),
    Indices(Vec<usize>),
}
impl SliceOrIndices {
    /// How many source points this entry contributes, without gathering them.
    #[inline]
    fn len(&self) -> usize {
        match self {
            SliceOrIndices::Slice(s, e) => e - s,
            SliceOrIndices::Indices(idx) => idx.len(),
        }
    }

    /// Extract this entry's fragment from `data`. When `shared` is `Some` it must
    /// be a copy of `data` over a span containing this entry; a contiguous slice
    /// then becomes a zero-copy view into that buffer. Strided fragments, and the
    /// `shared == None` case, copy out.
    #[inline]
    fn gather(&self, data: &[f64], shared: Option<&SharedSource>) -> SpectralData {
        match self {
            SliceOrIndices::Slice(s, e) => match shared {
                Some(src) => SpectralData {
                    buf: Arc::clone(&src.buf),
                    start: s - src.base,
                    len: e - s,
                },
                None => SpectralData::from_arc(Arc::from(&data[*s..*e])),
            },
            SliceOrIndices::Indices(idx) => {
                let vec: Vec<f64> = idx.iter().map(|&i| data[i]).collect();
                SpectralData::from_arc(Arc::from(vec))
            }
        }
    }
}

/// One key's source curve, copied over just the span the plan reaches, so every
/// contiguous fragment of that key can be a view into it instead of a copy of
/// its own.
///
/// The engine has to own this: the caller's buffer is borrowed for the length of
/// the call (and, from Python, is a numpy array the caller may write to
/// afterwards), while a fragment lives as long as the frame holding it. `base`
/// is the source offset the copy starts at -- frames rarely tile the whole input
/// curve, and copying the part of it no fragment reads is the one saving
/// available here that does not change what a fragment *is*.
struct SharedSource {
    buf: Arc<[f64]>,
    base: usize,
}

type DistributionPlan = Vec<(Arc<SpectralDataFrame>, SliceOrIndices)>;

/// Above this many fragments in one batch unweave -- frames in the plan times
/// keys in the batch -- the per-key loop goes to the rayon pool.
///
/// The count, not the byte volume, is the thing to gate on. The bytes are a
/// `memcpy` per key and that is bound by DRAM bandwidth, which does not
/// parallelise: 512 keys x 100k points measures 39.7 / 35.1 / 36.5 ms on 1 / 4
/// / 32 threads, 22 GB/s of traffic either way. What *does* scale is the
/// per-fragment work -- gather, hash, frame lock -- and that goes as the
/// fragment count. Measured across the bench grid the split is clean: every
/// configuration with 16384 fragments gains (1.06x to 2.11x) and every one with
/// 4096 or fewer loses (up to 2.2x at 10k points x 8 frames x 32 keys, where
/// the pool hand-off costs more than the whole call). 8192 sits between them.
///
/// Serial and parallel write the same bytes to the same places; this is a speed
/// knob, not a semantic one.
const UNWEAVE_PAR_FRAGMENTS: usize = 8192;

/// Top-level weaving store: fragment distribution plus re-assembly.
/// `generation` bumps on every structural change so caches and Python-side
/// snapshots can detect staleness; `invalidate_cache` clears the LRU plan cache.
pub struct OpticalWeaver {
    pub inner: OpticalCollection,
    distribution_cache: RwLock<LruCache<WlSig, (usize, Arc<[f64]>, DistributionPlan)>>,
    generation: AtomicUsize,
}

impl OpticalWeaver {
    /// Create a weaver with an LRU distribution-plan cache of `cache_size` grids.
    pub fn new(cache_size: usize) -> Self {
        OpticalWeaver {
            inner: OpticalCollection::new(),
            distribution_cache: RwLock::new(LruCache::new(
                NonZeroUsize::new(cache_size.max(1)).unwrap(),
            )),
            generation: AtomicUsize::new(0),
        }
    }

    /// Record a structural change (new frame or new key mapping).
    pub fn bump_generation(&self) {
        self.generation.fetch_add(1, Ordering::SeqCst);
    }
    /// Monotonic structural-change counter for staleness checks.
    pub fn generation(&self) -> usize {
        self.generation.load(Ordering::SeqCst)
    }

    /// Ingest one curve: convert units, route to the matching frame
    /// (creating it when the grid is new), and bump the generation on
    /// structural change.
    pub fn set_data(
        &self,
        key: OpticalKey,
        value: &[f64],
        wavelength: &[f64],
        input_spectral: Unit,
        input_intensity: Unit,
    ) -> Result<(), String> {
        if self
            .inner
            .set_data(key, value, wavelength, input_spectral, input_intensity)?
        {
            self.bump_generation();
        }
        Ok(())
    }

    /// Re-assemble the continuous `(wavelength, data)` curve for `key` by
    /// concatenating its fragments in grid order. Errors on unknown keys.
    pub fn get_weaved(&self, key: &OpticalKey) -> Result<(Vec<f64>, Vec<f64>), String> {
        let frames = self
            .inner
            .frames_for_key(key)
            .ok_or_else(|| "Key not found.".to_string())?;
        let mut fragments: SmallVec<[(f64, &[f64], SpectralData); 4]> = SmallVec::new();
        for frm in &frames {
            if let Some(data) = frm.get_data(key) {
                fragments.push((frm.wl_bounds().0, frm.wavelength(), data));
            }
        }
        if fragments.is_empty() {
            return Ok((Vec::new(), Vec::new()));
        }
        fragments.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

        let total: usize = fragments.iter().map(|f| f.1.len()).sum();
        let mut all_wl = Vec::with_capacity(total);
        let mut all_data = Vec::with_capacity(total);
        for (_, wl, data) in fragments {
            all_wl.extend_from_slice(wl);
            all_data.extend_from_slice(&data);
        }
        Ok((all_wl, all_data))
    }

    /// All keys re-assembled, grouped by woven grid: one entry per distinct
    /// grid with the curves sampled on it.
    pub fn get_weaved_collections(&self) -> Vec<(Vec<f64>, AHashMap<OpticalKey, Vec<f64>>)> {
        let mut groups: AHashMap<WlSig, (Vec<f64>, AHashMap<OpticalKey, Vec<f64>>)> = AHashMap::new();
        for key in self.inner.keys() {
            if let Ok((wl, data)) = self.get_weaved(&key) {
                if wl.is_empty() {
                    continue;
                }
                let sig = wl_signature(&wl);
                let entry = groups.entry(sig).or_insert_with(|| (wl.clone(), AHashMap::new()));
                entry.1.insert(key, data);
            }
        }
        groups.into_values().collect()
    }

    /// Distribute one long curve across frames per the cached plan.
    /// Returns the number of fragments written.
    /// Explain a short gather before it reaches the frame as a bare count
    /// mismatch. `unweave` distributes by **exact** wavelength match, so a
    /// target curve that merely spans the same range -- a denser or coarser
    /// linspace over the same interval, say -- lands on a handful of shared
    /// points and nothing says why. Names what was covered and what the
    /// contract is.
    fn check_coverage(frm: &SpectralDataFrame, got: usize) -> Result<(), String> {
        let want = frm.wavelength().len();
        if got == want {
            return Ok(());
        }
        let (lo, hi) = frm.wl_bounds();
        Err(format!(
            "unweave: the supplied grid carries {got} of the {want} wavelengths in the frame spanning [{lo}, {hi}] nm. A target curve is distributed by exact wavelength match (1e-12), not by interpolation, so it must contain every point of each frame it overlaps -- take the grid from get_weaved() rather than building a fresh linspace over the same range."
        ))
    }

    /// Everything `unweave` validates that depends only on the plan and the
    /// grid, not on which curve is being distributed: each fragment has to
    /// cover its frame's grid exactly.
    ///
    /// Hoisting it out of the per-key loop is what makes a batch all-or-nothing.
    /// Checked inside the loop it would reject the batch after an arbitrary
    /// prefix of it had already been written -- and once the loop runs in
    /// parallel, "arbitrary" would mean *scheduling-dependent*.
    fn check_plan(plan: &DistributionPlan) -> Result<(), String> {
        for (frm, indices) in plan {
            Self::check_coverage(frm, indices.len())?;
        }
        Ok(())
    }

    /// Reject a curve that is not one value per wavelength, before the plan
    /// indexes it.
    ///
    /// The plan addresses the source by position, so a short curve would index
    /// out of bounds -- a panic, and from a rayon worker at that. The Python
    /// wrapper checks this for a single `unweave` but not for a batch, which is
    /// exactly the path that reaches here with many curves at once.
    fn check_length(key: &OpticalKey, got: usize, want: usize) -> Result<(), String> {
        if got == want {
            return Ok(());
        }
        Err(format!(
            "unweave: the curve for key ({}, {}, {}) has {got} value(s) for a {want}-point grid; a target curve must carry one value per wavelength.",
            key.wavelength, key.data_type, key.polarisation
        ))
    }

    /// The half-open source span every contiguous fragment of `plan` lies in,
    /// or `None` when the plan has none (an all-strided plan gathers
    /// element-wise and never reads a shared buffer).
    fn shared_span(plan: &DistributionPlan) -> Option<(usize, usize)> {
        let mut lo = usize::MAX;
        let mut hi = 0usize;
        for (_, indices) in plan {
            if let SliceOrIndices::Slice(s, e) = indices {
                lo = lo.min(*s);
                hi = hi.max(*e);
            }
        }
        (lo < hi).then_some((lo, hi))
    }

    /// Distribute one curve over `plan` and register the frames it is new to.
    /// Returns the fragments written.
    ///
    /// Infallible in practice once `check_plan` and `check_length` have passed
    /// -- the only error `set_data` can still return is the length mismatch
    /// `check_plan` already ruled out -- but it is not the place to assert that,
    /// so the `Result` stays.
    fn unweave_one(
        &self,
        plan: &DistributionPlan,
        span: Option<(usize, usize)>,
        key: &OpticalKey,
        full_data: &[f64],
    ) -> Result<usize, String> {
        let shared = span.map(|(lo, hi)| SharedSource {
            buf: Arc::from(&full_data[lo..hi]),
            base: lo,
        });
        let mut fresh: SmallVec<[Arc<SpectralDataFrame>; 2]> = SmallVec::new();
        for (frm, indices) in plan {
            let subset = indices.gather(full_data, shared.as_ref());
            if frm.set_data(key.clone(), subset, None)? {
                fresh.push(Arc::clone(frm));
            }
        }
        if !fresh.is_empty() {
            self.inner.map_frames_to_key(key, &fresh);
        }
        Ok(plan.len())
    }

    pub fn unweave(
        &self,
        key: OpticalKey,
        full_wavelength: &[f64],
        full_data: &[f64],
    ) -> Result<usize, String> {
        let plan = self.resolve_plan(full_wavelength)?;
        Self::check_plan(&plan)?;
        Self::check_length(&key, full_data.len(), full_wavelength.len())?;
        self.unweave_one(&plan, Self::shared_span(&plan), &key, full_data)
    }

    /// Distribute many curves sharing one grid, reusing a single plan.
    /// Returns the total fragments written.
    ///
    /// The per-key work is a copy of the source span plus one map insert per
    /// frame; it shares nothing across keys but the plan, which is read-only
    /// here, so above [`UNWEAVE_PAR_FRAGMENTS`] fragments it runs on the rayon
    /// pool. What each key writes is independent of the others -- distinct
    /// map entries, in frames that already exist -- and each key is handled by
    /// one thread from start to finish, so the frame order recorded under a key
    /// is the plan's either way. The keys arrive in an `AHashMap`'s iteration
    /// order, which is already not the caller's, so nothing observable was
    /// ordered by the loop to begin with.
    pub fn unweave_collection(
        &self,
        common_wavelength: &[f64],
        data_batch: AHashMap<OpticalKey, &[f64]>,
    ) -> Result<usize, String> {
        if data_batch.is_empty() {
            return Ok(0);
        }
        let plan = self.resolve_plan(common_wavelength)?;
        Self::check_plan(&plan)?;
        let items: Vec<(&OpticalKey, &[f64])> =
            data_batch.iter().map(|(k, d)| (k, *d)).collect();
        for (key, full_data) in &items {
            Self::check_length(key, full_data.len(), common_wavelength.len())?;
        }

        let span = Self::shared_span(&plan);
        let fragments = plan.len().saturating_mul(items.len());
        let written: Vec<Result<usize, String>> =
            if items.len() > 1 && fragments >= UNWEAVE_PAR_FRAGMENTS {
                items
                    .par_iter()
                    .map(|(key, full_data)| self.unweave_one(&plan, span, key, full_data))
                    .collect()
            } else {
                items
                    .iter()
                    .map(|(key, full_data)| self.unweave_one(&plan, span, key, full_data))
                    .collect()
            };

        let mut total = 0usize;
        for count in written {
            total += count?;
        }
        Ok(total)
    }

    /// Drop all cached distribution plans (grids re-resolve on next use).
    pub fn invalidate_cache(&self) {
        self.distribution_cache.write().clear();
    }

    /// LRU-cached distribution plan for `full_wavelength` (hit or rebuild).
    fn resolve_plan(&self, full_wavelength: &[f64]) -> Result<DistributionPlan, String> {
        let sig = wl_signature(full_wavelength);
        let current_gen = self.generation.load(Ordering::SeqCst);
        {
            let mut cache = self.distribution_cache.write();
            if let Some((cached_gen, cached_wl, plan)) = cache.get(&sig)
                && *cached_gen == current_gen && wl_bits_eq(cached_wl, full_wavelength) {
                return Ok(plan.clone());
            }
            cache.pop(&sig);
        }
        let plan = self.build_distribution_plan(full_wavelength)?;
        self.distribution_cache.write().put(
            sig,
            (current_gen, Arc::from(full_wavelength), plan.clone()),
        );
        Ok(plan)
    }

    /// Fragment `full_wavelength` across existing frames (contiguous slices
    /// where grids overlap, strided gathers elsewhere).
    fn build_distribution_plan(&self, full_wavelength: &[f64]) -> Result<DistributionPlan, String> {
        let mut plan = Vec::new();
        if full_wavelength.is_empty() {
            return Ok(plan);
        }
        let (fw_min, fw_max) = (
            full_wavelength[0],
            full_wavelength[full_wavelength.len() - 1],
        );

        for frm in self.inner.frames_snapshot() {
            let frame_wl = frm.wavelength();
            let (f_min, f_max) = frm.wl_bounds();
            if f_min > fw_max || f_max < fw_min {
                continue;
            }

            let idx_lo = full_wavelength.partition_point(|&x| x < f_min);
            let idx_hi = full_wavelength.partition_point(|&x| x <= f_max);
            let candidate = &full_wavelength[idx_lo..idx_hi];

            if candidate.len() == frame_wl.len() && candidate == frame_wl {
                plan.push((frm, SliceOrIndices::Slice(idx_lo, idx_hi)));
                continue;
            }

            // Two-pointer sweep of two sorted grids (O(n+m)): `fw` advances
            // monotonically as we walk the (strictly increasing) frame grid,
            // collecting the exact-match positions. Same sweep shape as the merit
            // interpolation loop, minus the value math (this aligns indices only).
            let mut indices = Vec::with_capacity(frame_wl.len());
            let mut fw = idx_lo;
            for &target in frame_wl {
                while fw < full_wavelength.len() && full_wavelength[fw] < target {
                    fw += 1;
                }
                if fw < full_wavelength.len() && (full_wavelength[fw] - target).abs() < 1e-12 {
                    indices.push(fw);
                    fw += 1; // exact match consumed; next target is strictly greater
                }
            }
            if !indices.is_empty() {
                if indices
                    .iter()
                    .enumerate()
                    .all(|(i, &idx)| i == 0 || idx == indices[i - 1] + 1)
                {
                    plan.push((
                        frm,
                        SliceOrIndices::Slice(indices[0], indices[indices.len() - 1] + 1),
                    ));
                } else {
                    plan.push((frm, SliceOrIndices::Indices(indices)));
                }
            }
        }
        Ok(plan)
    }
}


// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    fn key(wl: f64, pol: &str) -> OpticalKey {
        OpticalKey::from((wl, "R".to_string(), pol.to_string()))
    }

    /// A weaver whose frames tile `[0, n)` in `n / frames` blocks, seeded so the
    /// frames exist before anything is unwoven onto them.
    fn tiled(n: usize, frames: usize) -> (OpticalWeaver, Vec<f64>) {
        let grid: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let w = OpticalWeaver::new(16);
        let block = n / frames;
        let seed = key(0.0, "seed");
        for f in 0..frames {
            let sub = &grid[f * block..(f + 1) * block];
            w.set_data(seed.clone(), &vec![0.0; sub.len()], sub, Unit::NM, Unit::RAW)
                .unwrap();
        }
        (w, grid)
    }

    fn curve(grid: &[f64], scale: f64) -> Vec<f64> {
        grid.iter().map(|&x| x * scale + 1.0).collect()
    }

    // -- the length guard --------------------------------------------------

    #[test]
    fn a_short_curve_in_a_batch_is_an_error_not_a_panic() {
        let (w, grid) = tiled(64, 4);
        let k = key(500.0, "s");
        let short = vec![0.0; 32];
        let mut batch = AHashMap::new();
        batch.insert(k, short.as_slice());
        let err = w.unweave_collection(&grid, batch).unwrap_err();
        assert!(err.contains("32 value(s) for a 64-point grid"), "{err}");
        assert!(err.contains("one value per wavelength"), "{err}");
    }

    #[test]
    fn the_length_error_names_the_key_it_came_from() {
        let (w, grid) = tiled(64, 4);
        let mut batch = AHashMap::new();
        let short = vec![0.0; 10];
        batch.insert(key(633.0, "p"), short.as_slice());
        let err = w.unweave_collection(&grid, batch).unwrap_err();
        assert!(err.contains("(633, R, p)"), "{err}");
    }

    #[test]
    fn a_long_curve_is_rejected_too() {
        // Indexing would succeed; the curve is still not one value per
        // wavelength, and silently ignoring the tail is how a caller's off-by-one
        // grid goes unnoticed.
        let (w, grid) = tiled(64, 4);
        let long = vec![0.0; 65];
        assert!(w.unweave(key(500.0, "s"), &grid, &long).is_err());
    }

    // -- the batch is all-or-nothing ---------------------------------------

    #[test]
    fn a_grid_that_misses_part_of_a_frame_writes_nothing() {
        // The plan is checked before any key is distributed, so a bad grid
        // cannot leave an arbitrary prefix of the batch written.
        let (w, grid) = tiled(64, 4);
        let every_other: Vec<f64> = grid.iter().step_by(2).copied().collect();
        let data = curve(&every_other, 1.0);
        let mut batch = AHashMap::new();
        for pol in ["s", "p"] {
            batch.insert(key(500.0, pol), data.as_slice());
        }
        assert!(w.unweave_collection(&every_other, batch).is_err());
        assert!(!w.inner.contains_key(&key(500.0, "s")));
        assert!(!w.inner.contains_key(&key(500.0, "p")));
    }

    // -- batch == per-key ---------------------------------------------------

    #[test]
    fn a_batch_writes_exactly_what_the_per_key_calls_would() {
        let (batched, grid) = tiled(600, 6);
        let (singly, _) = tiled(600, 6);
        let curves: Vec<Vec<f64>> = (0..8).map(|i| curve(&grid, i as f64 + 0.5)).collect();
        let keys: Vec<OpticalKey> = (0..8).map(|i| key(400.0 + i as f64, "s")).collect();

        let mut batch = AHashMap::new();
        for (k, c) in keys.iter().zip(&curves) {
            batch.insert(k.clone(), c.as_slice());
        }
        assert_eq!(batched.unweave_collection(&grid, batch).unwrap(), 8 * 6);
        for (k, c) in keys.iter().zip(&curves) {
            singly.unweave(k.clone(), &grid, c).unwrap();
        }
        for (k, c) in keys.iter().zip(&curves) {
            let (bw, bd) = batched.get_weaved(k).unwrap();
            let (sw, sd) = singly.get_weaved(k).unwrap();
            assert_eq!(bd, sd, "data differs for {}", k.wavelength);
            assert_eq!(bw, sw);
            assert_eq!(bd, *c);
        }
    }

    #[test]
    fn the_parallel_path_writes_what_the_serial_path_writes() {
        // 64 frames x 256 keys = 16384 fragments, comfortably over
        // UNWEAVE_PAR_FRAGMENTS; the 4-frame twin stays serial. Same answers.
        let n = 1024;
        let (par, grid) = tiled(n, 64);
        let (ser, _) = tiled(n, 4);
        let curves: Vec<Vec<f64>> = (0..256).map(|i| curve(&grid, i as f64)).collect();
        let keys: Vec<OpticalKey> = (0..256).map(|i| key(i as f64, "s")).collect();

        let mut a = AHashMap::new();
        let mut b = AHashMap::new();
        for (k, c) in keys.iter().zip(&curves) {
            a.insert(k.clone(), c.as_slice());
            b.insert(k.clone(), c.as_slice());
        }
        const { assert!(64 * 256 >= UNWEAVE_PAR_FRAGMENTS) };
        const { assert!(4 * 256 < UNWEAVE_PAR_FRAGMENTS) };
        assert_eq!(par.unweave_collection(&grid, a).unwrap(), 256 * 64);
        assert_eq!(ser.unweave_collection(&grid, b).unwrap(), 256 * 4);
        for (k, c) in keys.iter().zip(&curves) {
            assert_eq!(par.get_weaved(k).unwrap().1, *c);
            assert_eq!(ser.get_weaved(k).unwrap().1, *c);
        }
    }

    #[test]
    fn every_key_records_its_frames_in_plan_order() {
        // `get_converted` hands fragments back in the order they were recorded,
        // so that order must not depend on how the batch was scheduled.
        let n = 1024;
        let (w, grid) = tiled(n, 64);
        let curves: Vec<Vec<f64>> = (0..256).map(|i| curve(&grid, i as f64)).collect();
        let mut batch = AHashMap::new();
        let keys: Vec<OpticalKey> = (0..256).map(|i| key(i as f64, "s")).collect();
        for (k, c) in keys.iter().zip(&curves) {
            batch.insert(k.clone(), c.as_slice());
        }
        w.unweave_collection(&grid, batch).unwrap();
        let reference: Vec<usize> = w
            .inner
            .frames_for_key(&keys[0])
            .unwrap()
            .iter()
            .map(|f| f.uid)
            .collect();
        assert_eq!(reference.len(), 64);
        for k in &keys[1..] {
            let uids: Vec<usize> = w
                .inner
                .frames_for_key(k)
                .unwrap()
                .iter()
                .map(|f| f.uid)
                .collect();
            assert_eq!(uids, reference);
        }
    }

    // -- the span-limited shared source -------------------------------------

    #[test]
    fn a_plan_reaching_only_the_tail_of_the_curve_still_gathers_the_tail() {
        // Frames cover [300, 600) of a 600-point grid, so the shared copy starts
        // at 300 and every fragment index has to be rebased onto it. Getting the
        // offset wrong reads the head of the curve and nothing else notices.
        let grid: Vec<f64> = (0..600).map(|i| i as f64).collect();
        let w = OpticalWeaver::new(8);
        let seed = key(0.0, "seed");
        for (lo, hi) in [(300, 450), (450, 600)] {
            w.set_data(seed.clone(), &vec![0.0; hi - lo], &grid[lo..hi], Unit::NM, Unit::RAW)
                .unwrap();
        }
        let c = curve(&grid, 3.0);
        let k = key(500.0, "s");
        let mut batch = AHashMap::new();
        batch.insert(k.clone(), c.as_slice());
        assert_eq!(w.unweave_collection(&grid, batch).unwrap(), 2);
        let (wl, data) = w.get_weaved(&k).unwrap();
        assert_eq!(wl, grid[300..]);
        assert_eq!(data, c[300..]);
    }

    #[test]
    fn shared_span_is_none_when_no_fragment_is_contiguous() {
        // An all-strided plan gathers element-wise and must not materialise a
        // copy it will never read.
        let grid: Vec<f64> = (0..200).map(|i| i as f64).collect();
        let strided: Vec<f64> = grid.iter().step_by(2).copied().collect();
        let w = OpticalWeaver::new(8);
        w.set_data(key(0.0, "seed"), &vec![0.0; strided.len()], &strided, Unit::NM, Unit::RAW)
            .unwrap();
        let c = curve(&grid, 2.0);
        let k = key(500.0, "s");
        assert_eq!(w.unweave(k.clone(), &grid, &c).unwrap(), 1);
        let (_, data) = w.get_weaved(&k).unwrap();
        assert_eq!(data, c.iter().step_by(2).copied().collect::<Vec<f64>>());
    }

    // -- the key map --------------------------------------------------------

    #[test]
    fn mapping_the_same_frames_twice_adds_nothing() {
        let (w, grid) = tiled(64, 4);
        let c = curve(&grid, 1.0);
        let k = key(500.0, "s");
        for _ in 0..3 {
            let mut batch = AHashMap::new();
            batch.insert(k.clone(), c.as_slice());
            assert_eq!(w.unweave_collection(&grid, batch).unwrap(), 4);
        }
        assert_eq!(w.inner.frames_for_key(&k).unwrap().len(), 4);
        assert_eq!(w.inner.len_keys(), 2); // the seed key and this one
    }

    #[test]
    fn map_frames_to_key_reports_only_what_was_new() {
        let (w, _) = tiled(64, 4);
        let frames = w.inner.frames_snapshot();
        let k = key(500.0, "s");
        assert_eq!(w.inner.map_frames_to_key(&k, &frames), 4);
        assert_eq!(w.inner.map_frames_to_key(&k, &frames), 0);
        assert!(!w.inner.map_frame_to_key(&k, &frames[0]));
    }

    #[test]
    fn an_empty_batch_is_a_no_op() {
        let (w, grid) = tiled(64, 4);
        assert_eq!(w.unweave_collection(&grid, AHashMap::new()).unwrap(), 0);
        assert_eq!(w.inner.len_keys(), 1);
    }

    #[test]
    fn overwriting_a_key_does_not_re_count_it_as_new() {
        let (w, grid) = tiled(64, 4);
        let a = curve(&grid, 1.0);
        let b = curve(&grid, 7.0);
        let k = key(500.0, "s");
        for c in [&a, &b] {
            let mut batch = AHashMap::new();
            batch.insert(k.clone(), c.as_slice());
            w.unweave_collection(&grid, batch).unwrap();
        }
        assert_eq!(w.get_weaved(&k).unwrap().1, b);
    }
}
