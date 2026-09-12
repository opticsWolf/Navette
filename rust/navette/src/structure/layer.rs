// SPDX-License-Identifier: LGPL-3.0-or-later
//! One design layer: material name + geometry + flags.
//!
//! Mirrors `navette.structure.models.Layer` exactly (field names, defaults,
//! refinement rule, state keys). Differences are deliberate and tested:
//! - `sub_layer_count` is derived ([`Layer::sub_layer_count`]), never stored.
//! - `set_properties` returns its warnings instead of emitting them (the
//!   Python boundary re-emits via `warnings.warn`).
//! - Enum fields are typed; raw-int coercion is fail-closed
//!   (`try_from_i32`, never a default).
//!
//! All lengths in nanometres.

use std::collections::BTreeMap;
use std::fmt;

use serde::ser::SerializeMap;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_json::Value;

use crate::structure::enums::{LayerType, RoughnessType};
use crate::structure::validation::ValidationIssue;
use crate::structure::version::{SCHEMA_VERSION, check_schema_version};

/// Serialize a fieldless enum as its wire int (solver contract, not name).
pub(crate) fn ser_int<T: Copy + Into<i32>, S: Serializer>(v: &T, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_i32((*v).into())
}

macro_rules! impl_serde_int {
    ($t:ty) => {
        impl From<$t> for i32 {
            fn from(v: $t) -> i32 {
                v as i32
            }
        }
        impl Serialize for $t {
            fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                crate::structure::layer::ser_int(self, s)
            }
        }
        impl<'de> Deserialize<'de> for $t {
            fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                let v = i32::deserialize(d)?;
                Self::try_from_i32(v).map_err(serde::de::Error::custom)
            }
        }
    };
}

impl_serde_int!(crate::structure::enums::ErrorType);
impl_serde_int!(crate::structure::enums::RoughnessType);
impl_serde_int!(crate::structure::enums::ErrorMask);
impl_serde_int!(crate::structure::enums::LayerMask);
impl_serde_int!(crate::structure::enums::OptMask);
impl_serde_int!(crate::structure::enums::LayerType);
impl_serde_int!(crate::structure::enums::BlockKind);

/// One physical film: material (unresolved name), thickness [nm],
/// coherence/roughness/grading/interface flags, optimizer flags.
#[derive(Debug, Clone, PartialEq)]
pub struct Layer {
    /// Governing material name (resolved via provider at expansion).
    pub material: String,
    /// Film thickness [nm].
    pub thickness: f64,
    /// Coherent (`false` = incoherent intensity treatment).
    pub coherent: bool,
    /// Graded film (split into sub-layers for the solver).
    pub inhomogen: bool,
    /// Roughness form factor (solver contract).
    pub rough_type: RoughnessType,
    /// Grading strength driving the sub-layer refinement.
    pub inh_delta: f64,
    /// Roughness rms sigma [nm].
    pub roughness: f64,
    /// Header interface slice emitted.
    pub interface: bool,
    /// Interface slice width [nm].
    pub interface_thickness: f64,
    /// Thickness optimizable.
    pub optimize: bool,
    /// Needle-insertion site.
    pub needle: bool,
    /// Design role (ambient/film/substrate markers delimit stacks).
    pub layer_type: LayerType,
}

impl Default for Layer {
    /// Python `Layer()` defaults: 1 nm unnamed film, coherent, optimizable.
    fn default() -> Self {
        Self {
            material: String::new(),
            thickness: 1.0,
            coherent: true,
            inhomogen: false,
            rough_type: RoughnessType::None,
            inh_delta: 0.1,
            roughness: 0.0,
            interface: false,
            interface_thickness: 0.0,
            optimize: true,
            needle: true,
            layer_type: LayerType::Film,
        }
    }
}

impl Layer {
    /// Design film: `Layer::default()` with material + thickness set.
    pub fn film(thickness: f64, material: impl Into<String>) -> Self {
        Self {
            thickness,
            material: material.into(),
            ..Self::default()
        }
    }

    /// Solver sub-layer count (Python `_refine_layer_count`, transliterated
    /// exactly: `int(ceil(t^0.4) * factor) + 1`).
    ///
    /// NOTE: `powf` may differ from NumPy's power by 1 ulp; `ceil` at exact
    /// integer boundaries would then disagree. The differential suite pins
    /// counts over randomized thicknesses — any boundary divergence fails
    /// loudly there, not silently here.
    pub fn sub_layer_count(&self) -> u32 {
        if self.inhomogen && self.thickness > 0.0 {
            let factor = 1.0 + (self.inh_delta / 0.1) * 0.5;
            (self.thickness.powf(0.4).ceil() * factor) as u32 + 1
        } else {
            1
        }
    }

    /// Numeric sanity of the authored properties, as findings rather than a
    /// refusal — the caller decides what to do with them.
    ///
    /// This is the SINGLE SOURCE OF TRUTH for what a layer's numbers may be.
    /// Three doors use it: the PyO3 constructor and setters raise on any
    /// `Error` (so a bad value fails where it was written, not three calls
    /// later), [`crate::structure::structure::Structure::validate`] collects
    /// them, and the synthesis assembler checks each film it is handed. The
    /// flat-array solver surface is deliberately *not* a door — `ScatterMatrix`
    /// and the native `Solver` take raw arrays and stay permissive, as R3.1
    /// promised.
    ///
    /// `label` prefixes every message so each door can name itself:
    /// `Structure::validate` passes `"Layer 3 (TiO2)"` and keeps the exact
    /// wording it has always had (pinned by `validate_collects_like_python`);
    /// the constructor passes `"Layer"`.
    ///
    /// Non-finite values are refused rather than corrected throughout. There is
    /// no value a NaN was meant to be, and a NaN anywhere in the geometry turns
    /// every output NaN with nothing to point at.
    pub fn property_issues(&self, label: &str) -> Vec<ValidationIssue> {
        let mut issues = Vec::new();
        let bad = |m: String| ValidationIssue::error(format!("{label}: {m}"));
        let note = |m: String| ValidationIssue::warning(format!("{label}: {m}"));

        if !self.thickness.is_finite() {
            issues.push(bad(format!("Non-finite thickness {} nm.", self.thickness)));
        } else if self.thickness < 0.0 {
            issues.push(bad(format!("Negative thickness {} nm.", self.thickness)));
        }

        if !self.roughness.is_finite() {
            issues.push(bad(format!("Non-finite roughness {} nm.", self.roughness)));
        } else if self.roughness < 0.0 {
            // Not corrected to |sigma|: every roughness form factor squares it, so
            // a negative sigma is byte-identical to its positive twin and a sign
            // slip would never surface. Refusing is the only way the caller hears
            // about it. `Structure::validate` has refused this since it existed;
            // this moves the refusal to where the value is written.
            issues.push(bad(format!("Negative roughness {} nm.", self.roughness)));
        }

        if !self.interface_thickness.is_finite() {
            issues.push(bad(format!(
                "Non-finite interface thickness {} nm.",
                self.interface_thickness
            )));
        } else if self.interface_thickness < 0.0 {
            issues.push(bad(format!(
                "Negative interface thickness {} nm.",
                self.interface_thickness
            )));
        }
        // Overhang is LEGAL but suspicious: advisory, never blocking.
        if self.interface && self.interface_thickness >= self.thickness {
            issues.push(note(format!(
                "Interface thickness ({}) >= layer thickness ({}); clamped at expansion.",
                self.interface_thickness, self.thickness
            )));
        }

        // Grading strength. Expansion ramps the whole complex index by factors
        // `1 - d/2 ..= 1 + d/2`, so the physical window is [0, 2).
        if !self.inh_delta.is_finite() {
            issues.push(bad(format!("Non-finite inh_delta {}.", self.inh_delta)));
        } else if self.inh_delta < 0.0 {
            // Measured before this gate: a negative delta drives the refinement
            // factor negative, `as u32` saturates it to 0, and `sub_layer_count`
            // returns 1 — the grading is dropped and the caller gets a flat film
            // with no indication. Writing -0.2 to mean "ramp the other way" is the
            // obvious way to hit this.
            issues.push(bad(format!(
                "Negative inh_delta {}. Grading strength is a magnitude, not a \
         direction: a negative value collapses the sub-layer count to 1, so \
         the profile is dropped and the film solves as homogeneous.",
                self.inh_delta
            )));
        } else if self.inh_delta >= 2.0 {
            // Measured before this gate: inh_delta = 2.5 on a 2.35 + 0.05i film
            // put rows of n = -0.59, k = -0.0125 into the solver. Negative k is
            // gain — the layer amplifies. There is no correction worth making;
            // clamping to 1.99 only buys an index that grazes zero instead.
            issues.push(bad(format!(
                "inh_delta {} is outside [0, 2). The profile scales the complex index \
         by 1 - d/2 ..= 1 + d/2, so d >= 2 drives Re(n) through zero and flips \
         the sign of k -- which is optical gain, not a graded film.",
                self.inh_delta
            )));
        } else if self.inhomogen && self.inh_delta == 0.0 {
            issues.push(note(format!(
                "Graded layer with inh_delta 0 expands to {} identical sub-layers; \
         set inhomogen = false instead.",
                self.sub_layer_count()
            )));
        }

        issues
    }

    /// Per-layer status vector indexed by `LayerMask` (ACTIVE always 1).
    pub fn mask(&self) -> [i32; 4] {
        [
            1,
            i32::from(self.coherent),
            i32::from(self.inhomogen),
            i32::from(self.rough_type != RoughnessType::None),
        ]
    }

    /// `(material, thickness)` pair (Python `__call__`).
    pub fn as_pair(&self) -> (&str, f64) {
        (self.material.as_str(), self.thickness)
    }

    /// Bulk-set known properties; unknown/read-only/bad-enum keys become
    /// returned warnings (Python emits them via `warnings.warn`).
    /// `sub_layer_count` is derived and read-only here.
    pub fn set_properties(&mut self, props: &BTreeMap<String, Value>) -> Vec<ValidationIssue> {
        let mut warnings = Vec::new();
        for (key, value) in props {
            let bad =
                |msg: String| ValidationIssue::warning(format!("Layer.set_properties: {msg}"));
            match key.as_str() {
                "material" => match value.as_str() {
                    Some(s) => self.material = s.to_string(),
                    None => warnings.push(bad("ignoring non-string 'material'.".to_string())),
                },
                "thickness" => match value.as_f64() {
                    Some(v) => self.thickness = v,
                    None => warnings.push(bad("ignoring non-numeric 'thickness'.".to_string())),
                },
                "coherent" => match value.as_bool() {
                    Some(v) => self.coherent = v,
                    None => warnings.push(bad("ignoring non-bool 'coherent'.".to_string())),
                },
                "inhomogen" => match value.as_bool() {
                    Some(v) => self.inhomogen = v,
                    None => warnings.push(bad("ignoring non-bool 'inhomogen'.".to_string())),
                },
                "inh_delta" => match value.as_f64() {
                    Some(v) => self.inh_delta = v,
                    None => warnings.push(bad("ignoring non-numeric 'inh_delta'.".to_string())),
                },
                "roughness" => match value.as_f64() {
                    Some(v) => self.roughness = v,
                    None => warnings.push(bad("ignoring non-numeric 'roughness'.".to_string())),
                },
                "interface" => match value.as_bool() {
                    Some(v) => self.interface = v,
                    None => warnings.push(bad("ignoring non-bool 'interface'.".to_string())),
                },
                "interface_thickness" => match value.as_f64() {
                    Some(v) => self.interface_thickness = v,
                    None => warnings.push(bad(
                        "ignoring non-numeric 'interface_thickness'.".to_string()
                    )),
                },
                "optimize" => match value.as_bool() {
                    Some(v) => self.optimize = v,
                    None => warnings.push(bad("ignoring non-bool 'optimize'.".to_string())),
                },
                "needle" => match value.as_bool() {
                    Some(v) => self.needle = v,
                    None => warnings.push(bad("ignoring non-bool 'needle'.".to_string())),
                },
                "rough_type" => match value.as_i64().map(|v| v as i32) {
                    Some(v) => match RoughnessType::try_from_i32(v) {
                        Ok(e) => self.rough_type = e,
                        Err(_) => warnings.push(bad(format!("unknown rough_type {v}; ignoring."))),
                    },
                    None => warnings.push(bad("ignoring non-integer 'rough_type'.".to_string())),
                },
                "layer_type" => match value.as_i64().map(|v| v as i32) {
                    Some(v) => match LayerType::try_from_i32(v) {
                        Ok(e) => self.layer_type = e,
                        Err(_) => warnings.push(bad(format!("unknown layer_type {v}; ignoring."))),
                    },
                    None => warnings.push(bad("ignoring non-integer 'layer_type'.".to_string())),
                },
                other => warnings.push(bad(format!("ignoring unknown attribute '{other}'."))),
            }
        }
        warnings
    }
}

impl fmt::Display for Layer {
    /// Python `__repr__` format, verbatim.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Layer(mat='{}', d={:.2}nm, rough={:.2}nm, opt={})",
            self.material,
            self.thickness,
            self.roughness,
            if self.optimize { "True" } else { "False" }
        )
    }
}

impl Serialize for Layer {
    /// Python `get_state` key-for-key (`material_name`, int enums, version).
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut m = s.serialize_map(Some(13))?;
        m.serialize_entry("schema_version", &SCHEMA_VERSION)?;
        m.serialize_entry("thickness", &self.thickness)?;
        m.serialize_entry("material_name", &self.material)?;
        m.serialize_entry("coherent", &self.coherent)?;
        m.serialize_entry("inhomogen", &self.inhomogen)?;
        m.serialize_entry("inh_delta", &self.inh_delta)?;
        m.serialize_entry("rough_type", &self.rough_type)?;
        m.serialize_entry("roughness", &self.roughness)?;
        m.serialize_entry("interface", &self.interface)?;
        m.serialize_entry("interface_thickness", &self.interface_thickness)?;
        m.serialize_entry("optimize", &self.optimize)?;
        m.serialize_entry("needle", &self.needle)?;
        m.serialize_entry("layer_type", &self.layer_type)?;
        m.end()
    }
}

impl<'de> Deserialize<'de> for Layer {
    /// Python `from_state`: version-checked, unknown keys ignored.
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let map: BTreeMap<String, Value> = BTreeMap::deserialize(d)?;
        let found = map
            .get("schema_version")
            .and_then(|v| v.as_u64())
            .map(|v| v as u32);
        check_schema_version(found, "Layer").map_err(serde::de::Error::custom)?;
        let get = |k: &str| map.get(k).cloned().unwrap_or(Value::Null);
        let req_f = |k: &str| {
            get(k).as_f64().ok_or_else(|| {
                serde::de::Error::custom(format!("Layer: '{k}' missing/non-numeric"))
            })
        };
        let req_b = |k: &str| {
            get(k)
                .as_bool()
                .ok_or_else(|| serde::de::Error::custom(format!("Layer: '{k}' missing/non-bool")))
        };
        Ok(Self {
            material: map
                .get("material_name")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    serde::de::Error::custom("Layer: 'material_name' missing/non-string")
                })?
                .to_string(),
            thickness: req_f("thickness")?,
            coherent: req_b("coherent")?,
            inhomogen: req_b("inhomogen")?,
            rough_type: match get("rough_type").as_i64().map(|v| v as i32) {
                Some(v) => RoughnessType::try_from_i32(v).map_err(serde::de::Error::custom)?,
                None => {
                    return Err(serde::de::Error::custom(
                        "Layer: 'rough_type' missing/non-integer",
                    ));
                }
            },
            inh_delta: req_f("inh_delta")?,
            roughness: req_f("roughness")?,
            interface: req_b("interface")?,
            interface_thickness: req_f("interface_thickness")?,
            optimize: req_b("optimize")?,
            needle: req_b("needle")?,
            layer_type: match get("layer_type").as_i64().map(|v| v as i32) {
                Some(v) => LayerType::try_from_i32(v).map_err(serde::de::Error::custom)?,
                None => {
                    return Err(serde::de::Error::custom(
                        "Layer: 'layer_type' missing/non-integer",
                    ));
                }
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// Oracle twins (values captured from Python `Layer` live).
    #[test]
    fn refinement_counts_match_python() {
        let graded = |t: f64, d: f64| Layer {
            thickness: t,
            inh_delta: d,
            inhomogen: true,
            ..Layer::default()
        };
        assert_eq!(graded(50.0, 0.2).sub_layer_count(), 11);
        assert_eq!(graded(100.0, 0.1).sub_layer_count(), 11);
        assert_eq!(graded(57.0, 0.3).sub_layer_count(), 16);
        assert_eq!(graded(0.0, 0.1).sub_layer_count(), 1);
        assert_eq!(graded(1.0, 0.1).sub_layer_count(), 2);
        assert_eq!(graded(1000.0, 0.5).sub_layer_count(), 57);
        assert_eq!(graded(8.0, 0.0).sub_layer_count(), 4);
        assert_eq!(Layer::film(50.0, "TiO2").sub_layer_count(), 1);
    }

    #[test]
    fn mask_repr_pair_match_python() {
        let mut l = Layer::film(50.0, "TiO2");
        l.rough_type = RoughnessType::Step;
        l.roughness = 1.5;
        l.inhomogen = true;
        l.interface = true;
        assert_eq!(l.mask(), [1, 1, 1, 1]);
        assert_eq!(l.as_pair(), ("TiO2", 50.0));
        assert_eq!(
            l.to_string(),
            "Layer(mat='TiO2', d=50.00nm, rough=1.50nm, opt=True)"
        );
        let flat = Layer::film(50.0, "TiO2");
        assert_eq!(flat.mask(), [1, 1, 0, 0]);
    }

    #[test]
    fn state_round_trip_key_for_key() {
        let mut l = Layer::film(50.0, "TiO2");
        l.rough_type = RoughnessType::Step;
        l.roughness = 1.5;
        l.inhomogen = true;
        l.interface = true;
        let v = serde_json::to_value(&l).unwrap();
        let mut keys: Vec<_> = v.as_object().unwrap().keys().cloned().collect();
        keys.sort();
        assert_eq!(
            keys,
            [
                "coherent",
                "inh_delta",
                "inhomogen",
                "interface",
                "interface_thickness",
                "layer_type",
                "material_name",
                "needle",
                "optimize",
                "rough_type",
                "roughness",
                "schema_version",
                "thickness"
            ]
        );
        assert_eq!(v["rough_type"], json!(2));
        assert_eq!(v["layer_type"], json!(1));
        assert_eq!(v["schema_version"], json!(1));
        // Unknown keys ignored; version enforced.
        let mut with_bogus = v.clone();
        with_bogus["bogus"] = json!(1);
        let back: Layer = serde_json::from_value(with_bogus).unwrap();
        assert_eq!(back, l);
        let mut nover = v.clone();
        nover.as_object_mut().unwrap().remove("schema_version");
        assert!(serde_json::from_value::<Layer>(nover).is_err());
        let mut future = v.clone();
        future["schema_version"] = json!(999);
        assert!(serde_json::from_value::<Layer>(future).is_err());
        let mut bad_enum = v.clone();
        bad_enum["rough_type"] = json!(9);
        assert!(serde_json::from_value::<Layer>(bad_enum).is_err());
    }

    #[test]
    fn set_properties_warns_like_python() {
        let mut l = Layer::default();
        let mut props = BTreeMap::new();
        props.insert("thickness".to_string(), json!(25.0));
        props.insert("bogus".to_string(), json!(1));
        props.insert("rough_type".to_string(), json!(9));
        props.insert("sub_layer_count".to_string(), json!(99));
        let ws = l.set_properties(&props);
        assert_eq!(l.thickness, 25.0);
        assert_eq!(ws.len(), 3);
        assert!(ws.iter().all(|w| !w.is_error()));
    }

    #[test]
    fn property_issues_flags_each_bad_number_once() {
        let probe = |f: &dyn Fn(&mut Layer)| {
            let mut l = Layer::film(100.0, "TiO2");
            f(&mut l);
            l.property_issues("Layer")
        };
        let errs = |v: &[ValidationIssue]| -> Vec<String> {
            v.iter()
                .filter(|i| i.is_error())
                .map(|i| i.message.clone())
                .collect()
        };

        assert!(errs(&probe(&|l| l.thickness = -5.0))[0].contains("Negative thickness"));
        assert!(errs(&probe(&|l| l.thickness = f64::NAN))[0].contains("Non-finite thickness"));
        assert!(errs(&probe(&|l| l.roughness = -20.0))[0].contains("Negative roughness"));
        assert!(errs(&probe(&|l| l.interface_thickness = -2.0))[0].contains("Negative interface"));
        assert!(errs(&probe(&|l| l.inh_delta = -0.2))[0].contains("Negative inh_delta"));
        assert!(errs(&probe(&|l| l.inh_delta = 2.5))[0].contains("outside [0, 2)"));
        assert!(errs(&probe(&|l| l.inh_delta = f64::INFINITY))[0].contains("Non-finite inh_delta"));

        // 2.0 is out, 1.999 is in: the window is half-open on purpose (d = 2
        // drives Re(n) to zero).
        assert_eq!(errs(&probe(&|l| l.inh_delta = 2.0)).len(), 1);
        assert!(errs(&probe(&|l| l.inh_delta = 1.999)).is_empty());

        // Advisory, never blocking.
        let overhang = probe(&|l| {
            l.interface = true;
            l.interface_thickness = 200.0;
        });
        assert!(errs(&overhang).is_empty());
        assert_eq!(overhang.len(), 1);
        let flat_grade = probe(&|l| {
            l.inhomogen = true;
            l.inh_delta = 0.0;
        });
        assert!(errs(&flat_grade).is_empty());
        assert_eq!(flat_grade.len(), 1);

        // Ordinary layers, and the zero-thickness half-spaces, say nothing.
        assert!(probe(&|_| {}).is_empty());
        assert!(probe(&|l| l.thickness = 0.0).is_empty());
        assert!(
            probe(&|l| {
                l.inhomogen = true;
                l.inh_delta = 0.2;
            })
            .is_empty()
        );
    }

    /// These messages cross into Python and land on a cp1252 console. An em
    /// dash there raises `UnicodeEncodeError` on a bare Windows terminal, so
    /// the rule is pinned here rather than left to review.
    #[test]
    fn every_issue_message_is_ascii() {
        let mut l = Layer::film(f64::NAN, "TiO2");
        l.roughness = f64::NAN;
        l.interface_thickness = -1.0;
        l.interface = true;
        l.inh_delta = -1.0;
        let mut all: Vec<ValidationIssue> = l.property_issues("Layer");
        for d in [2.5, 0.0, f64::NAN] {
            let mut g = Layer::film(100.0, "TiO2");
            g.inhomogen = true;
            g.inh_delta = d;
            all.extend(g.property_issues("Layer"));
        }
        assert!(
            all.len() >= 6,
            "the probe stopped covering the messages: {}",
            all.len()
        );
        for i in &all {
            assert!(i.message.is_ascii(), "non-ASCII in {:?}", i.message);
        }
    }

    #[test]
    fn defaults_match_python_ctor() {
        let l = Layer::default();
        assert_eq!(l.thickness, 1.0);
        assert_eq!(l.material, "");
        assert!(l.coherent && l.optimize && l.needle);
        assert!(!l.inhomogen && !l.interface);
        assert_eq!(l.inh_delta, 0.1);
        assert_eq!(l.layer_type, LayerType::Film);
        assert_eq!(l.sub_layer_count(), 1);
    }
}
