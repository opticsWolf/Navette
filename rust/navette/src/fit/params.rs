// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::params — what a fit varies: addresses, bounds, scales, transforms.
//!
//! A [`ParamAddr`] names one scalar in a [`FitModel`](super::FitModel) and
//! has a textual form (configs, reports, selectors):
//!
//! | Kind | Address |
//! |---|---|
//! | thickness | `layer[main:3].thickness`, `layer[main:#TiO2_2].thickness` |
//! | grading | `layer[...].inh_delta`, `.inh_rate`, `.inh_cap` |
//! | roughness / interface | `layer[...].roughness`, `.interface_thickness` |
//! | gradient | `layer[...].gradient.f_start`, `.f_end`, `.rate`, `.f_min`, `.f_max` |
//! | group | `group[TiO2].n_factor`, `.k_factor`, `.thick_factor`, `.thick_summand`, `.inh_delta_summand`, `.roughness_summand`, `.interface_summand` |
//! | material | `material[TiO2].A`, `material[TiO2].osc[1][0]`, `material[porous].host.B1` |
//!
//! `layer[S:i]` is row `i` (0-based, ambient included) of structure `S`;
//! `layer[S:#M_k]` is the `k`-th (1-based) row of material `M` in it. A
//! material address walks the spec's plain-JSON params; a nested spec
//! (`host`, `inclusion`, `bottom`, `top`) is entered by its own keys, so
//! `host.B1` reaches `params.host.params.B1`.
//!
//! The solver never sees a parameter's own value: it sees
//! `u = (t(x) - t(x0)) / scale`, `t` the identity or `ln`. FD steps, damping
//! and tolerances then mean the same thing for a 100 nm thickness, a 0.01
//! Sellmeier `C1` and a 1e3 Urbach `alpha0`.

use std::fmt;

/// Which row of a structure a layer address names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LayerRef {
    /// Row index, 0-based, ambient row included.
    Index(usize),
    /// The `nth` (1-based) row whose material is `material`.
    Material { material: String, nth: usize },
}

/// A writable scalar on a [`Layer`](crate::structure::Layer).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LayerField {
    Thickness,
    InhDelta,
    InhRate,
    InhCap,
    Roughness,
    InterfaceThickness,
    GradientFStart,
    GradientFEnd,
    GradientRate,
    GradientFMin,
    GradientFMax,
}

const LAYER_FIELDS: [(&str, LayerField); 11] = [
    ("thickness", LayerField::Thickness),
    ("inh_delta", LayerField::InhDelta),
    ("inh_rate", LayerField::InhRate),
    ("inh_cap", LayerField::InhCap),
    ("roughness", LayerField::Roughness),
    ("interface_thickness", LayerField::InterfaceThickness),
    ("gradient.f_start", LayerField::GradientFStart),
    ("gradient.f_end", LayerField::GradientFEnd),
    ("gradient.rate", LayerField::GradientRate),
    ("gradient.f_min", LayerField::GradientFMin),
    ("gradient.f_max", LayerField::GradientFMax),
];

/// A writable scalar on a [`Group`](crate::structure::Group).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupField {
    NFactor,
    KFactor,
    ThickFactor,
    ThickSummand,
    InhDeltaSummand,
    RoughnessSummand,
    InterfaceSummand,
}

const GROUP_FIELDS: [(&str, GroupField); 7] = [
    ("n_factor", GroupField::NFactor),
    ("k_factor", GroupField::KFactor),
    ("thick_factor", GroupField::ThickFactor),
    ("thick_summand", GroupField::ThickSummand),
    ("inh_delta_summand", GroupField::InhDeltaSummand),
    ("roughness_summand", GroupField::RoughnessSummand),
    ("interface_summand", GroupField::InterfaceSummand),
];

/// One step of a material address below the material.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PathSeg {
    Key(String),
    Index(usize),
}

/// One scalar in a fit model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParamAddr {
    Layer {
        structure: String,
        layer: LayerRef,
        field: LayerField,
    },
    Group {
        group: String,
        field: GroupField,
    },
    Material {
        material: String,
        path: Vec<PathSeg>,
    },
}

/// What kind of quantity an address names (default bounds, selectors).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParamKind {
    Thickness,
    Grading,
    Roughness,
    Interface,
    Gradient,
    Group,
    Material,
}

impl ParamAddr {
    /// Parse the textual form (module docs).
    pub fn parse(text: &str) -> Result<Self, String> {
        let bad = |why: &str| format!("parameter address {text:?}: {why}");
        let (head, rest) = text
            .split_once('[')
            .ok_or_else(|| bad("expected layer[...], group[...] or material[...]"))?;
        let (inner, tail) = rest.split_once(']').ok_or_else(|| bad("unclosed '['"))?;
        if inner.is_empty() {
            return Err(bad("empty name in brackets"));
        }
        let field = tail.strip_prefix('.');
        match head {
            "layer" => {
                let (structure, row) = inner
                    .split_once(':')
                    .ok_or_else(|| bad("a layer is named as layer[structure:row]"))?;
                let layer = if let Some(spec) = row.strip_prefix('#') {
                    let (material, nth) = spec
                        .rsplit_once('_')
                        .ok_or_else(|| bad("a material row is named as #Material_k"))?;
                    let nth: usize = nth
                        .parse()
                        .ok()
                        .filter(|k| *k >= 1)
                        .ok_or_else(|| bad("the k in #Material_k counts from 1"))?;
                    LayerRef::Material {
                        material: material.to_string(),
                        nth,
                    }
                } else {
                    LayerRef::Index(row.parse().map_err(|_| bad("row is not an index"))?)
                };
                let field = field.ok_or_else(|| bad("missing .field"))?;
                let field = LAYER_FIELDS
                    .iter()
                    .find(|(name, _)| *name == field)
                    .map(|(_, f)| *f)
                    .ok_or_else(|| {
                        bad(&format!(
                            "unknown layer field; known: {}",
                            LAYER_FIELDS.map(|(n, _)| n).join(", ")
                        ))
                    })?;
                Ok(ParamAddr::Layer {
                    structure: structure.to_string(),
                    layer,
                    field,
                })
            }
            "group" => {
                let field = field.ok_or_else(|| bad("missing .field"))?;
                let field = GROUP_FIELDS
                    .iter()
                    .find(|(name, _)| *name == field)
                    .map(|(_, f)| *f)
                    .ok_or_else(|| {
                        bad(&format!(
                            "unknown group field; known: {}",
                            GROUP_FIELDS.map(|(n, _)| n).join(", ")
                        ))
                    })?;
                Ok(ParamAddr::Group {
                    group: inner.to_string(),
                    field,
                })
            }
            "material" => {
                let path = parse_path(tail).map_err(|e| bad(&e))?;
                if path.is_empty() {
                    return Err(bad("missing .param"));
                }
                if !matches!(path[0], PathSeg::Key(_)) {
                    return Err(bad("a material path starts with a parameter name"));
                }
                Ok(ParamAddr::Material {
                    material: inner.to_string(),
                    path,
                })
            }
            _ => Err(bad("expected layer[...], group[...] or material[...]")),
        }
    }

    pub fn kind(&self) -> ParamKind {
        match self {
            ParamAddr::Layer { field, .. } => match field {
                LayerField::Thickness => ParamKind::Thickness,
                LayerField::InhDelta | LayerField::InhRate | LayerField::InhCap => {
                    ParamKind::Grading
                }
                LayerField::Roughness => ParamKind::Roughness,
                LayerField::InterfaceThickness => ParamKind::Interface,
                _ => ParamKind::Gradient,
            },
            ParamAddr::Group { .. } => ParamKind::Group,
            ParamAddr::Material { .. } => ParamKind::Material,
        }
    }
}

/// `.a.b[1][0]` → segments; empty input → empty path.
fn parse_path(mut s: &str) -> Result<Vec<PathSeg>, String> {
    let mut out = Vec::new();
    while !s.is_empty() {
        if let Some(rest) = s.strip_prefix('.') {
            let end = rest.find(['.', '[']).unwrap_or(rest.len());
            if end == 0 {
                return Err("empty name in path".into());
            }
            out.push(PathSeg::Key(rest[..end].to_string()));
            s = &rest[end..];
        } else if let Some(rest) = s.strip_prefix('[') {
            let (idx, after) = rest.split_once(']').ok_or("unclosed '[' in path")?;
            out.push(PathSeg::Index(
                idx.parse()
                    .map_err(|_| format!("{idx:?} is not an index"))?,
            ));
            s = after;
        } else {
            return Err(format!("unexpected {s:?} in path"));
        }
    }
    Ok(out)
}

impl fmt::Display for ParamAddr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParamAddr::Layer {
                structure,
                layer,
                field,
            } => {
                let name = LAYER_FIELDS.iter().find(|(_, x)| x == field).unwrap().0;
                match layer {
                    LayerRef::Index(i) => write!(f, "layer[{structure}:{i}].{name}"),
                    LayerRef::Material { material, nth } => {
                        write!(f, "layer[{structure}:#{material}_{nth}].{name}")
                    }
                }
            }
            ParamAddr::Group { group, field } => {
                let name = GROUP_FIELDS.iter().find(|(_, x)| x == field).unwrap().0;
                write!(f, "group[{group}].{name}")
            }
            ParamAddr::Material { material, path } => {
                write!(f, "material[{material}]")?;
                for seg in path {
                    match seg {
                        PathSeg::Key(k) => write!(f, ".{k}")?,
                        PathSeg::Index(i) => write!(f, "[{i}]")?,
                    }
                }
                Ok(())
            }
        }
    }
}

/// How the solver coordinate relates to the parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Transform {
    #[default]
    Linear,
    /// `t(x) = ln x`: for strictly positive quantities spanning decades
    /// (oscillator widths, `alpha0`, group factors).
    Log,
}

impl Transform {
    pub fn forward(self, x: f64) -> f64 {
        match self {
            Transform::Linear => x,
            Transform::Log => x.ln(),
        }
    }

    pub fn inverse(self, t: f64) -> f64 {
        match self {
            Transform::Linear => t,
            Transform::Log => t.exp(),
        }
    }

    /// `dx/dt` at `x`.
    pub fn slope(self, x: f64) -> f64 {
        match self {
            Transform::Linear => 1.0,
            Transform::Log => x,
        }
    }
}

/// One fitted (or held) scalar. Unset fields take the defaults of its
/// kind when the problem is built; an unset start reads the model.
#[derive(Debug, Clone, PartialEq)]
pub struct Parameter {
    pub addr: ParamAddr,
    /// Start value; `None` reads the model.
    pub value: Option<f64>,
    pub lb: Option<f64>,
    pub ub: Option<f64>,
    /// Typical size of a change, in the transformed coordinate.
    pub scale: Option<f64>,
    pub transform: Transform,
    /// Written once, never varied.
    pub fixed: bool,
    /// Free-form labels for selectors (phase 3).
    pub tags: Vec<String>,
}

impl Parameter {
    pub fn new(addr: ParamAddr) -> Self {
        Self {
            addr,
            value: None,
            lb: None,
            ub: None,
            scale: None,
            transform: Transform::Linear,
            fixed: false,
            tags: Vec::new(),
        }
    }

    /// [`Parameter::new`] from the textual address.
    pub fn parse(addr: &str) -> Result<Self, String> {
        ParamAddr::parse(addr).map(Self::new)
    }

    pub fn start(mut self, value: f64) -> Self {
        self.value = Some(value);
        self
    }

    pub fn bounds(mut self, lb: f64, ub: f64) -> Self {
        self.lb = Some(lb);
        self.ub = Some(ub);
        self
    }

    pub fn scale(mut self, scale: f64) -> Self {
        self.scale = Some(scale);
        self
    }

    pub fn log(mut self) -> Self {
        self.transform = Transform::Log;
        self
    }

    pub fn fixed(mut self) -> Self {
        self.fixed = true;
        self
    }

    pub fn tag(mut self, tag: impl Into<String>) -> Self {
        self.tags.push(tag.into());
        self
    }
}

/// Default `[lb, ub]` of an address's kind: the gates `Layer::validate`,
/// `InhMode::issues`, `GradientSpec::issues` and `Group::validate` enforce.
/// Material parameters are unbounded until phase 2's per-model table.
pub(crate) fn default_bounds(addr: &ParamAddr, x0: f64) -> (f64, f64) {
    const OPEN: f64 = 1e-9;
    match addr {
        ParamAddr::Layer { field, .. } => match field {
            LayerField::Thickness => (0.0, (10.0 * x0).max(1000.0)),
            LayerField::InhDelta => (-2.0 + OPEN, 2.0 - OPEN),
            LayerField::InhCap
            | LayerField::GradientFStart
            | LayerField::GradientFEnd
            | LayerField::GradientFMin
            | LayerField::GradientFMax => (0.0, 1.0),
            LayerField::Roughness | LayerField::InterfaceThickness => (0.0, f64::INFINITY),
            LayerField::InhRate | LayerField::GradientRate => (f64::NEG_INFINITY, f64::INFINITY),
        },
        ParamAddr::Group { field, .. } => match field {
            GroupField::NFactor | GroupField::KFactor | GroupField::ThickFactor => {
                (f64::MIN_POSITIVE, f64::INFINITY)
            }
            _ => (f64::NEG_INFINITY, f64::INFINITY),
        },
        ParamAddr::Material { .. } => (f64::NEG_INFINITY, f64::INFINITY),
    }
}

/// Default scale: 1 % of the start in the linear coordinate (at least
/// 1e-3, so a zero start still moves), 0.01 in the log coordinate.
pub(crate) fn default_scale(transform: Transform, x0: f64) -> f64 {
    match transform {
        Transform::Linear => (0.01 * x0.abs()).max(1e-3),
        Transform::Log => 0.01,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addresses_round_trip_through_text() {
        for text in [
            "layer[main:3].thickness",
            "layer[main:#TiO2_2].thickness",
            "layer[coat:1].gradient.f_end",
            "layer[coat:0].inh_cap",
            "group[TiO2].n_factor",
            "group[SiO2].interface_summand",
            "material[TiO2].A",
            "material[TiO2].osc[1][0]",
            "material[porous].host.B1",
        ] {
            let addr = ParamAddr::parse(text).unwrap();
            assert_eq!(addr.to_string(), text);
        }
    }

    #[test]
    fn material_row_names_keep_underscores() {
        let addr = ParamAddr::parse("layer[s:#Ta2O5_hi_3].thickness").unwrap();
        let ParamAddr::Layer { layer, .. } = addr else {
            panic!()
        };
        assert_eq!(
            layer,
            LayerRef::Material {
                material: "Ta2O5_hi".into(),
                nth: 3
            }
        );
    }

    #[test]
    fn malformed_addresses_name_themselves() {
        for (text, why) in [
            ("layer[main].thickness", "layer[structure:row]"),
            ("layer[main:x].thickness", "not an index"),
            ("layer[main:#TiO2_0].thickness", "counts from 1"),
            ("layer[main:1].colour", "unknown layer field"),
            ("group[H].n", "unknown group field"),
            ("material[TiO2]", "missing .param"),
            ("material[TiO2].osc[x]", "not an index"),
            ("stack[1].thickness", "expected layer"),
            ("layer[main:1", "unclosed"),
        ] {
            let err = ParamAddr::parse(text).unwrap_err();
            assert!(err.contains(text) && err.contains(why), "{err}");
        }
    }

    #[test]
    fn log_transform_inverts() {
        let t = Transform::Log;
        assert!((t.inverse(t.forward(37.5)) - 37.5).abs() < 1e-12);
        assert_eq!(t.slope(2.0), 2.0);
    }
}
