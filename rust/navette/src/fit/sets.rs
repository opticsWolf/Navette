// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::sets — selectors and parameter sets (plan §3.3.4, §3.3.5).
//!
//! A selector picks parameters out of a problem's declared list:
//!
//! | Text | Picks |
//! |---|---|
//! | `*` | every parameter |
//! | `kind=thickness` | by kind: `thickness`, `grading`, `roughness`, `interface`, `gradient`, `group`, `material` |
//! | `tag=H` | parameters carrying the tag |
//! | `set:dispersion_H` | the members of a set |
//! | anything else | a glob on the address: `material[TiO2].*`, `group[*].n_factor`; `*` is any run, `?` one character, everything else literal |
//!
//! A **parameter set** is a named union of selectors: the unit a strategy
//! orders. Built-in sets need no declaration: `thicknesses`, `grading`,
//! `roughness`, `interfaces`, `gradients`, `groups`, `materials` (one per
//! kind), and `material[X]` / `group[G]` for every material and group a
//! parameter names. A user set of the same name wins.

use serde::{Deserialize, Serialize};

use super::params::{ParamKind, Parameter};

/// A named union of selectors.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParamSet {
    pub name: String,
    pub members: Vec<String>,
}

impl ParamSet {
    pub fn new(name: impl Into<String>, members: &[&str]) -> Self {
        Self {
            name: name.into(),
            members: members.iter().map(|s| s.to_string()).collect(),
        }
    }
}

/// One parsed selector.
#[derive(Debug, Clone, PartialEq)]
pub enum Selector {
    All,
    Kind(ParamKind),
    Tag(String),
    Set(String),
    Glob(String),
}

impl Selector {
    pub fn parse(text: &str) -> Result<Self, String> {
        let text = text.trim();
        if text.is_empty() {
            return Err("empty selector".into());
        }
        if text == "*" {
            return Ok(Selector::All);
        }
        if let Some(kind) = text.strip_prefix("kind=") {
            return parse_kind(kind).map(Selector::Kind);
        }
        if let Some(tag) = text.strip_prefix("tag=") {
            return Ok(Selector::Tag(tag.to_string()));
        }
        if let Some(set) = text.strip_prefix("set:") {
            return Ok(Selector::Set(set.to_string()));
        }
        Ok(Selector::Glob(text.to_string()))
    }
}

fn parse_kind(kind: &str) -> Result<ParamKind, String> {
    Ok(match kind {
        "thickness" => ParamKind::Thickness,
        "grading" => ParamKind::Grading,
        "roughness" => ParamKind::Roughness,
        "interface" => ParamKind::Interface,
        "gradient" => ParamKind::Gradient,
        "group" => ParamKind::Group,
        "material" => ParamKind::Material,
        other => {
            return Err(format!(
                "unknown kind {other:?}; one of thickness, grading, roughness, interface, \
                 gradient, group, material"
            ));
        }
    })
}

/// `*` any run, `?` one character, the rest literal.
fn glob(pattern: &[char], text: &[char]) -> bool {
    match pattern.split_first() {
        None => text.is_empty(),
        Some(('*', rest)) => (0..=text.len()).any(|i| glob(rest, &text[i..])),
        Some(('?', rest)) => !text.is_empty() && glob(rest, &text[1..]),
        Some((c, rest)) => text.first() == Some(c) && glob(rest, &text[1..]),
    }
}

/// The built-in sets' selectors, or `None`.
fn builtin(name: &str) -> Option<Vec<Selector>> {
    let kind = match name {
        "thicknesses" => Some(ParamKind::Thickness),
        "grading" => Some(ParamKind::Grading),
        "roughness" => Some(ParamKind::Roughness),
        "interfaces" => Some(ParamKind::Interface),
        "gradients" => Some(ParamKind::Gradient),
        "groups" => Some(ParamKind::Group),
        "materials" => Some(ParamKind::Material),
        _ => None,
    };
    if let Some(k) = kind {
        return Some(vec![Selector::Kind(k)]);
    }
    if (name.starts_with("material[") || name.starts_with("group[")) && name.ends_with(']') {
        return Some(vec![Selector::Glob(format!("{name}.*"))]);
    }
    None
}

/// Resolve selectors against `parameters`: one flag per parameter.
pub fn select(
    selectors: &[String],
    parameters: &[Parameter],
    sets: &[ParamSet],
) -> Result<Vec<bool>, String> {
    let mut picked = vec![false; parameters.len()];
    for text in selectors {
        let sel = Selector::parse(text)?;
        let hits = apply(&sel, parameters, sets, &mut Vec::new())?;
        if !hits.iter().any(|h| *h) {
            return Err(format!("selector {text:?} matches no parameter"));
        }
        for (p, h) in picked.iter_mut().zip(hits) {
            *p |= h;
        }
    }
    Ok(picked)
}

fn apply(
    sel: &Selector,
    parameters: &[Parameter],
    sets: &[ParamSet],
    visiting: &mut Vec<String>,
) -> Result<Vec<bool>, String> {
    Ok(match sel {
        Selector::All => vec![true; parameters.len()],
        Selector::Kind(k) => parameters.iter().map(|p| p.addr.kind() == *k).collect(),
        Selector::Tag(t) => parameters.iter().map(|p| p.tags.contains(t)).collect(),
        Selector::Glob(g) => {
            let pat: Vec<char> = g.chars().collect();
            parameters
                .iter()
                .map(|p| glob(&pat, &p.addr.to_string().chars().collect::<Vec<_>>()))
                .collect()
        }
        Selector::Set(name) => {
            if visiting.contains(name) {
                return Err(format!("set {name:?} contains itself"));
            }
            let members = match sets.iter().find(|s| s.name == *name) {
                Some(s) => s
                    .members
                    .iter()
                    .map(|m| Selector::parse(m))
                    .collect::<Result<Vec<_>, _>>()?,
                None => builtin(name).ok_or_else(|| {
                    format!(
                        "no set {name:?}; declared {:?}, built-in thicknesses, grading, roughness, \
                         interfaces, gradients, groups, materials, material[X], group[G]",
                        sets.iter().map(|s| &s.name).collect::<Vec<_>>()
                    )
                })?,
            };
            visiting.push(name.clone());
            let mut out = vec![false; parameters.len()];
            for m in &members {
                for (o, h) in out.iter_mut().zip(apply(m, parameters, sets, visiting)?) {
                    *o |= h;
                }
            }
            visiting.pop();
            out
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> Vec<Parameter> {
        [
            "layer[main:1].thickness",
            "layer[main:2].thickness",
            "group[TiO2].n_factor",
            "group[SiO2].n_factor",
            "material[TiO2].A",
            "material[TiO2].osc[0][1]",
        ]
        .iter()
        .map(|a| Parameter::parse(a).unwrap())
        .collect()
    }

    fn picked(sel: &[&str], sets: &[ParamSet]) -> Vec<usize> {
        let sel: Vec<String> = sel.iter().map(|s| s.to_string()).collect();
        select(&sel, &params(), sets)
            .unwrap()
            .iter()
            .enumerate()
            .filter(|(_, p)| **p)
            .map(|(i, _)| i)
            .collect()
    }

    #[test]
    fn selectors_pick_by_kind_glob_tag_and_set() {
        assert_eq!(picked(&["kind=thickness"], &[]), [0, 1]);
        assert_eq!(picked(&["group[*].n_factor"], &[]), [2, 3]);
        assert_eq!(picked(&["material[TiO2].*"], &[]), [4, 5]);
        assert_eq!(picked(&["*"], &[]).len(), 6);
        assert_eq!(picked(&["set:material[TiO2]"], &[]), [4, 5]);
        assert_eq!(
            picked(&["set:thicknesses", "group[TiO2].*"], &[]),
            [0, 1, 2]
        );
        let h = ParamSet::new("H_all", &["material[TiO2].*", "group[TiO2].*"]);
        assert_eq!(picked(&["set:H_all"], &[h]), [2, 4, 5]);
        let mut p = params();
        p[1] = p[1].clone().tag("L");
        let flags = select(&["tag=L".into()], &p, &[]).unwrap();
        assert_eq!(flags, [false, true, false, false, false, false]);
    }

    #[test]
    fn bad_selectors_say_why() {
        let p = params();
        let err = |s: &str, sets: &[ParamSet]| select(&[s.to_string()], &p, sets).unwrap_err();
        assert!(err("kind=colour", &[]).contains("unknown kind"));
        assert!(err("material[SiO2].*", &[]).contains("matches no parameter"));
        assert!(err("set:nope", &[]).contains("no set \"nope\""));
        let loop_ = ParamSet::new("a", &["set:a"]);
        assert!(err("set:a", &[loop_]).contains("contains itself"));
    }
}
