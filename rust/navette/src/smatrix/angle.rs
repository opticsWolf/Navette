// SPDX-License-Identifier: LGPL-3.0-or-later
//! `Radians`: the one type an incidence angle has inside Navette.
//!
//! Every calculation works in radians. A degree value exists only where a
//! person types one in or reads one out, so the only ways in are the two
//! named, range-checked constructors - there is no `From<f64>`, and a bare
//! `f64` cannot become an angle by accident. `repr(transparent)` makes the
//! type free at run time: `&[Radians]` is the memory of `&[f64]`.
//!
//! The range is [0, π/2] (normal to grazing incidence). The check is what
//! catches the likeliest mistake, degrees handed over as radians: 30 rad is
//! out of range.

use std::f64::consts::FRAC_PI_2;

#[repr(transparent)]
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct Radians(f64);

impl Radians {
    /// Normal incidence (0 in either unit).
    pub const ZERO: Radians = Radians(0.0);

    /// An angle given in radians, in [0, π/2].
    pub fn new(rad: f64) -> Result<Self, String> {
        if (0.0..=FRAC_PI_2).contains(&rad) {
            Ok(Radians(rad))
        } else {
            Err(format!(
                "incidence angle {rad} rad is outside [0, pi/2] (0 to 90 degrees)"
            ))
        }
    }

    /// An angle given in degrees, in [0, 90]. One formula for the whole
    /// code base: `deg.to_radians()`, i.e. `deg * (π / 180)`.
    pub fn from_degrees(deg: f64) -> Result<Self, String> {
        if (0.0..=90.0).contains(&deg) {
            Ok(Radians(deg.to_radians()))
        } else {
            Err(format!("incidence angle {deg} deg is outside [0, 90]"))
        }
    }

    /// The angle whose sine is `sin_theta` (`asin`), for rebuilding an
    /// angle axis from the solver's precomputed sines.
    pub fn from_sin(sin_theta: f64) -> Result<Self, String> {
        Self::new(sin_theta.asin())
    }

    /// `from_sin` without the check, for sines the solver already holds:
    /// they were built from a checked grid, and a kernel that only labels
    /// or projects with them has no error path to report through.
    #[inline]
    pub(crate) fn from_checked_sin(sin_theta: f64) -> Self {
        Radians(sin_theta.asin())
    }

    #[inline(always)]
    pub fn get(self) -> f64 {
        self.0
    }

    #[inline(always)]
    pub fn sin(self) -> f64 {
        self.0.sin()
    }

    #[inline(always)]
    pub fn cos(self) -> f64 {
        self.0.cos()
    }

    /// For display only (reports, warnings, Python-facing labels).
    #[inline]
    pub fn to_degrees(self) -> f64 {
        self.0.to_degrees()
    }
}

/// A grid given in radians.
pub fn radians_slice(rad: &[f64]) -> Result<Vec<Radians>, String> {
    rad.iter().map(|&a| Radians::new(a)).collect()
}

/// A grid given in degrees.
pub fn degrees_slice(deg: &[f64]) -> Result<Vec<Radians>, String> {
    deg.iter().map(|&a| Radians::from_degrees(a)).collect()
}

/// The sines the solver kernels take.
pub fn sines(angles: &[Radians]) -> Vec<f64> {
    angles.iter().map(|a| a.sin()).collect()
}

/// The grid in degrees, for display and for the degree-labelled doors.
pub fn to_degrees_vec(angles: &[Radians]) -> Vec<f64> {
    angles.iter().map(|a| a.to_degrees()).collect()
}

/// Test shorthand: an angle literal in degrees.
#[cfg(test)]
pub(crate) fn deg(d: f64) -> Radians {
    Radians::from_degrees(d).unwrap()
}

/// Test shorthand: an angle grid literal in degrees.
#[cfg(test)]
pub(crate) fn degs(d: &[f64]) -> Vec<Radians> {
    degrees_slice(d).unwrap()
}

/// serde adapter: one angle stored in JSON as degrees.
pub mod serde_degrees {
    use super::Radians;
    use serde::{Deserialize, Deserializer, Serializer};

    pub fn serialize<S: Serializer>(a: &Radians, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_f64(a.to_degrees())
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Radians, D::Error> {
        Radians::from_degrees(f64::deserialize(d)?).map_err(serde::de::Error::custom)
    }
}

/// serde adapter: an angle grid stored in JSON as degrees.
pub mod serde_degrees_vec {
    use super::Radians;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S: Serializer>(a: &[Radians], s: S) -> Result<S::Ok, S::Error> {
        super::to_degrees_vec(a).serialize(s)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Vec<Radians>, D::Error> {
        super::degrees_slice(&Vec::<f64>::deserialize(d)?).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_range_is_normal_to_grazing_in_either_unit() {
        assert!(Radians::new(0.0).is_ok() && Radians::new(FRAC_PI_2).is_ok());
        assert!(Radians::from_degrees(0.0).is_ok() && Radians::from_degrees(90.0).is_ok());
        for bad in [-1e-12, FRAC_PI_2 + 1e-12, f64::NAN, f64::INFINITY] {
            assert!(Radians::new(bad).is_err(), "{bad}");
        }
        for bad in [-0.5, 90.5, f64::NAN] {
            assert!(Radians::from_degrees(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn degrees_given_as_radians_are_refused() {
        let err = radians_slice(&[0.0, 30.0, 60.0]).unwrap_err();
        assert!(err.contains("30 rad"), "{err}");
    }

    #[test]
    fn from_degrees_is_to_radians_bit_for_bit() {
        for i in 0..=9000 {
            let d = i as f64 / 100.0;
            assert_eq!(
                Radians::from_degrees(d).unwrap().get().to_bits(),
                d.to_radians().to_bits()
            );
        }
        assert_eq!(Radians::from_degrees(90.0).unwrap().get(), FRAC_PI_2);
    }

    #[test]
    fn a_slice_of_radians_is_a_slice_of_f64() {
        assert_eq!(std::mem::size_of::<Radians>(), std::mem::size_of::<f64>());
        assert_eq!(std::mem::align_of::<Radians>(), std::mem::align_of::<f64>());
    }
}
