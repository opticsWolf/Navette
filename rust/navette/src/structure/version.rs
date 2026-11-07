// SPDX-License-Identifier: LGPL-3.0-or-later
//! State-schema version gate.
//!
//! F1.4 widens the point gate to a readable RANGE: every state this code
//! can still read truthfully is accepted, everything else is refused with
//! a reason that names the direction of the mismatch. Mirrors
//! `navette.structure.types.SCHEMA_VERSION` / `MIN_READABLE_SCHEMA_VERSION`.
//!
//! Why a range (F1.4): the first additive feature keys (`gradient`,
//! schema v2) are written with defaults that reconstruct on read, so a v1
//! state written by a pre-F1.4 build is still *truthfully* readable - the
//! reader reconstructs what the older writer left implicit. Refusing it
//! would discard user data for no misread. The range closes again at the
//! top: a state naming a version above this build's writer was written by
//! NEWER code whose keys/meanings this build cannot know - that refusal
//! is a data-integrity refusal, not a stale-state refusal, and the two
//! get different messages so a user can tell them apart.

/// Current state-schema version (structure states + config states).
///
/// v3 flips the grading sign to the literature's: `+inh_delta` now raises
/// the index toward the outer surface. A v1/v2 state still reads
/// truthfully - the reader negates its grading keys (see
/// [`legacy_grading_sign`]) - so the readable range stays open below.
pub const SCHEMA_VERSION: u32 = 3;

/// First schema whose grading keys carry the literature sign.
pub const GRADING_SIGN_VERSION: u32 = 3;

/// Whether a state of version `found` wrote its grading keys
/// (`inh_delta`, the `RateCapped` rate, `inh_delta_summand`, the
/// `inh_delta` error means) with the pre-v3 sign, so the reader must
/// negate them to keep the profile they described.
pub(crate) fn legacy_grading_sign(found: Option<u32>) -> bool {
    found.is_some_and(|v| v < GRADING_SIGN_VERSION)
}

/// Negate a legacy grading value; zero stays `+0.0`, so a v2 state with
/// no grading re-serializes without a `-0.0`.
pub(crate) fn flip_grading(x: f64) -> f64 {
    if x == 0.0 { x } else { -x }
}

/// Oldest state version this build still reads truthfully. Bumped only
/// when an older schema becomes UNREADABLE (a key whose meaning changed
/// without a new name); additive keys lower it never.
pub const MIN_READABLE_SCHEMA_VERSION: u32 = 1;

/// Refuse states outside this build's readable range.
///
/// Two distinct refusals, and each names its direction: below the range
/// is *stale* (written by an older build that knew less than we do -
/// refuse, and say so); above the range is *newer* (written by a build
/// that knows more than we do - refuse, and say a newer build wrote
/// it, because the remedy is upgrading, not hand-editing the file).
pub fn check_schema_version(found: Option<u32>, what: &str) -> Result<(), String> {
    match found {
        None => Err(format!(
            "{what}: missing schema_version tag (malformed state)."
        )),
        Some(v) if (MIN_READABLE_SCHEMA_VERSION..=SCHEMA_VERSION).contains(&v) => Ok(()),
        Some(v) if v > SCHEMA_VERSION => Err(format!(
            "{what}: schema_version {v} unsupported \
       (this code reads up to {SCHEMA_VERSION}); refusing a state \
       written by a newer build - upgrade navette to read it."
        )),
        Some(v) => Err(format!(
            "{what}: schema_version {v} unsupported \
       (this code reads up to {SCHEMA_VERSION}); refusing a stale state."
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_gate_reads_a_range() {
        assert_eq!(SCHEMA_VERSION, 3);
        assert_eq!(MIN_READABLE_SCHEMA_VERSION, 1);
        // The whole readable range is accepted - the F1.4 contract.
        assert!(check_schema_version(Some(1), "Layer").is_ok());
        assert!(check_schema_version(Some(2), "Layer").is_ok());
        assert!(check_schema_version(Some(3), "Layer").is_ok());
        assert!(check_schema_version(None, "Layer").is_err());
        // Below the range: stale.
        let v0 = check_schema_version(Some(0), "Layer").unwrap_err();
        assert!(v0.contains("stale"), "{v0}");
        // Above the range: newer build, named as such.
        let v4 = check_schema_version(Some(4), "Layer").unwrap_err();
        assert!(v4.contains("newer build"), "{v4}");
        let v1000 = check_schema_version(Some(1000), "Layer").unwrap_err();
        assert!(v1000.contains("newer build"), "{v1000}");
    }

    #[test]
    fn only_pre_v3_states_carry_the_old_grading_sign() {
        assert!(legacy_grading_sign(Some(1)));
        assert!(legacy_grading_sign(Some(2)));
        assert!(!legacy_grading_sign(Some(3)));
        assert!(!legacy_grading_sign(None));
        assert_eq!(flip_grading(0.2), -0.2);
        assert!(flip_grading(0.0).is_sign_positive());
    }
}
