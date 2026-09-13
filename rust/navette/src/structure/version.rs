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
pub const SCHEMA_VERSION: u32 = 2;

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
        assert_eq!(SCHEMA_VERSION, 2);
        assert_eq!(MIN_READABLE_SCHEMA_VERSION, 1);
        // The whole readable range is accepted - the F1.4 contract.
        assert!(check_schema_version(Some(1), "Layer").is_ok());
        assert!(check_schema_version(Some(2), "Layer").is_ok());
        assert!(check_schema_version(None, "Layer").is_err());
        // Below the range: stale.
        let v0 = check_schema_version(Some(0), "Layer").unwrap_err();
        assert!(v0.contains("stale"), "{v0}");
        // Above the range: newer build, named as such.
        let v3 = check_schema_version(Some(3), "Layer").unwrap_err();
        assert!(v3.contains("newer build"), "{v3}");
        let v1000 = check_schema_version(Some(1000), "Layer").unwrap_err();
        assert!(v1000.contains("newer build"), "{v1000}");
    }
}
