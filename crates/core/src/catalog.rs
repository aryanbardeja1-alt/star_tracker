//! The star catalogue: ICRS/J2000 unit vectors with V magnitudes.
//!
//! `data/catalog.bin` is produced by the `build-catalog` tool and is loaded
//! here from bytes, so the same code serves the native CLI and wasm-bindgen
//! (via `include_bytes!` in the WASM crate). This module performs no I/O.

use crate::math::Vec3;
use crate::{Error, Result};
use serde::{Deserialize, Serialize};

/// Tolerance accepted on `|unit| - 1` when decoding a catalogue.
const UNIT_TOLERANCE: f64 = 1e-9;

/// One catalogue entry.
///
/// A blend -- two or more stars closer together than the blend radius used by
/// `build-catalog` -- is stored as a single `Star` carrying the combined flux
/// and the flux-weighted direction. `id` and `name` are then those of the
/// brightest component.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Star {
    /// Hipparcos (HIP) catalogue number.
    pub id: u32,
    /// Direction in the inertial frame (ICRS, epoch J2000.0) as a unit vector,
    /// `r = (cos d cos a, cos d sin a, sin d)`.
    pub unit: [f64; 3],
    /// Johnson V magnitude; for a blend, the combined magnitude.
    pub mag: f32,
    /// IAU-approved proper name, where the star has one.
    pub name: Option<String>,
}

impl Star {
    /// The star's direction as an inertial (ICRS/J2000) unit vector.
    pub fn direction(&self) -> Vec3 {
        Vec3::new(self.unit[0], self.unit[1], self.unit[2])
    }
}

/// A loaded star catalogue.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Catalog {
    /// Entries ordered by increasing magnitude, ties broken by `id`, so that
    /// any magnitude-limited subset is a prefix of this slice.
    pub stars: Vec<Star>,
}

impl Catalog {
    /// Decodes a catalogue from the bincode bytes written by `build-catalog`.
    ///
    /// Validates that the catalogue is non-empty and that every direction is a
    /// finite unit vector, so a stale or truncated `catalog.bin` fails loudly
    /// here instead of silently skewing every attitude downstream.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let (catalog, _consumed): (Self, usize) =
            bincode::serde::decode_from_slice(bytes, bincode::config::standard())
                .map_err(|e| Error::CatalogDecode(e.to_string()))?;

        if catalog.stars.is_empty() {
            return Err(Error::CatalogEmpty);
        }
        for star in &catalog.stars {
            let norm = star.direction().norm();
            if !norm.is_finite() || (norm - 1.0).abs() > UNIT_TOLERANCE {
                return Err(Error::CatalogNotUnit { id: star.id, norm });
            }
        }
        Ok(catalog)
    }

    /// Builds a lookup from catalogue id to catalogue index.
    ///
    /// The catalogue is ordered by magnitude, not by id, so recall statistics
    /// -- which have to ask whether a truth star was in the database at all --
    /// would otherwise need a linear scan per star.
    pub fn id_index(&self) -> IdIndex {
        let mut sorted: Vec<(u32, u32)> = self
            .stars
            .iter()
            .enumerate()
            .map(|(index, star)| (star.id, index as u32))
            .collect();
        sorted.sort_unstable_by_key(|(id, _)| *id);
        IdIndex { sorted }
    }

    /// Encodes the catalogue to bincode bytes.
    ///
    /// Lives here rather than in `build-catalog` so that encode and decode
    /// cannot drift apart in configuration.
    pub fn to_bytes(&self) -> Result<Vec<u8>> {
        bincode::serde::encode_to_vec(self, bincode::config::standard())
            .map_err(|e| Error::CatalogEncode(e.to_string()))
    }
}

/// A lookup from catalogue id to catalogue index, built by
/// [`Catalog::id_index`].
#[derive(Clone, Debug)]
pub struct IdIndex {
    /// `(id, index)` ordered by id.
    sorted: Vec<(u32, u32)>,
}

impl IdIndex {
    /// Catalogue index of the star with this id, if the catalogue holds it.
    pub fn get(&self, id: u32) -> Option<usize> {
        self.sorted
            .binary_search_by_key(&id, |(key, _)| *key)
            .ok()
            .map(|at| self.sorted[at].1 as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{angle_between, unit_to_radec};
    use std::f64::consts::PI;

    /// One arcsecond in radians.
    const ARCSEC: f64 = PI / (180.0 * 3600.0);

    /// The committed catalogue, embedded at compile time so the test needs no
    /// I/O and exercises exactly the bytes the WASM crate will carry.
    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

    fn catalog() -> Catalog {
        Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode")
    }

    /// The blend radius `build-catalog` applies: 6 px at the default 20 deg
    /// field over 1024 px, i.e. `6 * FOV / W`.
    const BLEND_RADIUS: f64 = 6.0 * (20.0 * PI / 180.0) / 1024.0;

    #[test]
    fn catalog_bin_is_within_the_size_budget() {
        assert!(
            CATALOG_BIN.len() < 1_000_000,
            "catalog.bin is {} bytes, over the 1 MB budget",
            CATALOG_BIN.len()
        );
    }

    #[test]
    fn star_count_is_plausible() {
        let n = catalog().stars.len();
        assert!(
            (8000..=9500).contains(&n),
            "expected roughly 9k stars at V <= 6.5, got {n}"
        );
    }

    #[test]
    fn every_entry_is_sane() {
        let catalog = catalog();
        for star in &catalog.stars {
            assert!(
                star.mag <= 6.5,
                "star {} has V = {} above the 6.5 cut",
                star.id,
                star.mag
            );
            assert!(star.mag.is_finite());
            assert!(star.id > 0, "HIP numbers start at 1");
            // `from_bytes` already checked the norm; confirm the components
            // themselves are finite rather than cancelling to a unit length.
            assert!(star.unit.iter().all(|c| c.is_finite()));
        }
    }

    #[test]
    fn ids_are_unique() {
        let mut ids: Vec<u32> = catalog().stars.iter().map(|s| s.id).collect();
        ids.sort_unstable();
        let before = ids.len();
        ids.dedup();
        assert_eq!(before, ids.len(), "duplicate HIP ids in the catalogue");
    }

    #[test]
    fn stars_are_sorted_brightest_first() {
        let catalog = catalog();
        for pair in catalog.stars.windows(2) {
            let (a, b) = (&pair[0], &pair[1]);
            assert!(
                (a.mag, a.id) <= (b.mag, b.id),
                "order broken at {} ({}) before {} ({})",
                a.id,
                a.mag,
                b.id,
                b.mag
            );
        }
    }

    /// The post-merge invariant: nothing closer than the blend radius survives
    /// as two separate entries, or the identification stage would be handed a
    /// pair it can never resolve.
    #[test]
    fn no_pair_is_closer_than_the_blend_radius() {
        let catalog = catalog();
        let mut by_z: Vec<&Star> = catalog.stars.iter().collect();
        by_z.sort_by(|a, b| a.unit[2].total_cmp(&b.unit[2]));

        let mut closest = f64::INFINITY;
        for (i, a) in by_z.iter().enumerate() {
            for b in &by_z[i + 1..] {
                // |sin d_b - sin d_a| <= |d_b - d_a| <= angle, so once the z
                // gap exceeds the radius no later star can be inside it.
                if b.unit[2] - a.unit[2] > BLEND_RADIUS {
                    break;
                }
                closest = closest.min(angle_between(&a.direction(), &b.direction()));
            }
        }
        assert!(
            closest > BLEND_RADIUS,
            "closest surviving pair is {:.1}\" but the blend radius is {:.1}\"",
            closest / ARCSEC,
            BLEND_RADIUS / ARCSEC
        );
    }

    /// Acceptance check for Phase 1: five bright stars against their published
    /// ICRS J2000.0 positions.
    ///
    /// References are VizieR's own proper-motion-propagated `_RAJ2000` /
    /// `_DEJ2000` for Hipparcos I/239, each of which reproduces the standard
    /// published position (Sirius 06h45m08.92s -16d42'58.0", Vega
    /// 18h36m56.34s +38d47'01.3"). They are independent of the propagation
    /// `build-catalog` performs, which works in vector form instead.
    #[test]
    fn bright_stars_match_known_positions_within_one_arcsecond() {
        // (HIP, name, RA J2000 deg, Dec J2000 deg)
        const REFERENCES: [(u32, &str, f64, f64); 5] = [
            (32349, "Sirius", 101.2871553865, -16.7161158193),
            (91262, "Vega", 279.2347351065, 38.7836917958),
            (69673, "Arcturus", 213.9153001021, 19.1824102958),
            (27989, "Betelgeuse", 88.7929385961, 7.4070627358),
            (11767, "Polaris", 37.9545153525, 89.2641095074),
        ];

        let catalog = catalog();
        for (hip, name, ra_deg, dec_deg) in REFERENCES {
            let star = catalog
                .stars
                .iter()
                .find(|s| s.id == hip)
                .unwrap_or_else(|| panic!("HIP {hip} ({name}) missing from the catalogue"));

            let expected = crate::math::radec_to_unit(ra_deg.to_radians(), dec_deg.to_radians());
            let separation = angle_between(&star.direction(), &expected);
            assert!(
                separation < ARCSEC,
                "HIP {hip} ({name}) is {:.3}\" from its published position",
                separation / ARCSEC
            );

            // The stored position must also round-trip through RA/Dec.
            let (ra, dec) = unit_to_radec(&star.direction());
            assert!((ra.to_degrees() - ra_deg).abs() < 1e-3);
            assert!((dec.to_degrees() - dec_deg).abs() < 1e-3);
        }
    }

    #[test]
    fn proper_names_are_populated() {
        let catalog = catalog();
        let named = catalog.stars.iter().filter(|s| s.name.is_some()).count();
        assert!(
            named > 300,
            "expected 400-odd IAU names to survive the V <= 6.5 cut, got {named}"
        );

        for (hip, want) in [(32349, "Sirius"), (91262, "Vega"), (11767, "Polaris")] {
            let star = catalog
                .stars
                .iter()
                .find(|s| s.id == hip)
                .unwrap_or_else(|| panic!("HIP {hip} missing"));
            assert_eq!(star.name.as_deref(), Some(want));
        }
    }

    #[test]
    fn round_trips_through_bytes() {
        let catalog = catalog();
        let bytes = catalog.to_bytes().expect("encode");
        let back = Catalog::from_bytes(&bytes).expect("decode");
        assert_eq!(catalog, back);
        assert_eq!(bytes.as_slice(), CATALOG_BIN, "encoding is not stable");
    }

    #[test]
    fn rejects_an_empty_catalogue() {
        let empty = Catalog::default();
        let bytes = empty.to_bytes().expect("encode");
        assert!(matches!(
            Catalog::from_bytes(&bytes),
            Err(Error::CatalogEmpty)
        ));
    }

    #[test]
    fn rejects_a_non_unit_vector() {
        let bad = Catalog {
            stars: vec![Star {
                id: 1,
                unit: [1.0, 1.0, 0.0],
                mag: 3.0,
                name: None,
            }],
        };
        let bytes = bad.to_bytes().expect("encode");
        assert!(matches!(
            Catalog::from_bytes(&bytes),
            Err(Error::CatalogNotUnit { id: 1, .. })
        ));
    }

    #[test]
    fn rejects_garbage_bytes() {
        assert!(matches!(
            Catalog::from_bytes(&[0xff; 8]),
            Err(Error::CatalogDecode(_))
        ));
    }
}
