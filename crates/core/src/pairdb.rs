//! Pair database with a k-vector index.
//!
//! Every catalogue pair whose separation could fit inside the field is stored
//! once, sorted by angle, alongside a k-vector index (Mortari) that turns a
//! "which pairs have a separation near this angle?" query into O(1) arithmetic
//! plus a slice.
//!
//! Indices in a [`Pair`] are **catalogue** indices, never indices into a frame's
//! observed stars. Because the catalogue is sorted brightest first, the stars
//! with `V <= db_mag_cut` are exactly its leading prefix, so a pair index is
//! also a direct index into `Catalog::stars`.
//!
//! Angles are radians. Separations are stored as `f32` as docs/SPEC.md specifies,
//! which at a 0.5 rad separation resolves to about 0.006 arcsec -- three orders
//! finer than the search tolerance -- while halving the size of a table that
//! holds the better part of a million entries.

use crate::catalog::Catalog;
use crate::math::{Vec3, angle_between};
use crate::simulate::SimConfig;
use crate::{Error, Result};

/// Margin added to the field diagonal when deciding which pairs to keep.
///
/// A pair exactly on the diagonal can still be measured slightly wider than
/// the ideal geometry predicts, once the focal-length error and distortion the
/// presets inject are in play, so the table has to reach a little past it.
const DIAGONAL_MARGIN_DEG: f64 = 1.0;

/// A pair is indexed by `u16`, so this is the most database stars possible.
const MAX_DB_STARS: usize = u16::MAX as usize + 1;

/// One catalogue pair and its angular separation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pair {
    /// Angular separation in radians.
    pub angle: f32,
    /// Catalogue index of the first star; always less than `j`.
    pub i: u16,
    /// Catalogue index of the second star.
    pub j: u16,
}

impl Pair {
    /// The pair's two catalogue indices as `usize`.
    pub fn indices(&self) -> (usize, usize) {
        (usize::from(self.i), usize::from(self.j))
    }

    /// Whether `index` is one of this pair's two stars.
    pub fn contains(&self, index: u16) -> bool {
        self.i == index || self.j == index
    }

    /// The pair's other star, given one of them.
    pub fn other(&self, index: u16) -> Option<u16> {
        if self.i == index {
            Some(self.j)
        } else if self.j == index {
            Some(self.i)
        } else {
            None
        }
    }
}

/// All catalogue pairs within the field, sorted by separation, plus the
/// k-vector that indexes them.
#[derive(Clone, Debug)]
pub struct PairDb {
    /// Pairs in ascending order of `angle`.
    pairs: Vec<Pair>,
    /// `k[i]` is how many pairs have an angle below `q + m*i`. Monotonic, with
    /// `k[0] == 0` and `k[n-1] == n`.
    k: Vec<u32>,
    /// Slope of the k-vector's index line, radians per entry.
    m: f64,
    /// Intercept of that line, radians.
    q: f64,
    /// How many catalogue stars went into the database.
    star_count: usize,
    /// Largest separation the table holds.
    max_angle: f64,
}

impl PairDb {
    /// Builds the database from the leading `V <= cfg.db_mag_cut` prefix of
    /// `catalog`, keeping separations out to the field diagonal plus a margin.
    ///
    /// Fails only if the magnitude cut admits more stars than a `u16` index can
    /// address.
    pub fn build(catalog: &Catalog, cfg: &SimConfig) -> Result<Self> {
        // The catalogue is sorted by magnitude, so the cut is a prefix.
        let star_count = catalog
            .stars
            .partition_point(|star| star.mag <= cfg.db_mag_cut);
        if star_count > MAX_DB_STARS {
            return Err(Error::TooManyDatabaseStars { count: star_count });
        }
        let max_angle =
            2.0 * cfg.nominal_camera().half_diagonal_fov() + DIAGONAL_MARGIN_DEG.to_radians();

        let pairs = collect_pairs(catalog, star_count, max_angle);
        let (k, m, q) = build_k_vector(&pairs);

        Ok(Self {
            pairs,
            k,
            m,
            q,
            star_count,
            max_angle,
        })
    }

    /// Pairs whose separation lies within `eps` of `angle`, as a slice.
    ///
    /// O(1): the k-vector converts the angle bounds straight into array bounds,
    /// which are then nudged onto the exact range. Both inputs are radians, and
    /// the result is ordered by angle like the table itself.
    pub fn query(&self, angle: f64, eps: f64) -> &[Pair] {
        let (low, high) = (angle - eps.abs(), angle + eps.abs());
        let n = self.pairs.len();
        if n == 0 || high < self.first_angle() || low > self.last_angle() {
            return &[];
        }

        // The index line is increasing, so the bracketing entries of `k` give a
        // range that is guaranteed to contain the answer.
        let (mut start, mut end) = if self.m > 0.0 {
            let lo_index = ((low - self.q) / self.m).floor();
            let hi_index = ((high - self.q) / self.m).ceil();
            (
                self.k[clamp_index(lo_index, n)] as usize,
                self.k[clamp_index(hi_index, n)] as usize,
            )
        } else {
            // Every separation is identical; the whole table is the candidate.
            (0, n)
        };

        // Nudge to the exact bounds. The index line models a uniform density,
        // so on this distribution each loop runs a couple of times at most.
        while start > 0 && f64::from(self.pairs[start - 1].angle) >= low {
            start -= 1;
        }
        while start < n && f64::from(self.pairs[start].angle) < low {
            start += 1;
        }
        end = end.max(start);
        while end < n && f64::from(self.pairs[end].angle) <= high {
            end += 1;
        }
        while end > start && f64::from(self.pairs[end - 1].angle) > high {
            end -= 1;
        }
        &self.pairs[start..end]
    }

    /// Brute-force equivalent of [`PairDb::query`], as a test oracle.
    ///
    /// Scans the whole table, so it is O(n) and exists only to prove the
    /// k-vector path returns exactly the same pairs.
    pub fn query_slow(&self, angle: f64, eps: f64) -> Vec<Pair> {
        let (low, high) = (angle - eps.abs(), angle + eps.abs());
        self.pairs
            .iter()
            .filter(|pair| {
                let a = f64::from(pair.angle);
                a >= low && a <= high
            })
            .copied()
            .collect()
    }

    /// Every pair, ascending by separation.
    pub fn pairs(&self) -> &[Pair] {
        &self.pairs
    }

    /// How many catalogue stars the database covers. Pair indices are below it.
    pub fn star_count(&self) -> usize {
        self.star_count
    }

    /// Largest separation stored, radians.
    pub fn max_angle(&self) -> f64 {
        self.max_angle
    }

    /// The k-vector index, for tests and diagnostics.
    pub fn k_vector(&self) -> &[u32] {
        &self.k
    }

    /// Bytes of heap the table occupies, for start-up reporting.
    pub fn heap_bytes(&self) -> usize {
        self.pairs.len() * size_of::<Pair>() + self.k.len() * size_of::<u32>()
    }

    fn first_angle(&self) -> f64 {
        self.pairs.first().map_or(0.0, |p| f64::from(p.angle))
    }

    fn last_angle(&self) -> f64 {
        self.pairs.last().map_or(0.0, |p| f64::from(p.angle))
    }
}

/// Clamps a line index to a valid `k` subscript.
fn clamp_index(value: f64, n: usize) -> usize {
    if value <= 0.0 {
        0
    } else if value >= (n - 1) as f64 {
        n - 1
    } else {
        value as usize
    }
}

/// Collects every pair of the first `star_count` catalogue stars separated by
/// no more than `max_angle`, sorted by separation.
fn collect_pairs(catalog: &Catalog, star_count: usize, max_angle: f64) -> Vec<Pair> {
    // Walking the stars in order of z lets the inner loop stop early: for unit
    // vectors |sin d_b - sin d_a| <= |d_b - d_a| <= angle, so once the z gap
    // exceeds the limit no later star can be within it. That alone removes
    // more than half the candidate pairs.
    let mut order: Vec<u16> = (0..star_count as u16).collect();
    order.sort_unstable_by(|&a, &b| {
        catalog.stars[usize::from(a)].unit[2].total_cmp(&catalog.stars[usize::from(b)].unit[2])
    });
    let directions: Vec<Vec3> = order
        .iter()
        .map(|&index| catalog.stars[usize::from(index)].direction())
        .collect();

    // Comparing dot products rejects a candidate without an atan2; the real
    // angle is only computed for pairs that survive.
    let cos_min = max_angle.cos();

    // Pre-sizing matters: the table runs to the better part of a million
    // entries, and growing into it would copy megabytes several times over.
    // For a uniform sky the count is n(n-1)/2 * (1 - cos)/2; the real sky is
    // clustered, so allow a wide margin over that.
    let n = star_count as f64;
    let estimate = n * (n - 1.0) / 2.0 * (1.0 - cos_min) / 2.0 * 1.3;
    let mut pairs: Vec<Pair> = Vec::with_capacity(estimate as usize);

    for (a, first) in directions.iter().enumerate() {
        for (offset, second) in directions[a + 1..].iter().enumerate() {
            if second.z - first.z > max_angle {
                break;
            }
            if first.dot(second) < cos_min {
                continue;
            }
            // Store catalogue indices, not positions in the z-ordered scan.
            let (i, j) = {
                let (p, q) = (order[a], order[a + 1 + offset]);
                if p < q { (p, q) } else { (q, p) }
            };
            pairs.push(Pair {
                angle: angle_between(first, second) as f32,
                i,
                j,
            });
        }
    }

    pairs.sort_unstable_by(|a, b| {
        a.angle
            .total_cmp(&b.angle)
            .then(a.i.cmp(&b.i))
            .then(a.j.cmp(&b.j))
    });
    pairs
}

/// Builds the k-vector index over pairs already sorted by angle.
///
/// The index line runs from just below the smallest separation to just above
/// the largest, and `k[i]` counts the pairs below the line at `i`. Because the
/// line only increases, the whole vector is filled in one pass.
fn build_k_vector(pairs: &[Pair]) -> (Vec<u32>, f64, f64) {
    let n = pairs.len();
    if n < 2 {
        return (vec![0; n], 0.0, 0.0);
    }

    let first = f64::from(pairs[0].angle);
    let last = f64::from(pairs[n - 1].angle);
    let span = last - first;
    if span <= 0.0 {
        // Every separation is identical, so no line can separate them.
        return (vec![0; n], 0.0, first);
    }

    // The offset guarantees the line starts strictly below the first angle and
    // ends strictly above the last, so k[0] is 0 and k[n-1] is n.
    let delta = span * 1e-9;
    let m = (span + 2.0 * delta) / (n - 1) as f64;
    let q = first - delta;

    let mut k = Vec::with_capacity(n);
    let mut below = 0usize;
    for index in 0..n {
        let level = q + m * index as f64;
        while below < n && f64::from(pairs[below].angle) < level {
            below += 1;
        }
        k.push(below as u32);
    }
    (k, m, q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::Star;
    use crate::simulate::Preset;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use std::collections::HashSet;
    use std::f64::consts::PI;

    /// One arcsecond in radians.
    const ARCSEC: f64 = PI / (180.0 * 3600.0);

    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

    fn catalog() -> Catalog {
        Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode")
    }

    /// The database the solver actually uses: everything to `db_mag_cut`.
    fn full(catalog: &Catalog) -> PairDb {
        PairDb::build(catalog, &SimConfig::preset(Preset::Nominal)).expect("build")
    }

    /// A much smaller database, for the structural tests. The full one holds
    /// 850k pairs, and building it repeatedly would dominate a debug test run.
    fn small(catalog: &Catalog) -> (PairDb, SimConfig) {
        let cfg = SimConfig {
            db_mag_cut: 4.5,
            ..SimConfig::preset(Preset::Nominal)
        };
        let db = PairDb::build(catalog, &cfg).expect("build");
        (db, cfg)
    }

    fn rng(seed: u64) -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(seed)
    }

    // --- structure ---

    #[test]
    fn covers_exactly_the_magnitude_cut() {
        let catalog = catalog();
        let (db, cfg) = small(&catalog);

        assert!(db.star_count() > 0);
        // The cut is a prefix of the catalogue, so the star just past the end
        // must be fainter than the cut and the last one inside it brighter.
        assert!(catalog.stars[db.star_count() - 1].mag <= cfg.db_mag_cut);
        assert!(catalog.stars[db.star_count()].mag > cfg.db_mag_cut);
    }

    #[test]
    fn max_angle_is_the_field_diagonal_plus_a_margin() {
        let catalog = catalog();
        let (db, cfg) = small(&catalog);
        let diagonal = 2.0 * cfg.nominal_camera().half_diagonal_fov();
        assert!((db.max_angle() - diagonal - 1.0f64.to_radians()).abs() < 1e-12);
        // Roughly 28 degrees for the default camera.
        assert!((28.0..30.0).contains(&db.max_angle().to_degrees()));
    }

    #[test]
    fn pairs_are_sorted_and_within_range() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        assert!(!db.pairs().is_empty());

        for window in db.pairs().windows(2) {
            assert!(
                window[0].angle <= window[1].angle,
                "table is not sorted by angle"
            );
        }
        for pair in db.pairs() {
            assert!(f64::from(pair.angle) <= db.max_angle());
            assert!(pair.angle >= 0.0);
            assert!(pair.i < pair.j, "indices must be ordered within a pair");
            assert!(usize::from(pair.j) < db.star_count());
        }
    }

    /// docs/SPEC.md pitfall 6: the table must be keyed by catalogue index, not by
    /// any per-frame ordering. Recomputing each separation straight from the
    /// catalogue is what proves it.
    #[test]
    fn indices_are_catalogue_indices() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let mut r = rng(0xC0FFEE);

        for _ in 0..2000 {
            let pair = db.pairs()[r.random_range(0..db.pairs().len())];
            let (i, j) = pair.indices();
            let expected =
                angle_between(&catalog.stars[i].direction(), &catalog.stars[j].direction());
            // Stored as f32, so compare at f32 resolution.
            assert!(
                (f64::from(pair.angle) - expected).abs() < 1e-6,
                "pair {i}-{j} stored {} but the catalogue says {expected}",
                pair.angle
            );
        }
    }

    #[test]
    fn pairs_are_unique() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let mut seen = HashSet::with_capacity(db.pairs().len());
        for pair in db.pairs() {
            assert!(
                seen.insert((pair.i, pair.j)),
                "pair {}-{} appears twice",
                pair.i,
                pair.j
            );
        }
    }

    /// Nothing inside the field may be missing, or identification would fail
    /// on a frame that happens to contain that pair.
    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "quadratic in the star count, run with --release"
    )]
    fn every_pair_within_range_is_present() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let stored: HashSet<(u16, u16)> = db.pairs().iter().map(|p| (p.i, p.j)).collect();

        let mut expected = 0usize;
        for i in 0..db.star_count() {
            for j in (i + 1)..db.star_count() {
                let angle =
                    angle_between(&catalog.stars[i].direction(), &catalog.stars[j].direction());
                if angle <= db.max_angle() {
                    expected += 1;
                    assert!(
                        stored.contains(&(i as u16, j as u16)),
                        "pair {i}-{j} at {angle} rad is missing from the table"
                    );
                }
            }
        }
        assert_eq!(expected, db.pairs().len(), "table holds extra pairs");
    }

    // --- the k-vector ---

    #[test]
    fn k_vector_is_monotonic_and_spans_the_table() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let k = db.k_vector();
        let n = db.pairs().len();

        assert_eq!(k.len(), n);
        assert_eq!(k[0], 0, "the line must start below the first separation");
        assert_eq!(
            k[n - 1] as usize,
            n,
            "the line must end above the last separation"
        );
        for window in k.windows(2) {
            assert!(window[0] <= window[1], "k-vector must not decrease");
        }
    }

    /// The point of the index is that it brackets tightly. On this
    /// distribution the density varies only about twofold across the range, so
    /// a query should not be handed many more candidates than it keeps.
    #[test]
    fn k_vector_brackets_tightly() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let mut r = rng(4);
        let eps = 20.0 * ARCSEC;

        let mut worst_ratio = 0.0f64;
        for _ in 0..500 {
            let angle = 0.02 + r.random::<f64>() * (db.max_angle() - 0.04);
            let fast = db.query(angle, eps).len();
            let slow = db.query_slow(angle, eps).len();
            assert_eq!(fast, slow);
            if slow > 0 {
                worst_ratio = worst_ratio.max(fast as f64 / slow as f64);
            }
        }
        // `query` returns the exact set, so this is 1 by construction; the
        // assertion guards against a future change that returns candidates.
        assert!((worst_ratio - 1.0).abs() < 1e-12);
    }

    // --- acceptance: query against the brute-force oracle ---

    /// Phase 4 acceptance: `query` equals `query_slow` on 10k random queries.
    ///
    /// Run against the full database, the one the solver uses, and across a
    /// spread of tolerances including degenerate ones. `query_slow` scans all
    /// 850k pairs per call, so this is a release-only test.
    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "10k linear scans of 850k pairs, run with --release"
    )]
    fn query_matches_brute_force_on_10k_random_queries() {
        let catalog = catalog();
        let db = full(&catalog);
        let mut r = rng(0x10_000);

        let mut nonempty = 0usize;
        for attempt in 0..10_000 {
            // Mostly inside the table, sometimes outside it, with tolerances
            // from zero to far wider than any real search.
            let angle = match attempt % 10 {
                0 => -0.1 + r.random::<f64>() * 0.1,
                1 => db.max_angle() + r.random::<f64>() * 0.1,
                _ => r.random::<f64>() * db.max_angle(),
            };
            let eps = match attempt % 7 {
                0 => 0.0,
                1 => 1.0 * ARCSEC,
                2 => 0.05,
                _ => r.random::<f64>() * 60.0 * ARCSEC,
            };

            let fast = db.query(angle, eps);
            let slow = db.query_slow(angle, eps);
            assert_eq!(
                fast,
                slow.as_slice(),
                "mismatch at angle {angle} eps {eps}: {} vs {} pairs",
                fast.len(),
                slow.len()
            );
            if !fast.is_empty() {
                nonempty += 1;
            }
        }
        // A test where everything came back empty would prove nothing. Two in
        // ten queries deliberately sit outside the table and one in seven asks
        // for a zero tolerance, so around 6900 hits is the expected yield.
        assert!(
            nonempty > 6000,
            "only {nonempty} of 10000 queries returned anything"
        );
    }

    #[test]
    fn query_results_lie_inside_the_tolerance() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let mut r = rng(11);
        for _ in 0..500 {
            let angle = r.random::<f64>() * db.max_angle();
            let eps = r.random::<f64>() * 30.0 * ARCSEC;
            for pair in db.query(angle, eps) {
                assert!((f64::from(pair.angle) - angle).abs() <= eps);
            }
        }
    }

    #[test]
    fn query_handles_the_edges() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let pairs = db.pairs();
        let (first, last) = (
            f64::from(pairs[0].angle),
            f64::from(pairs[pairs.len() - 1].angle),
        );

        // Well outside the table in either direction.
        assert!(db.query(-1.0, 0.1).is_empty());
        assert!(db.query(db.max_angle() + 1.0, 0.1).is_empty());

        // Exactly on the smallest and largest separations, with no tolerance.
        assert!(!db.query(first, 0.0).is_empty(), "first pair unreachable");
        assert!(!db.query(last, 0.0).is_empty(), "last pair unreachable");
        assert_eq!(db.query(first, 0.0), db.query_slow(first, 0.0).as_slice());
        assert_eq!(db.query(last, 0.0), db.query_slow(last, 0.0).as_slice());

        // A tolerance wide enough to take everything.
        assert_eq!(db.query(last / 2.0, last).len(), pairs.len());

        // A negative tolerance is treated as its magnitude.
        assert_eq!(db.query(first, -0.01), db.query(first, 0.01));
    }

    // --- helpers on Pair ---

    #[test]
    fn pair_accessors_agree_with_the_fields() {
        let pair = Pair {
            angle: 0.25,
            i: 7,
            j: 19,
        };
        assert_eq!(pair.indices(), (7, 19));
        assert!(pair.contains(7) && pair.contains(19) && !pair.contains(8));
        assert_eq!(pair.other(7), Some(19));
        assert_eq!(pair.other(19), Some(7));
        assert_eq!(pair.other(8), None);
    }

    // --- determinism and limits ---

    #[test]
    fn two_builds_are_identical() {
        let catalog = catalog();
        let (first, cfg) = small(&catalog);
        let second = PairDb::build(&catalog, &cfg).expect("build");
        assert_eq!(first.pairs(), second.pairs());
        assert_eq!(first.k_vector(), second.k_vector());
        assert_eq!(first.star_count(), second.star_count());
    }

    #[test]
    fn a_tighter_cut_gives_a_subset() {
        let catalog = catalog();
        let bright_cfg = SimConfig {
            db_mag_cut: 3.5,
            ..SimConfig::preset(Preset::Nominal)
        };
        let bright = PairDb::build(&catalog, &bright_cfg).expect("build");
        let (wider, _) = small(&catalog);

        assert!(bright.star_count() < wider.star_count());
        assert!(bright.pairs().len() < wider.pairs().len());

        let wide: HashSet<(u16, u16)> = wider.pairs().iter().map(|p| (p.i, p.j)).collect();
        for pair in bright.pairs() {
            assert!(
                wide.contains(&(pair.i, pair.j)),
                "pair {}-{} is absent from the wider database",
                pair.i,
                pair.j
            );
        }
    }

    /// Pair indices are `u16`, so a cut admitting more stars than that has to
    /// be refused rather than silently wrapping.
    #[test]
    fn refuses_more_stars_than_a_u16_index_can_address() {
        // The guard runs before any pair work, so this stays cheap.
        let oversized = Catalog {
            stars: (0..=MAX_DB_STARS)
                .map(|n| Star {
                    id: n as u32 + 1,
                    unit: [1.0, 0.0, 0.0],
                    mag: 1.0,
                    name: None,
                })
                .collect(),
        };
        let result = PairDb::build(&oversized, &SimConfig::preset(Preset::Nominal));
        assert!(matches!(
            result,
            Err(Error::TooManyDatabaseStars { count }) if count == MAX_DB_STARS + 1
        ));
    }

    #[test]
    fn heap_use_is_reported() {
        let catalog = catalog();
        let (db, _) = small(&catalog);
        let expected = db.pairs().len() * 8 + db.k_vector().len() * 4;
        assert_eq!(db.heap_bytes(), expected);
        assert_eq!(size_of::<Pair>(), 8, "Pair must stay 8 bytes");
    }
}
