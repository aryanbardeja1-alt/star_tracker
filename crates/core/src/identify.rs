//! Star identification by the Pyramid algorithm (Mortari et al. 2004).
//!
//! Given unit vectors of the observed stars in the camera frame, the pair
//! database and the catalogue, this works out which catalogue star each
//! observation is, with no prior attitude knowledge.
//!
//! The shape of it, following docs/SPEC.md:
//!
//! 1. For a triple of observed stars, three pair-database lookups give
//!    candidate catalogue pairs; a triangle candidate is an assignment whose
//!    catalogue indices close the loop.
//! 2. A fourth observed star confirms the triangle, checking all three of its
//!    angles against it.
//! 3. Triples are visited in Mortari's order, which varies the index *gaps* in
//!    the outer loops, so a single false star cannot sit in every early try.
//! 4. A confirmed pyramid gives a provisional attitude.
//! 5. That attitude identifies the remaining observations by nearest
//!    catalogue neighbour.
//!
//! All angles are radians; observed vectors are camera frame, catalogue vectors
//! inertial, related by `b = A r`.

use crate::Clock;
use crate::attitude::solve_wahba;
use crate::camera::CameraModel;
use crate::catalog::Catalog;
use crate::math::{Mat3, Vec3, angle_between};
use crate::pairdb::PairDb;
use crate::simulate::SimConfig;

/// A pyramid needs four stars.
const PYRAMID_STARS: usize = 4;

/// Triangle assignments examined per triple before moving on.
///
/// A triple normally yields zero or one, so this only bites on a degenerate
/// field where the tolerance admits a crowd of look-alikes, and there it stops
/// one triple from eating the whole iteration budget.
const MAX_CANDIDATES_PER_TRIPLE: u32 = 32;

/// One identified star.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StarMatch {
    /// Index into the observed vectors handed to [`identify`].
    pub observed: u16,
    /// Index into `Catalog::stars`.
    pub catalog_index: u16,
    /// Catalogue (HIP) number of that star -- the same star as
    /// `catalog_index`, in the form the ground-truth check compares against.
    pub id: u32,
    /// Whether this star came from the confirmed pyramid rather than the
    /// attitude-assisted sweep. The self-check requires at least four of these.
    pub from_pyramid: bool,
}

/// Counters and timing from one identification attempt.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Diagnostics {
    /// Observed-star triples examined.
    pub tries: u32,
    /// Pair-database queries issued.
    pub pair_queries: u32,
    /// Pairs returned by those queries, summed.
    pub pairs_examined: u64,
    /// Triangle assignments whose catalogue indices closed the loop.
    pub triangle_candidates: u32,
    /// Triangles a fourth star confirmed.
    pub confirmed_pyramids: u32,
    /// Confirmed pyramids thrown away because they explained too little of the
    /// frame. A coincidence can confirm four stars; it cannot predict ten.
    pub discarded_pyramids: u32,
    /// Stars added by the attitude-assisted sweep.
    pub swept_matches: u32,
    /// Whether the iteration cap was reached.
    pub hit_try_cap: bool,
    /// Wall time for the attempt, nanoseconds, from the host's clock.
    pub elapsed_ns: u64,
}

/// What an identification attempt found.
#[derive(Clone, Debug)]
pub struct Solution {
    /// Every identified star, the four pyramid members first.
    pub matches: Vec<StarMatch>,
    /// Attitude from the matched stars, with `b = A r`.
    pub attitude: Mat3,
}

impl Solution {
    /// How many matches came from the confirmed pyramid.
    pub fn pyramid_count(&self) -> usize {
        self.matches.iter().filter(|m| m.from_pyramid).count()
    }
}

/// The outcome of one identification attempt.
///
/// The counters come back whether or not anything was found, because a trial
/// that identified nothing is exactly the one whose counters are interesting.
#[derive(Clone, Debug)]
pub struct Identification {
    /// The solution, when a pyramid was confirmed.
    pub solution: Option<Solution>,
    /// Counters and timing.
    pub diagnostics: Diagnostics,
}

/// Angular tolerances derived from the configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    /// Constant part of the tolerance on an angle *between two* observed
    /// stars, radians: `k_sigma * sqrt(2) * sigma_px / f + calib_margin`.
    pub pair_base: f64,
    /// Part of the pair tolerance proportional to the separation itself.
    ///
    /// A focal length off by `d` stretches every observed angle by about `d`,
    /// so a pair's angular error grows with the angle. Without this term the
    /// tolerance covers only close pairs and the pyramid, which needs a spread
    /// triangle, never closes.
    pub pair_per_rad: f64,
    /// Tolerance on a *single* direction, radians. The same expression without
    /// the `sqrt(2)`, since only one centroid contributes.
    pub vector: f64,
}

impl Tolerance {
    /// Derives both tolerances from the configuration, per docs/SPEC.md:
    /// `eps = k_sigma * sigma_angle + calib_margin` with
    /// `sigma_angle = sqrt(2) * sigma_centroid_px / f`.
    pub fn from_config(cfg: &SimConfig) -> Self {
        let per_vector = cfg.centroid_sigma_px / cfg.nominal_camera().f;
        let margin = cfg.calib_margin_arcsec * std::f64::consts::PI / (180.0 * 3600.0);
        Self {
            pair_base: cfg.id_k_sigma * std::f64::consts::SQRT_2 * per_vector + margin,
            pair_per_rad: cfg.id_focal_tolerance_frac,
            vector: cfg.id_k_sigma * per_vector + margin,
        }
    }

    /// Tolerance for an observed separation of `theta` radians.
    pub fn pair(&self, theta: f64) -> f64 {
        self.pair_base + self.pair_per_rad * theta
    }
}

/// Scratch buffers reused across frames.
#[derive(Clone, Debug, Default)]
pub struct Workspace {
    /// Angles between observed stars, row-major and square.
    angles: Vec<f64>,
    /// Stamped partner maps for the two secondary pair lookups. `stamp` says
    /// which generation an entry belongs to, so only touched entries are
    /// cleared.
    stamp_ik: Vec<u32>,
    partners_ik: Vec<Vec<u16>>,
    stamp_jk: Vec<u32>,
    partners_jk: Vec<Vec<u16>>,
    generation: u32,
    /// Catalogue stars inside the provisional field, with their predicted
    /// camera-frame directions.
    in_field: Vec<(u16, Vec3)>,
    /// Result buffer.
    matches: Vec<StarMatch>,
}

impl Workspace {
    /// An empty workspace; buffers size themselves on first use.
    pub fn new() -> Self {
        Self::default()
    }

    fn prepare(&mut self, observed: usize, db_stars: usize) {
        self.angles.clear();
        self.angles.resize(observed * observed, 0.0);
        self.stamp_ik.resize(db_stars, 0);
        self.stamp_jk.resize(db_stars, 0);
        self.partners_ik.resize(db_stars, Vec::new());
        self.partners_jk.resize(db_stars, Vec::new());
        self.matches.clear();
        self.in_field.clear();
    }
}

/// Identifies the observed stars.
///
/// `observed` are camera-frame unit vectors, conventionally ordered brightest
/// first so the Mortari sweep reaches the most reliable stars soonest.
pub fn identify(
    observed: &[Vec3],
    catalog: &Catalog,
    db: &PairDb,
    cfg: &SimConfig,
    clock: &dyn Clock,
    ws: &mut Workspace,
) -> Identification {
    let started = clock.now_ns();
    let mut diagnostics = Diagnostics::default();
    let solution = solve(observed, catalog, db, cfg, ws, &mut diagnostics);
    diagnostics.elapsed_ns = clock.now_ns().saturating_sub(started);
    Identification {
        solution,
        diagnostics,
    }
}

/// The identification proper, with timing left to [`identify`].
fn solve(
    observed: &[Vec3],
    catalog: &Catalog,
    db: &PairDb,
    cfg: &SimConfig,
    ws: &mut Workspace,
    diagnostics: &mut Diagnostics,
) -> Option<Solution> {
    let n = observed.len();
    if n < PYRAMID_STARS {
        return None;
    }

    let base = Tolerance::from_config(cfg);
    ws.prepare(n, db.star_count());

    // Every observed angle, once.
    for a in 0..n {
        for b in (a + 1)..n {
            let angle = angle_between(&observed[a], &observed[b]);
            ws.angles[a * n + b] = angle;
            ws.angles[b * n + a] = angle;
        }
    }

    // Two passes over the tolerance, tight before wide.
    //
    // The wide pass exists for a camera whose focal length is off, where a
    // pair's angular error grows with the separation and a constant tolerance
    // finds nothing. But it is not free: a wider window admits wrong pairs, and
    // with them triangles that a fourth star can still confirm by coincidence.
    // Charging every frame for that would make the common case both slower and
    // less reliable, so a frame only pays for the wide pass if the tight one
    // cannot solve it.
    let passes: &[f64] = if base.pair_per_rad > 0.0 {
        &[0.0, base.pair_per_rad]
    } else {
        &[0.0]
    };
    // One budget across both passes: `id_max_tries` caps the triples an
    // identification will try, and splitting it per pass would quietly double
    // that. A tight pass over the default fifteen centroids spends 455 of the
    // thousand, so the wide pass still has room for a full second walk.
    let mut budget = cfg.id_max_tries;
    for &per_rad in passes {
        let tolerance = Tolerance {
            pair_per_rad: per_rad,
            ..base
        };

        // Walk the triples once, across however many pyramids get confirmed: a
        // confirmed pyramid that turns out not to describe the frame must not end
        // the search, or one coincidence costs the whole solution.
        let mut triples = MortariTriples::new(n);
        while let Some(pyramid) = find_pyramid(
            &mut triples,
            &mut budget,
            observed,
            catalog,
            db,
            &tolerance,
            ws,
            diagnostics,
        ) {
            // A provisional attitude from the four confirmed stars, then the sweep.
            let Some(provisional) = attitude_from(&pyramid, observed, catalog) else {
                continue;
            };
            ws.matches.clear();
            for (observed_index, catalog_index) in pyramid {
                ws.matches.push(StarMatch {
                    observed: observed_index,
                    catalog_index,
                    id: catalog.stars[usize::from(catalog_index)].id,
                    from_pyramid: true,
                });
            }

            // Take the workspace's match list out so the sweep can borrow both it
            // and the workspace's own buffers.
            let mut matches = std::mem::take(&mut ws.matches);
            let added = extend_matches(
                observed,
                &provisional,
                &cfg.nominal_camera(),
                catalog,
                db.star_count(),
                tolerance.vector,
                &mut matches,
                ws,
            );

            // Does this pyramid explain the rest of the frame? A right one does:
            // its attitude puts catalogue stars on the other centroids, usually ten
            // or more of them, while a wrong one predicts nothing at all -- every
            // false pyramid measured swept up exactly zero. So the bar is one star
            // beyond the four, which is the weakest test that separates them, and
            // it only applies when the frame has another star to offer: a frame of
            // four has nothing left to predict and must not be failed for it.
            if added == 0 && n > PYRAMID_STARS {
                diagnostics.discarded_pyramids += 1;
                matches.clear();
                ws.matches = matches;
                continue;
            }

            diagnostics.swept_matches += added as u32;
            let attitude = attitude_over(&matches, observed, catalog).unwrap_or(provisional);

            return Some(Solution { matches, attitude });
        }
    }
    None
}

/// Pair examinations one identification may spend before giving up.
///
/// `id_max_tries` caps the *triples* tried, which bounds the work only while
/// every triple costs about the same. The wide tolerance pass breaks that: its
/// windows are several times broader, so a triple there examines several times
/// the pairs, and a frame that cannot be solved spends the whole budget at the
/// higher rate. Measured over 200 trials a solved frame examines at most 15.2
/// million pairs on `hard` and 103 thousand on `nominal`, while a hopeless
/// `brutal` frame runs to 387 million. This bounds the damage at roughly twice
/// what the hardest solvable frame needs. It is deterministic, so it does not
/// disturb native/WASM parity.
const MAX_PAIRS_EXAMINED: u64 = 32_000_000;

/// The four `(observed index, catalogue index)` pairs of a confirmed pyramid.
type Pyramid = [(u16, u16); PYRAMID_STARS];

/// Searches observed triples in Mortari's order for a confirmable pyramid.
#[allow(clippy::too_many_arguments)]
fn find_pyramid(
    triples: &mut MortariTriples,
    budget: &mut u32,
    observed: &[Vec3],
    catalog: &Catalog,
    db: &PairDb,
    tolerance: &Tolerance,
    ws: &mut Workspace,
    diagnostics: &mut Diagnostics,
) -> Option<Pyramid> {
    for (i, j, k) in triples.by_ref() {
        if *budget == 0 || diagnostics.pairs_examined >= MAX_PAIRS_EXAMINED {
            diagnostics.hit_try_cap = true;
            return None;
        }
        *budget -= 1;
        diagnostics.tries += 1;

        if let Some(pyramid) =
            try_triple(i, j, k, observed, catalog, db, tolerance, ws, diagnostics)
        {
            return Some(pyramid);
        }
    }
    None
}

/// Mortari's ordering over triples of observed stars.
///
/// The outer loops walk the index *gaps* rather than the indices, so
/// consecutive triples do not share all their stars: after `(0,1,2)` come
/// `(1,2,3)` and `(2,3,4)` before anything returns to index 0. A plain nested
/// `i<j<k` sweep would instead put index 0 in its whole first run, letting one
/// false star there poison every early attempt.
///
/// Yields each triple `i < j < k` exactly once, so `C(n, 3)` in total.
struct MortariTriples {
    n: usize,
    gap_j: usize,
    gap_k: usize,
    i: usize,
}

impl MortariTriples {
    fn new(n: usize) -> Self {
        Self {
            n,
            gap_j: 1,
            gap_k: 1,
            i: 0,
        }
    }
}

impl Iterator for MortariTriples {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while self.gap_j + self.gap_k < self.n {
            if self.i + self.gap_j + self.gap_k < self.n {
                let i = self.i;
                self.i += 1;
                let j = i + self.gap_j;
                return Some((i, j, j + self.gap_k));
            }
            // This pair of gaps is exhausted; widen the inner one, and the
            // outer one once the inner can grow no further.
            self.i = 0;
            self.gap_k += 1;
            if self.gap_j + self.gap_k >= self.n {
                self.gap_j += 1;
                self.gap_k = 1;
            }
        }
        None
    }
}

/// Tries one observed triple: close a triangle, then confirm it.
#[allow(clippy::too_many_arguments)]
fn try_triple(
    i: usize,
    j: usize,
    k: usize,
    observed: &[Vec3],
    catalog: &Catalog,
    db: &PairDb,
    tolerance: &Tolerance,
    ws: &mut Workspace,
    diagnostics: &mut Diagnostics,
) -> Option<Pyramid> {
    let n = observed.len();
    let angle = |a: usize, b: usize| ws.angles[a * n + b];
    let (theta_ij, theta_ik, theta_jk) = (angle(i, j), angle(i, k), angle(j, k));

    let pairs_ij = db.query(theta_ij, tolerance.pair(theta_ij));
    let pairs_ik = db.query(theta_ik, tolerance.pair(theta_ik));
    let pairs_jk = db.query(theta_jk, tolerance.pair(theta_jk));
    diagnostics.pair_queries += 3;
    diagnostics.pairs_examined += (pairs_ij.len() + pairs_ik.len() + pairs_jk.len()) as u64;
    if pairs_ij.is_empty() || pairs_ik.is_empty() || pairs_jk.is_empty() {
        return None;
    }

    // Index the two secondary lists by star, so closing the loop is a lookup
    // rather than a scan. Both maps are generation-stamped, so only the
    // entries this triple touches need clearing.
    ws.generation = ws.generation.wrapping_add(1);
    let generation = ws.generation;
    for pair in pairs_ik {
        push_partner(
            &mut ws.stamp_ik,
            &mut ws.partners_ik,
            generation,
            pair.i,
            pair.j,
        );
        push_partner(
            &mut ws.stamp_ik,
            &mut ws.partners_ik,
            generation,
            pair.j,
            pair.i,
        );
    }
    for pair in pairs_jk {
        push_partner(
            &mut ws.stamp_jk,
            &mut ws.partners_jk,
            generation,
            pair.i,
            pair.j,
        );
        push_partner(
            &mut ws.stamp_jk,
            &mut ws.partners_jk,
            generation,
            pair.j,
            pair.i,
        );
    }

    // A triangle candidate is an assignment (I, J, K) present in all three
    // lists. Both orientations of each i-j pair have to be tried, because the
    // triangle is not symmetric: theta_ik and theta_jk differ.
    let mut examined = 0u32;
    let mut confirmed: Option<Pyramid> = None;
    for pair in pairs_ij {
        for (star_i, star_j) in [(pair.i, pair.j), (pair.j, pair.i)] {
            let candidates_k = partners_of(&ws.stamp_ik, &ws.partners_ik, generation, star_i);
            if candidates_k.is_empty() {
                continue;
            }
            let from_j = partners_of(&ws.stamp_jk, &ws.partners_jk, generation, star_j);
            if from_j.is_empty() {
                continue;
            }

            for &star_k in candidates_k {
                if !from_j.contains(&star_k) {
                    continue;
                }
                diagnostics.triangle_candidates += 1;
                examined += 1;
                if examined > MAX_CANDIDATES_PER_TRIPLE {
                    return None;
                }

                let triangle = [(i, star_i), (j, star_j), (k, star_k)];
                if let Some(fourth) = confirm_with_fourth(
                    &triangle,
                    observed,
                    catalog,
                    db,
                    tolerance,
                    ws,
                    diagnostics,
                ) {
                    diagnostics.confirmed_pyramids += 1;
                    let pyramid = [
                        (i as u16, star_i),
                        (j as u16, star_j),
                        (k as u16, star_k),
                        fourth,
                    ];
                    if confirmed.is_some_and(|existing| existing != pyramid) {
                        // Two different readings of the same triple both
                        // confirm: the field is ambiguous at this tolerance, so
                        // claiming either would be a guess.
                        return None;
                    }
                    confirmed = Some(pyramid);
                }
            }
        }
    }
    confirmed
}

/// Looks for a fourth observed star that confirms a triangle.
///
/// All three of its angles to the triangle must agree, which is what makes a
/// pyramid far harder to satisfy by accident than a triangle alone.
#[allow(clippy::too_many_arguments)]
fn confirm_with_fourth(
    triangle: &[(usize, u16); 3],
    observed: &[Vec3],
    catalog: &Catalog,
    db: &PairDb,
    tolerance: &Tolerance,
    ws: &Workspace,
    diagnostics: &mut Diagnostics,
) -> Option<(u16, u16)> {
    let n = observed.len();
    let angle = |a: usize, b: usize| ws.angles[a * n + b];
    let [(i, star_i), (j, star_j), (k, star_k)] = *triangle;

    for r in 0..n {
        if r == i || r == j || r == k {
            continue;
        }
        // Candidates for the fourth catalogue star come from its separation to
        // the triangle's first star; the other two angles are then checked
        // straight against the catalogue.
        let candidates = db.query(angle(r, i), tolerance.pair(angle(r, i)));
        diagnostics.pair_queries += 1;
        diagnostics.pairs_examined += candidates.len() as u64;

        for pair in candidates {
            let Some(star_r) = pair.other(star_i) else {
                continue;
            };
            if star_r == star_j || star_r == star_k {
                continue;
            }
            let direction = catalog.stars[usize::from(star_r)].direction();
            let to_j = angle_between(&direction, &catalog.stars[usize::from(star_j)].direction());
            let to_k = angle_between(&direction, &catalog.stars[usize::from(star_k)].direction());
            if (to_j - angle(r, j)).abs() <= tolerance.pair(angle(r, j))
                && (to_k - angle(r, k)).abs() <= tolerance.pair(angle(r, k))
            {
                return Some((r as u16, star_r));
            }
        }
    }
    None
}

/// Adds `partner` to `star`'s list, resetting the list if it is stale.
fn push_partner(
    stamp: &mut [u32],
    partners: &mut [Vec<u16>],
    generation: u32,
    star: u16,
    partner: u16,
) {
    let index = usize::from(star);
    if stamp[index] != generation {
        stamp[index] = generation;
        partners[index].clear();
    }
    partners[index].push(partner);
}

/// This generation's partners of `star`, or empty if it has none.
fn partners_of<'a>(
    stamp: &[u32],
    partners: &'a [Vec<u16>],
    generation: u32,
    star: u16,
) -> &'a [u16] {
    let index = usize::from(star);
    if stamp[index] == generation {
        &partners[index]
    } else {
        &[]
    }
}

/// Attitude from a confirmed pyramid, every star weighted equally.
fn attitude_from(pyramid: &Pyramid, observed: &[Vec3], catalog: &Catalog) -> Option<Mat3> {
    let mut pairs = [(Vec3::zeros(), Vec3::zeros(), 0.0); PYRAMID_STARS];
    for (slot, &(observed_index, catalog_index)) in pairs.iter_mut().zip(pyramid.iter()) {
        *slot = (
            observed[usize::from(observed_index)],
            catalog.stars[usize::from(catalog_index)].direction(),
            1.0,
        );
    }
    solve_wahba(&pairs)
}

/// Extends `matches` with the observations an attitude can account for.
///
/// Given an attitude and the camera it was solved with, each unmatched
/// observation is paired with the nearest catalogue star inside the field, so
/// long as that star is within `tolerance` and no runner-up is -- an ambiguous
/// observation is left alone rather than guessed at. Catalogue stars are
/// reduced to the field once, so each observation searches a few dozen rather
/// than a few thousand.
///
/// Public because the calibration step has to redo it. With an uncorrected
/// focal length the sweep only reaches stars near the field centre: on the
/// `hard` preset that leaves 4 to 7 matches spanning a tenth of the field
/// radially, far too few and too clustered to fit a focal length from. Once a
/// first fit has improved the camera, running the sweep again picks up the rest
/// and the next fit is properly conditioned.
#[allow(clippy::too_many_arguments)]
pub fn extend_matches(
    observed: &[Vec3],
    attitude: &Mat3,
    camera: &CameraModel,
    catalog: &Catalog,
    db_star_count: usize,
    tolerance: f64,
    matches: &mut Vec<StarMatch>,
    ws: &mut Workspace,
) -> usize {
    collect_field(
        attitude,
        camera,
        catalog,
        db_star_count,
        tolerance,
        &mut ws.in_field,
    );
    match_field(observed, &ws.in_field, catalog, tolerance, matches)
}

/// Collects the database stars an attitude puts in the field.
///
/// Each entry is `(catalogue index, direction in the camera frame)`, with
/// `b = A r`. The reach is the field's half-diagonal widened by `tolerance`,
/// in radians, so a star just outside the corner is still a candidate.
pub fn collect_field(
    attitude: &Mat3,
    camera: &CameraModel,
    catalog: &Catalog,
    db_star_count: usize,
    tolerance: f64,
    out: &mut Vec<(u16, Vec3)>,
) {
    // The boresight's inertial direction is the third row of A, since b = A r.
    let boresight = attitude.row(2).transpose();
    let cos_reach = (camera.half_diagonal_fov() + tolerance).cos();

    out.clear();
    for index in 0..db_star_count {
        let direction = catalog.stars[index].direction();
        if boresight.dot(&direction) >= cos_reach {
            out.push((index as u16, attitude * direction));
        }
    }
}

/// Matches observed directions against an already-prepared field.
///
/// `field` is `(catalogue index, camera-frame direction)` as
/// [`collect_field`] produces. `tolerance` is radians. Appends to `matches`
/// and returns how many were added. Tracking prepares its own field and calls
/// this directly, which is why the two halves are separable.
pub fn match_field(
    observed: &[Vec3],
    field: &[(u16, Vec3)],
    catalog: &Catalog,
    tolerance: f64,
    matches: &mut Vec<StarMatch>,
) -> usize {
    // A star further than the tolerance can be neither the match nor an
    // ambiguity, so reject it on a dot product and leave the angle -- which
    // costs an atan2 -- for the few that survive. The guard band makes the
    // rejection provably looser than the test it stands in for, so this is an
    // optimisation and not a change of behaviour.
    let cos_reject = (tolerance * 1.001 + 1e-12).cos();

    let mut added = 0usize;
    for (observed_index, observation) in observed.iter().enumerate() {
        let observed_index = observed_index as u16;
        if matches.iter().any(|m| m.observed == observed_index) {
            continue;
        }

        // Nearest in-field catalogue star, and the runner-up, so an ambiguous
        // observation can be left alone rather than guessed at.
        let mut best: Option<(f64, u16)> = None;
        let mut second = f64::INFINITY;
        for &(catalog_index, predicted) in field {
            if observation.dot(&predicted) < cos_reject {
                continue;
            }
            let separation = angle_between(observation, &predicted);
            match best {
                Some((best_separation, _)) if separation >= best_separation => {
                    second = second.min(separation);
                }
                _ => {
                    if let Some((previous, _)) = best {
                        second = second.min(previous);
                    }
                    best = Some((separation, catalog_index));
                }
            }
        }

        let Some((separation, catalog_index)) = best else {
            continue;
        };
        if separation > tolerance || second <= tolerance {
            continue;
        }
        if matches.iter().any(|m| m.catalog_index == catalog_index) {
            continue;
        }

        matches.push(StarMatch {
            observed: observed_index,
            catalog_index,
            id: catalog.stars[usize::from(catalog_index)].id,
            from_pyramid: false,
        });
        added += 1;
    }
    added
}

/// An equally weighted Wahba solution over a match list.
///
/// The identifier's own attitude is provisional; the brightness-weighted
/// solution in `attitude` is what the pipeline reports.
/// Solves Wahba over a match list, with `b = A r` and brightness weights.
///
/// Returns `None` when the matches are too few or degenerate to fix a rotation.
pub fn attitude_over(matches: &[StarMatch], observed: &[Vec3], catalog: &Catalog) -> Option<Mat3> {
    let pairs: Vec<(Vec3, Vec3, f64)> = matches
        .iter()
        .map(|m| {
            (
                observed[usize::from(m.observed)],
                catalog.stars[usize::from(m.catalog_index)].direction(),
                1.0,
            )
        })
        .collect();
    solve_wahba(&pairs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NoClock;
    use crate::math::{quat_to_matrix, rotation_vector, splitmix64};
    use crate::simulate::Preset;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use rand_distr::StandardNormal;
    use std::collections::HashSet;
    use std::f64::consts::PI;

    /// One arcsecond in radians.
    const ARCSEC: f64 = PI / (180.0 * 3600.0);

    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

    fn catalog() -> Catalog {
        Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode")
    }

    /// A uniformly distributed attitude, built the way the simulator does.
    fn random_attitude(rng: &mut ChaCha8Rng) -> Mat3 {
        quat_to_matrix([
            rng.sample::<f64, _>(StandardNormal),
            rng.sample::<f64, _>(StandardNormal),
            rng.sample::<f64, _>(StandardNormal),
            rng.sample::<f64, _>(StandardNormal),
        ])
    }

    /// Perturbs a unit vector by a Gaussian angular error of `sigma` radians,
    /// applied in the two directions tangent to it.
    fn jitter(v: &Vec3, sigma: f64, rng: &mut ChaCha8Rng) -> Vec3 {
        let helper = if v.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let east = v.cross(&helper).normalize();
        let north = v.cross(&east);
        let a = rng.sample::<f64, _>(StandardNormal) * sigma;
        let b = rng.sample::<f64, _>(StandardNormal) * sigma;
        (v + east * a + north * b).normalize()
    }

    /// One vector-level frame: no image, just directions.
    struct Frame {
        observed: Vec<Vec3>,
        /// Catalogue index per observation, or `None` for an injected false
        /// vector.
        truth: Vec<Option<u16>>,
        attitude: Mat3,
    }

    /// Builds a frame of exact in-field catalogue directions, optionally
    /// perturbed, with `false_count` random directions mixed in.
    ///
    /// The false vectors are *inserted* at random positions rather than
    /// appended: a flux-ordered centroid list interleaves them with the real
    /// stars, and appending them would leave the first triple Mortari tries
    /// always clean, which would not test the ordering at all.
    fn frame(
        seed: u64,
        cfg: &SimConfig,
        catalog: &Catalog,
        db: &PairDb,
        sigma: f64,
        false_count: usize,
    ) -> Frame {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let attitude = random_attitude(&mut rng);
        let camera = cfg.nominal_camera();

        let mut visible: Vec<(f32, u16, Vec3)> = Vec::new();
        for index in 0..db.star_count() {
            let star = &catalog.stars[index];
            let body = attitude * star.direction();
            if let Some(pixel) = camera.project(&body)
                && camera.contains(pixel[0], pixel[1])
            {
                visible.push((star.mag, index as u16, body));
            }
        }
        visible.sort_by(|a, b| a.0.total_cmp(&b.0));
        visible.truncate(cfg.max_centroids.saturating_sub(false_count));

        let mut observed = Vec::new();
        let mut truth = Vec::new();
        for (_, index, body) in &visible {
            observed.push(if sigma > 0.0 {
                jitter(body, sigma, &mut rng)
            } else {
                *body
            });
            truth.push(Some(*index));
        }
        for _ in 0..false_count {
            let x = rng.random::<f64>() * (cfg.width as f64 - 1.0);
            let y = rng.random::<f64>() * (cfg.height as f64 - 1.0);
            let at = rng.random_range(0..=observed.len());
            observed.insert(at, camera.unproject(x, y));
            truth.insert(at, None);
        }
        Frame {
            observed,
            truth,
            attitude,
        }
    }

    /// Per-vector angular sigma the configuration implies.
    fn nominal_sigma(cfg: &SimConfig) -> f64 {
        cfg.centroid_sigma_px / cfg.nominal_camera().f
    }

    /// How a batch of vector-level trials turned out.
    #[derive(Default)]
    struct Tally {
        trials: u32,
        correct: u32,
        no_solution: u32,
        wrong: u32,
    }

    impl Tally {
        fn fraction_correct(&self) -> f64 {
            f64::from(self.correct) / f64::from(self.trials)
        }
    }

    /// Runs `trials` vector-level trials and counts the outcomes. A trial is
    /// correct when a solution came back, every claimed star is the right one,
    /// and the pyramid's four stars are all present.
    fn run_trials(trials: u32, sigma: f64, false_count: usize) -> Tally {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();
        let mut tally = Tally {
            trials,
            ..Default::default()
        };

        for trial in 0..u64::from(trials) {
            let f = frame(
                splitmix64(0xBEEF ^ trial),
                &cfg,
                &catalog,
                &db,
                sigma,
                false_count,
            );
            let outcome = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut ws);
            match outcome.solution {
                None => tally.no_solution += 1,
                Some(solution) => {
                    let all_right = solution
                        .matches
                        .iter()
                        .all(|m| f.truth[usize::from(m.observed)] == Some(m.catalog_index));
                    if all_right && solution.pyramid_count() == PYRAMID_STARS {
                        tally.correct += 1;
                    } else {
                        tally.wrong += 1;
                    }
                }
            }
        }
        tally
    }

    // --- the Mortari ordering ---

    #[test]
    fn mortari_order_covers_every_triple_exactly_once() {
        for n in 3..=16usize {
            let triples: Vec<_> = MortariTriples::new(n).collect();
            let expected = n * (n - 1) * (n - 2) / 6;
            assert_eq!(triples.len(), expected, "n = {n}");

            let unique: HashSet<_> = triples.iter().copied().collect();
            assert_eq!(unique.len(), expected, "n = {n} repeated a triple");
            for &(i, j, k) in &triples {
                assert!(i < j && j < k && k < n, "n = {n} gave ({i},{j},{k})");
            }
        }
        // Too few stars for a triple at all.
        assert_eq!(MortariTriples::new(2).count(), 0);
        assert_eq!(MortariTriples::new(0).count(), 0);
    }

    /// The property that makes the ordering worth having: a star at any single
    /// index must be absent from one of the first few triples, so one false
    /// star cannot block every early attempt.
    #[test]
    fn mortari_order_spreads_the_indices() {
        let n = 15;
        let first: Vec<_> = MortariTriples::new(n).take(4).collect();
        assert_eq!(first[0], (0, 1, 2));
        // The next tries move off index 0 entirely.
        assert_eq!(first[1], (1, 2, 3));
        assert_eq!(first[2], (2, 3, 4));
        for index in 0..n {
            assert!(
                first
                    .iter()
                    .any(|&(i, j, k)| i != index && j != index && k != index),
                "index {index} appears in all of the first four triples"
            );
        }
    }

    // --- tolerance ---

    #[test]
    fn tolerance_follows_the_documented_formula() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let tolerance = Tolerance::from_config(&cfg);

        let per_vector = cfg.centroid_sigma_px / cfg.nominal_camera().f;
        let margin = cfg.calib_margin_arcsec * ARCSEC;
        let expected_base = cfg.id_k_sigma * 2.0f64.sqrt() * per_vector + margin;
        assert!((tolerance.pair_base - expected_base).abs() < 1e-15);
        assert!((tolerance.vector - (cfg.id_k_sigma * per_vector + margin)).abs() < 1e-15);
        // The pair tolerance is the looser of the two, since two centroids
        // contribute to an angle between stars.
        assert!(tolerance.pair_base > tolerance.vector);
        // About 50 arcsec for the defaults, at zero separation.
        assert!((40.0..70.0).contains(&(tolerance.pair(0.0) / ARCSEC)));
        assert!((tolerance.pair(0.0) - tolerance.pair_base).abs() < 1e-15);

        // And it grows with the separation, by the configured fraction of it.
        // A focal length off by that fraction stretches an observed angle by
        // about the same, which is what this term has to cover.
        assert!((tolerance.pair_per_rad - cfg.id_focal_tolerance_frac).abs() < 1e-15);
        for degrees in [1.0, 10.0, 20.0] {
            let theta = degrees * PI / 180.0;
            let expected = expected_base + cfg.id_focal_tolerance_frac * theta;
            assert!((tolerance.pair(theta) - expected).abs() < 1e-15);
            assert!(tolerance.pair(theta) > tolerance.pair_base);
        }
        // Measured on the hard preset: 0.2% of focal error costs about 8
        // arcseconds of pair error per degree of separation, so the default
        // fraction has to buy at least that much.
        let per_degree = (tolerance.pair(PI / 180.0) - tolerance.pair_base) / ARCSEC;
        assert!(per_degree >= 8.0, "only {per_degree:.1}\" per degree");
    }

    // --- acceptance ---

    /// Phase 5 acceptance: 100% correct over 1000 noiseless vector-level
    /// trials.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "1000 trials, run with --release")]
    fn noiseless_vectors_identify_perfectly() {
        let tally = run_trials(1000, 0.0, 0);
        assert_eq!(
            tally.correct, tally.trials,
            "{} wrong and {} unsolved out of {}",
            tally.wrong, tally.no_solution, tally.trials
        );
    }

    /// Phase 5 acceptance: at least 99.5% correct with nominal noise.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "1000 trials, run with --release")]
    fn nominal_noise_identifies_at_least_995_per_cent() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let tally = run_trials(1000, nominal_sigma(&cfg), 0);
        assert!(
            tally.fraction_correct() >= 0.995,
            "{:.3}% correct ({} wrong, {} unsolved)",
            100.0 * tally.fraction_correct(),
            tally.wrong,
            tally.no_solution
        );
    }

    /// Phase 5 acceptance: at least 98% correct with three false vectors.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "1000 trials, run with --release")]
    fn three_false_vectors_still_identify() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let tally = run_trials(1000, nominal_sigma(&cfg), 3);
        assert!(
            tally.fraction_correct() >= 0.98,
            "{:.3}% correct ({} wrong, {} unsolved)",
            100.0 * tally.fraction_correct(),
            tally.wrong,
            tally.no_solution
        );
    }

    /// Degrading gracefully matters as much as the headline rate: the solver
    /// must not fall off a cliff as false stars pile up.
    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "several batches of trials, run with --release"
    )]
    fn accuracy_degrades_gradually_with_more_false_vectors() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let sigma = nominal_sigma(&cfg);
        let six = run_trials(300, sigma, 6);
        let nine = run_trials(300, sigma, 9);
        assert!(
            six.fraction_correct() >= 0.97,
            "six false vectors gave {:.2}%",
            100.0 * six.fraction_correct()
        );
        assert!(
            nine.fraction_correct() >= 0.95,
            "nine false vectors gave {:.2}%",
            100.0 * nine.fraction_correct()
        );
    }

    // --- what the solution contains ---

    #[test]
    fn a_solution_is_internally_consistent() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();

        for trial in 0..20u64 {
            let f = frame(splitmix64(trial), &cfg, &catalog, &db, 0.0, 0);
            let outcome = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut ws);
            let solution = outcome.solution.expect("noiseless frames must solve");

            assert_eq!(solution.pyramid_count(), PYRAMID_STARS);
            let mut observed_seen = HashSet::new();
            let mut catalog_seen = HashSet::new();
            for m in &solution.matches {
                assert!(usize::from(m.observed) < f.observed.len());
                assert!(usize::from(m.catalog_index) < db.star_count());
                // The id and the index must name the same star.
                assert_eq!(m.id, catalog.stars[usize::from(m.catalog_index)].id);
                assert!(
                    observed_seen.insert(m.observed),
                    "observation matched twice"
                );
                assert!(
                    catalog_seen.insert(m.catalog_index),
                    "catalogue star reused"
                );
            }
            // The four pyramid members come first.
            assert!(
                solution.matches[..PYRAMID_STARS]
                    .iter()
                    .all(|m| m.from_pyramid)
            );
        }
    }

    /// With exact vectors and correct identifications, the attitude should come
    /// back essentially exactly. Phase 6 characterises it under noise.
    #[test]
    fn noiseless_attitude_is_recovered_exactly() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();

        let mut worst = 0.0f64;
        for trial in 0..50u64 {
            let f = frame(splitmix64(0xA77 ^ trial), &cfg, &catalog, &db, 0.0, 0);
            let solution = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut ws)
                .solution
                .expect("noiseless frames must solve");
            let error = rotation_vector(&(solution.attitude * f.attitude.transpose())).norm();
            worst = worst.max(error);
        }
        assert!(
            worst < 1e-9,
            "worst noiseless attitude error {:.3e} rad ({:.4} arcsec)",
            worst,
            worst / ARCSEC
        );
    }

    /// Attitude error under noise, against what the geometry predicts.
    ///
    /// Cross-boresight error should sit near `sigma / sqrt(N)`. Roll is far
    /// worse in a narrow field, because the stars give it only a short lever
    /// arm: dividing by `sin(half-diagonal)` -- about 0.24 at 28 degrees
    /// diagonal -- is the expected penalty. Medians are used rather than a
    /// worst case, since roll has a long tail: measured over 3000 trials the
    /// median total is 9.8 arcsec but the maximum reaches 61.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "300 trials, run with --release")]
    fn nominal_noise_attitude_error_matches_the_predicted_scaling() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();
        let sigma = nominal_sigma(&cfg);

        let mut cross = Vec::new();
        let mut roll = Vec::new();
        let mut stars = 0usize;
        for trial in 0..300u64 {
            let f = frame(splitmix64(0xB88 ^ trial), &cfg, &catalog, &db, sigma, 0);
            let solution = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut ws)
                .solution
                .expect("must solve");
            let phi = rotation_vector(&(solution.attitude * f.attitude.transpose()));
            cross.push(phi.x.hypot(phi.y));
            roll.push(phi.z.abs());
            stars += solution.matches.len();
        }

        let median = |values: &mut Vec<f64>| {
            values.sort_by(f64::total_cmp);
            values[values.len() / 2]
        };
        let cross_median = median(&mut cross);
        let roll_median = median(&mut roll);

        let mean_stars = stars as f64 / 300.0;
        let predicted_cross = sigma / mean_stars.sqrt();
        let predicted_roll = predicted_cross / cfg.nominal_camera().half_diagonal_fov().sin();

        for (label, measured, predicted) in [
            ("cross-boresight", cross_median, predicted_cross),
            ("roll", roll_median, predicted_roll),
        ] {
            let ratio = measured / predicted;
            assert!(
                (0.5..2.0).contains(&ratio),
                "{label} median {:.2} arcsec against a predicted {:.2} \
                 (ratio {ratio:.2}, wanted 0.5 to 2)",
                measured / ARCSEC,
                predicted / ARCSEC
            );
        }
        // Roll really is the dominant term at this field of view.
        assert!(roll_median > 2.0 * cross_median);
    }

    // --- diagnostics, caps and edge cases ---

    #[test]
    fn diagnostics_come_back_even_when_nothing_is_found() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();

        // Three vectors cannot make a pyramid.
        let outcome = identify(
            &[Vec3::z(), Vec3::x(), Vec3::y()],
            &catalog,
            &db,
            &cfg,
            &NoClock,
            &mut ws,
        );
        assert!(outcome.solution.is_none());
        assert_eq!(outcome.diagnostics.tries, 0);

        // A solvable frame reports real counters.
        let f = frame(1, &cfg, &catalog, &db, 0.0, 0);
        let outcome = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut ws);
        assert!(outcome.solution.is_some());
        assert!(outcome.diagnostics.tries >= 1);
        assert!(outcome.diagnostics.pair_queries >= 3);
        assert!(outcome.diagnostics.pairs_examined > 0);
        assert_eq!(outcome.diagnostics.confirmed_pyramids, 1);
        assert!(!outcome.diagnostics.hit_try_cap);
    }

    #[test]
    fn the_iteration_cap_is_honoured() {
        let catalog = catalog();
        let cfg = SimConfig {
            id_max_tries: 5,
            ..SimConfig::preset(Preset::Nominal)
        };
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();

        // Fifteen directions that are not a star field: the search will run
        // out of tries rather than out of triples.
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let camera = cfg.nominal_camera();
        let observed: Vec<Vec3> = (0..15)
            .map(|_| camera.unproject(rng.random::<f64>() * 1023.0, rng.random::<f64>() * 1023.0))
            .collect();

        let outcome = identify(&observed, &catalog, &db, &cfg, &NoClock, &mut ws);
        assert!(outcome.diagnostics.tries <= 5);
        if outcome.solution.is_none() {
            assert!(outcome.diagnostics.hit_try_cap);
        }
    }

    /// A clock is only read twice per call, so the elapsed figure must follow
    /// whatever the host supplies.
    #[test]
    fn elapsed_time_comes_from_the_supplied_clock() {
        struct Fixed(std::cell::Cell<u64>);
        impl Clock for Fixed {
            fn now_ns(&self) -> u64 {
                let value = self.0.get();
                self.0.set(value + 1_234);
                value
            }
        }

        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();
        let f = frame(2, &cfg, &catalog, &db, 0.0, 0);

        let clock = Fixed(std::cell::Cell::new(0));
        let outcome = identify(&f.observed, &catalog, &db, &cfg, &clock, &mut ws);
        assert_eq!(outcome.diagnostics.elapsed_ns, 1_234);

        // And `NoClock` reports nothing rather than something wrong.
        let outcome = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut ws);
        assert_eq!(outcome.diagnostics.elapsed_ns, 0);
    }

    #[test]
    fn repeated_calls_agree_and_a_reused_workspace_changes_nothing() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let f = frame(77, &cfg, &catalog, &db, nominal_sigma(&cfg), 2);

        let mut fresh = Workspace::new();
        let first = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut fresh)
            .solution
            .expect("solves");

        let mut reused = Workspace::new();
        for seed in 0..5u64 {
            let other = frame(seed, &cfg, &catalog, &db, 0.0, 0);
            identify(&other.observed, &catalog, &db, &cfg, &NoClock, &mut reused);
        }
        let second = identify(&f.observed, &catalog, &db, &cfg, &NoClock, &mut reused)
            .solution
            .expect("solves");

        assert_eq!(first.matches, second.matches);
        assert_eq!(first.attitude, second.attitude);
    }

    #[test]
    fn an_empty_observation_list_is_handled() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();
        assert!(
            identify(&[], &catalog, &db, &cfg, &NoClock, &mut ws)
                .solution
                .is_none()
        );
    }
}
