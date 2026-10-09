//! Verification: two independent layers, then an outcome.
//!
//! The **self-check** is truth-blind and is what a real tracker would do with
//! no ground truth available: reproject the catalogue through the estimated
//! attitude, count how many observed centroids that explains, and measure the
//! residual. It decides whether the solution is CONFIDENT.
//!
//! The **ground-truth check** uses the simulator's truth to ask whether the
//! solution was actually right. It reads only [`SimFrame`] truth and the
//! solver's output, and never touches identification internals -- otherwise it
//! could not catch a mistake the identifier made.
//!
//! Keeping the two apart is the whole point: the self-check is the mechanism a
//! flight tracker ships, and the ground-truth check is how we find out whether
//! that mechanism is trustworthy.

use crate::camera::CameraModel;
use crate::catalog::{Catalog, IdIndex};
use crate::centroid::Centroid;
use crate::identify::{Solution, StarMatch};
use crate::math::{Mat3, Vec3, rotation_vector};
use crate::simulate::{SimConfig, SimFrame};

/// A pyramid's four stars: the floor the self-check insists on.
const PYRAMID_STARS: usize = 4;

/// How a trial ended, as docs/SPEC.md tabulates it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Outcome {
    /// CONFIDENT, every claimed id right, and the attitude error below
    /// `err_threshold`.
    Correct,
    /// CONFIDENT but with a wrong id or a large error. The worst case, and the
    /// one whose target rate is zero.
    WrongConfident,
    /// A solution was found and the self-check refused it.
    Rejected,
    /// Identification found nothing at all.
    NoSolution,
}

impl Outcome {
    /// Every outcome, in the order reports list them.
    pub const ALL: [Outcome; 4] = [
        Outcome::Correct,
        Outcome::WrongConfident,
        Outcome::Rejected,
        Outcome::NoSolution,
    ];

    /// The name used in reports, CSV and JSON.
    pub fn as_str(self) -> &'static str {
        match self {
            Outcome::Correct => "CORRECT",
            Outcome::WrongConfident => "WRONG_CONFIDENT",
            Outcome::Rejected => "REJECTED",
            Outcome::NoSolution => "NO_SOLUTION",
        }
    }
}

impl std::fmt::Display for Outcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// What the truth-blind self-check concluded, and the numbers behind it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SelfCheck {
    /// Observed centroids explained by a reprojected catalogue star.
    pub matched: usize,
    /// RMS of those reprojection residuals, pixels.
    pub residual_rms_px: f64,
    /// How many matches the pyramid itself confirmed.
    pub pyramid_stars: usize,
    /// Claims whose own catalogue star reprojects onto the centroid it was
    /// claimed for.
    pub claims_explained: usize,
    /// Claims made in total.
    pub claims: usize,
    /// Every condition held.
    pub confident: bool,
}

/// What the ground-truth check found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TruthCheck {
    /// Catalogue ids the solver claimed.
    pub claimed: usize,
    /// How many of those claims were right.
    pub correct: usize,
    /// Visible, rendered truth stars that the database could have identified:
    /// the denominator for recall.
    pub available: usize,
    /// Total attitude error, radians.
    pub total_error: f64,
    /// Cross-boresight component, radians.
    pub cross_error: f64,
    /// Roll component, radians.
    pub roll_error: f64,
}

impl TruthCheck {
    /// Whether every claim was right.
    pub fn all_ids_right(&self) -> bool {
        self.claimed > 0 && self.correct == self.claimed
    }
}

/// Scratch buffers reused across trials.
#[derive(Clone, Debug, Default)]
pub struct Workspace {
    /// Catalogue stars inside the estimated field, as predicted pixels.
    predicted: Vec<[f64; 2]>,
}

impl Workspace {
    /// An empty workspace; buffers size themselves on first use.
    pub fn new() -> Self {
        Self::default()
    }
}

/// Squared pixel distance between a predicted pixel and a centroid.
///
/// These comparisons sit in loops that run tens of thousands of times a frame,
/// and `hypot` is a careful, slow libm routine -- on wasm32 a software one,
/// where it dominated verification at seventeen times the native cost. Squares
/// order exactly as distances do, so comparing them against a squared threshold
/// decides the same way without the square root.
fn distance_squared(pixel: [f64; 2], x: f64, y: f64) -> f64 {
    let dx = pixel[0] - x;
    let dy = pixel[1] - y;
    dx * dx + dy * dy
}

/// The truth-blind self-check.
///
/// Reprojects every database catalogue star through `attitude`, counts the
/// observed centroids that lands on, and measures the residual RMS. A solution
/// is CONFIDENT only when at least `verify_min_matches` centroids are
/// explained, the residual RMS is within `verify_max_rms_px`, at least four
/// stars came from the pyramid, and **every** claim reprojects onto the
/// centroid it was made for.
///
/// `db_star_count` is how many of the catalogue's leading entries the solver's
/// database covers; stars past it are invisible to the solver and must not
/// count towards its own confidence.
///
/// `camera` is whatever the solver believes it has, which after calibration is
/// the nominal camera with a fitted focal length. Using the nominal one here
/// instead would charge the solution for a focal error it has already measured
/// and corrected.
#[allow(clippy::too_many_arguments)]
pub fn self_check(
    attitude: &Mat3,
    matches: &[StarMatch],
    centroids: &[Centroid],
    catalog: &Catalog,
    camera: &CameraModel,
    db_star_count: usize,
    cfg: &SimConfig,
    ws: &mut Workspace,
) -> SelfCheck {
    // The boresight's inertial direction is the third row of A, since b = A r.
    let boresight = attitude.row(2).transpose();
    let cos_reach = camera.half_diagonal_fov().cos();

    ws.predicted.clear();
    for star in catalog.stars.iter().take(db_star_count) {
        let direction = star.direction();
        if boresight.dot(&direction) < cos_reach {
            continue;
        }
        if let Some(pixel) = camera.project(&(attitude * direction))
            && camera.contains(pixel[0], pixel[1])
        {
            ws.predicted.push(pixel);
        }
    }

    let match_limit_squared = cfg.verify_match_px * cfg.verify_match_px;
    let mut matched = 0usize;
    let mut sum_squares = 0.0;
    for centroid in centroids {
        let mut nearest = f64::INFINITY;
        for pixel in &ws.predicted {
            nearest = nearest.min(distance_squared(*pixel, centroid.x, centroid.y));
        }
        if nearest <= match_limit_squared {
            matched += 1;
            // Already a square, which is what the running sum wants.
            sum_squares += nearest;
        }
    }

    let residual_rms_px = if matched > 0 {
        (sum_squares / matched as f64).sqrt()
    } else {
        f64::INFINITY
    };
    let pyramid_stars = matches.iter().filter(|m| m.from_pyramid).count();

    // Each claim, checked on its own terms: does the star the solver named
    // actually reproject onto the centroid it named it for? The residual RMS
    // alone will not catch this. One wrong id among fifteen leaves the RMS at
    // 0.30 px, under the 0.5 px limit, because the fourteen good matches
    // dominate it -- which is exactly how a wrong solution slipped through as
    // CONFIDENT on the golden run.
    let mut claims_explained = 0usize;
    for claim in matches {
        let Some(centroid) = centroids.get(usize::from(claim.observed)) else {
            continue;
        };
        let Some(star) = catalog.stars.get(usize::from(claim.catalog_index)) else {
            continue;
        };
        if let Some(pixel) = camera.project(&(attitude * star.direction()))
            && distance_squared(pixel, centroid.x, centroid.y) <= match_limit_squared
        {
            claims_explained += 1;
        }
    }

    SelfCheck {
        matched,
        residual_rms_px,
        pyramid_stars,
        claims_explained,
        claims: matches.len(),
        confident: matched >= cfg.verify_min_matches
            && residual_rms_px <= cfg.verify_max_rms_px
            && pyramid_stars >= PYRAMID_STARS
            && claims_explained == matches.len(),
    }
}

/// The ground-truth check.
///
/// A claim is right when the truth star carrying the claimed id was actually
/// rendered and its true position lies within `truth_match_px` of the centroid
/// the claim was made about. Asking it that way round avoids the ambiguity of
/// a nearest-truth-star search, and treats a blended centroid that sits well
/// away from the star it names as the wrong claim it is.
///
/// Reads only `frame`'s truth and the solver's output.
pub fn truth_check(
    attitude: &Mat3,
    matches: &[StarMatch],
    centroids: &[Centroid],
    frame: &SimFrame,
    ids: &IdIndex,
    db_star_count: usize,
    cfg: &SimConfig,
) -> TruthCheck {
    let truth_limit_squared = cfg.truth_match_px * cfg.truth_match_px;
    let mut correct = 0usize;
    for claim in matches {
        let Some(centroid) = centroids.get(usize::from(claim.observed)) else {
            continue;
        };
        let right = frame.truth.iter().any(|star| {
            !star.dropped
                && star.id == claim.id
                && distance_squared(star.pixel, centroid.x, centroid.y) <= truth_limit_squared
        });
        if right {
            correct += 1;
        }
    }

    // Recall's denominator: rendered truth stars the database actually holds.
    let available = frame
        .rendered()
        .filter(|star| ids.get(star.id).is_some_and(|index| index < db_star_count))
        .count();

    let phi = rotation_vector(&(attitude * frame.a_true.transpose()));
    TruthCheck {
        claimed: matches.len(),
        correct,
        available,
        total_error: phi.norm(),
        cross_error: phi.x.hypot(phi.y),
        roll_error: phi.z.abs(),
    }
}

/// Classifies a trial from its two checks.
///
/// `solution` is `None` when identification found nothing. Note that a
/// CONFIDENT solution with a large error is WRONG_CONFIDENT, not REJECTED: the
/// self-check having passed is exactly what makes it the dangerous case.
pub fn classify(
    solution: Option<&Solution>,
    self_check: Option<&SelfCheck>,
    truth: Option<&TruthCheck>,
    cfg: &SimConfig,
) -> Outcome {
    let (Some(_), Some(self_check), Some(truth)) = (solution, self_check, truth) else {
        return Outcome::NoSolution;
    };
    if !self_check.confident {
        return Outcome::Rejected;
    }
    let threshold = cfg.err_threshold_arcsec * std::f64::consts::PI / (180.0 * 3600.0);
    if truth.all_ids_right() && truth.total_error <= threshold {
        Outcome::Correct
    } else {
        Outcome::WrongConfident
    }
}

/// Camera-frame directions for a set of centroids, through `camera`.
///
/// The one conversion from measured pixels to the unit vectors every later
/// stage works in.
pub fn directions(centroids: &[Centroid], camera: &CameraModel, out: &mut Vec<Vec3>) {
    out.clear();
    out.extend(
        centroids
            .iter()
            .map(|centroid| camera.unproject(centroid.x, centroid.y)),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NoClock;
    use crate::camera::CameraModel;
    use crate::centroid;
    use crate::identify;
    use crate::math::splitmix64;
    use crate::pairdb::PairDb;
    use crate::simulate::{self, Preset};

    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

    fn catalog() -> Catalog {
        Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode")
    }

    /// Everything a solved frame produces, for the checks to work on.
    struct Solved {
        frame: SimFrame,
        centroids: Vec<Centroid>,
        solution: crate::identify::Solution,
        camera: CameraModel,
    }

    /// Runs the pipeline up to identification on one seed.
    fn solve(seed: u64, preset: Preset, catalog: &Catalog, db: &PairDb) -> Option<Solved> {
        let cfg = SimConfig::preset(preset);
        let camera = cfg.nominal_camera();
        let frame = simulate::simulate(seed, &cfg, catalog, &mut simulate::Workspace::new());
        let mut centroid_ws = centroid::Workspace::new();
        let detections = centroid::detect(
            &frame.image,
            frame.width,
            frame.height,
            &cfg,
            &mut centroid_ws,
        );
        let centroids = centroid::brightest(detections, &cfg).to_vec();
        let mut rays = Vec::new();
        directions(&centroids, &camera, &mut rays);
        let solution = identify::identify(
            &rays,
            catalog,
            db,
            &cfg,
            &NoClock,
            &mut identify::Workspace::new(),
        )
        .solution?;
        Some(Solved {
            frame,
            centroids,
            solution,
            camera,
        })
    }

    // --- outcome ---

    #[test]
    fn outcome_names_are_stable() {
        assert_eq!(Outcome::Correct.as_str(), "CORRECT");
        assert_eq!(Outcome::WrongConfident.as_str(), "WRONG_CONFIDENT");
        assert_eq!(Outcome::Rejected.as_str(), "REJECTED");
        assert_eq!(Outcome::NoSolution.as_str(), "NO_SOLUTION");
        assert_eq!(Outcome::ALL.len(), 4);
        for outcome in Outcome::ALL {
            assert_eq!(outcome.to_string(), outcome.as_str());
        }
    }

    // --- classification ---

    #[test]
    fn classification_covers_every_outcome() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let threshold = cfg.err_threshold_arcsec * std::f64::consts::PI / (180.0 * 3600.0);

        let confident = SelfCheck {
            matched: 10,
            residual_rms_px: 0.05,
            pyramid_stars: 4,
            claims_explained: 10,
            claims: 10,
            confident: true,
        };
        let refused = SelfCheck {
            confident: false,
            ..confident
        };
        let right = TruthCheck {
            claimed: 10,
            correct: 10,
            available: 20,
            total_error: threshold / 2.0,
            cross_error: 0.0,
            roll_error: 0.0,
        };
        let wrong_id = TruthCheck {
            correct: 9,
            ..right
        };
        let wild = TruthCheck {
            total_error: threshold * 10.0,
            ..right
        };

        // A solution object is only needed to show one exists.
        let solution = crate::identify::Solution {
            matches: Vec::new(),
            attitude: Mat3::identity(),
        };
        let at = |s: &SelfCheck, t: &TruthCheck| classify(Some(&solution), Some(s), Some(t), &cfg);

        assert_eq!(at(&confident, &right), Outcome::Correct);
        assert_eq!(at(&confident, &wrong_id), Outcome::WrongConfident);
        assert_eq!(at(&confident, &wild), Outcome::WrongConfident);
        assert_eq!(at(&refused, &right), Outcome::Rejected);
        // Refused beats wrong: the self-check caught it, which is the point.
        assert_eq!(at(&refused, &wrong_id), Outcome::Rejected);
        assert_eq!(classify(None, None, None, &cfg), Outcome::NoSolution);
    }

    // --- the self-check ---

    #[test]
    fn a_good_solution_passes_the_self_check() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let mut ws = Workspace::new();

        for trial in 0..20u64 {
            let solved = solve(splitmix64(trial), Preset::Nominal, &catalog, &db)
                .expect("nominal frames solve");
            let check = self_check(
                &solved.solution.attitude,
                &solved.solution.matches,
                &solved.centroids,
                &catalog,
                &solved.camera,
                db.star_count(),
                &cfg,
                &mut ws,
            );
            assert!(check.confident, "{check:?}");
            assert!(check.matched >= cfg.verify_min_matches);
            assert!(check.residual_rms_px <= cfg.verify_max_rms_px);
            assert_eq!(check.claims_explained, check.claims);
            assert!(check.pyramid_stars >= 4);
        }
    }

    #[test]
    fn a_rotated_attitude_fails_the_self_check() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let solved = solve(7, Preset::Nominal, &catalog, &db).expect("solves");

        // Ten degrees off: the catalogue no longer explains the frame.
        let spoiled = crate::math::quat_to_matrix([
            (5.0f64).to_radians().cos(),
            (5.0f64).to_radians().sin(),
            0.0,
            0.0,
        ]) * solved.solution.attitude;
        let check = self_check(
            &spoiled,
            &solved.solution.matches,
            &solved.centroids,
            &catalog,
            &solved.camera,
            db.star_count(),
            &cfg,
            &mut Workspace::new(),
        );
        assert!(!check.confident, "a 10 degree error passed: {check:?}");
    }

    /// The regression that drove the per-claim check: on the golden run one
    /// wrong id among fifteen left the residual RMS at 0.30 px, inside the
    /// 0.5 px limit, because the fourteen good matches dominated it.
    #[test]
    fn one_wrong_claim_fails_the_self_check_even_with_a_fine_residual() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let solved = solve(11, Preset::Nominal, &catalog, &db).expect("solves");

        // Point one claim at a different catalogue star, leaving the rest and
        // the attitude untouched.
        let mut spoiled = solved.solution.matches.clone();
        let last = spoiled.len() - 1;
        let wrong_index = (spoiled[last].catalog_index + 500) % db.star_count() as u16;
        spoiled[last].catalog_index = wrong_index;
        spoiled[last].id = catalog.stars[usize::from(wrong_index)].id;

        let check = self_check(
            &solved.solution.attitude,
            &spoiled,
            &solved.centroids,
            &catalog,
            &solved.camera,
            db.star_count(),
            &cfg,
            &mut Workspace::new(),
        );
        assert_eq!(check.claims_explained, check.claims - 1);
        assert!(
            check.residual_rms_px <= cfg.verify_max_rms_px,
            "the residual should still look fine, which is the whole problem"
        );
        assert!(!check.confident, "a wrong claim was accepted: {check:?}");
    }

    #[test]
    fn too_few_pyramid_stars_fails_the_self_check() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let solved = solve(3, Preset::Nominal, &catalog, &db).expect("solves");

        let mut demoted = solved.solution.matches.clone();
        for matched in &mut demoted {
            matched.from_pyramid = false;
        }
        let check = self_check(
            &solved.solution.attitude,
            &demoted,
            &solved.centroids,
            &catalog,
            &solved.camera,
            db.star_count(),
            &cfg,
            &mut Workspace::new(),
        );
        assert_eq!(check.pyramid_stars, 0);
        assert!(!check.confident);
    }

    // --- the ground-truth check ---

    #[test]
    fn truth_check_scores_a_good_solution_perfectly() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let ids = catalog.id_index();

        for trial in 0..20u64 {
            let solved =
                solve(splitmix64(0x7777 ^ trial), Preset::Nominal, &catalog, &db).expect("solves");
            let truth = truth_check(
                &solved.solution.attitude,
                &solved.solution.matches,
                &solved.centroids,
                &solved.frame,
                &ids,
                db.star_count(),
                &cfg,
            );
            assert_eq!(truth.claimed, solved.solution.matches.len());
            assert_eq!(truth.correct, truth.claimed, "a claim was wrong");
            assert!(truth.all_ids_right());
            // Recall's denominator must at least cover what was claimed.
            assert!(truth.available >= truth.correct);
            // A good solution is accurate.
            assert!(truth.total_error < 1e-3, "{} rad", truth.total_error);
            // Cross-boresight and roll must account for the total.
            let combined = truth.cross_error.hypot(truth.roll_error);
            assert!((combined - truth.total_error).abs() < 1e-12);
        }
    }

    #[test]
    fn truth_check_catches_a_wrong_claim() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Nominal);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let ids = catalog.id_index();
        let solved = solve(5, Preset::Nominal, &catalog, &db).expect("solves");

        let mut spoiled = solved.solution.matches.clone();
        spoiled[0].id = catalog.stars[4000].id;
        let truth = truth_check(
            &solved.solution.attitude,
            &spoiled,
            &solved.centroids,
            &solved.frame,
            &ids,
            db.star_count(),
            &cfg,
        );
        assert_eq!(truth.correct, truth.claimed - 1);
        assert!(!truth.all_ids_right());
    }

    #[test]
    fn an_empty_claim_list_is_not_all_right() {
        let empty = TruthCheck {
            claimed: 0,
            correct: 0,
            available: 10,
            total_error: 0.0,
            cross_error: 0.0,
            roll_error: 0.0,
        };
        assert!(!empty.all_ids_right());
    }

    // --- directions ---

    #[test]
    fn directions_round_trip_through_the_camera() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let camera = cfg.nominal_camera();
        let centroids = [
            Centroid {
                x: 100.0,
                y: 200.0,
                flux: 1.0,
                pixels: 9,
            },
            Centroid {
                x: 511.5,
                y: 511.5,
                flux: 2.0,
                pixels: 9,
            },
        ];
        let mut out = Vec::new();
        directions(&centroids, &camera, &mut out);
        assert_eq!(out.len(), 2);
        for (centroid, direction) in centroids.iter().zip(&out) {
            assert!((direction.norm() - 1.0).abs() < 1e-15);
            let back = camera.project(direction).expect("in front");
            assert!((back[0] - centroid.x).abs() < 1e-9);
            assert!((back[1] - centroid.y).abs() < 1e-9);
        }
        // The boresight centroid maps onto +z.
        assert!(out[1].z > 0.999999);

        // The buffer is reused, not appended to.
        directions(&centroids[..1], &camera, &mut out);
        assert_eq!(out.len(), 1);
    }
}
