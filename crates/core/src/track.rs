//! Tracking mode: identification from a prior attitude.
//!
//! Lost-in-space identification has to search the whole sky, which is what the
//! pyramid is for. Once an attitude is known the next frame is a far smaller
//! problem: the camera has barely moved, so every star it is about to see is
//! already predictable, and identification collapses to matching each observed
//! direction against the predicted field. No pair database, no triples, no
//! fourth-star confirmation.
//!
//! The prior here is the previous frame's solution, used as it stands rather
//! than propagated forward by a known rate. That is the weaker assumption --
//! a tracker that knows its own angular rate can predict better and search a
//! smaller window -- so the search window instead has to cover the whole slew
//! between frames.
//!
//! Frames and units: attitudes map inertial to camera, `b = A r`. All angles
//! are radians and all times are seconds, except the nanosecond counters the
//! host clock supplies.

use serde::{Deserialize, Serialize};

use crate::Clock;
use crate::bench::{Solver, Workspaces, trial_seed};
use crate::catalog::Catalog;
use crate::centroid;
use crate::identify::{self, Solution, StarMatch};
use crate::math::{Mat3, Vec3, matrix_to_quat, quat_to_matrix};
use crate::simulate::{self, SimConfig};
use crate::verify::{self, Outcome, SelfCheck};

/// Fewest matched stars a tracked frame needs before its attitude is fitted.
///
/// Two directions already fix a rotation, but a tracked solution is only as
/// trustworthy as the prior that seeded it, so this matches the four a pyramid
/// must confirm.
pub const MIN_TRACK_MATCHES: usize = 4;

/// How a slew is simulated and how wide the tracker searches.
///
/// Separate from [`SimConfig`], which describes one frame: these parameters
/// only mean anything for a sequence of them.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub struct TrackConfig {
    /// Frames in the sequence, the first of which is acquired lost-in-space.
    pub frames: usize,
    /// Body rate, degrees per second.
    pub rate_deg_s: f64,
    /// Time between frames, seconds.
    pub interval_s: f64,
    /// Body axis the slew turns about; normalised on use.
    ///
    /// The default is the camera's +x, so the field drifts along a column and
    /// new stars keep entering it. A slew about +z would spin about the
    /// boresight and keep the same stars in view, which tests much less.
    pub axis: [f64; 3],
    /// Added to the search window, arcseconds, for the prior's own error.
    pub search_margin_arcsec: f64,
    /// How far the boresight may drift before the predicted field is rebuilt,
    /// degrees.
    ///
    /// Rebuilding scans the whole database, which is the one expensive thing a
    /// tracked frame would otherwise do every time. The field moves slowly, so
    /// the list is built with this much slack around the boresight and reused
    /// until the boresight leaves it. Larger slack means rarer rebuilds but a
    /// longer list to match against on every frame.
    pub field_slack_deg: f64,
}

impl Default for TrackConfig {
    fn default() -> Self {
        Self {
            frames: 60,
            rate_deg_s: 0.5,
            interval_s: 0.1,
            axis: [1.0, 0.0, 0.0],
            search_margin_arcsec: 60.0,
            field_slack_deg: 1.0,
        }
    }
}

impl TrackConfig {
    /// Angle the body turns between two frames, radians.
    pub fn step_rad(&self) -> f64 {
        self.rate_deg_s.to_radians() * self.interval_s
    }

    /// The rotation between consecutive frames, in the body frame.
    pub fn step(&self) -> Mat3 {
        let axis = Vec3::new(self.axis[0], self.axis[1], self.axis[2]);
        let norm = axis.norm();
        if norm < f64::EPSILON {
            return Mat3::identity();
        }
        let half = self.step_rad() / 2.0;
        let unit = axis / norm;
        let sin = half.sin();
        quat_to_matrix([half.cos(), unit.x * sin, unit.y * sin, unit.z * sin])
    }

    /// Search window when only one past attitude is known, radians.
    ///
    /// It has to cover the whole inter-frame slew, because a single attitude
    /// says nothing about where the next frame will be.
    pub fn search_radius(&self, cfg: &SimConfig) -> f64 {
        self.step_rad() + self.aided_radius(cfg)
    }

    /// Search window once the rate is known, radians.
    ///
    /// Two past attitudes give the rotation between them, and applying it once
    /// more predicts the next frame. What is left to search for is then only
    /// the error in that prediction -- how much the rate changed, plus the
    /// noise in the attitudes it was derived from -- not the slew itself. That
    /// is the difference between a window of a few arcminutes and one of
    /// degrees, and a wide window is not merely slower: every observation in it
    /// has a rival, so the ambiguity check throws the match away.
    pub fn aided_radius(&self, cfg: &SimConfig) -> f64 {
        let tolerance = identify::Tolerance::from_config(cfg);
        let margin = self.search_margin_arcsec * std::f64::consts::PI / (180.0 * 3600.0);
        tolerance.vector + margin
    }
}

/// The predicted field, kept across frames.
///
/// Stars are held as inertial directions so the list survives the attitude
/// changing under it; only the rotation into the camera frame is redone each
/// frame, and that is over the few dozen stars in view rather than the whole
/// database.
#[derive(Debug, Default)]
pub struct FieldCache {
    /// `(catalogue index, inertial direction)` for stars around `built_at`.
    stars: Vec<(u16, Vec3)>,
    /// Boresight the list was built around, inertial. `None` until first built.
    built_at: Option<Vec3>,
    /// The same stars rotated into the camera frame, rebuilt every frame.
    projected: Vec<(u16, Vec3)>,
    /// How many times the list has been rebuilt; a tracking cost worth seeing.
    pub rebuilds: usize,
}

impl FieldCache {
    /// An empty cache; it fills on first use.
    pub fn new() -> Self {
        Self::default()
    }

    /// The field around `attitude`, in the camera frame, rebuilt if needed.
    ///
    /// `tolerance` and `slack` are radians.
    pub fn field(
        &mut self,
        attitude: &Mat3,
        camera: &crate::camera::CameraModel,
        catalog: &Catalog,
        db_star_count: usize,
        tolerance: f64,
        slack: f64,
    ) -> &[(u16, Vec3)] {
        // The boresight's inertial direction is the third row of A, since
        // b = A r.
        let boresight = attitude.row(2).transpose();
        let stale = match self.built_at {
            Some(built_at) => crate::math::angle_between(&boresight, &built_at) > slack,
            None => true,
        };
        if stale {
            let cos_reach = (camera.half_diagonal_fov() + tolerance + slack).cos();
            self.stars.clear();
            for index in 0..db_star_count {
                let direction = catalog.stars[index].direction();
                if boresight.dot(&direction) >= cos_reach {
                    self.stars.push((index as u16, direction));
                }
            }
            self.built_at = Some(boresight);
            self.rebuilds += 1;
        }

        self.projected.clear();
        self.projected.extend(
            self.stars
                .iter()
                .map(|&(index, direction)| (index, attitude * direction)),
        );
        &self.projected
    }
}

/// Identifies a frame from a prior attitude, searching only the predicted field.
///
/// `prior` maps inertial to camera and is where this frame is expected to be
/// pointing. `search_radius` is how far a star may be from where `prior` puts
/// it, in radians. Returns `None` when too few stars match to fit an attitude.
#[allow(clippy::too_many_arguments)]
pub fn track(
    prior: &Mat3,
    observed: &[Vec3],
    catalog: &Catalog,
    db_star_count: usize,
    search_radius: f64,
    slack: f64,
    cfg: &SimConfig,
    cache: &mut FieldCache,
    matches: &mut Vec<StarMatch>,
) -> Option<Solution> {
    let camera = cfg.nominal_camera();
    let field = cache.field(prior, &camera, catalog, db_star_count, search_radius, slack);
    matches.clear();
    identify::match_field(observed, field, catalog, search_radius, matches);
    if matches.len() < MIN_TRACK_MATCHES {
        return None;
    }
    let attitude = identify::attitude_over(matches, observed, catalog)?;
    Some(Solution {
        matches: matches.clone(),
        attitude,
    })
}

/// Whether a tracked solution may be called CONFIDENT.
///
/// The same bars the lost-in-space self-check applies, minus the one that asks
/// for four pyramid stars: a tracked solution has no pyramid, because its
/// confirmation is the prior it came from, which was itself checked. Every
/// other clause still holds, including that each claim must land on its own
/// centroid.
fn tracked_confidence(check: &SelfCheck, claims: usize, cfg: &SimConfig) -> bool {
    check.matched >= cfg.verify_min_matches
        && check.residual_rms_px <= cfg.verify_max_rms_px
        && check.claims_explained == claims
}

/// One frame of a tracking sequence.
#[derive(Clone, Copy, Debug)]
pub struct TrackStep {
    /// Index within the sequence; frame 0 is the lost-in-space acquisition.
    pub index: usize,
    /// How the tracked solution ended.
    pub outcome: Outcome,
    /// Total attitude error, radians. Meaningless without a solution.
    pub total_error: f64,
    /// Catalogue ids claimed, and how many were right.
    pub claimed: usize,
    /// How many claims were right.
    pub correct: usize,
    /// Centroids the self-check matched.
    pub matched: usize,
    /// Blobs the centroider found.
    pub detected: usize,
    /// Nanoseconds spent identifying this frame, tracked.
    pub track_ns: u64,
    /// Nanoseconds the lost-in-space search took on the very same directions.
    pub lost_ns: u64,
    /// Whether that lost-in-space search found anything.
    pub lost_solved: bool,
    /// True when the prior failed and the frame had to be acquired afresh.
    pub reacquired: bool,
    /// Whether the rate was known, so the prior could be propagated forward.
    pub aided: bool,
    /// Estimated attitude, scalar-first quaternion with `w >= 0`.
    pub quaternion: [f64; 4],
    /// The attitude the simulator actually used, same convention.
    pub truth_quaternion: [f64; 4],
}

/// What a sequence did, reduced to the figures worth reporting.
#[derive(Clone, Copy, Debug, Default)]
pub struct TrackReport {
    /// Frames in the sequence.
    pub frames: usize,
    /// Frames whose tracked solution came out CORRECT.
    pub correct: usize,
    /// Frames that were CONFIDENT but wrong. The number that must stay zero.
    pub wrong_confident: usize,
    /// Frames where the prior failed and a full search was needed.
    pub reacquisitions: usize,
    /// Median identification time while tracking, nanoseconds.
    pub track_ns: u64,
    /// Median identification time searching the whole sky, nanoseconds.
    pub lost_ns: u64,
    /// Median total attitude error, radians.
    pub median_error: f64,
}

impl TrackReport {
    /// How many times faster tracking identified a frame. Zero if untimed.
    pub fn speedup(&self) -> f64 {
        if self.track_ns == 0 {
            return 0.0;
        }
        self.lost_ns as f64 / self.track_ns as f64
    }
}

/// Runs a slew and tracks through it, timing both paths on every frame.
///
/// Frame 0 is acquired lost-in-space, as a real tracker must. From then on the
/// previous solution seeds the next frame, and the full search still runs
/// alongside purely to time it: that comparison is the point of the mode, so
/// the figures come from the same directions rather than from separate runs.
pub fn run_sequence(
    base_seed: u64,
    cfg: &SimConfig,
    track_cfg: &TrackConfig,
    solver: &Solver<'_>,
    clock: &dyn Clock,
    ws: &mut Workspaces,
) -> Vec<TrackStep> {
    let step = track_cfg.step();
    let radius = track_cfg.aided_radius(cfg);
    let camera = cfg.nominal_camera();

    let mut steps = Vec::with_capacity(track_cfg.frames);
    let mut last: Option<Mat3> = None;
    let mut before: Option<Mat3> = None;
    let mut truth = Mat3::identity();
    let mut matches = Vec::new();
    let mut cache = FieldCache::new();
    let slack = track_cfg.field_slack_deg.to_radians();

    for index in 0..track_cfg.frames {
        let seed = trial_seed(base_seed, index);
        let frame = if index == 0 {
            let frame = simulate::simulate(seed, cfg, solver.catalog, &mut ws.simulate);
            truth = frame.a_true;
            frame
        } else {
            truth = step * truth;
            simulate::simulate_at(truth, seed, cfg, solver.catalog, &mut ws.simulate)
        };

        let detections = centroid::detect(
            &frame.image,
            frame.width,
            frame.height,
            cfg,
            &mut ws.centroid,
        );
        let detected = detections.len();
        ws.centroids.clear();
        ws.centroids.extend_from_slice(detections);
        verify::directions(&ws.centroids, &camera, &mut ws.directions);
        let identified_from = centroid::brightest(&ws.centroids, cfg).len();

        // The full search, always, for the timing comparison.
        let started = clock.now_ns();
        let lost = identify::identify(
            &ws.directions[..identified_from],
            solver.catalog,
            solver.db,
            cfg,
            clock,
            &mut ws.identify,
        );
        let lost_ns = clock.now_ns().saturating_sub(started);

        // Where this frame is expected to be. Two past attitudes give the
        // rotation between them, and applying it once more predicts this one.
        //
        // One past attitude is not enough. The window would have to cover the
        // entire slew, and at any appreciable rate such a window is ambiguous:
        // every observation has a rival inside it, so the matches are thrown
        // away and what survives is worse than nothing. A rate needs two
        // samples, so the second frame of a sequence is searched in full like
        // the first, and tracking starts at the third.
        let (prior, aided) = match (before, last) {
            (Some(before), Some(last)) => (Some(last * before.transpose() * last), true),
            _ => (None, false),
        };

        let started = clock.now_ns();
        let tracked = prior.and_then(|prior| {
            track(
                &prior,
                &ws.directions,
                solver.catalog,
                solver.db.star_count(),
                radius,
                slack,
                cfg,
                &mut cache,
                &mut matches,
            )
        });
        let track_ns = clock.now_ns().saturating_sub(started);

        // The tracked solution first, the full search behind it. A prediction
        // that produced matches can still fail to calibrate, and when it does
        // the frame must fall back rather than be given up: a tracked solution
        // that is merely non-empty is not a solved frame.
        let mut solution = None;
        let mut calibrated = None;
        let mut used_prediction = false;
        for (at, candidate) in [tracked, lost.solution.clone()].into_iter().enumerate() {
            let Some(candidate) = candidate else {
                continue;
            };
            let fit = crate::attitude::solve_calibrated(
                &candidate.matches,
                &ws.centroids,
                solver.catalog,
                solver.db.star_count(),
                cfg,
                &mut ws.identify,
            );
            if let Some(fit) = fit {
                used_prediction = at == 0;
                solution = Some(candidate);
                calibrated = Some(fit);
                break;
            }
        }
        // A frame counts as reacquired when a prediction was available but the
        // full search is what actually solved it.
        let reacquired = prior.is_some() && !used_prediction;

        let (outcome, total_error, claimed, correct, matched) = match &calibrated {
            Some(calibrated) => {
                let mut check = verify::self_check(
                    &calibrated.attitude,
                    &calibrated.matches,
                    &ws.centroids,
                    solver.catalog,
                    &calibrated.camera,
                    solver.db.star_count(),
                    cfg,
                    &mut ws.verify,
                );
                check.confident = tracked_confidence(&check, calibrated.matches.len(), cfg);
                let truth_check = verify::truth_check(
                    &calibrated.attitude,
                    &calibrated.matches,
                    &ws.centroids,
                    &frame,
                    solver.ids,
                    solver.db.star_count(),
                    cfg,
                );
                let outcome =
                    verify::classify(solution.as_ref(), Some(&check), Some(&truth_check), cfg);
                (
                    outcome,
                    truth_check.total_error,
                    truth_check.claimed,
                    truth_check.correct,
                    check.matched,
                )
            }
            None => (Outcome::NoSolution, f64::NAN, 0, 0, 0),
        };

        // Only a solution this frame's own check trusted may seed the next one.
        match (&calibrated, outcome) {
            (Some(calibrated), Outcome::Correct | Outcome::WrongConfident) => {
                // Consecutive frames are one interval apart however each was
                // solved, so a frame that had to be reacquired still pairs with
                // the one before it for the rate estimate. Only a frame that
                // produced no attitude at all breaks the chain, which the other
                // arm handles.
                before = last;
                last = Some(calibrated.attitude);
            }
            _ => {
                before = None;
                last = None;
            }
        }

        steps.push(TrackStep {
            index,
            quaternion: calibrated
                .as_ref()
                .map_or([1.0, 0.0, 0.0, 0.0], |c| matrix_to_quat(&c.attitude)),
            truth_quaternion: matrix_to_quat(&truth),
            outcome,
            total_error,
            claimed,
            correct,
            matched,
            detected,
            track_ns,
            lost_ns,
            lost_solved: lost.solution.is_some(),
            reacquired,
            aided,
        });
    }
    steps
}

/// Reduces a sequence to its report.
pub fn summarise(steps: &[TrackStep]) -> TrackReport {
    let median = |mut values: Vec<u64>| -> u64 {
        if values.is_empty() {
            return 0;
        }
        values.sort_unstable();
        values[values.len() / 2]
    };

    // Frame 0 has no prior, so its timing belongs to acquisition, not tracking.
    let tracked: Vec<u64> = steps
        .iter()
        .filter(|step| step.index > 0 && !step.reacquired)
        .map(|step| step.track_ns)
        .collect();
    let lost: Vec<u64> = steps.iter().map(|step| step.lost_ns).collect();

    let mut errors: Vec<f64> = steps
        .iter()
        .filter(|step| step.total_error.is_finite())
        .map(|step| step.total_error)
        .collect();
    errors.sort_by(f64::total_cmp);

    TrackReport {
        frames: steps.len(),
        correct: steps
            .iter()
            .filter(|step| step.outcome == Outcome::Correct)
            .count(),
        wrong_confident: steps
            .iter()
            .filter(|step| step.outcome == Outcome::WrongConfident)
            .count(),
        reacquisitions: steps.iter().filter(|step| step.reacquired).count(),
        track_ns: median(tracked),
        lost_ns: median(lost),
        median_error: errors.get(errors.len() / 2).copied().unwrap_or(f64::NAN),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::Catalog;
    use crate::math::angle_between;
    use crate::pairdb::PairDb;
    use crate::simulate::Preset;
    use crate::{NoClock, StdClock};

    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");
    const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

    fn setup(preset: Preset) -> (Catalog, PairDb, crate::catalog::IdIndex, SimConfig) {
        let catalog = Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode");
        let cfg = SimConfig::preset(preset);
        let db = PairDb::build(&catalog, &cfg).expect("database builds");
        let ids = catalog.id_index();
        (catalog, db, ids, cfg)
    }

    #[test]
    fn the_step_follows_the_configured_rate() {
        let cfg = TrackConfig {
            rate_deg_s: 2.0,
            interval_s: 0.5,
            ..TrackConfig::default()
        };
        // One degree per frame, about +x.
        assert!((cfg.step_rad() - 1.0_f64.to_radians()).abs() < 1e-15);

        let step = cfg.step();
        let turned = angle_between(
            &(step * Vec3::new(0.0, 0.0, 1.0)),
            &Vec3::new(0.0, 0.0, 1.0),
        );
        assert!((turned - cfg.step_rad()).abs() < 1e-12, "turned {turned}");
        // A rotation about +x leaves +x alone.
        let axis = step * Vec3::new(1.0, 0.0, 0.0);
        assert!(angle_between(&axis, &Vec3::new(1.0, 0.0, 0.0)) < 1e-12);
    }

    /// A slew the tracker should hold from end to end.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "a full slew, run with --release")]
    fn tracking_holds_a_slew() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let track_cfg = TrackConfig {
            rate_deg_s: 5.0,
            ..TrackConfig::default()
        };
        let mut ws = Workspaces::new();
        let steps = run_sequence(42, &cfg, &track_cfg, &solver, &NoClock, &mut ws);
        let report = summarise(&steps);

        assert_eq!(report.frames, track_cfg.frames);
        assert_eq!(report.wrong_confident, 0, "{report:?}");
        assert_eq!(report.correct, track_cfg.frames, "{report:?}");
        assert_eq!(report.reacquisitions, 0, "{report:?}");
        // Measured 0.73"; the bound leaves room for the tail without being
        // loose enough to pass if tracking stopped working.
        assert!(
            report.median_error / ARCSEC <= 5.0,
            "median error {:.2}\"",
            report.median_error / ARCSEC
        );
    }

    /// The point of the mode: a predicted window is cheaper than the sky.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "timing, run with --release")]
    fn tracking_identifies_faster_than_the_full_search() {
        let clock = StdClock::default();
        for (preset, factor) in [(Preset::Nominal, 1.0), (Preset::Hard, 3.0)] {
            let (catalog, db, ids, cfg) = setup(preset);
            let solver = Solver::new(&catalog, &db, &ids);
            let mut ws = Workspaces::new();
            let steps = run_sequence(42, &cfg, &TrackConfig::default(), &solver, &clock, &mut ws);
            let report = summarise(&steps);
            // Measured 2.4-3.1x on nominal and 13.8x on hard, where the full
            // search has to widen its tolerance; the factors asserted are well
            // inside that so a loaded machine cannot make this flap.
            assert!(
                (report.track_ns as f64) * factor < report.lost_ns as f64,
                "{preset:?}: tracked {} ns against {} ns for the full search",
                report.track_ns,
                report.lost_ns,
            );
        }
    }

    /// Rate aiding is what makes a fast slew survivable.
    ///
    /// At 60 deg/s the field moves six degrees between frames. Searching that
    /// whole span would find a rival for every star and match nothing, so this
    /// only passes because the prior is propagated before it is used.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "a full slew, run with --release")]
    fn tracking_survives_a_fast_slew() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let track_cfg = TrackConfig {
            rate_deg_s: 60.0,
            ..TrackConfig::default()
        };
        let mut ws = Workspaces::new();
        let steps = run_sequence(7, &cfg, &track_cfg, &solver, &NoClock, &mut ws);
        let report = summarise(&steps);
        assert_eq!(report.wrong_confident, 0, "{report:?}");
        assert_eq!(report.correct, track_cfg.frames, "{report:?}");
        // The window never widens with the rate; it is the prediction that
        // absorbs it.
        assert!(track_cfg.aided_radius(&cfg) < track_cfg.step_rad() / 10.0);
    }

    /// The field is rebuilt only once the boresight has left the slack.
    #[test]
    fn the_field_outlives_small_motions() {
        let (catalog, db, _ids, cfg) = setup(Preset::Nominal);
        let camera = cfg.nominal_camera();
        let slack = 1.0_f64.to_radians();
        let tolerance = identify::Tolerance::from_config(&cfg).vector;

        let mut cache = FieldCache::new();
        let attitude = Mat3::identity();
        let first = cache
            .field(
                &attitude,
                &camera,
                &catalog,
                db.star_count(),
                tolerance,
                slack,
            )
            .len();
        assert!(first > 0, "the identity attitude should see stars");
        assert_eq!(cache.rebuilds, 1);

        // A tenth of the slack: reuse.
        let small = TrackConfig {
            rate_deg_s: 0.1,
            interval_s: 1.0,
            ..TrackConfig::default()
        }
        .step();
        let moved = small * attitude;
        cache.field(&moved, &camera, &catalog, db.star_count(), tolerance, slack);
        assert_eq!(cache.rebuilds, 1, "a 0.1 deg move should not rebuild");

        // Twice the slack: rebuild.
        let far = TrackConfig {
            rate_deg_s: 2.0,
            interval_s: 1.0,
            ..TrackConfig::default()
        }
        .step();
        let gone = far * attitude;
        cache.field(&gone, &camera, &catalog, db.star_count(), tolerance, slack);
        assert_eq!(cache.rebuilds, 2, "a 2 deg move should rebuild");
    }

    /// A prior pointing somewhere else must fail rather than invent a match.
    #[test]
    fn a_wrong_prior_matches_nothing() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let mut ws = Workspaces::new();

        // One honest frame, so the directions are real.
        let steps = run_sequence(
            42,
            &cfg,
            &TrackConfig {
                frames: 1,
                ..TrackConfig::default()
            },
            &solver,
            &NoClock,
            &mut ws,
        );
        assert_eq!(steps.len(), 1);
        let truth = quat_to_matrix(steps[0].truth_quaternion);

        // Thirty degrees away is a different patch of sky entirely.
        let wrong = TrackConfig {
            rate_deg_s: 30.0,
            interval_s: 1.0,
            ..TrackConfig::default()
        }
        .step()
            * truth;
        let mut cache = FieldCache::new();
        let mut matches = Vec::new();
        let solution = track(
            &wrong,
            &ws.directions,
            &catalog,
            db.star_count(),
            TrackConfig::default().aided_radius(&cfg),
            1.0_f64.to_radians(),
            &cfg,
            &mut cache,
            &mut matches,
        );
        assert!(
            solution.is_none(),
            "a 30 deg error should not produce a solution"
        );
    }
}
