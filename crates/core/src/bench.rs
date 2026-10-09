//! Benchmark harness: run one trial, or a thousand, and summarise them.
//!
//! A trial is the whole pipeline on one seed -- simulate, centroid, identify,
//! attitude, verify -- reduced to a [`TrialResult`]. [`run_benchmark`] rolls a
//! run of them into a [`BenchReport`].
//!
//! Per-trial seeds are `splitmix64(base_seed ^ index)`, as docs/SPEC.md
//! specifies, so any single trial can be reproduced on its own from the base
//! seed and its index -- which is what lets the web UI reopen a failure.

use crate::attitude::solve_calibrated;
use crate::catalog::{Catalog, IdIndex};
use crate::centroid::{self, Centroid};
use crate::identify::{self, Diagnostics};
use crate::math::{Mat3, Vec3, matrix_to_quat, splitmix64};
use crate::pairdb::PairDb;
use crate::simulate::{self, SimConfig};
use crate::verify::{self, Outcome, SelfCheck, TruthCheck};
use crate::{Clock, Error, Result};

/// Radians per arcsecond, for reporting.
const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// Everything a trial needs that does not change between trials.
pub struct Solver<'a> {
    /// The full catalogue, down to the simulator's limiting magnitude.
    pub catalog: &'a Catalog,
    /// The pair database, covering the catalogue's magnitude-limited prefix.
    pub db: &'a PairDb,
    /// Catalogue id to index, for recall statistics.
    pub ids: &'a IdIndex,
}

impl<'a> Solver<'a> {
    /// Bundles the three, checking they belong together.
    pub fn new(catalog: &'a Catalog, db: &'a PairDb, ids: &'a IdIndex) -> Self {
        Self { catalog, db, ids }
    }
}

/// Scratch buffers for every stage, reused across trials.
#[derive(Debug, Default)]
pub struct Workspaces {
    // Crate-visible rather than private: `track::run_sequence` drives the same
    // pipeline over a slew and must reuse these buffers rather than allocate
    // its own per frame.
    pub(crate) simulate: simulate::Workspace,
    pub(crate) centroid: centroid::Workspace,
    pub(crate) identify: identify::Workspace,
    pub(crate) verify: verify::Workspace,
    pub(crate) directions: Vec<Vec3>,
    pub(crate) centroids: Vec<Centroid>,
}

impl Workspaces {
    /// Empty buffers; they size themselves on first use.
    pub fn new() -> Self {
        Self::default()
    }
}

/// Nanoseconds spent in each stage of one trial.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StageTimings {
    /// Rendering the frame. Excluded from the per-frame solve budget.
    pub simulate_ns: u64,
    /// Background, thresholding, blobs and centroids.
    pub centroid_ns: u64,
    /// Pyramid identification, including its own sweep.
    pub identify_ns: u64,
    /// The weighted Wahba solution.
    pub attitude_ns: u64,
    /// Both verification layers.
    pub verify_ns: u64,
}

impl StageTimings {
    /// Everything except rendering: the figure the 5 ms budget applies to.
    pub fn solve_ns(&self) -> u64 {
        self.centroid_ns + self.identify_ns + self.attitude_ns + self.verify_ns
    }

    /// Every stage, rendering included.
    pub fn total_ns(&self) -> u64 {
        self.simulate_ns + self.solve_ns()
    }
}

/// One trial, reduced to numbers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrialResult {
    /// The derived seed this trial ran on; re-running it reproduces the frame.
    pub seed: u64,
    /// Index within the run.
    pub index: usize,
    /// How it ended.
    pub outcome: Outcome,
    /// Total attitude error, radians. Meaningless without a solution.
    pub total_error: f64,
    /// Cross-boresight error, radians.
    pub cross_error: f64,
    /// Roll error, radians.
    pub roll_error: f64,
    /// Catalogue ids the solver claimed.
    pub claimed: usize,
    /// How many of those were right.
    pub correct: usize,
    /// Rendered truth stars the database holds: recall's denominator.
    pub available: usize,
    /// Blobs the centroider found.
    pub detected: usize,
    /// Centroids handed to the identifier.
    pub used: usize,
    /// What the self-check saw.
    pub matched: usize,
    /// Reprojection residual RMS, pixels.
    pub residual_rms_px: f64,
    /// Estimated attitude as a quaternion, scalar-first with `w >= 0`.
    pub quaternion: [f64; 4],
    /// Fitted focal length over the nominal one. One when nothing was fitted.
    pub focal_scale: f64,
    /// Identification counters.
    pub diagnostics: Diagnostics,
    /// Per-stage timings.
    pub timings: StageTimings,
}

impl TrialResult {
    /// Whether a solution came back at all.
    pub fn solved(&self) -> bool {
        self.outcome != Outcome::NoSolution
    }

    /// Total error in arcseconds, for reporting.
    pub fn total_error_arcsec(&self) -> f64 {
        self.total_error / ARCSEC
    }
}

/// Runs one trial end to end.
///
/// `seed` is the fully derived per-trial seed. The `clock` times the stages;
/// pass `NoClock` to skip timing.
pub fn run_trial(
    seed: u64,
    index: usize,
    cfg: &SimConfig,
    solver: &Solver<'_>,
    clock: &dyn Clock,
    ws: &mut Workspaces,
) -> TrialResult {
    let mut timings = StageTimings::default();

    let started = clock.now_ns();
    let frame = simulate::simulate(seed, cfg, solver.catalog, &mut ws.simulate);
    timings.simulate_ns = clock.now_ns().saturating_sub(started);

    // Centroid once, and keep every detection.
    //
    // Identification sees only the brightest few, because its cost is
    // combinatorial in the number of stars. Everything downstream -- the focal
    // fit, the sweep that follows it, and both checks -- sees all of them,
    // because those are linear and more stars is strictly better: the roll
    // angle is the weakly constrained one in a narrow field and it improves
    // with both the count and the spread. The brightest are a prefix of the
    // flux-sorted list, so an index into them is already an index into the
    // whole list and the pyramid's indices stay valid.
    let started = clock.now_ns();
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
    verify::directions(&ws.centroids, &cfg.nominal_camera(), &mut ws.directions);
    let identified_from = centroid::brightest(&ws.centroids, cfg).len();
    timings.centroid_ns = clock.now_ns().saturating_sub(started);

    let identification = identify::identify(
        &ws.directions[..identified_from],
        solver.catalog,
        solver.db,
        cfg,
        clock,
        &mut ws.identify,
    );
    timings.identify_ns = identification.diagnostics.elapsed_ns;

    // The attitude stage proper: fit the focal length, then re-solve with
    // brightness weights through the camera that fit implies.
    let started = clock.now_ns();
    let calibrated = identification.solution.as_ref().and_then(|solution| {
        solve_calibrated(
            &solution.matches,
            &ws.centroids,
            solver.catalog,
            solver.db.star_count(),
            cfg,
            &mut ws.identify,
        )
    });
    timings.attitude_ns = clock.now_ns().saturating_sub(started);

    let started = clock.now_ns();
    let (self_check, truth) = match (&identification.solution, &calibrated) {
        (Some(_), Some(calibrated)) => (
            Some(verify::self_check(
                &calibrated.attitude,
                &calibrated.matches,
                &ws.centroids,
                solver.catalog,
                &calibrated.camera,
                solver.db.star_count(),
                cfg,
                &mut ws.verify,
            )),
            Some(verify::truth_check(
                &calibrated.attitude,
                &calibrated.matches,
                &ws.centroids,
                &frame,
                solver.ids,
                solver.db.star_count(),
                cfg,
            )),
        ),
        _ => (None, None),
    };
    let outcome = verify::classify(
        identification.solution.as_ref(),
        self_check.as_ref(),
        truth.as_ref(),
        cfg,
    );
    timings.verify_ns = clock.now_ns().saturating_sub(started);

    let blank_self = SelfCheck {
        matched: 0,
        residual_rms_px: f64::INFINITY,
        pyramid_stars: 0,
        claims_explained: 0,
        claims: 0,
        confident: false,
    };
    let blank_truth = TruthCheck {
        claimed: 0,
        correct: 0,
        available: 0,
        total_error: f64::NAN,
        cross_error: f64::NAN,
        roll_error: f64::NAN,
    };
    let self_check = self_check.unwrap_or(blank_self);
    let truth = truth.unwrap_or(blank_truth);

    TrialResult {
        seed,
        index,
        outcome,
        total_error: truth.total_error,
        cross_error: truth.cross_error,
        roll_error: truth.roll_error,
        claimed: truth.claimed,
        correct: truth.correct,
        available: truth.available,
        detected,
        used: identified_from,
        matched: self_check.matched,
        residual_rms_px: self_check.residual_rms_px,
        quaternion: calibrated
            .as_ref()
            .map_or([1.0, 0.0, 0.0, 0.0], |c| matrix_to_quat(&c.attitude)),
        focal_scale: calibrated.as_ref().map_or(1.0, |c| c.focal_scale),
        diagnostics: identification.diagnostics,
        timings,
    }
}

/// The per-trial seed for a run: `splitmix64(base_seed ^ index)`.
pub fn trial_seed(base_seed: u64, index: usize) -> u64 {
    splitmix64(base_seed ^ index as u64)
}

/// Median, 95th percentile and maximum of a sample.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Spread {
    /// Middle value.
    pub median: f64,
    /// 95th percentile.
    pub p95: f64,
    /// Largest value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// How many values went in.
    pub count: usize,
}

impl Spread {
    /// Summarises `values`, which is sorted in place. Empty input gives zeros.
    pub fn of(values: &mut [f64]) -> Self {
        if values.is_empty() {
            return Self::default();
        }
        values.sort_by(f64::total_cmp);
        let count = values.len();
        // The p95 index is clamped, so small samples report their largest.
        let p95_at = ((count as f64 * 0.95).ceil() as usize).min(count) - 1;
        Self {
            median: values[count / 2],
            p95: values[p95_at],
            max: values[count - 1],
            mean: values.iter().sum::<f64>() / count as f64,
            count,
        }
    }

    /// The same spread with every value divided by `divisor`.
    pub fn scaled(self, divisor: f64) -> Self {
        Self {
            median: self.median / divisor,
            p95: self.p95 / divisor,
            max: self.max / divisor,
            mean: self.mean / divisor,
            count: self.count,
        }
    }
}

/// Mean and 95th percentile for one stage, in nanoseconds.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StageSpread {
    /// Rendering.
    pub simulate: Spread,
    /// Centroiding, including the conversion to directions.
    pub centroid: Spread,
    /// Identification.
    pub identify: Spread,
    /// The weighted attitude solution.
    pub attitude: Spread,
    /// Verification.
    pub verify: Spread,
    /// Everything but rendering.
    pub solve: Spread,
}

/// A whole run, summarised.
#[derive(Clone, Debug, PartialEq)]
pub struct BenchReport {
    /// Base seed the run was derived from.
    pub base_seed: u64,
    /// How many trials ran.
    pub trials: usize,
    /// Trials per outcome, in [`Outcome::ALL`] order.
    pub counts: [usize; 4],
    /// Percentage of trials that were CORRECT: the Score.
    pub score: f64,
    /// Correct ids over claimed ids, across every trial that claimed any.
    pub id_precision: f64,
    /// Correct ids over the rendered truth stars the database holds.
    pub id_recall: f64,
    /// Total attitude error over solved trials, arcseconds.
    pub total_error_arcsec: Spread,
    /// Cross-boresight error, arcseconds.
    pub cross_error_arcsec: Spread,
    /// Roll error, arcseconds.
    pub roll_error_arcsec: Spread,
    /// Per-stage timings, nanoseconds.
    pub timings: StageSpread,
}

impl BenchReport {
    /// Trials with the given outcome.
    pub fn count(&self, outcome: Outcome) -> usize {
        self.counts[Outcome::ALL.iter().position(|o| *o == outcome).unwrap_or(0)]
    }

    /// Percentage of trials with the given outcome.
    pub fn percent(&self, outcome: Outcome) -> f64 {
        if self.trials == 0 {
            return 0.0;
        }
        100.0 * self.count(outcome) as f64 / self.trials as f64
    }

    /// The count that must be zero for a run to be acceptable.
    pub fn wrong_confident(&self) -> usize {
        self.count(Outcome::WrongConfident)
    }
}

/// Runs `trials` trials and summarises them, keeping every per-trial row.
///
/// Returns the rows alongside the report so a caller can write a CSV or open
/// one failing seed.
pub fn run_benchmark(
    base_seed: u64,
    trials: usize,
    cfg: &SimConfig,
    solver: &Solver<'_>,
    clock: &dyn Clock,
) -> Result<(BenchReport, Vec<TrialResult>)> {
    if trials == 0 {
        return Err(Error::EmptyBenchmark);
    }
    let mut ws = Workspaces::new();
    let mut rows = Vec::with_capacity(trials);
    for index in 0..trials {
        rows.push(run_trial(
            trial_seed(base_seed, index),
            index,
            cfg,
            solver,
            clock,
            &mut ws,
        ));
    }
    Ok((summarise(base_seed, &rows), rows))
}

/// Rolls per-trial rows into a report.
pub fn summarise(base_seed: u64, rows: &[TrialResult]) -> BenchReport {
    let trials = rows.len();
    let mut counts = [0usize; 4];
    for row in rows {
        if let Some(at) = Outcome::ALL.iter().position(|o| *o == row.outcome) {
            counts[at] += 1;
        }
    }

    let claimed: usize = rows.iter().map(|r| r.claimed).sum();
    let correct: usize = rows.iter().map(|r| r.correct).sum();
    let available: usize = rows.iter().map(|r| r.available).sum();

    // Errors are only meaningful where a solution exists.
    let mut total: Vec<f64> = rows
        .iter()
        .filter(|r| r.solved())
        .map(|r| r.total_error)
        .collect();
    let mut cross: Vec<f64> = rows
        .iter()
        .filter(|r| r.solved())
        .map(|r| r.cross_error)
        .collect();
    let mut roll: Vec<f64> = rows
        .iter()
        .filter(|r| r.solved())
        .map(|r| r.roll_error)
        .collect();

    let stage = |pick: fn(&StageTimings) -> u64| -> Spread {
        let mut values: Vec<f64> = rows.iter().map(|r| pick(&r.timings) as f64).collect();
        Spread::of(&mut values)
    };

    BenchReport {
        base_seed,
        trials,
        counts,
        score: if trials == 0 {
            0.0
        } else {
            100.0 * counts[0] as f64 / trials as f64
        },
        id_precision: if claimed == 0 {
            0.0
        } else {
            correct as f64 / claimed as f64
        },
        id_recall: if available == 0 {
            0.0
        } else {
            correct as f64 / available as f64
        },
        total_error_arcsec: Spread::of(&mut total).scaled(ARCSEC),
        cross_error_arcsec: Spread::of(&mut cross).scaled(ARCSEC),
        roll_error_arcsec: Spread::of(&mut roll).scaled(ARCSEC),
        timings: StageSpread {
            simulate: stage(|t| t.simulate_ns),
            centroid: stage(|t| t.centroid_ns),
            identify: stage(|t| t.identify_ns),
            attitude: stage(|t| t.attitude_ns),
            verify: stage(|t| t.verify_ns),
            solve: stage(StageTimings::solve_ns),
        },
    }
}

/// The CSV header matching [`csv_row`].
pub fn csv_header() -> &'static str {
    "index,seed,outcome,total_arcsec,cross_arcsec,roll_arcsec,claimed,correct,available,\
detected,used,matched,residual_rms_px,focal_scale,qw,qx,qy,qz,tries,pair_queries,\
simulate_ns,centroid_ns,identify_ns,attitude_ns,verify_ns"
}

/// One trial as a CSV row, in the order [`csv_header`] names.
pub fn csv_row(row: &TrialResult) -> String {
    let [qw, qx, qy, qz] = row.quaternion;
    format!(
        "{},{},{},{:.6},{:.6},{:.6},{},{},{},{},{},{},{:.6},\
{:.9},{:.9},{:.9},{:.9},{:.9},{},{},{},{},{},{},{}",
        row.index,
        row.seed,
        row.outcome,
        row.total_error / ARCSEC,
        row.cross_error / ARCSEC,
        row.roll_error / ARCSEC,
        row.claimed,
        row.correct,
        row.available,
        row.detected,
        row.used,
        row.matched,
        row.residual_rms_px,
        row.focal_scale,
        qw,
        qx,
        qy,
        qz,
        row.diagnostics.tries,
        row.diagnostics.pair_queries,
        row.timings.simulate_ns,
        row.timings.centroid_ns,
        row.timings.identify_ns,
        row.timings.attitude_ns,
        row.timings.verify_ns,
    )
}

/// The report as a human-readable block, as the CLI prints it.
pub fn format_report(report: &BenchReport, cfg: &SimConfig) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    let ms = |spread: Spread| (spread.mean / 1e6, spread.p95 / 1e6);

    let _ = writeln!(
        out,
        "trials {}  base seed {}  db cut V{:.1}",
        report.trials, report.base_seed, cfg.db_mag_cut
    );
    let _ = writeln!(out);
    let _ = writeln!(out, "  Score  {:.1}%   (CORRECT / total)", report.score);
    let _ = writeln!(
        out,
        "  WRONG_CONFIDENT  {}{}",
        report.wrong_confident(),
        if report.wrong_confident() == 0 {
            ""
        } else {
            "   <-- must be zero"
        }
    );
    let _ = writeln!(out);
    for outcome in Outcome::ALL {
        let _ = writeln!(
            out,
            "  {:<16} {:>6}  {:>6.2}%",
            outcome.as_str(),
            report.count(outcome),
            report.percent(outcome)
        );
    }
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "  ID precision {:.4}   recall {:.4}",
        report.id_precision, report.id_recall
    );
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "  {:<18} {:>10} {:>10} {:>10}",
        "error (arcsec)", "median", "p95", "max"
    );
    for (label, spread) in [
        ("total", report.total_error_arcsec),
        ("cross-boresight", report.cross_error_arcsec),
        ("roll", report.roll_error_arcsec),
    ] {
        let _ = writeln!(
            out,
            "  {:<18} {:>10.2} {:>10.2} {:>10.2}",
            label, spread.median, spread.p95, spread.max
        );
    }
    let _ = writeln!(out);
    let _ = writeln!(out, "  {:<18} {:>10} {:>10}", "stage (ms)", "mean", "p95");
    for (label, spread) in [
        ("simulate", report.timings.simulate),
        ("centroid", report.timings.centroid),
        ("identify", report.timings.identify),
        ("attitude", report.timings.attitude),
        ("verify", report.timings.verify),
        ("solve (no render)", report.timings.solve),
    ] {
        let (mean, p95) = ms(spread);
        let _ = writeln!(out, "  {:<18} {:>10.3} {:>10.3}", label, mean, p95);
    }
    out
}

/// The estimated boresight as right ascension, declination and roll, radians.
///
/// Reported by the CLI and the web UI next to the quaternion. `A` maps
/// inertial to camera, so the boresight's inertial direction is its third row,
/// and roll is the rotation about it.
pub fn boresight_radec_roll(attitude: &Mat3) -> (f64, f64, f64) {
    let boresight = attitude.row(2).transpose();
    let (ra, dec) = crate::math::unit_to_radec(&boresight);
    // Roll: where the camera's +x axis lands relative to the local north.
    let roll = attitude[(1, 2)].atan2(attitude[(0, 2)]);
    (ra, dec, roll)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NoClock;
    use crate::simulate::Preset;

    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

    fn catalog() -> Catalog {
        Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode")
    }

    /// Exact equality between two trial results, NaN included.
    ///
    /// An unsolved trial carries NaN error fields -- honestly, since the
    /// summary filters on `solved()` -- and `PartialEq` says NaN differs from
    /// itself, so it cannot express "these two runs produced the same thing".
    /// Comparing bit patterns can.
    fn identical(a: &TrialResult, b: &TrialResult) -> bool {
        let bits = |v: f64| v.to_bits();
        a.seed == b.seed
            && a.index == b.index
            && a.outcome == b.outcome
            && bits(a.total_error) == bits(b.total_error)
            && bits(a.cross_error) == bits(b.cross_error)
            && bits(a.roll_error) == bits(b.roll_error)
            && bits(a.residual_rms_px) == bits(b.residual_rms_px)
            && bits(a.focal_scale) == bits(b.focal_scale)
            && a.quaternion.map(bits) == b.quaternion.map(bits)
            && (
                a.claimed,
                a.correct,
                a.available,
                a.detected,
                a.used,
                a.matched,
            ) == (
                b.claimed,
                b.correct,
                b.available,
                b.detected,
                b.used,
                b.matched,
            )
            && a.diagnostics == b.diagnostics
    }

    /// Catalogue, database and id index for a preset.
    fn setup(preset: Preset) -> (Catalog, PairDb, IdIndex, SimConfig) {
        let catalog = catalog();
        let cfg = SimConfig::preset(preset);
        let db = PairDb::build(&catalog, &cfg).expect("build");
        let ids = catalog.id_index();
        (catalog, db, ids, cfg)
    }

    // --- seeds and reproducibility ---

    #[test]
    fn trial_seeds_follow_the_documented_rule() {
        for base in [0u64, 42, u64::MAX] {
            for index in [0usize, 1, 7, 999] {
                assert_eq!(
                    trial_seed(base, index),
                    splitmix64(base ^ index as u64),
                    "base {base} index {index}"
                );
            }
        }
        // Distinct indices give distinct seeds, which is what makes a run a
        // spread of frames rather than one frame repeated.
        let seeds: std::collections::HashSet<u64> =
            (0..1000).map(|index| trial_seed(42, index)).collect();
        assert_eq!(seeds.len(), 1000);
    }

    /// Phase 9 needs clicking a failed trial to reproduce it exactly, which
    /// means a trial taken out of a run must behave the same on its own.
    #[test]
    fn a_trial_reproduces_outside_its_run() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let (_, rows) = run_benchmark(42, 12, &cfg, &solver, &NoClock).expect("run");

        let mut ws = Workspaces::new();
        for row in &rows {
            let alone = run_trial(row.seed, row.index, &cfg, &solver, &NoClock, &mut ws);
            assert!(
                identical(&alone, row),
                "trial {} differs run-alone: {alone:?} vs {row:?}",
                row.index
            );
        }
    }

    #[test]
    fn a_reused_workspace_changes_nothing() {
        let (catalog, db, ids, cfg) = setup(Preset::Brutal);
        let solver = Solver::new(&catalog, &db, &ids);
        let seed = trial_seed(42, 3);

        let fresh = run_trial(seed, 3, &cfg, &solver, &NoClock, &mut Workspaces::new());
        let mut shared = Workspaces::new();
        for index in 100..105 {
            run_trial(
                trial_seed(42, index),
                index,
                &cfg,
                &solver,
                &NoClock,
                &mut shared,
            );
        }
        let reused = run_trial(seed, 3, &cfg, &solver, &NoClock, &mut shared);
        assert!(identical(&fresh, &reused), "{fresh:?} vs {reused:?}");
    }

    // --- one trial ---

    #[test]
    fn a_nominal_trial_comes_out_correct() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let mut ws = Workspaces::new();

        for index in 0..15usize {
            let row = run_trial(
                trial_seed(42, index),
                index,
                &cfg,
                &solver,
                &NoClock,
                &mut ws,
            );
            assert_eq!(row.outcome, Outcome::Correct, "trial {index}: {row:?}");
            assert_eq!(row.index, index);
            assert_eq!(row.seed, trial_seed(42, index));
            assert!(row.solved());
            assert_eq!(row.correct, row.claimed);
            assert!(row.claimed >= 4);
            assert!(row.detected >= row.used);
            assert!(row.used <= cfg.max_centroids);
            assert!(row.available >= row.correct);
            // The quaternion is normalised with a non-negative scalar part.
            let [w, x, y, z] = row.quaternion;
            assert!(w >= 0.0);
            assert!((w * w + x * x + y * y + z * z - 1.0).abs() < 1e-12);
            // Nothing to calibrate on this preset.
            assert!((row.focal_scale - 1.0).abs() < 1e-3);
        }
    }

    #[test]
    fn timings_are_recorded_when_a_clock_is_given() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let clock = crate::StdClock::new();
        let row = run_trial(
            trial_seed(42, 0),
            0,
            &cfg,
            &solver,
            &clock,
            &mut Workspaces::new(),
        );
        assert!(row.timings.simulate_ns > 0);
        assert!(row.timings.centroid_ns > 0);
        assert!(row.timings.identify_ns > 0);
        assert_eq!(
            row.timings.solve_ns(),
            row.timings.centroid_ns
                + row.timings.identify_ns
                + row.timings.attitude_ns
                + row.timings.verify_ns
        );
        assert_eq!(
            row.timings.total_ns(),
            row.timings.simulate_ns + row.timings.solve_ns()
        );

        // And `NoClock` leaves them at zero rather than at something wrong.
        let untimed = run_trial(
            trial_seed(42, 0),
            0,
            &cfg,
            &solver,
            &NoClock,
            &mut Workspaces::new(),
        );
        assert_eq!(untimed.timings, StageTimings::default());
    }

    // --- spreads ---

    #[test]
    fn spread_summarises_a_sample() {
        let mut values = vec![5.0, 1.0, 4.0, 2.0, 3.0];
        let spread = Spread::of(&mut values);
        assert_eq!(spread.count, 5);
        assert_eq!(spread.median, 3.0);
        assert_eq!(spread.max, 5.0);
        assert_eq!(spread.mean, 3.0);
        // Sorted in place, as documented.
        assert_eq!(values, [1.0, 2.0, 3.0, 4.0, 5.0]);

        // A single value is its own median, p95 and max.
        let single = Spread::of(&mut [7.0]);
        assert_eq!((single.median, single.p95, single.max), (7.0, 7.0, 7.0));

        // Empty input must not panic.
        assert_eq!(Spread::of(&mut []), Spread::default());

        // p95 of a hundred values is the 95th.
        let mut hundred: Vec<f64> = (1..=100).map(f64::from).collect();
        assert_eq!(Spread::of(&mut hundred).p95, 95.0);

        // Scaling divides every statistic.
        let scaled = Spread::of(&mut [2.0, 4.0, 6.0]).scaled(2.0);
        assert_eq!((scaled.median, scaled.max), (2.0, 3.0));
    }

    // --- reports ---

    #[test]
    fn a_report_accounts_for_every_trial() {
        let (catalog, db, ids, cfg) = setup(Preset::Hard);
        let solver = Solver::new(&catalog, &db, &ids);
        let (report, rows) = run_benchmark(9, 60, &cfg, &solver, &NoClock).expect("run");

        assert_eq!(report.trials, 60);
        assert_eq!(rows.len(), 60);
        assert_eq!(report.counts.iter().sum::<usize>(), 60);
        assert_eq!(report.base_seed, 9);

        // Percentages agree with the counts.
        let mut total_percent = 0.0;
        for outcome in Outcome::ALL {
            let count = report.count(outcome);
            assert_eq!(rows.iter().filter(|r| r.outcome == outcome).count(), count);
            total_percent += report.percent(outcome);
        }
        assert!((total_percent - 100.0).abs() < 1e-9);
        assert!((report.score - report.percent(Outcome::Correct)).abs() < 1e-12);
        assert_eq!(
            report.wrong_confident(),
            report.count(Outcome::WrongConfident)
        );

        // Precision and recall are fractions, and the error spreads only cover
        // the trials that produced a solution.
        assert!((0.0..=1.0).contains(&report.id_precision));
        assert!((0.0..=1.0).contains(&report.id_recall));
        let solved = rows.iter().filter(|r| r.solved()).count();
        assert_eq!(report.total_error_arcsec.count, solved);

        // Summarising the rows again gives the same report.
        assert_eq!(summarise(9, &rows), report);
    }

    #[test]
    fn a_benchmark_needs_at_least_one_trial() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        assert!(matches!(
            run_benchmark(1, 0, &cfg, &solver, &NoClock),
            Err(Error::EmptyBenchmark)
        ));
    }

    // --- output formats ---

    #[test]
    fn csv_rows_line_up_with_the_header() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let (_, rows) = run_benchmark(42, 4, &cfg, &solver, &NoClock).expect("run");

        let columns = csv_header().split(',').count();
        for row in &rows {
            let line = csv_row(row);
            assert_eq!(
                line.split(',').count(),
                columns,
                "row has the wrong field count: {line}"
            );
            assert!(!line.contains('\n'));
        }
        // The outcome goes in by name.
        assert!(csv_row(&rows[0]).contains(rows[0].outcome.as_str()));
    }

    #[test]
    fn the_report_formats_without_panicking() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let (report, _) = run_benchmark(42, 8, &cfg, &solver, &NoClock).expect("run");
        let text = format_report(&report, &cfg);
        for needle in [
            "Score",
            "WRONG_CONFIDENT",
            "CORRECT",
            "ID precision",
            "median",
            "solve (no render)",
        ] {
            assert!(text.contains(needle), "report is missing {needle:?}");
        }
    }

    #[test]
    fn boresight_angles_match_the_attitude() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let row = run_trial(
            trial_seed(42, 0),
            0,
            &cfg,
            &solver,
            &NoClock,
            &mut Workspaces::new(),
        );
        let attitude = crate::math::quat_to_matrix(row.quaternion);
        let (ra, dec, roll) = boresight_radec_roll(&attitude);
        assert!((0.0..std::f64::consts::TAU).contains(&ra));
        assert!(dec.abs() <= std::f64::consts::FRAC_PI_2);
        assert!(roll.abs() <= std::f64::consts::PI);

        // The reported direction is the boresight: the inertial vector that
        // maps onto the camera's +z.
        let boresight = crate::math::radec_to_unit(ra, dec);
        let camera_frame = attitude * boresight;
        assert!(camera_frame.z > 0.999999, "{camera_frame:?}");
    }

    // --- the golden run ---

    /// Phase 7 acceptance: seed 42, nominal, 1000 trials must give a Score of
    /// at least 99%, no WRONG_CONFIDENT at all, and a median total error within
    /// 10 arcsec. This is the run the Golden Benchmark Record tracks.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "1000 full frames, run with --release")]
    fn the_golden_run_meets_its_targets() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let (report, _) = run_benchmark(42, 1000, &cfg, &solver, &NoClock).expect("run");

        assert_eq!(
            report.wrong_confident(),
            0,
            "WRONG_CONFIDENT must be zero, got {}",
            report.wrong_confident()
        );
        assert!(
            report.score >= 99.0,
            "Score {:.2}% is under the 99% floor",
            report.score
        );
        assert!(
            report.total_error_arcsec.median <= 10.0,
            "median total error {:.2} arcsec is over the 10 arcsec ceiling",
            report.total_error_arcsec.median
        );
    }

    /// The easy preset should be flawless; anything less means something broke
    /// upstream of the difficulty knobs.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "300 full frames, run with --release")]
    fn the_easy_preset_is_flawless() {
        let (catalog, db, ids, cfg) = setup(Preset::Easy);
        let solver = Solver::new(&catalog, &db, &ids);
        let (report, _) = run_benchmark(42, 300, &cfg, &solver, &NoClock).expect("run");
        assert_eq!(report.wrong_confident(), 0);
        assert_eq!(report.count(Outcome::Correct), 300, "{report:?}");
    }

    /// Phase 7's second acceptance criterion: `hard` scores at least 95% with
    /// nothing confidently wrong.
    ///
    /// The preset renders through a camera whose focal length is 0.2% off and
    /// which has a little radial distortion, neither of which the solver models.
    /// Measured, that puts an observed pair angle out by about 8 arcseconds per
    /// degree of separation, so a constant identification tolerance of 50"
    /// admits nothing wider than 8 degrees and the pyramid, which needs a spread
    /// triangle, never closes.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "500 full frames, run with --release")]
    fn the_hard_preset_meets_its_target() {
        let (catalog, db, ids, cfg) = setup(Preset::Hard);
        let solver = Solver::new(&catalog, &db, &ids);
        const TRIALS: usize = 500;
        let (report, _) = run_benchmark(42, TRIALS, &cfg, &solver, &NoClock).expect("run");
        assert_eq!(report.wrong_confident(), 0, "{report:?}");
        let score = 100.0 * report.count(Outcome::Correct) as f64 / TRIALS as f64;
        assert!(
            score >= 95.0,
            "hard scored {score:.1}%, target 95%; {} no-solution, {} rejected",
            report.count(Outcome::NoSolution),
            report.count(Outcome::Rejected),
        );
    }

    /// The solve stages, rendering excluded, have a 5 ms budget per frame.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "timing, run with --release")]
    fn the_solve_stages_fit_their_budget() {
        let (catalog, db, ids, cfg) = setup(Preset::Nominal);
        let solver = Solver::new(&catalog, &db, &ids);
        let clock = crate::StdClock::new();
        let (report, _) = run_benchmark(42, 200, &cfg, &solver, &clock).expect("run");
        let mean_ms = report.timings.solve.mean / 1e6;
        assert!(
            mean_ms <= 5.0,
            "solve averaged {mean_ms:.3} ms against a 5 ms budget"
        );
    }
}
