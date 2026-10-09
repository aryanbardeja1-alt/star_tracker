//! wasm-bindgen wrapper around `tracker-core`.
//!
//! [`init`] loads the embedded catalogue and builds the pair database once;
//! everything else hangs off the [`Tracker`] it returns, so the 10 MB table is
//! built on start-up rather than per call.
//!
//! Configuration crosses as a JS value and may be either a preset name
//! (`"nominal"`) or a whole `SimConfig` object, which `serde-wasm-bindgen`
//! converts without going through a JSON string. Images cross as
//! `Uint16Array`, never as JSON: [`Frame::image`] hands over a copy, and
//! [`Frame::image_ptr`] with [`Frame::image_len`] let the caller build a view
//! straight over WASM memory when copying a megapixel would matter.
//!
//! Timing comes from `performance.now()`. `std::time::Instant` compiles for
//! this target and then traps at run time, which is why `tracker-core` takes a
//! [`Clock`] rather than reading one itself.

#![forbid(unsafe_code)]

use serde::Serialize;
use tracker_core::bench::{self, Solver, TrialResult, Workspaces};
use tracker_core::catalog::{Catalog, IdIndex};
use tracker_core::centroid::{self, Centroid};
use tracker_core::identify::{self, StarMatch};
use tracker_core::math::quat_to_matrix;
use tracker_core::pairdb::PairDb;
use tracker_core::simulate::{self, Preset, SimConfig, SimFrame};
use tracker_core::track::{self, TrackConfig};
use tracker_core::verify::{self, SelfCheck};
use tracker_core::{Clock, math};
use wasm_bindgen::prelude::*;

/// The committed catalogue, embedded so the module needs no second fetch.
const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = performance, js_name = now)]
    fn performance_now() -> f64;
}

/// A clock backed by `performance.now()`.
struct PerformanceClock;

impl Clock for PerformanceClock {
    fn now_ns(&self) -> u64 {
        // `performance.now()` is milliseconds as an f64; browsers clamp its
        // resolution, so a single fast stage can read as zero.
        (performance_now() * 1.0e6) as u64
    }
}

/// Crate version, so a page can confirm which build it loaded.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// SplitMix64 of `x`, re-exported from `tracker-core` as a determinism check.
#[wasm_bindgen]
pub fn splitmix64(x: u64) -> u64 {
    math::splitmix64(x)
}

/// Loads the catalogue and builds the pair database.
///
/// The returned handle owns both and should be kept for the life of the page.
#[wasm_bindgen]
pub fn init() -> Result<Tracker, JsError> {
    Tracker::new()
}

/// The solver, with its catalogue and pair database built.
#[wasm_bindgen]
pub struct Tracker {
    catalog: Catalog,
    ids: IdIndex,
    db: PairDb,
    /// The configuration the database was built for. Only `db_mag_cut` and the
    /// field of view affect it, so it is rebuilt only when those move.
    db_key: (f32, f64, usize),
    trial_ws: Workspaces,
    sim_ws: simulate::Workspace,
    centroid_ws: centroid::Workspace,
    identify_ws: identify::Workspace,
    verify_ws: verify::Workspace,
    centroids: Vec<Centroid>,
    directions: Vec<tracker_core::math::Vec3>,
}

#[wasm_bindgen]
impl Tracker {
    /// As [`init`].
    #[wasm_bindgen(constructor)]
    pub fn new() -> Result<Tracker, JsError> {
        let catalog = Catalog::from_bytes(CATALOG_BIN).map_err(to_js)?;
        let ids = catalog.id_index();
        let cfg = SimConfig::default();
        let db = PairDb::build(&catalog, &cfg).map_err(to_js)?;
        Ok(Tracker {
            catalog,
            ids,
            db,
            db_key: db_key(&cfg),
            trial_ws: Workspaces::new(),
            sim_ws: simulate::Workspace::new(),
            centroid_ws: centroid::Workspace::new(),
            identify_ws: identify::Workspace::new(),
            verify_ws: verify::Workspace::new(),
            centroids: Vec::new(),
            directions: Vec::new(),
        })
    }

    /// Catalogue entries loaded.
    #[wasm_bindgen(getter)]
    pub fn catalog_stars(&self) -> usize {
        self.catalog.stars.len()
    }

    /// Stars the pair database covers.
    #[wasm_bindgen(getter)]
    pub fn database_stars(&self) -> usize {
        self.db.star_count()
    }

    /// Pairs the database holds.
    #[wasm_bindgen(getter)]
    pub fn database_pairs(&self) -> usize {
        self.db.pairs().len()
    }

    /// The default configuration, as a JS object the UI can edit and hand back.
    #[wasm_bindgen]
    pub fn default_config(&self, preset: &str) -> Result<JsValue, JsError> {
        let preset: Preset = preset.parse().map_err(to_js)?;
        serde_wasm_bindgen::to_value(&SimConfig::preset(preset)).map_err(JsError::from)
    }

    /// Simulates one frame and returns it with its ground truth.
    ///
    /// `config` is a preset name or a `SimConfig` object.
    #[wasm_bindgen]
    pub fn simulate_frame(&mut self, seed: u64, config: JsValue) -> Result<Frame, JsError> {
        let cfg = read_config(&config)?;
        let frame = simulate::simulate(seed, &cfg, &self.catalog, &mut self.sim_ws);
        Ok(Frame::new(seed, frame))
    }

    /// Solves a frame with no access to truth, as a real tracker would.
    ///
    /// `image` is row-major `width * height` samples. Returns the centroids,
    /// the identifications, the attitude and the self-check.
    #[wasm_bindgen]
    pub fn solve_frame(&mut self, image: &[u16], config: JsValue) -> Result<JsValue, JsError> {
        let cfg = read_config(&config)?;
        self.ensure_database(&cfg)?;
        let clock = PerformanceClock;

        let detections =
            centroid::detect(image, cfg.width, cfg.height, &cfg, &mut self.centroid_ws);
        let detected = detections.len();
        // Every detection, with identification reading only the brightest
        // prefix of them; see `bench::run_trial`, which this mirrors.
        self.centroids.clear();
        self.centroids.extend_from_slice(detections);
        verify::directions(&self.centroids, &cfg.nominal_camera(), &mut self.directions);
        let identified_from = centroid::brightest(&self.centroids, &cfg).len();

        let identification = identify::identify(
            &self.directions[..identified_from],
            &self.catalog,
            &self.db,
            &cfg,
            &clock,
            &mut self.identify_ws,
        );

        let solved = identification.solution.as_ref().and_then(|solution| {
            tracker_core::attitude::solve_calibrated(
                &solution.matches,
                &self.centroids,
                &self.catalog,
                self.db.star_count(),
                &cfg,
                &mut self.identify_ws,
            )
        });

        let report = match &solved {
            Some(calibrated) => {
                let check = verify::self_check(
                    &calibrated.attitude,
                    &calibrated.matches,
                    &self.centroids,
                    &self.catalog,
                    &calibrated.camera,
                    self.db.star_count(),
                    &cfg,
                    &mut self.verify_ws,
                );
                SolveReport {
                    detected,
                    centroids: self.centroids.iter().map(CentroidDto::from).collect(),
                    matches: calibrated
                        .matches
                        .iter()
                        .map(|matched| MatchDto::new(matched, &self.catalog))
                        .collect(),
                    attitude: Some(AttitudeDto::new(
                        &calibrated.attitude,
                        calibrated.focal_scale,
                    )),
                    self_check: Some(SelfCheckDto::from(&check)),
                    identify_tries: identification.diagnostics.tries,
                    identify_ms: identification.diagnostics.elapsed_ns as f64 / 1.0e6,
                }
            }
            None => SolveReport {
                detected,
                centroids: self.centroids.iter().map(CentroidDto::from).collect(),
                matches: Vec::new(),
                attitude: None,
                self_check: None,
                identify_tries: identification.diagnostics.tries,
                identify_ms: identification.diagnostics.elapsed_ns as f64 / 1.0e6,
            },
        };
        serde_wasm_bindgen::to_value(&report).map_err(JsError::from)
    }

    /// Runs one scored trial, returning the result plus what a drawing needs.
    ///
    /// `seed` here is the *base* seed and `index` the trial number, so this
    /// reproduces exactly the trial a benchmark run at the same base seed
    /// produced: that is what lets the UI reopen a failure.
    #[wasm_bindgen]
    pub fn run_trial(
        &mut self,
        base_seed: u64,
        index: usize,
        config: JsValue,
    ) -> Result<JsValue, JsError> {
        let cfg = read_config(&config)?;
        self.ensure_database(&cfg)?;
        let seed = bench::trial_seed(base_seed, index);

        // The frame is re-simulated so the drawing data matches the trial, at
        // the cost of rendering twice; a single click can afford it.
        let row = {
            let solver = Solver::new(&self.catalog, &self.db, &self.ids);
            bench::run_trial(
                seed,
                index,
                &cfg,
                &solver,
                &PerformanceClock,
                &mut self.trial_ws,
            )
        };
        let frame = simulate::simulate(seed, &cfg, &self.catalog, &mut self.sim_ws);

        let detail = TrialDetail {
            trial: TrialRow::from(&row),
            truth: frame.truth.iter().map(TruthDto::from).collect(),
            false_stars: frame.false_stars.clone(),
            true_attitude: AttitudeDto::new(&frame.a_true, 1.0),
        };
        serde_wasm_bindgen::to_value(&detail).map_err(JsError::from)
    }

    /// Runs `count` trials starting at `start`, for a worker to accumulate.
    ///
    /// Chunking is what lets the page show progress and stay responsive; the
    /// rows are the same ones a single long run would produce.
    #[wasm_bindgen]
    pub fn run_benchmark_chunk(
        &mut self,
        base_seed: u64,
        start: usize,
        count: usize,
        config: JsValue,
    ) -> Result<JsValue, JsError> {
        let cfg = read_config(&config)?;
        self.ensure_database(&cfg)?;
        let solver = Solver::new(&self.catalog, &self.db, &self.ids);
        let clock = PerformanceClock;

        let mut rows = Vec::with_capacity(count);
        for index in start..start + count {
            let row = bench::run_trial(
                bench::trial_seed(base_seed, index),
                index,
                &cfg,
                &solver,
                &clock,
                &mut self.trial_ws,
            );
            rows.push(TrialRow::from(&row));
        }
        serde_wasm_bindgen::to_value(&rows).map_err(JsError::from)
    }

    /// Runs a slew and reports tracking against the full search on every frame.
    ///
    /// `config` is a preset name or a `SimConfig`; `track` is a `TrackConfig`
    /// object, or null for the defaults.
    #[wasm_bindgen]
    pub fn run_tracking(
        &mut self,
        base_seed: u64,
        config: JsValue,
        track: JsValue,
    ) -> Result<JsValue, JsError> {
        let cfg = read_config(&config)?;
        self.ensure_database(&cfg)?;
        let track_cfg: TrackConfig = if track.is_undefined() || track.is_null() {
            TrackConfig::default()
        } else {
            serde_wasm_bindgen::from_value(track).map_err(JsError::from)?
        };
        let solver = Solver::new(&self.catalog, &self.db, &self.ids);
        let steps = track::run_sequence(
            base_seed,
            &cfg,
            &track_cfg,
            &solver,
            &PerformanceClock,
            &mut self.trial_ws,
        );
        let report = track::summarise(&steps);
        let run = TrackRun {
            steps: steps.iter().map(TrackStepDto::from).collect(),
            frames: report.frames,
            correct: report.correct,
            wrong_confident: report.wrong_confident,
            reacquisitions: report.reacquisitions,
            track_ms: report.track_ns as f64 / 1.0e6,
            lost_ms: report.lost_ns as f64 / 1.0e6,
            speedup: report.speedup(),
            median_arcsec: report.median_error / ARCSEC,
            step_deg: track_cfg.step_rad().to_degrees(),
            window_arcsec: track_cfg.aided_radius(&cfg) / ARCSEC,
        };
        serde_wasm_bindgen::to_value(&run).map_err(JsError::from)
    }

    /// Rebuilds the pair database if the configuration changed what goes in it.
    fn ensure_database(&mut self, cfg: &SimConfig) -> Result<(), JsError> {
        let key = db_key(cfg);
        if key != self.db_key {
            self.db = PairDb::build(&self.catalog, cfg).map_err(to_js)?;
            self.db_key = key;
        }
        Ok(())
    }
}

/// What the pair database depends on: the magnitude cut and the field.
fn db_key(cfg: &SimConfig) -> (f32, f64, usize) {
    (cfg.db_mag_cut, cfg.fov_deg, cfg.width)
}

/// A simulated frame, holding its image on the WASM side.
#[wasm_bindgen]
pub struct Frame {
    seed: u64,
    frame: SimFrame,
}

#[wasm_bindgen]
impl Frame {
    fn new(seed: u64, frame: SimFrame) -> Self {
        Self { seed, frame }
    }

    /// The seed this frame was rendered from.
    #[wasm_bindgen(getter)]
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// Sensor width in pixels.
    #[wasm_bindgen(getter)]
    pub fn width(&self) -> usize {
        self.frame.width
    }

    /// Sensor height in pixels.
    #[wasm_bindgen(getter)]
    pub fn height(&self) -> usize {
        self.frame.height
    }

    /// The image as a `Uint16Array`, copied out of WASM memory.
    #[wasm_bindgen(getter)]
    pub fn image(&self) -> Vec<u16> {
        self.frame.image.clone()
    }

    /// Address of the image inside WASM memory.
    ///
    /// With [`Frame::image_len`] this lets the caller build a `Uint16Array`
    /// view over the buffer rather than copying a megapixel. The view is only
    /// valid while this `Frame` is alive and only until WASM memory grows, so
    /// read it out before allocating anything substantial.
    #[wasm_bindgen(getter)]
    pub fn image_ptr(&self) -> *const u16 {
        self.frame.image.as_ptr()
    }

    /// Length of the image, in samples.
    #[wasm_bindgen(getter)]
    pub fn image_len(&self) -> usize {
        self.frame.image.len()
    }

    /// The true attitude and the frame's truth stars.
    #[wasm_bindgen(getter)]
    pub fn truth(&self) -> Result<JsValue, JsError> {
        let truth = FrameTruth {
            attitude: AttitudeDto::new(&self.frame.a_true, 1.0),
            stars: self.frame.truth.iter().map(TruthDto::from).collect(),
            false_stars: self.frame.false_stars.clone(),
        };
        serde_wasm_bindgen::to_value(&truth).map_err(JsError::from)
    }
}

/// Reads a configuration from a preset name or a whole `SimConfig` object.
fn read_config(value: &JsValue) -> Result<SimConfig, JsError> {
    if let Some(name) = value.as_string() {
        let preset: Preset = name.parse().map_err(to_js)?;
        return Ok(SimConfig::preset(preset));
    }
    if value.is_undefined() || value.is_null() {
        return Ok(SimConfig::default());
    }
    serde_wasm_bindgen::from_value(value.clone()).map_err(JsError::from)
}

/// Turns a core error into something JS can read.
fn to_js(error: tracker_core::Error) -> JsError {
    JsError::new(&error.to_string())
}

// --- data handed to JavaScript -------------------------------------------

/// Radians per arcsecond.
const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// An attitude in the forms the UI shows.
#[derive(Serialize)]
struct AttitudeDto {
    /// Hamilton quaternion, scalar-first, `w >= 0`.
    quaternion: [f64; 4],
    /// Boresight right ascension, degrees.
    ra_deg: f64,
    /// Boresight declination, degrees.
    dec_deg: f64,
    /// Roll about the boresight, degrees.
    roll_deg: f64,
    /// Fitted focal length over the nominal one.
    focal_scale: f64,
}

impl AttitudeDto {
    fn new(attitude: &math::Mat3, focal_scale: f64) -> Self {
        let (ra, dec, roll) = bench::boresight_radec_roll(attitude);
        Self {
            quaternion: math::matrix_to_quat(attitude),
            ra_deg: ra.to_degrees(),
            dec_deg: dec.to_degrees(),
            roll_deg: roll.to_degrees(),
            focal_scale,
        }
    }
}

/// One detected blob.
#[derive(Serialize)]
struct CentroidDto {
    x: f64,
    y: f64,
    flux: f64,
    pixels: u32,
}

impl From<&Centroid> for CentroidDto {
    fn from(centroid: &Centroid) -> Self {
        Self {
            x: centroid.x,
            y: centroid.y,
            flux: centroid.flux,
            pixels: centroid.pixels,
        }
    }
}

/// One identification.
#[derive(Serialize)]
struct MatchDto {
    /// Index into the centroid list.
    observed: u16,
    /// Catalogue index.
    catalog_index: u16,
    /// Catalogue (HIP) number.
    id: u32,
    /// IAU proper name, for the 335 catalogue entries that have one.
    name: Option<String>,
    /// Whether the pyramid confirmed this star.
    from_pyramid: bool,
}

impl MatchDto {
    /// Builds the DTO, looking the proper name up in `catalog`.
    ///
    /// Not a `From` impl because the name is not on the match: it lives in the
    /// catalogue entry that `catalog_index` points at.
    fn new(matched: &StarMatch, catalog: &Catalog) -> Self {
        Self {
            observed: matched.observed,
            catalog_index: matched.catalog_index,
            id: matched.id,
            name: catalog
                .stars
                .get(usize::from(matched.catalog_index))
                .and_then(|star| star.name.clone()),
            from_pyramid: matched.from_pyramid,
        }
    }
}

/// The truth-blind self-check.
#[derive(Serialize)]
struct SelfCheckDto {
    matched: usize,
    residual_rms_px: f64,
    pyramid_stars: usize,
    claims_explained: usize,
    claims: usize,
    confident: bool,
}

impl From<&SelfCheck> for SelfCheckDto {
    fn from(check: &SelfCheck) -> Self {
        Self {
            matched: check.matched,
            residual_rms_px: check.residual_rms_px,
            pyramid_stars: check.pyramid_stars,
            claims_explained: check.claims_explained,
            claims: check.claims,
            confident: check.confident,
        }
    }
}

/// What [`Tracker::solve_frame`] returns.
#[derive(Serialize)]
struct SolveReport {
    detected: usize,
    centroids: Vec<CentroidDto>,
    matches: Vec<MatchDto>,
    attitude: Option<AttitudeDto>,
    self_check: Option<SelfCheckDto>,
    identify_tries: u32,
    identify_ms: f64,
}

/// One truth star.
#[derive(Serialize)]
struct TruthDto {
    id: u32,
    x: f64,
    y: f64,
    mag: f32,
    dropped: bool,
}

impl From<&simulate::TruthStar> for TruthDto {
    fn from(star: &simulate::TruthStar) -> Self {
        Self {
            id: star.id,
            x: star.pixel[0],
            y: star.pixel[1],
            mag: star.mag,
            dropped: star.dropped,
        }
    }
}

/// Ground truth for a frame.
#[derive(Serialize)]
struct FrameTruth {
    attitude: AttitudeDto,
    stars: Vec<TruthDto>,
    false_stars: Vec<[f64; 2]>,
}

/// One scored trial, in the units the UI reports.
#[derive(Serialize)]
struct TrialRow {
    index: usize,
    /// Derived seed, as a string: it does not fit a JS number.
    seed: String,
    outcome: String,
    total_arcsec: f64,
    cross_arcsec: f64,
    roll_arcsec: f64,
    claimed: usize,
    correct: usize,
    available: usize,
    detected: usize,
    used: usize,
    matched: usize,
    residual_rms_px: f64,
    focal_scale: f64,
    quaternion: [f64; 4],
    tries: u32,
    simulate_ms: f64,
    centroid_ms: f64,
    identify_ms: f64,
    attitude_ms: f64,
    verify_ms: f64,
}

impl From<&TrialResult> for TrialRow {
    fn from(row: &TrialResult) -> Self {
        Self {
            index: row.index,
            seed: row.seed.to_string(),
            outcome: row.outcome.as_str().to_string(),
            total_arcsec: row.total_error / ARCSEC,
            cross_arcsec: row.cross_error / ARCSEC,
            roll_arcsec: row.roll_error / ARCSEC,
            claimed: row.claimed,
            correct: row.correct,
            available: row.available,
            detected: row.detected,
            used: row.used,
            matched: row.matched,
            residual_rms_px: row.residual_rms_px,
            focal_scale: row.focal_scale,
            quaternion: row.quaternion,
            tries: row.diagnostics.tries,
            simulate_ms: row.timings.simulate_ns as f64 / 1.0e6,
            centroid_ms: row.timings.centroid_ns as f64 / 1.0e6,
            identify_ms: row.timings.identify_ns as f64 / 1.0e6,
            attitude_ms: row.timings.attitude_ns as f64 / 1.0e6,
            verify_ms: row.timings.verify_ns as f64 / 1.0e6,
        }
    }
}

/// A trial plus everything a drawing of it needs.
#[derive(Serialize)]
struct TrialDetail {
    trial: TrialRow,
    truth: Vec<TruthDto>,
    false_stars: Vec<[f64; 2]>,
    true_attitude: AttitudeDto,
}

/// One frame of a tracking sequence, as the UI reads it.
#[derive(Serialize)]
struct TrackStepDto {
    index: usize,
    outcome: String,
    aided: bool,
    reacquired: bool,
    detected: usize,
    claimed: usize,
    correct: usize,
    matched: usize,
    total_arcsec: f64,
    track_ms: f64,
    lost_ms: f64,
    /// Estimated boresight, degrees.
    ra_deg: f64,
    dec_deg: f64,
    /// True boresight, degrees.
    true_ra_deg: f64,
    true_dec_deg: f64,
}

impl From<&track::TrackStep> for TrackStepDto {
    fn from(step: &track::TrackStep) -> Self {
        let estimate = quat_to_matrix(step.quaternion);
        let truth = quat_to_matrix(step.truth_quaternion);
        let (ra, dec, _) = bench::boresight_radec_roll(&estimate);
        let (true_ra, true_dec, _) = bench::boresight_radec_roll(&truth);
        Self {
            index: step.index,
            outcome: step.outcome.as_str().to_string(),
            aided: step.aided,
            reacquired: step.reacquired,
            detected: step.detected,
            claimed: step.claimed,
            correct: step.correct,
            matched: step.matched,
            total_arcsec: step.total_error / ARCSEC,
            track_ms: step.track_ns as f64 / 1.0e6,
            lost_ms: step.lost_ns as f64 / 1.0e6,
            ra_deg: ra.to_degrees(),
            dec_deg: dec.to_degrees(),
            true_ra_deg: true_ra.to_degrees(),
            true_dec_deg: true_dec.to_degrees(),
        }
    }
}

/// A whole tracking run: the frames and the figures that summarise them.
#[derive(Serialize)]
struct TrackRun {
    steps: Vec<TrackStepDto>,
    frames: usize,
    correct: usize,
    wrong_confident: usize,
    reacquisitions: usize,
    track_ms: f64,
    lost_ms: f64,
    speedup: f64,
    median_arcsec: f64,
    step_deg: f64,
    window_arcsec: f64,
}
