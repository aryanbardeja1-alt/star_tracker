//! Image simulator: a seed and a [`SimConfig`] in, a rendered 12-bit frame
//! plus its ground truth out.
//!
//! The simulator is deliberately *not* a mirror of the solver's assumptions,
//! which is what keeps the benchmark from being an inverse crime:
//!
//! - stars are rendered down to [`SimConfig::sim_mag_limit`] (6.5), fainter
//!   than the database cut (6.0), so the frame contains stars the solver
//!   cannot know about;
//! - rendered magnitudes are jittered;
//! - the true focal length and radial distortion may differ from the nominal
//!   camera the solver uses, by amounts the preset controls;
//! - false stars, dropouts and sensor defects are injected.
//!
//! All randomness comes from a seeded `ChaCha8Rng`, consumed in a fixed order,
//! so the same seed gives a bit-identical frame on any platform.

use crate::camera::CameraModel;
use crate::catalog::Catalog;
use crate::math::{Mat3, quat_to_matrix};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rand_distr::StandardNormal;
use serde::{Deserialize, Serialize};
use std::f64::consts::SQRT_2;

/// Gaussian flux is integrated out to this many sigma.
///
/// Six rather than four: a finite window truncates the profile asymmetrically
/// once the star sits off a pixel centre, and at four sigma that alone biases
/// the recovered centroid by 1.4e-4 px, against the 0.02 px bias Phase 3 has
/// to measure. Six sigma drops it to 6e-9 px and loses 1e-9 of the flux, for
/// 169 pixels of work per star instead of 81 -- nothing beside the per-pixel
/// noise pass over the whole frame.
const PSF_WINDOW_SIGMAS: f64 = 6.0;

/// Below this mean, shot noise is drawn as an exact Poisson variate; above it
/// the Gaussian form is used. At a mean of 30 the Poisson skewness is 0.18, so
/// the approximation is good, and with the default 100 e- background every
/// pixel takes the cheap branch anyway.
const POISSON_EXACT_BELOW: f64 = 30.0;

/// Standard deviation of the wide blob the `brutal` preset adds, in pixels.
const BLOB_SIGMA_PX: f64 = 40.0;

/// Effective magnitude of that blob. With the default photometry this drives
/// the central pixels far past the full well and leaves a saturated disc
/// roughly 65 px in radius, which is what a sunlit body in the field looks
/// like; a merely bright star would not saturate at all.
const BLOB_MAGNITUDE: f64 = -8.0;

/// Difficulty presets, as tabulated in docs/SPEC.md.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Preset {
    /// No defects at all: the solver's assumptions hold exactly.
    Easy,
    /// The default operating point.
    Nominal,
    /// Noticeable noise, dropouts, false stars and a miscalibrated camera.
    Hard,
    /// Everything at once, plus a saturated bright object in the field.
    Brutal,
}

impl Preset {
    /// Every preset, in increasing difficulty.
    pub const ALL: [Preset; 4] = [Preset::Easy, Preset::Nominal, Preset::Hard, Preset::Brutal];

    /// The lowercase name used on the command line and in JSON.
    pub fn as_str(self) -> &'static str {
        match self {
            Preset::Easy => "easy",
            Preset::Nominal => "nominal",
            Preset::Hard => "hard",
            Preset::Brutal => "brutal",
        }
    }
}

impl std::fmt::Display for Preset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for Preset {
    type Err = crate::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Preset::ALL
            .into_iter()
            .find(|p| p.as_str() == s)
            .ok_or_else(|| crate::Error::UnknownPreset(s.to_string()))
    }
}

/// Every parameter of a simulated trial, in one serialisable struct shared by
/// the CLI and the web front end.
///
/// Photometric quantities are in electrons; angles in this struct are in
/// degrees because it is a user-facing boundary, and are converted on use.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SimConfig {
    // --- Sensor and optics ---
    /// Sensor width in pixels.
    pub width: usize,
    /// Sensor height in pixels.
    pub height: usize,
    /// Nominal full-width field of view, degrees. The solver assumes this.
    pub fov_deg: f64,
    /// Gaussian point-spread standard deviation, pixels.
    pub psf_sigma_px: f64,
    /// Output bit depth; 12 gives a 0..4095 range.
    pub bit_depth: u32,

    // --- Photometry, in electrons ---
    /// Electrons per second from a magnitude-0 star.
    pub f0_electrons_per_s: f64,
    /// Exposure time, seconds.
    pub exposure_s: f64,
    /// Sky and dark background per pixel, electrons.
    pub background_e: f64,
    /// Read noise, electrons RMS.
    pub read_noise_e: f64,
    /// Electrons per ADU. With `bit_depth` this fixes the full well, and
    /// therefore where saturation begins.
    pub gain_e_per_adu: f64,

    // --- Catalogue limits ---
    /// Faintest magnitude rendered into the image.
    pub sim_mag_limit: f32,
    /// Faintest magnitude the solver's database may contain. Kept here so one
    /// struct configures a whole trial; used from Phase 4 on.
    pub db_mag_cut: f32,
    /// Standard deviation of the per-star magnitude jitter, magnitudes.
    pub mag_jitter: f64,

    // --- Errors the solver is not told about ---
    /// Fractional error on the true focal length, e.g. 0.002 for 0.2%.
    pub focal_error_frac: f64,
    /// True radial distortion coefficient. The solver assumes zero.
    pub k1: f64,

    // --- Defects ---
    /// Inclusive range for the number of false stars injected.
    pub false_star_range: [u32; 2],
    /// Probability that a visible star is dropped from the image.
    pub dropout_prob: f64,
    /// Number of stuck-high hot pixels.
    pub hot_pixels: u32,
    /// Whether to add a wide saturated blob, standing in for the Moon or an
    /// Earth limb.
    pub bright_blob: bool,

    // --- Detection ---
    /// How many of the brightest centroids the solver keeps.
    pub max_centroids: usize,
    /// Stride of the sub-sampled grid the background and noise level are
    /// estimated from.
    pub background_stride: usize,
    /// Detection threshold, in noise sigmas above the background median.
    pub centroid_k_sigma: f64,
    /// Smallest blob accepted, in pixels. Two rejects single hot pixels.
    pub centroid_min_pixels: u32,
    /// Largest blob accepted, in pixels. This is what discards a sunlit body:
    /// the brightest real star spans about 60 px, the wide blob of the
    /// `brutal` preset nearly 20,000.
    pub centroid_max_pixels: u32,

    // --- Identification ---
    /// Centroid error the solver *assumes*, in pixels per axis, which sets the
    /// search tolerance. Phase 3 measured 0.011 px at magnitude 4 and about
    /// 0.04 px at the faint end, so the 0.1 default is deliberately generous.
    pub centroid_sigma_px: f64,
    /// Search tolerance in sigmas: the `k_sigma` of
    /// `eps = k_sigma * sigma_angle + calib_margin`.
    pub id_k_sigma: f64,
    /// Fixed addition to the search tolerance, arcseconds, standing in for
    /// calibration error the solver does not model.
    pub calib_margin_arcsec: f64,
    /// Fractional scale error the solver allows for when matching *pairs*.
    ///
    /// A focal length off by `d` stretches every observed angle by about `d`,
    /// so the error in a pair's angle grows with the angle itself and no
    /// constant tolerance can cover it: measured on the `hard` preset, a 0.2%
    /// focal error puts a pair out by 8 arcseconds per degree of separation,
    /// which at 20 degrees is 160 arcseconds against a 50 arcsecond constant.
    /// This is the solver's own assumption about its calibration, not knowledge
    /// of any preset, and it is deliberately the same for all of them.
    pub id_focal_tolerance_frac: f64,
    /// Cap on how many star triples the identifier will try before giving up.
    pub id_max_tries: u32,

    // --- Verification ---
    /// How close a reprojected catalogue star must fall to a centroid to count
    /// as matched by the self-check, in pixels.
    pub verify_match_px: f64,
    /// Fewest self-check matches a solution needs to be called CONFIDENT.
    pub verify_min_matches: usize,
    /// Largest reprojection residual RMS a CONFIDENT solution may have, in
    /// pixels.
    pub verify_max_rms_px: f64,
    /// Attitude error above which a solution counts as wrong rather than
    /// correct, in arcseconds.
    pub err_threshold_arcsec: f64,
    /// How close a centroid must be to the true position of the star whose id
    /// was claimed, for the ground-truth check to call that claim right.
    pub truth_match_px: f64,
}

impl Default for SimConfig {
    /// The `nominal` preset.
    fn default() -> Self {
        Self::preset(Preset::Nominal)
    }
}

impl SimConfig {
    /// The configuration for one of the tabulated presets.
    pub fn preset(preset: Preset) -> Self {
        // Shared baseline: the default camera from docs/SPEC.md. The photometry
        // is scaled so a magnitude 6.5 star -- the faintest rendered -- still
        // reaches a total signal-to-noise of about 36 on the nominal preset.
        let base = Self {
            width: 1024,
            height: 1024,
            fov_deg: 20.0,
            psf_sigma_px: 1.0,
            bit_depth: 12,
            f0_electrons_per_s: 1.0e7,
            exposure_s: 0.1,
            background_e: 100.0,
            read_noise_e: 10.0,
            gain_e_per_adu: 5.0,
            sim_mag_limit: 6.5,
            db_mag_cut: 6.0,
            mag_jitter: 0.1,
            focal_error_frac: 0.0,
            k1: 0.0,
            false_star_range: [0, 0],
            dropout_prob: 0.0,
            hot_pixels: 5,
            bright_blob: false,
            max_centroids: 15,
            background_stride: 8,
            centroid_k_sigma: 5.0,
            centroid_min_pixels: 2,
            centroid_max_pixels: 200,
            centroid_sigma_px: 0.1,
            id_k_sigma: 4.0,
            calib_margin_arcsec: 10.0,
            id_focal_tolerance_frac: 0.003,
            id_max_tries: 1000,
            verify_match_px: 1.0,
            verify_min_matches: 6,
            verify_max_rms_px: 0.5,
            err_threshold_arcsec: 60.0,
            truth_match_px: 1.5,
        };

        match preset {
            Preset::Easy => Self {
                read_noise_e: 5.0,
                ..base
            },
            Preset::Nominal => Self {
                false_star_range: [0, 1],
                dropout_prob: 0.05,
                ..base
            },
            Preset::Hard => Self {
                read_noise_e: 25.0,
                false_star_range: [0, 3],
                dropout_prob: 0.10,
                focal_error_frac: 0.002,
                k1: 0.01,
                ..base
            },
            Preset::Brutal => Self {
                read_noise_e: 50.0,
                false_star_range: [2, 5],
                dropout_prob: 0.20,
                focal_error_frac: 0.005,
                k1: 0.03,
                bright_blob: true,
                ..base
            },
        }
    }

    /// The camera the *solver* is given: nominal focal length, no distortion.
    pub fn nominal_camera(&self) -> CameraModel {
        CameraModel::from_fov(self.width, self.height, self.fov_deg.to_radians())
    }

    /// The camera that actually took the picture: focal length and distortion
    /// perturbed by the preset.
    pub fn true_camera(&self) -> CameraModel {
        let mut camera = self.nominal_camera();
        camera.f *= 1.0 + self.focal_error_frac;
        camera.k1 = self.k1;
        camera
    }

    /// Largest ADU value the sensor can output.
    pub fn saturation_adu(&self) -> u16 {
        // 12 bits -> 4095. Saturating arithmetic keeps absurd bit depths sane.
        ((1u32 << self.bit_depth.min(16)) - 1).min(u32::from(u16::MAX)) as u16
    }

    /// Total electrons collected from a star of the given magnitude:
    /// `F0 * 10^(-0.4 m) * t_exp`.
    pub fn flux_electrons(&self, mag: f64) -> f64 {
        self.f0_electrons_per_s * libm::pow(10.0, -0.4 * mag) * self.exposure_s
    }
}

/// One rendered star, as the simulator knows it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TruthStar {
    /// Catalogue (HIP) number.
    pub id: u32,
    /// True position on the sensor, in pixels, from the *true* camera.
    pub pixel: [f64; 2],
    /// Magnitude actually rendered, after jitter.
    pub mag: f32,
    /// Whether this star was dropped from the image rather than rendered.
    pub dropped: bool,
}

/// A simulated frame and everything true about it.
#[derive(Clone, Debug)]
pub struct SimFrame {
    /// The image, row-major, `width * height` entries, in ADU.
    pub image: Vec<u16>,
    /// Sensor width in pixels.
    pub width: usize,
    /// Sensor height in pixels.
    pub height: usize,
    /// True attitude, mapping inertial to camera: `b = A r`.
    pub a_true: Mat3,
    /// Every catalogue star that landed on the sensor, in catalogue order.
    pub truth: Vec<TruthStar>,
    /// Pixel positions of the injected false stars.
    pub false_stars: Vec<[f64; 2]>,
}

impl SimFrame {
    /// Truth stars that were actually rendered, i.e. not dropped.
    pub fn rendered(&self) -> impl Iterator<Item = &TruthStar> {
        self.truth.iter().filter(|s| !s.dropped)
    }
}

/// Scratch buffers reused across frames, so rendering allocates nothing per
/// star.
#[derive(Clone, Debug, Default)]
pub struct Workspace {
    /// Accumulated signal in electrons, one entry per pixel.
    electrons: Vec<f64>,
    /// Per-column flux fractions for the star currently being rendered. The
    /// pixel-integrated Gaussian is separable, so each column's fraction is
    /// computed once and reused down the rows.
    column_fractions: Vec<f64>,
}

impl Workspace {
    /// An empty workspace; buffers size themselves on first use.
    pub fn new() -> Self {
        Self::default()
    }

    fn reset(&mut self, pixels: usize) {
        self.electrons.clear();
        self.electrons.resize(pixels, 0.0);
    }
}

/// Simulates one frame.
///
/// `seed` is the fully derived per-trial seed; the benchmark computes it as
/// `splitmix64(base_seed ^ trial_index)`. The catalogue supplies the stars and
/// `ws` the scratch buffers.
pub fn simulate(seed: u64, cfg: &SimConfig, catalog: &Catalog, ws: &mut Workspace) -> SimFrame {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let a_true = random_attitude(&mut rng);
    render(a_true, &mut rng, cfg, catalog, ws)
}

/// Renders a frame at an attitude the caller supplies, rather than a drawn one.
///
/// Tracking needs this: along a slew each attitude follows from the last, not
/// from a seed. `a_true` maps inertial to camera (`b = A r`). The seed still
/// drives magnitude jitter, dropouts, false stars and all the noise, so two
/// frames at one attitude with different seeds differ exactly as two trials do.
pub fn simulate_at(
    a_true: Mat3,
    seed: u64,
    cfg: &SimConfig,
    catalog: &Catalog,
    ws: &mut Workspace,
) -> SimFrame {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    render(a_true, &mut rng, cfg, catalog, ws)
}

/// The rendering both entry points share, with the attitude already decided.
///
/// The draw order is what makes [`simulate`] reproducible, so the attitude is
/// taken from `rng` before this is called and everything else after.
fn render(
    a_true: Mat3,
    rng: &mut ChaCha8Rng,
    cfg: &SimConfig,
    catalog: &Catalog,
    ws: &mut Workspace,
) -> SimFrame {
    let pixels = cfg.width * cfg.height;
    ws.reset(pixels);

    let camera = cfg.true_camera();

    // Stars, in catalogue order so the RNG is consumed deterministically.
    let mut truth = Vec::new();
    for star in &catalog.stars {
        if star.mag > cfg.sim_mag_limit {
            continue;
        }
        let b = a_true * star.direction();
        let Some(pixel) = camera.project(&b) else {
            continue;
        };
        if !camera.contains(pixel[0], pixel[1]) {
            continue;
        }

        let jitter: f64 = rng.sample::<f64, _>(StandardNormal) * cfg.mag_jitter;
        let mag = f64::from(star.mag) + jitter;
        let dropped = rng.random::<f64>() < cfg.dropout_prob;

        if !dropped {
            render_gaussian(ws, cfg, pixel, cfg.flux_electrons(mag), cfg.psf_sigma_px);
        }
        truth.push(TruthStar {
            id: star.id,
            pixel,
            mag: mag as f32,
            dropped,
        });
    }

    // False stars: plausible but unrelated blobs, over the brighter half of
    // the rendered magnitude range so they genuinely compete for detection.
    let [low, high] = cfg.false_star_range;
    let count = if high > low {
        rng.random_range(low..=high)
    } else {
        low
    };
    let mut false_stars = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let x = rng.random::<f64>() * (cfg.width as f64 - 1.0);
        let y = rng.random::<f64>() * (cfg.height as f64 - 1.0);
        let mag = 3.0 + rng.random::<f64>() * (f64::from(cfg.sim_mag_limit) - 3.0);
        render_gaussian(ws, cfg, [x, y], cfg.flux_electrons(mag), cfg.psf_sigma_px);
        false_stars.push([x, y]);
    }

    if cfg.bright_blob {
        let x = rng.random::<f64>() * (cfg.width as f64 - 1.0);
        let y = rng.random::<f64>() * (cfg.height as f64 - 1.0);
        let flux = cfg.flux_electrons(BLOB_MAGNITUDE);
        render_gaussian(ws, cfg, [x, y], flux, BLOB_SIGMA_PX);
    }

    SimFrame {
        image: read_out(ws, cfg, rng),
        width: cfg.width,
        height: cfg.height,
        a_true,
        truth,
        false_stars,
    }
}

/// Renders a frame containing exactly the stars given as `(pixel, magnitude)`,
/// with no attitude or catalogue involved.
///
/// Same optics, noise, background, saturation and hot pixels as [`simulate`] --
/// it shares the same read-out path -- but the caller places the stars. This is
/// what the centroider is calibrated against, where a controlled sub-pixel
/// offset is the whole point of the measurement.
pub fn simulate_stars(
    seed: u64,
    cfg: &SimConfig,
    stars: &[([f64; 2], f64)],
    ws: &mut Workspace,
) -> Vec<u16> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    ws.reset(cfg.width * cfg.height);
    for &(pixel, mag) in stars {
        render_gaussian(ws, cfg, pixel, cfg.flux_electrons(mag), cfg.psf_sigma_px);
    }
    read_out(ws, cfg, &mut rng)
}

/// Turns accumulated electrons into a quantised frame: shot and read noise,
/// then the sensor's bit depth, then stuck-high hot pixels.
///
/// Shared by [`simulate`] and [`simulate_stars`] so there is exactly one noise
/// model in the crate.
fn read_out(ws: &Workspace, cfg: &SimConfig, rng: &mut ChaCha8Rng) -> Vec<u16> {
    let saturation = cfg.saturation_adu();
    let mut image = vec![0u16; ws.electrons.len()];

    for (out, signal) in image.iter_mut().zip(&ws.electrons) {
        let electrons = sample_pixel(*signal + cfg.background_e, cfg.read_noise_e, rng);
        let adu = (electrons / cfg.gain_e_per_adu).round();
        *out = if adu <= 0.0 {
            0
        } else if adu >= f64::from(saturation) {
            saturation
        } else {
            adu as u16
        };
    }

    let pixels = image.len();
    for _ in 0..cfg.hot_pixels {
        image[rng.random_range(0..pixels)] = saturation;
    }
    image
}

/// A uniformly distributed attitude.
///
/// Four independent N(0,1) samples normalised into a quaternion give a point
/// uniform on the unit 3-sphere, which is Haar-uniform over rotations.
fn random_attitude(rng: &mut ChaCha8Rng) -> Mat3 {
    let q = [
        rng.sample::<f64, _>(StandardNormal),
        rng.sample::<f64, _>(StandardNormal),
        rng.sample::<f64, _>(StandardNormal),
        rng.sample::<f64, _>(StandardNormal),
    ];
    // `quat_to_matrix` normalises, so the raw Gaussian quadruple is fine.
    quat_to_matrix(q)
}

/// Fraction of a Gaussian's flux landing in the pixel centred on `centre`.
///
/// Integrates across the pixel's unit extent with erf rather than sampling the
/// Gaussian at the pixel centre; sampling biases the recovered centroid, which
/// is exactly what Phase 3 measures.
fn pixel_fraction(centre: f64, mean: f64, sigma: f64) -> f64 {
    let scale = sigma * SQRT_2;
    let lo = (centre - 0.5 - mean) / scale;
    let hi = (centre + 0.5 - mean) / scale;
    0.5 * (libm::erf(hi) - libm::erf(lo))
}

/// Adds a pixel-integrated Gaussian of total flux `flux` electrons.
fn render_gaussian(ws: &mut Workspace, cfg: &SimConfig, centre: [f64; 2], flux: f64, sigma: f64) {
    let radius = (PSF_WINDOW_SIGMAS * sigma).ceil() as i64;
    let (centre_x, centre_y) = (centre[0], centre[1]);
    let col_lo = (centre_x.round() as i64 - radius).max(0);
    let col_hi = (centre_x.round() as i64 + radius).min(cfg.width as i64 - 1);
    let row_lo = (centre_y.round() as i64 - radius).max(0);
    let row_hi = (centre_y.round() as i64 + radius).min(cfg.height as i64 - 1);
    if col_lo > col_hi || row_lo > row_hi {
        return;
    }

    // The 2-D integral is separable, so each column's fraction is computed
    // once here instead of once per (row, column) pair.
    ws.column_fractions.clear();
    for col in col_lo..=col_hi {
        ws.column_fractions
            .push(pixel_fraction(col as f64, centre_x, sigma));
    }

    for row in row_lo..=row_hi {
        let row_fraction = pixel_fraction(row as f64, centre_y, sigma);
        if row_fraction == 0.0 {
            continue;
        }
        let offset = row as usize * cfg.width;
        for (col, column_fraction) in (col_lo..=col_hi).zip(&ws.column_fractions) {
            ws.electrons[offset + col as usize] += flux * row_fraction * column_fraction;
        }
    }
}

/// Draws one pixel's read-out in electrons: shot noise on `mean`, plus read
/// noise.
fn sample_pixel(mean: f64, read_noise: f64, rng: &mut ChaCha8Rng) -> f64 {
    if mean < POISSON_EXACT_BELOW {
        let counts = poisson_small(mean, rng);
        counts + read_noise * rng.sample::<f64, _>(StandardNormal)
    } else {
        // Shot and read noise are independent Gaussians here, so they combine
        // into a single draw.
        let sigma = (mean + read_noise * read_noise).sqrt();
        mean + sigma * rng.sample::<f64, _>(StandardNormal)
    }
}

/// Exact Poisson variate by Knuth's product method, for small means only.
///
/// Uses `libm::exp` so the threshold is identical on every target.
fn poisson_small(mean: f64, rng: &mut ChaCha8Rng) -> f64 {
    if mean <= 0.0 {
        return 0.0;
    }
    let limit = libm::exp(-mean);
    let mut product = 1.0;
    let mut count = 0u32;
    loop {
        product *= rng.random::<f64>();
        if product <= limit {
            return f64::from(count);
        }
        count += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::angle_between;
    use std::str::FromStr;

    const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

    fn catalog() -> Catalog {
        Catalog::from_bytes(CATALOG_BIN).expect("committed catalog.bin must decode")
    }

    fn frame(seed: u64, preset: Preset, catalog: &Catalog) -> SimFrame {
        simulate(
            seed,
            &SimConfig::preset(preset),
            catalog,
            &mut Workspace::new(),
        )
    }

    // --- presets and configuration ---

    #[test]
    fn presets_match_the_documented_table() {
        let easy = SimConfig::preset(Preset::Easy);
        let nominal = SimConfig::preset(Preset::Nominal);
        let hard = SimConfig::preset(Preset::Hard);
        let brutal = SimConfig::preset(Preset::Brutal);

        // Read noise rises monotonically: low, normal, high, very high.
        assert!(
            easy.read_noise_e < nominal.read_noise_e
                && nominal.read_noise_e < hard.read_noise_e
                && hard.read_noise_e < brutal.read_noise_e
        );
        // False stars and dropout rise too.
        assert_eq!(easy.false_star_range, [0, 0]);
        assert_eq!(nominal.false_star_range, [0, 1]);
        assert_eq!(hard.false_star_range, [0, 3]);
        assert_eq!(brutal.false_star_range, [2, 5]);
        assert_eq!(
            [
                easy.dropout_prob,
                nominal.dropout_prob,
                hard.dropout_prob,
                brutal.dropout_prob
            ],
            [0.0, 0.05, 0.10, 0.20]
        );
        // Only hard and brutal miscalibrate the camera.
        assert_eq!((easy.focal_error_frac, easy.k1), (0.0, 0.0));
        assert_eq!((nominal.focal_error_frac, nominal.k1), (0.0, 0.0));
        assert_eq!(hard.focal_error_frac, 0.002);
        assert_eq!(brutal.focal_error_frac, 0.005);
        assert!(hard.k1 > 0.0 && brutal.k1 > hard.k1);
        // Only brutal adds the bright object.
        assert!(!hard.bright_blob && brutal.bright_blob);
        // The simulator must see fainter than the database, or the benchmark
        // would be an inverse crime.
        assert!(nominal.sim_mag_limit > nominal.db_mag_cut);
    }

    #[test]
    fn default_is_the_nominal_preset() {
        assert_eq!(SimConfig::default(), SimConfig::preset(Preset::Nominal));
    }

    #[test]
    fn presets_round_trip_through_their_names() {
        for preset in Preset::ALL {
            assert_eq!(Preset::from_str(preset.as_str()).expect("parse"), preset);
            assert_eq!(preset.to_string(), preset.as_str());
        }
        assert!(matches!(
            Preset::from_str("medium"),
            Err(crate::Error::UnknownPreset(_))
        ));
    }

    #[test]
    fn config_serialises_round_trip() {
        // bincode rather than JSON: serde_json is not on the approved
        // dependency list, and this only needs to prove the derives work.
        let cfg = SimConfig::preset(Preset::Brutal);
        let bytes =
            bincode::serde::encode_to_vec(&cfg, bincode::config::standard()).expect("serialise");
        let (back, _): (SimConfig, usize) =
            bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
                .expect("deserialise");
        assert_eq!(back, cfg);
    }

    #[test]
    fn true_camera_differs_from_nominal_only_where_the_preset_says() {
        let nominal = SimConfig::preset(Preset::Nominal);
        assert_eq!(nominal.true_camera(), nominal.nominal_camera());

        let hard = SimConfig::preset(Preset::Hard);
        let (truth, assumed) = (hard.true_camera(), hard.nominal_camera());
        assert!((truth.f / assumed.f - 1.002).abs() < 1e-12);
        assert_eq!(truth.k1, 0.01);
        assert_eq!(assumed.k1, 0.0);
    }

    #[test]
    fn saturation_follows_the_bit_depth() {
        let mut cfg = SimConfig::default();
        assert_eq!(cfg.saturation_adu(), 4095);
        cfg.bit_depth = 8;
        assert_eq!(cfg.saturation_adu(), 255);
        cfg.bit_depth = 16;
        assert_eq!(cfg.saturation_adu(), 65535);
    }

    #[test]
    fn flux_follows_pogsons_law() {
        let cfg = SimConfig::default();
        // Five magnitudes is a factor of exactly 100 in flux.
        let bright = cfg.flux_electrons(1.0);
        let faint = cfg.flux_electrons(6.0);
        assert!((bright / faint - 100.0).abs() < 1e-9);
        // Magnitude 0 collects F0 * t_exp.
        assert!((cfg.flux_electrons(0.0) - cfg.f0_electrons_per_s * cfg.exposure_s).abs() < 1e-6);
    }

    // --- the pixel-integrated PSF ---

    #[test]
    fn psf_integrates_to_the_whole_flux() {
        // Summing the separable fractions over a wide window must recover
        // essentially all of the flux.
        for sigma in [0.5, 1.0, 2.5] {
            for offset in [0.0, 0.25, 0.5] {
                let mut total = 0.0;
                for pixel in -30..=30 {
                    total += pixel_fraction(f64::from(pixel), offset, sigma);
                }
                assert!(
                    (total - 1.0).abs() < 1e-9,
                    "sigma {sigma} offset {offset} integrated to {total}"
                );
            }
        }
    }

    /// The renderer's own intrinsic centroid bias, measured by driving
    /// `render_gaussian` and taking the first moment of what it wrote, rather
    /// than by reimplementing the integral in the test.
    ///
    /// Two effects set the floor here. The window is finite, so once the star
    /// sits off a pixel centre the profile is truncated asymmetrically; and
    /// binning a Gaussian into pixels displaces the discrete first moment by
    /// the usual aliasing term, of order `exp(-2*pi^2*sigma^2)`. At the
    /// default 1 px sigma with a 6 sigma window both are under 1e-8 px, five
    /// orders inside the 0.02 px bias Phase 3 must measure. The aliasing term
    /// is what punishes an undersampled PSF: at sigma = 0.5 px it alone is
    /// 2e-3 px.
    #[test]
    fn renderer_centroid_bias_is_far_inside_the_phase_3_budget() {
        let cfg = SimConfig {
            width: 64,
            height: 64,
            ..SimConfig::default()
        };
        let pixels = cfg.width * cfg.height;
        let mut ws = Workspace::new();

        let mut worst = 0.0f64;
        for step in 0..37 {
            // Sweep a full pixel diagonally, so both axes see every offset.
            let offset = f64::from(step) / 37.0;
            let centre = [32.0 + offset, 32.0 + offset * 0.5];

            ws.reset(pixels);
            render_gaussian(&mut ws, &cfg, centre, 1.0, cfg.psf_sigma_px);

            let (mut weight, mut mx, mut my) = (0.0, 0.0, 0.0);
            for (index, value) in ws.electrons.iter().enumerate() {
                let (col, row) = (index % cfg.width, index / cfg.width);
                weight += value;
                mx += value * col as f64;
                my += value * row as f64;
            }
            // The window keeps essentially all of the flux.
            assert!((weight - 1.0).abs() < 1e-8, "flux lost {}", 1.0 - weight);
            worst = worst
                .max((mx / weight - centre[0]).abs())
                .max((my / weight - centre[1]).abs());
        }
        assert!(worst < 1e-7, "worst renderer centroid bias {worst:e} px");
    }

    /// Phase 2 acceptance: 1000 random attitudes give a uniform boresight
    /// distribution, chi-square on an equal-area sky grid, p > 0.01.
    ///
    /// The grid uses 10 equal bands in `z` (equal `dz` is equal area on a
    /// sphere) times 10 equal right-ascension sectors, so all 100 cells have
    /// area `4*pi/100` and the expectation is 10 counts each. The critical
    /// value for 99 degrees of freedom at p = 0.01 is 134.6416, computed from
    /// the regularized incomplete gamma function.
    #[test]
    fn random_attitudes_point_uniformly_over_the_sky() {
        const BANDS: usize = 10;
        const SECTORS: usize = 10;
        const TRIALS: usize = 1000;
        const CHI2_CRITICAL: f64 = 134.6416;

        let mut counts = [0usize; BANDS * SECTORS];
        for trial in 0..TRIALS {
            let mut rng = ChaCha8Rng::seed_from_u64(crate::math::splitmix64(trial as u64));
            let a = random_attitude(&mut rng);
            // The inertial direction that maps onto the camera boresight
            // (+z) is the third row of A, since b = A r.
            let boresight = a.row(2).transpose();
            assert!((boresight.norm() - 1.0).abs() < 1e-12);

            let (ra, _) = crate::math::unit_to_radec(&boresight);
            let band = (((boresight.z + 1.0) / 2.0) * BANDS as f64) as usize;
            let sector = (ra / std::f64::consts::TAU * SECTORS as f64) as usize;
            counts[band.min(BANDS - 1) * SECTORS + sector.min(SECTORS - 1)] += 1;
        }

        let expected = TRIALS as f64 / (BANDS * SECTORS) as f64;
        let chi2: f64 = counts
            .iter()
            .map(|&observed| {
                let diff = observed as f64 - expected;
                diff * diff / expected
            })
            .sum();
        assert!(
            chi2 < CHI2_CRITICAL,
            "chi-square {chi2:.2} over {} cells exceeds the p = 0.01 critical \
             value {CHI2_CRITICAL}; boresights are not uniform",
            BANDS * SECTORS
        );
    }

    #[test]
    fn random_attitudes_are_proper_rotations() {
        for trial in 0..500u64 {
            let mut rng = ChaCha8Rng::seed_from_u64(trial);
            let a = random_attitude(&mut rng);
            assert!((a.determinant() - 1.0).abs() < 1e-12);
            assert!((a.transpose() * a - Mat3::identity()).norm() < 1e-12);
        }
    }

    // --- determinism ---

    /// Phase 2 acceptance: the same seed gives an identical image on repeated
    /// runs.
    #[test]
    fn the_same_seed_gives_an_identical_frame() {
        let catalog = catalog();
        for preset in Preset::ALL {
            let first = frame(4242, preset, &catalog);
            let second = frame(4242, preset, &catalog);
            assert_eq!(first.image, second.image, "{preset} image differs");
            assert_eq!(first.truth, second.truth, "{preset} truth differs");
            assert_eq!(first.false_stars, second.false_stars);
            assert_eq!(first.a_true, second.a_true);
        }
    }

    /// A reused workspace must not leak state between frames, or a benchmark
    /// would not reproduce a single trial run on its own.
    #[test]
    fn a_reused_workspace_gives_the_same_frame_as_a_fresh_one() {
        let catalog = catalog();
        let cfg = SimConfig::default();
        let mut shared = Workspace::new();

        // Warm the workspace on other seeds first.
        for seed in [1u64, 2, 3] {
            simulate(seed, &cfg, &catalog, &mut shared);
        }
        let reused = simulate(99, &cfg, &catalog, &mut shared);
        let fresh = simulate(99, &cfg, &catalog, &mut Workspace::new());
        assert_eq!(reused.image, fresh.image);
        assert_eq!(reused.truth, fresh.truth);
    }

    #[test]
    fn different_seeds_give_different_frames() {
        let catalog = catalog();
        let a = frame(1, Preset::Nominal, &catalog);
        let b = frame(2, Preset::Nominal, &catalog);
        assert_ne!(a.image, b.image);
        assert_ne!(a.a_true, b.a_true);
    }

    // --- frame contents ---

    #[test]
    fn frame_geometry_and_range_are_right() {
        let catalog = catalog();
        let cfg = SimConfig::default();
        let f = frame(7, Preset::Nominal, &catalog);
        assert_eq!(f.image.len(), cfg.width * cfg.height);
        assert_eq!((f.width, f.height), (cfg.width, cfg.height));
        assert!(f.image.iter().all(|&v| v <= cfg.saturation_adu()));
        // Hot pixels guarantee at least a few saturated entries.
        assert!(f.image.iter().any(|&v| v == cfg.saturation_adu()));
    }

    /// Every truth entry must sit on the sensor and reproject to its recorded
    /// pixel through the *true* camera and `A_true`.
    #[test]
    fn truth_positions_are_consistent_with_the_true_attitude() {
        let catalog = catalog();
        for preset in Preset::ALL {
            let cfg = SimConfig::preset(preset);
            let camera = cfg.true_camera();
            let f = frame(31337, preset, &catalog);
            assert!(!f.truth.is_empty(), "{preset} produced no visible stars");

            for star in &f.truth {
                assert!(camera.contains(star.pixel[0], star.pixel[1]));
                let entry = catalog
                    .stars
                    .iter()
                    .find(|s| s.id == star.id)
                    .expect("truth ids come from the catalogue");
                let expected = camera
                    .project(&(f.a_true * entry.direction()))
                    .expect("visible stars are in front");
                assert!((expected[0] - star.pixel[0]).abs() < 1e-9);
                assert!((expected[1] - star.pixel[1]).abs() < 1e-9);
                // Jitter is small, so the rendered magnitude stays close.
                assert!((f64::from(star.mag) - f64::from(entry.mag)).abs() < 1.0);
                assert!(f64::from(star.mag) <= f64::from(cfg.sim_mag_limit) + 1.0);
            }
        }
    }

    #[test]
    fn all_visible_stars_lie_within_the_half_diagonal_field() {
        let catalog = catalog();
        let cfg = SimConfig::default();
        let limit = cfg.true_camera().half_diagonal_fov();
        let f = frame(555, Preset::Nominal, &catalog);
        let boresight = f.a_true.row(2).transpose();
        for star in &f.truth {
            let entry = catalog.stars.iter().find(|s| s.id == star.id).expect("id");
            assert!(angle_between(&boresight, &entry.direction()) <= limit + 1e-12);
        }
    }

    #[test]
    fn easy_preset_drops_nothing_and_injects_nothing() {
        let catalog = catalog();
        for seed in 0..25u64 {
            let f = simulate(
                seed,
                &SimConfig::preset(Preset::Easy),
                &catalog,
                &mut Workspace::new(),
            );
            assert!(f.truth.iter().all(|s| !s.dropped));
            assert!(f.false_stars.is_empty());
        }
    }

    #[test]
    fn false_star_counts_stay_inside_the_preset_range() {
        let catalog = catalog();
        for preset in Preset::ALL {
            let [low, high] = SimConfig::preset(preset).false_star_range;
            let mut seen_low = false;
            let mut seen_high = false;
            for seed in 0..60u64 {
                let count = simulate(
                    seed,
                    &SimConfig::preset(preset),
                    &catalog,
                    &mut Workspace::new(),
                )
                .false_stars
                .len() as u32;
                assert!(
                    (low..=high).contains(&count),
                    "{preset} produced {count} false stars, outside {low}..={high}"
                );
                seen_low |= count == low;
                seen_high |= count == high;
            }
            // Over 60 seeds both ends of the range should appear.
            assert!(seen_low && seen_high, "{preset} never spanned its range");
        }
    }

    #[test]
    fn dropout_rate_is_close_to_the_configured_probability() {
        let catalog = catalog();
        let cfg = SimConfig::preset(Preset::Brutal);
        let (mut total, mut dropped) = (0usize, 0usize);
        for seed in 0..40u64 {
            let f = simulate(seed, &cfg, &catalog, &mut Workspace::new());
            total += f.truth.len();
            dropped += f.truth.iter().filter(|s| s.dropped).count();
        }
        let rate = dropped as f64 / total as f64;
        assert!(
            (rate - cfg.dropout_prob).abs() < 0.03,
            "dropout rate {rate:.4} against a configured {}",
            cfg.dropout_prob
        );
    }

    /// A dropped star must leave no trace in the image.
    #[test]
    fn dropped_stars_are_not_rendered() {
        let catalog = catalog();
        let mut cfg = SimConfig::preset(Preset::Easy);
        // Drop everything, and silence the other sources of signal.
        cfg.dropout_prob = 1.0;
        cfg.hot_pixels = 0;
        let dropped = simulate(17, &cfg, &catalog, &mut Workspace::new());
        assert!(dropped.truth.iter().all(|s| s.dropped));

        cfg.dropout_prob = 0.0;
        let kept = simulate(17, &cfg, &catalog, &mut Workspace::new());

        // With no stars rendered the frame is background only, so its
        // brightest pixel is far below a frame that does contain stars.
        let peak = |f: &SimFrame| *f.image.iter().max().unwrap_or(&0);
        assert!(
            peak(&dropped) < peak(&kept) / 2,
            "empty frame peaked at {} against {} with stars",
            peak(&dropped),
            peak(&kept)
        );
    }

    #[test]
    fn stars_actually_appear_where_truth_says_they_do() {
        let catalog = catalog();
        let mut cfg = SimConfig::preset(Preset::Easy);
        cfg.hot_pixels = 0;
        let f = simulate(2024, &cfg, &catalog, &mut Workspace::new());

        // Take the brightest rendered star and check the image is locally
        // bright at its position.
        let brightest = f
            .rendered()
            .min_by(|a, b| a.mag.total_cmp(&b.mag))
            .expect("some star was rendered");
        let col = brightest.pixel[0].round() as usize;
        let row = brightest.pixel[1].round() as usize;
        let here = f.image[row * f.width + col];

        let background_adu = cfg.background_e / cfg.gain_e_per_adu;
        assert!(
            f64::from(here) > 3.0 * background_adu,
            "pixel at the brightest star reads {here} ADU, background is {background_adu:.1}"
        );
    }

    #[test]
    fn the_bright_blob_saturates_a_wide_region() {
        let catalog = catalog();
        let mut cfg = SimConfig::preset(Preset::Brutal);
        cfg.hot_pixels = 0;
        let with_blob = simulate(8, &cfg, &catalog, &mut Workspace::new());
        cfg.bright_blob = false;
        let without = simulate(8, &cfg, &catalog, &mut Workspace::new());

        let saturated = |f: &SimFrame| {
            f.image
                .iter()
                .filter(|&&v| v == cfg.saturation_adu())
                .count()
        };
        // A wide saturated object, not a handful of saturated star cores.
        assert!(
            saturated(&with_blob) > saturated(&without) + 1000,
            "blob saturated {} pixels against {} without it",
            saturated(&with_blob),
            saturated(&without)
        );
    }

    #[test]
    fn noise_statistics_match_the_configuration() {
        let catalog = catalog();
        let mut cfg = SimConfig::preset(Preset::Easy);
        cfg.hot_pixels = 0;
        cfg.dropout_prob = 1.0; // background only, so statistics are clean
        let f = simulate(1234, &cfg, &catalog, &mut Workspace::new());

        let values: Vec<f64> = f.image.iter().map(|&v| f64::from(v)).collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance =
            values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / values.len() as f64;

        // Expected in ADU: background/gain, with variance
        // (background + read^2)/gain^2 plus 1/12 from rounding.
        let expected_mean = cfg.background_e / cfg.gain_e_per_adu;
        let expected_variance = (cfg.background_e + cfg.read_noise_e * cfg.read_noise_e)
            / (cfg.gain_e_per_adu * cfg.gain_e_per_adu)
            + 1.0 / 12.0;
        assert!(
            (mean - expected_mean).abs() < 0.1,
            "mean {mean:.3} ADU against an expected {expected_mean:.3}"
        );
        assert!(
            (variance / expected_variance - 1.0).abs() < 0.05,
            "variance {variance:.3} ADU^2 against an expected {expected_variance:.3}"
        );
    }

    #[test]
    fn poisson_branch_has_the_right_mean_and_variance() {
        // The small-mean branch must be an honest Poisson: mean = variance.
        let mut rng = ChaCha8Rng::seed_from_u64(0x5EED_0155);
        for mean in [0.5, 5.0, 20.0] {
            let samples: Vec<f64> = (0..20_000).map(|_| poisson_small(mean, &mut rng)).collect();
            let m = samples.iter().sum::<f64>() / samples.len() as f64;
            let v = samples.iter().map(|s| (s - m) * (s - m)).sum::<f64>() / samples.len() as f64;
            assert!((m / mean - 1.0).abs() < 0.05, "mean {m} against {mean}");
            assert!((v / mean - 1.0).abs() < 0.08, "variance {v} against {mean}");
            assert!(samples.iter().all(|s| s.fract() == 0.0 && *s >= 0.0));
        }
    }

    #[test]
    fn hot_pixel_count_is_respected() {
        let catalog = catalog();
        let mut cfg = SimConfig::preset(Preset::Easy);
        cfg.dropout_prob = 1.0;
        cfg.background_e = 0.0;
        cfg.read_noise_e = 0.0;
        cfg.hot_pixels = 7;
        let f = simulate(3, &cfg, &catalog, &mut Workspace::new());
        let saturated = f
            .image
            .iter()
            .filter(|&&v| v == cfg.saturation_adu())
            .count();
        // Distinct positions are overwhelmingly likely over 1M pixels.
        assert_eq!(saturated, 7);
    }
}
