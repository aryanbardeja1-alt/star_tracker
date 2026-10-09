//! Attitude determination: Wahba's problem solved by SVD.
//!
//! The rotation `A` maps inertial (ICRS/J2000) directions into the camera
//! frame: `b = A r`.
//!
//! [`solve_wahba`] is the solver itself, taking weighted direction pairs.
//! [`solve_weighted`] is the pipeline stage on top of it, turning identified
//! stars and their centroids into an attitude with `w = 1/sigma^2` weights
//! derived from brightness. [`solve_davenport`] is an independent route to the
//! same answer, kept as a cross-check.

use crate::camera::CameraModel;
use crate::catalog::Catalog;
use crate::centroid::Centroid;
use crate::identify::{self, StarMatch};
use crate::math::{Mat3, Vec3, quat_to_matrix};
use crate::simulate::SimConfig;
use nalgebra::Matrix4;

/// Solves Wahba's problem for the rotation `A` with `b = A r`.
///
/// `pairs` holds `(body vector b, inertial vector r, weight)`, all vectors unit
/// length and both in their own frame. Builds `B = sum w b r^T`, takes the SVD
/// `B = U S V^T` and returns `A = U diag(1, 1, det(U) det(V)) V^T`.
///
/// The determinant correction is mandatory, not cosmetic: without it the SVD
/// can return a reflection instead of a rotation, which still looks like a
/// plausible attitude and is wrong.
///
/// Returns `None` for fewer than two directions -- one direction cannot fix
/// roll -- or if the SVD fails to converge.
pub fn solve_wahba(pairs: &[(Vec3, Vec3, f64)]) -> Option<Mat3> {
    if pairs.len() < 2 {
        return None;
    }

    let mut profile = Mat3::zeros();
    for (body, inertial, weight) in pairs {
        profile += (body * inertial.transpose()) * *weight;
    }

    let svd = profile.svd(true, true);
    let u = svd.u?;
    let v_t = svd.v_t?;

    // det(V) equals det(V^T), so this is the det(U) det(V) of the formula.
    let sign = (u.determinant() * v_t.determinant()).signum();
    let correction = Mat3::from_diagonal(&Vec3::new(1.0, 1.0, sign));
    Some(u * correction * v_t)
}

/// Total flux, in ADU, at which the point-spread peak reaches full well.
///
/// A Gaussian places a fraction `1 / (2*pi*sigma^2)` of its flux in the
/// central pixel, so the peak saturates once the total passes this.
fn saturation_flux(cfg: &SimConfig) -> f64 {
    let peak_fraction = 1.0 / (std::f64::consts::TAU * cfg.psf_sigma_px * cfg.psf_sigma_px);
    f64::from(cfg.saturation_adu()) / peak_fraction
}

/// Weight for a centroid of the given background-subtracted flux, `w = 1/sigma^2`.
///
/// A shot-noise-limited centroid has `sigma ~ sigma_psf / sqrt(N)`, so
/// `w = 1/sigma^2` is proportional to the flux itself; the constant of
/// proportionality cancels, because scaling every weight scales `B` uniformly
/// and leaves the solution untouched. Measured across magnitude 3 to 6.5,
/// `sigma * sqrt(flux)` holds to within about 1.5x, which is the law this
/// relies on.
///
/// The flux is capped where the peak saturates. Past that point the measured
/// flux keeps rising while the centroid stops improving -- a magnitude 1 star
/// measured 12x the flux of a magnitude 4 one but scattered just as widely --
/// so an uncapped weight would trust it an order of magnitude too much.
pub fn centroid_weight(flux: f64, cfg: &SimConfig) -> f64 {
    flux.clamp(0.0, saturation_flux(cfg))
}

/// The attitude stage: a brightness-weighted Wahba solution from identified
/// stars.
///
/// Observed directions come from un-projecting each centroid through the
/// *nominal* camera, which is the only one the solver knows. Returns `None` if
/// fewer than two matches survive.
///
/// Against uniform weights, measured over 1498 nominal frames, this improves
/// the median attitude error from 1.59 to 1.20 arcsec and the 95th percentile
/// from 5.74 to 4.71.
pub fn solve_weighted(
    matches: &[StarMatch],
    centroids: &[Centroid],
    catalog: &Catalog,
    cfg: &SimConfig,
) -> Option<Mat3> {
    let camera = cfg.nominal_camera();
    let pairs: Vec<(Vec3, Vec3, f64)> = matches
        .iter()
        .filter_map(|matched| {
            let centroid = centroids.get(usize::from(matched.observed))?;
            let star = catalog.stars.get(usize::from(matched.catalog_index))?;
            Some((
                camera.unproject(centroid.x, centroid.y),
                star.direction(),
                centroid_weight(centroid.flux, cfg),
            ))
        })
        .collect();
    solve_wahba(&pairs)
}

/// Solves Wahba's problem by Davenport's q-method, as a cross-check on
/// [`solve_wahba`].
///
/// Builds Davenport's `K` matrix from the attitude profile and takes the
/// eigenvector of its largest eigenvalue, which is the optimal quaternion.
///
/// This is the exact eigensolution rather than QUEST. QUEST is the same method
/// with the largest eigenvalue approximated by Newton iteration on the
/// characteristic quartic -- faster, but carrying a non-convergence case and
/// needing the method of sequential rotations near 180 degrees. docs/SPEC.md wants
/// this "as a cross-check only", and for that purpose an exact, independent
/// solution with no failure modes of its own is worth more than a fast one.
pub fn solve_davenport(pairs: &[(Vec3, Vec3, f64)]) -> Option<Mat3> {
    if pairs.len() < 2 {
        return None;
    }

    let mut profile = Mat3::zeros();
    let mut z = Vec3::zeros();
    for (body, inertial, weight) in pairs {
        profile += (body * inertial.transpose()) * *weight;
        // `r x b`, not `b x r`. The opposite sign conjugates the quaternion,
        // which agrees with the SVD solution on symmetric configurations and
        // is 180 degrees out on general ones -- the error this cross-check
        // exists to catch, and did.
        z += inertial.cross(body) * *weight;
    }
    let sigma = profile.trace();
    let s = profile + profile.transpose();

    let mut k = Matrix4::zeros();
    k.fixed_view_mut::<3, 3>(0, 0)
        .copy_from(&(s - Mat3::identity() * sigma));
    k.fixed_view_mut::<3, 1>(0, 3).copy_from(&z);
    k.fixed_view_mut::<1, 3>(3, 0).copy_from(&z.transpose());
    k[(3, 3)] = sigma;

    let eigen = k.symmetric_eigen();
    let mut largest = 0;
    for index in 1..4 {
        if eigen.eigenvalues[index] > eigen.eigenvalues[largest] {
            largest = index;
        }
    }
    let v = eigen.eigenvectors.column(largest);
    // Davenport's eigenvector is vector-first; the crate's convention is
    // scalar-first.
    Some(quat_to_matrix([v[3], v[0], v[1], v[2]]))
}

/// Iterations of the focal-length fit. The pure scale converges in one; the
/// extra passes settle the residual that radial distortion leaves behind.
const FOCAL_ITERATIONS: usize = 3;

/// How much wider the bootstrap sweeps run than the final one.
///
/// While the focal length is still being fitted the outer stars are displaced
/// past the normal tolerance, and those are precisely the stars the fit needs.
const BOOTSTRAP_TOLERANCE_FACTOR: f64 = 10.0;

/// Matches closer to the boresight than this normalised radius are skipped
/// when fitting the focal length: the radius ratio is the measurement, and
/// near the axis it is a ratio of two small, noisy numbers.
const MIN_FIT_RADIUS: f64 = 0.02;

/// Largest focal-length correction the fit is allowed to apply, as a fraction
/// of the nominal value.
///
/// Ground calibration pins a real camera's focal length far better than this,
/// so a fit that wants more is not measuring a focal length -- it is absorbing
/// a bad identification. Left unbounded it reached 0.56 and 1.34 on `hard`
/// frames whose matches were wrong, and fitting a free parameter to bad data
/// made those solutions look self-consistent enough to be accepted.
const MAX_FOCAL_CORRECTION: f64 = 0.02;

/// An attitude together with the focal length the frame was actually taken at.
#[derive(Clone, Debug, PartialEq)]
pub struct Calibrated {
    /// Rotation with `b = A r`.
    pub attitude: Mat3,
    /// The camera as fitted: nominal in every respect but the focal length.
    pub camera: CameraModel,
    /// Fitted focal length over the nominal one. Exactly 1 when nothing was
    /// fitted.
    pub focal_scale: f64,
    /// The match list, which the re-sweeps may have grown beyond what
    /// identification first found.
    pub matches: Vec<StarMatch>,
}

/// Solves for attitude *and* focal length from identified stars.
///
/// The solver un-projects with the nominal focal length. Where the true one
/// differs by a factor `s`, every reconstructed direction sits further from the
/// boresight than it should by that same factor, because the normalised image
/// radius scales as `f_true / f_nominal`. So comparing each matched star's
/// measured radius with the radius its catalogue direction predicts recovers
/// `s` directly; the median over the matches makes that robust to one bad
/// match, and a few passes settle the part of the residual that distortion,
/// which is not a pure scale, leaves behind.
///
/// This matters because an uncalibrated focal length is not a small error. At
/// the `hard` preset's 0.2% it perturbs a 20 degree pair angle by 163 arcsec,
/// and leaves the attitude 54 arcsec out even when every identification is
/// right. Fitting it brings that to a few arcseconds.
pub fn solve_calibrated(
    matches: &[StarMatch],
    centroids: &[Centroid],
    catalog: &Catalog,
    db_star_count: usize,
    cfg: &SimConfig,
    ws: &mut identify::Workspace,
) -> Option<Calibrated> {
    let mut camera = cfg.nominal_camera();
    let nominal_focal = camera.f;
    let mut matches = matches.to_vec();
    let mut attitude = solve_at_focal(&matches, centroids, catalog, cfg, &camera)?;

    let tolerance = identify::Tolerance::from_config(cfg);
    let mut directions = Vec::with_capacity(centroids.len());
    let mut ratios = Vec::with_capacity(centroids.len());

    for _ in 0..FOCAL_ITERATIONS {
        ratios.clear();
        for matched in &matches {
            let Some(centroid) = centroids.get(usize::from(matched.observed)) else {
                continue;
            };
            let Some(star) = catalog.stars.get(usize::from(matched.catalog_index)) else {
                continue;
            };
            // Measured radius, in the normalised coordinates the current focal
            // length implies.
            let measured = ((centroid.x - camera.cx).hypot(centroid.y - camera.cy)) / camera.f;
            // Radius the catalogue direction predicts through this attitude.
            let predicted_direction = attitude * star.direction();
            if predicted_direction.z <= 0.0 {
                continue;
            }
            let predicted = (predicted_direction.x / predicted_direction.z)
                .hypot(predicted_direction.y / predicted_direction.z);
            if predicted < MIN_FIT_RADIUS {
                continue;
            }
            ratios.push(measured / predicted);
        }
        if ratios.len() < 2 {
            break;
        }

        ratios.sort_by(f64::total_cmp);
        let scale = ratios[ratios.len() / 2];
        // A measured radius larger than predicted means the focal length in
        // use is too short, so it grows by the same ratio -- but only within
        // what the camera's calibration could plausibly be out by.
        camera.f = (camera.f * scale).clamp(
            nominal_focal * (1.0 - MAX_FOCAL_CORRECTION),
            nominal_focal * (1.0 + MAX_FOCAL_CORRECTION),
        );
        attitude = solve_at_focal(&matches, centroids, catalog, cfg, &camera)?;

        // With a better camera the sweep reaches further out, so run it again:
        // the extra stars are what make the next fit well conditioned. This is
        // the half of the loop that breaks the deadlock between fit and sweep.
        //
        // The tolerance here is deliberately generous. While the camera is
        // still uncalibrated the outer stars sit well outside the normal
        // tolerance -- a residual 0.1% scale error displaces a star at the
        // field edge by 47 arcsec -- so holding the bootstrap to that tolerance
        // leaves the loop stuck exactly where it needs to reach. A wrong match
        // at this width is still very unlikely (a 384 arcsec circle holds 0.004
        // database stars on average) and the ambiguity check still applies; the
        // final pass below re-derives the match list at full strictness anyway.
        unproject_all(centroids, &camera, &mut directions);
        identify::extend_matches(
            &directions,
            &attitude,
            &camera,
            catalog,
            db_star_count,
            tolerance.vector * BOOTSTRAP_TOLERANCE_FACTOR,
            &mut matches,
            ws,
        );
        attitude = solve_at_focal(&matches, centroids, catalog, cfg, &camera)?;
    }

    // Final pass at the proper tolerance: rebuild the match list from the
    // pyramid alone, so nothing the calibrated camera cannot actually account
    // for survives into the reported solution.
    matches.retain(|matched| matched.from_pyramid);
    unproject_all(centroids, &camera, &mut directions);
    identify::extend_matches(
        &directions,
        &attitude,
        &camera,
        catalog,
        db_star_count,
        tolerance.vector,
        &mut matches,
        ws,
    );
    attitude = solve_at_focal(&matches, centroids, catalog, cfg, &camera)?;

    Some(Calibrated {
        attitude,
        camera,
        focal_scale: camera.f / nominal_focal,
        matches,
    })
}

/// Un-projects every centroid through `camera` into `out`.
fn unproject_all(centroids: &[Centroid], camera: &CameraModel, out: &mut Vec<Vec3>) {
    out.clear();
    out.extend(
        centroids
            .iter()
            .map(|centroid| camera.unproject(centroid.x, centroid.y)),
    );
}

/// A weighted Wahba solution with the centroids un-projected through `camera`.
fn solve_at_focal(
    matches: &[StarMatch],
    centroids: &[Centroid],
    catalog: &Catalog,
    cfg: &SimConfig,
    camera: &CameraModel,
) -> Option<Mat3> {
    let pairs: Vec<(Vec3, Vec3, f64)> = matches
        .iter()
        .filter_map(|matched| {
            let centroid = centroids.get(usize::from(matched.observed))?;
            let star = catalog.stars.get(usize::from(matched.catalog_index))?;
            Some((
                camera.unproject(centroid.x, centroid.y),
                star.direction(),
                centroid_weight(centroid.flux, cfg),
            ))
        })
        .collect();
    solve_wahba(&pairs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{angle_between, rotation_vector, splitmix64};
    use crate::simulate::Preset;
    use proptest::prelude::*;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use rand_distr::StandardNormal;
    use std::f64::consts::{PI, TAU};

    /// One arcsecond in radians.
    const ARCSEC: f64 = PI / (180.0 * 3600.0);

    fn random_attitude(rng: &mut ChaCha8Rng) -> Mat3 {
        quat_to_matrix([
            rng.sample::<f64, _>(StandardNormal),
            rng.sample::<f64, _>(StandardNormal),
            rng.sample::<f64, _>(StandardNormal),
            rng.sample::<f64, _>(StandardNormal),
        ])
    }

    fn random_unit(rng: &mut ChaCha8Rng) -> Vec3 {
        let ra = rng.random::<f64>() * TAU;
        let z = rng.random::<f64>() * 2.0 - 1.0;
        let rho = (1.0 - z * z).sqrt();
        Vec3::new(rho * ra.cos(), rho * ra.sin(), z)
    }

    /// Perturbs a unit vector by `sigma` radians in each tangential direction.
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

    /// `n` observations of a random attitude, spread over the whole sphere.
    fn observations(
        seed: u64,
        n: usize,
        sigma: f64,
        weighted: bool,
    ) -> (Mat3, Vec<(Vec3, Vec3, f64)>) {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let attitude = random_attitude(&mut rng);
        let pairs = (0..n)
            .map(|_| {
                let inertial = random_unit(&mut rng);
                let body = jitter(&(attitude * inertial), sigma, &mut rng);
                let weight = if weighted {
                    0.1 + rng.random::<f64>() * 10.0
                } else {
                    1.0
                };
                (body, inertial, weight)
            })
            .collect();
        (attitude, pairs)
    }

    // --- property tests, as docs/SPEC.md requires ---

    proptest! {
        /// Phase 6 acceptance: noiseless recovery error below 1e-9 rad.
        #[test]
        fn noiseless_recovery_is_exact(seed: u64, n in 2usize..24, weighted: bool) {
            let (truth, pairs) = observations(seed, n, 0.0, weighted);
            let solved = solve_wahba(&pairs).expect("two or more directions always solve");
            let error = rotation_vector(&(solved * truth.transpose())).norm();
            prop_assert!(
                error < 1e-9,
                "recovered attitude is {error:e} rad from the truth with {n} stars"
            );
        }

        /// Phase 6 acceptance: `det = +1` always.
        #[test]
        fn the_result_is_always_a_proper_rotation(seed: u64, n in 2usize..24, sigma in 0.0..0.2) {
            let (_, pairs) = observations(seed, n, sigma, true);
            let solved = solve_wahba(&pairs).expect("solves");
            prop_assert!((solved.determinant() - 1.0).abs() < 1e-12, "det = {}", solved.determinant());
            let orthogonality = (solved.transpose() * solved - Mat3::identity()).norm();
            prop_assert!(orthogonality < 1e-12, "not orthonormal: {orthogonality:e}");
        }

        /// Weights scale `B` uniformly when they are all equal, and cannot move
        /// a solution that already fits exactly.
        #[test]
        fn weights_do_not_disturb_an_exact_fit(seed: u64, n in 2usize..16, scale in 0.01..100.0) {
            let (truth, pairs) = observations(seed, n, 0.0, false);
            let scaled: Vec<_> = pairs.iter().map(|(b, r, w)| (*b, *r, w * scale)).collect();
            let solved = solve_wahba(&scaled).expect("solves");
            prop_assert!(rotation_vector(&(solved * truth.transpose())).norm() < 1e-9);
        }

        /// Davenport's q-method is an independent route to the same answer.
        ///
        /// "The same answer" is the same minimum of Wahba's cost, which is not
        /// always the same rotation. When every direction lies within a narrow
        /// cone the cost is nearly flat about the axis through them, and two
        /// solvers can settle a nanoradian apart while fitting the data equally
        /// well. Measured: two stars 1.76 degrees apart put the solvers 1.0e-9
        /// rad apart, the same pair at 6.8 degrees 4.3e-12, at 144.7 degrees
        /// 4.1e-15, and a third star takes it to 1.3e-15. So the cost is
        /// compared on every input, and the rotations only where the geometry
        /// determines one.
        #[test]
        fn davenport_agrees_with_the_svd_solution(seed: u64, n in 2usize..20, sigma in 0.0..0.05) {
            let (_, pairs) = observations(seed, n, sigma, true);
            let svd = solve_wahba(&pairs).expect("solves");
            let davenport = solve_davenport(&pairs).expect("solves");

            // Both claim to minimise the same sum; disagreeing about its value
            // would mean one of them is not at the minimum, whatever the
            // geometry. This is the part that holds unconditionally.
            let cost = |a: &Mat3| -> f64 {
                pairs.iter().map(|(b, r, w)| w * (b - a * r).norm_squared()).sum()
            };
            let (left, right) = (cost(&svd), cost(&davenport));
            prop_assert!(
                (left - right).abs() <= 1e-12 * (1.0 + left.abs()),
                "the solvers disagree about the cost: {left:e} against {right:e}"
            );

            // How well the directions pin a rotation down at all: two stars a
            // degree apart leave the roll about them almost free.
            let spread = pairs
                .iter()
                .flat_map(|(a, _, _)| pairs.iter().map(move |(b, _, _)| angle_between(a, b)))
                .fold(0.0_f64, f64::max);
            if spread > 10.0_f64.to_radians() {
                let difference = rotation_vector(&(svd * davenport.transpose())).norm();
                prop_assert!(
                    difference < 1e-9,
                    "the two solvers differ by {difference:e} rad with {n} stars                      spread over {:.1} degrees",
                    spread.to_degrees()
                );
            }
        }
    }

    // --- the determinant correction ---

    /// The correction is mandatory, not cosmetic. Body vectors that are a
    /// *reflection* of the inertial ones make the best orthogonal fit improper,
    /// so `det(U) det(V)` is -1 and the raw `U V^T` is a reflection. Measured
    /// on this configuration: `det(B) = -8.97` and `det(U V^T) = -1`.
    #[test]
    fn the_determinant_correction_turns_a_reflection_into_a_rotation() {
        let mut rng = ChaCha8Rng::seed_from_u64(99);
        let mirror = Mat3::from_diagonal(&Vec3::new(1.0, 1.0, -1.0));
        let pairs: Vec<(Vec3, Vec3, f64)> = (0..8)
            .map(|_| {
                let inertial = random_unit(&mut rng);
                (mirror * inertial, inertial, 1.0)
            })
            .collect();

        // Confirm the configuration really does provoke the correction.
        let mut profile = Mat3::zeros();
        for (body, inertial, weight) in &pairs {
            profile += (body * inertial.transpose()) * *weight;
        }
        let svd = profile.svd(true, true);
        let (u, v_t) = (svd.u.expect("u"), svd.v_t.expect("v_t"));
        assert!(profile.determinant() < 0.0);
        assert!((u.determinant() * v_t.determinant()) < 0.0);
        assert!(
            (u * v_t).determinant() < 0.0,
            "the uncorrected product should be a reflection"
        );

        // Both solvers must nonetheless return a rotation.
        for (label, solved) in [
            ("svd", solve_wahba(&pairs).expect("solves")),
            ("davenport", solve_davenport(&pairs).expect("solves")),
        ] {
            assert!(
                (solved.determinant() - 1.0).abs() < 1e-12,
                "{label} returned det {}",
                solved.determinant()
            );
        }
    }

    #[test]
    fn two_directions_are_enough_and_one_is_not() {
        let (truth, pairs) = observations(7, 2, 0.0, false);
        let solved = solve_wahba(&pairs).expect("two directions fix an attitude");
        assert!(rotation_vector(&(solved * truth.transpose())).norm() < 1e-9);
        assert!((solved.determinant() - 1.0).abs() < 1e-12);

        // One direction leaves roll free, so there is nothing to return.
        assert!(solve_wahba(&pairs[..1]).is_none());
        assert!(solve_wahba(&[]).is_none());
        assert!(solve_davenport(&pairs[..1]).is_none());
    }

    // --- error scaling under noise ---

    /// Phase 6 acceptance: the error under noise matches the predicted
    /// `sigma / sqrt(N)` scaling within a factor of two.
    ///
    /// The exact prediction for directions spread over the sphere is
    /// `sigma * sqrt(3 / (2N))` per axis, so the ratio to a bare
    /// `sigma / sqrt(N)` is `sqrt(1.5) = 1.225`; the total over three axes is
    /// `3 / sqrt(2) = 2.121` times it. Measured over 4000 repetitions per N
    /// those ratios came out at 1.21 to 1.46 and 2.10 to 2.53, approaching the
    /// analytic values as N grows. This asserts both the 1/sqrt(N) shape and
    /// the coefficient.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "statistical, run with --release")]
    fn error_under_noise_follows_the_predicted_scaling() {
        const SIGMA: f64 = 1e-4;
        const REPEATS: u64 = 1500;

        let mut scaled = Vec::new();
        for n in [4usize, 8, 16, 32, 64] {
            let mut sum_squares = 0.0;
            for repeat in 0..REPEATS {
                let seed = splitmix64((n as u64) << 40 | repeat);
                let (truth, pairs) = observations(seed, n, SIGMA, false);
                let solved = solve_wahba(&pairs).expect("solves");
                let phi = rotation_vector(&(solved * truth.transpose()));
                // Per-axis, so the comparison is against sigma/sqrt(N).
                sum_squares += phi.norm_squared() / 3.0;
            }
            let rms = (sum_squares / REPEATS as f64).sqrt();
            let predicted = SIGMA / (n as f64).sqrt();
            let ratio = rms / predicted;
            assert!(
                (0.5..2.0).contains(&ratio),
                "N = {n}: per-axis RMS {rms:e} against a predicted {predicted:e} (ratio {ratio:.3})"
            );
            scaled.push(rms * (n as f64).sqrt());
        }

        // The 1/sqrt(N) shape itself: rms * sqrt(N) must be flat.
        let smallest = scaled.iter().copied().fold(f64::INFINITY, f64::min);
        let largest = scaled.iter().copied().fold(0.0, f64::max);
        assert!(
            largest / smallest < 1.25,
            "rms * sqrt(N) drifted from {smallest:e} to {largest:e}, so the scaling is not 1/sqrt(N)"
        );
    }

    // --- brightness weighting ---

    #[test]
    fn centroid_weight_follows_flux_until_the_peak_saturates() {
        let cfg = SimConfig::preset(Preset::Nominal);
        let cap = saturation_flux(&cfg);
        // The default camera saturates a Gaussian peak at about 25,700 ADU of
        // total flux.
        assert!(
            (25_000.0..26_500.0).contains(&cap),
            "saturation flux {cap:.0} ADU"
        );

        // Below the cap the weight is the flux itself.
        for flux in [1.0, 100.0, 5_000.0, 20_000.0] {
            assert_eq!(centroid_weight(flux, &cfg), flux);
        }
        // Above it the weight stops growing.
        assert_eq!(centroid_weight(cap * 10.0, &cfg), cap);
        // And a negative flux cannot produce a negative weight.
        assert_eq!(centroid_weight(-5.0, &cfg), 0.0);
    }

    /// Weighting by brightness must help when the directions really do differ
    /// in quality: noise is scaled down for the stars given more weight.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "statistical, run with --release")]
    fn brightness_weighting_beats_uniform_weighting() {
        const REPEATS: u64 = 2000;
        let mut weighted_sum = 0.0;
        let mut uniform_sum = 0.0;

        for repeat in 0..REPEATS {
            let mut rng = ChaCha8Rng::seed_from_u64(splitmix64(0x5151 ^ repeat));
            let truth = random_attitude(&mut rng);

            // Half the stars are ten times better than the other half, so the
            // optimal weights differ by a hundred.
            let mut pairs = Vec::new();
            for index in 0..12 {
                let good = index % 2 == 0;
                let sigma = if good { 1e-5 } else { 1e-4 };
                let inertial = random_unit(&mut rng);
                let body = jitter(&(truth * inertial), sigma, &mut rng);
                pairs.push((body, inertial, 1.0 / (sigma * sigma)));
            }

            let weighted = solve_wahba(&pairs).expect("solves");
            let uniform: Vec<_> = pairs.iter().map(|(b, r, _)| (*b, *r, 1.0)).collect();
            let uniform = solve_wahba(&uniform).expect("solves");

            weighted_sum += rotation_vector(&(weighted * truth.transpose())).norm_squared();
            uniform_sum += rotation_vector(&(uniform * truth.transpose())).norm_squared();
        }

        let weighted_rms = (weighted_sum / REPEATS as f64).sqrt();
        let uniform_rms = (uniform_sum / REPEATS as f64).sqrt();
        assert!(
            weighted_rms < 0.5 * uniform_rms,
            "weighted RMS {:.4} arcsec against uniform {:.4}; weighting bought nothing",
            weighted_rms / ARCSEC,
            uniform_rms / ARCSEC
        );
    }
}
