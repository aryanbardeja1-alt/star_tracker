//! Vector, quaternion and rotation helpers.
//!
//! Every angle here is in radians and every float is `f64`, per docs/SPEC.md;
//! conversion to degrees or arcseconds happens only at the UI/CLI boundary.
//!
//! A rotation matrix called `A` maps inertial (ICRS/J2000) vectors into the
//! camera frame: `b = A r`. Quaternions are Hamilton convention, scalar-first
//! `[w, x, y, z]`, normalised with `w >= 0`.

use nalgebra::{Matrix3, Quaternion, Rotation3, UnitQuaternion, Vector3};
use std::f64::consts::TAU;

/// A 3-vector. Each function's doc comment names the frame it belongs to.
pub type Vec3 = Vector3<f64>;

/// A 3x3 matrix. Rotation matrices follow the `b = A r` convention above.
pub type Mat3 = Matrix3<f64>;

/// Angle between two vectors, in radians, in `[0, pi]`.
///
/// Frame-agnostic: both vectors must simply be in the same frame. Computed as
/// `atan2(|a x b|, a . b)`, which stays accurate for small angles where
/// `acos(a . b)` loses most of its significant digits. The inputs need not be
/// normalised.
pub fn angle_between(a: &Vec3, b: &Vec3) -> f64 {
    a.cross(b).norm().atan2(a.dot(b))
}

/// Inertial (ICRS/J2000) unit vector from right ascension and declination, radians.
///
/// `r = (cos d cos a, cos d sin a, sin d)`.
pub fn radec_to_unit(ra: f64, dec: f64) -> Vec3 {
    let (sin_ra, cos_ra) = ra.sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();
    Vec3::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
}

/// Right ascension and declination, radians, from an inertial (ICRS/J2000) vector.
///
/// Returns `(ra, dec)` with `ra` wrapped into `[0, TAU)` and `dec` in
/// `[-pi/2, pi/2]`. The input need not be normalised.
pub fn unit_to_radec(v: &Vec3) -> (f64, f64) {
    let mut ra = v.y.atan2(v.x);
    if ra < 0.0 {
        ra += TAU;
    }
    // A tiny negative angle can round up to exactly TAU when TAU is added.
    if ra >= TAU {
        ra = 0.0;
    }
    let dec = v.z.atan2(v.x.hypot(v.y));
    (ra, dec)
}

/// Rotation matrix from a Hamilton, scalar-first quaternion `[w, x, y, z]`.
///
/// The quaternion is normalised first, so an unnormalised input is fine; a
/// zero quaternion is not a rotation and yields NaN. The result `A` encodes
/// the same rotation, i.e. `b = A r` for the quaternion's own frames.
pub fn quat_to_matrix(q: [f64; 4]) -> Mat3 {
    let unit = UnitQuaternion::from_quaternion(Quaternion::new(q[0], q[1], q[2], q[3]));
    unit.to_rotation_matrix().into_inner()
}

/// Hamilton, scalar-first quaternion `[w, x, y, z]` from a rotation matrix.
///
/// `a` must be a proper rotation (`det = +1`, orthonormal). The result is
/// normalised and has `w >= 0`, which picks one of the two quaternions
/// encoding the rotation so that comparisons and serialisation are
/// unambiguous.
pub fn matrix_to_quat(a: &Mat3) -> [f64; 4] {
    let unit = UnitQuaternion::from_rotation_matrix(&Rotation3::from_matrix_unchecked(*a));
    let q = unit.into_inner();
    let q = if q.w < 0.0 { -q } else { q };
    [q.w, q.i, q.j, q.k]
}

/// Rotation vector (axis times angle, radians) of a rotation matrix.
///
/// The matrix logarithm: `|phi|` is the rotation angle in `[0, pi]` and `phi`
/// points along the rotation axis, expressed in the frame `a` rotates into.
/// Routed through the quaternion so it keeps full precision at the
/// arcsecond-scale angles that attitude-error reporting needs.
///
/// For the attitude error of docs/SPEC.md, pass `A_est * A_true.transpose()`; the
/// result is then a rotation vector in the camera frame, whose components are
/// cross-boresight (`x`, `y`) and roll (`z`).
pub fn rotation_vector(a: &Mat3) -> Vec3 {
    let q = matrix_to_quat(a);
    let axis = Vec3::new(q[1], q[2], q[3]);
    let sin_half = axis.norm();
    if sin_half == 0.0 {
        return Vec3::zeros();
    }
    // 2 * atan2(s, w) / s stays well conditioned as s -> 0, tending to 2 / w.
    axis * (2.0 * sin_half.atan2(q[0]) / sin_half)
}

/// SplitMix64 bit mixer.
///
/// Pure wrapping `u64` arithmetic, so native and WASM agree bit for bit. The
/// seed for trial `i` of a benchmark is `splitmix64(base_seed ^ i)`.
pub fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use std::f64::consts::{FRAC_PI_2, PI};

    /// One arcsecond in radians -- the scale attitude error is reported at.
    const ARCSEC: f64 = PI / (180.0 * 3600.0);

    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(0xA5A5_0F0F_1234_5678)
    }

    /// Uniform point on the sphere, from a uniform azimuth and a uniform `z`.
    fn random_unit(r: &mut ChaCha8Rng) -> Vec3 {
        let ra = r.random::<f64>() * TAU;
        let z = r.random::<f64>() * 2.0 - 1.0;
        let rho = (1.0 - z * z).sqrt();
        Vec3::new(rho * ra.cos(), rho * ra.sin(), z)
    }

    /// Rotation about a uniform axis by an angle in `[0, max_angle)`.
    fn random_rotation(r: &mut ChaCha8Rng, max_angle: f64) -> Mat3 {
        let axis = random_unit(r);
        let angle = r.random::<f64>() * max_angle;
        UnitQuaternion::from_scaled_axis(axis * angle)
            .to_rotation_matrix()
            .into_inner()
    }

    #[test]
    fn angle_between_known_values() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert!((angle_between(&x, &y) - FRAC_PI_2).abs() < 1e-15);
        assert!(angle_between(&x, &x).abs() < 1e-15);
        assert!((angle_between(&x, &-x) - PI).abs() < 1e-15);
    }

    #[test]
    fn angle_between_ignores_magnitude() {
        let mut r = rng();
        for _ in 0..200 {
            let a = random_unit(&mut r);
            let b = random_unit(&mut r);
            let scaled = angle_between(&(a * 1234.5), &(b * 1e-6));
            assert!((scaled - angle_between(&a, &b)).abs() < 1e-14);
        }
    }

    /// Pitfall 3 of docs/SPEC.md, pinned. The atan2 form holds a flat absolute
    /// error of about 1.1e-16 rad at every angle, because what limits it is
    /// the storage of `b` as a unit vector rather than the formula. The
    /// `acos(a . b)` form degrades as the angle shrinks -- measured at 3e-8 rad
    /// for a 1e-9 rad separation, i.e. 30 times larger than the angle itself --
    /// so no relative bound is achievable with it.
    #[test]
    fn angle_between_is_accurate_at_tiny_angles() {
        let mut r = rng();
        for exponent in 0..13 {
            let theta = 10f64.powi(-exponent);
            let a = random_unit(&mut r);
            // Turn `a` about an axis perpendicular to it, so the angle swept
            // is exactly `theta`.
            let perp = random_unit(&mut r).cross(&a).normalize();
            let b = UnitQuaternion::from_scaled_axis(perp * theta).transform_vector(&a);

            let err = (angle_between(&a, &b) - theta).abs();
            assert!(err <= 1e-15, "theta={theta:e} absolute error={err:e}");

            // The same quantity via acos, to keep the contrast measured rather
            // than asserted in a comment.
            let acos_err = (a.dot(&b).clamp(-1.0, 1.0).acos() - theta).abs();
            if theta <= 1e-6 {
                assert!(
                    acos_err > 100.0 * err,
                    "acos should be far worse at theta={theta:e}: \
                     atan2={err:e} acos={acos_err:e}"
                );
            }
        }
    }
    #[test]
    fn quat_to_matrix_is_a_proper_rotation() {
        let mut r = rng();
        for _ in 0..1000 {
            let a = random_rotation(&mut r, TAU);
            assert!((a.determinant() - 1.0).abs() < 1e-14);
            assert!((a.transpose() * a - Mat3::identity()).norm() < 1e-14);
        }
    }

    #[test]
    fn quat_matrix_round_trip() {
        let mut r = rng();
        for _ in 0..1000 {
            let axis = random_unit(&mut r);
            // An angle in [0, pi) keeps w > 0, so components compare directly.
            let angle = r.random::<f64>() * PI;
            let unit = UnitQuaternion::from_scaled_axis(axis * angle);
            let q = [unit.w, unit.i, unit.j, unit.k];
            let back = matrix_to_quat(&quat_to_matrix(q));
            for k in 0..4 {
                assert!((back[k] - q[k]).abs() < 1e-14, "{back:?} vs {q:?}");
            }
        }
    }

    #[test]
    fn matrix_to_quat_normalises_and_fixes_the_sign() {
        let mut r = rng();
        for _ in 0..1000 {
            // Angles past pi are where the raw conversion can return w < 0.
            let a = random_rotation(&mut r, TAU);
            let q = matrix_to_quat(&a);
            assert!(q[0] >= 0.0, "w must be non-negative, got {q:?}");
            let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
            assert!((norm - 1.0).abs() < 1e-14);
            // Flipping the sign must not have changed the rotation itself.
            assert!((quat_to_matrix(q) - a).norm() < 1e-14);
        }
    }

    #[test]
    fn quat_to_matrix_normalises_its_input() {
        let q = [0.5, 0.5, 0.5, 0.5];
        let scaled = [1.5, 1.5, 1.5, 1.5];
        assert!((quat_to_matrix(scaled) - quat_to_matrix(q)).norm() < 1e-15);
    }

    #[test]
    fn rotation_vector_of_identity_is_zero() {
        assert_eq!(rotation_vector(&Mat3::identity()), Vec3::zeros());
    }

    #[test]
    fn rotation_vector_recovers_axis_times_angle() {
        let mut r = rng();
        for _ in 0..1000 {
            let axis = random_unit(&mut r);
            // (0, pi) is the range the logarithm can represent.
            let angle = 1e-6 + r.random::<f64>() * (PI - 2e-6);
            let a = UnitQuaternion::from_scaled_axis(axis * angle)
                .to_rotation_matrix()
                .into_inner();
            let phi = rotation_vector(&a);
            assert!((phi.norm() - angle).abs() < 1e-13);
            assert!((phi - axis * angle).norm() < 1e-13);
        }
    }

    #[test]
    fn rotation_vector_is_accurate_at_arcsecond_scale() {
        let mut r = rng();
        for exponent in 0..7 {
            let angle = ARCSEC * 10f64.powi(-exponent);
            let axis = random_unit(&mut r);
            let a = UnitQuaternion::from_scaled_axis(axis * angle)
                .to_rotation_matrix()
                .into_inner();
            let phi = rotation_vector(&a);
            assert!(
                (phi.norm() - angle).abs() <= 1e-9 * angle,
                "angle={angle:e} measured={:e}",
                phi.norm()
            );
            assert!((phi - axis * angle).norm() <= 1e-9 * angle);
        }
    }

    /// The error metric of docs/SPEC.md end to end: a known perturbation of a
    /// true attitude must come back out of `A_est * A_true^T`.
    #[test]
    fn rotation_vector_measures_attitude_error() {
        let mut r = rng();
        for _ in 0..200 {
            let a_true = random_rotation(&mut r, TAU);
            let delta = Vec3::new(3.0, -4.0, 12.0) * ARCSEC; // 13 arcsec total
            let a_est = UnitQuaternion::from_scaled_axis(delta)
                .to_rotation_matrix()
                .into_inner()
                * a_true;
            let phi = rotation_vector(&(a_est * a_true.transpose()));
            assert!((phi - delta).norm() < 1e-15);
            assert!((phi.norm() / ARCSEC - 13.0).abs() < 1e-9);
        }
    }

    #[test]
    fn radec_round_trip() {
        let mut r = rng();
        for _ in 0..1000 {
            let v = random_unit(&mut r);
            let (ra, dec) = unit_to_radec(&v);
            assert!((0.0..TAU).contains(&ra));
            assert!(dec.abs() <= FRAC_PI_2);
            assert!((radec_to_unit(ra, dec) - v).norm() < 1e-15);

            // Magnitude must not matter.
            let (ra2, dec2) = unit_to_radec(&(v * 42.0));
            assert!((ra2 - ra).abs() < 1e-15 && (dec2 - dec).abs() < 1e-15);
        }
    }

    #[test]
    fn radec_known_values_and_wrapping() {
        for (ra, dec, expected) in [
            (0.0, 0.0, Vec3::new(1.0, 0.0, 0.0)),
            (FRAC_PI_2, 0.0, Vec3::new(0.0, 1.0, 0.0)),
            (0.0, FRAC_PI_2, Vec3::new(0.0, 0.0, 1.0)),
            (3.0 * FRAC_PI_2, 0.0, Vec3::new(0.0, -1.0, 0.0)),
        ] {
            assert!((radec_to_unit(ra, dec) - expected).norm() < 1e-15);
            let (back, _) = unit_to_radec(&expected);
            assert!((back - ra).abs() < 1e-15, "{back} vs {ra}");
        }

        // Just short of RA = 0 must wrap into [0, TAU), not land on TAU.
        let (ra, _) = unit_to_radec(&Vec3::new(1.0, -1e-18, 0.0));
        assert!((0.0..TAU).contains(&ra), "ra={ra}");
    }

    /// Pinned against the published SplitMix64 stream, whose first output for
    /// state 0 is 0xE220A8397B1DCDAF. Native and WASM must both reproduce
    /// these, so treat any change here as a determinism regression.
    #[test]
    fn splitmix64_matches_the_reference_stream() {
        assert_eq!(splitmix64(0), 0xE220_A839_7B1D_CDAF);
        assert_eq!(splitmix64(0x9E37_79B9_7F4A_7C15), 0x6E78_9E6A_A1B9_65F4);
        assert_eq!(splitmix64(0x3C6E_F372_FE94_F82A), 0x06C4_5D18_8009_454F);
    }

    #[test]
    fn splitmix64_is_deterministic_and_injective_over_a_sample() {
        let mut seen = std::collections::HashSet::new();
        for i in 0..10_000u64 {
            let h = splitmix64(i);
            assert_eq!(h, splitmix64(i), "not deterministic at {i}");
            assert!(seen.insert(h), "collision at {i}");
        }
        // Adjacent seeds must not stay adjacent: that is what makes
        // splitmix64(base_seed ^ i) usable as a per-trial seed.
        for i in 0..1000u64 {
            let gap = splitmix64(i).wrapping_sub(splitmix64(i + 1));
            assert!(gap > (1u64 << 20), "seeds {i} and {} too close", i + 1);
        }
    }
}
