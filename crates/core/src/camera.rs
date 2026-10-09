//! Pinhole camera model with optional radial distortion.
//!
//! Camera frame (OpenCV convention, right-handed): `+z` is the boresight, `+x`
//! points along increasing column (right in the image) and `+y` along
//! increasing row (down the image). Pixel centres sit at integer `(col, row)`
//! and the principal point defaults to `((W-1)/2, (H-1)/2)`.
//!
//! Attitude maps inertial to camera: `b = A r`.

use crate::math::Vec3;

/// A pinhole camera: sensor size, focal length in pixels, principal point and
/// one radial distortion coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CameraModel {
    /// Sensor width in pixels.
    pub width: usize,
    /// Sensor height in pixels.
    pub height: usize,
    /// Focal length in pixels.
    pub f: f64,
    /// Principal point column.
    pub cx: f64,
    /// Principal point row.
    pub cy: f64,
    /// Second-order radial distortion coefficient, applied to *normalised*
    /// image coordinates. Zero for an ideal pinhole.
    pub k1: f64,
}

impl CameraModel {
    /// Builds an ideal camera (`k1 = 0`) from its full-width field of view.
    ///
    /// `fov` is in radians across the full sensor width, so
    /// `f = (W/2) / tan(FOV/2)`; the default 1024 px at 20 deg gives
    /// f = 2904.1 px, as docs/SPEC.md states.
    pub fn from_fov(width: usize, height: usize, fov: f64) -> Self {
        let f = (width as f64 / 2.0) / (fov / 2.0).tan();
        Self {
            width,
            height,
            f,
            cx: (width as f64 - 1.0) / 2.0,
            cy: (height as f64 - 1.0) / 2.0,
            k1: 0.0,
        }
    }

    /// Projects a camera-frame direction to pixel coordinates.
    ///
    /// Applies `k1`, so this is the *physical* forward model: the simulator
    /// uses it with the true focal length and distortion. Returns `None` when
    /// `b.z <= 0`, where the projection is undefined (the direction is at or
    /// behind the focal plane).
    ///
    /// `x = cx + f*bx/bz`, `y = cy + f*by/bz`, with the normalised radius
    /// scaled by `1 + k1*r^2` before the focal length is applied.
    pub fn project(&self, b: &Vec3) -> Option<[f64; 2]> {
        if b.z <= 0.0 {
            return None;
        }
        let xn = b.x / b.z;
        let yn = b.y / b.z;
        let scale = 1.0 + self.k1 * (xn * xn + yn * yn);
        Some([self.cx + self.f * xn * scale, self.cy + self.f * yn * scale])
    }

    /// The pinhole inverse: `b = normalize(x - cx, y - cy, f)`.
    ///
    /// Deliberately ignores `k1`. This is the model the *solver* uses, and the
    /// simulator's distortion is an error the solver is never told about, so
    /// this is the exact inverse of [`CameraModel::project`] only when
    /// `k1 == 0`. See `no_pair_is_closer...`-style tests below for the
    /// measured residual when it is not.
    pub fn unproject(&self, x: f64, y: f64) -> Vec3 {
        Vec3::new(x - self.cx, y - self.cy, self.f).normalize()
    }

    /// Half the diagonal field of view, radians.
    ///
    /// The largest angle from the boresight at which a star can still land on
    /// the sensor, which is what bounds the catalogue search and the pair
    /// database's separation range.
    pub fn half_diagonal_fov(&self) -> f64 {
        let corner = self.cx.hypot(self.cy);
        (corner / self.f).atan()
    }

    /// True when `(x, y)` falls within the sensor's physical area.
    ///
    /// Pixel centres are at integers, so pixel `0` spans `[-0.5, 0.5)` and the
    /// sensor spans `[-0.5, W-0.5)`.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= -0.5 && x < self.width as f64 - 0.5 && y >= -0.5 && y < self.height as f64 - 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::angle_between;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use std::f64::consts::PI;

    fn nominal() -> CameraModel {
        CameraModel::from_fov(1024, 1024, 20.0 * PI / 180.0)
    }

    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(0x0CAD_EA5E)
    }

    #[test]
    fn focal_length_matches_the_documented_default() {
        let cam = nominal();
        // 512/tan(10 deg); docs/SPEC.md quotes this as "f ~ 2904 px".
        assert!((cam.f - 2_903.696_292).abs() < 1e-5, "f = {}", cam.f);
        assert_eq!((cam.cx, cam.cy), (511.5, 511.5));
    }

    #[test]
    fn boresight_lands_on_the_principal_point() {
        let cam = nominal();
        let p = cam.project(&Vec3::new(0.0, 0.0, 1.0)).expect("in front");
        assert!((p[0] - cam.cx).abs() < 1e-12 && (p[1] - cam.cy).abs() < 1e-12);
    }

    /// `+x` must be right (increasing column) and `+y` down (increasing row).
    #[test]
    fn axes_point_the_documented_way() {
        let cam = nominal();
        let right = cam.project(&Vec3::new(0.1, 0.0, 1.0)).expect("in front");
        let down = cam.project(&Vec3::new(0.0, 0.1, 1.0)).expect("in front");
        assert!(right[0] > cam.cx && (right[1] - cam.cy).abs() < 1e-12);
        assert!(down[1] > cam.cy && (down[0] - cam.cx).abs() < 1e-12);
    }

    #[test]
    fn rejects_directions_at_or_behind_the_focal_plane() {
        let cam = nominal();
        assert!(cam.project(&Vec3::new(0.0, 0.0, -1.0)).is_none());
        assert!(cam.project(&Vec3::new(1.0, 0.0, 0.0)).is_none());
    }

    /// Phase 2 acceptance: the projection round-trip error is below 1e-9 px.
    #[test]
    fn pixel_round_trip_is_exact_for_an_ideal_pinhole() {
        let cam = nominal();
        let mut r = rng();
        let mut worst = 0.0f64;
        for _ in 0..10_000 {
            let x = r.random::<f64>() * (cam.width as f64 - 1.0);
            let y = r.random::<f64>() * (cam.height as f64 - 1.0);
            let back = cam
                .project(&cam.unproject(x, y))
                .expect("unprojected rays point forward");
            worst = worst.max((back[0] - x).abs().max((back[1] - y).abs()));
        }
        assert!(worst < 1e-9, "worst pixel round-trip error {worst:e} px");
    }

    /// The other direction: direction -> pixel -> direction must return the
    /// same ray.
    #[test]
    fn direction_round_trip_is_exact_for_an_ideal_pinhole() {
        let cam = nominal();
        let mut r = rng();
        let mut worst = 0.0f64;
        for _ in 0..10_000 {
            // Any direction inside the field, pointing forward.
            let xn = (r.random::<f64>() - 0.5) * 2.0 * 0.18;
            let yn = (r.random::<f64>() - 0.5) * 2.0 * 0.18;
            let b = Vec3::new(xn, yn, 1.0).normalize();
            let p = cam.project(&b).expect("in front");
            worst = worst.max(angle_between(&cam.unproject(p[0], p[1]), &b));
        }
        assert!(worst < 1e-15, "worst direction round-trip {worst:e} rad");
    }

    /// `unproject` is the solver's model and ignores `k1`, so with distortion
    /// present the round-trip deliberately does not close. This pins how big
    /// that unmodelled error is, since it is the error budget the identifier
    /// has to absorb on the `hard` and `brutal` presets.
    #[test]
    fn distortion_is_the_unmodelled_error_it_is_meant_to_be() {
        let mut cam = nominal();
        let corner = cam.unproject(0.0, 0.0);

        for (k1, lower, upper) in [(0.01, 0.3, 0.6), (0.03, 1.0, 1.8)] {
            cam.k1 = k1;
            let p = cam.project(&corner).expect("in front");
            let shift = ((p[0] - 0.0).powi(2) + (p[1] - 0.0).powi(2)).sqrt();
            assert!(
                (lower..upper).contains(&shift),
                "k1 = {k1} shifts the corner {shift:.3} px, expected {lower}..{upper}"
            );
        }
    }

    #[test]
    fn half_diagonal_fov_is_consistent_with_the_corner() {
        let cam = nominal();
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let corner = cam.unproject(0.0, 0.0);
        let measured = angle_between(&axis, &corner);
        assert!((measured - cam.half_diagonal_fov()).abs() < 1e-12);
        // 20 deg across the width implies a little under 28 deg diagonally.
        let full_diagonal = 2.0 * cam.half_diagonal_fov().to_degrees();
        assert!(
            (27.0..29.0).contains(&full_diagonal),
            "diagonal FOV {full_diagonal} deg"
        );
    }

    #[test]
    fn contains_covers_the_physical_sensor_area() {
        let cam = nominal();
        assert!(cam.contains(0.0, 0.0));
        assert!(cam.contains(-0.5, -0.5));
        assert!(cam.contains(1023.0, 1023.0));
        assert!(!cam.contains(-0.51, 0.0));
        assert!(!cam.contains(1023.5, 0.0));
        assert!(!cam.contains(0.0, 1023.5));
    }

    #[test]
    fn non_square_sensors_keep_their_principal_point_centred() {
        let cam = CameraModel::from_fov(640, 480, 0.3);
        assert_eq!((cam.cx, cam.cy), (319.5, 239.5));
        let p = cam.project(&Vec3::new(0.0, 0.0, 1.0)).expect("in front");
        assert!((p[0] - 319.5).abs() < 1e-12 && (p[1] - 239.5).abs() < 1e-12);
    }
}
