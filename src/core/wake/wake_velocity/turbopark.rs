//! TurbOPark Wake Velocity Model
//!
//! Based on the TurbOPark model from Ørsted.
//! Uses pre-computed overlap lookup table

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array1, Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

const OVERLAP_DIST_N: usize = 100;
const OVERLAP_RADIUS_N: usize = 100;

#[derive(Debug, Clone)]
pub struct TurbOParkVelocityDeficit {
    pub base: BaseModel,
    pub a: Float,
    pub sigma_max_rel: Float,
    overlap_gauss_interp: OverlapInterpolator,
}

impl TurbOParkVelocityDeficit {
    pub fn new(a: Float, sigma_max_rel: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("a".to_string(), a);
        params.insert("sigma_max_rel".to_string(), sigma_max_rel);

        Self {
            base: BaseModel::new(params, "turbopark"),
            a,
            sigma_max_rel,
            overlap_gauss_interp: OverlapInterpolator::new(),
        }
    }

    pub fn default() -> Self {
        Self::new(0.04, 4.0)
    }
}

impl Default for TurbOParkVelocityDeficit {
    fn default() -> Self {
        Self::default()
    }
}

#[derive(Debug, Clone)]
pub struct OverlapInterpolator {
    dist: Vec<Float>,
    radius_down: Vec<Float>,
}

impl OverlapInterpolator {
    fn new() -> Self {
        Self {
            dist: super::turbopark_lookup_matlab::dist().to_vec(),
            radius_down: super::turbopark_lookup_matlab::radius_down().to_vec(),
        }
    }

    fn interpolate(&self, normalized_r: Float, normalized_radius_down: Float) -> Float {
        if normalized_r < 0.0 || normalized_radius_down <= 0.0 {
            return 0.0;
        }

        let i = Self::find_index(&self.dist, normalized_r);
        let j = Self::find_index(&self.radius_down, normalized_radius_down);

        let (i0, i1, f_r) = if i + 1 < self.dist.len() {
            (
                i,
                i + 1,
                (normalized_r - self.dist[i]) / (self.dist[i + 1] - self.dist[i]),
            )
        } else {
            (i.saturating_sub(1), i, 0.0)
        };

        let (j0, j1, f_rd) = if j + 1 < self.radius_down.len() {
            (
                j,
                j + 1,
                (normalized_radius_down - self.radius_down[j])
                    / (self.radius_down[j + 1] - self.radius_down[j]),
            )
        } else {
            (j.saturating_sub(1), j, 0.0)
        };

        let i0 = i0.min(super::turbopark_lookup_matlab::DIST_N - 1);
        let i1 = i1.min(super::turbopark_lookup_matlab::DIST_N - 1);
        let j0 = j0.min(super::turbopark_lookup_matlab::RADIUS_DOWN_N - 1);
        let j1 = j1.min(super::turbopark_lookup_matlab::RADIUS_DOWN_N - 1);

        let v00 = super::turbopark_lookup_matlab::overlap_at(i0, j0);
        let v01 = super::turbopark_lookup_matlab::overlap_at(i0, j1);
        let v10 = super::turbopark_lookup_matlab::overlap_at(i1, j0);
        let v11 = super::turbopark_lookup_matlab::overlap_at(i1, j1);

        let v0 = v00 * (1.0 - f_rd) + v01 * f_rd;
        let v1 = v10 * (1.0 - f_rd) + v11 * f_rd;

        v0 * (1.0 - f_r) + v1 * f_r
    }

    #[inline]
    fn find_index(arr: &[Float], value: Float) -> usize {
        arr.iter()
            .enumerate()
            .take(arr.len() - 1)
            .take_while(|(_, &v)| value >= v)
            .count()
            .min(arr.len() - 1)
    }
}

#[allow(dead_code)]
fn create_overlap_lookup_table() -> (Vec<Float>, Vec<Float>, Array2) {
    let max_dist = 10.0;
    let max_radius_down = 20.0;

    let dist: Vec<Float> = (0..OVERLAP_DIST_N)
        .map(|i| max_dist * i as Float / (OVERLAP_DIST_N - 1) as Float)
        .collect();

    let radius_down: Vec<Float> = (0..OVERLAP_RADIUS_N)
        .map(|i| max_radius_down * i as Float / (OVERLAP_RADIUS_N - 1) as Float)
        .collect();

    let mut overlap_gauss = Array::zeros((OVERLAP_DIST_N, OVERLAP_RADIUS_N));

    for (i, &d) in dist.iter().enumerate() {
        for (j, &rd) in radius_down.iter().enumerate() {
            overlap_gauss[[i, j]] = calculate_overlap_integral(d, rd);
        }
    }

    (dist, radius_down, overlap_gauss)
}

#[allow(dead_code)]
fn calculate_overlap_integral(dist: Float, radius_down: Float) -> Float {
    if radius_down <= 0.0 {
        return (-(dist.powi(2)) / 2.0).exp();
    }

    if dist > 10.0 {
        return (-(dist.powi(2)) / 2.0).exp();
    }

    let n_points = 100;
    let dr = radius_down / (n_points as Float - 1.0);
    let mut integral = 0.0;

    for k in 0..n_points {
        let r = k as Float * dr;
        // 2D polar integral: the theta integral of exp(dist*r*cos(theta))
        // evaluates to 2*pi*I_0(dist*r) (modified Bessel function)
        let decay = (-((r.powi(2) + dist.powi(2)) / 2.0)).exp()
            * 2.0
            * std::f64::consts::PI as Float
            * bessel_i0(dist * r);
        let integrand = r * decay;
        let weight = if k == 0 || k == n_points - 1 {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * integrand;
    }

    integral *= dr / 3.0;
    let area = std::f64::consts::PI as Float * radius_down.powi(2);
    (integral / area).max(0.0).min(1.0)
}

#[inline]
fn bessel_i0(x: Float) -> Float {
    let mut sum: Float = 1.0;
    let mut term: Float = 1.0;
    let mut k = 1;
    while term.abs() > 1e-12 * sum.abs() && k < 50 {
        term *= (x * x / 4.0) / (k as f64 * k as f64);
        sum += term;
        k += 1;
    }
    sum
}

#[inline]
fn calculate_epsilon(ct: Float) -> Float {
    // Python (turbopark.py:106-108): 0.25 * sqrt(0.5*(1+sqrt(1-Ct))/sqrt(1-Ct))
    // with NO scalar cap — the np.min(..., 3, keepdims=True) reduces over a
    // length-1 axis (the solver passes Cts as (n_findex, n_turbines, 1, 1)),
    // so it is a no-op (Bug #32). Ct >= 1 degenerates to epsilon -> inf, which
    // then yields C -> 0 downstream, matching Python.
    let sqrt_one_minus_ct = (1.0 - ct).max(0.0).sqrt();
    if sqrt_one_minus_ct < 1e-10 {
        return Float::INFINITY;
    }
    let ratio = 0.5 * (1.0 + sqrt_one_minus_ct) / sqrt_one_minus_ct;
    0.25 * ratio.sqrt()
}

fn characteristic_wake_width(x_dist: Float, ti: Float, ct: Float, a: Float) -> Float {
    let c1 = 1.5;
    let c2 = 0.8;

    let alpha = ti * c1;
    let beta = (c2 * ti / ct.sqrt()).max(1e-10);

    let term1 = ((alpha + beta * x_dist).powi(2) + 1.0).sqrt();
    let term2 = (1.0 + alpha.powi(2)).sqrt();

    let dw = a * ti / beta
        * (term1
            - term2
            - ((term1 + 1.0) * alpha / ((term2 + 1.0) * (alpha + beta * x_dist))).ln());

    dw.max(0.0)
}

#[inline]
fn calculate_peak_deficit(ct: Float, sigma_over_d: Float) -> Float {
    if sigma_over_d <= 0.0 {
        return 0.0;
    }
    let val = 1.0 - ct / (8.0 * sigma_over_d.powi(2));
    if val <= 0.0 {
        return 1.0;
    }
    1.0 - val.sqrt()
}

/// Target-centric deficit (port of Python turbopark_solver's model call).
///
/// Computes the deficit experienced by turbine `target_index` from all
/// upstream turbines (including mirrored/image wakes), each evaluated at the
/// source's grid points, combined via sqrt-sum-of-squares.
#[allow(clippy::too_many_arguments)]
pub fn turbopark_velocity_deficit(
    turbulence_intensity: &Array4,
    thrust_coefficients: &Array2,
    rotor_diameters: &Array1,
    grid_x: &Array4,
    grid_y: &Array4,
    grid_z: &Array4,
    deflection_field: &Array4,
    target_index: usize,
    a: Float,
    sigma_max_rel: Float,
) -> Array4 {
    let shape = grid_x.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];
    let n_y = shape[2];
    let n_z = shape[3];

    let mut sum_sq = Array4::zeros((n_findex, n_turbines, n_y, n_z));
    let interp = OverlapInterpolator::new();

    // Python base.py NUM_EPS = 0.001 (Bug #32)
    const NUM_EPS: Float = 0.001;

    for fi in 0..n_findex {
        // Target turbine's mean positions (x is constant per turbine)
        let x_i = grid_x[[fi, target_index, 0, 0]];
        let mut y_i = 0.0;
        let mut z_i = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                y_i += grid_y[[fi, target_index, iy, iz]];
                z_i += grid_z[[fi, target_index, iy, iz]];
            }
        }
        let n_pts = (n_y * n_z) as f64;
        y_i /= n_pts;
        z_i /= n_pts;

        let d_i = rotor_diameters[target_index];

        for m in 0..n_turbines {
            let d_m = rotor_diameters[m];
            let ct_m = thrust_coefficients[[fi, m]];

            for iy in 0..n_y {
                for iz in 0..n_z {
                    let x_m = grid_x[[fi, m, iy, iz]];
                    let y_m = grid_y[[fi, m, iy, iz]];
                    let z_m = grid_z[[fi, m, iy, iz]];

                    // Normalized distance from the target to the source's grid
                    // point; positive only for sources upstream of the target
                    let x_dist = (x_i - x_m) / d_m;
                    if x_dist <= NUM_EPS {
                        continue;
                    }

                    // Per-point TI of source m (Python passes the evolved
                    // turbine_turbulence_intensity field (n_findex, n_turbines,
                    // gy, gz); characteristic_wake_width broadcasts pointwise)
                    let ti_m = turbulence_intensity[[fi, m, iy, iz]];
                    let dw = characteristic_wake_width(x_dist, ti_m, ct_m, a);
                    let eps = calculate_epsilon(ct_m);
                    let sigma = d_m * (eps + dw);

                    let defl = deflection_field[[fi, m, iy, iz]];
                    let dy = y_i - (y_m + defl);
                    let dz = z_i - z_m;
                    let r_dist = (dy * dy + dz * dz).sqrt();

                    let c_m = calculate_peak_deficit(ct_m, sigma / d_m);

                    let is_overlapping = sigma_max_rel * sigma / 2.0 + d_i / 2.0 > r_dist;
                    if !is_overlapping {
                        continue;
                    }

                    let overlap_real = interp.interpolate(r_dist / sigma, d_i / 2.0 / sigma);
                    let dz_image = z_i + z_m;
                    let r_dist_image = (dy * dy + dz_image * dz_image).sqrt();
                    let overlap_image = interp.interpolate(r_dist_image / sigma, d_i / 2.0 / sigma);

                    let d_real = c_m * overlap_real;
                    let d_image = c_m * overlap_image;

                    let cur = sum_sq[[fi, target_index, iy, iz]];
                    sum_sq[[fi, target_index, iy, iz]] = cur + d_real * d_real + d_image * d_image;
                }
            }
        }
    }

    sum_sq.mapv(|v| v.sqrt())
}

impl VelocityModel for TurbOParkVelocityDeficit {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn clone_box(&self) -> Box<dyn VelocityModel> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn function(
        &self,
        x: Array4,
        y: Array4,
        z: Array4,
        _axial_induction: Array2,
        deflection_field: Array4,
        _yaw_angle: Array2,
        turbulence_intensity: Array4,
        thrust_coefficient: Array2,
        hub_height: Float,
        rotor_diameter: Float,
        turbine_index: usize,
        _model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let n_y = shape[2];
        let n_z = shape[3];

        let r0 = rotor_diameter / 2.0;
        let mut velocity_deficit = Array::zeros((n_findex, n_turbines, n_y, n_z));

        for fi in 0..n_findex {
            let mut x_wake_source = 0.0;
            let mut y_wake_source = 0.0;
            let mut count = 0.0;
            for iy in 0..n_y {
                for iz in 0..n_z {
                    x_wake_source += x[[fi, turbine_index, iy, iz]];
                    y_wake_source += y[[fi, turbine_index, iy, iz]];
                    count += 1.0;
                }
            }
            x_wake_source /= count;
            y_wake_source /= count;

            // Turbine upstream of reference point (x < 0) doesn't generate a wake
            if x_wake_source < 0.0 {
                continue;
            }

            let epsilon = calculate_epsilon(thrust_coefficient[[fi, turbine_index]]);

            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        if x_point <= x_wake_source {
                            continue;
                        }

                        // Per-point TI (source's slice of the evolving 4D field)
                        let ti_pt = turbulence_intensity[[fi, turbine_index, iy, iz]];
                        let c = calculate_peak_deficit(
                            thrust_coefficient[[fi, turbine_index]],
                            self.a * ti_pt + epsilon,
                        );

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        // Wake center with deflection evaluated at this point
                        // (Python: y_i - (y + deflection_field))
                        let wake_center_y = y_wake_source + deflection_field[[fi, ti, iy, iz]];
                        let dy = y_point - wake_center_y;
                        let dz = z_point - hub_height;
                        let r = (dy.powi(2) + dz.powi(2)).sqrt();

                        let x_dist = (x_point - x_wake_source) / rotor_diameter;
                        let dw = characteristic_wake_width(
                            x_dist,
                            ti_pt,
                            thrust_coefficient[[fi, turbine_index]],
                            self.a,
                        );
                        let sigma_downstream = (self.a * ti_pt + epsilon + dw).max(epsilon);

                        let normalized_r = r / (sigma_downstream * rotor_diameter);
                        let normalized_radius_down = r0 / (sigma_downstream * rotor_diameter);
                        let overlap = self
                            .overlap_gauss_interp
                            .interpolate(normalized_r, normalized_radius_down);

                        velocity_deficit[[fi, ti, iy, iz]] = (c * overlap).max(0.0);
                    }
                }
            }
        }

        Ok(velocity_deficit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;

    #[test]
    fn test_turbopark_creation() {
        let park = TurbOParkVelocityDeficit::new(0.04, 4.0);
        assert_eq!(park.a, 0.04);
        assert_eq!(park.sigma_max_rel, 4.0);
    }

    #[test]
    fn test_turbopark_default() {
        let park = TurbOParkVelocityDeficit::default();
        assert_eq!(park.a, 0.04);
        assert_eq!(park.sigma_max_rel, 4.0);
    }

    #[test]
    fn test_turbopark_prepare_function() {
        let park = TurbOParkVelocityDeficit::default();
        let result = park.prepare_function(&fake_grid(), &fake_flow_field());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn test_characteristic_wake_width() {
        let width = characteristic_wake_width(5.0, 0.06, 0.8, 0.04);
        assert!(width > 0.0);
    }

    #[test]
    fn test_epsilon() {
        let eps = calculate_epsilon(0.8);
        assert!(eps > 0.0 && eps < 1.0);
    }

    #[test]
    fn test_peak_deficit() {
        let c = calculate_peak_deficit(0.8, 0.5);
        assert!(c > 0.0 && c < 1.0);
    }

    #[test]
    fn test_overlap_interpolation() {
        let interp = OverlapInterpolator::new();
        let center = interp.interpolate(0.0, 1.0);
        let far = interp.interpolate(5.0, 1.0);
        // Center should have higher overlap than edge
        assert!(center > 0.0);
        // Far should have lower overlap than near
        assert!(far < center);
    }

    #[test]
    fn test_calculate_overlap_integral() {
        // At zero distance with rotor radius
        let overlap_zero = calculate_overlap_integral(0.0, 1.0);
        assert!(overlap_zero > 0.0);
        // At large distance, overlap should be small
        let overlap_large = calculate_overlap_integral(10.0, 1.0);
        assert!(overlap_large < overlap_zero);
    }

    fn fake_grid() -> impl crate::core::Grid {
        struct FakeGrid;
        impl crate::core::Grid for FakeGrid {
            fn n_turbines(&self) -> usize {
                1
            }
            fn n_findex(&self) -> usize {
                1
            }
            fn x_sorted(&self) -> &Array4 {
                panic!()
            }
            fn y_sorted(&self) -> &Array4 {
                panic!()
            }
            fn z_sorted(&self) -> &Array4 {
                panic!()
            }
            fn x_sorted_inertial_frame(&self) -> &Array4 {
                panic!()
            }
            fn y_sorted_inertial_frame(&self) -> &Array4 {
                panic!()
            }
            fn z_sorted_inertial_frame(&self) -> &Array4 {
                panic!()
            }
            fn cubature_weights(&self) -> Option<&Array2> {
                None
            }
            fn average_method(&self) -> crate::core::AveragingMethod {
                crate::core::AveragingMethod::CubicMean
            }
            fn sorted_indices(&self) -> &Array2 {
                static INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
                INDICES.get_or_init(|| Array2::zeros((1, 1)))
            }
            fn sorted_coord_indices(&self) -> &Array2 {
                static COORD_INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
                COORD_INDICES.get_or_init(|| Array2::zeros((1, 1)))
            }
            fn resolution(&self) -> usize {
                1
            }
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }
        }
        FakeGrid
    }

    fn fake_flow_field() -> crate::core::FlowField {
        let empty_4d = ndarray::Array4::zeros((0, 0, 0, 0));
        crate::core::FlowField {
            wind_speeds: Array1::from_vec(vec![8.0]),
            wind_directions: Array1::from_vec(vec![270.0]),
            turbulence_intensities: Array1::from_vec(vec![0.06]),
            air_density: 1.225,
            wind_shear: 0.12,
            wind_veer: 0.0,
            reference_wind_height: 90.0,
            n_findex: 1,
            u_initial_sorted: empty_4d.clone(),
            v_initial_sorted: empty_4d.clone(),
            w_initial_sorted: empty_4d.clone(),
            u_sorted: empty_4d.clone(),
            v_sorted: empty_4d.clone(),
            w_sorted: empty_4d.clone(),
            u: empty_4d.clone(),
            v: empty_4d.clone(),
            w: empty_4d.clone(),
            turbulence_intensity_field: empty_4d.clone(),
            dudz_initial_sorted: empty_4d.clone(),
            turbulence_intensity_field_sorted: empty_4d,
        }
    }
}
