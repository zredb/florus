//! Cumulative Gauss Curl Wake Velocity Model
//!
//! Implements the cumulative-curl wake model for capturing deep array effects
//! and impacts to wake steering.
//!
//! Based on Bay et al. (2023) "Addressing deep array effects and impacts to wake
//! steering with the cumulative-curl wake model" - Wind Energy Science
//!
//! This model extends the Gaussian wake model with:
//! - Wake expansion continuation (WEC) for cumulative wake effects
//! - Curl terms for yawed turbine wake deformation

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array1, Array2, Array4, ArrayD, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Cumulative Gauss Curl wake velocity model
#[derive(Debug, Clone)]
pub struct CumulativeCurlVelocityDeficit {
    pub base: BaseModel,
    // WEC (Wake Expansion Continuation) parameters
    pub wec_factor: Float,
    // Near-field parameters
    pub a_s: Float,
    pub b_s: Float,
    pub c_s1: Float,
    pub c_s2: Float,
    // Far-field parameters
    pub a_f: Float,
    pub b_f: Float,
    pub c_f: Float,
}

impl CumulativeCurlVelocityDeficit {
    pub fn new(
        wec_factor: Float,
        a_s: Float,
        b_s: Float,
        c_s1: Float,
        c_s2: Float,
        a_f: Float,
        b_f: Float,
        c_f: Float,
    ) -> Self {
        let mut params = HashMap::new();
        params.insert("wec_factor".to_string(), wec_factor);
        params.insert("a_s".to_string(), a_s);
        params.insert("b_s".to_string(), b_s);
        params.insert("c_s1".to_string(), c_s1);
        params.insert("c_s2".to_string(), c_s2);
        params.insert("a_f".to_string(), a_f);
        params.insert("b_f".to_string(), b_f);
        params.insert("c_f".to_string(), c_f);

        Self {
            base: BaseModel::new(params, "cumulative_gauss_curl"),
            wec_factor,
            a_s,
            b_s,
            c_s1,
            c_s2,
            a_f,
            b_f,
            c_f,
        }
    }

    pub fn default() -> Self {
        Self::new(
            1.0,          // wec_factor
            0.179367259,  // a_s
            0.0118889215, // b_s
            0.0563691592, // c_s1
            0.13290157,   // c_s2
            3.11,         // a_f
            -0.68,        // b_f
            2.41,         // c_f
        )
    }
}

impl Default for CumulativeCurlVelocityDeficit {
    fn default() -> Self {
        Self::default()
    }
}

/// Calculate wake expansion (port of Python cumulative_gauss_curl.wake_expansion).
/// Returns sigma in rotor-diameter units.
fn calculate_wake_spread(
    x: Float,
    ti: Float,
    ct: Float,
    rotor_diameter: Float,
    params: &CumulativeCurlVelocityDeficit,
) -> Float {
    let _ = rotor_diameter;

    let beta = 0.5 * (1.0 + (1.0 - ct).sqrt()) / (1.0 - ct).sqrt();
    let k = params.a_s * ti + params.b_s;
    let eps = (params.c_s1 * ct + params.c_s2) * beta.sqrt();

    let x_tilde = x.abs() / rotor_diameter;
    k * x_tilde + eps
}

/// Calculate Gaussian wake deficit
fn gaussian_deficit(r: Float, sigma: Float, ct: Float) -> Float {
    if sigma <= 0.0 {
        return 0.0;
    }

    // Gaussian profile
    let exponent = -r.powi(2) / (2.0 * sigma.powi(2));
    let profile = exponent.exp();

    // Maximum deficit based on actuator disk theory
    let deficit_max = 0.5 * (1.0 - (1.0 - ct).sqrt());

    deficit_max * profile
}

/// Cumulative curl deficit (port of Python cumulative_gauss_curl.function).
///
/// Blondel super-Gaussian deficit with cumulative effects: the wake evolution
/// from previous turbines (m < source_index - 1) accumulates into sum_lbda over
/// the full grid and corrects the peak deficit. The deficit is accumulated into
/// turb_u_wake (velocity units) scaled by the source turbine's rotor-averaged
/// inflow velocity. Returns the updated (turb_u_wake, Ctmp).
#[allow(clippy::too_many_arguments)]
pub fn cc_velocity_deficit(
    params: &CumulativeCurlVelocityDeficit,
    source_index: usize,
    grid_x: &Array4,
    grid_y: &Array4,
    grid_z: &Array4,
    u_i: &Array4,
    u_initial: &Array4,
    deflection_field: &Array4,
    yaw_angle: &Array1,
    turbulence_intensity: &Array4,
    turb_cts: &Array2,
    rotor_diameters: &Array2,
    mut turb_u_wake: Array4,
    mut ctmp: ArrayD,
) -> (Array4, ArrayD) {
    let shape = grid_x.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];
    let n_y = shape[2];
    let n_z = shape[3];
    let ii = source_index;

    let wake_expansion = |delta_x_m: Float, ct_m: Float, ti_m: Float, d_m: Float| -> Float {
        let sqrt_omc = (1.0 - ct_m).max(1e-10).sqrt();
        let beta = 0.5 * (1.0 + sqrt_omc) / sqrt_omc;
        let k = params.a_s * ti_m + params.b_s;
        let eps = (params.c_s1 * ct_m + params.c_s2) * beta.sqrt();
        k * (delta_x_m.abs() / d_m) + eps
    };

    // Per-(findex, turbine) mean coordinates
    let mut x_mean = Array2::zeros((n_findex, n_turbines));
    let mut y_mean = Array2::zeros((n_findex, n_turbines));
    let mut z_mean = Array2::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for m in 0..n_turbines {
            let mut ys = 0.0;
            let mut zs = 0.0;
            for iy in 0..n_y {
                for iz in 0..n_z {
                    ys += grid_y[[fi, m, iy, iz]];
                    zs += grid_z[[fi, m, iy, iz]];
                }
            }
            let n_pts = (n_y * n_z) as f64;
            x_mean[[fi, m]] = grid_x[[fi, m, 0, 0]];
            y_mean[[fi, m]] = ys / n_pts;
            z_mean[[fi, m]] = zs / n_pts;
        }
    }

    // The source turbine's rotor-averaged inflow velocity (cbrt mean)
    let mut src_avg_vel = Array1::zeros(n_findex);
    for fi in 0..n_findex {
        let mut v_cubed = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                v_cubed += u_i[[fi, 0, iy, iz]].powi(3);
            }
        }
        src_avg_vel[fi] = (v_cubed / (n_y * n_z) as f64).powf(1.0 / 3.0);
    }

    for fi in 0..n_findex {
        let ct_i = turb_cts[[fi, ii]];
        let yaw_i = yaw_angle[fi];
        let d_i = rotor_diameters[[fi, ii]];
        let x_i_loc = x_mean[[fi, ii]];
        let y_i_loc = y_mean[[fi, ii]];
        let z_i_loc = z_mean[[fi, ii]];

        // Cumulative loop: previous turbines m < ii - 1, over the full grid
        let mut sum_lbda = Array4::zeros((n_findex, n_turbines, n_y, n_z));
        for m in 0..ii.saturating_sub(1) {
            let ct_m = turb_cts[[fi, m]];
            let d_m = rotor_diameters[[fi, m]];
            let y_m_loc = y_mean[[fi, m]];
            let z_m_loc = z_mean[[fi, m]];

            for t in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_t = grid_x[[fi, t, iy, iz]];
                        // Per-point TI: source slices of the evolving 4D field
                        // (Python: turbine_ti[:, ii:ii+1] / [:, m:m+1],
                        // (n_findex, 1, gy, gz), broadcast over targets)
                        let ti_i = turbulence_intensity[[fi, ii, iy, iz]];
                        let ti_m = turbulence_intensity[[fi, m, iy, iz]];
                        let sigma_n = wake_expansion(x_t - x_i_loc, ct_i, ti_i, d_i);
                        let sigma_i = wake_expansion(x_t - x_mean[[fi, m]], ct_m, ti_m, d_m);

                        let s_i = sigma_n * sigma_n + sigma_i * sigma_i;
                        let defl = deflection_field[[fi, t, iy, iz]];
                        let y_diff = y_i_loc - y_m_loc - defl;
                        let z_diff = z_i_loc - z_m_loc;
                        let lbda = params.wec_factor * sigma_i * sigma_i / s_i
                            * (-(y_diff * y_diff) / (2.0 * s_i)).exp()
                            * (-(z_diff * z_diff) / (2.0 * s_i)).exp();

                        let ctmp_m = ctmp[[m, fi, t, iy, iz]];
                        let u_init = u_initial[[fi, t, iy, iz]].max(1e-10);
                        sum_lbda[[fi, t, iy, iz]] += lbda * (ctmp_m / u_init);
                    }
                }
            }
        }

        // Blondel super-Gaussian deficit with cumulative correction, full grid
        for t in 0..n_turbines {
            for iy in 0..n_y {
                for iz in 0..n_z {
                    let x_pt = grid_x[[fi, t, iy, iz]];
                    let y_pt = grid_y[[fi, t, iy, iz]];
                    let z_pt = grid_z[[fi, t, iy, iz]];

                    let delta_x = x_pt - x_i_loc;
                    let x_tilde = delta_x.abs() / d_i;
                    let defl = deflection_field[[fi, t, iy, iz]];
                    let r_tilde =
                        ((y_pt - y_i_loc - defl).powi(2) + (z_pt - z_i_loc).powi(2)).sqrt() / d_i;

                    let ti_i = turbulence_intensity[[fi, ii, iy, iz]];
                    let sigma_n = wake_expansion(delta_x, ct_i, ti_i, d_i);

                    let n_exp = params.a_f * (params.b_f * x_tilde).exp() + params.c_f;
                    let a1 = 2.0_f64.powf(2.0 / n_exp - 1.0);
                    let a2 = 2.0_f64.powf(4.0 / n_exp - 2.0);

                    let slb = sum_lbda[[fi, t, iy, iz]];
                    let gamma_val = gamma_fn(2.0 / n_exp);
                    let mut tmp = a2
                        - ((n_exp * ct_i) * yaw_i.to_radians().cos()
                            / (16.0
                                * gamma_val
                                * sigma_n.signum()
                                * sigma_n.abs().powf(4.0 / n_exp)
                                * (1.0 - slb).powi(2)));
                    if tmp < 0.0 {
                        tmp = 0.0;
                    }

                    let c_val = (a1 - tmp.sqrt()) * (1.0 - slb);
                    ctmp[[ii, fi, t, iy, iz]] = c_val;

                    let y_r = y_pt - y_i_loc;
                    let x_r = y_r * yaw_i.to_radians().tan() + x_i_loc;

                    let vel_def = c_val
                        * (-r_tilde.powf(n_exp) / (2.0 * sigma_n * sigma_n)).exp()
                        * if x_pt - x_r >= 0.1 { 1.0 } else { 0.0 };

                    turb_u_wake[[fi, t, iy, iz]] += src_avg_vel[fi] * vel_def;
                }
            }
        }
    }

    (turb_u_wake, ctmp)
}

/// Gamma function via the Lanczos approximation (g=7, 9 coefficients)
fn gamma_fn(x: Float) -> Float {
    const G: Float = 7.0;
    const COEF: [Float; 9] = [
        0.9999999999998099,
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI as Float
            / ((std::f64::consts::PI as Float * x).sin() * gamma_fn(1.0 - x))
    } else {
        let z = x - 1.0;
        let t = z + G + 0.5;
        let mut a = COEF[0];
        for (i, &c) in COEF.iter().enumerate().skip(1) {
            a += c / (z + i as Float);
        }
        (2.0 * std::f64::consts::PI as Float).sqrt() * t.powf(z + 0.5) * (-t).exp() * a
    }
}

impl VelocityModel for CumulativeCurlVelocityDeficit {
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
        hub_heights: &Array2,
        rotor_diameters: &Array2,
        turbine_index: usize,
        _model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let n_y = shape[2];
        let n_z = shape[3];

        let mut velocity_deficit = Array::zeros((n_findex, n_turbines, n_y, n_z));

        for fi in 0..n_findex {
            let hub_height = hub_heights[[fi, 0]];
            let rotor_diameter = rotor_diameters[[fi, 0]];
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

            // Convert yaw angle to radians
            let yaw_rad = _yaw_angle[[fi, turbine_index]].to_radians();
            let ct_i = thrust_coefficient[[fi, turbine_index]];

            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        if x_point <= x_wake_source {
                            continue;
                        }

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        // Wake center with deflection evaluated at this point
                        let wake_center_y = y_wake_source + deflection_field[[fi, ti, iy, iz]];
                        let dy = y_point - wake_center_y;
                        let dz = z_point - hub_height;

                        // Radial distance from wake center
                        let r = (dy.powi(2) + dz.powi(2)).sqrt();

                        // Per-point TI (source's slice of the evolving 4D field)
                        let ti_i = turbulence_intensity[[fi, turbine_index, iy, iz]];

                        // Wake expansion in D-units, converted to meters
                        let downstream_x = x_point - x_wake_source;
                        let sigma =
                            calculate_wake_spread(downstream_x, ti_i, ct_i, rotor_diameter, self)
                                * rotor_diameter;

                        // Apply curl correction for yawed turbines
                        let curl_correction = self.calculate_curl_correction(
                            downstream_x,
                            yaw_rad,
                            ct_i,
                            rotor_diameter,
                            dz,
                        );

                        let sigma_effective = sigma * curl_correction;

                        // Calculate deficit
                        let deficit = gaussian_deficit(r, sigma_effective, ct_i);

                        // Wake curl causes velocity to increase at edges
                        let curl_deficit = deficit * curl_correction;

                        velocity_deficit[[fi, ti, iy, iz]] = curl_deficit;
                    }
                }
            }
        }

        Ok(velocity_deficit)
    }
}

impl CumulativeCurlVelocityDeficit {
    /// Calculate curl correction factor for yawed turbines
    fn calculate_curl_correction(
        &self,
        x: Float,
        yaw_rad: Float,
        ct: Float,
        rotor_diameter: Float,
        dz: Float,
    ) -> Float {
        if yaw_rad.abs() < 0.01 {
            return 1.0;
        }

        let x_d = x / rotor_diameter;

        // Curl effect parameters
        // The curl effect is stronger near the turbine and decreases downstream
        let curl_strength = 0.5 * (1.0 - (-x_d / 5.0).exp());

        // Yaw-induced curl
        let yaw_curl = yaw_rad * ct.sqrt() * curl_strength;

        // Apply curl correction - increases deficit on one side, decreases on other
        let curl_factor = 1.0 + yaw_curl * (dz / (x + 1.0));

        curl_factor.max(0.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;

    #[test]
    fn test_cumulative_curl_creation() {
        let curl =
            CumulativeCurlVelocityDeficit::new(1.0, 0.179, 0.012, 0.056, 0.133, 3.11, -0.68, 2.41);
        assert_eq!(curl.wec_factor, 1.0);
        assert_eq!(curl.a_s, 0.179);
    }

    #[test]
    fn test_cumulative_curl_default() {
        let curl = CumulativeCurlVelocityDeficit::default();
        assert_eq!(curl.wec_factor, 1.0);
        assert_eq!(curl.a_f, 3.11);
        assert_eq!(curl.c_f, 2.41);
    }

    #[test]
    fn test_cumulative_curl_prepare_function() {
        let curl = CumulativeCurlVelocityDeficit::default();
        let result = curl.prepare_function(&fake_grid(), &fake_flow_field());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn test_wake_spread() {
        let params = CumulativeCurlVelocityDeficit::default();
        let sigma = calculate_wake_spread(100.0, 0.06, 0.8, 126.0, &params);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_gaussian_deficit() {
        let deficit = gaussian_deficit(10.0, 20.0, 0.8);
        assert!(deficit > 0.0);
        assert!(deficit < 0.5);
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
