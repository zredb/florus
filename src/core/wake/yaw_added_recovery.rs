//! Yaw-added recovery (GCH) transverse velocity and turbulence mixing
//!
//! Ports of `floris/floris/core/wake_deflection/gauss.py`:
//! `gamma`, `calculate_transverse_velocity`, `yaw_added_turbulence_mixing`.
//! These implement the GCH wake-added turbulence recovery mechanism.

use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::{cosd, sind};
use ndarray::Array;
use std::f64::consts::PI;

/// Numerical epsilon to prevent divide-by-zero, mirrors `BaseModel.NUM_EPS`.
const NUM_EPS: Float = 0.001;

/// Vortex circulation strength.
///
/// Mirrors `gamma(D, velocity, Uinf, Ct, scale=1.0)` in gauss.py.
/// `scale * (pi / 8) * D * velocity * Uinf * Ct` (the `cosd(yaw)` factor is
/// folded into `Ct` upstream, as noted in the Python source).
pub fn gamma(d: Float, velocity: Float, uinf: Float, ct: Float, scale: Float) -> Float {
    scale * (PI / 8.0) * d * velocity * uinf * ct
}

/// Compute transverse (spanwise `V`) and vertical (`W`) velocity components
/// induced by the vortices at the current turbine, for all downstream points.
///
/// Port of `calculate_transverse_velocity(...)`. All grid arrays are
/// `(n_findex, n_turbines, n_y, n_z)`. The per-turbine scalar parameters
/// (`yaw`, `ct`, `tsr`, `axial_induction`) are provided per findex for the
/// current turbine `i`.
pub fn calculate_transverse_velocity(
    u_i: &Array4,
    u_initial: &Array4,
    dudz_initial: &Array4,
    delta_x: &Array4,
    delta_y: &Array4,
    z: &Array4,
    rotor_diameters: &Array2,
    hub_heights: &Array2,
    yaw: &Array1,
    ct_i: &Array1,
    tsr_i: &Array1,
    axial_induction_i: &Array1,
    wind_shear: Float,
    scale: Float,
) -> (Array4, Array4) {
    let shape = delta_x.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];
    let n_y = shape[2];
    let n_z = shape[3];

    let mut v = Array::zeros((n_findex, n_turbines, n_y, n_z));
    let mut w = Array::zeros((n_findex, n_turbines, n_y, n_z));

    let eps_gain = 0.2;

    for fi in 0..n_findex {
        let d = rotor_diameters[[fi, 0]];
        let hh = hub_heights[[fi, 0]];
        let eps = eps_gain * d;

        let yaw_f = yaw[fi];
        let ct_f = ct_i[fi];
        let tsr_f = tsr_i[fi];
        let ai_f = axial_induction_i[fi];

        // Uinf = mean of u_initial over (n_turbines, n_y, n_z) for this findex
        let mut sum = 0.0;
        let mut count = 0.0;
        for ti in 0..n_turbines {
            for iy in 0..n_y {
                for iz in 0..n_z {
                    sum += u_initial[[fi, ti, iy, iz]];
                    count += 1.0;
                }
            }
        }
        let uinf = sum / count;

        let vel_top = ((hh + d / 2.0) / hh).powf(wind_shear);
        let gamma_top = sind(yaw_f) * cosd(yaw_f) * gamma(d, vel_top, uinf, ct_f, scale);
        let vel_bottom = ((hh - d / 2.0) / hh).powf(wind_shear);
        let gamma_bottom =
            -1.0 * sind(yaw_f) * cosd(yaw_f) * gamma(d, vel_bottom, uinf, ct_f, scale);

        // turbine_average_velocity = cbrt(mean(u_i**3, axis=(2, 3)))
        let mut sum3 = 0.0;
        let mut count3 = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                let uu = u_i[[fi, 0, iy, iz]];
                sum3 += uu * uu * uu;
                count3 += 1.0;
            }
        }
        let turbine_average_velocity = (sum3 / count3).powf(1.0 / 3.0);
        let gamma_wake_rotation =
            0.25 * 2.0 * PI * d * (ai_f - ai_f * ai_f) * turbine_average_velocity / tsr_f;

        // Mixing-length based eddy viscosity for downstream decay
        let lmda = d / 8.0;
        let kappa = 0.41;

        for ti in 0..n_turbines {
            for iy in 0..n_y {
                for iz in 0..n_z {
                    let dx = delta_x[[fi, ti, iy, iz]];
                    let dy = delta_y[[fi, ti, iy, iz]];
                    let zz = z[[fi, ti, iy, iz]];

                    let lm = kappa * zz / (1.0 + kappa * zz / lmda);
                    let nu = lm * lm * dudz_initial[[fi, ti, iy, iz]].abs();
                    let decay = eps * eps / (4.0 * nu * dx / uinf + eps * eps);

                    let ylocs = dy + NUM_EPS;

                    // top vortex
                    let zt = zz - (hh + d / 2.0) + NUM_EPS;
                    let rt2 = ylocs * ylocs + zt * zt;
                    let core = 1.0 - (-rt2 / (eps * eps)).exp();
                    let v1 = (gamma_top * zt) / (2.0 * PI * rt2) * core * decay;
                    let w1 = (-1.0 * gamma_top * ylocs) / (2.0 * PI * rt2) * core * decay;

                    // bottom vortex
                    let zb = zz - (hh - d / 2.0) + NUM_EPS;
                    let rb2 = ylocs * ylocs + zb * zb;
                    let coreb = 1.0 - (-rb2 / (eps * eps)).exp();
                    let v2 = (gamma_bottom * zb) / (2.0 * PI * rb2) * coreb * decay;
                    let w2 = (-1.0 * gamma_bottom * ylocs) / (2.0 * PI * rb2) * coreb * decay;

                    // wake rotation vortex
                    let zc = zz - hh + NUM_EPS;
                    let rc2 = ylocs * ylocs + zc * zc;
                    let corec = 1.0 - (-rc2 / (eps * eps)).exp();
                    let v5 = (gamma_wake_rotation * zc) / (2.0 * PI * rc2) * corec * decay;
                    let w5 =
                        (-1.0 * gamma_wake_rotation * ylocs) / (2.0 * PI * rc2) * corec * decay;

                    // ground mirror - top vortex
                    let ztb = zz + (hh + d / 2.0) + NUM_EPS;
                    let r_tb2 = ylocs * ylocs + ztb * ztb;
                    let coretb = 1.0 - (-r_tb2 / (eps * eps)).exp();
                    let v3 = (-1.0 * gamma_top * ztb) / (2.0 * PI * r_tb2) * coretb * decay;
                    let w3 = (gamma_top * ylocs) / (2.0 * PI * r_tb2) * coretb * decay;

                    // ground mirror - bottom vortex
                    let zbb = zz + (hh - d / 2.0) + NUM_EPS;
                    let r_bb2 = ylocs * ylocs + zbb * zbb;
                    let corebb = 1.0 - (-r_bb2 / (eps * eps)).exp();
                    let v4 = (-1.0 * gamma_bottom * zbb) / (2.0 * PI * r_bb2) * corebb * decay;
                    let w4 = (gamma_bottom * ylocs) / (2.0 * PI * r_bb2) * corebb * decay;

                    // ground mirror - wake rotation vortex
                    let zcb = zz + hh + NUM_EPS;
                    let r_cb2 = ylocs * ylocs + zcb * zcb;
                    let corecb = 1.0 - (-r_cb2 / (eps * eps)).exp();
                    let v6 =
                        (-1.0 * gamma_wake_rotation * zcb) / (2.0 * PI * r_cb2) * corecb * decay;
                    let w6 = (gamma_wake_rotation * ylocs) / (2.0 * PI * r_cb2) * corecb * decay;

                    let mut vv = v1 + v2 + v3 + v4 + v5 + v6;
                    let mut ww = w1 + w2 + w3 + w4 + w5 + w6;

                    // No spanwise/vertical velocity upstream of the turbine
                    if dx < 0.0 {
                        vv = 0.0;
                        ww = 0.0;
                    }
                    if ww < 0.0 {
                        ww = 0.0;
                    }

                    v[[fi, ti, iy, iz]] = vv;
                    w[[fi, ti, iy, iz]] = ww;
                }
            }
        }
    }

    (v, w)
}

/// Compute the effective yaw angle added by secondary steering.
///
/// Port of `wake_added_yaw(...)` in `wake_deflection/gauss.py`. Given the
/// spanwise velocity induced by a yawed turbine, returns the yaw angle that
/// would have produced the same average spanwise velocity.
///
/// * `u_i` / `v_i`: current-turbine rotor grid (u, v), shape `[findex, 1, gy, gz]`
/// * `u_initial`: full turbine-grid freestream, shape `[findex, n_turbines, gy, gz]`
/// * `delta_y`: `turbine_grid.y_sorted[:, i:i+1] - y_i`, shape `[findex, 1, gy, gz]`
/// * `z_i`: `turbine_grid.z_sorted[:, i:i+1]`, shape `[findex, 1, gy, gz]`
/// * `ct_i` / `tip_speed_ratio` / `axial_induction_i`: per-findex `[findex, 1, 1, 1]`
///
/// Returns the added yaw per findex, shape `(n_findex,)`.
pub fn wake_added_yaw(
    u_i: &Array4,
    v_i: &Array4,
    u_initial: &Array4,
    delta_y: &Array4,
    z_i: &Array4,
    rotor_diameters: &Array2,
    hub_heights: &Array2,
    ct_i: &Array4,
    tip_speed_ratio: &Array4,
    axial_induction_i: &Array4,
    wind_shear: Float,
    scale: Float,
) -> Array1 {
    let n_findex = u_i.shape()[0];
    let n_y = u_i.shape()[2];
    let n_z = u_i.shape()[3];

    let eps_gain = 0.2;

    let mut y = Array1::zeros(n_findex);

    for fi in 0..n_findex {
        let d = rotor_diameters[[fi, 0]];
        let hh = hub_heights[[fi, 0]];
        let eps = eps_gain * d;

        let ct_f = ct_i[[fi, 0, 0, 0]];
        let tsr_f = tip_speed_ratio[[fi, 0, 0, 0]];
        let ai_f = axial_induction_i[[fi, 0, 0, 0]];

        // avg_v = mean(v_i, axis=(2, 3))  -- scalar for this findex
        let mut sum_v = 0.0;
        let mut cnt_v = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                sum_v += v_i[[fi, 0, iy, iz]];
                cnt_v += 1.0;
            }
        }
        let avg_v = sum_v / cnt_v;

        // Uinf = mean(u_initial, axis=(1, 2, 3))
        let mut sum_u = 0.0;
        let mut cnt_u = 0.0;
        for ti in 0..u_initial.shape()[1] {
            for iy in 0..n_y {
                for iz in 0..n_z {
                    sum_u += u_initial[[fi, ti, iy, iz]];
                    cnt_u += 1.0;
                }
            }
        }
        let uinf = sum_u / cnt_u;

        let vel_top = ((hh + d / 2.0) / hh).powf(wind_shear);
        let gamma_top = gamma(d, vel_top, uinf, ct_f, scale);

        let vel_bottom = ((hh - d / 2.0) / hh).powf(wind_shear);
        let gamma_bottom = -1.0 * gamma(d, vel_bottom, uinf, ct_f, scale);

        // turbine_average_velocity = cbrt(mean(u_i**3, axis=(2, 3)))
        let mut sum3 = 0.0;
        let mut cnt3 = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                let uu = u_i[[fi, 0, iy, iz]];
                sum3 += uu * uu * uu;
                cnt3 += 1.0;
            }
        }
        let turbine_average_velocity = (sum3 / cnt3).powf(1.0 / 3.0);
        let gamma_wake_rotation =
            0.25 * 2.0 * PI * d * (ai_f - ai_f * ai_f) * turbine_average_velocity / tsr_f;

        // Compute v_top, v_bottom, v_core averaged over (iy, iz)
        let mut sum_vtop = 0.0;
        let mut sum_vbottom = 0.0;
        let mut sum_vcore = 0.0;
        let mut cnt2 = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                let ylocs = delta_y[[fi, 0, iy, iz]] + NUM_EPS;

                // top vortex
                let zt = z_i[[fi, 0, iy, iz]] - (hh + d / 2.0) + NUM_EPS;
                let rt2 = ylocs * ylocs + zt * zt;
                let core_top = 1.0 - (-rt2 / (eps * eps)).exp();
                sum_vtop += (gamma_top * zt) / (2.0 * PI * rt2) * core_top;

                // bottom vortex
                let zb = z_i[[fi, 0, iy, iz]] - (hh - d / 2.0) + NUM_EPS;
                let rb2 = ylocs * ylocs + zb * zb;
                let core_bottom = 1.0 - (-rb2 / (eps * eps)).exp();
                sum_vbottom += (gamma_bottom * zb) / (2.0 * PI * rb2) * core_bottom;

                // wake rotation vortex
                let zc = z_i[[fi, 0, iy, iz]] - hh + NUM_EPS;
                let rc2 = ylocs * ylocs + zc * zc;
                let core_core = 1.0 - (-rc2 / (eps * eps)).exp();
                sum_vcore += (gamma_wake_rotation * zc) / (2.0 * PI * rc2) * core_core;

                cnt2 += 1.0;
            }
        }
        let v_top = sum_vtop / cnt2;
        let v_bottom = sum_vbottom / cnt2;
        let v_core = sum_vcore / cnt2;

        // Cap the effective yaw values between -45 and 45 degrees
        let mut val = 2.0 * (avg_v - v_core) / (v_top + v_bottom);
        if val < -1.0 {
            val = -1.0;
        }
        if val > 1.0 {
            val = 1.0;
        }
        y[fi] = (0.5 * val.asin()).to_degrees();
    }

    y
}

/// Compute the turbulence intensity added by wake mixing (yaw-added recovery).
///
/// Port of `yaw_added_turbulence_mixing(...)`. Returns `I_mixing` per findex
/// `(n_findex,)`. `i_i` is the per-findex ambient TI for the current turbine.
pub fn yaw_added_turbulence_mixing(
    u_i: &Array4,
    i_i: &Array1,
    v_i: &Array4,
    w_i: &Array4,
    turb_v_i: &Array4,
    turb_w_i: &Array4,
) -> Array1 {
    let n_findex = u_i.shape()[0];
    let n_y = u_i.shape()[2];
    let n_z = u_i.shape()[3];
    let mut i_mixing = Array::zeros(n_findex);

    let two_thirds = 2.0 / 3.0;
    for fi in 0..n_findex {
        // average_u_i = cbrt(mean(u_i**3, axis=(1, 2, 3)))
        let mut sum3 = 0.0;
        let mut cnt = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                let uu = u_i[[fi, 0, iy, iz]];
                sum3 += uu * uu * uu;
                cnt += 1.0;
            }
        }
        let average_u_i = (sum3 / cnt).powf(1.0 / 3.0);
        let ii = i_i[fi];

        // Convert ambient TI to TKE (eq 24)
        let k = (average_u_i * ii).powi(2) / two_thirds;
        let u_term = (2.0 * k).sqrt();

        // v_term/w_term = mean(v_i + turb_v_i, axis=(1, 2, 3))
        let mut sum_v = 0.0;
        let mut sum_w = 0.0;
        let mut cnt2 = 0.0;
        for iy in 0..n_y {
            for iz in 0..n_z {
                sum_v += v_i[[fi, 0, iy, iz]] + turb_v_i[[fi, 0, iy, iz]];
                sum_w += w_i[[fi, 0, iy, iz]] + turb_w_i[[fi, 0, iy, iz]];
                cnt2 += 1.0;
            }
        }
        let v_term = sum_v / cnt2;
        let w_term = sum_w / cnt2;

        // New TKE (eq 23)
        let k_total = 0.5 * (u_term * u_term + v_term * v_term + w_term * w_term);

        // Convert TKE back to TI
        let i_total = ((two_thirds * k_total).sqrt()) / average_u_i;

        // Remove ambient from total TI, leaving only the TI due to mixing
        i_mixing[fi] = i_total - ii;
    }

    i_mixing
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gamma_sign_and_magnitude() {
        // gamma should be positive for positive inputs; scales linearly with D and Ct
        let g1 = gamma(100.0, 10.0, 8.0, 0.8, 1.0);
        let g2 = gamma(200.0, 10.0, 8.0, 0.8, 1.0);
        assert!((g2 - 2.0 * g1).abs() < 1e-9);
        let g3 = gamma(100.0, 10.0, 8.0, 0.8, 2.0);
        assert!((g3 - 2.0 * g1).abs() < 1e-9);
    }

    #[test]
    fn test_yaw_added_turbulence_mixing_zero_when_no_crossflow() {
        // With zero cross-flow velocities and ambient TI, I_mixing should be ~0
        let n_findex = 1;
        let n_y = 3;
        let n_z = 3;
        let u_i = Array::from_elem((n_findex, 1, n_y, n_z), 8.0);
        let i_i = Array::from_vec(vec![0.06]);
        let v_i = Array::zeros((n_findex, 1, n_y, n_z));
        let w_i = Array::zeros((n_findex, 1, n_y, n_z));
        let turb_v_i = Array::zeros((n_findex, 1, n_y, n_z));
        let turb_w_i = Array::zeros((n_findex, 1, n_y, n_z));
        let mixing = yaw_added_turbulence_mixing(&u_i, &i_i, &v_i, &w_i, &turb_v_i, &turb_w_i);
        assert!(mixing[0].abs() < 1e-9);
    }
}
