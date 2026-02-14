use crate::types::{Array1, Array2, Float};
use crate::floris_model::FlorisModel;
use ndarray::{Array, Axis};
use std::f64::consts::PI;

pub const POWER_SETPOINT_DEFAULT: Float = 5_000_000.0;
pub const POWER_SETPOINT_DISABLED: Float = 0.0;

/// Compute the turbine 'load turbulence intensity' (LTI) for the current layout.
///
/// LTI represents turbulence intensity used in load calculations and follows the
/// method of computing wake added turbulence described in Annex E of the IEC 61400-1 Ed. 4
/// standard. In principle this can be the same as the turbulence models used in the wake
/// velocity and deflection models within FLORIS, but for consistency with the IEC standard
/// is computed separately here.
///
/// # Arguments
///
/// * `fmodel` - FlorisModel object
/// * `ambient_lti` - Ambient 'load' turbulence intensity for each findex
/// * `wake_slope` - Wake slope, lateral expansion per unit downstream distance (default: 0.3)
/// * `max_dist_d` - Maximum distance downstream of a turbine beyond which wake
///   added turbulence is assumed to be zero, in rotor diameters (default: 10.0)
///
/// # Returns
///
/// Array of load turbulence intensity for each findex and turbine
///
/// # Errors
///
/// Returns error if FlorisModel has not been run yet
pub fn compute_lti(
    fmodel: &FlorisModel,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
) -> anyhow::Result<Array2> {
    if !fmodel.state.initialized {
        anyhow::bail!("FlorisModel must be run before computing load turbulence intensity");
    }

    let grid = fmodel.grid.as_ref()
        .ok_or_else(|| anyhow::anyhow!("Grid not initialized"))?;

    let n_findex = fmodel.flow_field.n_findex;
    let n_turbines = fmodel.farm.rotor_diameters.len();

    if ambient_lti.len() != n_findex {
        anyhow::bail!(
            "ambient_lti must have length n_findex ({}), got {}",
            n_findex,
            ambient_lti.len()
        );
    }

    let d = fmodel.farm.rotor_diameters[0];
    let sorted_indices = grid.sorted_indices();
    let x_sorted = grid.x_sorted();
    let y_sorted = grid.y_sorted();

    let mut lti = Array::from_shape_vec((n_findex, n_turbines), vec![0.0; n_findex * n_turbines])?;

    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            lti[[fi, ti]] = ambient_lti[fi];
        }
    }

    let cts = fmodel.get_turbine_thrust_coefficients();
    let ambient_wind_speeds = &fmodel.flow_field.wind_speeds;

    let mut x_sorted_mean = Array::zeros((n_findex, n_turbines));
    let mut y_sorted_mean = Array::zeros((n_findex, n_turbines));

    for fi in 0..n_findex {
        let grid_res = grid.resolution();
        for ti in 0..n_turbines {
            let mut sum_x = 0.0;
            let mut sum_y = 0.0;
            for i in 0..grid_res {
                for j in 0..grid_res {
                    sum_x += x_sorted[[fi, ti, i, j]];
                    sum_y += y_sorted[[fi, ti, i, j]];
                }
            }
            let n_points = (grid_res * grid_res) as Float;
            x_sorted_mean[[fi, ti]] = sum_x / n_points;
            y_sorted_mean[[fi, ti]] = sum_y / n_points;
        }
    }

    let mut ct_sorted = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let sorted_idx = sorted_indices[[fi, ti]] as usize;
            ct_sorted[[fi, ti]] = cts[[fi, sorted_idx]];
        }
    }

    for t in 0..n_turbines {
        for fi in 0..n_findex {
            let x_t = x_sorted_mean[[fi, t]];
            let y_t = y_sorted_mean[[fi, t]];
            let ct_t = ct_sorted[[fi, t]];
            let ambient_ws = ambient_wind_speeds[fi];

            for tj in 0..n_turbines {
                if tj == t {
                    continue;
                }

                let x_j = x_sorted_mean[[fi, tj]];
                let y_j = y_sorted_mean[[fi, tj]];
                let dx = x_j - x_t;
                let dy = y_j - y_t;

                if dx <= 0.0 {
                    continue;
                }

                let wake_width = d + wake_slope * dx;
                if dy.abs() > wake_width {
                    continue;
                }

                let distance = (dx * dx + dy * dy).sqrt();
                if distance >= d * max_dist_d {
                    continue;
                }

                let ws_std_wake_add = ambient_ws / (1.5 + 0.8 * (distance / d) / ct_t.sqrt());
                let lti_update = ((ws_std_wake_add.powi(2) + (ambient_lti[fi] * ambient_ws).powi(2)).sqrt()) / ambient_ws;

                let sorted_idx = tj;
                if lti_update > lti[[fi, sorted_idx]] {
                    lti[[fi, sorted_idx]] = lti_update;
                }
            }
        }
    }

    let mut lti_unsorted = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let sorted_idx = sorted_indices[[fi, ti]] as usize;
            lti_unsorted[[fi, ti]] = lti[[fi, sorted_idx]];
        }
    }

    Ok(lti_unsorted)
}

pub fn compute_turbine_voc(
    fmodel: &FlorisModel,
    a: Float,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
    exp_ws_std: Float,
    exp_thrust: Float,
) -> anyhow::Result<Array2> {
    let n_findex = fmodel.flow_field.n_findex;
    let n_turbines = fmodel.farm.rotor_diameters.len();

    let ambient_wind_speeds = Array::from_shape_vec((n_findex, n_turbines), {
        let mut v = Vec::with_capacity(n_findex * n_turbines);
        for &ws in &fmodel.flow_field.wind_speeds {
            for _ in 0..n_turbines {
                v.push(ws);
            }
        }
        v
    })?;

    let d = fmodel.farm.rotor_diameters[0];
    let area = PI * (d / 2.0).powi(2);

    let cts = fmodel.get_turbine_thrust_coefficients();
    let air_density = fmodel.flow_field.air_density;

    let mut thrust = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let ws = ambient_wind_speeds[[fi, ti]];
            let ct = cts[[fi, ti]];
            thrust[[fi, ti]] = 0.5 * air_density * area * ct * ws.powi(2);
        }
    }

    let load_ti = compute_lti(fmodel, ambient_lti, wake_slope, max_dist_d)?;

    let ws_std = &ambient_wind_speeds * &load_ti;

    let mut voc = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            voc[[fi, ti]] = a * ws_std[[fi, ti]].powf(exp_ws_std) * thrust[[fi, ti]].powf(exp_thrust);
        }
    }

    Ok(voc)
}

pub fn compute_farm_voc(
    fmodel: &FlorisModel,
    a: Float,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
    exp_ws_std: Float,
    exp_thrust: Float,
) -> anyhow::Result<Array1> {
    let turbine_voc = compute_turbine_voc(
        fmodel,
        a,
        ambient_lti,
        wake_slope,
        max_dist_d,
        exp_ws_std,
        exp_thrust,
    )?;

    let farm_voc = turbine_voc.sum_axis(Axis(1));
    Ok(farm_voc)
}

pub fn compute_farm_revenue(fmodel: &FlorisModel) -> anyhow::Result<Array1> {
    if !fmodel.state.initialized {
        anyhow::bail!("FlorisModel must be run before computing farm revenue");
    }

    let farm_power = fmodel.get_farm_power();
    let revenue = farm_power.mapv(|p| p * 1.0);

    Ok(revenue)
}

pub fn compute_net_revenue(
    fmodel: &FlorisModel,
    _a: Float,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
    exp_ws_std: Float,
    exp_thrust: Float,
) -> anyhow::Result<Array1> {
    let revenue = compute_farm_revenue(fmodel)?;

    let farm_voc = compute_farm_voc(
        fmodel,
        1.0,
        ambient_lti,
        wake_slope,
        max_dist_d,
        exp_ws_std,
        exp_thrust,
    )?;

    Ok(&revenue - &farm_voc)
}

pub fn find_a_to_satisfy_rev_voc_ratio(
    fmodel: &FlorisModel,
    target_rev_voc_ratio: Float,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
    exp_ws_std: Float,
    exp_thrust: Float,
) -> anyhow::Result<Float> {
    let farm_revenue = compute_farm_revenue(fmodel)?;

    let farm_voc = compute_farm_voc(
        fmodel,
        1.0,
        ambient_lti,
        wake_slope,
        max_dist_d,
        exp_ws_std,
        exp_thrust,
    )?;

    let total_revenue: Float = farm_revenue.sum();
    let total_voc: Float = farm_voc.sum();

    Ok((total_revenue / total_voc) / target_rev_voc_ratio)
}

pub fn find_a_to_satisfy_target_voc_per_mw(
    fmodel: &FlorisModel,
    target_voc_per_mw_findex: Float,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
    exp_ws_std: Float,
    exp_thrust: Float,
) -> anyhow::Result<Float> {
    if !fmodel.state.initialized {
        anyhow::bail!("FlorisModel must be run before finding A for target cost/MW/findex");
    }

    let farm_power = fmodel.get_farm_power();

    let farm_voc = compute_farm_voc(
        fmodel,
        1.0,
        ambient_lti,
        wake_slope,
        max_dist_d,
        exp_ws_std,
        exp_thrust,
    )?;

    let total_power: Float = farm_power.sum();
    let total_voc: Float = farm_voc.sum();

    Ok(1e-6 * target_voc_per_mw_findex / (total_voc / total_power))
}

pub fn optimize_power_setpoints(
    fmodel: &mut FlorisModel,
    a: Float,
    ambient_lti: &[Float],
    wake_slope: Float,
    max_dist_d: Float,
    exp_ws_std: Float,
    exp_thrust: Float,
    power_setpoint_initial: Option<&Array2>,
    power_setpoint_levels: &[Float],
) -> anyhow::Result<(Array2, Array1)> {
    let operation_model = fmodel.get_operation_model();
    if operation_model != "mixed" && operation_model != "simple-derating" {
        anyhow::bail!(
            "Operation model must include derating (e.g., 'mixed' or 'simple-derating'), got '{}'",
            operation_model
        );
    }

    if fmodel.farm.turbine_type.len() > 1 {
        anyhow::bail!("Only one turbine type is currently supported for optimization");
    }

    let n_findex = fmodel.flow_field.n_findex;
    let n_turbines = fmodel.farm.rotor_diameters.len();

    let power_setpoint_initial = if let Some(initial) = power_setpoint_initial {
        initial.clone()
    } else {
        let max_power = fmodel.farm.turbine_map[0].turbine_type.power_curve().values[fmodel.farm.turbine_map[0].turbine_type.power_curve().values.len() - 1] * 1000.0;
        Array::from_elem((n_findex, n_turbines), max_power)
    };

    let mut power_setpoint_test = power_setpoint_initial.clone();
    let mut power_setpoint_opt = power_setpoint_initial.clone();

    let grid = fmodel.grid.as_ref()
        .ok_or_else(|| anyhow::anyhow!("Grid not initialized"))?;

    // Clone sorted_indices before mutable borrow of fmodel
    let sorted_indices = grid.sorted_indices().clone();
    let _x_sorted = grid.x_sorted();
    let _y_sorted = grid.y_sorted();

    // Drop all references to grid before mutable borrow
    let _ = grid;

    fmodel.farm.set_power_setpoints(power_setpoint_initial.clone());
    fmodel.run()?;

    let mut net_revenue_opt = compute_net_revenue(
        fmodel,
        a,
        ambient_lti,
        wake_slope,
        max_dist_d,
        exp_ws_std,
        exp_thrust,
    )?;

    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let sorted_idx = sorted_indices[[fi, ti]] as usize;

            for derating_level in power_setpoint_levels {
                power_setpoint_test[[fi, sorted_idx]] = *derating_level;

                fmodel.farm.set_power_setpoints(power_setpoint_test.clone());
                fmodel.run()?;

                let test_net_revenue = compute_net_revenue(
                    fmodel,
                    a,
                    ambient_lti,
                    wake_slope,
                    max_dist_d,
                    exp_ws_std,
                    exp_thrust,
                )?;

                let mut update_mask = Array::from_elem(n_findex, false);
                for f_idx in 0..n_findex {
                    update_mask[f_idx] = test_net_revenue[f_idx] > net_revenue_opt[f_idx];
                }

                for f_idx in 0..n_findex {
                    if !update_mask[f_idx] {
                        power_setpoint_test[[f_idx, sorted_idx]] = power_setpoint_opt[[f_idx, sorted_idx]];
                    }
                }

                for f_idx in 0..n_findex {
                    if update_mask[f_idx] {
                        power_setpoint_opt[[f_idx, sorted_idx]] = power_setpoint_test[[f_idx, sorted_idx]];
                        net_revenue_opt[f_idx] = test_net_revenue[f_idx];
                    }
                }
            }
        }
    }

    Ok((power_setpoint_opt, net_revenue_opt))
}
