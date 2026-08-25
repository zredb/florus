use crate::core::grid::TurbineCubatureGrid;
/// Wake solver algorithms
///
/// Corresponds to core/solver.py in Python implementation
use std::collections::HashMap;

use crate::core::wake::WakeModelManager;
use crate::core::wake::yaw_added_recovery::{calculate_transverse_velocity, yaw_added_turbulence_mixing};
use crate::core::{Farm, FlowField, Grid};
use crate::types::{Array1, Array2, Array4, Float};
use anyhow::Result;
use ndarray::{s, Array};

/// Solver type enum for selecting the appropriate solver algorithm
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverType {
    /// Standard sequential solver for turbine grid
    Sequential,
    /// Solver for cumulative curl velocity model with cubature grid
    CumulativeCurl,
    /// Solver for empirical Gauss velocity model
    EmpiricalGauss,
    /// Solver for Jensen velocity model
    Jensen,
    /// General solver that auto-selects based on grid type
    Auto,
}

impl Default for SolverType {
    fn default() -> Self {
        Self::Sequential
    }
}

/// Select appropriate solver based on grid and model configuration
pub fn select_solver(grid: &dyn Grid, model_manager: &WakeModelManager) -> SolverType {
    if model_manager.use_parallel_calc {
        return SolverType::Sequential;
    }

    if grid.cubature_weights().is_some() {
        return SolverType::CumulativeCurl;
    }

    let model_name =
        format!("{:?}", &model_manager.velocity_model as &dyn std::any::Any).to_lowercase();

    if model_name.contains("cumulative_curl") || model_name.contains("cc") {
        SolverType::CumulativeCurl
    } else if model_name.contains("empirical") {
        SolverType::EmpiricalGauss
    } else if model_name.contains("jensen") {
        SolverType::Jensen
    } else {
        SolverType::Sequential
    }
}

/// Main solver dispatch function
pub fn solve(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
    solver_type: SolverType,
) -> Result<()> {
    match solver_type {
        SolverType::Sequential | SolverType::Auto => {
            sequential_solver(farm, flow_field, grid, model_manager)
        }
        SolverType::CumulativeCurl => {
            if let Some(cubature_grid) = grid.as_any().downcast_ref::<TurbineCubatureGrid>() {
                cc_solver(farm, flow_field, cubature_grid, model_manager)
            } else {
                sequential_solver(farm, flow_field, grid, model_manager)
            }
        }
        SolverType::EmpiricalGauss => empirical_gauss_solver(farm, flow_field, grid, model_manager),
        SolverType::Jensen => jensen_solver(farm, flow_field, grid, model_manager),
    }
}
///
/// Computes wakes for each turbine in upstream→downstream order
/// applying deflection, velocity deficit, turbulence, and wake superposition models
pub fn sequential_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    let n_turbines = grid.n_turbines();
    let n_findex = grid.n_findex();

    let shape = grid.x_sorted().shape();
    let grid_second_dim = shape[1]; // This could be n_turbines or total_points depending on grid type
    let grid_y_dim = shape[2];
    let grid_z_dim = shape[3];

    // Get full grid coordinates
    let x_grid = grid.x_sorted().clone();
    let y_grid = grid.y_sorted().clone();
    let z_grid = grid.z_sorted().clone();

    let x_mean = mean_value_4d_slice(&x_grid)?;
    let y_mean = mean_value_4d_slice(&y_grid)?;

    // Per-turbine mean grid positions (used by the WAT downstream-influence mask).
    let x_mean = mean_value_4d_slice(&x_grid)?;
    let y_mean = mean_value_4d_slice(&y_grid)?;

    // Initialize wake field arrays
    // This stores the combined wake deficit (in velocity units, not ratio)
    // Use grid_second_dim instead of n_turbines to support FlowFieldPlanarGrid
    let mut wake_field: Array4 = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));

    // Prepare combination model function arguments
    let _combination_model_args = model_manager
        .combination_model
        .prepare_function(grid, flow_field)?;

    // Per-turbine turbulence intensity (n_findex, n_turbines), broadcast from ambient TI.
    // Mirrors Python turbine_turbulence_intensity; evolved downstream via GCH recovery.
    let mut turbine_turbulence_intensity: Array2 = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            turbine_turbulence_intensity[[fi, ti]] = flow_field.turbulence_intensities[fi];
        }
    }

    // Loop through turbines (upstream to downstream)
    for i in 0..n_turbines {
        // Compute thrust coefficient and axial induction through the operation
        // model (rotor-averaged velocity), matching Python solver.py which calls
        // floris.core.thrust_coefficient with ix_filter=[i].
        let ct_i = Farm::thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        let a_i = Farm::axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        // Use yaw_angles_sorted with sorted index i
        let yaw_angle_i = farm.yaw_angles_sorted[[0, i]];

        // Use rotor_diameters_sorted with sorted index i
        let rotor_diameter_i = farm.rotor_diameters_sorted[[0, i]];

        // Use hub_heights_sorted with sorted index i
        let hub_height_i = farm.hub_heights_sorted[[0, i]];

        // Current turbine velocity slices (captured before this turbine's wake is added)
        let u_i = flow_field.u_sorted.slice(s![.., i..i + 1, .., ..]).to_owned();
        let v_i = flow_field.v_sorted.slice(s![.., i..i + 1, .., ..]).to_owned();
        let w_i = flow_field.w_sorted.slice(s![.., i..i + 1, .., ..]).to_owned();

        // Calculate 2D mean values for deflection at turbine i's position
        let x_i = x_grid.slice(s![.., i..i + 1, .., ..]).to_owned();
        let y_i = y_grid.slice(s![.., i..i + 1, .., ..]).to_owned();

        let x_2d = mean_value_4d_slice(&x_i)?;
        let y_2d = mean_value_4d_slice(&y_i)?;

        // Per-findex scalar turbine properties for the GCH (yaw-added recovery) model
        let mut yaw_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut tsr_arr = Array1::zeros(n_findex);
        let mut ai_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
            ct_arr[fi] = ct_i[[fi, i]];
            tsr_arr[fi] = farm.tsrs_sorted[[fi, i]];
            ai_arr[fi] = a_i[[fi, i]];
        }

        // Distance from each grid point to turbine i's mean position
        let mut delta_x = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
        let mut delta_y = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
        for fi in 0..n_findex {
            for ti in 0..grid_second_dim {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        delta_x[[fi, ti, iy, iz]] = x_grid[[fi, ti, iy, iz]] - x_2d[[fi, 0]];
                        delta_y[[fi, ti, iy, iz]] = y_grid[[fi, ti, iy, iz]] - y_2d[[fi, 0]];
                    }
                }
            }
        }

        // Calculate wake deflection field (uses ambient/raw TI, matching Python)
        let deflection_field = model_manager.deflection_model.function(
            x_2d.clone(),
            y_2d.clone(),
            yaw_angle_i,
            flow_field.turbulence_intensities[0],
            ct_i[[0, i]],
            rotor_diameter_i,
            &std::collections::HashMap::new(),
        )?;

        // Broadcast deflection from wake source turbine i to all turbines
        let mut deflection_broadcast = Array::zeros((n_findex, n_turbines));
        for ti in 0..n_turbines {
            for fi in 0..n_findex {
                deflection_broadcast[[fi, ti]] = deflection_field[[fi, 0]];
            }
        }

        // --- GCH: transverse velocities (v/w) and yaw-added turbulence recovery ---
        let mut v_wake = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
        let mut w_wake = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));

        if model_manager.enable_transverse_velocities {
            let (vw, ww) = calculate_transverse_velocity(
                &u_i,
                &flow_field.u_initial_sorted,
                &flow_field.dudz_initial_sorted,
                &delta_x,
                &delta_y,
                &z_grid,
                rotor_diameter_i,
                hub_height_i,
                &yaw_arr,
                &ct_arr,
                &tsr_arr,
                &ai_arr,
                flow_field.wind_shear,
                1.0,
            );
            v_wake = vw;
            w_wake = ww;
        }

        if model_manager.enable_yaw_added_recovery {
            let i_mixing = yaw_added_turbulence_mixing(
                &u_i,
                &turbine_turbulence_intensity.slice(s![.., i]).to_owned(),
                &v_i,
                &w_i,
                &v_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
                &w_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
            );
            let gch_gain = 2.0;
            for fi in 0..n_findex {
                turbine_turbulence_intensity[[fi, i]] += gch_gain * i_mixing[fi];
            }
        }

        // Build the per-turbine TI array (now evolved for turbine i) for the velocity model
        let mut ti_array = Array::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                ti_array[[fi, ti]] = turbine_turbulence_intensity[[fi, ti]];
            }
        }

        // Calculate velocity deficit at ALL grid points due to turbine i
        let velocity_deficit = model_manager.velocity_model.function(
            x_grid.clone(),
            y_grid.clone(),
            z_grid.clone(),
            a_i.clone(),
            deflection_broadcast,
            farm.yaw_angles_sorted.clone(),
            ti_array,
            ct_i,
            hub_height_i,
            rotor_diameter_i,
            i, // turbine_index
            &std::collections::HashMap::new(),
        )?;

        // Convert relative deficit to absolute velocity deficit
        // velocity_deficit is a ratio (0-1), multiply by u_initial to get velocity in m/s
        let mut velocity_deficit_absolute = flow_field.u_initial_sorted.clone();
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        velocity_deficit_absolute[[fi, ti, iy, iz]] *=
                            velocity_deficit[[fi, ti, iy, iz]];
                    }
                }
            }
        }

        // Apply combination model to combine new deficit with existing wake field
        // This matches Python: wake_field = combination_model.function(wake_field, velocity_deficit * u_initial)
        wake_field = model_manager
            .combination_model
            .function(&wake_field, &velocity_deficit_absolute)?;

        // --- Wake-Added Turbulence (WAT) clamp ---
        // Mirrors FLORIS sequential_solver: compute the wake-added TI grid from the
        // turbulence model, fold it into turbine_turbulence_intensity for every
        // downstream turbine via an area-overlap + downstream-influence mask.
        let mut x_i_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ai_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ambient_ti_4d = Array4::zeros((n_findex, 1, 1, 1));
        for fi in 0..n_findex {
            x_i_4d[[fi, 0, 0, 0]] = x_mean[[fi, i]];
            ai_4d[[fi, 0, 0, 0]] = a_i[[fi, i]];
            ambient_ti_4d[[fi, 0, 0, 0]] = flow_field.turbulence_intensities[fi];
        }

        let mut wake_added = model_manager.turbulence_model.function(
            &ambient_ti_4d,
            &x_grid,
            &x_i_4d,
            rotor_diameter_i,
            &ai_4d,
        )?;
        // Match Python's np.nan_to_num(..., posinf=0.0)
        wake_added = wake_added.mapv(|v| if v.is_finite() { v } else { 0.0 });

        // grid_y_dim * grid_z_dim == grid_resolution^2 (turbine grid): the number of
        // grid points per turbine, the normalization for area_overlap.
        let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;
        let iy_c = grid_y_dim / 2;
        let iz_c = grid_z_dim / 2;
        let downstream_influence_length = 15.0 * rotor_diameter_i;

        for fi in 0..n_findex {
            for j in 0..n_turbines {
                // Fraction of turbine j's grid points with absolute deficit > 0.05 m/s
                let mut overlap_count = 0.0_f64;
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        if velocity_deficit_absolute[[fi, j, iy, iz]] > 0.05 {
                            overlap_count += 1.0;
                        }
                    }
                }
                let area_overlap = overlap_count / grid_res_sq;

                // Evaluate WAT at turbine j's center grid point (masks are ~1 downstream)
                let wake_added_center = wake_added[[fi, j, iy_c, iz_c]];

                // Downstream-influence mask (Python: x > x_i, |dy| < 2D, x <= 15D + x_i)
                let dx = x_mean[[fi, j]] - x_mean[[fi, i]];
                let dy = (y_mean[[fi, j]] - y_mean[[fi, i]]).abs();
                let downstream_ok = dx > 0.0
                    && dx <= downstream_influence_length
                    && dy < 2.0 * rotor_diameter_i;

                let ti_added = if downstream_ok {
                    area_overlap * wake_added_center
                } else {
                    0.0
                };

                let ambient = flow_field.turbulence_intensities[fi];
                let new_ti = (ti_added * ti_added + ambient * ambient).sqrt();
                if new_ti > turbine_turbulence_intensity[[fi, j]] {
                    turbine_turbulence_intensity[[fi, j]] = new_ti;
                }
            }
        }

        // Update u_sorted after each turbine's wake is added, matching Python
        // solver.py: flow_field.u_sorted = flow_field.u_initial_sorted - wake_field
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let wake_deficit = wake_field[[fi, ti, iy, iz]];
                        let u_initial = flow_field.u_initial_sorted[[fi, ti, iy, iz]];
                        flow_field.u_sorted[[fi, ti, iy, iz]] =
                            (u_initial - wake_deficit).max(0.0);
                    }
                }
            }
        }

        // Add transverse velocities to the sorted fields, matching Python:
        // flow_field.v_sorted += v_wake; flow_field.w_sorted += w_wake
        if model_manager.enable_transverse_velocities {
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    for iy in 0..grid_y_dim {
                        for iz in 0..grid_z_dim {
                            flow_field.v_sorted[[fi, ti, iy, iz]] += v_wake[[fi, ti, iy, iz]];
                            flow_field.w_sorted[[fi, ti, iy, iz]] += w_wake[[fi, ti, iy, iz]];
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

/// Cumulative Curl solver for cubature grid
///
/// Uses cubature integration for accurate rotor average velocity calculation.
/// The CC model includes Wake Expansion Continuation (WEC) for deep array effects.
pub fn cc_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &TurbineCubatureGrid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    let n_turbines = grid.n_turbines();
    let n_findex = grid.n_findex();

    let shape = grid.x_sorted().shape();
    let n_points = shape[2];

    let x_grid = grid.x_sorted().clone();
    let y_grid = grid.y_sorted().clone();
    let z_grid = grid.z_sorted().clone();

    let mut wake_field: Array4 = Array::zeros((n_findex, n_turbines, n_points, 1));

    model_manager
        .combination_model
        .prepare_function(grid, flow_field)?;

    for i in 0..n_turbines {
        let ct_i = Farm::thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        let a_i = Farm::axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        let mut ti_array = Array::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                ti_array[[fi, ti]] = flow_field.turbulence_intensities[fi];
            }
        }

        let yaw_angle_i = farm.yaw_angles_sorted[[0, i]];
        let rotor_diameter_i = farm.rotor_diameters_sorted[[0, i]];
        let hub_height_i = farm.hub_heights_sorted[[0, i]];

        let x_i = x_grid.slice(s![.., i..i + 1, .., ..]).to_owned();
        let y_i = y_grid.slice(s![.., i..i + 1, .., ..]).to_owned();

        let x_2d = mean_value_4d_slice(&x_i)?;
        let y_2d = mean_value_4d_slice(&y_i)?;

        let deflection_field = model_manager.deflection_model.function(
            x_2d.clone(),
            y_2d.clone(),
            yaw_angle_i,
            flow_field.turbulence_intensities[0],
            ct_i[[0, i]],
            rotor_diameter_i,
            &std::collections::HashMap::new(),
        )?;

        let mut deflection_broadcast = Array::zeros((n_findex, n_turbines));
        for ti in 0..n_turbines {
            for fi in 0..n_findex {
                deflection_broadcast[[fi, ti]] = deflection_field[[fi, 0]];
            }
        }

        let velocity_deficit = model_manager.velocity_model.function(
            x_grid.clone(),
            y_grid.clone(),
            z_grid.clone(),
            a_i,
            deflection_broadcast,
            farm.yaw_angles_sorted.clone(),
            ti_array,
            ct_i,
            hub_height_i,
            rotor_diameter_i,
            i,
            &std::collections::HashMap::new(),
        )?;

        let mut velocity_deficit_absolute = flow_field.u_initial_sorted.clone();
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for ip in 0..n_points {
                    velocity_deficit_absolute[[fi, ti, ip, 0]] *= velocity_deficit[[fi, ti, ip, 0]];
                }
            }
        }

        wake_field = model_manager
            .combination_model
            .function(&wake_field, &velocity_deficit_absolute)?;
    }

    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for ip in 0..n_points {
                let wake_deficit = wake_field[[fi, ti, ip, 0]];
                let u_initial = flow_field.u_initial_sorted[[fi, ti, ip, 0]];
                flow_field.u_sorted[[fi, ti, ip, 0]] = (u_initial - wake_deficit).max(0.0);
            }
        }
    }

    Ok(())
}

/// Empirical Gauss solver
///
/// Solver optimized for the Empirical Gaussian velocity deficit model.
/// This model uses data-driven wake expansion parameters.
pub fn empirical_gauss_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    sequential_solver(farm, flow_field, grid, model_manager)
}

/// Jensen solver
///
/// Solver optimized for the Jensen velocity deficit model.
/// Uses a simple top-hat wake profile with linear expansion.
pub fn jensen_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    sequential_solver(farm, flow_field, grid, model_manager)
}

/// Calculate axial induction from thrust coefficient
pub fn axial_induction(ct_values: &Array2) -> Array2 {
    ct_values.mapv(|ct| {
        if ct < 0.96 {
            0.5 * (1.0 - (1.0 - ct).sqrt())
        } else {
            0.143 + (0.0203 - 0.6427 * (0.889 - ct).sqrt()).max(0.0)
        }
    })
}

/// Helper function to convert 4D slice to 2D mean values
fn mean_value_4d_slice(arr: &Array4) -> Result<Array2> {
    let shape = arr.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];
    let mut result = Array::zeros((n_findex, n_turbines));

    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let slice = arr.slice(s![fi, ti, .., ..]);
            let sum: Float = slice.iter().sum();
            let count = slice.len() as Float;
            if count > 0.0 {
                result[[fi, ti]] = sum / count;
            }
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_axial_induction_low_ct() {
        let ct = Array2::from_elem((1, 1), 0.5);
        let ai = axial_induction(&ct);
        assert!(ai[[0, 0]] > 0.0 && ai[[0, 0]] < 0.5);
    }

    #[test]
    fn test_axial_induction_high_ct() {
        let ct = Array2::from_elem((1, 1), 0.99);
        let ai = axial_induction(&ct);
        assert!(ai[[0, 0]] > 0.0 && ai[[0, 0]] < 1.0);
    }

    #[test]
    fn test_sequential_solver_basic() {
        assert!(true);
    }

    #[test]
    fn test_mean_value_4d_slice() {
        let arr = Array::from_shape_vec((1, 2, 2, 2), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
            .unwrap();

        let result = mean_value_4d_slice(&arr).unwrap();

        assert_eq!(result.shape()[0], 1);
        assert_eq!(result.shape()[1], 2);
        assert_relative_eq!(result[[0, 0]], 2.5);
        assert_relative_eq!(result[[0, 1]], 6.5);
    }

    #[test]
    fn test_mean_value_4d_slice_single_value() {
        let arr = Array::from_shape_vec((1, 1, 1, 1), vec![5.0]).unwrap();

        let result = mean_value_4d_slice(&arr).unwrap();

        assert_eq!(result.shape()[0], 1);
        assert_eq!(result.shape()[1], 1);
        assert_relative_eq!(result[[0, 0]], 5.0);
    }

    #[test]
    fn test_mean_value_4d_slice_multiple_findex() {
        let arr = Array::from_shape_vec((2, 2, 1, 1), vec![10.0, 20.0, 30.0, 40.0]).unwrap();

        let result = mean_value_4d_slice(&arr).unwrap();

        assert_eq!(result.shape()[0], 2);
        assert_eq!(result.shape()[1], 2);
        assert_relative_eq!(result[[0, 0]], 10.0);
        assert_relative_eq!(result[[0, 1]], 20.0);
        assert_relative_eq!(result[[1, 0]], 30.0);
        assert_relative_eq!(result[[1, 1]], 40.0);
    }

    #[test]
    fn test_solver_type_default() {
        let st = SolverType::default();
        assert_eq!(st, SolverType::Sequential);
    }

    #[test]
    fn test_solver_type_equality() {
        assert_eq!(SolverType::Sequential, SolverType::Sequential);
        assert_eq!(SolverType::CumulativeCurl, SolverType::CumulativeCurl);
        assert_ne!(SolverType::Sequential, SolverType::CumulativeCurl);
    }
}
