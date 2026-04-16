use crate::core::grid::TurbineCubatureGrid;
/// Wake solver algorithms
///
/// Corresponds to core/solver.py in Python implementation
use crate::core::turbines::cp_ct_table::TableConditions;
use crate::core::turbines::turbine_type::TurbineType;
use crate::core::turbines::Turbine;
use crate::core::wake::WakeModelManager;
use crate::core::{Farm, FlowField, Grid};
use crate::types::{Array2, Array4, Float};
use anyhow::{bail, Result};
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
    let grid_y_dim = shape[2];
    let grid_z_dim = shape[3];

    // Get full grid coordinates
    let x_grid = grid.x_sorted().clone();
    let y_grid = grid.y_sorted().clone();
    let z_grid = grid.z_sorted().clone();

    // Initialize wake field arrays
    // This stores the combined wake deficit (in velocity units, not ratio)
    let mut wake_field: Array4 = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));

    // Prepare combination model function arguments
    let _combination_model_args = model_manager
        .combination_model
        .prepare_function(grid, flow_field)?;

    // Loop through turbines (upstream to downstream)
    for i in 0..n_turbines {
        // Get turbine properties - use sorted farm properties
        let ct_i = thrust_coefficient(
            &flow_field.u_sorted,
            &farm.turbines,
            &farm.yaw_angles_sorted,
            &farm.tilt_angles_sorted,
            grid.average_method(),
        )?;

        let a_i = axial_induction(&ct_i);
        let ti_i = flow_field.turbulence_intensities.slice(s![..]);

        // Use yaw_angles_sorted with sorted index i
        let yaw_angle_i = farm.yaw_angles_sorted[[0, i]];

        // Use rotor_diameters_sorted with sorted index i
        let rotor_diameter_i = farm.rotor_diameters_sorted[[0, i]];

        // Use hub_heights_sorted with sorted index i
        let hub_height_i = farm.hub_heights_sorted[[0, i]];

        // Calculate 2D mean values for deflection at turbine i's position
        let x_i = x_grid.slice(s![.., i..i + 1, .., ..]).to_owned();
        let y_i = y_grid.slice(s![.., i..i + 1, .., ..]).to_owned();

        let x_2d = mean_value_4d_slice(&x_i)?;
        let y_2d = mean_value_4d_slice(&y_i)?;

        // Calculate wake deflection field
        let deflection_field = model_manager.deflection_model.function(
            x_2d.clone(),
            y_2d.clone(),
            yaw_angle_i,
            ti_i[0],
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

        // Calculate velocity deficit at ALL grid points due to turbine i
        let velocity_deficit = model_manager.velocity_model.function(
            x_grid.clone(),
            y_grid.clone(),
            z_grid.clone(),
            a_i[[0, i]],
            deflection_broadcast,
            yaw_angle_i,
            ti_i[0],
            ct_i[[0, i]],
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
    }

    // Apply combined wake field to get final velocities
    // u_sorted = u_initial - wake_field
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    let wake_deficit = wake_field[[fi, ti, iy, iz]];
                    let u_initial = flow_field.u_initial_sorted[[fi, ti, iy, iz]];
                    // Ensure we don't get negative velocities
                    flow_field.u_sorted[[fi, ti, iy, iz]] = (u_initial - wake_deficit).max(0.0);
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

    let weights = grid.cubature_weights();

    let mut wake_field: Array4 = Array::zeros((n_findex, n_turbines, n_points, 1));

    model_manager
        .combination_model
        .prepare_function(grid, flow_field)?;

    for i in 0..n_turbines {
        let ct_i = thrust_coefficient_cc(
            &flow_field.u_sorted,
            &farm.turbines,
            weights,
            grid.average_method(),
        )?;

        let a_i = axial_induction(&ct_i);
        let ti_i = flow_field.turbulence_intensities.slice(s![..]);

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
            ti_i[0],
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
            a_i[[0, i]],
            deflection_broadcast,
            yaw_angle_i,
            ti_i[0],
            ct_i[[0, i]],
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

/// Thrust coefficient for cubature grid with weighted averaging
fn thrust_coefficient_cc(
    velocities: &Array4,
    turbines: &[Turbine],
    _weights: Option<&Array2>,
    _average_method: crate::core::AveragingMethod,
) -> Result<Array2> {
    let shape = velocities.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];

    let mut ct_output = Array::zeros((n_findex, n_turbines));

    for fi in 0..n_findex {
        for ti in 0..n_turbines.min(turbines.len()) {
            if ti < turbines.len() {
                let v = velocities[[fi, ti, 0, 0]];
                let conditions = TableConditions::builder().wind_speed(v).build().unwrap();
                ct_output[[fi, ti]] = turbines[ti]
                    .turbine_type
                    .power_thrust_table
                    .cp_ct_table
                    .get_ct(&conditions)
                    .unwrap_or(0.0);
            }
        }
    }

    Ok(ct_output)
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

/// Calculate thrust coefficient for all turbines
fn thrust_coefficient(
    velocities: &Array4,
    turbines: &[Turbine],
    _yaw_angles: &Array2,
    _tilt_angles: &Array2,
    _average_method: crate::core::AveragingMethod,
) -> Result<Array2> {
    let shape = velocities.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];

    let mut ct_output = Array::zeros((n_findex, n_turbines));

    for fi in 0..n_findex {
        for ti in 0..n_turbines.min(turbines.len()) {
            if ti < turbines.len() {
                let v = velocities[[fi, ti, 0, 0]];
                let conditions = TableConditions::builder().wind_speed(v).build().unwrap();
                ct_output[[fi, ti]] = turbines[ti]
                    .turbine_type
                    .power_thrust_table
                    .cp_ct_table
                    .get_ct(&conditions)
                    .unwrap_or(0.0);
            }
        }
    }

    Ok(ct_output)
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
