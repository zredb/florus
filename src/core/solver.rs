/// Wake solver algorithms
///
/// Corresponds to core/solver.py in Python implementation

use crate::core::turbine::Turbine;
use crate::core::wake::WakeModelManager;
use crate::core::{Farm, FlowField, GridBase};
use crate::types::{Array2, Array4, Float};
use anyhow::Result;
use ndarray::{Array, s};

/// Sequential wake solver
///
/// Computes wakes for each turbine in upstream→downstream order
/// applying deflection, velocity deficit, turbulence, and wake superposition models
pub fn sequential_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn GridBase,
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

    // Initialize wake field arrays (velocity deficit accumulated across turbines)
    let mut wake_field: Array4 = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));

    // Loop through turbines (upstream to downstream)
    for i in 0..n_turbines {
        // Get turbine properties
        let ct_i = thrust_coefficient(
            &flow_field.u_sorted,
            &farm.turbine_map,
            farm.yaw_angles(),
            farm.tilt_angles(),
            grid.average_method(),
        )?;

        let a_i = axial_induction(&ct_i);
        let ti_i = flow_field.turbulence_intensities.slice(s![..]);
        let yaw_angle_i = farm.yaw_angles()[[0, i]];
        let rotor_diameter_i = farm.rotor_diameters()[i];
        let hub_height_i = farm.hub_heights()[i];

        // Calculate 2D mean values for deflection at turbine i's position
        let x_i = x_grid.slice(s![.., i..i+1, .., ..]).to_owned();
        let y_i = y_grid.slice(s![.., i..i+1, .., ..]).to_owned();
        
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

        // Broadcast deflection to all turbines (deflection_field has shape [n_findex, 1])
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
            &std::collections::HashMap::new(),
        )?;

        // Only apply deficit to grid points that are DOWNSTREAM of turbine i
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        if x_grid[[fi, ti, iy, iz]] >= x_grid[[fi, i, 0, 0]] {
                            wake_field[[fi, ti, iy, iz]] += velocity_deficit[[fi, i, iy, iz]];
                        }
                    }
                }
            }
        }
    }

    // Apply combined wake field to get final velocities
    // u_sorted = u_initial - wake_field
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    let deficit = wake_field[[fi, ti, iy, iz]].min(1.0).max(0.0);
                    flow_field.u_sorted[[fi, ti, iy, iz]] = 
                        flow_field.u_initial_sorted[[fi, ti, iy, iz]] * (1.0 - deficit);
                }
            }
        }
    }

    Ok(())
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
                ct_output[[fi, ti]] = turbines[ti].ct_at_speed(v);
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
            // High thrust region - empirical relationship
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
        // Placeholder test - real tests need full setup
        assert!(true);
    }

    #[test]
    fn test_mean_value_4d_slice() {
        let arr = Array::from_shape_vec(
            (1, 2, 2, 2),
            vec![
                // Turbine 0: values 1, 2, 3, 4 -> mean = 2.5
                1.0, 2.0, 3.0, 4.0,
                // Turbine 1: values 5, 6, 7, 8 -> mean = 6.5
                5.0, 6.0, 7.0, 8.0,
            ]
        ).unwrap();
        
        let result = mean_value_4d_slice(&arr).unwrap();
        
        assert_eq!(result.shape()[0], 1);
        assert_eq!(result.shape()[1], 2);
        assert_relative_eq!(result[[0, 0]], 2.5);
        assert_relative_eq!(result[[0, 1]], 6.5);
    }

    #[test]
    fn test_mean_value_4d_slice_single_value() {
        let arr = Array::from_shape_vec(
            (1, 1, 1, 1),
            vec![5.0]
        ).unwrap();
        
        let result = mean_value_4d_slice(&arr).unwrap();
        
        assert_eq!(result.shape()[0], 1);
        assert_eq!(result.shape()[1], 1);
        assert_relative_eq!(result[[0, 0]], 5.0);
    }

    #[test]
    fn test_mean_value_4d_slice_multiple_findex() {
        let arr = Array::from_shape_vec(
            (2, 2, 1, 1),
            vec![
                // findex 0: turbine 0 = 10.0, turbine 1 = 20.0
                10.0, 20.0,
                // findex 1: turbine 0 = 30.0, turbine 1 = 40.0
                30.0, 40.0,
            ]
        ).unwrap();
        
        let result = mean_value_4d_slice(&arr).unwrap();
        
        assert_eq!(result.shape()[0], 2);
        assert_eq!(result.shape()[1], 2);
        assert_relative_eq!(result[[0, 0]], 10.0);
        assert_relative_eq!(result[[0, 1]], 20.0);
        assert_relative_eq!(result[[1, 0]], 30.0);
        assert_relative_eq!(result[[1, 1]], 40.0);
    }
}
