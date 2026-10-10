use crate::core::grid::{TurbineCubatureGrid, TurbineGrid};
/// Wake solver algorithms
///
/// Corresponds to core/solver.py in Python implementation
use crate::core::wake::yaw_added_recovery::{
    calculate_transverse_velocity, wake_added_yaw, yaw_added_turbulence_mixing,
};
use crate::core::wake::WakeModelManager;
use crate::core::{Farm, FlowField, Grid};
use crate::types::{Array1, Array2, Array3, Array4, ArrayD, Float};
use crate::utilities::cosd;
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

        // Use rotor_diameters_sorted with sorted index i
        let rotor_diameter_i = farm.rotor_diameters_sorted[[0, i]];

        // Use hub_heights_sorted with sorted index i
        let hub_height_i = farm.hub_heights_sorted[[0, i]];

        // Current turbine velocity slices (captured before this turbine's wake is added)
        let u_i = flow_field
            .u_sorted
            .slice(s![.., i..i + 1, .., ..])
            .to_owned();
        let v_i = flow_field
            .v_sorted
            .slice(s![.., i..i + 1, .., ..])
            .to_owned();
        let w_i = flow_field
            .w_sorted
            .slice(s![.., i..i + 1, .., ..])
            .to_owned();

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

        // Effective yaw per findex including secondary steering
        // (Python: effective_yaw_i += added_yaw, per findex)
        let mut effective_yaw_arr = yaw_arr.clone();
        if model_manager.enable_secondary_steering {
            let z_i_block = z_grid.slice(s![.., i..i + 1, .., ..]).to_owned();
            let mut delta_y_ss = Array::zeros((n_findex, 1, grid_y_dim, grid_z_dim));
            for fi in 0..n_findex {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        delta_y_ss[[fi, 0, iy, iz]] = y_grid[[fi, i, iy, iz]] - y_2d[[fi, 0]];
                    }
                }
            }
            let ct_4d = ct_arr.clone().into_shape_with_order((n_findex, 1, 1, 1))?;
            let tsr_4d = tsr_arr.clone().into_shape_with_order((n_findex, 1, 1, 1))?;
            let ai_4d = ai_arr.clone().into_shape_with_order((n_findex, 1, 1, 1))?;
            let added_yaw = wake_added_yaw(
                &u_i,
                &v_i,
                &flow_field.u_initial_sorted,
                &delta_y_ss,
                &z_i_block,
                rotor_diameter_i,
                hub_height_i,
                &ct_4d,
                &tsr_4d,
                &ai_4d,
                flow_field.wind_shear,
                1.0,
            );
            for fi in 0..n_findex {
                effective_yaw_arr[fi] += added_yaw[fi];
            }
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

        // Calculate the deflection field of source turbine i at every grid
        // point (Python: per-findex effective yaw/TI/Ct, full-field evaluation)
        let ti_source: Array1 = (0..n_findex)
            .map(|fi| turbine_turbulence_intensity[[fi, i]])
            .collect();
        let deflection_field = model_manager.deflection_model.function(
            &x_2d,
            &y_2d,
            &effective_yaw_arr,
            &ti_source,
            &ct_arr,
            rotor_diameter_i,
            &x_grid,
            &std::collections::HashMap::new(),
        )?;

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
            deflection_field,
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
                let downstream_ok =
                    dx > 0.0 && dx <= downstream_influence_length && dy < 2.0 * rotor_diameter_i;

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
                        flow_field.u_sorted[[fi, ti, iy, iz]] = (u_initial - wake_deficit).max(0.0);
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

    // Store the evolved turbine TI (constant per rotor block) so downstream
    // consumers (e.g. full_flow_sequential_solver) can read the averaged TI.
    flow_field.turbulence_intensity_field_sorted =
        Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    flow_field.turbulence_intensity_field_sorted[[fi, ti, iy, iz]] =
                        turbine_turbulence_intensity[[fi, ti]];
                }
            }
        }
    }

    Ok(())
}

/// TurbOPark solver (port of Python turbopark_solver).
///
/// The turbopark model is target-centric: for each wake source turbine i
/// (upstream to downstream), the model computes the deficit experienced by
/// turbine i from all upstream turbines (including mirrored wakes), and the
/// solver accumulates it via the combination model.
pub fn turbopark_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    if model_manager.enable_secondary_steering {
        anyhow::bail!("Secondary steering not available for this model.");
    }
    if model_manager.enable_transverse_velocities {
        anyhow::bail!("Transverse velocities not used in this model.");
    }
    if model_manager.enable_yaw_added_recovery {
        anyhow::bail!("Yaw added recovery not used in this model.");
    }

    let tp = model_manager
        .velocity_model
        .as_any()
        .downcast_ref::<crate::core::wake_velocity::TurbOParkVelocityDeficit>()
        .ok_or_else(|| anyhow::anyhow!("turbopark_solver requires TurbOParkVelocityDeficit"))?;

    let n_turbines = grid.n_turbines();
    let n_findex = grid.n_findex();
    let shape = grid.x_sorted().shape();
    let grid_second_dim = shape[1];
    let grid_y_dim = shape[2];
    let grid_z_dim = shape[3];

    let mut wake_field: Array4 = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));

    // Per-turbine TI (n_findex, n_turbines), starting at ambient
    let mut turbine_turbulence_intensity: Array2 = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            turbine_turbulence_intensity[[fi, ti]] = flow_field.turbulence_intensities[fi];
        }
    }

    for i in 0..n_turbines {
        // All turbines' thrust coefficients from the current u_sorted (matching
        // Python: recomputed each iteration as the wakes evolve)
        let cts_all = Farm::thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        // Deflection field at turbine i due to each upstream source:
        // (n_findex, n_turbines, gy, gz), zeros for the unyawed case
        let mut deflection_field =
            Array4::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));

        if farm.yaw_angles_sorted.iter().any(|&v| v != 0.0) {
            log::warn!(
                "WARNING: Deflection with the TurbOPark model has not been fully validated."
            );
            for ii in 0..i {
                let x_ii = grid.x_sorted().slice(s![.., ii..ii + 1, .., ..]).to_owned();
                let y_ii = grid.y_sorted().slice(s![.., ii..ii + 1, .., ..]).to_owned();
                let x_2d = mean_value_4d_slice(&x_ii)?;
                let y_2d = mean_value_4d_slice(&y_ii)?;

                // Per-findex yaw/TI/Ct of upstream source ii (Python turbopark_solver)
                let mut yaw_ii_arr = Array1::zeros(n_findex);
                let mut ti_ii_arr = Array1::zeros(n_findex);
                let mut ct_ii_arr = Array1::zeros(n_findex);
                for fi in 0..n_findex {
                    yaw_ii_arr[fi] = farm.yaw_angles_sorted[[fi, ii]];
                    ti_ii_arr[fi] = turbine_turbulence_intensity[[fi, ii]];
                    ct_ii_arr[fi] = cts_all[[fi, ii]];
                }
                let d_ii = farm.rotor_diameters_sorted[[0, ii]];

                let field_ii = model_manager.deflection_model.function(
                    &x_2d,
                    &y_2d,
                    &yaw_ii_arr,
                    &ti_ii_arr,
                    &ct_ii_arr,
                    d_ii,
                    grid.x_sorted(),
                    &std::collections::HashMap::new(),
                )?;

                // Python: deflection_field[:, ii] = deflection_field_ii[:, i]
                // — source ii's deflection evaluated at TARGET i's grid points
                for fi in 0..n_findex {
                    for iy in 0..grid_y_dim {
                        for iz in 0..grid_z_dim {
                            deflection_field[[fi, ii, iy, iz]] = field_ii[[fi, i, iy, iz]];
                        }
                    }
                }
            }
        }

        // Target-centric turbopark deficit at turbine i
        let velocity_deficit = crate::core::wake_velocity::turbopark::turbopark_velocity_deficit(
            &flow_field.turbulence_intensities,
            &cts_all,
            &farm.rotor_diameters.clone(),
            grid.x_sorted(),
            grid.y_sorted(),
            grid.z_sorted(),
            &deflection_field,
            i,
            tp.a,
            tp.sigma_max_rel,
        );

        // Convert relative deficit to absolute, accumulate via the combination model
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
        wake_field = model_manager
            .combination_model
            .function(&wake_field, &velocity_deficit_absolute)?;

        // --- Wake-Added Turbulence (WAT) clamp, matching Python turbopark_solver ---
        let x_mean = mean_value_4d_slice(&grid.x_sorted().clone())?;
        let y_mean = mean_value_4d_slice(&grid.y_sorted().clone())?;

        let ai_all = Farm::axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        let x_i_4d_slice = grid.x_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
        let x_2d = mean_value_4d_slice(&x_i_4d_slice)?;

        let mut x_i_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ai_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ambient_ti_4d = Array4::zeros((n_findex, 1, 1, 1));
        for fi in 0..n_findex {
            x_i_4d[[fi, 0, 0, 0]] = x_2d[[fi, 0]];
            ai_4d[[fi, 0, 0, 0]] = ai_all[[fi, i]];
            ambient_ti_4d[[fi, 0, 0, 0]] = flow_field.turbulence_intensities[fi];
        }

        let rotor_diameter_i = farm.rotor_diameters_sorted[[0, i]];
        let mut wake_added = model_manager.turbulence_model.function(
            &ambient_ti_4d,
            grid.x_sorted(),
            &x_i_4d,
            rotor_diameter_i,
            &ai_4d,
        )?;
        wake_added = wake_added.mapv(|v| if v.is_finite() { v } else { 0.0 });

        let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;
        let iy_c = grid_y_dim / 2;
        let iz_c = grid_z_dim / 2;
        let downstream_influence_length = 15.0 * rotor_diameter_i;

        for fi in 0..n_findex {
            for j in 0..n_turbines {
                let mut overlap_count = 0.0_f64;
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        if velocity_deficit_absolute[[fi, j, iy, iz]] > 0.05 {
                            overlap_count += 1.0;
                        }
                    }
                }
                let area_overlap = overlap_count / grid_res_sq;
                let wake_added_center = wake_added[[fi, j, iy_c, iz_c]];

                let dx = x_mean[[fi, j]] - x_mean[[fi, i]];
                let dy = (y_mean[[fi, j]] - y_mean[[fi, i]]).abs();
                let downstream_ok =
                    dx > 0.0 && dx <= downstream_influence_length && dy < 2.0 * rotor_diameter_i;

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

        // u_sorted = u_initial_sorted - wake_field
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        flow_field.u_sorted[[fi, ti, iy, iz]] = flow_field.u_initial_sorted
                            [[fi, ti, iy, iz]]
                            - wake_field[[fi, ti, iy, iz]];
                    }
                }
            }
        }
    }

    flow_field.turbulence_intensity_field_sorted =
        Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    flow_field.turbulence_intensity_field_sorted[[fi, ti, iy, iz]] =
                        turbine_turbulence_intensity[[fi, ti]];
                }
            }
        }
    }

    Ok(())
}

/// Full-flow sequential solver for visualization grids (planar/cross sections).
///
/// Port of `full_flow_sequential_solver(...)` in `core/solver.py`. Solves the
/// wake on a fine turbine grid first, then evaluates each turbine's wake over
/// the full flow-field grid (e.g. a `FlowFieldPlanarGrid`) using a synthetic
/// grid so the velocity model's source turbine index resolves correctly.
pub fn full_flow_sequential_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    let n_turbines = farm.n_turbines();
    let n_findex = grid.n_findex();
    let grid_shape = grid.x_sorted().shape();
    let n_points = grid_shape[1];
    let grid_y_dim = grid_shape[2];
    let grid_z_dim = grid_shape[3];

    // Build a turbine grid to resolve per-turbine rotor-averaged quantities.
    let turbine_grid = TurbineGrid::new(
        farm.coordinates(),
        farm.rotor_diameters().clone(),
        flow_field.wind_directions.clone(),
        3,
    )?;

    let mut turbine_grid_farm = farm.clone();
    let mut turbine_grid_flow_field = flow_field.clone();
    turbine_grid_farm.expand_farm_properties(n_findex, turbine_grid.sorted_coord_indices());
    turbine_grid_flow_field.initialize_flow_field(
        turbine_grid.grid_shape(),
        turbine_grid.z_sorted(),
        &turbine_grid.hub_heights(),
    );
    turbine_grid_farm.initialize(turbine_grid.sorted_indices());
    sequential_solver(
        &turbine_grid_farm,
        &mut turbine_grid_flow_field,
        &turbine_grid,
        model_manager,
    )?;

    // Source positions on the turbine grid (mean of each turbine's rotor block).
    let source_x = mean_value_4d_slice(&turbine_grid.x_sorted().clone())?; // (n_findex, n_turbines)
    let source_y = mean_value_4d_slice(&turbine_grid.y_sorted().clone())?;
    let source_z = mean_value_4d_slice(&turbine_grid.z_sorted().clone())?;

    // Full-flow grid coordinates.
    let planar_x = grid.x_sorted().clone();
    let planar_y = grid.y_sorted().clone();
    let planar_z = grid.z_sorted().clone();

    // Prepare model-specific velocity parameters once for the full planar grid,
    // matching FLORIS full_flow_sequential_solver.
    let velocity_model_args = model_manager
        .velocity_model
        .prepare_function(grid, flow_field)?;

    // Wake and transverse velocity fields over the full grid.
    let mut wake_field: Array4 = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
    let mut v_wake: Array4 = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
    let mut w_wake: Array4 = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));

    // Turbulence intensity field over the full grid, initialized to ambient.
    let mut turbulence_intensity_field = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
    let mut ambient_turbulence_intensities =
        Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
    for fi in 0..n_findex {
        for p in 0..n_points {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    ambient_turbulence_intensities[[fi, p, iy, iz]] =
                        flow_field.turbulence_intensities[fi];
                    turbulence_intensity_field[[fi, p, iy, iz]] =
                        flow_field.turbulence_intensities[fi];
                }
            }
        }
    }

    for i in 0..n_turbines {
        // Per-turbine quantities from the turbine-grid solve.
        let ct_all = Farm::thrust_coefficient(
            &turbine_grid_flow_field.u_sorted,
            turbine_grid_flow_field.air_density,
            &turbine_grid_farm.turbines,
            Some(&turbine_grid_farm.yaw_angles_sorted),
            Some(&turbine_grid_farm.tilt_angles_sorted),
            turbine_grid.average_method(),
        )?;
        let a_all = Farm::axial_induction(
            &turbine_grid_flow_field.u_sorted,
            turbine_grid_flow_field.air_density,
            &turbine_grid_farm.turbines,
            Some(&turbine_grid_farm.yaw_angles_sorted),
            Some(&turbine_grid_farm.tilt_angles_sorted),
            turbine_grid.average_method(),
        )?;

        let a_i_scalar = a_all[[0, i]];

        // Per-findex scalars for the current turbine.
        let mut yaw_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut tsr_arr = Array1::zeros(n_findex);
        let mut ai_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = turbine_grid_farm.yaw_angles_sorted[[fi, i]];
            ct_arr[fi] = ct_all[[fi, i]];
            tsr_arr[fi] = turbine_grid_farm.tsrs_sorted[[fi, i]];
            ai_arr[fi] = a_all[[fi, i]];
        }

        // Rotor-averaged turbulence intensity for this turbine (constant per block).
        let ti_avg = mean_value_4d_slice(
            &turbine_grid_flow_field
                .turbulence_intensity_field_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned(),
        )?;
        let mut ti_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            ti_arr[fi] = ti_avg[[fi, 0]];
        }

        let rotor_diameter_i = turbine_grid_farm.rotor_diameters_sorted[[0, i]];
        let hub_height_i = turbine_grid_farm.hub_heights_sorted[[0, i]];

        // Effective yaw per findex including secondary steering.
        let mut effective_yaw_arr = yaw_arr.clone();
        if model_manager.enable_secondary_steering {
            let u_i = turbine_grid_flow_field
                .u_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let v_i = turbine_grid_flow_field
                .v_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let source_y_4d = source_y
                .slice(s![.., i..i + 1])
                .insert_axis(ndarray::Axis(2))
                .insert_axis(ndarray::Axis(3));
            let delta_y = &turbine_grid.y_sorted().slice(s![.., i..i + 1, .., ..]) - &source_y_4d;
            let z_i_block = turbine_grid
                .z_sorted()
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let ct_4d = ct_arr.clone().into_shape_with_order((n_findex, 1, 1, 1))?;
            let tsr_4d = tsr_arr.clone().into_shape_with_order((n_findex, 1, 1, 1))?;
            let ai_4d = ai_arr.clone().into_shape_with_order((n_findex, 1, 1, 1))?;
            let added_yaw = wake_added_yaw(
                &u_i,
                &v_i,
                &turbine_grid_flow_field.u_initial_sorted,
                &delta_y,
                &z_i_block,
                rotor_diameter_i,
                hub_height_i,
                &ct_4d,
                &tsr_4d,
                &ai_4d,
                flow_field.wind_shear,
                1.0,
            );
            for fi in 0..n_findex {
                effective_yaw_arr[fi] += added_yaw[fi];
            }
        }

        // Source position for this turbine.
        let x_i = source_x.slice(s![.., i..i + 1]).to_owned(); // (n_findex, 1)
        let y_i = source_y.slice(s![.., i..i + 1]).to_owned();

        // Distance from the source turbine to every full-grid point.
        let mut delta_x = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
        let mut delta_y = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
        for fi in 0..n_findex {
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        delta_x[[fi, p, iy, iz]] = planar_x[[fi, p, iy, iz]] - x_i[[fi, 0]];
                        delta_y[[fi, p, iy, iz]] = planar_y[[fi, p, iy, iz]] - y_i[[fi, 0]];
                    }
                }
            }
        }

        // Transverse velocities over the full grid.
        if model_manager.enable_transverse_velocities {
            let u_i = turbine_grid_flow_field
                .u_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let (vw, ww) = calculate_transverse_velocity(
                &u_i,
                &flow_field.u_initial_sorted,
                &flow_field.dudz_initial_sorted,
                &delta_x,
                &delta_y,
                &planar_z,
                rotor_diameter_i,
                hub_height_i,
                &yaw_arr,
                &ct_arr,
                &tsr_arr,
                &ai_arr,
                flow_field.wind_shear,
                1.0,
            );
            for fi in 0..n_findex {
                for p in 0..n_points {
                    for iy in 0..grid_y_dim {
                        for iz in 0..grid_z_dim {
                            v_wake[[fi, p, iy, iz]] += vw[[fi, p, iy, iz]];
                            w_wake[[fi, p, iy, iz]] += ww[[fi, p, iy, iz]];
                        }
                    }
                }
            }
        }

        // Build the synthetic grid: first n_turbines blocks are the source
        // turbines, the rest are the full-flow field points.
        let total = n_turbines + n_points;
        let mut x_full = Array::zeros((n_findex, total, 1, 1));
        let mut y_full = Array::zeros((n_findex, total, 1, 1));
        let mut z_full = Array::zeros((n_findex, total, 1, 1));
        for fi in 0..n_findex {
            for b in 0..n_turbines {
                x_full[[fi, b, 0, 0]] = source_x[[fi, b]];
                y_full[[fi, b, 0, 0]] = source_y[[fi, b]];
                z_full[[fi, b, 0, 0]] = source_z[[fi, b]];
            }
            for p in 0..n_points {
                x_full[[fi, n_turbines + p, 0, 0]] = planar_x[[fi, p, 0, 0]];
                y_full[[fi, n_turbines + p, 0, 0]] = planar_y[[fi, p, 0, 0]];
                z_full[[fi, n_turbines + p, 0, 0]] = planar_z[[fi, p, 0, 0]];
            }
        }

        // Deflection field of source i evaluated at every point of the
        // synthetic grid (Python: deflection on the full flow grid).
        let deflection_field = model_manager.deflection_model.function(
            &x_i,
            &y_i,
            &effective_yaw_arr,
            &ti_arr,
            &ct_arr,
            rotor_diameter_i,
            &x_full,
            &std::collections::HashMap::new(),
        )?;

        // Per-turbine broadcast arrays (only column `i` is read by the model).
        let mut a_arr = Array::zeros((n_findex, total));
        let mut yaw_bc = Array::zeros((n_findex, total));
        let mut ti_bc = Array::zeros((n_findex, total));
        let mut ct_bc = Array::zeros((n_findex, total));
        for fi in 0..n_findex {
            a_arr[[fi, i]] = a_i_scalar;
            yaw_bc[[fi, i]] = yaw_arr[fi];
            ti_bc[[fi, i]] = ti_arr[fi];
            ct_bc[[fi, i]] = ct_arr[fi];
        }

        let velocity_deficit = model_manager.velocity_model.function(
            x_full,
            y_full,
            z_full,
            a_arr,
            deflection_field,
            yaw_bc,
            ti_bc,
            ct_bc,
            hub_height_i,
            rotor_diameter_i,
            i,
            &velocity_model_args,
        )?;

        // Relative deficit over the field points (drop the source blocks).
        let mut velocity_deficit_absolute =
            Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));
        for fi in 0..n_findex {
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        velocity_deficit_absolute[[fi, p, iy, iz]] = flow_field.u_initial_sorted
                            [[fi, p, iy, iz]]
                            * velocity_deficit[[fi, n_turbines + p, iy, iz]];
                    }
                }
            }
        }

        wake_field = model_manager
            .combination_model
            .function(&wake_field, &velocity_deficit_absolute)?;

        // Wake-added turbulence over the full grid.
        let mut x_i_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ai_4d = Array4::zeros((n_findex, 1, 1, 1));
        for fi in 0..n_findex {
            x_i_4d[[fi, 0, 0, 0]] = x_i[[fi, 0]];
            ai_4d[[fi, 0, 0, 0]] = a_i_scalar;
        }
        let mut wake_added_ti = model_manager.turbulence_model.function(
            &ambient_turbulence_intensities,
            &planar_x,
            &x_i_4d,
            rotor_diameter_i,
            &ai_4d,
        )?;
        wake_added_ti = wake_added_ti.mapv(|v| if v.is_finite() { v } else { 0.0 });

        let downstream_influence_length = 15.0 * rotor_diameter_i;
        for fi in 0..n_findex {
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let area_overlap = if velocity_deficit_absolute[[fi, p, iy, iz]] > 0.05 {
                            1.0
                        } else {
                            0.0
                        };
                        let downstream = planar_x[[fi, p, iy, iz]] > x_i[[fi, 0]]
                            && (y_i[[fi, 0]] - planar_y[[fi, p, iy, iz]]).abs()
                                < 2.0 * rotor_diameter_i
                            && planar_x[[fi, p, iy, iz]]
                                <= downstream_influence_length + x_i[[fi, 0]];
                        let ti_added = if downstream {
                            area_overlap * wake_added_ti[[fi, p, iy, iz]]
                        } else {
                            0.0
                        };
                        let ambient = ambient_turbulence_intensities[[fi, p, iy, iz]];
                        let new_ti = (ti_added * ti_added + ambient * ambient).sqrt();
                        if new_ti > turbulence_intensity_field[[fi, p, iy, iz]] {
                            turbulence_intensity_field[[fi, p, iy, iz]] = new_ti;
                        }
                    }
                }
            }
        }

        // Update the full-flow sorted fields.
        for fi in 0..n_findex {
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        flow_field.u_sorted[[fi, p, iy, iz]] = (flow_field.u_initial_sorted
                            [[fi, p, iy, iz]]
                            - wake_field[[fi, p, iy, iz]])
                        .max(0.0);
                    }
                }
            }
        }
        if model_manager.enable_transverse_velocities {
            for fi in 0..n_findex {
                for p in 0..n_points {
                    for iy in 0..grid_y_dim {
                        for iz in 0..grid_z_dim {
                            flow_field.v_sorted[[fi, p, iy, iz]] += v_wake[[fi, p, iy, iz]];
                            flow_field.w_sorted[[fi, p, iy, iz]] += w_wake[[fi, p, iy, iz]];
                        }
                    }
                }
            }
        }
    }

    flow_field.turbulence_intensity_field_sorted = turbulence_intensity_field;

    Ok(())
}

/// Cumulative Curl solver (port of Python cc_solver).
///
/// The CC model maintains cumulative state (turb_u_wake, Ctmp): each turbine's
/// rotor inflow is updated to u_initial - turb_u_wake, Ct/AI are computed from
/// those inflow velocities, and the model accumulates the wake velocity.
/// The final u_sorted is the inflow field.
pub fn cc_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    let cc = model_manager
        .velocity_model
        .as_any()
        .downcast_ref::<crate::core::wake_velocity::CumulativeCurlVelocityDeficit>()
        .ok_or_else(|| anyhow::anyhow!("cc_solver requires CumulativeCurlVelocityDeficit"))?;

    let n_turbines = grid.n_turbines();
    let n_findex = grid.n_findex();
    let shape = grid.x_sorted().shape();
    let grid_y_dim = shape[2];
    let grid_z_dim = shape[3];

    let mut turb_u_wake: Array4 = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
    let mut ctmp: ArrayD = ArrayD::from_elem(
        ndarray::IxDyn(&[n_turbines, n_findex, n_turbines, grid_y_dim, grid_z_dim]),
        0.0,
    );
    let mut turb_inflow_field = flow_field.u_initial_sorted.clone();

    let mut turbine_turbulence_intensity: Array2 = Array::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            turbine_turbulence_intensity[[fi, ti]] = flow_field.turbulence_intensities[fi];
        }
    }

    let x_mean_all = mean_value_4d_slice(&grid.x_sorted().clone())?;
    let y_mean_all = mean_value_4d_slice(&grid.y_sorted().clone())?;

    for i in 0..n_turbines {
        let d_i = farm.rotor_diameters_sorted[[0, i]];

        let x_i_slice = grid.x_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
        let y_i_slice = grid.y_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
        let x_i_2d = mean_value_4d_slice(&x_i_slice)?;
        let y_i_2d = mean_value_4d_slice(&y_i_slice)?;

        // The inflow at turbine i's rotor = u_initial - turb_u_wake
        for fi in 0..n_findex {
            for t in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let at_rotor = (grid.x_sorted()[[fi, t, iy, iz]] - x_i_2d[[fi, 0]]).abs()
                            < 0.01
                            && (grid.y_sorted()[[fi, t, iy, iz]] - y_i_2d[[fi, 0]]).abs()
                                < 0.51 * d_i;
                        if at_rotor {
                            turb_inflow_field[[fi, t, iy, iz]] = flow_field.u_initial_sorted
                                [[fi, t, iy, iz]]
                                - turb_u_wake[[fi, t, iy, iz]];
                        }
                    }
                }
            }
        }

        // Per-turbine average inflow velocity (cubic-mean, matching Python average_velocity)
        let mut turb_avg_vels = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for t in 0..n_turbines {
                let mut v_cubed_sum = 0.0;
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        v_cubed_sum += turb_inflow_field[[fi, t, iy, iz]].powi(3);
                    }
                }
                turb_avg_vels[[fi, t]] =
                    (v_cubed_sum / (grid_y_dim * grid_z_dim) as f64).powf(1.0 / 3.0);
            }
        }

        // All turbines' Ct/AI from the inflow velocities (1x1 rotor grid)
        let mut avg_vels_4d = Array::zeros((n_findex, n_turbines, 1, 1));
        for fi in 0..n_findex {
            for t in 0..n_turbines {
                avg_vels_4d[[fi, t, 0, 0]] = turb_avg_vels[[fi, t]];
            }
        }
        let turb_cts_2d = Farm::thrust_coefficient(
            &avg_vels_4d,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;
        let turb_ais_2d = Farm::axial_induction(
            &avg_vels_4d,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        // The source's inflow (n_findex, 1, gy, gz) and yaw per findex
        let u_i = turb_inflow_field.slice(s![.., i..i + 1, .., ..]).to_owned();
        let mut yaw_i_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_i_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
        }

        // Per-findex scalar arrays for the source turbine (ct, tsr, ai, yaw)
        let hub_height_i = farm.hub_heights_sorted[[0, i]];
        let tsr_i = farm.tsrs_sorted[[0, i]];
        let mut yaw_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut tsr_arr = Array1::zeros(n_findex);
        let mut ai_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
            ct_arr[fi] = turb_cts_2d[[fi, i]];
            tsr_arr[fi] = tsr_i;
            ai_arr[fi] = turb_ais_2d[[fi, i]];
        }

        // GCH: transverse velocities (v/w) for the source turbine (Python cc: scale=2.0)
        let mut v_wake = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
        let mut w_wake = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
        if model_manager.enable_transverse_velocities {
            let mut delta_x = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
            let mut delta_y = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    for iy in 0..grid_y_dim {
                        for iz in 0..grid_z_dim {
                            delta_x[[fi, ti, iy, iz]] =
                                grid.x_sorted()[[fi, ti, iy, iz]] - x_i_2d[[fi, 0]];
                            delta_y[[fi, ti, iy, iz]] =
                                grid.y_sorted()[[fi, ti, iy, iz]] - y_i_2d[[fi, 0]];
                        }
                    }
                }
            }
            let (vw, ww) = calculate_transverse_velocity(
                &u_i,
                &flow_field.u_initial_sorted,
                &flow_field.dudz_initial_sorted,
                &delta_x,
                &delta_y,
                &grid.z_sorted().clone(),
                d_i,
                hub_height_i,
                &yaw_arr,
                &ct_arr,
                &tsr_arr,
                &ai_arr,
                flow_field.wind_shear,
                2.0,
            );
            v_wake = vw;
            w_wake = ww;
        }

        // GCH: yaw-added recovery (Python cc: gch_gain=1.0)
        if model_manager.enable_yaw_added_recovery {
            let v_i = flow_field
                .v_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let w_i = flow_field
                .w_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let i_mixing = yaw_added_turbulence_mixing(
                &u_i,
                &turbine_turbulence_intensity.slice(s![.., i]).to_owned(),
                &v_i,
                &w_i,
                &v_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
                &w_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
            );
            let gch_gain = 1.0;
            for fi in 0..n_findex {
                turbine_turbulence_intensity[[fi, i]] += gch_gain * i_mixing[fi];
            }
        }

        // Effective yaw per findex including secondary steering (Python cc: scale=2.0)
        let mut effective_yaw_arr = yaw_i_arr.clone();
        if model_manager.enable_secondary_steering {
            let v_i = flow_field
                .v_sorted
                .slice(s![.., i..i + 1, .., ..])
                .to_owned();
            let y_i_block = grid.y_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
            let z_i_block = grid.z_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
            let hub_height_i = farm.hub_heights_sorted[[0, i]];
            let tsr_i = farm.tsrs_sorted[[0, i]];

            let mut ct_4d = Array4::zeros((n_findex, 1, 1, 1));
            let mut tsr_4d = Array4::zeros((n_findex, 1, 1, 1));
            let mut ai_4d = Array4::zeros((n_findex, 1, 1, 1));
            for fi in 0..n_findex {
                ct_4d[[fi, 0, 0, 0]] = turb_cts_2d[[fi, i]];
                tsr_4d[[fi, 0, 0, 0]] = tsr_i;
                ai_4d[[fi, 0, 0, 0]] = turb_ais_2d[[fi, i]];
            }

            let mut y_i_4d = Array4::zeros((n_findex, 1, grid_y_dim, grid_z_dim));
            for fi in 0..n_findex {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        y_i_4d[[fi, 0, iy, iz]] = y_i_2d[[fi, 0]];
                    }
                }
            }
            let delta_y = &y_i_block - &y_i_4d;

            let added_yaw = wake_added_yaw(
                &u_i,
                &v_i,
                &flow_field.u_initial_sorted,
                &delta_y,
                &z_i_block,
                d_i,
                hub_height_i,
                &ct_4d,
                &tsr_4d,
                &ai_4d,
                flow_field.wind_shear,
                2.0,
            );
            for fi in 0..n_findex {
                effective_yaw_arr[fi] += added_yaw[fi];
            }
        }

        // The deflection field of source i at every grid point
        // (Python cc: per-findex effective yaw/TI/Ct, full-field evaluation)
        let ti_source: Array1 = (0..n_findex)
            .map(|fi| turbine_turbulence_intensity[[fi, i]])
            .collect();
        let ct_source: Array1 = (0..n_findex).map(|fi| turb_cts_2d[[fi, i]]).collect();
        let deflection_field = model_manager.deflection_model.function(
            &x_i_2d,
            &y_i_2d,
            &effective_yaw_arr,
            &ti_source,
            &ct_source,
            d_i,
            grid.x_sorted(),
            &std::collections::HashMap::new(),
        )?;

        // The cumulative CC model call
        let (new_turb_u_wake, new_ctmp) =
            crate::core::wake_velocity::cumulative_gauss_curl::cc_velocity_deficit(
                cc,
                i,
                grid.x_sorted(),
                grid.y_sorted(),
                grid.z_sorted(),
                &u_i,
                &flow_field.u_initial_sorted,
                &deflection_field,
                &yaw_i_arr,
                &turbine_turbulence_intensity,
                &turb_cts_2d,
                &farm.rotor_diameters.clone(),
                turb_u_wake,
                ctmp,
            );
        turb_u_wake = new_turb_u_wake;
        ctmp = new_ctmp;

        // --- WAT (matching Python cc_solver: the inverted area_overlap) ---
        let mut x_i_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ai_4d = Array4::zeros((n_findex, 1, 1, 1));
        let mut ambient_ti_4d = Array4::zeros((n_findex, 1, 1, 1));
        for fi in 0..n_findex {
            x_i_4d[[fi, 0, 0, 0]] = x_i_2d[[fi, 0]];
            ai_4d[[fi, 0, 0, 0]] = turb_ais_2d[[fi, i]];
            ambient_ti_4d[[fi, 0, 0, 0]] = flow_field.turbulence_intensities[fi];
        }

        let mut wake_added = model_manager.turbulence_model.function(
            &ambient_ti_4d,
            grid.x_sorted(),
            &x_i_4d,
            d_i,
            &ai_4d,
        )?;
        wake_added = wake_added.mapv(|v| if v.is_finite() { v } else { 0.0 });

        let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;
        let iy_c = grid_y_dim / 2;
        let iz_c = grid_z_dim / 2;
        let downstream_influence_length = 15.0 * d_i;

        for fi in 0..n_findex {
            for j in 0..n_turbines {
                // area_overlap = 1 - fraction of turbine j's rotor with wake <= 0.05
                let mut unwaked_count = 0.0_f64;
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        if turb_u_wake[[fi, j, iy, iz]] <= 0.05 {
                            unwaked_count += 1.0;
                        }
                    }
                }
                let area_overlap = 1.0 - unwaked_count / grid_res_sq;
                let wake_added_center = wake_added[[fi, j, iy_c, iz_c]];

                let dx = x_mean_all[[fi, j]] - x_i_2d[[fi, 0]];
                let dy = (y_i_2d[[fi, 0]] - y_mean_all[[fi, j]]).abs();
                let downstream_ok = dx > 0.0 && dx <= downstream_influence_length && dy < 2.0 * d_i;

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

        // Add transverse velocities to the flow field (Python cc: v_sorted += v_wake)
        flow_field.v_sorted += &v_wake;
        flow_field.w_sorted += &w_wake;
    }

    // The final u_sorted is the inflow field (matching Python)
    flow_field.u_sorted = turb_inflow_field;

    // Store the TI field (constant per rotor block)
    flow_field.turbulence_intensity_field_sorted =
        Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    flow_field.turbulence_intensity_field_sorted[[fi, ti, iy, iz]] =
                        turbine_turbulence_intensity[[fi, ti]];
                }
            }
        }
    }

    Ok(())
}

/// Empirical Gauss solver (port of Python empirical_gauss_solver).
///
/// For each wake source turbine (upstream to downstream), computes the deficit
/// it adds to every downstream turbine and accumulates it via the combination
/// model. The wake-induced mixing factor is built up iteratively: each source
/// contributes `area_overlap * a_i / downstream_D^2` to downstream turbines'
/// mixing terms, and the RSS over upstream contributions feeds the velocity
/// model's wake width via `mixing_gain_velocity`. Returns nothing; the mixing
/// field is an internal quantity of this solver (Python returns it for the
/// full-flow solver, which is not ported).
pub fn empirical_gauss_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    if model_manager.enable_secondary_steering {
        anyhow::bail!("Secondary steering not available for this model.");
    }
    if model_manager.enable_transverse_velocities {
        anyhow::bail!("Transverse velocities not used in this model.");
    }

    // The mixing architecture is defined for turbine grids (Python's full-flow
    // solver runs the turbine-grid solver on a copy first).
    if grid.as_any().downcast_ref::<TurbineGrid>().is_none() {
        return sequential_solver(farm, flow_field, grid, model_manager);
    }

    let eg = model_manager
        .velocity_model
        .as_any()
        .downcast_ref::<crate::core::wake_velocity::EmpiricalGaussVelocityDeficit>()
        .ok_or_else(|| {
            anyhow::anyhow!("empirical_gauss_solver requires EmpiricalGaussVelocityDeficit")
        })?;
    let defl_eg = model_manager
        .deflection_model
        .as_any()
        .downcast_ref::<crate::core::wake_deflection::EmpiricalGaussVelocityDeflection>()
        .ok_or_else(|| {
            anyhow::anyhow!("empirical_gauss_solver requires EmpiricalGaussVelocityDeflection")
        })?;

    let n_turbines = grid.n_turbines();
    let n_findex = grid.n_findex();
    let shape = grid.x_sorted().shape();
    let grid_y_dim = shape[2];
    let grid_z_dim = shape[3];

    let mut wake_field: Array4 = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));

    // Mean grid positions per turbine (n_findex, n_turbines)
    let x_mean_all = mean_value_4d_slice(&grid.x_sorted().clone())?;

    // Downstream distance in rotor diameters:
    // downstream_distance_D[fi, j, k] = max((x_j - x_k) / D_j, 0.1)
    let mut downstream_distance_d: Array3 = Array::zeros((n_findex, n_turbines, n_turbines));
    for fi in 0..n_findex {
        for j in 0..n_turbines {
            let d_j = farm.rotor_diameters_sorted[[fi, j]];
            for k in 0..n_turbines {
                downstream_distance_d[[fi, j, k]] =
                    ((x_mean_all[[fi, j]] - x_mean_all[[fi, k]]) / d_j).max(0.1);
            }
        }
    }

    // Initialize the mixing factor with the atmospheric TI gain on the diagonal
    if model_manager.atmospheric_ti_gain != 0.0 {
        log::warn!(
            "Running wake_induced_mixing with mixing contributions from the atmospheric \
             turbulence intensity has not been vetted. Set atmospheric_ti_gain=0."
        );
    }
    let mut mixing_factor: Array3 = Array::zeros((n_findex, n_turbines, n_turbines));
    for fi in 0..n_findex {
        for t in 0..n_turbines {
            mixing_factor[[fi, t, t]] =
                model_manager.atmospheric_ti_gain * flow_field.turbulence_intensities[fi];
        }
    }

    let wind_veer = flow_field.wind_veer;
    let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;

    for i in 0..n_turbines {
        // Ct / AI from the current u_sorted (recomputed as wakes evolve,
        // matching Python thrust_coefficient/axial_induction with ix_filter=[i])
        let ct_all = Farm::thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;
        let a_all = Farm::axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            &farm.turbines,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            grid.average_method(),
        )?;

        // Per-findex source-turbine quantities. Python computes tilt via
        // calculate_tilt_for_eff_velocities, which for non-floating turbines
        // (no tilt interpolation) returns tilt_angles_sorted unchanged.
        let mut yaw_arr = Array1::zeros(n_findex);
        let mut tilt_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut a_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
            tilt_arr[fi] = farm.tilt_angles_sorted[[fi, i]];
            ct_arr[fi] = ct_all[[fi, i]];
            a_arr[fi] = a_all[[fi, i]];
        }

        // Total wake-induced mixing for turbine i: RSS over all contributions
        let mut mixing_i = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            let mut sum_sq = 0.0;
            for k in 0..n_turbines {
                sum_sq += mixing_factor[[fi, i, k]] * mixing_factor[[fi, i, k]];
            }
            mixing_i[fi] = sum_sq.sqrt();
        }

        // Deflection field due to turbine i at every grid point
        // (Python: x_normalized = (x - x_i) evaluated on the full grid)
        let d_i = farm.rotor_diameters_sorted[[0, i]];
        let n_pts = n_turbines * grid_y_dim * grid_z_dim;
        let mut x_downstream = Array2::zeros((n_findex, n_pts));
        for fi in 0..n_findex {
            let mut idx = 0;
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        x_downstream[[fi, idx]] =
                            grid.x_sorted()[[fi, ti, iy, iz]] - x_mean_all[[fi, i]];
                        idx += 1;
                    }
                }
            }
        }
        let (dy_flat, dz_flat) =
            defl_eg.deflection_yz(&x_downstream, &yaw_arr, &tilt_arr, &mixing_i, &ct_arr, d_i);
        let deflection_y =
            dy_flat.into_shape_with_order((n_findex, n_turbines, grid_y_dim, grid_z_dim))?;
        let deflection_z =
            dz_flat.into_shape_with_order((n_findex, n_turbines, grid_y_dim, grid_z_dim))?;

        // Velocity deficit at all grid points due to turbine i (with mixing)
        let velocity_deficit = eg.velocity_deficit_emgauss(
            grid.x_sorted(),
            grid.y_sorted(),
            grid.z_sorted(),
            &deflection_y,
            &deflection_z,
            &yaw_arr,
            &tilt_arr,
            &ct_arr,
            &mixing_i,
            farm.hub_heights_sorted[[0, i]],
            d_i,
            i,
            wind_veer,
        );

        // Convert relative deficit to absolute velocity deficit
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

        // wake_field = combination_model.function(wake_field, deficit * u_initial)
        wake_field = model_manager
            .combination_model
            .function(&wake_field, &velocity_deficit_absolute)?;

        // Wake overlap: fraction of each receiving turbine's grid points with
        // an absolute deficit > 0.05 m/s
        let mut area_overlap: Array2 = Array::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for j in 0..n_turbines {
                let mut overlap_count = 0.0;
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        if velocity_deficit_absolute[[fi, j, iy, iz]] > 0.05 {
                            overlap_count += 1.0;
                        }
                    }
                }
                area_overlap[[fi, j]] = overlap_count / grid_res_sq;
            }
        }

        // WakeInducedMixing: turbine i contributes a_i / downstream_D^2 to
        // every turbine's mixing term
        for fi in 0..n_findex {
            for j in 0..n_turbines {
                let dist_sq = downstream_distance_d[[fi, j, i]] * downstream_distance_d[[fi, j, i]];
                mixing_factor[[fi, j, i]] += area_overlap[[fi, j]] * a_arr[fi] / dist_sq;
            }
        }

        // Influence of yawing on the wake (yaw_added_wake_mixing)
        if model_manager.enable_yaw_added_recovery {
            for fi in 0..n_findex {
                for j in 0..n_turbines {
                    let dist_sq =
                        downstream_distance_d[[fi, j, i]] * downstream_distance_d[[fi, j, i]];
                    let yaw_term =
                        a_arr[fi] * defl_eg.yaw_added_mixing_gain * (1.0 - cosd(yaw_arr[fi]))
                            / dist_sq;
                    mixing_factor[[fi, j, i]] += area_overlap[[fi, j]] * yaw_term;
                }
            }
        }

        // Influence of awc on the turbine's own wake (awc_added_wake_mixing)
        if model_manager.enable_active_wake_mixing {
            for fi in 0..n_findex {
                if farm.awc_modes_sorted[[fi, i]] == "helix" {
                    let amp = farm.awc_amplitudes_sorted[[fi, i]];
                    mixing_factor[[fi, i, i]] +=
                        amp.powf(eg.awc_wake_exp) / eg.awc_wake_denominator;
                }
            }
        }

        // Update u_sorted after each turbine's wake is added
        // (Python: flow_field.u_sorted = flow_field.u_initial_sorted - wake_field)
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let wake_deficit = wake_field[[fi, ti, iy, iz]];
                        let u_initial = flow_field.u_initial_sorted[[fi, ti, iy, iz]];
                        flow_field.u_sorted[[fi, ti, iy, iz]] = (u_initial - wake_deficit).max(0.0);
                    }
                }
            }
        }
    }

    // Store the TI field (constant per rotor block)
    flow_field.turbulence_intensity_field_sorted =
        Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    flow_field.turbulence_intensity_field_sorted[[fi, ti, iy, iz]] =
                        flow_field.turbulence_intensities[fi];
                }
            }
        }
    }

    Ok(())
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
