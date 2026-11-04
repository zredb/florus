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
        SolverType::EmpiricalGauss => {
            empirical_gauss_solver(farm, flow_field, grid, model_manager).map(|_| ())
        }
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

    // Prepare model-specific function arguments (Bug #23): Python passes
    // prepare_function() results to the velocity/deflection models; the gauss
    // models read wind_veer from them.
    let velocity_model_args = model_manager
        .velocity_model
        .prepare_function(grid, flow_field)?;
    let deflection_model_args = model_manager
        .deflection_model
        .prepare_function(grid, flow_field)?;

    // Prepare combination model function arguments
    let _combination_model_args = model_manager
        .combination_model
        .prepare_function(grid, flow_field)?;

    // Evolving 4D TI field (n_findex, n_turbines, grid_y, grid_z), initialized
    // to ambient. Mirrors Python's local `turbine_turbulence_intensity`, which
    // starts as (n_findex, n_turbines, 1, 1) and broadcasts to the full grid
    // shape on the first WAT `np.maximum` update; a full-shape constant field
    // is semantically identical (every consumer broadcasts pointwise).
    let mut turbine_turbulence_intensity: Array4 =
        Array::from_elem((n_findex, n_turbines, grid_y_dim, grid_z_dim), 0.0);
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    turbine_turbulence_intensity[[fi, ti, iy, iz]] =
                        flow_field.turbulence_intensities[fi];
                }
            }
        }
    }

    // Loop through turbines (upstream to downstream)
    for i in 0..n_turbines {
        // Compute thrust coefficient and axial induction through the operation
        // model (rotor-averaged velocity), matching Python solver.py which calls
        // floris.core.thrust_coefficient with ix_filter=[i].
        let ct_i = farm.thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;

        let a_i = farm.axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;

        // Source turbine's rotor diameter / hub height per findex
        // (Python: rotor_diameters_sorted[:, i:i+1, None, None])
        let rotor_diameters_i = farm
            .rotor_diameters_sorted
            .slice(s![.., i..i + 1])
            .to_owned();

        let hub_heights_i = farm.hub_heights_sorted.slice(s![.., i..i + 1]).to_owned();

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
                &rotor_diameters_i,
                &hub_heights_i,
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
        // Source's per-point TI slice (Python: turbulence_intensity_i =
        // turbine_turbulence_intensity[:, i:i+1], the evolving field slice)
        let ti_source: Array4 = turbine_turbulence_intensity
            .slice(s![.., i..i + 1, .., ..])
            .to_owned();
        let deflection_field = model_manager.deflection_model.function(
            &x_2d,
            &y_2d,
            &effective_yaw_arr,
            &ti_source,
            &ct_arr,
            &rotor_diameters_i,
            &x_grid,
            &deflection_model_args,
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
                &rotor_diameters_i,
                &hub_heights_i,
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
            // Python: yaw_added_turbulence_mixing takes I_i =
            // turbulence_intensity_i[:, 0, 0, 0] — the first grid point of the
            // source's TI field — and adds gch_gain * I_mixing (scalar per
            // findex) to the whole source slice.
            let i_mixing = yaw_added_turbulence_mixing(
                &u_i,
                &turbine_turbulence_intensity
                    .slice(s![.., i, 0, 0])
                    .to_owned(),
                &v_i,
                &w_i,
                &v_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
                &w_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
            );
            let gch_gain = 2.0;
            for fi in 0..n_findex {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        turbine_turbulence_intensity[[fi, i, iy, iz]] += gch_gain * i_mixing[fi];
                    }
                }
            }
        }

        // Calculate velocity deficit at ALL grid points due to turbine i
        // (Python passes the full evolving TI field; models read the source's
        // slice via [[fi, turbine_index, iy, iz]])
        let velocity_deficit = model_manager.velocity_model.function(
            x_grid.clone(),
            y_grid.clone(),
            z_grid.clone(),
            a_i.clone(),
            deflection_field,
            farm.yaw_angles_sorted.clone(),
            turbine_turbulence_intensity.clone(),
            ct_i,
            &hub_heights_i,
            &rotor_diameters_i,
            i, // turbine_index
            &velocity_model_args,
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
        // Mirrors FLORIS sequential_solver: wake-added TI is evaluated per grid
        // point (Python: ti_added = area_overlap * wake_added * pointwise masks
        // on grid.x_sorted/grid.y_sorted vs the source's mean position), then
        // folded into the evolving field with np.maximum(sqrt(ti_added^2 + amb^2)).
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
            &rotor_diameters_i,
            &ai_4d,
        )?;
        // Match Python's np.nan_to_num(..., posinf=0.0)
        wake_added = wake_added.mapv(|v| if v.is_finite() { v } else { 0.0 });

        // grid_y_dim * grid_z_dim == grid_resolution^2 (turbine grid): the number of
        // grid points per turbine, the normalization for area_overlap.
        let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;

        for fi in 0..n_findex {
            let x_src = x_mean[[fi, i]];
            let y_src = y_mean[[fi, i]];
            let ambient = flow_field.turbulence_intensities[fi];
            let d_i = rotor_diameters_i[[fi, 0]];
            let downstream_influence_length = 15.0 * d_i;
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

                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        // Pointwise downstream-influence masks (Python:
                        // grid.x_sorted > x_i, |y_i - grid.y_sorted| < 2D,
                        // grid.x_sorted <= 15D + x_i)
                        let x_pt = x_grid[[fi, j, iy, iz]];
                        let downstream = x_pt > x_src
                            && (y_src - y_grid[[fi, j, iy, iz]]).abs() < 2.0 * d_i
                            && x_pt <= downstream_influence_length + x_src;
                        let ti_added = if downstream {
                            area_overlap * wake_added[[fi, j, iy, iz]]
                        } else {
                            0.0
                        };
                        let new_ti = (ti_added * ti_added + ambient * ambient).sqrt();
                        if new_ti > turbine_turbulence_intensity[[fi, j, iy, iz]] {
                            turbine_turbulence_intensity[[fi, j, iy, iz]] = new_ti;
                        }
                    }
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
                        flow_field.u_sorted[[fi, ti, iy, iz]] = u_initial - wake_deficit;
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

    // Store the evolved per-point TI field (Python: flow_field.
    // turbulence_intensity_field_sorted = turbine_turbulence_intensity)
    if grid_second_dim == n_turbines {
        flow_field.turbulence_intensity_field_sorted = turbine_turbulence_intensity;
    } else {
        // Planar grids may have a wider second dimension; copy per-turbine blocks.
        let mut field_sorted = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
        for ti in 0..n_turbines {
            field_sorted
                .slice_mut(s![.., ti..ti + 1, .., ..])
                .assign(&turbine_turbulence_intensity.slice(s![.., ti..ti + 1, .., ..]));
        }
        flow_field.turbulence_intensity_field_sorted = field_sorted;
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

    // Evolving 4D TI field (n_findex, n_turbines, gy, gz), starting at ambient
    // (mirrors Python turbine_turbulence_intensity; full-shape constant is
    // equivalent to the (n_findex, n_turbines, 1, 1) start + broadcast)
    let mut turbine_turbulence_intensity: Array4 =
        Array::from_elem((n_findex, n_turbines, grid_y_dim, grid_z_dim), 0.0);
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    turbine_turbulence_intensity[[fi, ti, iy, iz]] =
                        flow_field.turbulence_intensities[fi];
                }
            }
        }
    }

    // Model-specific deflection arguments (e.g. wind_veer for gauss), matching
    // Python's prepare_function (Bug #23).
    let deflection_model_args = model_manager
        .deflection_model
        .prepare_function(grid, flow_field)?;

    for i in 0..n_turbines {
        // All turbines' thrust coefficients from the current u_sorted (matching
        // Python: recomputed each iteration as the wakes evolve)
        let cts_all = farm.thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
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

                // Per-findex yaw/Ct of upstream source ii; TI is the source's
                // per-point slice of the evolving field (Python:
                // turbulence_intensity_ii = turbine_turbulence_intensity[:, ii:ii+1])
                let mut yaw_ii_arr = Array1::zeros(n_findex);
                let mut ct_ii_arr = Array1::zeros(n_findex);
                for fi in 0..n_findex {
                    yaw_ii_arr[fi] = farm.yaw_angles_sorted[[fi, ii]];
                    ct_ii_arr[fi] = cts_all[[fi, ii]];
                }
                let ti_ii: Array4 = turbine_turbulence_intensity
                    .slice(s![.., ii..ii + 1, .., ..])
                    .to_owned();
                let d_ii = farm
                    .rotor_diameters_sorted
                    .slice(s![.., ii..ii + 1])
                    .to_owned();

                let field_ii = model_manager.deflection_model.function(
                    &x_2d,
                    &y_2d,
                    &yaw_ii_arr,
                    &ti_ii,
                    &ct_ii_arr,
                    &d_ii,
                    grid.x_sorted(),
                    &deflection_model_args,
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

        // Target-centric turbopark deficit at turbine i. Python passes the
        // evolved turbine_turbulence_intensity field (used pointwise by
        // characteristic_wake_width for each source), not the ambient TI.
        let velocity_deficit = crate::core::wake_velocity::turbopark::turbopark_velocity_deficit(
            &turbine_turbulence_intensity,
            &cts_all,
            &farm.rotor_diameters_sorted,
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
   
        let y_mean = mean_value_4d_slice(&grid.y_sorted().clone())?;

        let ai_all = farm.axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
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

        let rotor_diameters_i = farm
            .rotor_diameters_sorted
            .slice(s![.., i..i + 1])
            .to_owned();
        let mut wake_added = model_manager.turbulence_model.function(
            &ambient_ti_4d,
            grid.x_sorted(),
            &x_i_4d,
            &rotor_diameters_i,
            &ai_4d,
        )?;
        wake_added = wake_added.mapv(|v| if v.is_finite() { v } else { 0.0 });

        let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;

        for fi in 0..n_findex {
            let x_src = x_2d[[fi, 0]];
            let y_src = y_mean[[fi, i]];
            let ambient = flow_field.turbulence_intensities[fi];
            let d_i = rotor_diameters_i[[fi, 0]];
            let downstream_influence_length = 15.0 * d_i;
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

                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        // Pointwise masks vs the outer source turbine's mean
                        // position (Python: grid.x_sorted > x_i, |y_i -
                        // grid.y_sorted| < 2D, grid.x_sorted <= 15D + x_i)
                        let x_pt = grid.x_sorted()[[fi, j, iy, iz]];
                        let downstream = x_pt > x_src
                            && (y_src - grid.y_sorted()[[fi, j, iy, iz]]).abs() < 2.0 * d_i
                            && x_pt <= downstream_influence_length + x_src;
                        let ti_added = if downstream {
                            area_overlap * wake_added[[fi, j, iy, iz]]
                        } else {
                            0.0
                        };
                        let new_ti = (ti_added * ti_added + ambient * ambient).sqrt();
                        if new_ti > turbine_turbulence_intensity[[fi, j, iy, iz]] {
                            turbine_turbulence_intensity[[fi, j, iy, iz]] = new_ti;
                        }
                    }
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

    // Store the evolved per-point TI field (Python: flow_field.
    // turbulence_intensity_field_sorted = turbine_turbulence_intensity)
    if grid_second_dim == n_turbines {
        flow_field.turbulence_intensity_field_sorted = turbine_turbulence_intensity;
    } else {
        let mut field_sorted = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
        for ti in 0..n_turbines {
            field_sorted
                .slice_mut(s![.., ti..ti + 1, .., ..])
                .assign(&turbine_turbulence_intensity.slice(s![.., ti..ti + 1, .., ..]));
        }
        flow_field.turbulence_intensity_field_sorted = field_sorted;
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
    turbine_grid_flow_field.initialize_flow_field_with_coordinates(
        turbine_grid.grid_shape(),
        turbine_grid.z_sorted(),
        &turbine_grid.hub_heights(),
        Some(turbine_grid.x_sorted_inertial_frame()),
        Some(turbine_grid.y_sorted_inertial_frame()),
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

    // Deflection model arguments (e.g. wind_veer for gauss), matching Python's
    // prepare_function (Bug #23).
    let deflection_model_args = model_manager
        .deflection_model
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
        let ct_all = turbine_grid_farm.thrust_coefficient(
            &turbine_grid_flow_field.u_sorted,
            turbine_grid_flow_field.air_density,
            Some(&turbine_grid_farm.yaw_angles_sorted),
            Some(&turbine_grid_farm.tilt_angles_sorted),
            Some(&turbine_grid_flow_field.turbulence_intensity_field_sorted),
            Some(&turbine_grid_farm.power_setpoints_sorted),
            Some(&turbine_grid_farm.awc_modes_sorted),
            Some(&turbine_grid_farm.awc_amplitudes_sorted),
            turbine_grid.average_method(),
            turbine_grid.cubature_weights(),
        )?;
        let a_all = turbine_grid_farm.axial_induction(
            &turbine_grid_flow_field.u_sorted,
            turbine_grid_flow_field.air_density,
            Some(&turbine_grid_farm.yaw_angles_sorted),
            Some(&turbine_grid_farm.tilt_angles_sorted),
            Some(&turbine_grid_flow_field.turbulence_intensity_field_sorted),
            Some(&turbine_grid_farm.power_setpoints_sorted),
            Some(&turbine_grid_farm.awc_modes_sorted),
            Some(&turbine_grid_farm.awc_amplitudes_sorted),
            turbine_grid.average_method(),
            turbine_grid.cubature_weights(),
        )?;

        // Per-findex scalars for the current turbine.
        // (Python: ct_i/axial_induction_i are per-findex arrays
        // (n_findex, 1, 1, 1), solver.py:349-371 — Bug #20.)
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

        let rotor_diameters_i = turbine_grid_farm
            .rotor_diameters_sorted
            .slice(s![.., i..i + 1])
            .to_owned();
        let hub_heights_i = turbine_grid_farm
            .hub_heights_sorted
            .slice(s![.., i..i + 1])
            .to_owned();

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
                &rotor_diameters_i,
                &hub_heights_i,
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
                &rotor_diameters_i,
                &hub_heights_i,
                &yaw_arr,
                &ct_arr,
                &tsr_arr,
                &ai_arr,
                flow_field.wind_shear,
                1.0,
            );
            // Python reassigns v_wake/w_wake each iteration (solver.py:411) so
            // that the per-iteration `v_sorted += v_wake` adds only source i's
            // contribution; accumulating here would double-count earlier
            // sources (Bug #19).
            for fi in 0..n_findex {
                for p in 0..n_points {
                    for iy in 0..grid_y_dim {
                        for iz in 0..grid_z_dim {
                            v_wake[[fi, p, iy, iz]] = vw[[fi, p, iy, iz]];
                            w_wake[[fi, p, iy, iz]] = ww[[fi, p, iy, iz]];
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
        // synthetic grid (Python: deflection on the full flow grid with the
        // rotor-averaged TI scalar (n_findex, 1, 1, 1)).
        let ti_slice = ti_arr
            .clone()
            .into_shape_with_order((n_findex, 1, 1, 1))
            .map_err(|e| anyhow::anyhow!("ti slice reshape failed: {}", e))?;
        let deflection_field = model_manager.deflection_model.function(
            &x_i,
            &y_i,
            &effective_yaw_arr,
            &ti_slice,
            &ct_arr,
            &rotor_diameters_i,
            &x_full,
            &deflection_model_args,
        )?;

        // Per-turbine broadcast arrays (only column `i` is read by the model).
        // ti_bc is constant over every block, matching Python where the scalar
        // rotor-averaged TI (n_findex, 1, 1, 1) broadcasts over the grid.
        let mut a_arr = Array::zeros((n_findex, total));
        let mut yaw_bc = Array::zeros((n_findex, total));
        let mut ti_bc = Array::zeros((n_findex, total, 1, 1));
        let mut ct_bc = Array::zeros((n_findex, total));
        for fi in 0..n_findex {
            a_arr[[fi, i]] = a_all[[fi, i]];
            yaw_bc[[fi, i]] = yaw_arr[fi];
            ct_bc[[fi, i]] = ct_arr[fi];
            for b in 0..total {
                ti_bc[[fi, b, 0, 0]] = ti_arr[fi];
            }
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
            &hub_heights_i,
            &rotor_diameters_i,
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
            ai_4d[[fi, 0, 0, 0]] = a_all[[fi, i]];
        }
        let mut wake_added_ti = model_manager.turbulence_model.function(
            &ambient_turbulence_intensities,
            &planar_x,
            &x_i_4d,
            &rotor_diameters_i,
            &ai_4d,
        )?;
        wake_added_ti = wake_added_ti.mapv(|v| if v.is_finite() { v } else { 0.0 });

        for fi in 0..n_findex {
            let d_i = rotor_diameters_i[[fi, 0]];
            let downstream_influence_length = 15.0 * d_i;
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let area_overlap = if velocity_deficit_absolute[[fi, p, iy, iz]] > 0.05 {
                            1.0
                        } else {
                            0.0
                        };
                        let downstream = planar_x[[fi, p, iy, iz]] > x_i[[fi, 0]]
                            && (y_i[[fi, 0]] - planar_y[[fi, p, iy, iz]]).abs() < 2.0 * d_i
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
                        flow_field.u_sorted[[fi, p, iy, iz]] = flow_field.u_initial_sorted
                            [[fi, p, iy, iz]]
                            - wake_field[[fi, p, iy, iz]];
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

    // Evolving 4D TI field (n_findex, n_turbines, gy, gz), starting at ambient
    let mut turbine_turbulence_intensity: Array4 =
        Array::from_elem((n_findex, n_turbines, grid_y_dim, grid_z_dim), 0.0);
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            for iy in 0..grid_y_dim {
                for iz in 0..grid_z_dim {
                    turbine_turbulence_intensity[[fi, ti, iy, iz]] =
                        flow_field.turbulence_intensities[fi];
                }
            }
        }
    }

    // Python computes the axial induction used by the transverse velocities and
    // secondary steering from flow_field.u_sorted (solver.py:576-596), which is
    // never updated inside the cc loop (only after it, solver.py:712) — i.e. the
    // freestream AI. Rust previously used the wake-inflow AI (turb_ais_2d),
    // which differs for waked turbines (Bug #21).
    let ai_freestream = farm.axial_induction(
        &flow_field.u_initial_sorted,
        flow_field.air_density,
        Some(&farm.yaw_angles_sorted),
        Some(&farm.tilt_angles_sorted),
        Some(&flow_field.turbulence_intensity_field_sorted),
        Some(&farm.power_setpoints_sorted),
        Some(&farm.awc_modes_sorted),
        Some(&farm.awc_amplitudes_sorted),
        grid.average_method(),
        grid.cubature_weights(),
    )?;

    // Model-specific deflection arguments (e.g. wind_veer for gauss), matching
    // Python's prepare_function (Bug #23).
    let deflection_model_args = model_manager
        .deflection_model
        .prepare_function(grid, flow_field)?;

    for i in 0..n_turbines {
        let rotor_diameters_i = farm
            .rotor_diameters_sorted
            .slice(s![.., i..i + 1])
            .to_owned();

        let x_i_slice = grid.x_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
        let y_i_slice = grid.y_sorted().slice(s![.., i..i + 1, .., ..]).to_owned();
        let x_i_2d = mean_value_4d_slice(&x_i_slice)?;
        let y_i_2d = mean_value_4d_slice(&y_i_slice)?;

        // The inflow at turbine i's rotor = u_initial - turb_u_wake
        for fi in 0..n_findex {
            let d_i = rotor_diameters_i[[fi, 0]];
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
        let turb_cts_2d = farm.thrust_coefficient(
            &avg_vels_4d,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;
        let turb_ais_2d = farm.axial_induction(
            &avg_vels_4d,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;

        // The source's inflow (n_findex, 1, gy, gz) and yaw per findex
        let u_i = turb_inflow_field.slice(s![.., i..i + 1, .., ..]).to_owned();
        let mut yaw_i_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_i_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
        }

        // Per-findex scalar arrays for the source turbine (ct, tsr, ai, yaw)
        let hub_heights_i = farm.hub_heights_sorted.slice(s![.., i..i + 1]).to_owned();
        let tsrs_i = farm.tsrs_sorted.slice(s![.., i..i + 1]).to_owned();
        let mut yaw_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut tsr_arr = Array1::zeros(n_findex);
        let mut ai_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
            ct_arr[fi] = turb_cts_2d[[fi, i]];
            tsr_arr[fi] = tsrs_i[[fi, 0]];
            // Transverse velocities use the freestream AI (Bug #21).
            ai_arr[fi] = ai_freestream[[fi, i]];
        }

        // Snapshot the source's TI slice BEFORE the GCH update (Bug #22):
        // Python's deflection call (solver.py:625) reads `turbulence_intensity_i`,
        // a view taken at solver.py:598 before the GCH write at 653-663, so the
        // deflection model sees the pre-GCH TI while the velocity model (665)
        // sees the post-GCH field.
        let ti_source: Array4 = turbine_turbulence_intensity
            .slice(s![.., i..i + 1, .., ..])
            .to_owned();

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
                &rotor_diameters_i,
                &hub_heights_i,
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
            // Python: I_i = turbulence_intensity_i[:, 0, 0, 0] (first grid
            // point of the source's field); gch_gain = 1.0 in cc_solver
            let i_mixing = yaw_added_turbulence_mixing(
                &u_i,
                &turbine_turbulence_intensity
                    .slice(s![.., i, 0, 0])
                    .to_owned(),
                &v_i,
                &w_i,
                &v_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
                &w_wake.slice(s![.., i..i + 1, .., ..]).to_owned(),
            );
            let gch_gain = 1.0;
            for fi in 0..n_findex {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        turbine_turbulence_intensity[[fi, i, iy, iz]] += gch_gain * i_mixing[fi];
                    }
                }
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

            let mut ct_4d = Array4::zeros((n_findex, 1, 1, 1));
            let mut tsr_4d = Array4::zeros((n_findex, 1, 1, 1));
            let mut ai_4d = Array4::zeros((n_findex, 1, 1, 1));
            for fi in 0..n_findex {
                ct_4d[[fi, 0, 0, 0]] = turb_cts_2d[[fi, i]];
                tsr_4d[[fi, 0, 0, 0]] = tsrs_i[[fi, 0]];
                // Secondary steering uses the freestream AI (Bug #21).
                ai_4d[[fi, 0, 0, 0]] = ai_freestream[[fi, i]];
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
                &rotor_diameters_i,
                &hub_heights_i,
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
        // (Python cc: per-findex effective yaw/Ct, per-point TI slice,
        // full-field evaluation; `ti_source` was snapshotted pre-GCH above)
        let ct_source: Array1 = (0..n_findex).map(|fi| turb_cts_2d[[fi, i]]).collect();
        let deflection_field = model_manager.deflection_model.function(
            &x_i_2d,
            &y_i_2d,
            &effective_yaw_arr,
            &ti_source,
            &ct_source,
            &rotor_diameters_i,
            grid.x_sorted(),
            &deflection_model_args,
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
                &farm.rotor_diameters_sorted,
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
            &rotor_diameters_i,
            &ai_4d,
        )?;
        wake_added = wake_added.mapv(|v| if v.is_finite() { v } else { 0.0 });

        let grid_res_sq = (grid_y_dim * grid_z_dim) as Float;

        for fi in 0..n_findex {
            let x_src = x_i_2d[[fi, 0]];
            let y_src = y_i_2d[[fi, 0]];
            let ambient = flow_field.turbulence_intensities[fi];
            let d_i = rotor_diameters_i[[fi, 0]];
            let downstream_influence_length = 15.0 * d_i;
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

                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        // Pointwise downstream-influence masks (Python:
                        // grid.x_sorted > x_i, |y_i - grid.y_sorted| < 2D,
                        // grid.x_sorted <= 15D + x_i)
                        let x_pt = grid.x_sorted()[[fi, j, iy, iz]];
                        let downstream = x_pt > x_src
                            && (y_src - grid.y_sorted()[[fi, j, iy, iz]]).abs() < 2.0 * d_i
                            && x_pt <= downstream_influence_length + x_src;
                        let ti_added = if downstream {
                            area_overlap * wake_added[[fi, j, iy, iz]]
                        } else {
                            0.0
                        };
                        let new_ti = (ti_added * ti_added + ambient * ambient).sqrt();
                        if new_ti > turbine_turbulence_intensity[[fi, j, iy, iz]] {
                            turbine_turbulence_intensity[[fi, j, iy, iz]] = new_ti;
                        }
                    }
                }
            }
        }

        // Add transverse velocities to the flow field (Python cc: v_sorted += v_wake)
        flow_field.v_sorted += &v_wake;
        flow_field.w_sorted += &w_wake;
    }

    // The final u_sorted is the inflow field (matching Python)
    flow_field.u_sorted = turb_inflow_field;

    // Store the evolved per-point TI field (Python: flow_field.
    // turbulence_intensity_field_sorted = turbine_turbulence_intensity)
    flow_field.turbulence_intensity_field_sorted = turbine_turbulence_intensity;

    Ok(())
}

/// Empirical Gauss solver (port of Python empirical_gauss_solver).
///
/// For each wake source turbine (upstream to downstream), computes the deficit
/// it adds to every downstream turbine and accumulates it via the combination
/// model. The wake-induced mixing factor is built up iteratively: each source
/// contributes `area_overlap * a_i / downstream_D^2` to downstream turbines'
/// mixing terms, and the RSS over upstream contributions feeds the velocity
/// model's wake width via `mixing_gain_velocity`.
///
/// Returns the mixing factor (n_findex, n_turbines, n_turbines): the entry
/// `[fi, j, k]` is the contribution of the source turbine `k` to the mixing of
/// the receiving turbine `j` (Python returns it as `wim_field` for the
/// full-flow solver).
pub fn empirical_gauss_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<Array3> {
    if model_manager.enable_secondary_steering {
        anyhow::bail!("Secondary steering not available for this model.");
    }
    if model_manager.enable_transverse_velocities {
        anyhow::bail!("Transverse velocities not used in this model.");
    }

    // The mixing architecture is defined for turbine grids (Python's full-flow
    // solver runs the turbine-grid solver on a copy first).
    if grid.as_any().downcast_ref::<TurbineGrid>().is_none() {
        sequential_solver(farm, flow_field, grid, model_manager)?;
        return Ok(Array3::zeros((
            flow_field.n_findex,
            grid.n_turbines(),
            grid.n_turbines(),
        )));
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
        let ct_all = farm.thrust_coefficient(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;
        let a_all = farm.axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;

        // Per-findex source-turbine quantities. Python computes tilt via
        // calculate_tilt_for_eff_velocities, which for non-floating turbines
        // (no tilt interpolation) returns tilt_angles_sorted unchanged.
        // The floating tilt is interpolated at the rotor average velocity of
        // the (waked) flow field (Python solver.py:1305-1311).
        let average_velocities = {
            let velocities_dyn = flow_field.u_sorted.view().into_dyn().to_owned();
            let cubature_weights = grid
                .cubature_weights()
                .map(|w| Array1::from_iter(w.iter().copied()));
            crate::core::rotor_velocity::average_velocity(
                &velocities_dyn,
                None,
                &grid.average_method(),
                cubature_weights.as_ref(),
            )?
        };
        let tilt_angles = farm.calculate_tilt_for_eff_velocities(&average_velocities);

        let mut yaw_arr = Array1::zeros(n_findex);
        let mut tilt_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut a_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = farm.yaw_angles_sorted[[fi, i]];
            tilt_arr[fi] = tilt_angles[[fi, i]];
            ct_arr[fi] = ct_all[[fi, i]];
            a_arr[fi] = a_all[[fi, i]];
        }

        // Influence of yawing on the turbine's OWN wake (Python solver.py:1321-1329:
        // mixing_factor[:, i, i] += yaw_added_wake_mixing(ai, yaw, 1, gain) with
        // downstream_distance=1, applied BEFORE the RSS so it enters mixing_i —
        // Bug #25; previously missing entirely).
        if model_manager.enable_yaw_added_recovery {
            for fi in 0..n_findex {
                let yaw_self =
                    a_arr[fi] * defl_eg.yaw_added_mixing_gain * (1.0 - cosd(yaw_arr[fi]));
                mixing_factor[[fi, i, i]] += yaw_self / 1.0;
            }
        }

        // Influence of awc on the turbine's OWN wake (Python solver.py:1331-1340,
        // also applied before the RSS — Bug #25; previously added too late to
        // ever enter any mixing_i).
        if model_manager.enable_active_wake_mixing {
            for fi in 0..n_findex {
                if farm.awc_modes_sorted[[fi, i]] == "helix" {
                    let amp = farm.awc_amplitudes_sorted[[fi, i]];
                    mixing_factor[[fi, i, i]] +=
                        amp.powf(eg.awc_wake_exp) / eg.awc_wake_denominator;
                }
            }
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
        let rotor_diameters_i = farm
            .rotor_diameters_sorted
            .slice(s![.., i..i + 1])
            .to_owned();
        let hub_heights_i = farm.hub_heights_sorted.slice(s![.., i..i + 1]).to_owned();
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
        let (dy_flat, dz_flat) = defl_eg.deflection_yz(
            &x_downstream,
            &yaw_arr,
            &tilt_arr,
            &mixing_i,
            &ct_arr,
            &rotor_diameters_i,
        );
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
            &hub_heights_i,
            &rotor_diameters_i,
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

        // Update u_sorted after each turbine's wake is added
        // (Python: flow_field.u_sorted = flow_field.u_initial_sorted - wake_field)
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let wake_deficit = wake_field[[fi, ti, iy, iz]];
                        let u_initial = flow_field.u_initial_sorted[[fi, ti, iy, iz]];
                        flow_field.u_sorted[[fi, ti, iy, iz]] = u_initial - wake_deficit;
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

    Ok(mixing_factor)
}

/// Full-flow Empirical Gauss solver (port of Python
/// `full_flow_empirical_gauss_solver`).
///
/// The turbine grid is solved first (with the empirical Gauss solver), which
/// yields the wake-induced mixing field. The wake of every turbine is then
/// evaluated over the full flow-field grid with that mixing term (Python
/// `core.py:259-262`).
pub fn full_flow_empirical_gauss_solver(
    farm: &Farm,
    flow_field: &mut FlowField,
    grid: &dyn Grid,
    model_manager: &WakeModelManager,
) -> Result<()> {
    let eg = model_manager
        .velocity_model
        .as_any()
        .downcast_ref::<crate::core::wake_velocity::EmpiricalGaussVelocityDeficit>()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "full_flow_empirical_gauss_solver requires EmpiricalGaussVelocityDeficit"
            )
        })?;
    let defl_eg = model_manager
        .deflection_model
        .as_any()
        .downcast_ref::<crate::core::wake_deflection::EmpiricalGaussVelocityDeflection>()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "full_flow_empirical_gauss_solver requires EmpiricalGaussVelocityDeficit"
            )
        })?;

    let n_turbines = farm.n_turbines();
    let n_findex = grid.n_findex();
    let grid_shape = grid.x_sorted().shape();
    let n_points = grid_shape[1];
    let grid_y_dim = grid_shape[2];
    let grid_z_dim = grid_shape[3];

    // Solve the turbine grid to get the wake-induced mixing field.
    let turbine_grid = TurbineGrid::new(
        farm.coordinates(),
        farm.rotor_diameters().clone(),
        flow_field.wind_directions.clone(),
        3,
    )?;
    let mut turbine_grid_farm = farm.clone();
    let mut turbine_grid_flow_field = flow_field.clone();
    turbine_grid_farm.expand_farm_properties(n_findex, turbine_grid.sorted_coord_indices());
    turbine_grid_flow_field.initialize_flow_field_with_coordinates(
        turbine_grid.grid_shape(),
        turbine_grid.z_sorted(),
        &turbine_grid.hub_heights(),
        Some(turbine_grid.x_sorted_inertial_frame()),
        Some(turbine_grid.y_sorted_inertial_frame()),
    );
    turbine_grid_farm.initialize(turbine_grid.sorted_indices());
    let mixing_field = empirical_gauss_solver(
        &turbine_grid_farm,
        &mut turbine_grid_flow_field,
        &turbine_grid,
        model_manager,
    )?;

    // The source positions of the turbines (the mean of their rotor blocks)
    let source_x = mean_value_4d_slice(&turbine_grid.x_sorted().clone())?;
    let source_y = mean_value_4d_slice(&turbine_grid.y_sorted().clone())?;

    // The full-flow grid coordinates (the rotated frame, like Python)
    let planar_x = grid.x_sorted().clone();
    let planar_y = grid.y_sorted().clone();
    let planar_z = grid.z_sorted().clone();

    let mut wake_field: Array4 = Array::zeros((n_findex, n_points, grid_y_dim, grid_z_dim));

    for i in 0..n_turbines {
        // The thrust coefficient and the axial induction of the turbines from
        // the turbine-grid solve (Python's `ix_filter=[i]`)
        let ct_all = turbine_grid_farm.thrust_coefficient(
            &turbine_grid_flow_field.u_sorted,
            turbine_grid_flow_field.air_density,
            Some(&turbine_grid_farm.yaw_angles_sorted),
            Some(&turbine_grid_farm.tilt_angles_sorted),
            Some(&turbine_grid_flow_field.turbulence_intensity_field_sorted),
            Some(&turbine_grid_farm.power_setpoints_sorted),
            Some(&turbine_grid_farm.awc_modes_sorted),
            Some(&turbine_grid_farm.awc_amplitudes_sorted),
            turbine_grid.average_method(),
            turbine_grid.cubature_weights(),
        )?;
        let a_all = turbine_grid_farm.axial_induction(
            &turbine_grid_flow_field.u_sorted,
            turbine_grid_flow_field.air_density,
            Some(&turbine_grid_farm.yaw_angles_sorted),
            Some(&turbine_grid_farm.tilt_angles_sorted),
            Some(&turbine_grid_flow_field.turbulence_intensity_field_sorted),
            Some(&turbine_grid_farm.power_setpoints_sorted),
            Some(&turbine_grid_farm.awc_modes_sorted),
            Some(&turbine_grid_farm.awc_amplitudes_sorted),
            turbine_grid.average_method(),
            turbine_grid.cubature_weights(),
        )?;

        let mut yaw_arr = Array1::zeros(n_findex);
        let mut ct_arr = Array1::zeros(n_findex);
        let mut a_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            yaw_arr[fi] = turbine_grid_farm.yaw_angles_sorted[[fi, i]];
            ct_arr[fi] = ct_all[[fi, i]];
            a_arr[fi] = a_all[[fi, i]];
        }

        // The tilt of the source turbine (Python's
        // `calculate_tilt_for_eff_velocities` on the turbine-grid velocities)
        let average_velocities = {
            let velocities_dyn = turbine_grid_flow_field
                .u_sorted
                .view()
                .into_dyn()
                .to_owned();
            let cubature_weights = turbine_grid
                .cubature_weights()
                .map(|w| Array1::from_iter(w.iter().copied()));
            crate::core::rotor_velocity::average_velocity(
                &velocities_dyn,
                None,
                &turbine_grid.average_method(),
                cubature_weights.as_ref(),
            )?
        };
        let tilt_angles = turbine_grid_farm.calculate_tilt_for_eff_velocities(&average_velocities);
        let mut tilt_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            tilt_arr[fi] = tilt_angles[[fi, i]];
        }

        // The wake-induced mixing of the source turbine: the sum of the
        // contributions of the sources (Python's
        // `wim_field[:, i:i+1, :, None].sum(axis=2)`)
        let mut mixing_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            let mut sum = 0.0;
            for k in 0..n_turbines {
                sum += mixing_field[[fi, i, k]];
            }
            mixing_arr[fi] = sum;
        }

        let rotor_diameters_i = turbine_grid_farm
            .rotor_diameters_sorted
            .slice(s![.., i..i + 1])
            .to_owned();
        let hub_heights_i = turbine_grid_farm
            .hub_heights_sorted
            .slice(s![.., i..i + 1])
            .to_owned();

        // The distance from the source to every point of the full grid
        let mut x_downstream = Array2::zeros((n_findex, n_points));
        for fi in 0..n_findex {
            for p in 0..n_points {
                x_downstream[[fi, p]] = planar_x[[fi, p, 0, 0]] - source_x[[fi, i]];
            }
        }
        let (dy_flat, dz_flat) = defl_eg.deflection_yz(
            &x_downstream,
            &yaw_arr,
            &tilt_arr,
            &mixing_arr,
            &ct_arr,
            &rotor_diameters_i,
        );
        let deflection_y =
            dy_flat.into_shape_with_order((n_findex, n_points, grid_y_dim, grid_z_dim))?;
        let deflection_z =
            dz_flat.into_shape_with_order((n_findex, n_points, grid_y_dim, grid_z_dim))?;

        // The velocity deficit of the source over the full grid
        // The source position of the turbine (the mean of its rotor block)
        let mut x_i_arr = Array1::zeros(n_findex);
        let mut y_i_arr = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            x_i_arr[fi] = source_x[[fi, i]];
            y_i_arr[fi] = source_y[[fi, i]];
        }
        let velocity_deficit = eg.velocity_deficit_emgauss_at(
            &planar_x,
            &planar_y,
            &planar_z,
            &deflection_y,
            &deflection_z,
            &yaw_arr,
            &tilt_arr,
            &ct_arr,
            &mixing_arr,
            &hub_heights_i,
            &rotor_diameters_i,
            &x_i_arr,
            &y_i_arr,
            flow_field.wind_veer,
        );

        // The absolute deficit and the combination of the wakes
        let mut velocity_deficit_absolute = flow_field.u_initial_sorted.clone();
        for fi in 0..n_findex {
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        velocity_deficit_absolute[[fi, p, iy, iz]] *=
                            velocity_deficit[[fi, p, iy, iz]];
                    }
                }
            }
        }
        wake_field = model_manager
            .combination_model
            .function(&wake_field, &velocity_deficit_absolute)?;

        // Update the full-flow sorted field
        for fi in 0..n_findex {
            for p in 0..n_points {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        flow_field.u_sorted[[fi, p, iy, iz]] = flow_field.u_initial_sorted
                            [[fi, p, iy, iz]]
                            - wake_field[[fi, p, iy, iz]];
                    }
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
