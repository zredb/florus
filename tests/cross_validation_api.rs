mod common;

use florus::core::grid::TurbineGrid;
use florus::core::solver::sequential_solver;
use florus::core::turbines::cp_ct_table::TableConditions;
use florus::floris_config::{
    CombinationModelConfig, DeflectionModelConfig, FarmConfig, FlorisConfig, FlowFieldConfig,
    LoggingConfig, SolverConfig, TurbulenceModelConfig, VelocityModelConfig, WakeConfig,
    WakeModelStringsConfig,
};
use florus::floris_model::FlorisModel;
use florus::types::{Array2, Float};
use ndarray::Array;
use serde::Deserialize;

// ============================================================
// Golden value struct for air_density_effect (non-standard format)
// ============================================================

#[derive(Debug, Deserialize)]
struct AirDensityGolden {
    power_air_density_1225: Vec<Vec<f64>>,
    power_air_density_10: Vec<Vec<f64>>,
    ratio: f64,
}

// ============================================================
// Config builder matching Python FLORIS API golden value generator
// ============================================================

const API_WIND_DIRECTIONS: [f64; 1] = [270.0];
const API_WIND_SPEEDS_8: [f64; 1] = [8.0];
const API_WIND_SPEEDS_10: [f64; 1] = [10.0];
const API_TURBULENCE_INTENSITIES: [f64; 1] = [0.06];

// 2-turbine layout for single_condition_power, yaw_effect, no_wake_power
const LAYOUT_2T_X: [f64; 2] = [0.0, 630.0];
const LAYOUT_2T_Y: [f64; 2] = [0.0, 0.0];

// 3-turbine layout for air_density_effect (matches Python default)
const LAYOUT_3T_X: [f64; 3] = [0.0, 630.0, 1260.0];
const LAYOUT_3T_Y: [f64; 3] = [0.0, 0.0, 0.0];

fn create_api_config_2t(wind_speeds: Vec<f64>, air_density: f64) -> FlorisConfig {
    let turbine_types = vec!["nrel_5MW".to_string(); 2];

    FlorisConfig {
        name: "api_cross_validation".to_string(),
        description: Some("API cross-validation test config".to_string()),
        floris_version: "v4".to_string(),
        logging: LoggingConfig::default(),
        solver: SolverConfig::default(),
        farm: FarmConfig {
            layout_x: LAYOUT_2T_X.to_vec(),
            layout_y: LAYOUT_2T_Y.to_vec(),
            turbine_type: turbine_types,
        },
        flow_field: FlowFieldConfig {
            wind_speeds,
            wind_directions: API_WIND_DIRECTIONS.to_vec(),
            turbulence_intensities: API_TURBULENCE_INTENSITIES.to_vec(),
            wind_shear: 0.12,
            wind_veer: 0.0,
            air_density,
            reference_wind_height: 90.0,
            multidim_conditions: None,
            heterogeneous_inflow_config: None,
        },
        wake: WakeConfig {
            model_strings: WakeModelStringsConfig {
                velocity_model: VelocityModelConfig::Gauss {
                    alpha: 0.58,
                    beta: 0.077,
                    ka: 0.38,
                    kb: 0.004,
                },
                deflection_model: DeflectionModelConfig::Gauss {
                    alpha: 0.58,
                    beta: 0.077,
                    ad: 0.0,
                    bd: 0.0,
                    dm: 1.0,
                    ka: 0.38,
                    kb: 0.004,
                },
                combination_model: CombinationModelConfig::FLS,
                turbulence_model: TurbulenceModelConfig::CrespoHernandez {
                    initial: 0.1,
                    constant: 0.5,
                    ai: 0.8,
                    downstream: -0.32,
                },
            },
            enable_secondary_steering: true,
            enable_yaw_added_recovery: true,
            enable_transverse_velocities: true,
            enable_active_wake_mixing: false,
            enable_wake_mixing: false,
            use_parallel_calc: false,
        },
        turbine_library: "turbine_library".to_string(),
    }
}

fn create_api_config_3t(wind_speeds: Vec<f64>, air_density: f64) -> FlorisConfig {
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    FlorisConfig {
        name: "api_cross_validation_3t".to_string(),
        description: Some("API cross-validation test config 3 turbines".to_string()),
        floris_version: "v4".to_string(),
        logging: LoggingConfig::default(),
        solver: SolverConfig::default(),
        farm: FarmConfig {
            layout_x: LAYOUT_3T_X.to_vec(),
            layout_y: LAYOUT_3T_Y.to_vec(),
            turbine_type: turbine_types,
        },
        flow_field: FlowFieldConfig {
            wind_speeds,
            wind_directions: API_WIND_DIRECTIONS.to_vec(),
            turbulence_intensities: API_TURBULENCE_INTENSITIES.to_vec(),
            wind_shear: 0.12,
            wind_veer: 0.0,
            air_density,
            reference_wind_height: 90.0,
            multidim_conditions: None,
            heterogeneous_inflow_config: None,
        },
        wake: WakeConfig {
            model_strings: WakeModelStringsConfig {
                velocity_model: VelocityModelConfig::Gauss {
                    alpha: 0.58,
                    beta: 0.077,
                    ka: 0.38,
                    kb: 0.004,
                },
                deflection_model: DeflectionModelConfig::Gauss {
                    alpha: 0.58,
                    beta: 0.077,
                    ad: 0.0,
                    bd: 0.0,
                    dm: 1.0,
                    ka: 0.38,
                    kb: 0.004,
                },
                combination_model: CombinationModelConfig::FLS,
                turbulence_model: TurbulenceModelConfig::CrespoHernandez {
                    initial: 0.1,
                    constant: 0.5,
                    ai: 0.8,
                    downstream: -0.32,
                },
            },
            enable_secondary_steering: true,
            enable_yaw_added_recovery: true,
            enable_transverse_velocities: true,
            enable_active_wake_mixing: false,
            enable_wake_mixing: false,
            use_parallel_calc: false,
        },
        turbine_library: "turbine_library".to_string(),
    }
}

// ============================================================
// Manual solver power computation (matching regression test approach)
// ============================================================

fn compute_powers_manual(config: &FlorisConfig, yaw_angles: Option<&Array2>) -> Vec<Vec<f64>> {
    let mut model = FlorisModel::from_config(config.clone()).expect("Failed to create FlorisModel");

    if let Some(yaw) = yaw_angles {
        model
            .set_yaw_angles(yaw.clone())
            .expect("Failed to set yaw angles");
    }

    model.initialize_grid().expect("init grid failed");
    model
        .initialize_flow_field()
        .expect("init flow field failed");

    let grid = model.core().grid.as_ref().unwrap();
    let turbine_grid = grid.as_any().downcast_ref::<TurbineGrid>().unwrap();

    let farm = model.core().farm.clone();
    let mut flow_field = model.core().flow_field.clone();
    let wake = model.core().wake.clone();

    sequential_solver(&farm, &mut flow_field, turbine_grid, &wake).expect("solver failed");

    let n_findex = flow_field.u_sorted.shape()[0];
    let n_turbines = farm.turbines.len();
    let mut powers: Vec<Vec<f64>> = Vec::new();

    for fi in 0..n_findex {
        let mut row: Vec<f64> = Vec::new();
        for ti in 0..n_turbines {
            let mut v_cubed_sum = 0.0;
            let grid_points = flow_field.u_sorted.shape()[2] * flow_field.u_sorted.shape()[3];
            for iy in 0..flow_field.u_sorted.shape()[2] {
                for iz in 0..flow_field.u_sorted.shape()[3] {
                    v_cubed_sum += flow_field.u_sorted[[fi, ti, iy, iz]].powi(3);
                }
            }
            let rotor_avg_velocity = (v_cubed_sum / grid_points as f64).powf(1.0 / 3.0);

            let turbine = &farm.turbines[ti];
            let ref_air_density = turbine
                .turbine_type
                .power_thrust_table
                .ref_air_density
                .unwrap_or(1.225);
            let density_factor = (flow_field.air_density / ref_air_density).powf(1.0 / 3.0);
            let mut rotor_effective_velocity = rotor_avg_velocity * density_factor;

            let cosine_loss_exponent_yaw = turbine
                .turbine_type
                .power_thrust_table
                .cosine_loss_exponent_yaw
                .unwrap_or(1.88);
            let yaw = farm.yaw_angles[[fi, ti]];
            let yaw_correction = cosd(yaw).powf(cosine_loss_exponent_yaw / 3.0);
            rotor_effective_velocity *= yaw_correction;

            let cosine_loss_exponent_tilt = turbine
                .turbine_type
                .power_thrust_table
                .cosine_loss_exponent_tilt
                .unwrap_or(1.88);
            let ref_tilt = turbine
                .turbine_type
                .power_thrust_table
                .ref_tilt
                .unwrap_or(5.0);
            let tilt = farm.tilt_angles[[fi, ti]];
            let tilt_correction =
                (cosd(tilt) / cosd(ref_tilt)).powf(cosine_loss_exponent_tilt / 3.0);
            rotor_effective_velocity *= tilt_correction;

            let mut conditions = TableConditions::default();
            conditions.wind_speed = rotor_effective_velocity;
            let power_kw = turbine
                .turbine_type
                .power_thrust_table
                .cp_ct_table
                .get_cp(&conditions)
                .unwrap_or(0.0);
            let power_val = power_kw * 1000.0;

            row.push(power_val);
        }
        powers.push(row);
    }

    powers
}

fn cosd(degrees: f64) -> f64 {
    (degrees.to_radians()).cos()
}

// ============================================================
// Tolerance check — 30% relative or 250000W absolute
// ============================================================

fn power_close(test: f64, golden: f64) -> bool {
    let abs_diff = (test - golden).abs();
    let rel_diff = if golden.abs() > 1e-10 {
        abs_diff / golden.abs()
    } else {
        abs_diff
    };
    abs_diff <= 250000.0 || rel_diff <= 0.30
}

fn assert_powers_close(test: &[Vec<f64>], golden: &[Vec<f64>], name: &str) {
    assert_eq!(test.len(), golden.len(), "{}: findex count mismatch", name);
    for (fi, (t_row, g_row)) in test.iter().zip(golden.iter()).enumerate() {
        assert_eq!(
            t_row.len(),
            g_row.len(),
            "{}[{}]: turbine count mismatch",
            name,
            fi
        );
        for (ti, (t_val, g_val)) in t_row.iter().zip(g_row.iter()).enumerate() {
            assert!(
                power_close(*t_val, *g_val),
                "{}[{}][{}]: power mismatch: got {:.2}, expected {:.2}, diff {:.2} ({:.4}%)",
                name,
                fi,
                ti,
                t_val,
                g_val,
                (t_val - g_val).abs(),
                if g_val.abs() > 1e-10 {
                    (t_val - g_val).abs() / g_val.abs() * 100.0
                } else {
                    0.0
                }
            );
        }
    }
}

// ============================================================
// Test: single_condition_power
// Config: 2 turbines, 1 findex, 8m/s@270°, TI=0.06, Gauss model
// ============================================================

#[test]
fn test_single_condition_power() {
    let config = create_api_config_2t(API_WIND_SPEEDS_8.to_vec(), 1.225);

    let powers = compute_powers_manual(&config, None);
    let golden = common::load_api_golden("single_condition_power");
    assert_powers_close(&powers, &golden.turbine_powers, "single_condition_power");
}

// ============================================================
// Test: yaw_effect
// Config: 2 turbines, 1 findex, 8m/s@270°, TI=0.06, yaw=[10°, 0°]
// ============================================================

#[test]
fn test_yaw_effect() {
    let config = create_api_config_2t(API_WIND_SPEEDS_8.to_vec(), 1.225);

    let yaw = Array2::from_shape_fn((1, 2), |(fi, ti)| if ti == 0 { 10.0 } else { 0.0 });

    let powers = compute_powers_manual(&config, Some(&yaw));
    let golden = common::load_api_golden("yaw_effect");
    assert_powers_close(&powers, &golden.turbine_powers, "yaw_effect");
}

// ============================================================
// Test: air_density_effect
// Config: 3 turbines, 1 findex, 10m/s@270°, TI=0.06
// Two runs: air_density=1.225 and air_density=1.0
// Also verifies the power ratio between the two
// ============================================================

#[test]
fn test_air_density_effect() {
    let config_1225 = create_api_config_3t(API_WIND_SPEEDS_10.to_vec(), 1.225);
    let config_10 = create_api_config_3t(API_WIND_SPEEDS_10.to_vec(), 1.0);

    let powers_1225 = compute_powers_manual(&config_1225, None);
    let powers_10 = compute_powers_manual(&config_10, None);

    let golden: AirDensityGolden = serde_json::from_str(
        &std::fs::read_to_string("tests/golden_values/floris_model_api/air_density_effect.json")
            .unwrap(),
    )
    .unwrap();

    assert_powers_close(
        &powers_1225,
        &golden.power_air_density_1225,
        "air_density_1225",
    );
    assert_powers_close(&powers_10, &golden.power_air_density_10, "air_density_10");

    // Verify ratio: power at rho=1.0 / power at rho=1.225 ≈ 0.822
    let total_1225: f64 = powers_1225[0].iter().sum();
    let total_10: f64 = powers_10[0].iter().sum();
    let ratio = total_10 / total_1225;

    assert!(
        (ratio - golden.ratio).abs() < 0.05,
        "air_density ratio: got {:.4}, expected {:.4}",
        ratio,
        golden.ratio
    );
}

// ============================================================
// Test: no_wake_power
// Config: 2 turbines, 1 findex, 8m/s@270°, TI=0.06, no wake
// Both turbines should produce equal power
// ============================================================

#[test]
fn test_no_wake_power() {
    let config = create_api_config_2t(API_WIND_SPEEDS_8.to_vec(), 1.225);

    let mut model = FlorisModel::from_config(config.clone()).expect("Failed to create FlorisModel");

    model.initialize_grid().expect("init grid failed");
    model
        .initialize_flow_field()
        .expect("init flow field failed");

    let grid = model.core().grid.as_ref().unwrap();
    let turbine_grid = grid.as_any().downcast_ref::<TurbineGrid>().unwrap();

    let farm = model.core().farm.clone();
    let mut flow_field = model.core().flow_field.clone();
    let wake = model.core().wake.clone();

    // Run solver but then zero out wake effects for no_wake test
    // Actually, we should use the same approach as Python's run_no_wake
    // which skips the solver entirely and just uses the initial flow field

    // Re-initialize flow field to get initial state (no wake)
    model
        .initialize_flow_field()
        .expect("init flow field failed");
    let flow_field_no_wake = model.core().flow_field.clone();

    let n_findex = flow_field_no_wake.u_sorted.shape()[0];
    let n_turbines = farm.turbines.len();
    let mut powers: Vec<Vec<f64>> = Vec::new();

    for fi in 0..n_findex {
        let mut row: Vec<f64> = Vec::new();
        for ti in 0..n_turbines {
            let mut v_cubed_sum = 0.0;
            let grid_points =
                flow_field_no_wake.u_sorted.shape()[2] * flow_field_no_wake.u_sorted.shape()[3];
            for iy in 0..flow_field_no_wake.u_sorted.shape()[2] {
                for iz in 0..flow_field_no_wake.u_sorted.shape()[3] {
                    v_cubed_sum += flow_field_no_wake.u_sorted[[fi, ti, iy, iz]].powi(3);
                }
            }
            let rotor_avg_velocity = (v_cubed_sum / grid_points as f64).powf(1.0 / 3.0);

            let turbine = &farm.turbines[ti];
            let ref_air_density = turbine
                .turbine_type
                .power_thrust_table
                .ref_air_density
                .unwrap_or(1.225);
            let density_factor = (flow_field_no_wake.air_density / ref_air_density).powf(1.0 / 3.0);
            let rotor_effective_velocity = rotor_avg_velocity * density_factor;

            let mut conditions = TableConditions::default();
            conditions.wind_speed = rotor_effective_velocity;
            let power_kw = turbine
                .turbine_type
                .power_thrust_table
                .cp_ct_table
                .get_cp(&conditions)
                .unwrap_or(0.0);
            let power_val = power_kw * 1000.0;

            row.push(power_val);
        }
        powers.push(row);
    }

    let golden = common::load_api_golden("no_wake_power");
    assert_powers_close(&powers, &golden.turbine_powers, "no_wake_power");

    // Verify both turbines produce equal power (no wake)
    assert_eq!(powers[0].len(), 2, "expected 2 turbines");
    let diff = (powers[0][0] - powers[0][1]).abs();
    assert!(
        diff < 100.0,
        "no_wake: turbines should produce equal power, diff={:.2}",
        diff
    );
}
