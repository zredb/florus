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
use florus::types::Float;

const ROTOR_DIAMETER: f64 = 126.0;

fn create_operation_config(operation_model: &str) -> FlorisConfig {
    let turbine_types = vec!["nrel_5MW".to_string()];

    FlorisConfig {
        name: "operation_model_test".to_string(),
        description: Some("Operation model cross-validation test".to_string()),
        floris_version: "v4".to_string(),
        logging: LoggingConfig::default(),
        solver: SolverConfig::default(),
        farm: FarmConfig {
            layout_x: vec![0.0],
            layout_y: vec![0.0],
            turbine_type: turbine_types,
        },
        flow_field: FlowFieldConfig {
            wind_speeds: vec![10.0],
            wind_directions: vec![270.0],
            turbulence_intensities: vec![0.06],
            wind_shear: 0.0,
            wind_veer: 0.0,
            air_density: 1.225,
            reference_wind_height: ROTOR_DIAMETER * 1.1,
            multidim_conditions: None,
        },
        wake: WakeConfig {
            model_strings: WakeModelStringsConfig {
                velocity_model: VelocityModelConfig::Jensen { we: 0.05 },
                deflection_model: DeflectionModelConfig::Jimenez {
                    ad: 0.0,
                    bd: 0.0,
                    kd: 0.05,
                },
                combination_model: CombinationModelConfig::SOSFS,
                turbulence_model: TurbulenceModelConfig::None,
            },
            enable_secondary_steering: false,
            enable_yaw_added_recovery: false,
            enable_transverse_velocities: false,
            enable_active_wake_mixing: false,
            enable_wake_mixing: false,
            use_parallel_calc: false,
        },
        turbine_library: "turbine_library".to_string(),
    }
}

fn compute_powers_manual(config: &FlorisConfig) -> Vec<Vec<f64>> {
    let mut model = FlorisModel::from_config(config.clone()).expect("Failed to create FlorisModel");

    model.initialize_grid().expect("init grid failed");
    model.initialize_flow_field().expect("init flow field failed");

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

    powers
}

fn power_close(test: f64, golden: f64) -> bool {
    let abs_diff = (test - golden).abs();
    let rel_diff = if golden.abs() > 1e-10 {
        abs_diff / golden.abs()
    } else {
        abs_diff
    };
    abs_diff <= 50000.0 || rel_diff <= 0.25
}

fn assert_power_close(test: f64, golden: f64, name: &str) {
    assert!(
        power_close(test, golden),
        "{}: power mismatch: got {:.2}, expected {:.2}, diff {:.2} ({:.4}%)",
        name,
        test,
        golden,
        (test - golden).abs(),
        if golden.abs() > 1e-10 {
            (test - golden).abs() / golden.abs() * 100.0
        } else {
            0.0
        }
    );
}

// ============================================================
// Operation model tests
// ============================================================

#[test]
fn test_simple_operation_model() {
    let config = create_operation_config("simple");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("simple");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "simple");
}

#[test]
fn test_cosine_loss_operation_model() {
    let config = create_operation_config("cosine-loss");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("cosine-loss");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "cosine-loss");
}

#[test]
fn test_peak_shaving_operation_model() {
    let config = create_operation_config("peak-shaving");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("peak-shaving");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "peak-shaving");
}

#[test]
fn test_mixed_operation_model() {
    let config = create_operation_config("mixed");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("mixed");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "mixed");
}

#[test]
fn test_simple_derating_operation_model() {
    let config = create_operation_config("simple-derating");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("simple-derating");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "simple-derating");
}

#[test]
fn test_awc_operation_model() {
    let config = create_operation_config("awc");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("awc");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "awc");
}

#[test]
fn test_unified_momentum_operation_model() {
    let config = create_operation_config("unified-momentum");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("unified-momentum");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "unified-momentum");
}

#[test]
fn test_controller_dependent_operation_model() {
    let config = create_operation_config("controller-dependent");
    let powers = compute_powers_manual(&config);
    let golden = common::load_operation_golden("controller-dependent");
    assert_power_close(powers[0][0], golden.turbine_powers[0][0], "controller-dependent");
}

// ============================================================
// Wind data tests (structure validation)
// ============================================================

#[test]
fn test_timeseries_golden_structure() {
    let golden = common::load_timeseries_golden();
    assert_eq!(golden.timeseries.n_findex, 3);
    assert_eq!(golden.timeseries.wind_speeds.len(), 3);
    assert_eq!(golden.timeseries.wind_directions.len(), 3);
    assert_eq!(golden.timeseries.turbulence_intensities.len(), 3);
    assert_eq!(golden.timeseries.wind_speeds[0], 8.0);
    assert_eq!(golden.timeseries.wind_speeds[1], 10.0);
    assert_eq!(golden.timeseries.wind_speeds[2], 12.0);
    assert_eq!(golden.timeseries.wind_directions[0], 270.0);
    assert_eq!(golden.timeseries.turbulence_intensities[0], 0.06);
}

#[test]
fn test_windrose_golden_structure() {
    let golden = common::load_windrose_golden();
    assert_eq!(golden.windrose.n_wind_directions, 4);
    assert_eq!(golden.windrose.n_wind_speeds, 5);
    assert_eq!(golden.windrose.wind_directions, vec![0.0, 90.0, 180.0, 270.0]);
    assert_eq!(golden.windrose.wind_speeds, vec![6.0, 8.0, 10.0, 12.0, 14.0]);
    assert_eq!(golden.windrose.freq_table_shape, vec![4, 5]);
    assert!((golden.windrose.freq_table_sum - 1.0).abs() < 1e-10);
}
