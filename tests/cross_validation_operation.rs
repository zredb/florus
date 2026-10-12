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
        solver: SolverConfig {
            // Match the Python golden generator (TURBINE_GRID_RESOLUTION = 2).
            turbine_grid_points: 2,
            ..SolverConfig::default()
        },
        farm: FarmConfig {
            layout_x: vec![0.0],
            layout_y: vec![0.0],
            turbine_type: turbine_types,
        },
        flow_field: FlowFieldConfig {
            wind_speeds: vec![10.0],
            wind_directions: vec![270.0],
            turbulence_intensities: vec![0.06],
            // Match the Python golden generator (generate_golden_values.py
            // create_floris_model): wind_shear=0.12, reference_wind_height=90.
            // Previously shear=0 / ref_height=138.6 — a config mismatch that
            // the previous loose tolerance masked (Bug #35).
            wind_shear: 0.12,
            wind_veer: 0.0,
            air_density: 1.225,
            reference_wind_height: 90.0,
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

/// Run the sequential solver with the requested operation model ACTUALLY
/// enabled, then compute the power THROUGH the operation model via
/// `Farm::power`.
///
/// Previously `create_operation_config(name)` never used the model name (all
/// tests ran the turbine YAML's default model) and power was computed with an
/// inline Cp lookup that bypassed the operation model entirely — so a broken
/// peak-shaving implementation still "passed" against per-model Python
/// goldens under a 25% tolerance (Bug #35).
fn compute_powers_model(
    config: &FlorisConfig,
    operation_model: &str,
    yaw: Option<florus::types::Array2>,
    setpoints: Option<florus::types::Array2>,
    awc: Option<(String, florus::types::Float)>,
) -> Vec<Vec<f64>> {
    let mut model = FlorisModel::from_config(config.clone()).expect("Failed to create FlorisModel");

    // Match the Python golden generator: fmodel.set_operation_model(name)
    model
        .set_operation_model(operation_model)
        .expect("Failed to set operation model");

    if let Some(yaw) = yaw {
        model.set_yaw_angles(yaw).expect("Failed to set yaw angles");
    }

    model.initialize_grid().expect("init grid failed");
    model
        .initialize_flow_field()
        .expect("init flow field failed");

    let grid = model.core().grid.as_ref().unwrap();
    let turbine_grid = grid.as_any().downcast_ref::<TurbineGrid>().unwrap();

    let mut farm = model.core().farm.clone();
    if let Some(sp) = setpoints {
        farm.set_power_setpoints(sp);
    }
    if let Some((mode, amplitude)) = awc {
        // Match the Python generator: fmodel.set(awc_modes=..., awc_amplitudes=...)
        let modes = ndarray::Array2::from_elem((1, 1), mode);
        farm.awc_modes = modes.clone();
        farm.awc_modes_sorted = modes;
        let amps = florus::types::Array2::from_elem((1, 1), amplitude);
        farm.awc_amplitudes = amps.clone();
        farm.awc_amplitudes_sorted = amps;
    }
    let mut flow_field = model.core().flow_field.clone();
    let wake = model.core().wake.clone();

    sequential_solver(&farm, &mut flow_field, turbine_grid, &wake).expect("solver failed");

    // Model-aware power (Python: core power() dispatches through the
    // operation model; the inline Cp lookup bypassed it).
    let powers = farm
        .power(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            Some(&flow_field.turbulence_intensity_field_sorted),
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            turbine_grid.average_method,
        )
        .expect("Farm::power failed");

    powers.outer_iter().map(|row| row.to_vec()).collect()
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

/// Tight assertion for fully-ported operation models.
fn assert_power_exact(test: f64, golden: f64, name: &str) {
    let rel = if golden.abs() > 1e-10 {
        (test - golden).abs() / golden.abs()
    } else {
        (test - golden).abs()
    };
    assert!(
        rel <= 1e-5,
        "{}: power mismatch: got {:.4}, expected {:.4} (rel {:.3e})",
        name,
        test,
        golden,
        rel
    );
}

fn create_operation_config_ti(operation_model: &str, ti: Float) -> FlorisConfig {
    let mut config = create_operation_config(operation_model);
    config.flow_field.turbulence_intensities = vec![ti];
    config
}

// ============================================================
// Operation model tests
// ============================================================

#[test]
fn test_simple_operation_model() {
    let config = create_operation_config("simple");
    let powers = compute_powers_model(&config, "simple", None, None, None);
    let golden = common::load_operation_golden("simple");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "simple");
}

#[test]
fn test_cosine_loss_operation_model() {
    // The Python golden generator sets yaw 10 deg for the cosine-loss case.
    let config = create_operation_config("cosine-loss");
    let yaw = florus::types::Array2::from_shape_vec((1, 1), vec![10.0]).unwrap();
    let powers = compute_powers_model(&config, "cosine-loss", Some(yaw), None, None);
    let golden = common::load_operation_golden("cosine-loss");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "cosine-loss");
}

#[test]
fn test_peak_shaving_operation_model() {
    // TI=0.06 is below the nrel_5MW peak_shaving_TI_threshold (0.1): no
    // shaving is applied, so this checks the passthrough path only.
    let config = create_operation_config("peak-shaving");
    let powers = compute_powers_model(&config, "peak-shaving", None, None, None);
    let golden = common::load_operation_golden("peak-shaving");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "peak-shaving");
}

#[test]
fn test_peak_shaving_high_ti_operation_model() {
    // TI=0.15 exceeds the nrel_5MW peak_shaving_TI_threshold (0.1), so the
    // max-allowable-thrust cap and the ps_ct/ps_ai power fraction are actually
    // exercised — the default TI=0.06 case cannot detect a broken
    // implementation (Bug #35).
    let config = create_operation_config_ti("peak-shaving", 0.15);
    let powers = compute_powers_model(&config, "peak-shaving", None, None, None);
    let golden = common::load_operation_golden("peak_shaving_high_ti");
    assert_power_exact(
        powers[0][0],
        golden.turbine_powers[0][0],
        "peak_shaving_high_ti",
    );
}

#[test]
fn test_mixed_operation_model() {
    // Default yaw/setpoints: every turbine falls in the "neither" bucket and
    // the mixed model reduces to SimpleTurbine.
    let config = create_operation_config("mixed");
    let powers = compute_powers_model(&config, "mixed", None, None, None);
    let golden = common::load_operation_golden("mixed");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "mixed");
}

#[test]
fn test_simple_derating_operation_model() {
    // Default setpoints (1e12): no derating, passthrough path.
    let config = create_operation_config("simple-derating");
    let powers = compute_powers_model(&config, "simple-derating", None, None, None);
    let golden = common::load_operation_golden("simple-derating");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "simple-derating");
}

#[test]
fn test_derating_3mw_operation_model() {
    // 3 MW setpoint at 10 m/s: derating is active, so the linear
    // power-fraction Ct scaling and the min(base, setpoint) power clamp are
    // actually exercised (Bug #29 discriminative case).
    let config = create_operation_config("simple-derating");
    let setpoints = florus::types::Array2::from_shape_vec((1, 1), vec![3e6]).unwrap();
    let powers = compute_powers_model(&config, "simple-derating", None, Some(setpoints), None);
    let golden = common::load_operation_golden("derating_3MW");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "derating_3MW");
}

#[test]
fn test_awc_operation_model() {
    // Baseline mode (no excitation): the AWC model reduces to the simple
    // power/Ct lookup. Python golden = simple golden in this case.
    let config = create_operation_config("awc");
    let powers = compute_powers_model(&config, "awc", None, None, None);
    let golden = common::load_operation_golden("awc");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "awc");
}

#[test]
fn test_awc_helix_operation_model() {
    // Active helix excitation (2.5 deg): the AWC_model correction
    // `base * (1 - (b + c*base) * amp^a)` is actually applied to both power
    // and Ct — the baseline case cannot discriminate a placeholder (Bug #37).
    let config = create_operation_config("awc");
    let powers = compute_powers_model(&config, "awc", None, None, Some(("helix".to_string(), 2.5)));
    let golden = common::load_operation_golden("awc_helix");
    assert_power_exact(powers[0][0], golden.turbine_powers[0][0], "awc_helix");
}

#[test]
fn test_unified_momentum_operation_model() {
    // At zero yaw the UMM velocity-correction ratio is 1, so the model
    // reduces to the plain table lookup at the rotor-effective velocity.
    let config = create_operation_config("unified-momentum");
    let powers = compute_powers_model(&config, "unified-momentum", None, None, None);
    let golden = common::load_operation_golden("unified-momentum");
    assert_power_exact(
        powers[0][0],
        golden.turbine_powers[0][0],
        "unified-momentum",
    );
}

#[test]
fn test_unified_momentum_yawed_operation_model() {
    // 20 deg yaw: the Heck fixed-point iteration and the yawed rotor
    // velocity correction are actually exercised — the zero-yaw case cannot
    // discriminate a placeholder (Bug #37).
    let config = create_operation_config("unified-momentum");
    let yaw = florus::types::Array2::from_shape_vec((1, 1), vec![20.0]).unwrap();
    let powers = compute_powers_model(&config, "unified-momentum", Some(yaw), None, None);
    let golden = common::load_operation_golden("unified_momentum_yawed");
    assert_power_exact(
        powers[0][0],
        golden.turbine_powers[0][0],
        "unified_momentum_yawed",
    );
}

#[test]
fn test_controller_dependent_operation_model() {
    // Full TUM model port: BEM closed forms + control trajectory + cubic
    // not-a-knot Cp/Ct surface interpolation, validated against Python at
    // tight tolerance (Bug #37).
    let config = create_operation_config("controller-dependent");
    let powers = compute_powers_model(&config, "controller-dependent", None, None, None);
    let golden = common::load_operation_golden("controller-dependent");
    assert_power_exact(
        powers[0][0],
        golden.turbine_powers[0][0],
        "controller-dependent",
    );
}

#[test]
fn test_controller_dependent_yawed_operation_model() {
    // 15 deg yaw: the MU-dependent BEM terms, the cosd(yaw) control-trajectory
    // terms and the yawed/unyawed Cp ratio are actually exercised (Bug #37).
    let config = create_operation_config("controller-dependent");
    let yaw = florus::types::Array2::from_shape_vec((1, 1), vec![15.0]).unwrap();
    let powers = compute_powers_model(&config, "controller-dependent", Some(yaw), None, None);
    let golden = common::load_operation_golden("controller_dependent_yawed");
    assert_power_exact(
        powers[0][0],
        golden.turbine_powers[0][0],
        "controller_dependent_yawed",
    );
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
    assert_eq!(
        golden.windrose.wind_directions,
        vec![0.0, 90.0, 180.0, 270.0]
    );
    assert_eq!(
        golden.windrose.wind_speeds,
        vec![6.0, 8.0, 10.0, 12.0, 14.0]
    );
    assert_eq!(golden.windrose.freq_table_shape, vec![4, 5]);
    assert!((golden.windrose.freq_table_sum - 1.0).abs() < 1e-10);
}
