//! Cross-validation of the post-solve query/sample APIs against Python
//! FLORIS goldens (tests/golden_values/query_apis/*.json):
//!
//! - `get_turbine_powers` / `get_turbine_thrust_coefficients` /
//!   `get_turbine_ais` routed through the operation model over the SOLVED
//!   rotor velocities (previously: inline Cp estimate, freestream Ct table
//!   lookup, heuristic AI).
//! - `sample_flow_at_points` / `sample_ti_at_points`: wake solve ON the
//!   sample points for ALL findices (previously: first-findex
//!   nearest-neighbor on the turbine grid and a constant-ambient-TI stub).
//! - The floating-tilt Cp/Ct correction (`correct_cp_ct_for_tilt` +
//!   `floating_tilt_table`) through the cosine-loss power/thrust paths and
//!   the solver's Ct/AI chain (nrel_5MW_floating test turbine).
//! - The stale-state guard on `get_turbine_powers` (Python raises
//!   RuntimeError from `_get_turbine_powers` unless the model has run).

mod common;

use common::load_query_apis_golden;
use florus::floris_config::{
    CombinationModelConfig, DeflectionModelConfig, FarmConfig, FlorisConfig, FlowFieldConfig,
    LoggingConfig, SolverConfig, TurbulenceModelConfig, VelocityModelConfig, WakeConfig,
    WakeModelStringsConfig,
};
use florus::floris_model::FlorisModel;
use florus::types::Array2;

const WIND_DIRECTIONS: [f64; 16] = [
    270.0, 270.0, 270.0, 270.0, 360.0, 360.0, 360.0, 360.0, 285.0, 285.0, 285.0, 285.0, 315.0,
    315.0, 315.0, 315.0,
];

const WIND_SPEEDS: [f64; 16] = [
    8.0, 9.0, 10.0, 11.0, 8.0, 9.0, 10.0, 11.0, 8.0, 9.0, 10.0, 11.0, 8.0, 9.0, 10.0, 11.0,
];

const TURBULENCE_INTENSITIES: [f64; 16] = [
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
];

const ROTOR_DIAMETER: f64 = 126.0;
const X_COORDS: [f64; 3] = [0.0, 5.0 * ROTOR_DIAMETER, 10.0 * ROTOR_DIAMETER];
const Y_COORDS: [f64; 3] = [0.0, 0.0, 0.0];

fn create_query_apis_config(turbine_type: &str) -> FlorisConfig {
    FlorisConfig {
        name: "query_apis_cross_validation".to_string(),
        description: Some("Query APIs cross-validation test config".to_string()),
        floris_version: "v4".to_string(),
        logging: LoggingConfig::default(),
        solver: SolverConfig {
            turbine_grid_points: 2,
            ..SolverConfig::default()
        },
        farm: FarmConfig {
            layout_x: X_COORDS.to_vec(),
            layout_y: Y_COORDS.to_vec(),
            turbine_type: vec![turbine_type.to_string(); 3],
        },
        flow_field: FlowFieldConfig {
            wind_speeds: WIND_SPEEDS.to_vec(),
            wind_directions: WIND_DIRECTIONS.to_vec(),
            turbulence_intensities: TURBULENCE_INTENSITIES.to_vec(),
            wind_shear: 0.12,
            wind_veer: 0.0,
            air_density: 1.225,
            reference_wind_height: 90.0,
            multidim_conditions: None,
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

/// Yaw angles matching the Python golden: upstream turbine yawed 5 degrees.
fn yaw_angles() -> Array2 {
    let mut yaw = Array2::zeros((16, 3));
    for fi in 0..16 {
        yaw[[fi, 0]] = 5.0;
    }
    yaw
}

fn assert_2d_close(
    test: &florus::types::Array2,
    baseline: &[Vec<f64>],
    description: &str,
    rel_tol: f64,
    abs_tol: f64,
) {
    assert_eq!(
        test.shape(),
        [
            baseline.len(),
            baseline.first().map(|r| r.len()).unwrap_or(0)
        ],
        "{}: shape mismatch: got {:?}, expected [{}, {}]",
        description,
        test.shape(),
        baseline.len(),
        baseline.first().map(|r| r.len()).unwrap_or(0),
    );
    for (fi, row) in baseline.iter().enumerate() {
        for (ti, &b) in row.iter().enumerate() {
            let t = test[[fi, ti]];
            let diff = (t - b).abs();
            let rel = if b.abs() > 1e-10 {
                diff / b.abs()
            } else {
                diff
            };
            assert!(
                rel <= rel_tol || diff <= abs_tol,
                "{}[{}][{}]: got {}, expected {}, rel diff {} (tol {} / {})",
                description,
                fi,
                ti,
                t,
                b,
                rel,
                rel_tol,
                abs_tol
            );
        }
    }
}

fn assert_nested_close(
    test: &[Vec<f64>],
    baseline: &[Vec<f64>],
    description: &str,
    rel_tol: f64,
    abs_tol: f64,
) {
    for (fi, (t_row, b_row)) in test.iter().zip(baseline.iter()).enumerate() {
        for (pi, (&t, &b)) in t_row.iter().zip(b_row.iter()).enumerate() {
            let diff = (t - b).abs();
            let rel = if b.abs() > 1e-10 {
                diff / b.abs()
            } else {
                diff
            };
            assert!(
                rel <= rel_tol || diff <= abs_tol,
                "{}[{}][{}]: got {}, expected {}, rel diff {} (tol {} / {})",
                description,
                fi,
                pi,
                t,
                b,
                rel,
                rel_tol,
                abs_tol
            );
        }
    }
}

/// Python golden `gauss_yawed_query_apis.json`: operation-model power/Ct/AI
/// with yaw, and the on-points wake solve for velocity and TI across ALL 16
/// findices.
#[test]
fn test_query_apis_gauss_yawed_golden() {
    let mut model = FlorisModel::from_config(create_query_apis_config("nrel_5MW")).expect("model");
    model.set_yaw_angles(yaw_angles()).expect("set_yaw_angles");
    model.run().expect("run");

    let golden = load_query_apis_golden("gauss_yawed_query_apis");

    // Powers through the operation-model route (fixed-body, yawed).
    let powers = model.get_turbine_powers().expect("get_turbine_powers");
    assert_2d_close(&powers, &golden.turbine_powers, "powers", 1e-3, 500.0);

    // Ct from the solved rotor velocities (previously the freestream table).
    let cts = model.get_turbine_thrust_coefficients();
    assert_2d_close(
        &cts,
        golden.turbine_thrust_coefficients.as_ref().unwrap(),
        "thrust coefficients",
        1e-4,
        1e-4,
    );

    // Axial induction through the operation model (previously a heuristic).
    let ais = model.get_turbine_ais().expect("get_turbine_ais");
    assert_2d_close(
        &ais,
        golden.turbine_axial_inductions.as_ref().unwrap(),
        "axial inductions",
        1e-4,
        1e-4,
    );

    // Sampling: wake solved ON the points for ALL findices.
    let px = golden.sample_points_x.as_ref().unwrap();
    let py = golden.sample_points_y.as_ref().unwrap();
    let pz = golden.sample_points_z.as_ref().unwrap();
    let sampled = model
        .sample_flow_at_points(px, py, pz)
        .expect("sample_flow_at_points");
    let sampled_rows: Vec<Vec<f64>> = sampled.outer_iter().map(|r| r.to_vec()).collect();
    assert_nested_close(
        &sampled_rows,
        golden.sample_velocities.as_ref().unwrap(),
        "sample velocities",
        1e-3,
        0.05,
    );

    // Sampled TI carries the wake-added turbulence (the pre-fix stub returned
    // the constant ambient TI).
    let sampled_ti = model
        .sample_ti_at_points(px, py, pz)
        .expect("sample_ti_at_points");
    let ti_rows: Vec<Vec<f64>> = sampled_ti.outer_iter().map(|r| r.to_vec()).collect();
    assert_nested_close(
        &ti_rows,
        golden.sample_tis.as_ref().unwrap(),
        "sample TI",
        1e-3,
        1e-3,
    );
}

/// Python golden `floating_tilt.json` (nrel_5MW_floating: same curves as
/// nrel_5MW plus a wind-speed-dependent platform tilt with
/// `correct_cp_ct_for_tilt: true`). The tilt correction enters the
/// cosine-loss power/thrust paths and the solver's Ct/AI chain; a missing
/// wiring reproduces the fixed-body values, which differ by ~0.25% on the
/// upstream turbine (well above these tolerances).
#[test]
fn test_query_apis_floating_tilt_golden() {
    let mut model =
        FlorisModel::from_config(create_query_apis_config("nrel_5MW_floating")).expect("model");
    model.run().expect("run");

    let golden = load_query_apis_golden("floating_tilt");

    let powers = model.get_turbine_powers().expect("get_turbine_powers");
    assert_2d_close(
        &powers,
        &golden.turbine_powers,
        "floating powers",
        1e-3,
        500.0,
    );

    let cts = model.get_turbine_thrust_coefficients();
    assert_2d_close(
        &cts,
        golden.turbine_thrust_coefficients.as_ref().unwrap(),
        "floating thrust coefficients",
        1e-4,
        1e-4,
    );

    let ais = model.get_turbine_ais().expect("get_turbine_ais");
    assert_2d_close(
        &ais,
        golden.turbine_axial_inductions.as_ref().unwrap(),
        "floating axial inductions",
        1e-4,
        1e-4,
    );
}

/// `get_turbine_powers` must refuse stale state, mirroring Python's
/// RuntimeError from `_get_turbine_powers` (Python floris_model.py:545-548).
/// This is what prevents silently reading the pre-`set_operation` wake
/// field after `set_operation(disable_turbines=...)`.
#[test]
fn test_query_apis_state_check() {
    let mut model = FlorisModel::from_config(create_query_apis_config("nrel_5MW")).expect("model");

    let err = model
        .get_turbine_powers()
        .expect_err("query before run() must fail");
    assert!(
        err.to_string().contains("without first running"),
        "unexpected error message: {err:#}"
    );

    // After run() the query succeeds; after a state-invalidating
    // set_operation it must fail again until run() is called.
    model.run().expect("run");
    model
        .get_turbine_powers()
        .expect("query after run() must succeed");

    let disabled = ndarray::Array2::from_elem((16, 3), false);
    model
        .set_operation(None, None, None, None, None, Some(disabled))
        .expect("set_operation");
    let err = model
        .get_turbine_powers()
        .expect_err("query after set_operation without run() must fail");
    assert!(
        err.to_string().contains("without first running"),
        "unexpected error message: {err:#}"
    );

    model.run().expect("run after set_operation");
    model
        .get_turbine_powers()
        .expect("query after re-run must succeed");
}

/// The sampling APIs must reject models whose full-flow solver variant is
/// not available, matching Python's explicit refusal for turbopark.
#[test]
fn test_query_apis_sampling_rejects_turbopark() {
    let mut config = create_query_apis_config("nrel_5MW");
    config.wake.model_strings.velocity_model = VelocityModelConfig::Turbopark {
        kstar: 0.05,
        cstar: 1.5,
    };
    // Python disables the GCH extras for turbopark (generate_golden_values.py
    // lines 228-231); run() rejects them otherwise.
    config.wake.enable_secondary_steering = false;
    config.wake.enable_yaw_added_recovery = false;
    config.wake.enable_transverse_velocities = false;
    let mut model = FlorisModel::from_config(config).expect("model");
    model.run().expect("run");

    let err = model
        .sample_flow_at_points(&[630.0], &[0.0], &[90.0])
        .expect_err("turbopark sampling must fail");
    assert!(
        err.to_string().contains("turbopark"),
        "unexpected error message: {err:#}"
    );
}
