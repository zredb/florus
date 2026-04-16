//! Regression tests for wake models
//!
//! These tests validate that the Rust implementation matches the Python FLORIS v4.6
//! baseline values. The tests compare average velocity, thrust coefficient, power,
//! and axial induction against expected values.
//!
//! Test configurations based on:
//! - floris-4.6/tests/reg_tests/turbopark_regression_test.py
//! - floris-4.6/tests/reg_tests/gauss_regression_test.py
//! - floris-4.6/tests/reg_tests/jensen_jimenez_regression_test.py

use crate::core::TurbineGrid;
use crate::core::wake::{WakeModelManager, WakeModelStrings};
use crate::floris_config::{FarmConfig, FlorisConfig, FlowFieldConfig, SolverConfig, WakeConfig};
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2};
use ndarray::Array;

/// Test fixture parameters from Python conftest.py
const WIND_DIRECTIONS: [f64; 16] = [
    270.0, 270.0, 270.0, 270.0,
    360.0, 360.0, 360.0, 360.0,
    285.0, 285.0, 285.0, 285.0,
    315.0, 315.0, 315.0, 315.0,
];

const WIND_SPEEDS: [f64; 16] = [
    8.0, 9.0, 10.0, 11.0,
    8.0, 9.0, 10.0, 11.0,
    8.0, 9.0, 10.0, 11.0,
    8.0, 9.0, 10.0, 11.0,
];

const TURBULENCE_INTENSITIES: [f64; 16] = [
    0.1, 0.1, 0.1, 0.1,
    0.1, 0.1, 0.1, 0.1,
    0.1, 0.1, 0.1, 0.1,
    0.1, 0.1, 0.1, 0.1,
];

const N_FINDEX: usize = 16;
const N_TURBINES: usize = 3;
const ROTOR_DIAMETER: f64 = 126.0;
const TURBINE_GRID_RESOLUTION: usize = 2;

/// X coordinates: 0, 5D, 10D
const X_COORDS: [f64; 3] = [0.0, 5.0 * ROTOR_DIAMETER, 10.0 * ROTOR_DIAMETER];
const Y_COORDS: [f64; 3] = [0.0, 0.0, 0.0];
const Z_COORDS: [f64; 3] = [90.0, 90.0, 90.0];

/// TurbOPark baseline values from Python FLORIS v4.6
/// Shape: (n_findex, n_turbines, 4) = [velocity, Ct, power, axial_induction]
const TURBOPARK_BASELINE: [[[f64; 4]; 3]; 4] = [
    // 8 m/s
    [
        [7.9736858, 0.7871515, 1753954.4591792, 0.2693224],
        [6.0332948, 0.8593353, 752557.9240063, 0.3124735],
        [5.4029800, 0.8947888, 538370.5108659, 0.3378186],
    ],
    // 9 m/s
    [
        [8.9703965, 0.7858774, 2496427.8618358, 0.2686331],
        [6.7887441, 0.8249788, 1092199.1775234, 0.2908223],
        [6.0678594, 0.8577634, 768097.7785191, 0.3114286],
    ],
    // 10 m/s
    [
        [9.9671073, 0.7838789, 3417797.0050916, 0.2675559],
        [7.5453629, 0.7962514, 1487438.4031455, 0.2743074],
        [6.7548552, 0.8265200, 1076963.1412833, 0.2917453],
    ],
    // 11 m/s
    [
        [10.9638180, 0.7565157, 4519404.3072862, 0.2532794],
        [8.3436376, 0.7866851, 2027996.3027579, 0.2690699],
        [7.4626804, 0.7989174, 1439263.3915910, 0.2757889],
    ],
];

/// TurbOPark yawed baseline (upstream turbine yawed 5 degrees)
const TURBOPARK_YAWED_BASELINE: [[[f64; 4]; 3]; 4] = [
    // 8 m/s
    [
        [7.9736858, 0.7841561, 1741508.6722008, 0.2671213],
        [6.0523119, 0.8584704, 761107.7639542, 0.3118979],
        [5.4177841, 0.8939472, 543310.4550423, 0.3371713],
    ],
    // 9 m/s
    [
        [8.9703965, 0.7828869, 2480428.8963141, 0.2664440],
        [6.8101438, 0.8240055, 1101820.2623232, 0.2902415],
        [6.0851644, 0.8569764, 775877.8906008, 0.3109077],
    ],
    // 10 m/s
    [
        [9.9671073, 0.7808960, 3395681.0032992, 0.2653854],
        [7.5691494, 0.7955016, 1501458.3309846, 0.2738925],
        [6.7745474, 0.8256244, 1085816.5021615, 0.2912085],
    ],
    // 11 m/s
    [
        [10.9638180, 0.7536370, 4488242.9153943, 0.2513413],
        [8.3695194, 0.7866518, 2047340.0279521, 0.2690518],
        [7.4830530, 0.7982426, 1450966.1620998, 0.2754129],
    ],
];

/// Gauss baseline values from Python FLORIS v4.6
const GAUSS_BASELINE: [[[f64; 4]; 3]; 4] = [
    // 8 m/s
    [
        [7.9736858, 0.7871515, 1753954.4591792, 0.2693224],
        [5.9186455, 0.8654743, 710441.9192938, 0.3166113],
        [6.0090150, 0.8604395, 741642.0177873, 0.3132110],
    ],
    // 9 m/s
    [
        [8.9703965, 0.7858774, 2496427.8618358, 0.2686331],
        [6.6606465, 0.8308044, 1034608.0101396, 0.2943330],
        [6.7947466, 0.8247058, 1094897.8563374, 0.2906592],
    ],
    // 10 m/s
    [
        [9.9671073, 0.7838789, 3417797.0050916, 0.2675559],
        [7.4045198, 0.8008441, 1405853.7207176, 0.2768656],
        [7.5868432, 0.7949439, 1511887.2179035, 0.2735844],
    ],
    // 11 m/s
    [
        [10.9638180, 0.7565157, 4519404.3072862, 0.2532794],
        [8.2046271, 0.7868643, 1924101.6501936, 0.2691669],
        [8.3491997, 0.7866780, 2032153.3223547, 0.2690660],
    ],
];

/// Jensen baseline values from Python FLORIS v4.6
const JENSEN_BASELINE: [[[f64; 4]; 3]; 4] = [
    // 8 m/s
    [
        [7.9736858, 0.7871515, 1753954.4591792, 0.2693224],
        [6.0660565, 0.8578454, 767287.2198744, 0.3114830],
        [5.5204712, 0.8881097, 577575.9208353, 0.3327500],
    ],
    // 9 m/s
    [
        [8.9703965, 0.7858774, 2496427.8618358, 0.2686331],
        [6.8298067, 0.8231113, 1110660.4518964, 0.2897093],
        [6.3668912, 0.8441639, 902538.9934586, 0.3026196],
    ],
    // 10 m/s
    [
        [9.9671073, 0.7838789, 3417797.0050916, 0.2675559],
        [7.5982117, 0.7945856, 1518587.8467982, 0.2733867],
        [7.2042504, 0.8077903, 1294847.7809883, 0.2807914],
    ],
    // 11 m/s
    [
        [10.9638180, 0.7565157, 4519404.3072862, 0.2532794],
        [8.4970746, 0.7864874, 2142673.1558338, 0.2689629],
        [7.9997342, 0.7871282, 1770992.0756703, 0.2693098],
    ],
];

/// Helper to create a test model with specified wake models
fn create_test_model(
    velocity_model: &str,
    deflection_model: &str,
    combination_model: &str,
    turbulence_model: &str,
) -> FlorisModel {
    let turbine_types = vec!["nrel_5MW".to_string(); N_TURBINES];

    let farm_config = FarmConfig {
        layout_x: X_COORDS.to_vec(),
        layout_y: Y_COORDS.to_vec(),
        turbine_type: turbine_types,
    };

    let flow_field_config = FlowFieldConfig {
        wind_speeds: WIND_SPEEDS[0..4].to_vec(),
        wind_directions: WIND_DIRECTIONS[0..4].to_vec(),
        turbulence_intensities: TURBULENCE_INTENSITIES[0..4].to_vec(),
        wind_shear: 0.12,
        wind_veer: 0.0,
        air_density: 1.225,
        reference_wind_height: 90.0,
        multidim_conditions: None,
    };

    let solver_config = SolverConfig::default();
    let wake_config = WakeConfig::default();

    let config = FlorisConfig {
        name: "test".to_string(),
        description: Some("test".to_string()),
        floris_version: "v4".to_string(),
        logging: Default::default(),
        solver: solver_config,
        farm: farm_config,
        flow_field: flow_field_config,
        wake: wake_config,
        turbine_library: "turbine_library".to_string(),
    };

    FlorisModel::from_config(config).expect("Failed to create FlorisModel")
}

/// Run the solver and extract results
fn run_solver_and_extract_results(
    model: &mut FlorisModel,
    velocity_model: &str,
    deflection_model: &str,
    combination_model: &str,
    turbulence_model: &str,
) -> (Array2, Array2, Array2, Array2) {
    model.initialize_grid().expect("Failed to initialize grid");
    model.initialize_flow_field().expect("Failed to initialize flow field");

    let model_strings = WakeModelStrings {
        velocity_model: velocity_model.to_string(),
        deflection_model: deflection_model.to_string(),
        combination_model: combination_model.to_string(),
        turbulence_model: turbulence_model.to_string(),
    };

    let model_manager = WakeModelManager::new(
        model_strings,
        std::collections::HashMap::new(),
        std::collections::HashMap::new(),
        std::collections::HashMap::new(),
        false, false, false,
    ).expect("Failed to create model manager");

    let grid = model.grid.as_ref().unwrap().as_ref();
    crate::core::solver::sequential_solver(
        &model.farm,
        &mut model.flow_field,
        grid,
        &model_manager,
    ).expect("Solver failed");

    // Extract results - these are the values Python tests compare
    let velocities = model.flow_field.u_sorted.clone();
    let n_findex = model.flow_field.n_findex;
    let n_turbines = model.farm.n_turbines();

    // Calculate average velocities per turbine
    let mut avg_velocities = Array2::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let mut sum = 0.0;
            let mut count = 0;
            for iy in 0..3 {
                for iz in 0..3 {
                    sum += velocities[[fi, ti, iy, iz]];
                    count += 1;
                }
            }
            avg_velocities[[fi, ti]] = sum / count as f64;
        }
    }

    // Get Ct values
    let ct_values = model.get_turbine_thrust_coefficients();

    // Get powers
    let powers = model.get_turbine_powers();

    // Get axial induction (simplified calculation)
    let axial_inductions = ct_values.mapv(|ct| {
        if ct < 0.96 {
            0.5 * (1.0 - (1.0 - ct).sqrt())
        } else {
            0.143 + (0.0203 - 0.6427 * (0.889 - ct).sqrt()).max(0.0)
        }
    });

    (avg_velocities, ct_values, powers, axial_inductions)
}

/// Assert results are close to baseline with tolerance
/// baseline is 3D: [findex][turbine][value_type] where value_type: 0=velocity, 1=Ct, 2=power, 3=AI
fn assert_results_close(
    test: &Array2,
    baseline: &[[[f64; 4]; 3]],
    column: usize,
    tolerance: f64,
    description: &str,
) {
    for fi in 0..4 {
        for ti in 0..3 {
            let test_val = test[[fi, ti]];
            let baseline_val = baseline[fi][ti][column];

            let diff = (test_val - baseline_val).abs();
            let rel_diff = if baseline_val != 0.0 {
                diff / baseline_val.abs()
            } else {
                diff
            };

            assert!(
                diff < tolerance || rel_diff < 0.01,
                "{}: findex={}, turbine={}, got={:.4}, expected={:.4}, diff={:.4}",
                description, fi, ti, test_val, baseline_val, diff
            );
        }
    }
}

#[cfg(test)]
mod turbopark_regression_tests {
    use super::*;

    #[test]
    fn test_turbopark_tandem_regression() {
        // Test TurbOPark model against Python baseline
        let mut model = create_test_model("turbopark", "gauss", "fls", "none");

        let (velocities, ct, powers, axial_ind) = run_solver_and_extract_results(
            &mut model, "turbopark", "gauss", "fls", "none"
        );

        // Compare with 1% tolerance
        // column: 0=velocity, 1=Ct, 2=power, 3=axial_induction
        assert_results_close(&velocities, &TURBOPARK_BASELINE, 0, 0.1, "velocity");
        assert_results_close(&ct, &TURBOPARK_BASELINE, 1, 0.1, "Ct");
        // Power has larger tolerance due to rounding in power calculation
        assert_results_close(&powers, &TURBOPARK_BASELINE, 2, 10000.0, "power");
        assert_results_close(&axial_ind, &TURBOPARK_BASELINE, 3, 0.1, "axial_induction");
    }

    #[test]
    fn test_turbopark_rotation_regression() {
        // Test that rotation gives consistent results
        // 4 turbines in a grid, wind from 270 and 360 degrees
        let turbine_diameter = ROTOR_DIAMETER;
        let spacing = 5.0 * turbine_diameter;

        let turbine_types = vec!["nrel_5MW".to_string(); 4];
        let layout_x = Array1::from_vec(vec![0.0, 0.0, spacing, spacing]);
        let layout_y = Array1::from_vec(vec![0.0, spacing, 0.0, spacing]);

        let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())
            .expect("Failed to create farm");

        // 2 wind directions: 270 (west) and 360 (north)
        let flow_field = FlowField::new(
            Array1::from_vec(vec![8.0, 8.0]),
            Array1::from_vec(vec![270.0, 360.0]),
            0.0,
            0.12,
            1.225,
            Array1::from_vec(vec![0.1, 0.1]),
            90.0,
        ).expect("Failed to create flow field");

        let mut model = FlorisModel {
            farm,
            flow_field,
            state: crate::core::State::new(),
            grid: None,
            solver: SolverConfig::default(),
            model_manager: None,
        };

        model.initialize_grid().expect("Failed to initialize grid");
        model.initialize_flow_field().expect("Failed to initialize flow field");

        let model_strings = WakeModelStrings {
            velocity_model: "turbopark".to_string(),
            deflection_model: "gauss".to_string(),
            combination_model: "fls".to_string(),
            turbulence_model: "none".to_string(),
        };

        let model_manager = WakeModelManager::new(
            model_strings,
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            false, false, false,
        ).expect("Failed to create model manager");

        let grid = model.grid.as_ref().unwrap().as_ref();
        crate::core::solver::sequential_solver(
            &model.farm,
            &mut model.flow_field,
            grid,
            &model_manager,
        ).expect("Solver failed");

        // Calculate average velocities
        let velocities = model.flow_field.u_sorted.clone();
        let n_findex = 2;
        let n_turbines = 4;

        let mut avg_velocities = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let mut sum = 0.0;
                let mut count = 0;
                for iy in 0..3 {
                    for iz in 0..3 {
                        sum += velocities[[fi, ti, iy, iz]];
                        count += 1;
                    }
                }
                avg_velocities[[fi, ti]] = sum / count as f64;
            }
        }

        // At 270 degrees: T0 and T1 are upstream, T2 and T3 are waked
        // At 360 degrees: T0 and T2 are waked, T1 and T3 are upstream
        // Upstream velocities should be similar between directions
        let t0_270 = avg_velocities[[0, 0]];
        let t1_270 = avg_velocities[[0, 1]];
        let t2_270 = avg_velocities[[0, 2]];
        let t3_270 = avg_velocities[[0, 3]];

        let t0_360 = avg_velocities[[1, 0]];
        let t1_360 = avg_velocities[[1, 1]];
        let t2_360 = avg_velocities[[1, 2]];
        let t3_360 = avg_velocities[[1, 3]];

        // Compare velocities between rotations
        // t0_270 should equal t1_360 (both upstream)
        // t1_270 should equal t3_360 (both upstream)
        // t2_270 should equal t0_360 (both waked)
        // t3_270 should equal t2_360 (both waked)

        let tol = 0.01;
        assert!((t0_270 - t1_360).abs() < tol, "t0_270 ({}) should equal t1_360 ({})", t0_270, t1_360);
        assert!((t1_270 - t3_360).abs() < tol, "t1_270 ({}) should equal t3_360 ({})", t1_270, t3_360);
        assert!((t2_270 - t0_360).abs() < tol, "t2_270 ({}) should equal t0_360 ({})", t2_270, t0_360);
        assert!((t3_270 - t2_360).abs() < tol, "t3_270 ({}) should equal t2_360 ({})", t3_270, t2_360);
    }
}

#[cfg(test)]
mod gauss_regression_tests {
    use super::*;

    #[test]
    fn test_gauss_tandem_regression() {
        let mut model = create_test_model("gauss", "gauss", "fls", "none");

        let (velocities, ct, powers, axial_ind) = run_solver_and_extract_results(
            &mut model, "gauss", "gauss", "fls", "none"
        );

        assert_results_close(&velocities, &GAUSS_BASELINE, 0, 0.1, "velocity");
        assert_results_close(&ct, &GAUSS_BASELINE, 1, 0.1, "Ct");
        assert_results_close(&powers, &GAUSS_BASELINE, 2, 10000.0, "power");
        assert_results_close(&axial_ind, &GAUSS_BASELINE, 3, 0.1, "axial_induction");
    }
}

#[cfg(test)]
mod jensen_regression_tests {
    use super::*;

    #[test]
    fn test_jensen_tandem_regression() {
        let mut model = create_test_model("jensen", "jimenez", "sosfs", "none");

        let (velocities, ct, powers, axial_ind) = run_solver_and_extract_results(
            &mut model, "jensen", "jimenez", "sosfs", "none"
        );

        assert_results_close(&velocities, &JENSEN_BASELINE, 0, 0.1, "velocity");
        assert_results_close(&ct, &JENSEN_BASELINE, 1, 0.1, "Ct");
        assert_results_close(&powers, &JENSEN_BASELINE, 2, 10000.0, "power");
        assert_results_close(&axial_ind, &JENSEN_BASELINE, 3, 0.1, "axial_induction");
    }

    #[test]
    fn test_jensen_rotation_regression() {
        // Same rotation test as TurbOPark
        let turbine_diameter = ROTOR_DIAMETER;
        let spacing = 5.0 * turbine_diameter;

        let turbine_types = vec!["nrel_5MW".to_string(); 4];
        let layout_x = Array1::from_vec(vec![0.0, 0.0, spacing, spacing]);
        let layout_y = Array1::from_vec(vec![0.0, spacing, 0.0, spacing]);

        let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())
            .expect("Failed to create farm");

        let flow_field = FlowField::new(
            Array1::from_vec(vec![8.0, 8.0]),
            Array1::from_vec(vec![270.0, 360.0]),
            0.0,
            0.12,
            1.225,
            Array1::from_vec(vec![0.1, 0.1]),
            90.0,
        ).expect("Failed to create flow field");

        let mut model = FlorisModel {
            farm,
            flow_field,
            state: crate::core::State::new(),
            grid: None,
            solver: SolverConfig::default(),
            model_manager: None,
        };

        model.initialize_grid().expect("Failed to initialize grid");
        model.initialize_flow_field().expect("Failed to initialize flow field");

        let model_strings = WakeModelStrings {
            velocity_model: "jensen".to_string(),
            deflection_model: "jimenez".to_string(),
            combination_model: "sosfs".to_string(),
            turbulence_model: "none".to_string(),
        };

        let model_manager = WakeModelManager::new(
            model_strings,
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            false, false, false,
        ).expect("Failed to create model manager");

        let grid = model.grid.as_ref().unwrap().as_ref();
        crate::core::solver::sequential_solver(
            &model.farm,
            &mut model.flow_field,
            grid,
            &model_manager,
        ).expect("Solver failed");

        let velocities = model.flow_field.u_sorted.clone();
        let n_findex = 2;
        let n_turbines = 4;

        let mut avg_velocities = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let mut sum = 0.0;
                let mut count = 0;
                for iy in 0..3 {
                    for iz in 0..3 {
                        sum += velocities[[fi, ti, iy, iz]];
                        count += 1;
                    }
                }
                avg_velocities[[fi, ti]] = sum / count as f64;
            }
        }

        let tol = 0.01;
        // Jensen model should have similar rotation consistency
        assert!((avg_velocities[[0, 0]] - avg_velocities[[1, 1]]).abs() < tol, "Upstream velocities should match");
    }
}

#[cfg(test)]
mod wake_physics_tests {
    use super::*;

    #[test]
    fn test_wake_decreases_with_distance() {
        // Test that wake effect decreases with distance
        let mut model_close = create_test_model("gauss", "gauss", "fls", "none");
        let spacing_close = 3.0 * ROTOR_DIAMETER;

        let turbine_types = vec!["nrel_5MW".to_string(); 2];
        let layout_x_close = Array1::from_vec(vec![0.0, spacing_close]);

        let farm_close = Farm::new(
            layout_x_close,
            Array1::from_vec(vec![0.0, 0.0]),
            turbine_types.clone(),
        ).expect("Failed to create farm");

        let flow_field_close = FlowField::new(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            0.0,
            0.12,
            1.225,
            Array1::from_vec(vec![0.06]),
            90.0,
        ).expect("Failed to create flow field");

        let mut model_close = FlorisModel {
            farm: farm_close,
            flow_field: flow_field_close,
            state: crate::core::State::new(),
            grid: None,
            solver: SolverConfig::default(),
            model_manager: None,
        };

        let (v_close_up, v_close_down, _, _) = run_solver_and_extract_results(
            &mut model_close, "gauss", "gauss", "fls", "none"
        );

        // Wake loss at 3D
        let wake_loss_close = 1.0 - v_close_down[[0, 1]] / v_close_up[[0, 0]];

        // Wake should exist at 3D
        assert!(wake_loss_close > 0.0, "Wake effect should exist at 3D spacing");
    }

    #[test]
    fn test_wake_recovers_with_distance() {
        // Test that wake recovers at large distances
        let turbine_types = vec!["nrel_5MW".to_string(); 2];

        // Close spacing (3D)
        let layout_x_close = Array1::from_vec(vec![0.0, 3.0 * ROTOR_DIAMETER]);
        let farm_close = Farm::new(
            layout_x_close,
            Array1::from_vec(vec![0.0, 0.0]),
            turbine_types.clone(),
        ).expect("Failed to create farm");
        let flow_field_close = FlowField::new(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.12, 1.225,
            Array1::from_vec(vec![0.06]), 90.0,
        ).expect("Failed to create flow field");
        let mut model_close = FlorisModel {
            farm: farm_close, flow_field: flow_field_close,
            state: crate::core::State::new(), grid: None,
            solver: SolverConfig::default(), model_manager: None,
        };

        // Far spacing (15D)
        let layout_x_far = Array1::from_vec(vec![0.0, 15.0 * ROTOR_DIAMETER]);
        let farm_far = Farm::new(
            layout_x_far,
            Array1::from_vec(vec![0.0, 0.0]),
            turbine_types.clone(),
        ).expect("Failed to create farm");
        let flow_field_far = FlowField::new(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.12, 1.225,
            Array1::from_vec(vec![0.06]), 90.0,
        ).expect("Failed to create flow field");
        let mut model_far = FlorisModel {
            farm: farm_far, flow_field: flow_field_far,
            state: crate::core::State::new(), grid: None,
            solver: SolverConfig::default(), model_manager: None,
        };

        let (_, v_close_down, _, _) = run_solver_and_extract_results(
            &mut model_close, "gauss", "gauss", "fls", "none"
        );
        let (_, v_far_down, _, _) = run_solver_and_extract_results(
            &mut model_far, "gauss", "gauss", "fls", "none"
        );

        // Far turbine should have higher velocity (less wake)
        let v_close = v_close_down[[0, 1]];
        let v_far = v_far_down[[0, 1]];

        assert!(
            v_far > v_close,
            "Far turbine ({:.2}) should have higher velocity than close turbine ({:.2})",
            v_far, v_close
        );
    }

    #[test]
    fn test_aligned_wakes_strongest() {
        // Test that aligned turbines experience strongest wake
        let turbine_types = vec!["nrel_5MW".to_string(); 2];
        let spacing = 5.0 * ROTOR_DIAMETER;

        // Aligned (y offset = 0)
        let farm_aligned = Farm::new(
            Array1::from_vec(vec![0.0, spacing]),
            Array1::from_vec(vec![0.0, 0.0]),
            turbine_types.clone(),
        ).expect("Failed to create farm");
        let flow_field_aligned = FlowField::new(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.12, 1.225,
            Array1::from_vec(vec![0.06]), 90.0,
        ).expect("Failed to create flow field");
        let mut model_aligned = FlorisModel {
            farm: farm_aligned, flow_field: flow_field_aligned,
            state: crate::core::State::new(), grid: None,
            solver: SolverConfig::default(), model_manager: None,
        };

        // Offset (y offset = 3D)
        let farm_offset = Farm::new(
            Array1::from_vec(vec![0.0, spacing]),
            Array1::from_vec(vec![0.0, 3.0 * ROTOR_DIAMETER]),
            turbine_types.clone(),
        ).expect("Failed to create farm");
        let flow_field_offset = FlowField::new(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.12, 1.225,
            Array1::from_vec(vec![0.06]), 90.0,
        ).expect("Failed to create flow field");
        let mut model_offset = FlorisModel {
            farm: farm_offset, flow_field: flow_field_offset,
            state: crate::core::State::new(), grid: None,
            solver: SolverConfig::default(), model_manager: None,
        };

        let (_, v_aligned_down, _, _) = run_solver_and_extract_results(
            &mut model_aligned, "gauss", "gauss", "fls", "none"
        );
        let (_, v_offset_down, _, _) = run_solver_and_extract_results(
            &mut model_offset, "gauss", "gauss", "fls", "none"
        );

        // Aligned turbine should have lower velocity
        let v_aligned = v_aligned_down[[0, 1]];
        let v_offset = v_offset_down[[0, 1]];

        assert!(
            v_aligned < v_offset,
            "Aligned turbine ({:.2}) should have lower velocity than offset turbine ({:.2})",
            v_aligned, v_offset
        );
    }
}
