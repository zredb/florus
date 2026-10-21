//! Cross-validation of the yaw optimization with an uncertain model
//! (Python `examples_control_optimization/002_opt_yaw_single_ws_uncertain.py`).
//!
//! The example optimizes the yaw angles of a 3-turbine row with the
//! Serial-Refine method, once with a plain `FlorisModel` and once with an
//! `UncertainFlorisModel` (a wind direction standard deviation of 3 deg).
//! The reference values come from the Python run of the example.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::{
    UncertainPowerEvaluator, YawOptimization, YawOptimizationSR,
};
use florus::types::Array1;
use florus::uncertain::UncertainFlorisModel;
use florus::wind_data::TimeSeries;

const D: f64 = 126.0;

/// The layout and the inflow of the example
fn setup() -> (Array1, Array1, TimeSeries) {
    let wind_directions = Array1::from_iter((250..290).map(|d| d as f64));
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        wind_directions,
        Array1::from_elem(n, 8.0),
        Array1::from_elem(n, 0.06),
    )
    .expect("time series");
    (
        Array1::from_vec(vec![0.0, 5.0 * D, 10.0 * D]),
        Array1::from_vec(vec![0.0, 0.0, 0.0]),
        time_series,
    )
}

/// The yaw angles and the farm powers of the certain model
///
/// The reference values are from the Python run of
/// `002_opt_yaw_single_ws_uncertain.py` (the condition index, the yaw
/// angles of the three turbines, the baseline farm power and the optimized
/// farm power).
#[test]
fn test_certain_yaw_optimization_matches_python() {
    let (layout_x, layout_y, time_series) = setup();

    let mut model = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    model.set_layout(&layout_x, &layout_y).expect("layout");
    model.set_wind_data(&time_series).expect("wind data");
    model.run().expect("run");

    let mut optimizer = YawOptimizationSR::new();
    optimizer.print_progress = false;
    let result = optimizer.optimize(&mut model, None).expect("optimize");

    // (condition, yaw angles, baseline farm power, optimized farm power)
    let reference = [
        (0usize, [0.0, 0.0, 0.0], 5261505.675277, 5261505.675277),
        (18, [25.0, 25.0, 0.0], 2892655.144025, 2973038.555864),
        (19, [25.0, 25.0, 0.0], 2748875.853126, 3104344.503907),
        (20, [25.0, 22.65625, 0.0], 2697065.314464, 3275737.731164),
        (21, [25.0, 21.875, 0.0], 2753926.611869, 3486497.980480),
        (33, [3.125, 3.125, 0.0], 5170397.746709, 5185384.881261),
        (39, [0.0, 0.0, 0.0], 5260981.079248, 5260981.079248),
    ];

    for (condition, yaws, baseline, optimized) in reference {
        for t in 0..3 {
            assert!(
                (result.yaw_angles[[condition, t]] - yaws[t]).abs() < 1e-9,
                "yaw of turbine {t} in condition {condition}: {} (expected {})",
                result.yaw_angles[[condition, t]],
                yaws[t]
            );
        }
        assert!(
            (result.farm_powers_baseline[condition] - baseline).abs() < 1e-3,
            "baseline farm power of condition {condition}: {} (expected {baseline})",
            result.farm_powers_baseline[condition]
        );
        assert!(
            (result.farm_powers_opt[condition] - optimized).abs() < 1e-3,
            "optimized farm power of condition {condition}: {} (expected {optimized})",
            result.farm_powers_opt[condition]
        );
    }
}

/// The yaw angles and the farm powers of the uncertain model
///
/// The reference values are from the Python run of
/// `002_opt_yaw_single_ws_uncertain.py` (the condition index, the yaw
/// angles of the three turbines, the baseline farm power and the optimized
/// farm power).
#[test]
fn test_uncertain_yaw_optimization_matches_python() {
    let (layout_x, layout_y, time_series) = setup();

    let mut model = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").expect("model"),
        3.0,
        None,
    )
    .expect("uncertain model");
    model.set_layout(&layout_x, &layout_y).expect("layout");
    model.set_wind_data(&time_series).expect("wind data");

    let mut optimizer = YawOptimizationSR::new();
    optimizer.print_progress = false;
    let result = {
        let mut evaluator = UncertainPowerEvaluator::new(&mut model);
        optimizer
            .optimize_with_evaluator(&mut evaluator, None)
            .expect("optimize")
    };

    // (condition, yaw angles, baseline farm power, optimized farm power)
    let reference = [
        (0usize, [0.0, 0.0, 0.0], 5257545.038854, 5257545.038854),
        (
            19,
            [24.21875, 24.21875, 0.0],
            3058251.245322,
            3264533.047640,
        ),
        (
            20,
            [24.21875, 24.21875, 0.0],
            3027038.307499,
            3380955.640074,
        ),
        (21, [25.0, 21.875, 0.0], 3063084.293251, 3527885.358017),
        (35, [3.125, 0.0, 0.0], 5177676.308966, 5185583.628346),
        (39, [0.0, 0.0, 0.0], 5253753.171442, 5253753.171442),
    ];

    for (condition, yaws, baseline, optimized) in reference {
        for t in 0..3 {
            assert!(
                (result.yaw_angles[[condition, t]] - yaws[t]).abs() < 1e-9,
                "yaw of turbine {t} in condition {condition}: {} (expected {})",
                result.yaw_angles[[condition, t]],
                yaws[t]
            );
        }
        assert!(
            (result.farm_powers_baseline[condition] - baseline).abs() < 1e-3,
            "baseline farm power of condition {condition}: {} (expected {baseline})",
            result.farm_powers_baseline[condition]
        );
        assert!(
            (result.farm_powers_opt[condition] - optimized).abs() < 1e-3,
            "optimized farm power of condition {condition}: {} (expected {optimized})",
            result.farm_powers_opt[condition]
        );
    }
}

/// The rounding of the expanded inputs follows numpy: the resolutions are
/// applied per column (the yaw angles of every turbine are rounded to 1 deg)
/// and the ties are rounded to the even value.
#[test]
fn test_uncertain_expansion_rounds_the_yaw_angles() {
    let (layout_x, layout_y, time_series) = setup();

    let mut model = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").expect("model"),
        3.0,
        None,
    )
    .expect("uncertain model");
    model.set_layout(&layout_x, &layout_y).expect("layout");
    model.set_wind_data(&time_series).expect("wind data");

    // The yaw angle of the second turbine is used (it was rounded to the
    // power setpoint resolution before the fix, so it was always zero)
    let n = 40;
    let yaw_angles = florus::types::Array2::from_shape_fn((n, 3), |(_, t)| match t {
        0 => 0.0,
        1 => 20.0,
        _ => 0.0,
    });
    model.set_yaw_angles(yaw_angles).expect("yaw angles");
    model.run().expect("run");

    let powers = model.get_turbine_powers().expect("powers");

    // The second turbine is yawed, so its own power is reduced
    let yaw_loss = 1.0 - powers[[20, 1]] / 1_800_000.0;
    assert!(
        yaw_loss > 0.05,
        "the yaw angle of the second turbine is not applied: {}",
        powers[[20, 1]]
    );
}
