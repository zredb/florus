//! Regression tests for the query/optimization APIs exposed by the
//! FLORIS examples (Python cross-validation, 2026-10).

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use ndarray::Array2;

/// `set_yaw_angles` must invalidate the solved state so the next `run()`
/// re-sorts the yaw angles into the solver's sorted order.
///
/// Before the fix, `FlorisModel::set_yaw_angles` left the model "initialized",
/// so `run()` skipped `initialize_domain()` and the solver used the
/// previously sorted yaw array — applying the yaw angles to the WRONG
/// turbines whenever the wake sort order (by rotated x) differs from the
/// layout order. Python resets `core.state` to UNINITIALIZED in `set()`
/// (floris_model.py:398).
#[test]
fn test_set_yaw_angles_invalidates_state() {
    let mut model = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");

    // 2x2 grid, 5D spacing
    let d = 126.0;
    let spacing = 5.0 * d;
    model
        .set_layout(
            &Array1::from_vec(vec![0.0, 0.0, spacing, spacing]),
            &Array1::from_vec(vec![0.0, spacing, 0.0, spacing]),
        )
        .expect("set_layout");

    // Wind from the east (90 deg): the sorted order is (t0, t2, t1, t3) while
    // the layout order is (t0, t1, t2, t3)
    model
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![90.0]),
            Array1::from_vec(vec![0.06]),
        )
        .expect("set_wind_conditions");
    model.run().expect("run");

    let base = model.get_turbine_powers().expect("powers");

    // Yaw turbine 0 only (the upstream turbine in the sorted order)
    let mut yaw = Array2::zeros((1, 4));
    yaw[[0, 0]] = 25.0;
    model.set_yaw_angles(yaw).expect("set_yaw_angles");
    model.run().expect("run after yaw");

    let yawed = model.get_turbine_powers().expect("powers");

    let d0 = yawed[[0, 0]] - base[[0, 0]];
    let d2 = yawed[[0, 2]] - base[[0, 2]];

    assert!(
        d0 < 0.0,
        "yawed turbine 0 should lose power (cosine loss), got {d0}"
    );
    assert!(
        d2.abs() < 1.0,
        "turbine 2 (not yawed) must be nearly unchanged, got {d2}"
    );
}

/// The Serial-Refine yaw optimizer must find the same optimum as Python
/// (examples_control_optimization/003_opt_yaw_multiple_ws.py): a 3-turbine
/// row at 270 deg, 8 m/s, 25 deg maximum yaw.
#[test]
fn test_serial_refine_matches_python_optimum() {
    use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
    use florus::optimization::yaw_optimization::YawOptimization;

    let mut model = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    let d = 126.0;
    model
        .set_layout(
            &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .expect("set_layout");
    model
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .expect("set_wind_conditions");
    model.run().expect("run");

    let mut opt = YawOptimizationSR::new();
    let result = opt.optimize(&mut model, None).expect("optimize");

    // Python (examples_control_optimization/008): [25.0, 22.65625, 0.0]
    assert!(
        (result.yaw_angles[[0, 0]] - 25.0).abs() < 1e-9,
        "turbine 0 yaw: {}",
        result.yaw_angles[[0, 0]]
    );
    assert!(
        (result.yaw_angles[[0, 1]] - 22.65625).abs() < 1e-9,
        "turbine 1 yaw: {}",
        result.yaw_angles[[0, 1]]
    );
    assert_eq!(result.yaw_angles[[0, 2]], 0.0);
}

/// The Geometric yaw optimizer must reproduce Python's analytical result
/// (examples_control_optimization/008): [19.9952335557674, 19.9952335557674, 0].
#[test]
fn test_geometric_yaw_matches_python() {
    use florus::optimization::yaw_optimization::yaw_optimizer_geometric::YawOptimizationGeometric;
    use florus::optimization::yaw_optimization::YawOptimization;

    let mut model = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    let d = 126.0;
    model
        .set_layout(
            &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .expect("set_layout");
    model
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .expect("set_wind_conditions");
    model.run().expect("run");

    let mut opt = YawOptimizationGeometric::new();
    let result = opt.optimize(&mut model, None).expect("optimize");

    let expected = 19.9952335557674;
    assert!(
        (result.yaw_angles[[0, 0]] - expected).abs() < 1e-9,
        "turbine 0 yaw: {}",
        result.yaw_angles[[0, 0]]
    );
    assert!(
        (result.yaw_angles[[0, 1]] - expected).abs() < 1e-9,
        "turbine 1 yaw: {}",
        result.yaw_angles[[0, 1]]
    );
    assert_eq!(result.yaw_angles[[0, 2]], 0.0);
}
