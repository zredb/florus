//! Example: Optimizing yaw angles with disabled turbines
//!
//! This example demonstrates how to optimize yaw angles in FLORUS, when some
//! turbines are disabled. The example optimization is run using both
//! YawOptimizerSR and YawOptimizerGeometric, the two yaw optimizers that
//! support disabling turbines.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/008_optimize_yaw_with_disabled_turbines.py.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_geometric::YawOptimizationGeometric;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::YawOptimization;
use florus::types::{Array1, Array2, Float};
use ndarray::Array;

/// numpy-style 2D array print: " [[v v v]\n  [v v v]]"
fn np_print(yaw: &Array2) -> String {
    let mut s = String::from(" [");
    for (i, row) in yaw.outer_iter().enumerate() {
        if i > 0 {
            s.push_str("\n  ");
        }
        for (j, &v) in row.iter().enumerate() {
            if j > 0 {
                s.push(' ');
            }
            s.push_str(&format!("{v}"));
        }
    }
    s.push(']');
    s
}

fn main() -> anyhow::Result<()> {
    // Load a 3-turbine model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set wind conditions to be the same for two cases
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 8.0]),
        Array1::from_vec(vec![270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06]),
    )?;
    fmodel.run()?;

    // First run the case where all turbines are active and print results
    let mut yaw_opt = YawOptimizationSR::new();
    let df_opt = yaw_opt.optimize(&mut fmodel, None)?;
    println!(
        "Serial Refine optimized yaw angles (all turbines active) [deg]:\n{}",
        np_print(&df_opt.yaw_angles)
    );

    let mut yaw_opt = YawOptimizationGeometric::new();
    let df_opt = yaw_opt.optimize(&mut fmodel, None)?;
    println!(
        "\nGeometric optimized yaw angles (all turbines active) [deg]:\n{}",
        np_print(&df_opt.yaw_angles)
    );

    // Disable turbines (different pattern for each of the two cases)
    // First case: disable the middle turbine
    // Second case: disable the front turbine
    fmodel.set_operation_model("mixed")?;
    let disable = ndarray::Array2::<bool>::from_shape_fn((2, 3), |(i, j)| {
        (i == 0 && j == 1) || (i == 1 && j == 0)
    });
    fmodel.set_operation(None, None, None, None, None, Some(disable))?;
    fmodel.run()?;

    // Rerun optimizations and print results
    let mut yaw_opt = YawOptimizationSR::new();
    let df_opt = yaw_opt.optimize(&mut fmodel, None)?;
    println!(
        "\nSerial Refine optimized yaw angles (some turbines disabled) [deg]:\n{}",
        np_print(&df_opt.yaw_angles)
    );

    let mut yaw_opt = YawOptimizationGeometric::new();
    let df_opt = yaw_opt.optimize(&mut fmodel, None)?;
    println!(
        "\nGeometric optimized yaw angles (some turbines disabled) [deg]:\n{}",
        np_print(&df_opt.yaw_angles)
    );

    Ok(())
}
