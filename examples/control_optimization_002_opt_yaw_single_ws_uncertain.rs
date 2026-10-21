//! Example: Optimize yaw for a single wind speed and multiple wind directions.
//! Compare certain and uncertain results.
//!
//! Use the serial-refine method to optimize the yaw angles for a 3-turbine
//! wind farm.  In one case use the FlorisModel without uncertainty and in the
//! other use the UncertainFlorisModel with a wind direction standard
//! deviation of 3 degrees.  Compare the results.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/002_opt_yaw_single_ws_uncertain.py.

use florus::optimization::yaw_optimization::{
    UncertainPowerEvaluator, YawOptimization, YawOptimizationSR,
};
use florus::types::Array1;
use florus::uncertain::UncertainFlorisModel;
use florus::{FlorisModel, Result, TimeSeries};

fn main() -> Result<()> {
    // Load the floris model and uncertain floris model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut ufmodel = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml")?,
        3.0,
        None,
    )?;

    // Define an inflow that keeps wind speed and TI constant while sweeping
    // the wind directions
    let wind_directions = Array1::from_iter((250..290).map(|d| d as f64));
    let n_conditions = wind_directions.len();
    let wind_speeds = Array1::from_elem(n_conditions, 8.0);
    let turbulence_intensities = Array1::from_elem(n_conditions, 0.06);
    let time_series = TimeSeries::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        turbulence_intensities.clone(),
    )?;

    // Reinitialize as a 3-turbine farm using the above inflow
    let d = 126.0; // Rotor diameter for the NREL 5 MW
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    ufmodel.set_layout(&layout_x, &layout_y)?;
    ufmodel.set_wind_data(&time_series)?;

    // Initialize optimizer object and run optimization using the Serial-Refine
    // method
    println!("++++++++++CERTAIN++++++++++++");
    let mut yaw_opt = YawOptimizationSR::new();
    let result_certain = yaw_opt.optimize(&mut fmodel, None)?;

    // Repeat with uncertain model
    println!("++++++++++UNCERTAIN++++++++++++");
    let mut yaw_opt_u = YawOptimizationSR::new();
    let result_uncertain = {
        let mut evaluator = UncertainPowerEvaluator::new(&mut ufmodel);
        yaw_opt_u.optimize_with_evaluator(&mut evaluator, None)?
    };

    // Show the yaw and turbine results
    println!("certain: wind direction, yaw angles, baseline farm power, optimized farm power");
    for i in 0..n_conditions {
        println!(
            "certain {i}: {:.0}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6}",
            wind_directions[i],
            result_certain.yaw_angles[[i, 0]],
            result_certain.yaw_angles[[i, 1]],
            result_certain.yaw_angles[[i, 2]],
            result_certain.farm_powers_baseline[i],
            result_certain.farm_powers_opt[i],
        );
    }

    println!("uncertain: wind direction, yaw angles, baseline farm power, optimized farm power");
    for i in 0..n_conditions {
        println!(
            "uncertain {i}: {:.0}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6}",
            wind_directions[i],
            result_uncertain.yaw_angles[[i, 0]],
            result_uncertain.yaw_angles[[i, 1]],
            result_uncertain.yaw_angles[[i, 2]],
            result_uncertain.farm_powers_baseline[i],
            result_uncertain.farm_powers_opt[i],
        );
    }

    Ok(())
}
