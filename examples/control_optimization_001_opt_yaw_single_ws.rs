//! Example: Optimize the yaw angles for a single wind speed
//!
//! This example demonstrates how to perform a yaw optimization for a sweep of
//! wind directions with a constant wind speed, using the Serial-Refine method.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/001_opt_yaw_single_ws.py: a 3-turbine row on
//! gch.yaml, the wind directions from 0 to 357 degrees in 3 degree steps, a wind
//! speed of 8 m/s and TI = 0.06.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::YawOptimization;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A 3-turbine row
    let d = 126.0; // The rotor diameter of the NREL 5 MW
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // The wind direction is swept while the wind speed and the turbulence
    // intensity are held constant
    let wind_directions: Vec<f64> = (0..120).map(|i| 3.0 * i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    // The Serial-Refine optimizer (Python defaults: yaw in [0, 25] deg,
    // Ny_passes = [5, 4]).
    let mut yaw_opt = YawOptimizationSR::new();
    let result = yaw_opt.optimize(&mut fmodel, None)?;

    println!("Wind direction [deg], yaw t0 [deg], yaw t1 [deg], yaw t2 [deg], farm power opt [W], farm power baseline [W]");
    for i in 0..n {
        println!(
            "{:.0}, {:.3}, {:.3}, {:.3}, {:.1}, {:.1}",
            wind_directions[i],
            result.yaw_angles[[i, 0]],
            result.yaw_angles[[i, 1]],
            result.yaw_angles[[i, 2]],
            result.farm_powers_opt[i],
            result.farm_powers_baseline[i]
        );
    }

    // The baseline and the optimized farm power
    println!(
        "Baseline farm power: {:.3} W, optimized farm power: {:.3} W ({:.3} %)",
        result.baseline_power, result.optimized_power, result.improvement_percentage
    );

    Ok(())
}
