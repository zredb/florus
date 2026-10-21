//! Example: Sweeping the model variables
//!
//! This example demonstrates sweeping:
//!     1) The wind speeds
//!     2) The wind directions
//!     3) The turbulence intensities
//!     4) The yaw angles
//!     5) The power setpoints
//!     6) The disabled turbines
//!
//! Scenario aligned 1:1 with examples/007_sweeping_variables.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A two-turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 126.0 * 5.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    // 1) Sweep the wind speeds (the wind direction and the turbulence
    // intensity are constant)
    let wind_speeds: Vec<f64> = (0..50).map(|i| 5.0 + 0.1 * i as f64).collect();
    let n = wind_speeds.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Wind speed sweep [m/s], upstream [kW], downstream [kW]:");
    for i in 0..n {
        println!(
            "{:.1}, {:.1}, {:.1}",
            wind_speeds[i],
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0
        );
    }

    // 2) Sweep the wind directions
    let wind_directions: Vec<f64> = (250..290).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Wind direction sweep [deg], upstream [kW], downstream [kW]:");
    for i in 0..n {
        println!(
            "{:.0}, {:.1}, {:.1}",
            wind_directions[i],
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0
        );
    }

    // 3) Sweep the turbulence intensities
    let turbulence_intensities: Vec<f64> = (0..17).map(|i| 0.03 + 0.01 * i as f64).collect();
    let n = turbulence_intensities.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(turbulence_intensities.clone()),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Turbulence intensity sweep, upstream [kW], downstream [kW]:");
    for i in 0..n {
        println!(
            "{:.2}, {:.1}, {:.1}",
            turbulence_intensities[i],
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0
        );
    }

    // 4) Sweep the upstream yaw angle
    let n_yaw = 100;
    let yaw_angles_upstream: Vec<f64> = (0..n_yaw)
        .map(|i| -30.0 + 60.0 * i as f64 / (n_yaw - 1) as f64)
        .collect();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_yaw]),
        Array1::from_vec(vec![8.0; n_yaw]),
        Array1::from_vec(vec![0.06; n_yaw]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    let mut yaw_angles = Array2::zeros((n_yaw, 2));
    for (i, yaw) in yaw_angles_upstream.iter().enumerate() {
        yaw_angles[[i, 0]] = *yaw;
    }
    fmodel.set_yaw_angles(yaw_angles)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Upstream yaw angle sweep [deg], upstream [kW], downstream [kW]:");
    for i in 0..n_yaw {
        println!(
            "{:.3}, {:.1}, {:.1}",
            yaw_angles_upstream[i],
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0
        );
    }

    // 5) Sweep the upstream power rating (the derating needs the
    // simple-derating operation model)
    fmodel.reset_operation();
    fmodel.set_operation_model("simple-derating")?;

    let rated_power = 5.0e6;
    let n_derating_levels = 150;
    let upstream_power_setpoint: Vec<f64> = (0..n_derating_levels)
        .map(|i| rated_power * 0.5 * i as f64 / (n_derating_levels - 1) as f64)
        .collect();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_derating_levels]),
        Array1::from_vec(vec![8.0; n_derating_levels]),
        Array1::from_vec(vec![0.06; n_derating_levels]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    let mut power_setpoints = Array2::from_elem((n_derating_levels, 2), rated_power);
    for (i, setpoint) in upstream_power_setpoint.iter().enumerate() {
        power_setpoints[[i, 0]] = *setpoint;
    }
    fmodel.set_operation(None, Some(power_setpoints), None, None, None, None)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Upstream power setpoint sweep [kW], upstream [kW], downstream [kW]:");
    for i in 0..n_derating_levels {
        println!(
            "{:.1}, {:.1}, {:.1}",
            upstream_power_setpoint[i] / 1000.0,
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0
        );
    }

    // 6) Sweep through the disabling turbine combinations
    fmodel.reset_operation();
    let disable_combinations = [[false, false], [true, false], [false, true], [true, true]];
    let names = ["None", "T0", "T1", "T0 & T1"];
    let n_combinations = disable_combinations.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_combinations]),
        Array1::from_vec(vec![8.0; n_combinations]),
        Array1::from_vec(vec![0.06; n_combinations]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((n_combinations, 2), false);
    for (i, combination) in disable_combinations.iter().enumerate() {
        disable_turbines[[i, 0]] = combination[0];
        disable_turbines[[i, 1]] = combination[1];
    }
    fmodel.set_operation(None, None, None, None, None, Some(disable_turbines))?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Disabled turbine combinations, upstream [kW], downstream [kW]:");
    for i in 0..n_combinations {
        println!(
            "{}, {:.1}, {:.1}",
            names[i],
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0
        );
    }

    Ok(())
}
