//! Example: Heterogeneous speed-up by wind direction and wind speed
//!
//! This example builds a `HeterogeneousMap` whose speed-ups are defined for
//! two wind directions and two wind speeds, applies it to a set of wind
//! conditions and reports the multipliers that are selected for each
//! condition together with the resulting turbine powers.
//!
//! Scenario aligned 1:1 with
//! examples_heterogeneous/003_heterogeneous_speedup_by_wd_and_ws.py.

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::HeterogeneousMap;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    // The speed-ups imply no heterogeneity for the first wind direction
    // (270 deg) and heterogeneity for the second wind direction (90 deg),
    // where the speed-ups depend on the wind speed.
    let heterogeneous_map = HeterogeneousMap::new(
        Array1::from_vec(vec![0.0, 0.0, 250.0, 500.0, 500.0]),
        Array1::from_vec(vec![0.0, 500.0, 250.0, 0.0, 500.0]),
        Array2::from_shape_vec(
            (4, 5),
            vec![
                1.0, 1.0, 1.0, 1.0, 1.0, // (270 deg, 5 m/s)
                1.0, 1.0, 1.0, 1.0, 1.0, // (270 deg, 10 m/s)
                1.5, 1.0, 1.25, 1.5, 1.0, // (90 deg, 5 m/s)
                1.0, 1.5, 1.25, 1.0, 1.5, // (90 deg, 10 m/s)
            ],
        )?,
        None,
        Some(Array1::from_vec(vec![270.0, 270.0, 90.0, 90.0])),
        Some(Array1::from_vec(vec![5.0, 10.0, 5.0, 10.0])),
        florus::core::models::InterpMethod::Linear,
    )?;

    // Initialize FlorisModel
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change the layout to a 2 turbine layout within the heterogeneous domain
    fmodel.set_layout(
        &Array1::from_vec(vec![200.0, 200.0]),
        &Array1::from_vec(vec![50.0, 450.0]),
    )?;

    // Three wind conditions with the heterogeneous map applied
    let wind_directions = Array1::from_vec(vec![275.0, 95.0, 75.0]);
    let wind_speeds = Array1::from_vec(vec![7.0, 6.2, 8.0]);
    fmodel.set_wind_conditions(
        wind_speeds.clone(),
        wind_directions.clone(),
        Array1::from_vec(vec![0.06; 3]),
    )?;

    // The multipliers selected for each condition
    let config = heterogeneous_map
        .get_heterogeneous_inflow_config(wind_directions.clone(), wind_speeds.clone())?;
    for fi in 0..config.speed_multipliers.shape()[0] {
        let multipliers: Vec<String> = (0..config.x.len())
            .map(|p| format!("{:.2}", config.speed_multipliers[[fi, p]]))
            .collect();
        println!(
            "Condition {} (wd {:.1} deg, ws {:.1} m/s): multipliers [{}]",
            fi,
            wind_directions[fi],
            wind_speeds[fi],
            multipliers.join(", ")
        );
    }

    // Apply the heterogeneous inflow and run the simulation
    fmodel.set_heterogeneous_inflow_config(config)?;
    fmodel.run()?;

    let powers = fmodel.get_turbine_powers()? / 1000.0;
    println!("Turbine powers [kW]:");
    for t in 0..powers.shape()[1] {
        let row: Vec<String> = (0..powers.shape()[0])
            .map(|i| format!("{:.1}", powers[[i, t]]))
            .collect();
        println!("  turbine {t}: {}", row.join(" "));
    }

    Ok(())
}
