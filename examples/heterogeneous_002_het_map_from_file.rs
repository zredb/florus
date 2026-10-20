//! Example: Heterogeneous inflow with wind data
//!
//! This example sets up a heterogeneous inflow condition for a set of wind
//! conditions, first with the per-condition speed multipliers and then with
//! a `HeterogeneousMap` that selects the multipliers by wind direction.
//!
//! Scenario aligned 1:1 with
//! examples_heterogeneous/002_heterogeneous_using_wind_data.py.

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::{HeterogeneousInflowConfig, HeterogeneousMap};
use florus::types::{Array1, Array2};
use ndarray::Array;

fn main() -> anyhow::Result<()> {
    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change the layout to a 4 turbine layout in a box
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 0.0, 500.0, 500.0]),
        &Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0]),
    )?;

    // Four wind directions, constant wind speed and turbulence intensity
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0; 4]),
        Array1::from_vec(vec![269.0, 270.0, 271.0, 282.0]),
        Array1::from_vec(vec![0.06; 4]),
    )?;

    // The x/y locations of the heterogeneous inflow box
    let x_locs = Array1::from_vec(vec![-500.0, -500.0, 1000.0, 1000.0]);
    let y_locs = Array1::from_vec(vec![-500.0, 1000.0, -500.0, 1000.0]);

    // Per-condition speed multipliers: the first three conditions share the
    // same multipliers, the last one uses larger speed-ups
    let speed_multipliers = Array2::from_shape_vec(
        (4, 4),
        vec![
            1.0, 1.25, 1.0, 1.25, //
            1.0, 1.25, 1.0, 1.25, //
            1.0, 1.25, 1.0, 1.25, //
            1.0, 1.35, 1.0, 1.35,
        ],
    )?;

    fmodel.set_heterogeneous_inflow_config(HeterogeneousInflowConfig {
        x: x_locs.clone(),
        y: y_locs.clone(),
        z: None,
        wind_speeds: None,
        wind_directions: None,
        speed_multipliers: speed_multipliers.clone(),
    })?;
    fmodel.run()?;

    let powers = fmodel.get_turbine_powers()? / 1000.0;
    println!("Per-condition speed multipliers:");
    for t in 0..4 {
        let row: Vec<String> = (0..4).map(|i| format!("{:.1}", powers[[i, t]])).collect();
        println!("  turbine {t}: {}", row.join(" "));
    }

    // Repeat using a HeterogeneousMap that selects the multipliers by wind
    // direction (270 deg -> 1.25x, 280 deg -> 1.35x)
    let map_multipliers =
        Array2::from_shape_vec((2, 4), vec![1.0, 1.25, 1.0, 1.25, 1.0, 1.35, 1.0, 1.35])?;
    let heterogeneous_map = HeterogeneousMap::new(
        x_locs,
        y_locs,
        map_multipliers,
        None,
        Some(Array1::from_vec(vec![270.0, 280.0])),
        None,
        florus::core::models::InterpMethod::Linear,
    )?;

    let config = heterogeneous_map.get_heterogeneous_inflow_config(
        Array1::from_vec(vec![269.0, 270.0, 271.0, 282.0]),
        Array1::from_vec(vec![8.0; 4]),
    )?;
    fmodel.set_heterogeneous_inflow_config(config)?;
    fmodel.run()?;

    let powers = fmodel.get_turbine_powers()? / 1000.0;
    println!("HeterogeneousMap (by wind direction):");
    for t in 0..4 {
        let row: Vec<String> = (0..4).map(|i| format!("{:.1}", powers[[i, t]])).collect();
        println!("  turbine {t}: {}", row.join(" "));
    }

    Ok(())
}
