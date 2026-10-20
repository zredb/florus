//! Example: Multi-dimensional power/thrust coefficient with 2 Hs values
//!
//! This example follows the previous example but shows the effect of
//! changing the wave height `Hs`: the wind speeds are swept twice, once with
//! `Hs = 3.1` and once with `Hs = 1.0`, using the array-based
//! multidimensional functionality of FLORIS v4.6 (a different
//! multidimensional condition per wind condition).
//!
//! Scenario aligned 1:1 with
//! examples_multidim/002_multi_dimensional_cp_ct_2Hs.py.

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::MultidimConditions;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // Initialize FLORIS with the given input file
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct.yaml")?;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Conditions to evaluate: a sweep of wind speeds, twice
    let n_wind_speeds = 16;
    let mut wind_speeds: Vec<f64> = (0..n_wind_speeds)
        .map(|i| 5.0 + (20.0 - 5.0) * i as f64 / (n_wind_speeds - 1) as f64)
        .collect();
    let first_sweep = wind_speeds.clone();
    wind_speeds.extend_from_slice(&first_sweep);
    let n_conditions = wind_speeds.len();

    // The wave height is 3.1 for the first sweep and 1.0 for the second one
    let mut hs: Vec<f64> = vec![3.1; n_wind_speeds];
    hs.extend(std::iter::repeat(1.0).take(n_wind_speeds));

    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_conditions]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.06; n_conditions]),
    )?;
    time_series.multidim_conditions = Some(MultidimConditions {
        tp: Array1::from_vec(vec![2.5; n_conditions]),
        hs: Some(Array1::from_vec(hs)),
        ti: None,
    });

    // Set the wind data onto the FlorisModel
    fmodel.set_wind_data(&time_series)?;

    // Calculate wakes
    fmodel.run()?;

    // Collect the turbine powers in kW
    let turbine_powers = fmodel.get_turbine_powers()? / 1000.0;

    println!("Wind speeds [m/s]:");
    let speeds: Vec<String> = wind_speeds.iter().map(|v| format!("{v:.1}")).collect();
    println!("  {}", speeds.join(" "));
    for t in 0..3 {
        let powers: Vec<String> = (0..n_conditions)
            .map(|i| format!("{:.1}", turbine_powers[[i, t]]))
            .collect();
        println!("Turbine {t} power [kW]:");
        println!("  {}", powers.join(" "));
    }

    Ok(())
}
