//! Example: Multi-dimensional power/thrust coefficient with a TI dimension
//!
//! This example sweeps the wind speed and evaluates the turbine powers for
//! two turbulence intensities. The turbulence intensity is the
//! multidimensional condition that selects the power/thrust table; with
//! `TI = 0.10` the handler finds the nearest defined value (0.08) and uses
//! that data.
//!
//! Scenario aligned 1:1 with
//! examples_multidim/003_multi_dimensional_cp_ct_TI.py.

use florus::floris_config::MultiDimConditionValue;
use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    // Initialize FLORIS with the given input file
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct_TI.yaml")?;

    // Set both cases to 3 turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Use a sweep of wind speeds
    let wind_speeds: Vec<f64> = (0..150).map(|i| 5.0 + 0.1 * i as f64).collect();
    let n_conditions = wind_speeds.len();

    // In this case, we run() the fmodel multiple times, one for each
    // turbulence intensity. The multidim condition is a scalar, which is
    // broadcast over the wind speeds.
    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_conditions]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.06; n_conditions]),
    )?;
    time_series.multidim_conditions = None;
    fmodel.set_wind_data(&time_series)?;

    for ti in [0.06, 0.10] {
        // The multidimensional data handler finds the nearest defined value of
        // the turbulence intensity and uses that data
        let mut conditions: HashMap<String, MultiDimConditionValue> = HashMap::new();
        conditions.insert("TI".to_string(), MultiDimConditionValue::Scalar(ti));
        fmodel.set_multidim_conditions(conditions)?;

        fmodel.run()?;
        let turbine_powers = fmodel.get_turbine_powers()? / 1000.0;

        println!("TI = {ti:.2}:");
        for t in 0..3 {
            let powers: Vec<String> = (0..n_conditions)
                .map(|i| format!("{:.1}", turbine_powers[[i, t]]))
                .collect();
            println!("  turbine {t} power [kW]: {}", powers.join(" "));
        }
    }

    Ok(())
}
