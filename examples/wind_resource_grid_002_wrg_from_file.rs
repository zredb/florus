//! Example: Set a FLORIS model with a wind resource grid (WRG)
//!
//! This example applies a WRG to a FLORIS model: the wind roses of the WRG
//! are generated for each turbine of the layout and applied as the model's
//! wind data.
//!
//! Scenario aligned 1:1 with
//! examples_wind_resource_grid/002_set_floris_with_wrg.py.
//!

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::WindRoseWRG;

fn main() -> anyhow::Result<()> {
    // Read the WRG file
    let mut wind_rose_wrg = WindRoseWRG::new("examples/inputs/wrg_example.wrg", None, None, None)?;

    // Use 2 m/s steps and a fixed TI of 7% for all wind directions
    wind_rose_wrg.set_wind_speeds(Array1::from_iter((0..13).map(|i| i as f64 * 2.0)));
    wind_rose_wrg.set_ti_table(0.07);

    // Select a turbine layout: a line of turbines along the diagonals of the
    // WRG grid
    let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array1::from_vec(vec![0.0, 1000.0, 2000.0]);

    // Set the turbine layout in the WindRoseWRG object; the wind roses are
    // generated for each turbine
    wind_rose_wrg.set_layout(layout_x.clone(), layout_y.clone())?;

    // Apply the wind_rose_wrg to a FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.set_wind_data(&wind_rose_wrg)?;
    fmodel.run()?;

    // Run the model and get the expected turbine powers and farm power
    let expected_turbine_powers = fmodel.get_expected_turbine_powers(None)?;
    let expected_farm_power = fmodel.get_expected_farm_power(None, None)?;

    println!("=====================================");
    println!("Expected turbine powers: {:?}", expected_turbine_powers);
    println!("Expected farm power: {}", expected_farm_power);
    println!(
        "Sum of expected turbine powers: {}",
        expected_turbine_powers.sum()
    );

    Ok(())
}
