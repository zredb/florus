//! Example: Setting the model parameters
//!
//! This example demonstrates setting each of the following:
//!     1) The wind conditions
//!     2) The wind farm layout
//!     3) The control settings
//!
//! Scenario aligned 1:1 with examples/004_set.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};
use florus::wind_data::{TimeSeries, WindRose};

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // ---------------------------------------------------------------
    // The atmospheric conditions
    // ---------------------------------------------------------------

    // The wind directions, wind speeds and turbulence intensities
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06]),
    )?;

    // The same wind conditions with a TimeSeries object
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06]),
    )?;
    fmodel.set_wind_data(&time_series)?;

    // The same wind conditions with a WindRose object
    let wind_rose = WindRose::new(
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array2::from_elem((1, 3), 0.06),
        None,
        None,
        false,
        None,
        None,
    )?;
    fmodel.set_wind_data(&wind_rose)?;

    // The wind shear
    fmodel.set_wind_shear(0.2)?;

    // The air density
    fmodel.set_air_density(1.1)?;

    // The reference wind height (the height at which the wind speed is given)
    fmodel.set_reference_wind_height(92.0)?;

    // ---------------------------------------------------------------
    // The layout
    // ---------------------------------------------------------------

    // A two-turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    // ---------------------------------------------------------------
    // The control settings
    // ---------------------------------------------------------------

    // The yaw angles: the matrix must match (n_findex, n_turbines)
    fmodel.set_yaw_angles(Array2::from_shape_vec(
        (3, 2),
        vec![0.0, 0.0, 25.0, 0.0, 0.0, 0.0],
    )?)?;

    // Clear the control signals
    fmodel.reset_operation();

    // Change to the mixed operation model
    fmodel.set_operation_model("mixed")?;

    // Shut down the front turbine for the first two wind conditions
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((3, 2), false);
    disable_turbines[[0, 0]] = true;
    disable_turbines[[1, 0]] = true;
    fmodel.set_operation(None, None, None, None, None, Some(disable_turbines))?;

    // Derate the front turbine for the first two wind conditions
    let rated_power = 5.0e6;
    let power_setpoints = Array2::from_shape_vec(
        (3, 2),
        vec![
            rated_power * 0.3,
            rated_power,
            rated_power * 0.3,
            rated_power,
            rated_power,
            rated_power,
        ],
    )?;

    // The power setpoints are set after the disabled turbines, so they take
    // precedence (as in the Python example, which uses two set() calls)
    fmodel.set_operation(None, Some(power_setpoints), None, None, None, None)?;

    fmodel.run()?;

    let powers = fmodel.get_turbine_powers()?;
    println!("Turbine powers [kW]:");
    for fi in 0..powers.shape()[0] {
        println!(
            "  condition {fi}: {:.1} {:.1}",
            powers[[fi, 0]] / 1000.0,
            powers[[fi, 1]] / 1000.0
        );
    }

    Ok(())
}
