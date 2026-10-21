//! Example: Derating control with the simple-derating operation model
//!
//! This example demonstrates the simple-derating control model: the power
//! setpoint of every turbine derates its power output to at most the
//! setpoint. The setpoint is swept from 1 W to 6 MW and the resulting
//! turbine powers are reported.
//!
//! Scenario aligned 1:1 with
//! examples_control_types/001_derating_control.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change to the simple-derating model turbine
    fmodel.set_operation_model("simple-derating")?;

    // A simple two turbine layout with derating turbines
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    // Constant wind conditions over n_findex = N time steps
    let n = 50;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![10.0; n]),
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.run()?;
    let powers_orig = fmodel.get_turbine_powers()?;

    // Add a derating level to both turbines: the setpoint is swept from 1 W
    // to 6 MW
    let setpoints: Vec<f64> = (0..n)
        .map(|i| 1.0 + (6.0e6 - 1.0) * i as f64 / (n - 1) as f64)
        .collect();
    let power_setpoints = Array2::from_shape_vec((n, 2), {
        let mut values = Vec::with_capacity(2 * n);
        for s in &setpoints {
            values.push(*s);
            values.push(*s);
        }
        values
    })?;
    fmodel.set_operation(None, Some(power_setpoints), None, None, None, None)?;
    fmodel.run()?;
    let powers_derated = fmodel.get_turbine_powers()?;

    // The downstream turbine is only derated (the upstream one keeps its
    // setpoint at the default value)
    let mut power_setpoints_2 = Array2::from_elem((n, 2), f64::NAN);
    for (i, s) in setpoints.iter().enumerate() {
        power_setpoints_2[[i, 0]] = *s;
    }
    fmodel.set_operation(None, Some(power_setpoints_2), None, None, None, None)?;
    fmodel.run()?;
    let powers_avail_ds = fmodel.get_turbine_powers()?;

    println!("Setpoint [kW], upstream power [kW], downstream power [kW], derated upstream [kW], derated downstream [kW], available downstream [kW]");
    for i in 0..n {
        println!(
            "{:.1}, {:.1}, {:.1}, {:.1}, {:.1}, {:.1}",
            setpoints[i] / 1000.0,
            powers_orig[[i, 0]] / 1000.0,
            powers_orig[[i, 1]] / 1000.0,
            powers_derated[[i, 0]] / 1000.0,
            powers_derated[[i, 1]] / 1000.0,
            powers_avail_ds[[i, 1]] / 1000.0
        );
    }

    Ok(())
}
