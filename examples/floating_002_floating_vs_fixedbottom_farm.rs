//! Example: Floating vs fixed-bottom farm
//!
//! This example demonstrates the impact of floating on the turbine power and
//! the wake behaviour: a 20-turbine gridded wind farm is simulated with the
//! empirical Gaussian wake model, once with a fixed-bottom turbine and once
//! with a floating turbine whose tilt varies with the wind speed.
//!
//! Scenario aligned 1:1 with
//! examples_floating/002_floating_vs_fixedbottom_farm.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::WindRose;

fn main() -> anyhow::Result<()> {
    // The 5 x 4 grid of the farm
    let mut layout_x: Vec<f64> = Vec::new();
    let mut layout_y: Vec<f64> = Vec::new();
    for j in 0..4 {
        for i in 0..5 {
            layout_x.push(i as f64 * 630.0);
            layout_y.push(j as f64 * 630.0);
        }
    }

    let mut powers_fixed: Vec<Vec<f64>> = Vec::new();
    let mut powers_floating: Vec<Vec<f64>> = Vec::new();

    for (input_file, powers) in [
        ("examples/inputs/emgauss_fixed.yaml", &mut powers_fixed),
        (
            "examples/inputs/emgauss_floating.yaml",
            &mut powers_floating,
        ),
    ] {
        let mut fmodel = FlorisModel::from_file(input_file)?;
        fmodel.set_layout(
            &Array1::from_vec(layout_x.clone()),
            &Array1::from_vec(layout_y.clone()),
        )?;
        fmodel.set_wind_conditions(
            Array1::from_vec(vec![10.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        fmodel.run()?;

        let turbine_powers = fmodel.get_turbine_powers()?;
        let row: Vec<f64> = (0..turbine_powers.shape()[1])
            .map(|t| turbine_powers[[0, t]])
            .collect();
        powers.push(row);
    }

    // The difference between the floating and the fixed-bottom farm
    println!("Turbine positions and powers:");
    println!(
        "{:>8} {:>8} {:>14} {:>14} {:>12}",
        "x [m]", "y [m]", "fixed [kW]", "floating [kW]", "difference [kW]"
    );
    for t in 0..powers_fixed[0].len() {
        let difference = powers_floating[0][t] - powers_fixed[0][t];
        println!(
            "{:>8.0} {:>8.0} {:>14.1} {:>14.1} {:>12.1}",
            layout_x[t],
            layout_y[t],
            powers_fixed[0][t] / 1000.0,
            powers_floating[0][t] / 1000.0,
            difference / 1000.0
        );
    }

    let sum_fixed: f64 = powers_fixed[0].iter().sum();
    let sum_floating: f64 = powers_floating[0].iter().sum();
    println!("Farm power fixed bottom: {:.1} kW", sum_fixed / 1000.0);
    println!("Farm power floating:     {:.1} kW", sum_floating / 1000.0);
    println!(
        "Difference:              {:.1} kW ({:.2}%)",
        (sum_floating - sum_fixed) / 1000.0,
        100.0 * (sum_floating - sum_fixed) / sum_fixed
    );

    // Compute the AEP with the wind rose of the example
    let wind_rose = WindRose::read_csv_long(
        "examples/inputs/wind_rose.csv",
        "ws",
        "wd",
        Some("freq_val"),
        0.06.into(),
        ',',
    )?;

    let mut aep_fixed = 0.0;
    let mut aep_floating = 0.0;
    for (input_file, aep) in [
        ("examples/inputs/emgauss_fixed.yaml", &mut aep_fixed),
        ("examples/inputs/emgauss_floating.yaml", &mut aep_floating),
    ] {
        let mut fmodel = FlorisModel::from_file(input_file)?;
        fmodel.set_layout(
            &Array1::from_vec(layout_x.clone()),
            &Array1::from_vec(layout_y.clone()),
        )?;
        fmodel.set_wind_data(&wind_rose)?;
        fmodel.run()?;
        *aep = fmodel.get_farm_aep(&wind_rose, 8760.0)?;
    }

    println!("Farm AEP (fixed bottom): {:.3} GWh", aep_fixed / 1.0e9);
    println!("Farm AEP (floating): {:.3} GWh", aep_floating / 1.0e9);
    println!(
        "Floating AEP increase: {:.3} GWh ({:.2}%)",
        (aep_floating - aep_fixed) / 1.0e9,
        100.0 * (aep_floating - aep_fixed) / aep_fixed
    );

    Ok(())
}
