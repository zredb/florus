//! Example: Specify turbine power curve
//!
//! This example demonstrates how to specify a turbine model based on a power
//! and thrust curve for the wind turbine, as well as possible physical
//! parameters (which default to the parameters of the NREL 5MW reference
//! turbine).
//!
//! The example builds a turbine definition at runtime (Python
//! `floris.turbine_library.build_cosine_loss_turbine_dict`), replaces the
//! turbine of the model with it and compares the power of the model over a
//! sweep of the wind speed with the specified power curve. The comparison is
//! written as a scatter plot (like the Python example) and printed as a
//! table.
//!
//! Scenario aligned 1:1 with examples_turbine/003_specify_turbine_power_curve.py.

use florus::core::turbines::turbine_type::TurbineType;
use florus::floris_model::FlorisModel;
use florus::types::Array1;
use plotters::prelude::*;

fn main() -> anyhow::Result<()> {
    // Generate an example turbine power and thrust curve for use in the FLORIS model
    let powers_orig: [f64; 11] = [
        0.0, 30.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0,
    ];
    let wind_speeds: [f64; 11] = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
    let power_coefficients: Vec<f64> = std::iter::once(0.0)
        .chain(
            wind_speeds
                .iter()
                .zip(powers_orig.iter())
                .skip(1)
                .map(|(&ws, &power)| {
                    power
                        / (0.5 * 126.0_f64.powi(2) * std::f64::consts::PI / 4.0
                            * 1.225
                            * ws.powi(3))
                }),
        )
        .collect();
    let thrust_coefficients: [f64; 11] = [0.0, 0.9, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2];

    let turbine = TurbineType::build_cosine_loss(
        "example_turbine",
        &wind_speeds,
        None,
        Some(&power_coefficients),
        None,
        Some(&thrust_coefficients),
        Some(1.0),
        90.0,
        126.0,
        8.0,
        1.225,
        5.0,
        1.88,
        1.88,
    )?;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Replace the turbine(s) in the FLORIS model with the created one
    fmodel.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    fmodel.set_turbine_type(&turbine)?;
    fmodel.assign_hub_height_to_ref_height()?;

    // The wind conditions of the sweep
    let n = 100;
    let sweep_wind_speeds: Vec<f64> = (0..n)
        .map(|i| 1.0 + 14.0 * i as f64 / (n - 1) as f64)
        .collect();
    fmodel.set_wind_conditions(
        Array1::from_vec(sweep_wind_speeds.clone()),
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.run()?;

    let powers = fmodel.get_farm_power()?;

    // The powers of the specified power curve (the power coefficient of the
    // curve times the kinetic power of the wind at the reference air density)
    let specified_powers: Vec<f64> = power_coefficients
        .iter()
        .zip(wind_speeds.iter())
        .map(|(&cp, &ws)| {
            cp * 0.5 * 1.225 * 126.0_f64.powi(2) * std::f64::consts::PI / 4.0 * ws.powi(3) / 1000.0
        })
        .collect();

    println!("Specified points (the power curve of the turbine)");
    println!("wind speed [m/s], power [kW]");
    for (ws, power) in wind_speeds.iter().zip(specified_powers.iter()) {
        println!("{:.1}, {:.6}", ws, power);
    }

    println!();
    println!("Test points (the power of the turbine in the model)");
    println!("wind speed [m/s], power [kW]");
    for (ws, power) in sweep_wind_speeds.iter().zip(powers.iter()) {
        println!("{:.10}, {:.10}", ws, power / 1000.0);
    }

    // The scatter plot of the two series (Python's figure)
    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();
    let output_path = format!("{output_dir}/003_turbine_power_curve.png");

    let max_wind_speed = wind_speeds
        .iter()
        .copied()
        .chain(sweep_wind_speeds.iter().copied())
        .fold(0.0_f64, f64::max);
    let max_power = specified_powers
        .iter()
        .copied()
        .chain(powers.iter().map(|p| p / 1000.0))
        .fold(0.0_f64, f64::max);

    {
        let root = BitMapBackend::new(&output_path, (640, 480)).into_drawing_area();
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(&root)
            .caption(
                "Power curve of the specified turbine",
                ("sans-serif", 18).into_font(),
            )
            .margin(10)
            .x_label_area_size(50)
            .y_label_area_size(60)
            .build_cartesian_2d(0.0_f64..(max_wind_speed * 1.05), 0.0_f64..(max_power * 1.1))?;

        chart
            .configure_mesh()
            .x_desc("Wind speed [m/s]")
            .y_desc("Power [kW]")
            .label_style(("sans-serif", 12).into_font())
            .axis_desc_style(("sans-serif", 14).into_font())
            .draw()?;

        // The test points (Python: color "C0", s=5)
        chart
            .draw_series(
                sweep_wind_speeds
                    .iter()
                    .zip(powers.iter())
                    .map(|(&ws, &power)| Circle::new((ws, power / 1000.0), 3, BLUE.filled())),
            )?
            .label("Test points")
            .legend(|c| Circle::new(c, 4, BLUE.filled()));

        // The specified points (Python: color "red", s=20)
        chart
            .draw_series(
                wind_speeds
                    .iter()
                    .zip(specified_powers.iter())
                    .map(|(&ws, &power)| Circle::new((ws, power), 6, RED.filled())),
            )?
            .label("Specified points")
            .legend(|c| Circle::new(c, 6, RED.filled()));

        chart
            .configure_series_labels()
            .background_style(WHITE.filled())
            .border_style(BLACK)
            .label_font(("sans-serif", 12).into_font())
            .draw()?;
    }

    println!();
    println!("The plot was written to {output_path}");

    Ok(())
}
