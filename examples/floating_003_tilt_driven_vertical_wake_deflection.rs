//! Example: Tilt-driven vertical wake deflection
//!
//! This example initializes two FLORIS models: one with 5 degrees of tilt
//! (fixed across all wind speeds) and one with 15 degrees of tilt, and
//! compares the wakes and the powers of a five-turbine row.
//!
//! Scenario aligned 1:1 with
//! examples_floating/003_tilt_driven_vertical_wake_deflection.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::visualization::flow_visualization;

fn main() -> anyhow::Result<()> {
    let num_in_row = 5;
    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();

    for (input_file, tilt) in [
        ("examples/inputs/emgauss_floating_fixedtilt5.yaml", 5.0),
        ("examples/inputs/emgauss_floating_fixedtilt15.yaml", 15.0),
    ] {
        let mut fmodel = FlorisModel::from_file(input_file)?;
        let d = fmodel.core().farm.rotor_diameters()[0];

        // A row of five turbines
        fmodel.set_layout(
            &Array1::from_iter((0..num_in_row).map(|x| x as f64 * 5.0 * d)),
            &Array1::from_vec(vec![0.0; num_in_row]),
        )?;
        fmodel.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        fmodel.run()?;

        let powers = fmodel.get_turbine_powers()?;
        let turbine_powers: Vec<String> = (0..num_in_row)
            .map(|t| format!("{:.1}", powers[[0, t]] / 1000.0))
            .collect();
        println!(
            "Tilt {tilt:.0} deg, turbine powers [kW]: {}",
            turbine_powers.join(" ")
        );

        // The horizontal plane and the cross planes
        let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 100, None, None)?;
        flow_visualization::visualize_cut_plane(
            &horizontal_plane,
            format!("{output_dir}/003_horizontal_plane_tilt{}.png", tilt as i32),
            Some(4.0),
            Some(10.0),
            "viridis",
            true,
            &format!("Horizontal plane, {tilt:.0} deg tilt"),
        )?;

        for downstream_dist in [10.0, 1200.0, 2500.0] {
            let cross_plane =
                fmodel.calculate_cross_plane(100, 100, downstream_dist, None, None)?;
            flow_visualization::visualize_cut_plane(
                &cross_plane,
                format!(
                    "{output_dir}/003_cross_plane_tilt{}_{:.0}m.png",
                    tilt as i32, downstream_dist
                ),
                Some(4.0),
                Some(10.0),
                "viridis",
                true,
                &format!("Cross plane at {downstream_dist:.0} m, {tilt:.0} deg tilt"),
            )?;
        }
    }

    Ok(())
}
