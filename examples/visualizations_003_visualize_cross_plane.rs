//! Example: Visualize Cross Plane
//!
//! Demonstrates vertical plane visualization perpendicular to the wind direction.
//!
//! Corresponds to: examples_visualizations/003_visualize_cross_plane.py

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== Cross Plane Visualization ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0]),
        &ndarray::arr1(&[0.0]),
    )?;
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;
    fmodel.run()?;

    std::fs::create_dir_all("output/examples_visualizations")?;

    println!("Calculating cross-plane (perpendicular to wind at 500m downstream)...");
    match fmodel.calculate_cross_plane(100, 100, 500.0) {
        Ok(cut_plane) => {
            println!("  Cross-plane computed successfully");
            match visualize_cut_plane(
                &cut_plane,
                "output/examples_visualizations/cross_plane.png",
                Some(3.0),
                Some(9.0),
                "coolwarm",
                true,
                "Cross Plane at 500m",
            ) {
                Ok(_) => println!("  Saved to output/examples_visualizations/cross_plane.png"),
                Err(e) => eprintln!("  Visualization error: {}", e),
            }
        }
        Err(e) => eprintln!("  Cross-plane calculation error: {}", e),
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
