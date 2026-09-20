//! Example: Visualize Y Cut Plane
//!
//! Demonstrates vertical plane visualization along the wind direction.
//!
//! Corresponds to: examples_visualizations/002_visualize_y_cut_plane.py

use florus::visualization::*;
use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Y Cut Plane Visualization ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
    )?;
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;
    fmodel.run()?;

    std::fs::create_dir_all("output/examples_visualizations")?;

    println!("Calculating Y-plane (streamwise vertical slice)...");
    match fmodel.calculate_y_plane(200, 100, 0.0) {
        Ok(cut_plane) => {
            println!("  Y-plane computed successfully");
            match visualize_cut_plane(
                &cut_plane,
                "output/examples_visualizations/y_cut_plane.png",
                Some(3.0),
                Some(9.0),
                "coolwarm",
                true,
                "Y Cut Plane",
            ) {
                Ok(_) => println!("  Saved to output/examples_visualizations/y_cut_plane.png"),
                Err(e) => eprintln!("  Visualization error: {}", e),
            }
        }
        Err(e) => eprintln!("  Y-plane calculation error: {}", e),
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
