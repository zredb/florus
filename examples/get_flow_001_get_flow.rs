//! Example: Extract Wind Speed at Turbines
//!
//! This example demonstrates how to extract the wind speed at turbine points
//! from the FLORUS model. Both the u velocities and turbine average velocities
//! are extracted, then the turbine average is recalculated to show equivalence.
//!
//! Corresponds to: examples_get_flow/001_extract_wind_speed_at_turbines.py

use florus::{FlorisModel, Result};
use ndarray::Array;

fn main() -> Result<()> {
    println!("=== Extract Wind Speed at Turbines ===\n");

    // Initialize the FLORUS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Use default 3-turbine layout from gch.yaml
    println!("Using default 3-turbine layout from config...\n");

    // Calculate wake
    println!("Running simulation...");
    fmodel.run()?;

    // Collect the wind speed at all turbine points
    // u shape: [n_findex, n_turbines, grid_points_y, grid_points_z]
    let u_points = fmodel.core().flow_field.u.clone();

    println!("U points shape: {:?}", u_points.shape());
    println!("  Dimensions: [findex, turbines, grid_y, grid_z]");
    println!("  For turbine_grid_points=3, this is 1 x 3 x 3 x 3\n");

    // Get turbine average velocities
    let turbine_avg_velocities = fmodel.turbine_average_velocities();
    println!("Turbine average velocities shape: {:?}", turbine_avg_velocities.shape());
    println!("  Dimensions: [findex, turbines]\n");

    // Show that one is equivalent to the other following averaging
    println!("Verification:");
    println!("  Turbine average velocities are computed by taking the cube root");
    println!("  of the mean of cubed values across the rotor grid points.\n");

    // Recompute: cbrt(mean(u^3))
    let u_cubed = &u_points.mapv(|x| x.powi(3));
    
    // Mean over axes 2 and 3 (grid_y and grid_z)
    let mean_u_cubed_axis2 = u_cubed.mean_axis(ndarray::Axis(2)).unwrap();
    let mean_u_cubed = mean_u_cubed_axis2.mean_axis(ndarray::Axis(2)).unwrap();
    let recomputed = mean_u_cubed.mapv(|x: f64| x.cbrt());

    println!("  Original turbine_average_velocities:");
    for (i, val) in turbine_avg_velocities.iter().enumerate() {
        println!("    Turbine {}: {:.4} m/s", i, val);
    }

    println!("\n  Recomputed (cbrt(mean(u^3))):");
    for (i, val) in recomputed.iter().enumerate() {
        println!("    Turbine {}: {:.4} m/s", i, val);
    }

    // Verify they match
    let max_diff: f64 = turbine_avg_velocities
        .iter()
        .zip(recomputed.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, |max: f64, diff: f64| max.max(diff));

    println!("\n  Maximum difference: {:.2e} m/s", max_diff);
    if max_diff < 1e-10 {
        println!("  ✅ Verification passed! Values match perfectly.");
    } else {
        println!("  ⚠️  Warning: Values differ slightly (expected due to numerical precision)");
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
