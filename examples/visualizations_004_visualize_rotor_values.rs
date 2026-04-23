//! Example: Visualize Rotor Values (from examples_visualizations folder)
//!
//! Demonstrate visualizing the flow velocities at the rotor using plot_rotor_values.
//! This corresponds to `examples_visualizations/004_visualize_rotor_values.py` in Python FLORIS.

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== Rotor Values Visualization Example ===\n");
    println!("This example demonstrates visualizing flow velocities at turbine rotors.\n");

    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 2 turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    // Create output directory
    std::fs::create_dir_all("output/examples_visualizations")?;

    println!("Running simulation...");
    fmodel.run()?;

    // Get the velocity field from the core
    let core = fmodel.core();
    let u_field = &core.flow_field.u;
    
    println!("\nVelocity field dimensions: {:?}", u_field.dim());
    println!("  - findex: {}", u_field.dim().0);
    println!("  - turbines: {}", u_field.dim().1);
    println!("  - y grid points: {}", u_field.dim().2);
    println!("  - z grid points: {}", u_field.dim().3);

    // Plot the values at each rotor
    println!("\nCreating rotor plane visualization...");
    match plot_rotor_values(
        u_field,
        0,           // findex
        1,           // n_rows
        2,           // n_cols
        "output/examples_visualizations/rotor_values.png",
        "coolwarm",  // colormap
    ) {
        Ok(_) => {
            println!("   ✓ Saved to output/examples_visualizations/rotor_values.png");
            println!("   ℹ Shows velocity distribution across each turbine's rotor disk");
            println!("   ℹ Left: Turbine 1 (upstream), Right: Turbine 2 (downstream)");
        }
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    println!("\n=== Rotor Values Example Complete ===");
    println!("\nKey observations:");
    println!("  - Upstream turbine shows uniform velocity (~8 m/s)");
    println!("  - Downstream turbine shows wake deficit (lower velocities)");
    println!("  - Velocity patterns reveal wake structure and recovery");

    Ok(())
}
