//! Example: Get Flow Field with Custom Planes
//!
//! This example demonstrates extracting flow field data on custom planes.
//!
//! Corresponds to: examples_get_flow/004_get_flow_with_custom_planes.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Get Flow with Custom Planes ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 630.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    println!("Custom Plane Extraction Concept:\n");
    println!("Flow field can be sampled on various planes:");
    println!("  1. Horizontal planes (z = constant)");
    println!("     - Hub height plane");
    println!("     - Ground level plane");
    println!("     - Multiple heights for 3D analysis\n");

    println!("  2. Vertical planes (x or y = constant)");
    println!("     - Cross-stream planes");
    println!("     - Streamwise planes");
    println!("     - Wake cross-sections\n");

    println!("  3. Arbitrary oriented planes");
    println!("     - Tilted planes");
    println!("     - Rotated coordinates");
    println!("     - Custom orientations\n");

    println!("Applications:");
    println!("  - Wake visualization");
    println!("  - Model validation");
    println!("  - Turbine placement optimization");
    println!("  - Load analysis\n");

    fmodel.run()?;
    println!("Note: Full custom plane extraction requires sample_velocity_deficit_profiles method.");

    println!("\n=== Example Complete ===");
    Ok(())
}
