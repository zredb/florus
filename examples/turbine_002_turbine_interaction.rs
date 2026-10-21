//! Example: Multiple turbine types
//!
//! This example runs a farm with two different turbine types and extracts the
//! flow planes around them.
//!
//! Scenario aligned 1:1 with
//! examples_turbine/002_multiple_turbine_types.py.

use florus::floris_model::FlorisModel;

fn main() -> anyhow::Result<()> {
    // A farm with two different turbine types
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multiple_turbine_types.yaml")?;
    fmodel.run()?;

    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 100)?;
    let y_plane = fmodel.calculate_y_plane(200, 100, 0.0)?;
    let cross_plane = fmodel.calculate_cross_plane(100, 100, 500.0)?;

    let u_range = |plane: &florus::core::cut_plane::CutPlane| {
        let u = &plane.data.u;
        (
            u.iter().copied().fold(f64::INFINITY, f64::min),
            u.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    println!("Flow field u ranges [m/s]:");
    println!(
        "  horizontal:   {:.3} .. {:.3}",
        u_range(&horizontal_plane).0,
        u_range(&horizontal_plane).1
    );
    println!(
        "  y plane:      {:.3} .. {:.3}",
        u_range(&y_plane).0,
        u_range(&y_plane).1
    );
    println!(
        "  cross plane:  {:.3} .. {:.3}",
        u_range(&cross_plane).0,
        u_range(&cross_plane).1
    );

    // The turbine powers of the two turbine types
    let powers = fmodel.get_turbine_powers()?;
    println!("Turbine powers [kW]:");
    for t in 0..powers.shape()[1] {
        println!("  turbine {t}: {:.1}", powers[[0, t]] / 1000.0);
    }

    Ok(())
}
