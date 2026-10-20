//! Example: Gridded layout design
//!
//! This example shows a layout optimization that places as many turbines as
//! possible into a given boundary using a gridded layout pattern.
//!
//! Scenario aligned 1:1 with
//! examples_layout_optimization/004_generate_gridded_layout.py (5D spacing,
//! 1000 x 1000 m boundary).

use florus::floris_model::FlorisModel;
use florus::optimization::layout_optimization::layout_optimization_base::Boundary;
use florus::optimization::layout_optimization::LayoutOptimizationGridded;

fn main() -> anyhow::Result<()> {
    // Load the FLORIS model
    let fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // The boundaries for the turbines, specified as vertices
    let boundaries = Boundary::Polygon(vec![
        (0.0, 0.0),
        (0.0, 1000.0),
        (1000.0, 1000.0),
        (1000.0, 0.0),
        (0.0, 0.0),
    ]);

    // Set up the optimization object with 5D spacing
    // (results in spacing of 5 * 125.88 = 629.4 m)
    let mut layout_opt = LayoutOptimizationGridded::new(&fmodel, boundaries, Some(5.0), None)?;
    layout_opt.optimize();

    println!(
        "Gridded layout complete: {} turbines in bounds.",
        layout_opt.n_turbines_max
    );
    println!("x: {:?}", layout_opt.x_opt);
    println!("y: {:?}", layout_opt.y_opt);

    Ok(())
}
