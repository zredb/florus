//! Example: Heterogeneous inflow in 2D and 3D
//!
//! This example runs a model whose heterogeneous inflow is defined in the
//! input file (`gch_heterogeneous_inflow.yaml`) and inspects the resulting
//! flow field: the horizontal plane, the streamwise profile and the
//! spanwise profile downstream of the farm.
//!
//! Scenario aligned 1:1 with
//! examples_heterogeneous/004_heterogeneous_2d_and_3d.py.

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::HeterogeneousInflowConfig;
use florus::types::{Array1, Array2};

/// Minimum value of an array (ndarray is built without the `std` feature, so
/// `ArrayBase::min` is unavailable).
fn min_of(a: &Array1) -> f64 {
    a.iter().copied().fold(f64::INFINITY, f64::min)
}

/// Maximum value of an array.
fn max_of(a: &Array1) -> f64 {
    a.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Mean value of an array.
fn mean_of(a: &Array1) -> f64 {
    if a.is_empty() {
        return f64::NAN;
    }
    a.iter().sum::<f64>() / a.len() as f64
}

fn main() -> anyhow::Result<()> {
    // The heterogeneous inflow is defined in the input file
    let mut fmodel_2d = FlorisModel::from_file("examples/inputs/gch_heterogeneous_inflow.yaml")?;

    // Set shear to 0.0 to highlight the heterogeneous inflow
    fmodel_2d.set_wind_shear(0.0)?;

    // Run FLORIS and inspect the flow field
    fmodel_2d.run()?;

    let horizontal_plane = fmodel_2d.calculate_horizontal_plane(90.0, 200, 100, None, None)?;
    let y_plane = fmodel_2d.calculate_y_plane(200, 100, 0.0, None, None)?;
    let cross_plane = fmodel_2d.calculate_cross_plane(100, 100, 500.0, None, None)?;

    // Report the heterogeneous inflow effect on the horizontal plane
    let hp = &horizontal_plane.data;
    println!("2D heterogeneous inflow (from the input file):");
    println!("Horizontal plane at hub height:");
    println!(
        "  x1 range: {:.1} .. {:.1} m",
        min_of(&hp.x1),
        max_of(&hp.x1)
    );
    println!(
        "  x2 range: {:.1} .. {:.1} m",
        min_of(&hp.x2),
        max_of(&hp.x2)
    );
    println!(
        "  u range: {:.3} .. {:.3} m/s",
        min_of(&hp.u),
        max_of(&hp.u)
    );
    println!("  mean u: {:.3} m/s", mean_of(&hp.u));

    let yp = &y_plane.data;
    println!("Streamwise profile (y = 0):");
    println!(
        "  u range: {:.3} .. {:.3} m/s",
        min_of(&yp.u),
        max_of(&yp.u)
    );

    let cp = &cross_plane.data;
    println!("Spanwise profile (500 m downstream):");
    println!(
        "  u range: {:.3} .. {:.3} m/s",
        min_of(&cp.u),
        max_of(&cp.u)
    );

    // The 3-dimensional case: the speed-ups are defined at (x, y, z)
    // locations, so the multiplier varies with height.
    let mut fmodel_3d = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_3d.set_wind_shear(0.0)?;

    let speed_multipliers =
        Array2::from_shape_vec((1, 8), vec![1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0])?;
    fmodel_3d.set_heterogeneous_inflow_config(HeterogeneousInflowConfig {
        x: Array1::from_vec(vec![
            -300.0, -300.0, -300.0, -300.0, 2600.0, 2600.0, 2600.0, 2600.0,
        ]),
        y: Array1::from_vec(vec![
            -300.0, 300.0, -300.0, 300.0, -300.0, 300.0, -300.0, 300.0,
        ]),
        z: Some(Array1::from_vec(vec![
            540.0, 540.0, 0.0, 0.0, 540.0, 540.0, 0.0, 0.0,
        ])),
        wind_speeds: None,
        wind_directions: None,
        speed_multipliers,
    })?;
    fmodel_3d.run()?;

    let horizontal_plane_3d = fmodel_3d.calculate_horizontal_plane(90.0, 200, 100, None, None)?;
    let y_plane_3d = fmodel_3d.calculate_y_plane(200, 100, 0.0, None, None)?;
    let cross_plane_3d = fmodel_3d.calculate_cross_plane(100, 100, 500.0, None, None)?;

    let hp3 = &horizontal_plane_3d.data;
    println!("3D heterogeneous inflow:");
    println!("Horizontal plane at hub height:");
    println!(
        "  u range: {:.3} .. {:.3} m/s",
        min_of(&hp3.u),
        max_of(&hp3.u)
    );
    println!("  mean u: {:.3} m/s", mean_of(&hp3.u));

    let yp3 = &y_plane_3d.data;
    println!("Streamwise profile (y = 0):");
    println!(
        "  u range: {:.3} .. {:.3} m/s",
        min_of(&yp3.u),
        max_of(&yp3.u)
    );

    let cp3 = &cross_plane_3d.data;
    println!("Spanwise profile (500 m downstream):");
    println!(
        "  u range: {:.3} .. {:.3} m/s",
        min_of(&cp3.u),
        max_of(&cp3.u)
    );

    Ok(())
}
