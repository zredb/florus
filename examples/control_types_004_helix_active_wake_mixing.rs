//! Example: Helix active wake control
//!
//! This example runs a four-turbine farm with active wake control: the most
//! upstream turbine applies a helix yaw trajectory, which rotates its wake
//! and reduces the loads of the downstream turbines.
//!
//! Scenario aligned 1:1 with
//! examples_control_types/004_helix_active_wake_mixing.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    // A model with the awc-enabled turbines
    let mut fmodel = FlorisModel::from_file("examples/inputs/emgauss_helix.yaml")?;
    fmodel.set_operation_model("awc")?;

    // One wind condition, with the helix amplitude of the first turbine
    let n = 1;
    let awc_modes = ndarray::Array2::<String>::from_shape_vec(
        (n, 4),
        vec![
            "helix".to_string(),
            "baseline".to_string(),
            "baseline".to_string(),
            "baseline".to_string(),
        ],
    )?;
    let awc_amplitudes = Array2::from_shape_vec((n, 4), vec![2.5, 0.0, 0.0, 0.0])?;

    // A 4-turbine layout with a lateral offset of 3D and a streamwise offset
    // of 4D
    let d = 240.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 4.0 * d, 0.0, 4.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, -3.0 * d, -3.0 * d]),
    )?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.set_operation(
        None,
        None,
        Some(awc_modes),
        Some(awc_amplitudes),
        None,
        None,
    )?;

    fmodel.run()?;
    let turbine_powers = fmodel.get_turbine_powers()?;
    println!("Turbine powers [kW]:");
    for t in 0..turbine_powers.shape()[1] {
        println!("  turbine {t}: {:.1}", turbine_powers[[0, t]] / 1000.0);
    }

    // The flow fields of the baseline and the helix case
    let horizontal_plane = fmodel.calculate_horizontal_plane(150.0, 200, 100, None, None)?;
    let y_plane_baseline = fmodel.calculate_y_plane(200, 100, 0.0, None, None)?;
    let y_plane_helix = fmodel.calculate_y_plane(200, 100, -3.0 * d, None, None)?;
    let cross_plane = fmodel.calculate_cross_plane(100, 100, 4.0 * d, None, None)?;

    let u_range = |plane: &florus::core::cut_plane::CutPlane| {
        let u = &plane.data.u;
        (
            u.iter().copied().fold(f64::INFINITY, f64::min),
            u.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    println!("Flow field u ranges [m/s]:");
    println!(
        "  horizontal:        {:.3} .. {:.3}",
        u_range(&horizontal_plane).0,
        u_range(&horizontal_plane).1
    );
    println!(
        "  y plane baseline:  {:.3} .. {:.3}",
        u_range(&y_plane_baseline).0,
        u_range(&y_plane_baseline).1
    );
    println!(
        "  y plane helix:     {:.3} .. {:.3}",
        u_range(&y_plane_helix).0,
        u_range(&y_plane_helix).1
    );
    println!(
        "  cross plane:       {:.3} .. {:.3}",
        u_range(&cross_plane).0,
        u_range(&cross_plane).1
    );

    Ok(())
}
