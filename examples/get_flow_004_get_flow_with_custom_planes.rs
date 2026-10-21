//! Example: Plot velocity deficit profiles
//!
//! This example illustrates how to extract velocity deficit profiles at
//! several locations downstream of a turbine. The velocity deficit is defined
//! as `(homogeneous_wind_speed - u) / homogeneous_wind_speed`, where `u` is
//! the wake velocity obtained with a homogeneous inflow.
//!
//! Scenario aligned 1:1 with
//! examples_get_flow/004_plot_velocity_deficit_profiles.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    let d = 125.88; // The turbine diameter
    let hub_height = 90.0;
    let homogeneous_wind_speed = 8.0;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;

    // The profiles are sampled at three distances downstream of the turbine,
    // along lines parallel to the cross-stream direction
    let downstream_dists = [3.0 * d, 5.0 * d, 7.0 * d];

    let profiles = fmodel.sample_velocity_deficit_profiles(
        "cross-stream",
        Some(&downstream_dists),
        None,
        100,
        None,
        Some(homogeneous_wind_speed),
        None,
        0.0,
        0.0,
        Some(hub_height),
    )?;

    for (i, profile) in profiles.iter().enumerate() {
        println!(
            "Profile at {:.1} m downstream (x1/D = {:.3}), velocity deficit:",
            downstream_dists[i], profile.x1[0]
        );
        // Print every tenth point of the profile
        for j in (0..profile.velocity_deficit.len()).step_by(10) {
            println!(
                "  x2/D = {:+.3}, deficit = {:.6}",
                profile.x2[j], profile.velocity_deficit[j]
            );
        }
    }

    Ok(())
}
