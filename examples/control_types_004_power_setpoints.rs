//! Example: Helix Active Wake Mixing
//!
//! Demonstrates helix AWC on upstream turbines in a 4-turbine farm.
//!
//! Corresponds to: examples_control_types/004_helix_active_wake_mixing.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Helix Active Wake Mixing Example ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/emgauss_helix.yaml")?;
    fmodel.set_operation_model("awc")?;

    let d = 240.0;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 4.0 * d, 0.0, 4.0 * d]),
        &ndarray::arr1(&[0.0, 0.0, -3.0 * d, -3.0 * d]),
    )?;
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    // Single findex with helix on T0
    let mut awc_modes = ndarray::Array2::from_elem((1, 4), "baseline".to_string());
    awc_modes[[0, 0]] = "helix".to_string();
    let mut awc_amplitudes = ndarray::Array2::zeros((1, 4));
    awc_amplitudes[[0, 0]] = 2.5;

    fmodel.set_operation(None, None, Some(awc_modes), Some(awc_amplitudes), None, None)?;
    fmodel.run()?;

    let powers = fmodel.get_turbine_powers() / 1000.0;
    println!("With helix AWC (amplitude=2.5 deg on T0):");
    println!("  T0: {:.2} kW", powers[[0, 0]]);
    println!("  T1: {:.2} kW", powers[[0, 1]]);
    println!("  T2: {:.2} kW", powers[[0, 2]]);
    println!("  T3: {:.2} kW", powers[[0, 3]]);
    println!("  Farm: {:.2} kW", powers.row(0).sum() / 1000.0);

    // Sweep helix amplitudes
    let n = 10;
    println!("\nHelix amplitude sweep (T0):");
    println!("{:<20} | {:<15} | {:<15} | {:<15}", "Amplitude (deg)", "T0 Power (kW)", "Farm (kW)", "T1 Power (kW)");
    println!("{}", "-".repeat(68));

    for i in 0..n {
        let amp = i as f64 * 5.0 / (n - 1) as f64;
        let mut modes = ndarray::Array2::from_elem((1, 4), "baseline".to_string());
        modes[[0, 0]] = "helix".to_string();
        let mut amps = ndarray::Array2::zeros((1, 4));
        amps[[0, 0]] = amp;

        fmodel.set_operation(None, None, Some(modes), Some(amps), None, None)?;
        fmodel.run()?;
        let p = fmodel.get_turbine_powers() / 1000.0;
        println!(
            "{:<20.1} | {:<15.2} | {:<15.2} | {:<15.2}",
            amp, p[[0, 0]], p.row(0).sum(), p[[0, 1]]
        );
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
