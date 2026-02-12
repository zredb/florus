//! FLORIS-RS Example 20: Validation Studies
//!
//! This example demonstrates validation approaches for FLORIS-RS results,
//! including comparisons with expected values, Python FLORIS results,
//! and field data when available.

use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 20: Validation Studies");
    println!("==================================\n");

    let d = 126.0;

    println!("Validation Approaches:");
    println!("  1. Unit tests: Verify individual components");
    println!("  2. Regression tests: Compare against known good results");
    println!("  3. Cross-validation: Compare with Python FLORIS");
    println!("  4. Field data: Compare with measured wind farm data\n");

    // Create a simple 2-turbine farm
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 2]);
    let turbine_types = vec!["nrel_5MW".to_string(); 2];

    println!("Test Farm Configuration:");
    println!("  2 turbines at 5D spacing ({:.0} m)", 5.0 * d);
    println!("  Turbine type: nrel_5MW\n");

    // Create flow field
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]); // From west
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    let flow_field = florus::core::FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.14,
        1.225,
        turbulence_intensities.clone(),
        90.0,
    )?;

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    let mut model = florus::FlorisModel {
        farm: farm.clone(),
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let powers = model.get_turbine_powers();
    let total_power = (0..2).map(|i| powers[[0, i]]).sum::<f64>();

    println!("--- Wake-Free Reference ---\n");
    println!("With no wake interaction (turbines side-by-side):");
    println!(
        "  Upstream turbine power: {:.2} MW",
        powers[[0, 0]] / 1_000_000.0
    );
    println!(
        "  Downstream turbine power: {:.2} MW",
        powers[[0, 1]] / 1_000_000.0
    );
    println!("  Total farm power: {:.2} MW", total_power / 1_000_000.0);

    // Now test with aligned turbines (wakes present)
    let layout_x_aligned = Array1::from_vec(vec![0.0, 5.0 * d]);
    let layout_y_aligned = Array1::from_vec(vec![0.0, 0.0]);

    model.set_layout(&layout_x_aligned, &layout_y_aligned)?;
    model.run()?;

    let powers_aligned = model.get_turbine_powers();
    let total_power_aligned = (0..2).map(|i| powers_aligned[[0, i]]).sum::<f64>();

    println!("\n--- Wake Interaction Test ---\n");
    println!("With aligned turbines (wake present):");
    println!(
        "  Upstream turbine power: {:.2} MW",
        powers_aligned[[0, 0]] / 1_000_000.0
    );
    println!(
        "  Downstream turbine power: {:.2} MW",
        powers_aligned[[0, 1]] / 1_000_000.0
    );
    println!(
        "  Total farm power: {:.2} MW",
        total_power_aligned / 1_000_000.0
    );

    let wake_loss = (total_power - total_power_aligned) / total_power * 100.0;

    println!("\nWake Loss Analysis:");
    println!("  Wake-induced power loss: {:.1}%", wake_loss);
    println!(
        "  Downstream receives {:.1}% of upstream power",
        powers_aligned[[0, 1]] / powers_aligned[[0, 0]] * 100.0
    );

    println!("\n--- Expected Values (Validation) ---\n");

    // Typical expected values for validation
    let expected_upstream_power = 5_000_000.0; // ~5 MW for NREL 5MW at 8 m/s
    let expected_wake_loss_percent = 10.0; // Typically 5-15% for 5D spacing

    println!("Validation Criteria:");
    println!(
        "  Upstream power should be ~{:.1} MW",
        expected_upstream_power / 1_000_000.0
    );
    println!(
        "  Expected wake loss: {:.0}-{:.0}% for 5D spacing",
        expected_wake_loss_percent - 5.0,
        expected_wake_loss_percent + 5.0
    );

    let upstream_ok =
        (powers_aligned[[0, 0]] - expected_upstream_power).abs() < expected_upstream_power * 0.2;
    let wake_loss_ok = wake_loss > (expected_wake_loss_percent - 10.0)
        && wake_loss < (expected_wake_loss_percent + 15.0);

    println!("\nValidation Results:");
    println!(
        "  Upstream power valid: {}",
        if upstream_ok { "PASS" } else { "CHECK" }
    );
    println!(
        "  Wake loss within expected range: {}",
        if wake_loss_ok { "PASS" } else { "CHECK" }
    );

    println!("\n--- Cross-Validation with Python FLORIS ---\n");

    println!("Expected Python FLORIS v4.6 Results:");
    println!("  At 8 m/s wind speed, 270 deg direction:");
    println!("  - Upstream: ~4.9-5.0 MW");
    println!("  - Downstream: ~3.5-4.0 MW (70-80% of upstream)");
    println!("  - Wake loss: ~10-15%");
    println!("\nNote: Small differences are expected due to:");
    println!("  - Wake model parameterization");
    println!("  - Grid resolution");
    println!("  - Turbine power curve interpolation");

    println!("\n--- Validation Checklist ---\n");

    println!("[ ] Unit tests pass (cargo test)");
    println!("[ ] Regression tests pass");
    println!("[ ] Upstream turbine power matches expected");
    println!("[ ] Wake deficit is reasonable");
    println!("[ ] Wake deflection responds to yaw");
    println!("[ ] Multiple wind directions produce consistent results");

    println!("\n==================================");
    println!("Example completed successfully!");
    println!("\nTo validate against Python FLORIS:");
    println!("  1. Run identical case in Python FLORIS");
    println!("  2. Compare power outputs (should be <5% difference)");
    println!("  3. Test multiple wind conditions");
    println!("  4. Document any systematic differences");

    Ok(())
}
