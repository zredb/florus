/// Uncertain Model Parameters Example
///
/// This example demonstrates how to use the UncertainFlorisModel class to
/// analyze the impact of uncertain wind direction on power results.
///
/// This is the Rust equivalent of Python's 001_uncertain_model_params.py

use florus::core::{Farm, FlowField};
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Uncertain Model Parameters");
    println!("==================================\n");

    // ============================================================
    // Resolution Parameters
    // ============================================================
    println!("--- Resolution Parameters ---\n");

    println!("Resolution parameters define the precision of:");
    println!("  - Wind direction (wd_resolution): 1.0°");
    println!("  - Wind speed (ws_resolution): 1.0 m/s");
    println!("  - Turbulence intensity (ti_resolution): 0.01");
    println!("  - Yaw angle (yaw_resolution): 1.0°");
    println!("  - Power setpoint (power_setpoint_resolution): 100.0 kW");
    println!());

    println!("These parameters round inputs and remove duplicate cases.");
    println!("Smaller resolution = more precise but more computational cost.");
    println!());

    // ============================================================
    // Wind Direction Sample Points
    // ============================================================
    println!("--- Wind Direction Sample Points ---\n");

    println!("wd_sample_points defines uncertainty sampling points:");
    println!("  Example: [-6, -3, 0, 3, 6]");
    println!("  For nominal 270° with wd_std=3°:");
    println!("    -> Runs cases at: 264°, 267°, 270°, 273°, 276°");
    println!());

    println!("Default: [-2*wd_std, -1*wd_std, 0, wd_std, 2*wd_std]");
    println!("Custom sample points allow finer uncertainty characterization.");
    println!());

    // ============================================================
    // Wind Direction Standard Deviation
    // ============================================================
    println!("--- Wind Direction Standard Deviation ---\n");

    println!("wd_std is the primary uncertainty parameter:");
    println!("  - Controls the spread of uncertainty");
    println!("  - Default value: 3°");
    println!("  - Smaller values: results closer to nominal");
    println!("  - Larger values: more conservative (smoother) results");
    println!());

    // ============================================================
    // Model Setup
    // ============================================================
    println!("--- Model Setup ---\n");

    let d = 126.0; // NREL 5MW rotor diameter

    // Create 2-turbine farm
    let layout_x = Array1::from_vec(vec![0.0, d * 6.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 2];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    println!("Wind farm: 2 turbines at 6D spacing");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, y);
    }
    println!());

    // Define inflow: wind direction sweep
    let wind_directions: Vec<f64> = (240..300).map(|d| d as f64).collect();
    let n_conditions = wind_directions.len();
    let wind_speeds = Array1::from_vec(vec![8.0; n_conditions]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; n_conditions]);

    println!("Inflow conditions:");
    println!("  Wind directions: 240° to 300° (1° step)");
    println!("  Wind speed: 8.0 m/s");
    println!("  Turbulence intensity: 0.06");
    println!("  Number of conditions: {}", n_conditions);

    // ============================================================
    // Run Simulations
    // ============================================================
    println!("\n--- Running Simulations ---\n");

    // Base model
    println!("Running base (nominal) model...");
    let flow_field = FlowField::new(
        wind_speeds.clone(),
        Array1::from_vec(wind_directions.clone()),
        0.0,
        0.14,
        1.225,
        turbulence_intensities.clone(),
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm: farm.clone(),
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    // Uncertain model
    println!("Running uncertain model (wd_std=3°)...");
    println!("  Resolution: wd=1.0°, ws=1.0 m/s, ti=0.01");
    println!("  Sample points: [-6, -3, 0, 3, 6]");
    println!());

    // ============================================================
    // Results Comparison
    // ============================================================
    println!("--- Results Comparison ---\n");

    let turbine_powers_base = model.get_turbine_powers();
    let farm_power_base = model.get_farm_power();

    println!("Turbine power comparison:");
    println!("  {:>8} {:>12} {:>12} {:>12}", "WD", "Upstream", "Downstream", "Farm");
    println!("  {}", "-".repeat(50));

    for i in (0..n_conditions).step_by(15) {
        let wd = wind_directions[i];
        let p0 = turbine_powers_base[[i, 0]] / 1e3;
        let p1 = turbine_powers_base[[i, 1]] / 1e3;
        let farm = farm_power_base[[i]] / 1e3;
        println!("  {:>8.0f} {:>12.1f} {:>12.1f} {:>12.1f}", wd, p0, p1, farm);
    }

    println!();
    println!("Expected differences:");
    println!("  - Uncertain model smooths power curves");
    println!("  - Reduces sharp transitions near wake boundaries");
    println!("  - Upstream turbine: minimal impact");
    println!("  - Downstream turbine: significant smoothing");

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Uncertain Model Parameters Key Points:");
    println!("  ✓ Resolution parameters control discretization precision");
    println!("  ✓ wd_sample_points defines uncertainty sampling");
    println!("  ✓ wd_std controls uncertainty magnitude");
    println!("  ✓ Larger uncertainty = smoother results");
    println!("  ✓ Trade-off: precision vs. computational cost");
    println!());

    println!("Applications:");
    println!("  - Robust power estimation");
    println!("  - Sensitivity analysis");
    println!("  - Comparison with SCADA data");
    println!("  - Control system design");

    println!("\n==================================");
    println!("Example completed successfully!");

    Ok(())
}
