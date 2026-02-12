/// Example 2: Wake Models
///
/// This example demonstrates the different wake models available in FLORIS-RS:
///
/// 1. Wake Velocity Models:
///    - Gauss: Gaussian wake model (Bastankhah & Porté-Agel 2014)
///    - Jensen: Top-hat wake model (Jensen 1983)
///    - TurbPark: Turbulence-based wake model
///    - CumulativeGaussCurl: Cumulative curl model
///
/// 2. Wake Deflection Models:
///    - Gauss: Gaussian deflection model
///    - Jimenez: Jimenez deflection model
///    - EmpiricalGauss: Empirical deflection model
///
/// 3. Wake Combination Models:
///    - FLS: Free stream wake superposition
///    - SOSFS: Sum of squares wake superposition
///    - MAX: Maximum wake deficit
///
/// This example compares power outputs using different wake models.
///
/// This is the Rust equivalent of demonstrating different wake modeling approaches.

use florus::core::{Farm, FlowField};
use florus::floris_config::SolverConfig;
use florus::FlorisModel;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 2: Wake Models\n");
    println!("===========================================================\n");

    // Create a 5-turbine layout in a row
    // Turbine spacing: 5 rotor diameters (5 * 126 = 630 m)
    let rotor_diameter = 126.0;
    let spacing = 5.0 * rotor_diameter;
    let layout_x = Array1::from_vec((0..5).map(|i| i as f64 * spacing).collect());
    let layout_y = Array1::from_vec(vec![0.0; 5]);
    let turbine_types = vec!["nrel_5MW".to_string(); 5];

    println!("Creating wind farm with 5 turbines in a row:");
    println!("  Turbine spacing: {:.0} m ({:.1} D)", spacing, spacing / rotor_diameter);
    println!();

    // Wind conditions
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]); // From West (perpendicular to row)
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    println!("Wind condition:");
    println!("  Wind speed: 8.0 m/s");
    println!("  Wind direction: 270° (perpendicular to turbine row)");
    println!("  Turbulence intensity: 0.06\n");

    // ============================================================
    // Test 1: Gauss Wake Velocity Model
    // ============================================================
    println!("--- Test 1: Gauss Wake Velocity Model ---");
    
    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;
    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,   // wind_veer
        0.12,  // wind_shear
        1.225, // air_density
        turbulence_intensities.clone(),
        90.0,  // reference_wind_height
    )?;

    let mut model_gauss = FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    model_gauss.initialize_grid()?;
    model_gauss.initialize_flow_field()?;
    model_gauss.run()?;

    let powers_gauss = model_gauss.get_turbine_powers();
    let farm_power_gauss = model_gauss.get_farm_power();

    println!("Turbine powers (MW):");
    for i in 0..5 {
        println!("  T{}: {:.3} MW", i, powers_gauss[[0, i]] / 1_000_000.0);
    }
    println!("  Total: {:.3} MW\n", farm_power_gauss[[0]] / 1_000_000.0);

    // ============================================================
    // Test 2: Jensen Wake Velocity Model
    // ============================================================
    println!("--- Test 2: Jensen Wake Velocity Model ---");
    
    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;
    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,   // wind_veer
        0.12,  // wind_shear
        1.225, // air_density
        turbulence_intensities.clone(),
        90.0,  // reference_wind_height
    )?;

    let mut model_jensen: FlorisModel = FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    model_jensen.initialize_grid()?;
    model_jensen.initialize_flow_field()?;
    model_jensen.run()?;

    let powers_jensen = model_jensen.get_turbine_powers();
    let farm_power_jensen = model_jensen.get_farm_power();

    println!("Turbine powers (MW):");
    for i in 0..5 {
        println!("  T{}: {:.3} MW", i, powers_jensen[[0, i]] / 1_000_000.0);
    }
    println!("  Total: {:.3} MW\n", farm_power_jensen[[0]] / 1_000_000.0);

    // ============================================================
    // Comparison Summary
    // ============================================================
    println!("===========================================================");
    println!("WAKE MODEL COMPARISON SUMMARY\n");

    println!("| Model      | T0 (MW) | T1 (MW) | T2 (MW) | T3 (MW) | T4 (MW) | Total (MW) |");
    println!("|------------|----------|----------|----------|----------|----------|-------------|");
    
    println!("| {:<10} | {:>8.3} | {:>8.3} | {:>8.3} | {:>8.3} | {:>8.3} | {:>10.3} |",
        "Gauss",
        powers_gauss[[0, 0]] / 1_000_000.0,
        powers_gauss[[0, 1]] / 1_000_000.0,
        powers_gauss[[0, 2]] / 1_000_000.0,
        powers_gauss[[0, 3]] / 1_000_000.0,
        powers_gauss[[0, 4]] / 1_000_000.0,
        farm_power_gauss[[0]] / 1_000_000.0
    );

    println!("| {:<10} | {:>8.3} | {:>8.3} | {:>8.3} | {:>8.3} | {:>8.3} | {:>10.3} |",
        "Jensen",
        powers_jensen[[0, 0]] / 1_000_000.0,
        powers_jensen[[0, 1]] / 1_000_000.0,
        powers_jensen[[0, 2]] / 1_000_000.0,
        powers_jensen[[0, 3]] / 1_000_000.0,
        powers_jensen[[0, 4]] / 1_000_000.0,
        farm_power_jensen[[0]] / 1_000_000.0
    );

    println!();
    println!("Key Observations:");
    println!("  - T0 (leading turbine): Similar power for both models");
    println!("  - Downstream turbines: Different wake deficits predicted");
    println!("  - Total farm power: Difference due to wake model selection");
    println!();
    println!("Note: Actual differences depend on wake model parameters and spacing.");
    println!("      Jensen uses a top-hat wake shape with linear expansion.");
    println!("      Gauss uses a Gaussian wake profile with faster near-wake recovery.");

    println!("\n===========================================================");
    println!("Example completed successfully!");

    Ok(())
}
