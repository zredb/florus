/// Example 10: Different Turbine Types
///
/// FLORIS supports multiple turbine types with different power curves,
/// thrust coefficients, and performance characteristics. This example
/// demonstrates:
///
/// 1. Available turbine types in the library
/// 2. Power curve comparison
/// 3. Mixing turbine types in a farm
/// 4. Performance characteristics analysis
///
/// Turbine types available:
/// - nrel_5MW: NREL 5MW reference turbine
/// - iea_10MW: IEA 10MW offshore turbine
/// - iea_15MW: IEA 15MW offshore turbine

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 10: Different Turbine Types");
    println!("=========================================\n");

    let spacing = 7.0 * 126.0; // 7D spacing

    // ============================================================
    // Turbine Library Overview
    // ============================================================
    println!("--- Available Turbine Types ---\n");

    println!("Turbine Library:");
    println!("  1. nrel_5MW: NREL 5MW reference turbine");
    println!("     - Rated power: 5 MW");
    println!("     - Rotor diameter: 126 m");
    println!("     - Hub height: 90 m");
    println!("     - Designed for onshore applications");
    println!();

    println!("  2. iea_10MW: IEA 10MW offshore turbine");
    println!("     - Rated power: 10 MW");
    println!("     - Rotor diameter: 178 m");
    println!("     - Hub height: 119 m");
    println!("     - Designed for offshore applications");
    println!();

    println!("  3. iea_15MW: IEA 15MW offshore turbine");
    println!("     - Rated power: 15 MW");
    println!("     - Rotor diameter: 242 m");
    println!("     - Hub height: 150 m");
    println!("     - Next-generation offshore turbine");
    println!();

    // ============================================================
    // Single Turbine Power Curves
    // ============================================================
    println!("--- Single Turbine Power Curves ---\n");

    let turbine_types = vec![
        "nrel_5MW".to_string(),
        "iea_10MW".to_string(),
        "iea_15MW".to_string(),
    ];

    let wind_speeds: Vec<f64> = (3..26).map(|i| i as f64).collect();

    println!("Comparing power curves at different wind speeds:");
    println!("  {:>8}  {:>10}  {:>10}  {:>10}", "WS (m/s)", "5MW (kW)", "10MW (kW)", "15MW (kW)");
    println!("  {}", "-".repeat(45));

    for &ws in wind_speeds.iter().step_by(3) {
        let farm_5mw = Farm::new(
            Array1::from_vec(vec![0.0]),
            Array1::from_vec(vec![0.0]),
            vec!["nrel_5MW".to_string()],
        )?;
        let flow_field_5mw = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.14, 1.225,
            Array1::from_vec(vec![0.06]),
            90.0,
        )?;
        let mut model_5mw = florus::FlorisModel {
            farm: farm_5mw,
            flow_field: flow_field_5mw,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };
        model_5mw.initialize_grid()?;
        model_5mw.initialize_flow_field()?;
        model_5mw.run()?;
        let powers_5mw = model_5mw.get_turbine_powers();

        let farm_10mw = Farm::new(
            Array1::from_vec(vec![0.0]),
            Array1::from_vec(vec![0.0]),
            vec!["iea_10MW".to_string()],
        )?;
        let flow_field_10mw = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.14, 1.225,
            Array1::from_vec(vec![0.06]),
            119.0, // IEA 10MW hub height
        )?;
        let mut model_10mw = florus::FlorisModel {
            farm: farm_10mw,
            flow_field: flow_field_10mw,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };
        model_10mw.initialize_grid()?;
        model_10mw.initialize_flow_field()?;
        model_10mw.run()?;
        let powers_10mw = model_10mw.get_turbine_powers();

        let farm_15mw = Farm::new(
            Array1::from_vec(vec![0.0]),
            Array1::from_vec(vec![0.0]),
            vec!["iea_15MW".to_string()],
        )?;
        let flow_field_15mw = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.14, 1.225,
            Array1::from_vec(vec![0.06]),
            150.0, // IEA 15MW hub height
        )?;
        let mut model_15mw = florus::FlorisModel {
            farm: farm_15mw,
            flow_field: flow_field_15mw,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };
        model_15mw.initialize_grid()?;
        model_15mw.initialize_flow_field()?;
        model_15mw.run()?;
        let powers_15mw = model_15mw.get_turbine_powers();

        println!("  {:>8.1}  {:>10.0}  {:>10.0}  {:>10.0}",
                 ws,
                 powers_5mw[[0, 0]] / 1000.0,
                 powers_10mw[[0, 0]] / 1000.0,
                 powers_15mw[[0, 0]] / 1000.0);
    }

    // ============================================================
    // Mixed Turbine Farm
    // ============================================================
    println!("\n--- Mixed Turbine Farm ---\n");

    // Create a farm with different turbine types
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let mixed_turbines = vec![
        "nrel_5MW".to_string(),
        "iea_10MW".to_string(),
        "iea_15MW".to_string(),
    ];

    println!("Creating mixed turbine farm:");
    println!("  Turbine 0: nrel_5MW (5 MW)");
    println!("  Turbine 1: iea_10MW (10 MW)");
    println!("  Turbine 2: iea_15MW (15 MW)");
    println!("  Spacing: {:.0} m\n", spacing);

    let farm_mixed = Farm::new(layout_x.clone(), layout_y.clone(), mixed_turbines)?;

    let flow_field_mixed = florus::core::FlowField::new(
        Array1::from_vec(vec![12.0]), // Above rated for all
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model_mixed = florus::FlorisModel {
        farm: farm_mixed,
        flow_field: flow_field_mixed,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_mixed.initialize_grid()?;
    model_mixed.initialize_flow_field()?;
    model_mixed.run()?;

    let powers_mixed = model_mixed.get_turbine_powers();
    let farm_power_mixed: f64 = powers_mixed.iter().sum();

    println!("Mixed Farm Results (12 m/s wind):");
    for ti in 0..model_mixed.farm.n_turbines() {
        let turbine_name = match ti {
            0 => "nrel_5MW",
            1 => "iea_10MW",
            _ => "iea_15MW",
        };
        let power = powers_mixed[[0, ti]] / 1000.0;
        println!("  Turbine {} ({}): {:.1} kW", ti, turbine_name, power);
    }
    println!("  Total Farm Power: {:.2} MW\n", farm_power_mixed / 1_000_000.0);

    // ============================================================
    // Same Capacity Comparison
    // ============================================================
    println!("--- Same Capacity Comparison ---\n");

    // Compare 3×5MW vs 1×15MW
    let farm_3x5mw = Farm::new(
        Array1::from_vec(vec![0.0, spacing, 2.0 * spacing]),
        Array1::from_vec(vec![0.0, 0.0, 0.0]),
        vec!["nrel_5MW".to_string(); 3],
    )?;

    let flow_field_3x5mw = florus::core::FlowField::new(
        Array1::from_vec(vec![12.0]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model_3x5mw = florus::FlorisModel {
        farm: farm_3x5mw,
        flow_field: flow_field_3x5mw,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_3x5mw.initialize_grid()?;
    model_3x5mw.initialize_flow_field()?;
    model_3x5mw.run()?;

    let powers_3x5mw = model_3x5mw.get_turbine_powers();
    let farm_power_3x5mw: f64 = powers_3x5mw.iter().sum();

    // Single 15MW turbine
    let farm_1x15mw = Farm::new(
        Array1::from_vec(vec![0.0]),
        Array1::from_vec(vec![0.0]),
        vec!["iea_15MW".to_string()],
    )?;

    let flow_field_1x15mw = florus::core::FlowField::new(
        Array1::from_vec(vec![12.0]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        150.0,
    )?;

    let mut model_1x15mw = florus::FlorisModel {
        farm: farm_1x15mw,
        flow_field: flow_field_1x15mw,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_1x15mw.initialize_grid()?;
    model_1x15mw.initialize_flow_field()?;
    model_1x15mw.run()?;

    let powers_1x15mw = model_1x15mw.get_turbine_powers();
    let farm_power_1x15mw: f64 = powers_1x15mw.iter().sum();

    println!("Comparing same rated capacity (15 MW):");
    println!("  Configuration 1: 3 × nrel_5MW (5 MW each)");
    println!("    Total Power: {:.2} MW", farm_power_3x5mw / 1_000_000.0);
    println!("    Wake losses due to multiple turbines");
    println!();
    println!("  Configuration 2: 1 × iea_15MW (15 MW)");
    println!("    Total Power: {:.2} MW", farm_power_1x15mw / 1_000_000.0);
    println!("    No wake losses (single turbine)");
    println!();

    let wake_loss_percent = (1.0 - farm_power_3x5mw / (3.0 * farm_power_1x15mw)) * 100.0;
    println!("  Wake loss impact: {:.1}% of total capacity", wake_loss_percent);

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Turbine Selection Considerations:");
    println!("  1. Large turbines:");
    println!("     - Higher rated power per foundation");
    println!("     - Lower wakes per MW installed");
    println!("     - Higher hub heights access better wind");
    println!();

    println!("  2. Small turbines:");
    println!("     - More flexible layout options");
    println!("     - Lower transportation constraints");
    println!("     - More turbines = more wake interactions");
    println!();

    println!("  3. Mixing turbines:");
    println!("     - Can optimize for different wind regimes");
    println!("     - More complex farm management");
    println!("     - Consider maintenance logistics");
    println!();

    println!("Key Takeaways:");
    println!("  - Larger turbines generally produce more energy");
    println!("  - Wake losses scale with turbine count");
    println!("  - Hub height affects wind resource quality");
    println!("  - Rotor diameter affects swept area");
    println!("  - Optimal choice depends on site conditions");

    println!("\n=========================================");
    println!("Example completed successfully!");

    Ok(())
}
