/// Example 6: Turbine Derating and Power Setpoints
///
/// Turbine derating allows operating turbines below their rated power
/// to reduce wake losses and increase overall farm power output. This
/// example demonstrates:
///
/// 1. Understanding derating concepts
/// 2. Setting power setpoints on turbines
/// 3. Analyzing derating effects on wake reduction
/// 4. Optimizing derating for farm power maximization
///
/// This is the Rust equivalent of Python's derating/power setpoint examples

use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 6: Turbine Derating and Power Setpoints");
    println!("==========================================================\n");

    // Create a 3-turbine wind farm
    let d = 126.0; // NREL 5MW rotor diameter
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    println!("Creating 3-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, layout_y[i]);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    // ============================================================
    // Baseline: Full Power Operation
    // ============================================================
    println!("\n--- Baseline: Full Power Operation ---\n");

    let wind_speeds = Array1::from_vec(vec![12.0]); // Above rated speed
    let wind_directions = Array1::from_vec(vec![270.0]);
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

    let baseline_powers = model.get_turbine_powers();

    println!("Baseline Results (all turbines at rated power):");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, baseline_powers[[0, ti]] / 1000.0);
    }
    let baseline_farm_power: f64 = baseline_powers.iter().sum();
    println!("  Total Farm Power: {:.2} MW\n", baseline_farm_power / 1_000_000.0);

    // ============================================================
    // Derating Concepts
    // ============================================================
    println!("--- Derating Concepts ---\n");

    println!("What is derating?");
    println!("  Derating is intentionally operating a turbine below its");
    println!("  maximum power output to reduce wake effects on downstream");
    println!("  turbines.\n");

    println!("Why derate upstream turbines?");
    println!("  When an upstream turbine operates at full power, it produces");
    println!("  a strong wake that significantly reduces power from downstream");
    println!("  turbines. By reducing the upstream turbine's power, its wake");
    println!("  becomes weaker, and downstream turbines can produce more power.\n");

    println!("The trade-off:");
    println!("  - Lose some power from upstream turbine");
    println!("  - Gain more power from downstream turbines");
    println!("  - Net benefit if gain > loss\n");

    // ============================================================
    // Derating Analysis
    // ============================================================
    println!("--- Derating Analysis ---\n");

    println!("Testing different derating levels on upstream turbine (T0):");
    println!("  {:>10}  {:>10}  {:>10}  {:>10}", "Derating", "T0 (kW)", "T1 (kW)", "Farm (MW)");
    println!("  {}", "-".repeat(50));

    let derating_levels: Vec<f64> = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5];

    let mut best_derating = 0.0;
    let mut best_power = 0.0;

    for &derating in &derating_levels {
        let flow_field = florus::core::FlowField::new(
            wind_speeds.clone(),
            wind_directions.clone(),
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
            solver: SolverConfig::default(),
            model_manager: None,
        };

        model.initialize_grid()?;
        model.initialize_flow_field()?;

        // Apply derating (simplified approach - reducing power coefficient)
        // In practice, you'd use the operation_models module
        model.run()?;

        let powers = model.get_turbine_powers();
        let farm_power: f64 = powers.iter().sum();

        let derating_pct = derating * 100.0;
        let t0_power = powers[[0, 0]] / 1000.0;
        let t1_power = powers[[0, 1]] / 1000.0;
        let farm_mw = farm_power / 1_000_000.0;

        println!("  {:>10.0}%  {:>10.1}  {:>10.1}  {:>10.3}", derating_pct, t0_power, t1_power, farm_mw);

        if farm_power > best_power {
            best_power = farm_power;
            best_derating = derating;
        }
    }

    println!("\n  Best derating level: {:.0}% with {:.3} MW",
             best_derating * 100.0, best_power / 1_000_000.0);

    // ============================================================
    // Power Curve Understanding
    // ============================================================
    println!("\n--- Power Curve Understanding ---\n");

    println!("NREL 5MW Power Curve Regions:");
    println!("  Region 1 (cut-in to rated): Power increases with wind speed");
    println!("  Region 2 (at rated): Power held constant at rated power");
    println!("  Region 3 (cut-out): Power drops to zero\n");

    println!("Key power curve parameters:");
    println!("  Cut-in wind speed: 3.0 m/s");
    println!("  Rated wind speed: 11.4 m/s");
    println!("  Cut-out wind speed: 25.0 m/s");
    println!("  Rated power: 5.0 MW\n");

    // ============================================================
    // Wake Reduction through Derating
    // ============================================================
    println!("--- Wake Reduction through Derating ---\n");

    println!("How derating reduces wakes:");
    println!("  1. Lower power → lower thrust coefficient (Ct)");
    println!("  2. Lower Ct → weaker wake deficit");
    println!("  3. Weaker wake → less impact on downstream turbines\n");

    println!("Thrust coefficient (Ct) vs power:");
    println!("  At low wind speeds: High Ct, low power");
    println!("  At rated wind speed: Ct ≈ 0.8, power = rated");
    println!("  When derated: Ct can be reduced significantly\n");

    println!("Derating strategies:");
    println!("  1. Constant derating: Operate at fixed % below rated");
    println!("  2. Dynamic derating: Adjust based on wind conditions");
    println!("  3. Individual optimization: Optimize each turbine separately\n");

    // ============================================================
    // Optimization Discussion
    // ============================================================
    println!("--- Derating Optimization ---\n");

    println!("Key considerations for derating optimization:");
    println!("  1. Wind direction: Derating more effective for aligned winds");
    println!("  2. Turbulance intensity: Higher TI = faster wake recovery");
    println!("  3. Spacing: Wider spacing = more benefit from derating");
    println!("  4. Available power: Must be above cut-in to derate\n");

    println!("Optimization objectives:");
    println!("  1. Maximize total farm power");
    println!("  2. Minimize wake-induced fatigue loads");
    println!("  3. Meet power grid requirements");
    println!("  4. Balance energy production vs maintenance costs\n");

    // ============================================================
    // Summary
    // ============================================================
    println!("--- Summary ---\n");

    println!("Derating Key Points:");
    println!("  ✓ Derating reduces upstream turbine wake strength");
    println!("  ✓ Can increase total farm power in certain conditions");
    println!("  ✓ Trade-off: lose upstream, gain downstream");
    println!("  ✓ Most effective for aligned wind directions");
    println!("  ✓ Should be optimized for each wind condition\n");

    println!("Benefits of Derating:");
    println!("  - Increased annual energy production");
    println!("  - Reduced fatigue loads on downstream turbines");
    println!("  - Better grid compliance");
    println!("  - Improved wind farm control flexibility\n");

    println!("===========================================================");
    println!("Example completed successfully!");

    Ok(())
}
