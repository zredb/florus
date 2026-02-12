/// Example 5: Yaw Angle Optimization and Wake Steering
///
/// Yaw angle control is a key strategy for maximizing wind farm power output
/// by deflecting wakes away from downstream turbines. This example demonstrates:
///
/// 1. Setting yaw angles on turbines
/// 2. Analyzing yaw angle effects on power
/// 3. Wake steering concepts and benefits
/// 4. Basic yaw optimization approach
///
/// This is the Rust equivalent of Python's yaw optimization examples

use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 5: Yaw Angle Optimization and Wake Steering");
    println!("============================================================\n");

    // Create a 3-turbine wind farm
    let d = 126.0; // NREL 5MW rotor diameter
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    println!("Creating 3-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, layout_y[i]);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    // ============================================================
    // Baseline: No Yaw Misalignment
    // ============================================================
    println!("\n--- Baseline Configuration (0° yaw) ---\n");

    let wind_speeds = Array1::from_vec(vec![9.0]);
    let wind_directions = Array1::from_vec(vec![270.0]); // From West, aligned with turbines
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    let flow_field = florus::core::FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
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
    let baseline_farm_power: f64 = baseline_powers.iter().sum();

    println!("Baseline Results (0° yaw for all turbines):");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, baseline_powers[[0, ti]] / 1000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", baseline_farm_power / 1_000_000.0);

    // ============================================================
    // Yaw Angle Effects on Power
    // ============================================================
    println!("--- Yaw Angle Effects on Turbine Power ---\n");

    // Test different yaw angles
    let yaw_angles_to_test: Vec<f64> = vec![-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0];

    println!("Testing yaw angles from -25° to +25°...");
    println!("  {:>10}  {:>10}  {:>10}", "Yaw (°)", "T0 (kW)", "Farm (MW)");
    println!("  {}", "-".repeat(35));

    let mut best_yaw = 0.0;
    let mut best_power = 0.0;

    for &yaw in &yaw_angles_to_test {
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

        // Set yaw angle (same for all turbines)
        let yaw_array = Array2::from_elem((1, model.farm.n_turbines()), yaw);
        model.set_yaw_angles(yaw_array)?;

        model.run()?;

        let powers = model.get_turbine_powers();
        let farm_power: f64 = powers.iter().sum();

        println!("  {:>10.0}  {:>10.1}  {:>10.3}", yaw, powers[[0, 0]] / 1000.0, farm_power / 1_000_000.0);

        if farm_power > best_power {
            best_power = farm_power;
            best_yaw = yaw;
        }
    }

    println!("\n  Best yaw angle: {}° with {:.3} MW", best_yaw, best_power / 1_000_000.0);

    // ============================================================
    // Individual Yaw Angle Optimization
    // ============================================================
    println!("\n--- Individual Turbine Yaw Optimization ---\n");

    // Test yaw angle on individual turbines
    println!("Testing individual turbine yaw angles (0°, 10°, 20°):");
    println!("  {:>20}  {:>12}", "Configuration", "Farm (MW)");
    println!("  {}", "-".repeat(35));

    let yaw_configs = vec![
        ("All 0°", vec![0.0, 0.0, 0.0]),
        ("T0 10°, others 0°", vec![10.0, 0.0, 0.0]),
        ("T0 20°, others 0°", vec![20.0, 0.0, 0.0]),
        ("T0 10°, T1 5°", vec![10.0, 5.0, 0.0]),
        ("T0 -10°, others 0°", vec![-10.0, 0.0, 0.0]),
    ];

    for (name, yaws) in &yaw_configs {
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

        let yaw_array = Array2::from_elem((1, model.farm.n_turbines()), 0.0);
        let mut yaw_array = yaw_array;
        for (ti, &yaw) in yaws.iter().enumerate() {
            yaw_array[[0, ti]] = yaw;
        }
        model.set_yaw_angles(yaw_array)?;

        model.run()?;

        let powers = model.get_turbine_powers();
        let farm_power: f64 = powers.iter().sum();

        println!("  {:>20}  {:>12.3}", name, farm_power / 1_000_000.0);
    }

    // ============================================================
    // Wake Steering Analysis
    // ============================================================
    println!("\n--- Wake Steering Analysis ---\n");

    println!("Wake steering concept:");
    println!("  When an upstream turbine yaws, it deflects its wake away from");
    println!("  downstream turbines. This can reduce wake losses and increase");
    println!("  total farm power.\n");

    println!("Key parameters for wake steering:");
    println!("  - Yaw angle: Larger angles give more deflection");
    println!("  - Cosine loss: Power loss due to reduced rotor thrust");
    println!("  - Wake recovery: Deflected wakes recover faster");

    // Calculate theoretical wake deflection
    println!("\nTheoretical wake deflection calculation:");
    println!("  For a typical turbine with Ct=0.8, yaw angle γ:");
    println!("  - Wake deflection angle ≈ γ × kd (where kd ≈ 0.1-0.3)");
    println!("  - At 5D downstream distance with 20° yaw:");
    println!("  - Deflection ≈ 20° × 0.1 × 5 = 10 m\n");

    // ============================================================
    // Optimization Strategy Discussion
    // ============================================================
    println!("--- Yaw Optimization Strategies ---\n");

    println!("Common yaw optimization approaches:");
    println!("  1. Sequential Quadratic Programming (SQP)");
    println!("  2. Gradient-based optimization");
    println!("  3. Golden section search");
    println!("  4. Genetic algorithms\n");

    println!("Key considerations:");
    println!("  - Computational cost vs accuracy");
    println!("  - Real-time constraints");
    println!("  - Turbine fatigue loads");
    println!("  - Wind direction uncertainty\n");

    // ============================================================
    // Summary
    // ============================================================
    println!("--- Summary ---\n");

    println!("Yaw Angle Optimization Key Points:");
    println!("  ✓ Yaw angles control wake deflection direction");
    println!("  ✓ Positive yaw (nose right) deflects wake right");
    println!("  ✓ Negative yaw (nose left) deflects wake left");
    println!("  ✓ Cosine loss must be balanced against wake reduction");
    println!("  ✓ Optimal yaw depends on wind direction and layout\n");

    println!("Wake Steering Benefits:");
    println!("  - Reduces wake losses on downstream turbines");
    println!("  - Increases total farm energy production");
    println!("  - Particularly effective in offshore wind farms");
    println!("  - Can increase annual energy production by 1-5%\n");

    println!("============================================================");
    println!("Example completed successfully!");

    Ok(())
}
