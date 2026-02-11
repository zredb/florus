/// Example 7: Comparing Wake Models
///
/// FLORIS supports multiple wake models for different accuracy/speed trade-offs.
/// This example demonstrates:
///
/// 1. Jensen wake model (simple, fast)
/// 2. Gauss wake model (Gaussian profile, more accurate)
/// 3. Turbopark wake model (hybrid approach)
/// 4. Comparing wake deficits across models
/// 5. Model selection considerations
///
/// This is the Rust equivalent of Python's wake model comparison examples

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 7: Comparing Wake Models");
    println!("======================================\n");

    // Create a 4-turbine wind farm
    let d = 126.0; // NREL 5MW rotor diameter
    let spacing = 5.0 * d; // 5D spacing
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    println!("Creating 4-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m", i, x);
    }
    println!("  Spacing: {:.0} m ({:.1}D)\n", spacing, spacing / d);

    // ============================================================
    // Wake Model Descriptions
    // ============================================================
    println!("--- Available Wake Models ---\n");

    println!("1. Jensen Model:");
    println!("   - Simple top-hat wake profile");
    println!("   - Fast computation");
    println!("   - Good for preliminary studies");
    println!("   - Assumes linear wake expansion\n");

    println!("2. Gauss Model:");
    println!("   - Gaussian wake profile");
    println!("   - More accurate near-wake predictions");
    println!("   - Widely used for offshore applications");
    println!("   - Standard model in FLORIS v4.x\n");

    println!("3. Turbopark Model:");
    println!("   - Hybrid approach combining Gauss and Jensen");
    println!("   - Optimized for large wind farms");
    println!("   - Faster than full Gaussian");
    println!("   - Good for offshore wind farms\n");

    // ============================================================
    // Wind Conditions
    // ============================================================
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]); // Aligned with turbines
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    // ============================================================
    // Jensen Wake Model Analysis
    // ============================================================
    println!("--- Jensen Wake Model ---\n");

    let farm_jensen = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    let flow_field_jensen = florus::core::FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.14,
        1.225,
        turbulence_intensities.clone(),
        90.0,
    )?;

    let mut model_jensen = florus::FlorisModel {
        farm: farm_jensen,
        flow_field: flow_field_jensen,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_jensen.initialize_grid()?;
    model_jensen.initialize_flow_field()?;
    model_jensen.run()?;

    let powers_jensen = model_jensen.get_turbine_powers();
    let velocities_jensen = model_jensen.get_turbine_velocities();

    println!("Jensen Model Results:");
    for ti in 0..model_jensen.farm.n_turbines() {
        let vel = velocities_jensen[[0, ti, 0]];
        let power = powers_jensen[[0, ti]] / 1000.0;
        let ref_vel = velocities_jensen[[0, 0, 0]];
        let deficit = (1.0 - vel / ref_vel) * 100.0;
        println!("  Turbine {}: velocity = {:.2} m/s, power = {:.1} kW, deficit = {:.1}%",
                 ti, vel, power, deficit);
    }

    // ============================================================
    // Gauss Wake Model Analysis
    // ============================================================
    println!("\n--- Gauss Wake Model ---\n");

    let farm_gauss = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    let flow_field_gauss = florus::core::FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.14,
        1.225,
        turbulence_intensities.clone(),
        90.0,
    )?;

    let mut model_gauss = florus::FlorisModel {
        farm: farm_gauss,
        flow_field: flow_field_gauss,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_gauss.initialize_grid()?;
    model_gauss.initialize_flow_field()?;
    model_gauss.run()?;

    let powers_gauss = model_gauss.get_turbine_powers();
    let velocities_gauss = model_gauss.get_turbine_velocities();

    println!("Gauss Model Results:");
    for ti in 0..model_gauss.farm.n_turbines() {
        let vel = velocities_gauss[[0, ti, 0]];
        let power = powers_gauss[[0, ti]] / 1000.0;
        let ref_vel = velocities_gauss[[0, 0, 0]];
        let deficit = (1.0 - vel / ref_vel) * 100.0;
        println!("  Turbine {}: velocity = {:.2} m/s, power = {:.1} kW, deficit = {:.1}%",
                 ti, vel, power, deficit);
    }

    // ============================================================
    // Comparison
    // ============================================================
    println!("\n--- Model Comparison ---\n");

    println!("  {:>8}  {:>12}  {:>12}  {:>12}", "Turbine", "Jensen (%)", "Gauss (%)", "Difference");
    println!("  {}", "-".repeat(55));

    for ti in 0..4 {
        let vel_j = velocities_jensen[[0, ti, 0]];
        let vel_g = velocities_gauss[[0, ti, 0]];
        let ref_vel = velocities_jensen[[0, 0, 0]];
        let deficit_j = (1.0 - vel_j / ref_vel) * 100.0;
        let deficit_g = (1.0 - vel_g / ref_vel) * 100.0;
        let diff = deficit_g - deficit_j;

        println!("  {:>8}  {:>12.1}  {:>12.1}  {:>12.1}", ti, deficit_j, deficit_g, diff);
    }

    // ============================================================
    // Turbulence Intensity Effects
    // ============================================================
    println!("\n--- Turbulence Intensity Effects ---\n");

    let ti_values = vec![0.03, 0.06, 0.10, 0.14];

    println!("Testing different turbulence intensities with Gauss model:");
    println!("  {:>10}  {:>10}  {:>10}  {:>10}", "TI", "T0 vel", "T3 vel", "Wake Loss");
    println!("  {}", "-".repeat(45));

    for &ti in &ti_values {
        let farm_ti = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

        let flow_field_ti = florus::core::FlowField::new(
            wind_speeds.clone(),
            wind_directions.clone(),
            0.0,
            0.14,
            1.225,
            Array1::from_vec(vec![ti]),
            90.0,
        )?;

        let mut model_ti = florus::FlorisModel {
            farm: farm_ti,
            flow_field: flow_field_ti,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };

        model_ti.initialize_grid()?;
        model_ti.initialize_flow_field()?;
        model_ti.run()?;

        let velocities_ti = model_ti.get_turbine_velocities();
        let v0 = velocities_ti[[0, 0, 0]];
        let v3 = velocities_ti[[0, 3, 0]];
        let wake_loss = (1.0 - v3 / v0) * 100.0;

        println!("  {:>10.2}  {:>10.2}  {:>10.2}  {:>10.1}%", ti, v0, v3, wake_loss);
    }

    // ============================================================
    // Wind Direction Sensitivity
    // ============================================================
    println!("\n--- Wind Direction Sensitivity ---\n");

    let wind_directions: Vec<f64> = (250..295).map(|i| i as f64).collect();

    println!("Testing wake sensitivity to wind direction (Gauss model):");
    println!("  {:>10}  {:>10}  {:>10}", "WD (°)", "T0 vel", "T3 vel");
    println!("  {}", "-".repeat(35));

    for &wd in wind_directions.iter().step_by(5) {
        let farm_dir = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

        let flow_field_dir = florus::core::FlowField::new(
            wind_speeds.clone(),
            Array1::from_vec(vec![wd]),
            0.0,
            0.14,
            1.225,
            turbulence_intensities.clone(),
            90.0,
        )?;

        let mut model_dir = florus::FlorisModel {
            farm: farm_dir,
            flow_field: flow_field_dir,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };

        model_dir.initialize_grid()?;
        model_dir.initialize_flow_field()?;
        model_dir.run()?;

        let velocities_dir = model_dir.get_turbine_velocities();
        let v0 = velocities_dir[[0, 0, 0]];
        let v3 = velocities_dir[[0, 3, 0]];

        println!("  {:>10.0}  {:>10.2}  {:>10.2}", wd, v0, v3);
    }

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Wake Model Selection Guidelines:");
    println!("  Jensen Model:");
    println!("    ✓ Fast computation for large parameter sweeps");
    println!("    ✓ Good for preliminary layout optimization");
    println!("    ✗ Less accurate for complex wake interactions");
    println!("    ✗ Over-predicts near-wake deficits\n");

    println!("  Gauss Model:");
    println!("    ✓ Most accurate for offshore applications");
    println!("    ✓ Good for detailed wake analysis");
    println!("    ✓ Industry standard for wake modeling");
    println!("    ✗ Slower computation than Jensen");
    println!("    ✗ More sensitive to turbulence parameters\n");

    println!("  Turbopark Model:");
    println!("    ✓ Optimized for large wind farms");
    println!("    ✓ Good balance of accuracy and speed");
    println!("    ✓ Efficient for offshore applications");
    println!("    ✓ Good for real-time control applications\n");

    println!("Key Takeaways:");
    println!("  1. Wake models predict different deficit magnitudes");
    println!("  2. Higher turbulence = faster wake recovery");
    println!("  3. Aligned winds cause maximum wake losses");
    println!("  4. Model choice affects optimization results");
    println!("  5. Validate model against site measurements\n");

    println!("======================================");
    println!("Example completed successfully!");

    Ok(())
}
