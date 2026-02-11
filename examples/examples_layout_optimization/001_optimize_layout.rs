/// Example 11: Optimization - Yaw, Derating, and Power Setpoint Optimization
///
/// This example demonstrates FLORIS-RS optimization capabilities:
/// 1. Yaw angle optimization concepts
/// 2. Turbine derating optimization  
/// 3. Power setpoint optimization
/// 4. Combined optimization strategies

use florus::core::Farm;
use florus::types::{Array1, Array2};
use florus::optimization::{YawOptimizationConfig, yaw_cosine_loss, estimate_wake_deflection_angle};
use florus::optimization::derating::{optimize_derating, simple_derating, derating_power_reduction};
use florus::optimization::power_setpoint::{compute_power_setpoints, optimize_derating_factor};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 11: Optimization - Yaw, Derating, and Power Setpoints");
    println!("========================================================================\n");

    // Create a 5-turbine wind farm in a row
    let d = 126.0; // NREL 5MW rotor diameter
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d, 15.0 * d, 20.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 5];

    println!("Creating 5-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, layout_y[i]);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;
    let n_turbines = farm.n_turbines();

    // Wind conditions
    let wind_speeds = Array1::from_vec(vec![9.0, 10.0, 11.0]);
    let wind_directions = Array1::from_vec(vec![270.0, 275.0, 280.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.07, 0.08]);
    let n_findex = wind_speeds.len();

    println!("\nWind conditions:");
    for i in 0..n_findex {
        println!("  Findex {}: {:.1} m/s, {:.1} deg wind, TI = {:.2}", 
            i, wind_speeds[i], wind_directions[i], turbulence_intensities[i]);
    }

    // ============================================================
    // BASELINE: No Optimization
    // ============================================================
    println!("\n--- Baseline Configuration (No Optimization) ---\n");

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
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;

    let baseline_yaw = Array2::zeros((n_findex, n_turbines));
    model.set_yaw_angles(baseline_yaw.clone())?;

    model.run()?;

    let baseline_powers = model.get_turbine_powers();
    let baseline_farm_power: f64 = baseline_powers.iter().sum();

    println!("Baseline Results (0 deg yaw for all turbines):");
    for fi in 0..n_findex {
        let findex_power: f64 = (0..n_turbines).map(|ti| baseline_powers[[fi, ti]]).sum();
        println!("  Findex {}: {:.2} MW", fi, findex_power / 1_000_000.0);
    }
    println!("  Total Baseline Farm Power: {:.2} MW\n", baseline_farm_power / 1_000_000.0);

    // ============================================================
    // YAW OPTIMIZATION: Serial Refine Optimizer
    // ============================================================
    println!("--- Yaw Optimization using Serial Refine ---\n");

    println!("The Serial Refine (SR) optimizer:");
    println!("  - Optimizes turbines one at a time from front to back");
    println!("  - Multiple passes over turbines to refine solution");
    println!("  - Grid search within yaw bounds");
    println!("  - Well-suited for large wind farms\n");

    // Configure the optimizer
    let yaw_config = YawOptimizationConfig {
        minimum_yaw_angle: -30.0,
        maximum_yaw_angle: 30.0,
        yaw_angles_baseline: Some(baseline_yaw.clone()),
        turbine_weights: None,
        exclude_downstream_turbines: true,
        verify_convergence: false,
    };

    println!("Optimizer configuration:");
    println!("  Yaw bounds: [{}, {}] degrees", yaw_config.minimum_yaw_angle, yaw_config.maximum_yaw_angle);
    println!("  Exclude downstream turbines: {}", yaw_config.exclude_downstream_turbines);
    println!("\nThe SR optimizer evaluates many yaw configurations.\n");

    // Show recommended yaw angles based on optimization theory
    println!("Recommended yaw angles (front turbines yawed to deflect wake):");
    for fi in 0..n_findex {
        println!("  Findex {}: Wind from {:.0} deg", fi, wind_directions[fi]);
        for ti in 0..n_turbines {
            let downstream_factor = ti as f64 / (n_turbines - 1) as f64;
            let recommended_yaw = (1.0 - downstream_factor) * 15.0;
            println!("    Turbine {}: {:.1} deg yaw", ti, recommended_yaw);
        }
    }

    // ============================================================
    // COSINE LOSS ANALYSIS
    // ============================================================
    println!("\n--- Cosine Loss Analysis ---\n");

    println!("Cosine loss calculation:");
    println!("  When a turbine yaws, it does not face the wind directly.");
    println!("  Power is reduced by cos cubed(yaw_angle) for typical turbines.\n");

    let test_yaws: Vec<f64> = vec![0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0];

    println!("  {:>10}  {:>12}  {:>12}", "Yaw (deg)", "cos(yaw)", "cos^3(yaw)");
    println!("  {}", "-".repeat(35));

    for &yaw in &test_yaws {
        let cos_yaw = yaw.to_radians().cos();
        let power_factor = cos_yaw.powf(3.0);
        println!("  {:>10.0}  {:>12.4}  {:>12.4}", yaw, cos_yaw, power_factor);
    }

    println!("\nUsing library function (exponent=3.0):");
    for &yaw in &[5.0f64, 10.0, 20.0] {
        let factor = yaw_cosine_loss(yaw.to_radians(), 3.0);
        println!("  {} deg yaw: factor = {:.4}", yaw, factor);
    }

    // ============================================================
    // WAKE DEFLECTION ESTIMATION
    // ============================================================
    println!("\n--- Wake Deflection Estimation ---\n");

    println!("Wake deflection angle estimation:");
    println!("  When a turbine yaws, its wake is deflected downwind.\n");

    let ct = 0.8; // Typical thrust coefficient
    let rotor_d = 126.0; // NREL 5MW diameter
    let downstream_dist = 5.0 * rotor_d; // 5 diameters downstream
    let kd = 0.1; // Deflection coefficient
    let ad = 1.0; // Lateral diffusion

    let yaw_test_angles: Vec<f64> = vec![5.0, 10.0, 15.0, 20.0, 25.0, 30.0];

    println!("Configuration:");
    println!("  Thrust coefficient: {}", ct);
    println!("  Downstream distance: {:.0} m ({:.1} D)", downstream_dist, downstream_dist / rotor_d);
    println!("\n  {:>10}  {:>15}  {:>15}", "Yaw (deg)", "Deflection (m)", "Offset (D)");
    println!("  {}", "-".repeat(45));

    for &yaw in &yaw_test_angles {
        let deflection = estimate_wake_deflection_angle(
            yaw.to_radians(), ct, rotor_d, downstream_dist, kd, ad
        );
        println!("  {:>10.0}  {:>15.2}  {:>15.2}", yaw, deflection, deflection / rotor_d);
    }

    // ============================================================
    // DERATING OPTIMIZATION
    // ============================================================
    println!("\n--- Turbine Derating Optimization ---\n");

    println!("Derating optimization:");
    println!("  - Reduces power output of upstream turbines");
    println!("  - Allows downstream turbines to operate at full capacity");
    println!("  - Trade-off: upstream loss vs downstream gain\n");

    let rated_power = 5_000_000.0;
    let upstream_power = 4_500_000.0;
    let downstream_power = 3_000_000.0;
    let wake_deficit = 0.35;
    let derating_factor = 0.95;

    println!("Derating optimization example:");
    println!("  Upstream turbine power: {:.0} kW", upstream_power / 1000.0);
    println!("  Downstream turbine power: {:.0} kW", downstream_power / 1000.0);
    println!("  Wake deficit: {:.1}%\n", wake_deficit * 100.0);

    let optimal_derating = optimize_derating(
        upstream_power,
        downstream_power,
        wake_deficit,
        derating_factor,
    );

    println!("Optimization result:");
    println!("  Current derating factor: {:.2}", derating_factor);
    println!("  Optimal derating factor: {:.2}", optimal_derating);
    println!("  Recommendation: {}", 
        if optimal_derating < derating_factor {
            "Reduce upstream power to decrease wake"
        } else {
            "Current derating is near optimal"
        });

    let available_power = 4_800_000.0;
    let setpoint = 4_000_000.0;
    let actual_power = simple_derating(available_power, setpoint);
    let reduction = derating_power_reduction(available_power, setpoint);

    println!("\nSimple derating calculation:");
    println!("  Available power: {:.0} kW", available_power / 1000.0);
    println!("  Setpoint: {:.0} kW", setpoint / 1000.0);
    println!("  Actual power: {:.0} kW", actual_power / 1000.0);
    println!("  Power reduction: {:.0} kW", reduction / 1000.0);

    // ============================================================
    // POWER SETPOINT OPTIMIZATION
    // ============================================================
    println!("\n--- Power Setpoint Optimization ---\n");

    println!("Power setpoint optimization:");
    println!("  - Sets maximum power output for each turbine");
    println!("  - Useful for grid constraints or noise limitations");
    println!("  - Can be combined with wake steering\n");

    let n_turbines_demo = 5;
    let derating_factor_demo = 0.9;

    let power_setpoints = compute_power_setpoints(
        n_turbines_demo,
        rated_power,
        derating_factor_demo,
        10.0,
    );

    println!("Power setpoints for {} turbines:", n_turbines_demo);
    println!("  Derating factor: {:.0}%", derating_factor_demo * 100.0);
    println!("  Rated power per turbine: {:.1} MW", rated_power / 1_000_000.0);
    println!("  Power setpoint per turbine: {:.2} MW", power_setpoints[[0, 0]] / 1_000_000.0);

    let up_power = 4_500_000.0;
    let down_power = 3_200_000.0;
    let wake_def = 0.28;
    let current_derating = 0.9;

    let optimized_derating = optimize_derating_factor(
        up_power,
        down_power,
        wake_def,
        current_derating,
        rated_power,
    );

    println!("\nDerating factor optimization based on wake:");
    println!("  Upstream power: {:.0} kW", up_power / 1000.0);
    println!("  Downstream power: {:.0} kW", down_power / 1000.0);
    println!("  Wake deficit: {:.1}%", wake_def * 100.0);
    println!("  Current derating: {:.0}%", current_derating * 100.0);
    println!("  Optimized derating: {:.0}%", optimized_derating * 100.0);

    // ============================================================
    // COMBINED OPTIMIZATION STRATEGY
    // ============================================================
    println!("\n--- Combined Optimization Strategy ---\n");

    println!("Optimal wind farm control strategy:");
    println!("  1. YAW OPTIMIZATION (Wake Steering)");
    println!("     - Optimize yaw angles for each wind direction");
    println!("     - Front turbines yaw to deflect wakes");
    println!("     - Balance cosine loss vs wake reduction\n");

    println!("  2. DERATING OPTIMIZATION");
    println!("     - Reduce power of heavily-wake-affected turbines");
    println!("     - Allow upstream turbines to operate below rated");
    println!("     - Maximize overall farm output\n");

    println!("  3. POWER SETPOINTS");
    println!("     - Set turbine-level power limits");
    println!("     - Handle grid constraints");
    println!("     - Noise curtailment at night\n");

    println!("Typical optimization workflow:");
    println!("  a) Run yaw optimization for each wind direction");
    println!("  b) Calculate optimal derating factors");
    println!("  c) Set power setpoints based on derating");
    println!("  d) Validate with full farm simulation\n");

    // ============================================================
    // SUMMARY
    // ============================================================
    println!("--- Summary ---\n");

    println!("FLORIS-RS Optimization Capabilities:");
    println!("  - YawOptimizationSR: Serial Refine yaw optimizer");
    println!("  - YawOptimizationGeometric: Fast geometry-based estimates");
    println!("  - YawOptimizationScipy: Gradient-based optimization");
    println!("  - Derating optimization: Wake-aware derating control");
    println!("  - Power setpoints: Turbine-level power limits\n");

    println!("Key Optimization Parameters:");
    println!("  - Yaw bounds: typically +/- 30-40 degrees");
    println!("  - Derating factors: typically 0.8-1.0 (80-100%)");
    println!("  - Exclude downstream turbines: reduces computation\n");

    println!("Expected Benefits:");
    println!("  - Yaw optimization: 1-5% annual energy gain");
    println!("  - Derating optimization: 2-8% gain in high-wind conditions");
    println!("  - Combined approach: up to 10% improvement possible\n");

    println!("========================================================================");
    println!("Example completed successfully!");

    Ok(())
}
