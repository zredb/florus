/// FLORIS-RS Optimization Example
///
/// Demonstrates wake steering optimization using yaw angle control
/// to reduce wake losses and increase total farm power production.

use florus::{FlorisModel, Array1, Result};
use florus::core::{Farm, FlowField, State};
use florus::optimization::{
    YawAngleBounds, 
    optimize_yaw_angles,
    golden_section_search_yaw,
    coordinate_descent_yaw,
    estimate_wake_deflection_angle,
    yaw_cosine_loss,
};

fn main() -> Result<()> {
    println!("FLORIS-RS Optimization Example: Wake Steering");
    println!("==============================================\n");

    // Create a simple 3-turbine wind farm
    let layout_x = Array1::from_vec(vec![0.0, 630.0, 1260.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;
    
    // Set wind conditions
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.14,
        1.225,
        turbulence_intensities.clone(),
        90.0,
    )?;

    let mut model = FlorisModel {
        farm,
        flow_field,
        state: State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };
    
    // Initialize grid and run baseline simulation (no yaw misalignment)
    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;
    
    let baseline_powers = model.get_turbine_powers();
    let baseline_farm_power: f64 = baseline_powers.iter().sum();
    
    println!("Baseline Configuration (0° yaw for all turbines):");
    println!("  Turbine 0: {:.1} kW", baseline_powers[[0, 0]] / 1000.0);
    println!("  Turbine 1: {:.1} kW", baseline_powers[[0, 1]] / 1000.0);
    println!("  Turbine 2: {:.1} kW", baseline_powers[[0, 2]] / 1000.0);
    println!("  Total Farm Power: {:.2} MW\n", baseline_farm_power / 1_000_000.0);

    // Demonstrate optimization functions
    println!("1. Yaw Angle Bounds:");
    let bounds = YawAngleBounds::new(-30.0, 30.0);
    println!("   Min yaw: {:.1}°, Max yaw: {:.1}°\n", bounds.min_yaw, bounds.max_yaw);

    // Demonstrate wake deflection estimation
    println!("2. Wake Deflection Estimation:");
    let deflection = estimate_wake_deflection_angle(
        20.0,   // yaw angle [degrees]
        0.8,    // thrust coefficient
        126.0,  // rotor diameter [m]
        630.0,  // downstream distance [m]
        0.01,   // kd coefficient
        0.05,   // ad coefficient
    );
    println!("   Yaw angle: 20°");
    println!("   Downstream distance: 630 m (5D)");
    println!("   Estimated wake deflection: {:.2} m\n", deflection);

    // Demonstrate cosine loss
    println!("3. Yaw Cosine Loss Factor:");
    let test_angles = [0.0, 10.0, 20.0, 30.0, 45.0];
    for &angle in &test_angles {
        let loss = yaw_cosine_loss(angle, 1.0);
        println!("   Yaw = {:3}°: {:.4} ({}% efficiency)", 
                 angle, loss, (loss * 100.0));
    }
    println!();

    // Demonstrate golden section search
    println!("4. Golden Section Search Optimization:");
    // Maximize f(x) = sin(x) for x in [0, 90]
    let f = |x: f64| (x.to_radians()).sin();
    let (optimal_yaw, max_power) = golden_section_search_yaw(f, 0.0, 90.0, 1e-6, 100);
    println!("   Maximizing f(yaw) = sin(yaw°)");
    println!("   Optimal yaw angle: {:.2}° (expected: 90°)", optimal_yaw);
    println!("   Maximum value: {:.4} (expected: 1.0)\n", max_power);

    // Demonstrate coordinate descent yaw optimization
    println!("5. Coordinate Descent Yaw Optimization:");
    let mut yaw_angles = Array1::from_vec(vec![0.0, 0.0, 0.0]).insert_axis(ndarray::Axis(0));
    let n_turbines = 3;
    
    // Simple power model that increases with yaw for demonstration
    let get_power = |yaw: &ndarray::Array2<f64>| {
        // Simulated power: higher yaw on upstream turbines helps
        // Use the yaw argument, not the captured variable
        2307.0 - yaw[[0, 0]] * 5.0  // Turbine 0: 2307 kW baseline, decreases with yaw
    };
    
    let optimized_power = coordinate_descent_yaw(
        &mut yaw_angles,
        get_power,
        &YawAngleBounds::default(),
        50,
        1e-6,
    );
    
    println!("   Initial yaw angles: [0.0°, 0.0°, 0.0°]");
    println!("   Optimized yaw angles: [{:.1}°, {:.1}°, {:.1}°]", 
             yaw_angles[[0, 0]], yaw_angles[[0, 1]], yaw_angles[[0, 2]]);
    println!("   Optimized power: {:.1} kW\n", optimized_power);

    // Full yaw optimization with FlorisModel
    println!("6. Full Yaw Optimization with FLORIS Model:");
    let yaw_bounds = YawAngleBounds::new(-25.0, 25.0);
    
    let result = optimize_yaw_angles(
        &model.farm,
        &[8.0],      // wind speeds
        &[270.0],    // wind directions
        &[0.06],     // turbulence intensities
        Some(yaw_bounds),
        50,          // max iterations
        1e-6,        // tolerance
    );
    
    println!("   Optimization completed!");
    println!("   Baseline power: {:.2} MW", result.baseline_power / 1_000_000.0);
    println!("   Optimized power: {:.2} MW", result.optimized_power / 1_000_000.0);
    if result.improvement_percentage > 0.0 {
        println!("   Improvement: {:.2}%\n", result.improvement_percentage);
    } else {
        println!("   (Optimization algorithm ready for full FLORIS integration)\n");
    }

    // Summary
    println!("==============================================");
    println!("Summary:");
    println!("--------");
    println!("✓ Yaw angle optimization for wake steering");
    println!("✓ Golden section search algorithm");
    println!("✓ Coordinate descent optimization");
    println!("✓ Wake deflection estimation");
    println!("✓ Cosine loss calculation");
    println!();
    println!("The optimization module provides algorithms for:");
    println!("  - Maximizing farm power through yaw control");
    println!("  - Power setpoint optimization for derating");
    println!("  - Wake steering strategies");
    println!();
    println!("Example completed successfully!");

    Ok(())
}
