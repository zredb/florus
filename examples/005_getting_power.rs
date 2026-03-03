/// Example 2: Getting Turbine and Farm Power
///
/// After setting the FlorisModel and running, the next step is typically to get the power output
/// of the turbines. FLORIS-RS provides several methods for getting power:
///
/// 1. `get_turbine_powers()`: Returns the power output of each turbine in the farm for each findex
///     Shape: (n_findex, n_turbines)
/// 2. `get_farm_power()`: Returns the total power output of the farm for each findex
///     Shape: (n_findex,)
/// 3. `get_farm_aep()`: Calculates Annual Energy Production
///
/// This example demonstrates:
/// - Loading config from YAML file (matching Python FLORIS)
/// - Setting custom wind farm layout
/// - Wind direction sweeping using array-based conditions
/// - Getting turbine and farm power
/// - Computing wake losses by comparing with no-wake case
///
/// This is the Rust equivalent of Python's 005_getting_power.py
///

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 2: Getting Turbine and Farm Power");
    println!("=====================================================\n");

    // ============================================================
    // Load model from YAML config (matching Python FLORIS behavior)
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Create a 3-turbine farm
    // In Python: fmodel.set(layout_x=[0, 126 * 5, 126 * 10], layout_y=[0, 0, 0])
    let d = 126.0; // NREL 5MW rotor diameter
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    
    println!("Creating 3-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, layout_y[i]);
    }
    
    // Set custom layout
    model.set_layout(&layout_x, &layout_y)?;

    // ============================================================
    // Using array-based wind conditions (sweeping wind directions)
    // ============================================================
    println!("\n--- Wind Direction Sweep (250° to 290°) ---");

    // In Python:
    // wind_directions = np.arange(250, 290, 1.0)
    // time_series = TimeSeries(
    //     wind_directions=wind_directions, wind_speeds=9.9, turbulence_intensities=0.06
    // )

    // Create wind direction sweep
    let wind_directions: Vec<f64> = (250..290).map(|d| d as f64).collect();
    let wind_speeds = Array1::from_vec(vec![9.9; 40]); // Constant wind speed for all conditions
    let turbulence_intensities = Array1::from_vec(vec![0.06; 40]); // Constant TI

    println!("Simulating {} wind directions from {}° to {}",
             wind_directions.len(),
             wind_directions.first().unwrap(),
             wind_directions.last().unwrap());
    
    // Set wind conditions
    model.set_wind_conditions(
        wind_speeds,
        Array1::from_vec(wind_directions.clone()),
        turbulence_intensities,
    )?;

    // ============================================================
    // Run simulation
    // ============================================================
    model.run()?;

    // Get powers with wake
    let turbine_powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();

    println!("\nResults with wake:");
    println!("  Turbine powers shape: ({}, {})",
             turbine_powers.shape()[0], turbine_powers.shape()[1]);
    println!("  Farm power shape: ({},)", farm_power.shape()[0]);

    // ============================================================
    // Compute wake losses (compare with no-wake case)
    // ============================================================
    println!("\n--- Wake Loss Analysis ---");

    // For no-wake case, we would need to run with a special no-wake mode
    // For now, let's compute relative losses between turbines
    println!("\nFarm power by wind direction:");
    println!("  {:>6}  {:>10}", "WD (°)", "Farm Power (kW)");
    println!("  {}", "-".repeat(20));

    for (i, wd) in wind_directions.iter().enumerate().step_by(5) {
        println!("  {:>6.0}  {:>10.1}", wd, farm_power[[i]] / 1000.0);
    }

    // ============================================================
    // Analyze turbine power distribution
    // ============================================================
    println!("\n--- Turbine Power Analysis ---");

    // Find conditions with max/min farm power
    let mut max_power_idx = 0;
    let mut min_power_idx = 0;
    for i in 1..farm_power.len() {
        if farm_power[[i]] > farm_power[[max_power_idx]] {
            max_power_idx = i;
        }
        if farm_power[[i]] < farm_power[[min_power_idx]] {
            min_power_idx = i;
        }
    }

    println!("\nBest condition: {}° wind direction", wind_directions[max_power_idx]);
    println!("  Turbine powers (kW):");
    for ti in 0..model.farm.n_turbines() {
        println!("    Turbine {}: {:.1} kW", ti, turbine_powers[[max_power_idx, ti]] / 1000.0);
    }
    println!("  Total farm power: {:.1} kW", farm_power[[max_power_idx]] / 1000.0);

    println!("\nWorst condition: {}° wind direction", wind_directions[min_power_idx]);
    println!("  Turbine powers (kW):");
    for ti in 0..model.farm.n_turbines() {
        println!("    Turbine {}: {:.1} kW", ti, turbine_powers[[min_power_idx, ti]] / 1000.0);
    }
    println!("  Total farm power: {:.1} kW", farm_power[[min_power_idx]] / 1000.0);

    // ============================================================
    // Wake loss calculation (simplified)
    // ============================================================
    println!("\n--- Simplified Wake Loss Analysis ---");

    // For aligned wind (270°), downstream turbines should have lower power
    let aligned_idx = wind_directions.iter().position(|&d| (d - 270.0).abs() < 0.5).unwrap();

    println!("\nAt 270° (aligned wind):");
    let upstream_power = turbine_powers[[aligned_idx, 0]] / 1000.0;
    let middle_power = turbine_powers[[aligned_idx, 1]] / 1000.0;
    let downstream_power = turbine_powers[[aligned_idx, 2]] / 1000.0;

    println!("  Turbine 0 (upstream): {:.1} kW", upstream_power);
    println!("  Turbine 1 (middle):    {:.1} kW", middle_power);
    println!("  Turbine 2 (downstream): {:.1} kW", downstream_power);

    let wake_loss_1 = (1.0 - middle_power / upstream_power) * 100.0;
    let wake_loss_2 = (1.0 - downstream_power / upstream_power) * 100.0;

    println!("\n  Wake loss T0→T1: {:.1}%", wake_loss_1);
    println!("  Wake loss T0→T2: {:.1}%", wake_loss_2);
    println!("  Combined wake effect on T2: {:.1}%", wake_loss_2);

    // ============================================================
    // Comparison with offset wind direction (less wake interaction)
    // ============================================================
    println!("\n--- Wind Direction Effect on Wake Losses ---");

    let offset_angles = [250.0, 260.0, 270.0, 280.0, 290.0];

    println!("\n  {:>6}  {:>10}  {:>10}  {:>10}", "WD (°)", "T0 (kW)", "T1 (kW)", "T2 (kW)");
    println!("  {}", "-".repeat(45));

    for &wd in &offset_angles {
        if let Some(idx) = wind_directions.iter().position(|&d| (d - wd).abs() < 0.5) {
            let p0 = turbine_powers[[idx, 0]] / 1000.0;
            let p1 = turbine_powers[[idx, 1]] / 1000.0;
            let p2 = turbine_powers[[idx, 2]] / 1000.0;
            println!("  {:>6.0}  {:>10.1}  {:>10.1}  {:>10.1}", wd, p0, p1, p2);
        }
    }

    println!("\n=======================================================");
    println!("Example completed successfully!");

    Ok(())
}
