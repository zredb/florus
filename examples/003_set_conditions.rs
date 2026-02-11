/// Example 9: Wind Conditions and Time Series Analysis
///
/// This example demonstrates working with different wind data types
/// and time series analysis for wind farm simulations.
///
/// Topics covered:
/// 1. TimeSeries for sequential wind conditions
/// 2. WindRose for aggregated wind statistics
/// 3. Wind direction and speed distributions
/// 4. Turbulence intensity effects
/// 5. Batch simulation workflows

use florus::core::Farm;
use florus::types::Array1;
use florus::wind_data::{TimeSeries, WindData};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 9: Wind Conditions and Time Series");
    println!("============================================\n");

    let d = 126.0; // NREL 5MW rotor diameter

    // Create a 3-turbine farm
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    println!("Creating 3-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m", i, x);
    }

    // ============================================================
    // TimeSeries Analysis
    // ============================================================
    println!("\n--- Time Series Analysis ---\n");

    // Create a time series with varying conditions
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0, 270.0, 270.0, 280.0, 280.0, 280.0, 290.0, 290.0, 290.0]),
        Array1::from_vec(vec![8.0, 10.0, 12.0, 8.0, 10.0, 12.0, 8.0, 10.0, 12.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06]),
    )?;

    println!("TimeSeries created with {} conditions:", time_series.n_conditions());

    // Run simulation for each condition
    println!("\nSimulating power output for each condition:");
    println!("  {:>8}  {:>8}  {:>8}  {:>12}", "WD (°)", "WS (m/s)", "TI", "Farm (MW)");
    println!("  {}", "-".repeat(45));

    let mut total_energy = 0.0;

    for i in 0..time_series.n_conditions() {
        let wd = time_series.wind_directions[i];
        let ws = time_series.wind_speeds[i];
        let ti = time_series.turbulence_intensities[i];

        let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

        let flow_field = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![wd]),
            0.0,
            0.14,
            1.225,
            Array1::from_vec(vec![ti]),
            90.0,
        )?;

        let mut model = florus::FlorisModel {
            farm,
            flow_field,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };

        model.initialize_grid()?;
        model.initialize_flow_field()?;
        model.run()?;

        let powers = model.get_turbine_powers();
        let farm_power: f64 = powers.iter().sum();

        println!("  {:>8.0}  {:>8.1}  {:>8.2}  {:>10.3}", wd, ws, ti, farm_power / 1_000_000.0);

        total_energy += farm_power; // Assuming 1 hour per condition
    }

    println!("\nTotal energy (assuming 1 hour per condition): {:.2} MWh", total_energy / 1_000_000.0);

    // ============================================================
    // Wind Direction Distribution
    // ============================================================
    println!("\n--- Wind Direction Distribution ---\n");

    // Create a typical wind rose distribution
    let wind_directions: Vec<f64> = (0..36).map(|i| (i as f64) * 10.0).collect();
    let wind_speeds = vec![8.0, 10.0, 12.0]; // Typical speeds

    // Create frequency distribution (more wind from SW)
    let mut freq_2d = Vec::new();
    for &wd in &wind_directions {
        for &ws in &wind_speeds {
            let wd_factor = if wd >= 180.0 && wd <= 270.0 { 0.4 } else { 0.1 };
            let ws_factor = if ws >= 8.0 && ws <= 12.0 { 0.15 } else { 0.05 };
            freq_2d.push(wd_factor * ws_factor);
        }
    }

    // Normalize
    let sum: f64 = freq_2d.iter().sum();
    for f in freq_2d.iter_mut() {
        *f /= sum;
    }

    println!("Wind direction bins: {} ({}° to {}°)",
             wind_directions.len(),
             wind_directions.first().unwrap(),
             wind_directions.last().unwrap());
    println!("Wind speed bins: {} ({} m/s to {} m/s)",
             wind_speeds.len(),
             wind_speeds.first().unwrap(),
             wind_speeds.last().unwrap());

    // ============================================================
    // Turbulence Intensity Effects
    // ============================================================
    println!("\n--- Turbulence Intensity Effects ---\n");

    let ti_values = vec![0.03, 0.06, 0.08, 0.10, 0.12, 0.14];

    println!("Testing power output at different turbulence intensities:");
    println!("  {:>8}  {:>10}  {:>10}", "TI", "Farm (MW)", "Wake Loss (%)");
    println!("  {}", "-".repeat(35));

    let farm_base = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    let flow_field_base = florus::core::FlowField::new(
        Array1::from_vec(vec![10.0]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model_base = florus::FlorisModel {
        farm: farm_base,
        flow_field: flow_field_base,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_base.initialize_grid()?;
    model_base.initialize_flow_field()?;
    model_base.run()?;
    let powers_base = model_base.get_turbine_powers();
    let base_power: f64 = powers_base.iter().sum();

    for &ti in &ti_values {
        let farm_ti = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

        let flow_field_ti = florus::core::FlowField::new(
            Array1::from_vec(vec![10.0]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.14, 1.225,
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

        let powers_ti = model_ti.get_turbine_powers();
        let farm_power_ti: f64 = powers_ti.iter().sum();

        let wake_loss = (1.0 - farm_power_ti / base_power) * 100.0;

        println!("  {:>8.2}  {:>10.3}  {:>10.1}%", ti, farm_power_ti / 1_000_000.0, wake_loss);
    }

    // ============================================================
    // Wind Speed Power Curve
    // ============================================================
    println!("\n--- Wind Speed Power Curve ---\n");

    let wind_speeds: Vec<f64> = (3..26).map(|i| i as f64).collect();

    println!("Computing power curve (3-turbine farm, 270°):");
    println!("  {:>8}  {:>10}  {:>10}", "WS (m/s)", "T0 (MW)", "Farm (MW)");
    println!("  {}", "-".repeat(35));

    for &ws in wind_speeds.iter().step_by(2) {
        let farm_ws = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

        let flow_field_ws = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.14, 1.225,
            Array1::from_vec(vec![0.06]),
            90.0,
        )?;

        let mut model_ws = florus::FlorisModel {
            farm: farm_ws,
            flow_field: flow_field_ws,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };

        model_ws.initialize_grid()?;
        model_ws.initialize_flow_field()?;
        model_ws.run()?;

        let powers_ws = model_ws.get_turbine_powers();
        let farm_power_ws: f64 = powers_ws.iter().sum();

        println!("  {:>8.1}  {:>10.3}  {:>10.3}", ws, powers_ws[[0, 0]] / 1_000_000.0, farm_power_ws / 1_000_000.0);
    }

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Wind Data Types in FLORIS:");
    println!("  1. TimeSeries:");
    println!("     - Sequential wind measurements");
    println!("     - Each time step has WD, WS, TI");
    println!("     - Good for temporal analysis");
    println!();
    println!("  2. WindRose:");
    println!("     - Aggregated wind statistics");
    println!("     - Binned by direction and speed");
    println!("     - Includes frequency information");
    println!("     - Essential for AEP calculations");
    println!();

    println!("Key Parameters:");
    println!("  - Wind Speed: Affects power output (cubic relationship)");
    println!("  - Wind Direction: Determines wake interactions");
    println!("  - Turbulence Intensity: Affects wake recovery");
    println!("  - Wind Shear: Vertical speed profile");
    println!("  - Wind Veer: Direction change with height");

    println!("\n============================================");
    println!("Example completed successfully!");

    Ok(())
}
