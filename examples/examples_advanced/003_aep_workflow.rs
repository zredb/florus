use florus::core::Farm;
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;
use florus::aep::calculate_aep_from_time_series;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 19: Complete AEP Workflow");
    println!("====================================\n");

    let d = 126.0;
    let n_turbines = 10;
    let spacing = 7.0 * d;

    let layout_x: Vec<f64> = (0..n_turbines).map(|i| i as f64 * spacing).collect();
    let layout_y = vec![0.0; n_turbines];
    let turbine_types = vec!["nrel_5MW".to_string(); n_turbines];

    println!("Complete AEP Workflow:");
    println!("  Farm: {} turbines at {:.1}D spacing", n_turbines, spacing / d);
    println!("  Rated power: {:.0} MW\n", n_turbines as f64 * 5.0);

    let farm = Farm::new(
        Array1::from_vec(layout_x),
        Array1::from_vec(layout_y),
        turbine_types,
    )?;

    println!("--- Step 1: Define Wind Rose ---\n");

    let wind_directions = Array1::from_vec(vec![0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]);
    let wind_speeds = Array1::from_vec(vec![4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0]);

    let freq_table = Array2::from_shape_vec((12, 10), vec![
        0.01, 0.02, 0.03, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.00,
        0.01, 0.02, 0.03, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.00,
        0.02, 0.03, 0.04, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.00,
        0.02, 0.03, 0.04, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.00,
        0.01, 0.02, 0.03, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.00,
        0.01, 0.01, 0.02, 0.03, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00,
        0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.00, 0.00,
        0.01, 0.02, 0.03, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.00,
        0.02, 0.04, 0.06, 0.08, 0.06, 0.04, 0.02, 0.02, 0.01, 0.00,
        0.03, 0.05, 0.07, 0.08, 0.07, 0.05, 0.03, 0.02, 0.01, 0.00,
        0.02, 0.04, 0.05, 0.06, 0.05, 0.03, 0.02, 0.01, 0.01, 0.00,
        0.01, 0.02, 0.03, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.00,
    ])?;

    let ti_table = Array2::from_elem((12, 10), 0.06);

    let wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        ti_table,
        Some(freq_table),
        None,
        false,
        None,
        None,
    )?;

    println!("Wind Rose:");
    println!("  {} directions x {} speeds", wind_directions.len(), wind_speeds.len());

    let total_freq: f64 = wind_rose.freq_table.iter().sum();
    println!("  Frequency sum: {:.2}", total_freq);

    println!("\n--- Step 2: Calculate AEP from Time Series ---\n");

    // Create a simple time series from wind rose
    let wd_vec = wind_directions.iter().cloned().collect::<Vec<f64>>();
    let ws_vec = vec![8.0; 12];
    let ti_vec = vec![0.06; 12];

    let time_series = TimeSeries::new(
        Array1::from_vec(wd_vec),
        Array1::from_vec(ws_vec),
        Array1::from_vec(ti_vec),
    )?;

    let aep_result = calculate_aep_from_time_series(&farm, &time_series, None);

    println!("AEP Calculation Result:");
    println!("  Annual energy: {:.2} GWh", aep_result.total_energy_mwh / 1000.0);
    println!("  Conditions processed: {}", aep_result.conditions_processed);

    println!("\n--- Step 3: Sample Calculations ---\n");

    println!("Sample power calculations:");
    for i in [0, 4, 8] {
        let ws = wind_speeds[i];
        let flow_field = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.12, 1.225,
            Array1::from_vec(vec![0.06]),
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

        let powers = model.get_turbine_powers();
        let total: f64 = (0..n_turbines).map(|ti| powers[[0, ti]]).sum();
        println!("  {} m/s: {:.1} MW", ws, total / 1_000_000.0);
    }

    println!("\n--- Step 4: Summary ---\n");

    println!("AEP Workflow Summary:");
    println!("  1. Define wind conditions (WindRose or TimeSeries)");
    println!("  2. Create FLORIS model with farm layout");
    println!("  3. Calculate AEP using calculate_aep_from_time_series()");
    println!("  4. Analyze results (energy, capacity factor)");

    let rated_capacity = n_turbines as f64 * 5.0;
    println!("\nFarm Statistics:");
    println!("  Rated capacity: {:.1} MW", rated_capacity);
    println!("  Expected annual energy: {:.2} GWh", aep_result.total_energy_mwh / 1000.0);

    println!("\n====================================");
    println!("Example completed successfully!");
    Ok(())
}
