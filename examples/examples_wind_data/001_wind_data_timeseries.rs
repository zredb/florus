use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::types::Array1;
use florus::wind_data::{TimeSeries, WindData};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 13: Time Series Wind Data");
    println!("==========================================\n");

    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d, 15.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 4]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    println!("Creating 4-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m", i, x);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("\n--- Creating Time Series ---\n");

    let wind_directions: Vec<f64> = vec![
        270.0, 268.0, 265.0, 260.0, 255.0, 250.0,
        248.0, 250.0, 255.0, 260.0, 265.0, 270.0,
        272.0, 275.0, 278.0, 280.0, 282.0, 283.0,
        282.0, 280.0, 278.0, 275.0, 272.0, 270.0,
    ];

    let wind_speeds: Vec<f64> = vec![
        6.0, 5.5, 5.0, 5.0, 5.5, 6.0,
        7.0, 8.5, 10.0, 11.0, 11.5, 12.0,
        12.0, 11.5, 11.0, 10.0, 9.0, 8.0,
        7.5, 7.0, 6.5, 6.0, 6.0, 6.0,
    ];

    let turbulence_intensities: Vec<f64> = vec![
        0.10, 0.11, 0.12, 0.12, 0.11, 0.10,
        0.08, 0.07, 0.06, 0.06, 0.06, 0.05,
        0.05, 0.05, 0.06, 0.06, 0.07, 0.08,
        0.08, 0.09, 0.10, 0.10, 0.11, 0.11,
    ];

    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions),
        Array1::from_vec(wind_speeds),
        Array1::from_vec(turbulence_intensities),
    )?;

    println!("Time Series Configuration:");
    println!("  Number of time steps: {}", time_series.n_conditions());

    println!("\n--- Running Time Series Simulation ---\n");

    let n_times = time_series.n_conditions();
    let mut hourly_powers: Vec<f64> = Vec::new();

    println!("Simulating {} time steps...", n_times);

    for i in 0..n_times {
        let ws = time_series.wind_speeds[i];
        let wd = time_series.wind_directions[i];
        let ti = time_series.turbulence_intensities[i];

        let flow_field = florus::core::FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![wd]),
            0.0,
            0.12,
            1.225,
            Array1::from_vec(vec![ti]),
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

        let powers = model.get_turbine_powers();
        let total_power: f64 = (0..4).map(|ti| powers[[0, ti]]).sum();
        hourly_powers.push(total_power / 1_000_000.0);

        if i % 6 == 0 {
            println!("  Hour {}: {:.3} MW", i, total_power / 1_000_000.0);
        }
    }

    println!("\n--- Time Series Analysis ---\n");

    let max_power = hourly_powers.iter().cloned().fold(0.0 / 0.0, f64::max);
    let avg_power: f64 = hourly_powers.iter().sum::<f64>() / hourly_powers.len() as f64;
    let total_daily_energy = hourly_powers.iter().sum::<f64>();

    println!("Daily Power Statistics:");
    println!("  Maximum power: {:.3} MW", max_power);
    println!("  Average power: {:.3} MW", avg_power);
    println!("  Daily energy: {:.2} MWh", total_daily_energy);

    println!("\n--- Time Series with Values ---\n");

    let prices: Vec<f64> = vec![
        0.05, 0.05, 0.05, 0.05, 0.05, 0.06,
        0.08, 0.12, 0.15, 0.18, 0.20, 0.22,
        0.25, 0.28, 0.30, 0.28, 0.25, 0.22,
        0.18, 0.15, 0.12, 0.10, 0.08, 0.06,
    ];

    let time_series_priced = TimeSeries::with_values(
        Array1::from_vec(vec![270.0; 24]),
        Array1::from_vec(vec![8.0; 24]),
        Array1::from_vec(vec![0.06; 24]),
        Array1::from_vec(prices),
    )?;

    let total_revenue: f64 = hourly_powers.iter().zip(time_series_priced.values.iter())
        .map(|(p, v)| p * 1000.0 * v).sum();

    println!("Time Series with Electricity Prices:");
    println!("  Estimated daily revenue: ${:.2}", total_revenue);

    println!("\n==========================================");
    println!("Example completed successfully!");

    Ok(())
}
