use florus::core::Farm;
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 14: Wind Rose Analysis");
    println!("====================================\n");

    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 3]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    println!("Creating 3-turbine wind farm:");
    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    println!("\n--- Wind Rose Configuration ---\n");

    let wind_directions = Array1::from_vec(vec![
        0.0, 30.0, 60.0, 90.0, 120.0, 150.0,
        180.0, 210.0, 240.0, 270.0, 300.0, 330.0,
    ]);

    let wind_speeds = Array1::from_vec(vec![5.0, 7.5, 10.0, 12.5]);

    let freq_table = Array2::from_shape_vec((12, 4), vec![
        0.02, 0.03, 0.02, 0.01,
        0.03, 0.04, 0.02, 0.01,
        0.04, 0.05, 0.03, 0.01,
        0.03, 0.04, 0.02, 0.01,
        0.02, 0.03, 0.02, 0.01,
        0.01, 0.02, 0.01, 0.00,
        0.01, 0.02, 0.01, 0.00,
        0.02, 0.03, 0.02, 0.01,
        0.05, 0.08, 0.05, 0.02,
        0.08, 0.10, 0.06, 0.03,
        0.06, 0.08, 0.04, 0.02,
        0.03, 0.04, 0.02, 0.01,
    ])?;

    let ti_table = Array2::from_elem((12, 4), 0.06);

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

    println!("Wind Rose Configuration:");
    println!("  Wind directions: {} sectors", wind_directions.len());
    println!("  Wind speeds: {} bins", wind_speeds.len());
    println!("  Total sectors: {}", wind_directions.len() * wind_speeds.len());

    println!("\n--- Wind Direction Frequency ---\n");

    println!("{:>8} {:>10}", "Dir (deg)", "Freq");
    println!("{}", "-".repeat(20));
    for i in 0..wind_directions.len() {
        let freq: f64 = (0..wind_speeds.len()).map(|j| wind_rose.freq_table[[i, j]]).sum();
        println!("{:>8.0} {:>10.1}%", wind_directions[i], freq * 100.0);
    }

    println!("\n--- Running Simulations ---\n");

    let mut sector_powers: Vec<f64> = Vec::new();

    for i in 0..wind_directions.len() {
        let dominant_ws = 10.0;
        let ti = wind_rose.ti_table[[i, 2]];

        let flow_field = florus::core::FlowField::new(
            Array1::from_vec(vec![dominant_ws]),
            Array1::from_vec(vec![wind_directions[i]]),
            0.0, 0.12, 1.225,
            Array1::from_vec(vec![ti]),
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
        let total: f64 = (0..3).map(|ti| powers[[0, ti]]).sum();
        sector_powers.push(total / 1_000_000.0);
    }

    println!("Power by Wind Direction:");
    println!("{:>8} {:>12}", "Dir (deg)", "Power (MW)");
    println!("{}", "-".repeat(22));
    for i in 0..wind_directions.len() {
        println!("{:>8.0} {:>12.3}", wind_directions[i], sector_powers[i]);
    }

    let weighted_power: f64 = (0..wind_directions.len())
        .map(|i| {
            let freq: f64 = (0..wind_speeds.len()).map(|j| wind_rose.freq_table[[i, j]]).sum();
            freq * sector_powers[i]
        }).sum();

    println!("\n--- Summary ---\n");
    println!("Wind Rose Analysis:");
    println!("  Weighted average power: {:.3} MW", weighted_power);
    println!("  Dominant wind direction: 270 deg (onshore)");

    println!("\n====================================");
    println!("Example completed successfully!");
    Ok(())
}
