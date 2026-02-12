use florus::{FlorisModel, Array1, Result};
use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::wind_data::{TimeSeries, WindData};

fn main() -> Result<()> {
    println!("Testing set_wind_data method");
    println!("==============================\n");

    let layout_x = Array1::from_vec(vec![0.0, 630.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 2];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    let time_series = TimeSeries::new(
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.07, 0.08]),
    )?;

    let mut model = FlorisModel {
        farm,
        flow_field: florus::core::FlowField::new(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            0.0,
            0.14,
            1.225,
            Array1::from_vec(vec![0.06]),
            90.0,
        )?,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    println!("Using set_wind_data() with TimeSeries...");
    model.set_wind_data(&time_series)?;

    println!("Wind data set successfully!");
    println!("  n_findex: {}", model.flow_field.n_findex);
    println!("  n_conditions: {}", time_series.n_conditions());

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let powers = model.get_turbine_powers();
    println!("\nResults:");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, powers[[0, ti]] / 1000.0);
    }

    println!("\nTest completed successfully!");

    Ok(())
}
