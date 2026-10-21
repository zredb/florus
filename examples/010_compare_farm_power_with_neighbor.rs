//! Example: Compare the farm power with a neighboring farm
//!
//! This example computes the power of a wind farm with and without a
//! neighboring farm to the east: the turbines of the neighboring farm are
//! given a weight of zero, so that only the power of the base farm is
//! reported, while the wakes of the neighboring farm are still modeled.
//!
//! Scenario aligned 1:1 with examples/010_compare_farm_power_with_neighbor.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A four turbine farm
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, d * 6.0, 0.0, d * 6.0]),
        &Array1::from_vec(vec![0.0, 0.0, d * 3.0, d * 3.0]),
    )?;

    // A simple inflow with a single wind speed
    let wind_directions: Vec<f64> = (0..90).map(|i| 4.0 * i as f64).collect();
    let n = wind_directions.len();
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;

    fmodel.run()?;
    let farm_power_base = fmodel.get_farm_power()?;

    // Add a neighbor to the east
    let layout_x = Array1::from_vec(vec![
        0.0,
        d * 6.0,
        0.0,
        d * 6.0,
        d * 12.0,
        d * 15.0,
        d * 12.0,
        d * 15.0,
    ]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, d * 3.0, d * 3.0, 0.0, 0.0, d * 3.0, d * 3.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;

    // The weights exclude the neighboring farm from the power calculation
    let mut turbine_weights = Array2::zeros((n, 8));
    for fi in 0..n {
        for t in 0..4 {
            turbine_weights[[fi, t]] = 1.0;
        }
    }

    fmodel.run()?;
    let farm_power_neighbor = fmodel.get_farm_power_weighted(Some(&turbine_weights))?;

    println!(
        "Wind direction [deg], farm power [kW], farm power with the neighbor [kW], difference [kW]"
    );
    for i in 0..n {
        println!(
            "{:.0}, {:.1}, {:.1}, {:.1}",
            wind_directions[i],
            farm_power_base[i] / 1000.0,
            farm_power_neighbor[i] / 1000.0,
            (farm_power_neighbor[i] - farm_power_base[i]) / 1000.0
        );
    }

    Ok(())
}
