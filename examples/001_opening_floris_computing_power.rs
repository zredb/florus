use florus::{Array1, FlorisModel};
use florus::core::turbines::TurbineLibrary;
use ndarray::Axis;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 630.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 8.0, 10.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06]),
    )?;

    fmodel.run()?;

    let turbine_powers = fmodel.get_turbine_powers() / 1000.0;
    let farm_power = fmodel.get_farm_power() / 1000.0;

    println!("The turbine power matrix should be of dimensions 4 (n_findex) X 2 (n_turbines)");
    for (i, row) in turbine_powers.axis_iter(Axis(0)).enumerate() {
        println!("  findex {}: {:?}", i, row.as_slice().unwrap());
    }
    println!(
        "Shape: ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );

    println!("\nThe farm power should be a 1D array of length 4 (n_findex)");
    println!("Farm power (kW): {:?}", farm_power.as_slice().unwrap());
    println!("Shape: ({},)", farm_power.shape()[0]);

    Ok(())
}