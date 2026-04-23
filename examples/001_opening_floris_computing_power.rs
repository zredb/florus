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
    // Print the full array similar to Python's numpy output
    for i in 0..turbine_powers.shape()[0] {
        if i == 0 {
            print!("[");
        } else {
            print!(" ");
        }
        print!("[");
        for j in 0..turbine_powers.shape()[1] {
            if j > 0 {
                print!(" ");
            }
            print!("{:12.8}", turbine_powers[[i, j]]);
        }
        print!("]");
        if i < turbine_powers.shape()[0] - 1 {
            println!();
        } else {
            println!("]");
        }
    }
    println!("Shape:  ({}, {})", turbine_powers.shape()[0], turbine_powers.shape()[1]);

    println!("\nThe farm power should be a 1D array of length 4 (n_findex)");
    print!("[");
    for i in 0..farm_power.shape()[0] {
        if i > 0 {
            print!(" ");
        }
        print!("{:12.8}", farm_power[i]);
    }
    println!("]");
    println!("Shape:  ({},)", farm_power.shape()[0]);

    Ok(())
}