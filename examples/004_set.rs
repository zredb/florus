/// Example 4: Set
///
/// This example illustrates the use of the set method. The set method is used to
/// change the wind conditions, the wind farm layout, and the controls settings.
///
/// This example demonstrates setting each of the following:
///     1) Wind conditions
///     2) Wind farm layout
///     3) Controls settings

use florus::{Array1, Array2, FlorisModel};
use florus::core::turbines::TurbineLibrary;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    //////////////////////////////////////////////////
    // Atmospheric Conditions
    //////////////////////////////////////////////////

    // Change the wind directions, wind speeds, and turbulence intensities using arrays
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06]),
    )?;

    println!("Set wind conditions: 3 conditions (8, 9, 10 m/s at 270°)");
    println!("  n_findex: {}", fmodel.n_findex());

    // Set the wind shear
    fmodel.set_wind_shear(0.2)?;
    println!("Set wind shear: 0.2");

    // Set the air density
    fmodel.set_air_density(1.1)?;
    println!("Set air density: 1.1 kg/m³");

    // Set the reference wind height
    fmodel.set_reference_wind_height(92.0)?;
    println!("Set reference wind height: 92.0 m");

    //////////////////////////////////////////////////
    // Array Settings
    //////////////////////////////////////////////////

    // Changing the wind farm layout uses FLORIS' set method to a two-turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    println!("\nSet layout: 2 turbines at [0, 500]m");
    println!("  n_turbines: {}", fmodel.n_turbines());

    //////////////////////////////////////////////////
    // Controls Settings
    //////////////////////////////////////////////////

    // Changes to controls settings can be made using the set method
    // Note the dimension must match (n_findex, n_turbines)
    // Above we have n_findex = 3 and n_turbines = 2 so the matrix of yaw angles must be 3x2
    let yaw_angles = Array2::from_shape_vec(
        (3, 2),
        vec![
            0.0, 0.0,   // Condition 0: both turbines at 0°
            25.0, 0.0,  // Condition 1: front turbine at 25°, rear at 0°
            0.0, 0.0,   // Condition 2: both turbines at 0°
        ],
    ).expect("Failed to create yaw_angles array");
    println!("\nSet yaw angles: 3x2 matrix");
    println!("  Yaw angles shape: {:?}", yaw_angles.shape());
    fmodel.set_yaw_angles(yaw_angles)?;

    // Use the reset operation method to clear out control signals
    fmodel.reset_operation();
    println!("Reset operation (cleared control signals)");

    println!("\nExample 4 completed successfully!");

    Ok(())
}
