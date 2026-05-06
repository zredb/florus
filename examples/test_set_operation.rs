//! Example: Test set_operation method
//!
//! Demonstrate the use of set_operation method in FLORUS

use florus::{FlorisModel, Result};
use ndarray::Array2;

fn main() -> Result<()> {
    println!("=== FLORUS set_operation Method Test ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 3-turbine layout
    let layout_x = ndarray::Array::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = ndarray::Array::from_vec(vec![0.0, 0.0, 0.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;

    // Set wind conditions
    let wind_speeds = ndarray::Array::from_vec(vec![8.0]);
    let wind_directions = ndarray::Array::from_vec(vec![270.0]);
    fmodel.set(
        Some(wind_speeds),
        Some(wind_directions),
        None,  // wind_shear
        None,  // wind_veer
        None,  // reference_wind_height
        None,  // turbulence_intensities
        None,  // air_density
        None,  // layout_x
        None,  // layout_y
        None,  // yaw_angles
        None,  // power_setpoints
        None,  // awc_modes
        None,  // awc_amplitudes
        None,  // awc_frequencies
        None,  // disable_turbines
    )?;

    println!("Initial state:");
    println!("  Turbines: {}", fmodel.core.farm.n_turbines());
    println!("  Findex: {}", fmodel.core.flow_field.n_findex);
    println!("  Yaw angles: {:?}", fmodel.core.farm.yaw_angles);
    println!("  Power setpoints: {:?}", fmodel.core.farm.power_setpoints);
    println!("  AWC modes: {:?}", fmodel.core.farm.awc_modes);
    println!();

    // Test 1: Set yaw angles
    println!("Test 1: Setting yaw angles");
    let yaw_angles = Array2::from_shape_vec((1, 3), vec![0.0, 15.0, 30.0])?;
    fmodel.set_operation(
        Some(yaw_angles),
        None,  // power_setpoints
        None,  // awc_modes
        None,  // awc_amplitudes
        None,  // awc_frequencies
        None,  // disable_turbines
    )?;
    println!("  Yaw angles after set: {:?}", fmodel.core.farm.yaw_angles);
    println!();

    // Test 2: Set power setpoints
    println!("Test 2: Setting power setpoints");
    let power_setpoints = Array2::from_shape_vec((1, 3), vec![5000000.0, 4500000.0, 4000000.0])?;
    fmodel.set_operation(
        None,  // yaw_angles
        Some(power_setpoints),
        None,  // awc_modes
        None,  // awc_amplitudes
        None,  // awc_frequencies
        None,  // disable_turbines
    )?;
    println!("  Power setpoints after set: {:?}", fmodel.core.farm.power_setpoints);
    println!();

    // Test 3: Set AWC modes
    println!("Test 3: Setting AWC modes");
    use florus::types::NdArray2;
    let awc_modes = NdArray2::from_shape_vec((1, 3), vec!["baseline".to_string(), "helix".to_string(), "baseline".to_string()])?;
    fmodel.set_operation(
        None,  // yaw_angles
        None,  // power_setpoints
        Some(awc_modes),
        None,  // awc_amplitudes
        None,  // awc_frequencies
        None,  // disable_turbines
    )?;
    println!("  AWC modes after set: {:?}", fmodel.core.farm.awc_modes);
    println!();

    // Test 4: Set AWC amplitudes and frequencies
    println!("Test 4: Setting AWC amplitudes and frequencies");
    let awc_amplitudes = Array2::from_shape_vec((1, 3), vec![0.0, 5.0, 0.0])?;
    let awc_frequencies = Array2::from_shape_vec((1, 3), vec![0.0, 0.2, 0.0])?;
    fmodel.set_operation(
        None,  // yaw_angles
        None,  // power_setpoints
        None,  // awc_modes
        Some(awc_amplitudes),
        Some(awc_frequencies),
        None,  // disable_turbines
    )?;
    println!("  AWC amplitudes after set: {:?}", fmodel.core.farm.awc_amplitudes);
    println!("  AWC frequencies after set: {:?}", fmodel.core.farm.awc_frequencies);
    println!();

    // Test 5: Disable turbines
    println!("Test 5: Disabling turbines");
    let disable_turbines = NdArray2::from_shape_vec((1, 3), vec![false, true, false])?;
    fmodel.set_operation(
        None,  // yaw_angles
        None,  // power_setpoints
        None,  // awc_modes
        None,  // awc_amplitudes
        None,  // awc_frequencies
        Some(disable_turbines),
    )?;
    println!("  Yaw angles after disable: {:?}", fmodel.core.farm.yaw_angles);
    println!("  Power setpoints after disable: {:?}", fmodel.core.farm.power_setpoints);
    println!();

    // Test 6: Set all parameters at once
    println!("Test 6: Setting all parameters at once");
    let yaw_angles = Array2::from_shape_vec((1, 3), vec![5.0, 10.0, 15.0])?;
    let power_setpoints = Array2::from_shape_vec((1, 3), vec![4800000.0, 4600000.0, 4400000.0])?;
    let awc_modes = NdArray2::from_shape_vec((1, 3), vec!["baseline".to_string(); 3])?;
    let awc_amplitudes = Array2::from_shape_vec((1, 3), vec![2.0, 3.0, 4.0])?;
    let awc_frequencies = Array2::from_shape_vec((1, 3), vec![0.1, 0.15, 0.2])?;
    
    fmodel.set_operation(
        Some(yaw_angles),
        Some(power_setpoints),
        Some(awc_modes),
        Some(awc_amplitudes),
        Some(awc_frequencies),
        None,  // disable_turbines
    )?;
    println!("  Yaw angles: {:?}", fmodel.core.farm.yaw_angles);
    println!("  Power setpoints: {:?}", fmodel.core.farm.power_setpoints);
    println!("  AWC modes: {:?}", fmodel.core.farm.awc_modes);
    println!("  AWC amplitudes: {:?}", fmodel.core.farm.awc_amplitudes);
    println!("  AWC frequencies: {:?}", fmodel.core.farm.awc_frequencies);
    println!();

    println!("=== set_operation Method Test Complete! ===");

    Ok(())
}
