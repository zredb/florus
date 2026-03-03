/// Debug script to check velocities in Rust FLORIS
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    // Load model from YAML config
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Create a 3-turbine farm
    let layout_x = Array1::from_vec(vec![0.0, 630.0, 1260.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    
    // Set custom layout
    model.set_layout(&layout_x, &layout_y)?;

    // Set single wind condition
    let wind_speeds = Array1::from_vec(vec![9.9]);
    let wind_directions = Array1::from_vec(vec![270.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    model.set_wind_conditions(
        wind_speeds,
        wind_directions,
        turbulence_intensities,
    )?;

    // Run simulation
    model.run()?;

    // Get grid coordinates
    let grid = model.grid.as_ref().unwrap();
    let x_sorted = grid.x_sorted();
    let y_sorted = grid.y_sorted();
    let z_sorted = grid.z_sorted();
    
    println!("Rust x_sorted:");
    for ti in 0..3 {
        println!("  T{}: {:.2}", ti, x_sorted[[0, ti, 0, 0]]);
    }
    println!("\nRust y_sorted:");
    for ti in 0..3 {
        println!("  T{}: {:.2}", ti, y_sorted[[0, ti, 0, 0]]);
    }
    println!("\nRust z_sorted:");
    for ti in 0..3 {
        println!("  T{}: {:.2}", ti, z_sorted[[0, ti, 0, 0]]);
    }

    // Get velocities
    let velocities = &model.flow_field.u_sorted;
    let u_initial = &model.flow_field.u_initial_sorted;
    println!("\n\nRust FLORIS - Velocity field:");
    println!("  Shape: {:?}", velocities.shape());
    
    println!("\n  u_initial at T0:");
    for iy in 0..3 {
        for iz in 0..3 {
            print!("  {:.3}", u_initial[[0, 0, iy, iz]]);
        }
        println!();
    }
    
    println!("\n  u_initial at T1:");
    for iy in 0..3 {
        for iz in 0..3 {
            print!("  {:.3}", u_initial[[0, 1, iy, iz]]);
        }
        println!();
    }
    
    println!("\n  u_initial at T2:");
    for iy in 0..3 {
        for iz in 0..3 {
            print!("  {:.3}", u_initial[[0, 2, iy, iz]]);
        }
        println!();
    }

    println!("\n\n  Velocity at T0 (upstream):");
    for iy in 0..3 {
        for iz in 0..3 {
            print!("  {:.3}", velocities[[0, 0, iy, iz]]);
        }
        println!();
    }
    
    println!("\n  Velocity at T1 (middle):");
    for iy in 0..3 {
        for iz in 0..3 {
            print!("  {:.3}", velocities[[0, 1, iy, iz]]);
        }
        println!();
    }
    
    println!("\n  Velocity at T2 (downstream):");
    for iy in 0..3 {
        for iz in 0..3 {
            print!("  {:.3}", velocities[[0, 2, iy, iz]]);
        }
        println!();
    }

    // Get powers
    let powers = model.get_turbine_powers();
    println!("\n\nRust FLORIS - Power:");
    println!("  T0: {:.1} W", powers[[0, 0]]);
    println!("  T1: {:.1} W", powers[[0, 1]]);
    println!("  T2: {:.1} W", powers[[0, 2]]);
    println!("  Total: {:.1} W", powers.row(0).sum());

    Ok(())
}
