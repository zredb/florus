/// Example 10: Compare Farm Power with Neighboring Farm
///
/// This example demonstrates how to compare farm power production with and without
/// neighboring turbines. In Python FLORIS, this is done using turbine_weights to
/// exclude certain turbines from power calculations while still accounting for their
/// wake effects.
///
/// This Rust version demonstrates the concept by comparing a base farm with an
/// expanded farm that includes neighbors.

use florus::{Array1, FlorisModel};
use florus::core::turbines::TurbineLibrary;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    // Instantiate FLORIS using the GCH model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    println!("========== Farm Power Comparison with Neighbor ==========\n");

    // Define a 4 turbine farm (2x2 grid)
    let d = 126.0;
    println!("Step 1: Base farm with 4 turbines (2x2 grid)");
    let layout_x_base = vec![0.0, 6.0 * d, 0.0, 6.0 * d];
    let layout_y_base = vec![0.0, 0.0, 3.0 * d, 3.0 * d];
    
    fmodel.set_layout(
        &Array1::from_vec(layout_x_base.clone()),
        &Array1::from_vec(layout_y_base.clone()),
    )?;
    
    println!("  Layout:");
    for i in 0..4 {
        println!("    T{}: ({:.0}, {:.0}) m", i + 1, layout_x_base[i], layout_y_base[i]);
    }

    // Define wind conditions - sweep wind directions
    let wd_array: Vec<f64> = (0..360).step_by(4).map(|x| x as f64).collect();
    let n_wd = wd_array.len();
    let ws_array = vec![8.0; n_wd];
    let turbulence_intensities = vec![0.06; n_wd];
    
    println!("\n  Wind conditions: {} directions (0° to 356°, step 4°)", n_wd);
    println!("  Wind speed: 8.0 m/s (constant)");

    fmodel.set_wind_conditions(
        Array1::from_vec(ws_array.clone()),
        Array1::from_vec(wd_array.clone()),
        Array1::from_vec(turbulence_intensities.clone()),
    )?;

    // Calculate base case
    println!("\nRunning base case...");
    fmodel.run()?;
    let farm_power_base = fmodel.get_farm_power();
    
    println!("Base farm power at sample wind directions:");
    for &wd in &[0.0, 90.0, 180.0, 270.0] {
        let idx = (wd / 4.0) as usize;
        if idx < n_wd {
            println!(
                "  WD={:.0}°: Farm Power={:.0} kW",
                wd,
                farm_power_base[idx] / 1000.0
            );
        }
    }

    // Add a neighbor to the east (another 4 turbines)
    println!("\nStep 2: Adding neighboring farm to the east (4 more turbines)");
    let layout_x_neighbor = vec![
        0.0, 6.0 * d, 0.0, 6.0 * d,  // Original 4 turbines
        12.0 * d, 15.0 * d, 12.0 * d, 15.0 * d,  // 4 neighbor turbines to the east
    ];
    let layout_y_neighbor = vec![
        0.0, 0.0, 3.0 * d, 3.0 * d,  // Original y positions
        0.0, 0.0, 3.0 * d, 3.0 * d,  // Neighbor y positions (same)
    ];
    
    // Create a new model for the neighbor case to avoid state issues
    let mut fmodel_neighbor = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_neighbor.set_layout(
        &Array1::from_vec(layout_x_neighbor.clone()),
        &Array1::from_vec(layout_y_neighbor.clone()),
    )?;
    
    println!("  Total turbines: {}", layout_x_neighbor.len());
    println!("  Layout:");
    for i in 0..layout_x_neighbor.len() {
        let farm = if i < 4 { "Base" } else { "Neighbor" };
        println!(
            "    T{} ({}): ({:.0}, {:.0}) m",
            i + 1,
            farm,
            layout_x_neighbor[i],
            layout_y_neighbor[i]
        );
    }

    // Set the same wind conditions
    fmodel_neighbor.set_wind_conditions(
        Array1::from_vec(ws_array),
        Array1::from_vec(wd_array.clone()),
        Array1::from_vec(turbulence_intensities),
    )?;

    // Calculate with neighbor
    println!("\nRunning with neighboring farm...");
    fmodel_neighbor.run()?;
    
    // Get turbine powers for all 8 turbines
    let turbine_powers_all = fmodel_neighbor.get_turbine_powers();
    
    // Sum only the base farm turbines (first 4)
    let mut farm_power_with_neighbor = Array1::zeros(n_wd);
    for fi in 0..n_wd {
        for ti in 0..4 {
            farm_power_with_neighbor[fi] += turbine_powers_all[[fi, ti]];
        }
    }

    println!("Farm power (base turbines only) at sample wind directions:");
    for &wd in &[0.0, 90.0, 180.0, 270.0] {
        let idx = (wd / 4.0) as usize;
        if idx < n_wd {
            println!(
                "  WD={:.0}°: Farm Power={:.0} kW",
                wd,
                farm_power_with_neighbor[idx] / 1000.0
            );
        }
    }

    // Calculate and display the difference
    println!("\n========== Power Difference Analysis ==========");
    println!("Wind Dir | Base Farm | With Neighbor | Difference | % Change");
    println!("---------|-----------|---------------|------------|--------");
    
    for &wd in &[0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0] {
        let idx = (wd / 4.0) as usize;
        if idx < n_wd {
            let base = farm_power_base[idx] / 1000.0;
            let with_nbr = farm_power_with_neighbor[idx] / 1000.0;
            let diff = with_nbr - base;
            let pct_change = if base > 0.0 { 100.0 * diff / base } else { 0.0 };
            
            println!(
                "{:>7.0}° | {:>8.0} kW | {:>12.0} kW | {:>+9.0} kW | {:>+6.1}%",
                wd, base, with_nbr, diff, pct_change
            );
        }
    }

    // Find maximum impact
    let mut max_impact_wd = 0.0;
    let mut max_impact_pct = 0.0;
    
    for idx in 0..n_wd {
        let base = farm_power_base[idx];
        let with_nbr = farm_power_with_neighbor[idx];
        let diff = with_nbr - base;
        let pct_change = if base > 0.0 { 100.0 * diff.abs() / base } else { 0.0 };
        
        if pct_change > max_impact_pct {
            max_impact_pct = pct_change;
            max_impact_wd = wd_array[idx];
        }
    }
    
    println!("\nMaximum impact:");
    println!("  Wind direction: {:.0}°", max_impact_wd);
    println!("  Power change: {:.1}%", max_impact_pct);
    println!("\nNote: Negative changes indicate power loss due to neighbor wakes.");

    println!("\nExample 10 completed successfully!");
    println!("Note: Full turbine_weights feature will be implemented in future versions.");

    Ok(())
}
