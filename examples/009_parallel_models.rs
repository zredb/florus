/// Example 9: Parallel Models (Conceptual)
///
/// This example demonstrates the concept of parallel FLORIS models.
/// Note: Full ParFlorisModel implementation is not yet available in Rust.
/// In Python, ParFlorisModel uses multiprocessing to parallelize FLORIS calculations.
///
/// In Rust, parallelization can be achieved through:
/// - Rayon for data parallelism
/// - tokio for async parallelism
/// - std::thread for manual threading
///
/// This example shows how the same FlorisModel can be used in different contexts.

use florus::{Array1, FlorisModel};
use florus::core::turbines::TurbineLibrary;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    println!("========== Parallel Models Concept ==========\n");

    // Instantiate the FlorisModel
    let fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    println!("Created base FlorisModel from gch.yaml");

    // In Rust, we can create multiple independent models
    println!("\nCreating multiple independent models...");
    let mut fmodel1 = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut fmodel2 = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut fmodel3 = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    println!("Created 3 independent FlorisModel instances");

    // Set up a simple wind speed sweep
    let wind_speeds: Vec<f64> = (5..25).map(|x| x as f64).collect();
    let n_ws = wind_speeds.len();
    let wind_directions = vec![270.0; n_ws];
    let turbulence_intensities = vec![0.06; n_ws];

    // Configure all models with the same conditions
    println!("\nConfiguring models with wind speed sweep ({} conditions)...", n_ws);
    
    for fmodel in [&mut fmodel1, &mut fmodel2, &mut fmodel3] {
        fmodel.set_wind_conditions(
            Array1::from_vec(wind_speeds.clone()),
            Array1::from_vec(wind_directions.clone()),
            Array1::from_vec(turbulence_intensities.clone()),
        )?;
    }

    // Run models sequentially (in production, these could run in parallel)
    println!("\nRunning models...");
    
    println!("  Running model 1...");
    fmodel1.run()?;
    let powers1 = fmodel1.get_farm_power();
    
    println!("  Running model 2...");
    fmodel2.run()?;
    let powers2 = fmodel2.get_farm_power();
    
    println!("  Running model 3...");
    fmodel3.run()?;
    let powers3 = fmodel3.get_farm_power();

    // Verify results are identical
    println!("\nVerifying results consistency...");
    let all_close = powers1.iter().zip(powers2.iter()).zip(powers3.iter())
        .all(|((p1, p2), p3)| {
            (p1 - p2).abs() < 1e-6 && (p1 - p3).abs() < 1e-6
        });
    
    println!("  All models produce identical results: {}", all_close);

    // Show sample results
    println!("\nSample farm power results:");
    for i in [0, 5, 10, 15, 19] {
        if i < n_ws {
            println!(
                "  WS={:.0} m/s: Model1={:.0} kW, Model2={:.0} kW, Model3={:.0} kW",
                wind_speeds[i],
                powers1[i] / 1000.0,
                powers2[i] / 1000.0,
                powers3[i] / 1000.0
            );
        }
    }

    println!("\n========== Parallel Execution Notes ==========");
    println!("In a production environment, you could:");
    println!("  1. Use Rayon's par_iter() for data parallelism");
    println!("  2. Use tokio::spawn for async parallel execution");
    println!("  3. Use std::thread::spawn for manual threading");
    println!("\nExample parallel pattern with Rayon:");
    println!("  use rayon::prelude::*;");
    println!("  let results: Vec<_> = models.par_iter_mut()");
    println!("      .map(|m| {{ m.run()?; m.get_farm_power() }})");
    println!("      .collect();");

    println!("\nExample 9 completed successfully!");
    println!("Note: Full ParFlorisModel will be implemented in future versions.");

    Ok(())
}
