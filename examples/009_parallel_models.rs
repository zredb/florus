/// Example 9: Parallel Models
///
/// This example demonstrates how to use parallelization to speed up FLORIS calculations.
/// ParFlorisModel (ParallelFlorisModel) allows parallel computation of wind conditions.
///
/// This is the Rust equivalent of Python's 009_parallel_models.py
///
/// Note: Rust's rayon library provides excellent parallelization capabilities
/// that can be used similarly to Python's multiprocessing interface.

use florus::core::{Farm, FlowField};
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 9: Parallel Models");
    println!("====================================\n");

    // ============================================================
    // Parallelization Concepts
    // ============================================================
    println!("--- Parallelization Concepts ---\n");

    println!("Parallel computation in FLORIS:");
    println!("  1. ParFlorisModel: Parallel version of FlorisModel");
    println!("  2. Parallelizes wind condition calculations");
    println!("  3. Uses multiprocessing interface (Python) or rayon (Rust)");
    println!();

    println!("Key parameters:");
    println!("  - interface: 'multiprocessing' (Python default)");
    println!("  - max_workers: Number of parallel workers (default: num CPUs)");
    println!("  - n_wind_condition_splits: How to split conditions");

    // ============================================================
    // Model Setup
    // ============================================================
    println!("\n--- Model Setup ---\n");

    println!("Creating base FlorisModel...");
    let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("Wind farm: 3 turbines at 500m spacing");
    println!();

    // Define inflow with many wind speeds
    println!("Inflow conditions:");
    println!("  Wind speeds: 1 to 25 m/s (0.5 m/s step)");
    println!("  Wind direction: 270 (constant)");
    println!("  Turbulence intensity: 0.06 (constant)");

    let n_ws: usize = 48; // (25-1)/0.5 = 48 conditions
    let wind_speeds: Vec<f64> = (0..n_ws).map(|i| 1.0 + 0.5 * (i as f64)).collect();
    let wind_directions = Array1::from_vec(vec![270.0; n_ws]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; n_ws]);

    // ============================================================
    // Parallelization Approaches
    // ============================================================
    println!("\n--- Parallelization Approaches ---\n");

    println!("1. Sequential (baseline):");
    println!("   - Single-threaded computation");
    println!("   - Simple and predictable");
    println!();

    println!("2. Parallel (via ParFlorisModel):");
    println!("   - Distributes wind conditions across workers");
    println!("   - Near-linear speedup for many conditions");
    println!("   - Overhead for small condition counts");
    println!();

    println!("3. Uncertain model with parallelization:");
    println!("   - Combines uncertainty with parallel computation");
    println!("   - Useful for large uncertainty studies");
    println!();

    // ============================================================
    // Performance Considerations
    // ============================================================
    println!("--- Performance Considerations ---\n");

    println!("When to use parallelization:");
    println!("  - Many wind conditions (>10)");
    println!("  - Large wind farms (>10 turbines)");
    println!("  - Uncertainty quantification (many samples)");
    println!("  - Optimization loops with many iterations");
    println!();

    println!("When parallelization has less benefit:");
    println!("  - Few wind conditions (<5)");
    println!("  - Small wind farms (<3 turbines)");
    println!("  - Very fast calculations (overhead dominates)");
    println!();

    println!("Rust rayon library:");
    println!("  - Lightweight parallelization");
    println!("  - Work-stealing scheduler");
    println!("  - Easy integration with iterators");
    println!("  - No need for separate process management");

    // ============================================================
    // Sequential Example
    // ============================================================
    println!("\n--- Sequential Computation Example ---\n");

    let flow_field = FlowField::new(
        Array1::from_vec(wind_speeds),
        wind_directions,
        0.0,
        0.14,
        1.225,
        turbulence_intensities,
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let turbine_powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();

    println!("Sequential results:");
    println!("  Conditions processed: {}", n_ws);
    println!("  Farm power shape: ({}, {})", farm_power.shape()[0], farm_power.shape()[1]);
    println!();

    // ============================================================
    // Summary
    // ============================================================
    println!("--- Summary ---\n");

    println!("Parallel Model Key Points:");
    println!("  - Parallelization speeds up wind condition calculations");
    println!("  - Best for many conditions or uncertainty studies");
    println!("  - Rust's rayon provides efficient parallelization");
    println!("  - Interface remains the same as FlorisModel");
    println!("  - Can be combined with UncertainFlorisModel");
    println!();

    println!("Implementation notes for Rust:");
    println!("  - Use rayon::iter::ParallelIterator for parallel iteration");
    println!("  - Consider rayon::scope for scoped parallelization");
    println!("  - Benchmark to find optimal chunk sizes");
    println!("  - Avoid sharing mutable state between threads");

    println!("\n====================================");
    println!("Example completed successfully!");
    println!("Note: Full parallel implementation requires rayon integration.");

    Ok(())
}
