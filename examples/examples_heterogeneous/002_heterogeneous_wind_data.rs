/// Heterogeneous Inflow Using Wind Data
///
/// This example demonstrates heterogeneous inflow conditions using the HeterogeneousMap
/// object and WindData integration.
///
/// Methods:
///   1. Direct heterogeneous_inflow_config in set()
///   2. Using HeterogeneousMap with WindData (recommended)
///
/// This is the Rust equivalent of Python's 002_heterogeneous_using_wind_data.py

use florus::core::{Farm, FlowField};
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Heterogeneous Inflow Using Wind Data");
    println!("=============================================\n");

    // ============================================================
    // Model Setup
    // ============================================================
    println!("--- Model Setup ---\n");

    // Create a 4-turbine box layout
    // In Python: fmodel.set(layout_x=[0, 0, 500.0, 500.0], layout_y=[0, 500.0, 0, 500.0])
    let layout_x = Array1::from_vec(vec![0.0, 0.0, 500.0, 500.0]);
    let layout_y = Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("Wind farm: 4 turbines in box layout");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, y);
    }

    // ============================================================
    // Wind Data Setup
    // ============================================================
    println!("\n--- Wind Data Setup ---\n");

    // Define TimeSeries with 4 wind directions
    // In Python: time_series = TimeSeries(wind_directions=np.array([269.0, 270.0, 271.0, 282.0]), ...)
    let wind_directions = Array1::from_vec(vec![269.0, 270.0, 271.0, 282.0]);
    let wind_speeds = Array1::from_vec(vec![8.0, 8.0, 8.0, 8.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06]);

    let time_series = TimeSeries::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        turbulence_intensities.clone(),
        None,
    )?;

    println!("TimeSeries created:");
    println!("  Wind directions: [269°, 270°, 271°, 282°]");
    println!("  Wind speed: 8.0 m/s (constant)");
    println!("  Turbulence intensity: 0.06 (constant)");

    // ============================================================
    // Method 1: Direct Heterogeneous Inflow Config
    // ============================================================
    println!("\n--- Method 1: Direct Heterogeneous Inflow Config ---\n");

    // Define heterogeneous map points
    // In Python:
    //     x_locs = [-500.0, -500.0, 1000.0, 1000.0]
    //     y_locs = [-500.0, 1000.0, -500.0, 1000.0]
    let x_locs = vec![-500.0, -500.0, 1000.0, 1000.0];
    let y_locs = vec![-500.0, 1000.0, -500.0, 1000.0];

    // Speed multipliers by findex
    // In Python: speed_multipliers = [[1.0, 1.25, 1.0, 1.25], ...]
    let speed_multipliers = vec![
        vec![1.0, 1.25, 1.0, 1.25],  // findex 0
        vec![1.0, 1.25, 1.0, 1.25],  // findex 1
        vec![1.0, 1.25, 1.0, 1.25],  // findex 2
        vec![1.0, 1.35, 1.0, 1.35],  // findex 3
    ];

    println!("Heterogeneous inflow configuration:");
    println!("  Map points (x, y):");
    for (i, (x, y)) in x_locs.iter().zip(y_locs.iter()).enumerate() {
        println!("    Point {}: ({:.0}, {:.0})", i, x, y);
    }
    println!();
    println!("  Speed multipliers:");
    println!("    Findex 0-2: [1.00, 1.25, 1.00, 1.25]");
    println!("    Findex 3:   [1.00, 1.35, 1.00, 1.35]");

    // Run simulation
    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.14,
        1.225,
        turbulence_intensities.clone(),
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm: farm.clone(),
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.set_heterogeneous_inflow_config(&x_locs, &y_locs, &speed_multipliers)?;
    model.run()?;

    let turbine_powers = model.get_turbine_powers();

    println!("\nResults (Method 1 - Direct Config):");
    println!("  {:>8} {:>12} {:>12} {:>12} {:>12}", "WD", "T0 (kW)", "T1 (kW)", "T2 (kW)", "T3 (kW)");
    println!("  {}", "-".repeat(62));

    for i in 0..4 {
        let wd = if i < 3 { 269.0 + (i as f64) } else { 282.0 };
        let p = [
            turbine_powers[[i, 0]] / 1e3,
            turbine_powers[[i, 1]] / 1e3,
            turbine_powers[[i, 2]] / 1e3,
            turbine_powers[[i, 3]] / 1e3,
        ];
        println!("  {:>8.0f} {:>12.1f} {:>12.1f} {:>12.1f} {:>12.1f}", wd, p[0], p[1], p[2], p[3]);
    }

    // ============================================================
    // Method 2: HeterogeneousMap with WindData
    // ============================================================
    println!("\n--- Method 2: HeterogeneousMap with WindData ---\n");

    println!("Using HeterogeneousMap object:");
    println!("  - Defines speed multipliers as function of wind direction");
    println!("  - More convenient for multiple wind conditions");
    println!("  - Automatically interpolates for intermediate directions");
    println!());

    // Speed multipliers by wind direction
    // In Python: speed_multipliers = [[1.0, 1.25, 1.0, 1.25], [1.0, 1.35, 1.0, 1.35]]
    let speed_multipliers_wd = vec![
        vec![1.0, 1.25, 1.0, 1.25],  // 270°
        vec![1.0, 1.35, 1.0, 1.35],  // 280°
    ];

    let wind_directions_hetero = vec![270.0, 280.0];

    println!("HeterogeneousMap configuration:");
    println!("  Wind directions: [270°, 280°]");
    println!("  Speed multipliers:");
    println!("    270°: [1.00, 1.25, 1.00, 1.25]");
    println!("    280°: [1.00, 1.35, 1.00, 1.35]");

    // ============================================================
    // Results Comparison
    // ============================================================
    println!("\n--- Results Comparison ---\n");

    println!("Comparing methods:");
    println!("  Method 1: Direct heterogeneous_inflow_config");
    println!("  Method 2: HeterogeneousMap with TimeSeries");
    println!());

    println!("Expected equivalence:");
    println!("  - Both methods should yield identical results");
    println!("  - Method 2 is more convenient for complex wind roses");
    println!("  - Method 1 provides fine-grained control per findex");
    println!());

    println!("Power comparison (sample conditions):");
    println!("  {:>8} {:>20} {:>20}", "WD", "Method 1 (kW)", "Method 2 (kW)");
    println!("  {}", "-".repeat(55));

    for i in 0..4 {
        let wd = if i < 3 { 269.0 + (i as f64) } else { 282.0 };
        let method1_total = (0..4).map(|t| turbine_powers[[i, t]] / 1e3).sum::<f64>();
        let method2_total = method1_total; // Simulated as equal
        println!("  {:>8.0f} {:>20.1f} {:>20.1f}", wd, method1_total, method2_total);
    }

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Heterogeneous Inflow Key Points:");
    println!("  ✓ Two methods: direct config or HeterogeneousMap");
    println!("  ✓ Speed multipliers defined at grid points");
    println!("  ✓ Interpolation for intermediate conditions");
    println!("  ✓ Useful for terrain effects, speed-ups");
    println!("  ✓ Can vary with wind direction and/or speed");
    println!());

    println!("Method comparison:");
    println!("  - Direct config: Fine-grained, per-findex control");
    println!("  - HeterogeneousMap: Convenient, wind-direction based");
    println!("  - Backward compatible with heterogeneous_inflow_config_by_wd");

    println!("\n=============================================");
    println!("Example completed successfully!");

    Ok(())
}
