/// Disabling Turbines Example
///
/// This example demonstrates the ability of FLORIS to shut down some turbines
/// during a simulation using the "mixed" operation model.
///
/// This is the Rust equivalent of Python's 002_disable_turbines.py

use florus::core::{Farm, FlowField};
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Disabling Turbines");
    println!("============================\n");

    // ============================================================
    // Model Setup
    // ============================================================
    println!("--- Model Setup ---\n");

    println!("Initializing FLORIS with 'mixed' operation model...");
    println!("The 'mixed' model enables turbine disable and power derating.\n");

    // Create a 3-turbine aligned layout
    // In Python: layout = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
    let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("Wind farm: 3 turbines aligned at 500m spacing");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, y);
    }

    // ============================================================
    // Wind Conditions
    // ============================================================
    println!("\n--- Wind Conditions ---\n");

    // Two identical wind conditions (n_findex = 2)
    let wind_directions = Array1::from_vec(vec![270.0, 270.0]);
    let wind_speeds = Array1::from_vec(vec![8.0, 8.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.06]);

    println!("Condition 1: All turbines operating");
    println!("  Wind direction: 270°");
    println!("  Wind speed: 8.0 m/s");
    println!("  Turbulence intensity: 0.06");
    println!();
    println!("Condition 2: T0 and T1 disabled");
    println!("  Same wind conditions");
    println!("  Turbines 0 and 1 shut down");

    // ============================================================
    // Turbine Disable Configuration
    // ============================================================
    println!("\n--- Turbine Disable Configuration ---\n");

    println!("Disable configuration (n_findex x n_turbines):");
    println!("  [[false, false, false],   // Condition 1: All on");
    println!("   [true,  true,  false]]    // Condition 2: T0, T1 disabled");
    println!());

    // In Python: disable_turbines = np.array([[False, False, False], [True, True, False]])
    let disable_turbines = Array2::from_vec(vec![
        vec![false, false, false],
        vec![true, true, false],
    ]);

    // ============================================================
    // Simulation
    // ============================================================
    println!("--- Simulation ---\n");

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
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.set_disable_turbines(disable_turbines.clone())?;
    model.run()?;

    // ============================================================
    // Results
    // ============================================================
    println!("--- Results ---\n");

    let turbine_powers = model.get_turbine_powers();
    let effective_wind_speeds = model.get_turbine_average_velocities();

    println!("Turbine powers (kW):");
    println!("  {:>12} {:>12} {:>12}", "T0", "T1", "T2");
    println!("  {}", "-".repeat(42));

    for findex in 0..2 {
        let status = if findex == 0 { "All on" } else { "T0,T1 disabled" };
        let p0 = turbine_powers[[findex, 0]] / 1e3;
        let p1 = turbine_powers[[findex, 1]] / 1e3;
        let p2 = turbine_powers[[findex, 2]] / 1e3;
        println!("  Condition {}: {:>10.1f} {:>10.1f} {:>10.1f}  ({})", findex + 1, p0, p1, p2, status);
    }

    println!();
    println!("Effective wind speeds (m/s):");
    println!("  {:>12} {:>12} {:>12}", "T0", "T1", "T2");
    println!("  {}", "-".repeat(42));

    for findex in 0..2 {
        let status = if findex == 0 { "All on" } else { "T0,T1 disabled" };
        let w0 = effective_wind_speeds[[findex, 0]];
        let w1 = effective_wind_speeds[[findex, 1]];
        let w2 = effective_wind_speeds[[findex, 2]];
        println!("  Condition {}: {:>10.1f} {:>10.1f} {:>10.1f}  ({})", findex + 1, w0, w1, w2, status);
    }

    // ============================================================
    // Analysis
    // ============================================================
    println!("\n--- Analysis ---\n");

    println!("Effect of disabling upstream turbines:");
    println!("  1. Disabled turbines produce zero power");
    println!("  2. Downstream turbines see higher wind speeds (no upstream wake)");
    println!("  3. T2 effective wind speed increases when T0/T1 disabled");
    println!());

    let ws_increase = effective_wind_speeds[[1, 2]] - effective_wind_speeds[[0, 2]];
    let power_increase = (turbine_powers[[1, 2]] - turbine_powers[[0, 2]]) / 1e3;

    println!("Quantified effects on T2:");
    println!("  Wind speed increase: {:.1f} m/s", ws_increase);
    println!("  Power increase: {:.1f} kW", power_increase);
    println!("  Percentage increase: {:.1f}%", power_increase / turbine_powers[[0, 2]] * 1e3 * 100.0);

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Disabling Turbines Key Points:");
    println!("  ✓ Use 'mixed' operation model for disable capability");
    println!("  ✓ disable_turbines: 2D array (n_findex x n_turbines)");
    println!("  ✓ Disabled turbines produce zero power");
    println!("  ✓ Downstream turbines benefit from reduced wake");
    println!("  ✓ Useful for maintenance scheduling simulations");
    println!("  ✓ Can model partial farm outages");

    println!("\n============================");
    println!("Example completed successfully!");

    Ok(())
}
