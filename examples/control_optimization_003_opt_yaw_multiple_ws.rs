//! Example: Optimize yaw for multiple wind directions and multiple wind speeds.
//!
//! This example demonstrates how to perform a yaw optimization for multiple
//! wind directions and multiple wind speeds using the WindRose object.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/003_opt_yaw_multiple_ws.py: a 3-turbine row
//! on gch.yaml, a uniform wind rose (120 directions x 16 speeds, TI = 0.06),
//! and the Serial-Refine optimizer with Python defaults (Ny_passes = [5, 4],
//! yaw in [0, 25] deg). Prints the progress lines and the pandas-style
//! results table exactly like Python.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::YawOptimization;
use florus::types::{Array1, Array2, Float};
use florus::wind_data::WindRose;
use ndarray::Array;

/// Python pandas scientific float formatting: mantissa with 6 decimals,
/// exponent with sign and at least two digits ("1.095593e+05").
fn py_sci(v: f64) -> String {
    let s = format!("{v:.6e}"); // e.g. "1.095593e5"
    if let Some(pos) = s.find('e') {
        let (mantissa, exp) = s.split_at(pos);
        let digits: i32 = exp[1..].parse().unwrap();
        format!(
            "{mantissa}e{}{:02}",
            if digits < 0 { '-' } else { '+' },
            digits.abs()
        )
    } else {
        s
    }
}

/// Print the results table in the exact pandas truncation format recorded
/// from Python (head 5, ellipsis row, tail 5), with the middle columns
/// (turbulence_intensity, yaw_angles_opt) collapsed into "...".
fn print_results_table(
    wd: &Array1,
    ws: &Array1,
    farm_power_opt: &Array1,
    farm_power_baseline: &Array1,
) {
    let n = wd.len();
    // Column widths measured from the Python output
    const IDX_W: usize = 4; // index column (left-aligned)
    const WD_W: usize = 14; // wind_direction header width
    const WS_W: usize = 10; // wind_speed header width
    const FPO_W: usize = 14; // farm_power_opt header width
    const FPB_W: usize = 19; // farm_power_baseline header width

    let header = format!(
        "{:<IDX_W$}  {:>WD_W$}  {:>WS_W$}  ...  {:>FPO_W$} {:>FPB_W$}",
        "", "wind_direction", "wind_speed", "farm_power_opt", "farm_power_baseline"
    );
    println!("{header}");

    let mut print_row = |idx: &str, fi: usize| {
        let line = format!(
            "{:<IDX_W$}  {:>WD_W$.1}  {:>WS_W$.1}  ...  {:>FPO_W$} {:>FPB_W$}",
            idx,
            wd[fi],
            ws[fi],
            py_sci(farm_power_opt[fi]),
            py_sci(farm_power_baseline[fi])
        );
        println!("{line}");
    };

    for fi in 0..5.min(n) {
        print_row(&fi.to_string(), fi);
    }
    if n > 10 {
        println!(
            "{:<IDX_W$}  {:>WD_W$}  {:>WS_W$}  ...  {:>FPO_W$} {:>FPB_W$}",
            "...", "...", "...", "...", "..."
        );
    }
    for fi in n.saturating_sub(5)..n {
        print_row(&fi.to_string(), fi);
    }

    println!();
    println!("[{n} rows x 6 columns]");
}

fn main() -> anyhow::Result<()> {
    // Load the default example floris object
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Define a WindRose object with uniform TI and frequency table
    let wind_directions = Array::range(0.0, 360.0, 3.0);
    let wind_speeds = Array::range(2.0, 18.0, 1.0);
    let ti_table = Array2::from_elem((wind_directions.len(), wind_speeds.len()), 0.06 as Float);
    let wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        ti_table,
        None, // uniform frequency table
        None,
        false,
        None,
        None,
    )?;

    // Reinitialize as a 3-turbine farm with the full wind rose
    let d = 126.0; // Rotor diameter for the NREL 5 MW
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.set_wind_data(&wind_rose)?;
    fmodel.run()?;

    // Initialize optimizer object and run optimization using the
    // Serial-Refine method (Python defaults: yaw in [0, 25] deg,
    // Ny_passes = [5, 4], verify_convergence = False).
    let mut yaw_opt = YawOptimizationSR::new();
    let result = yaw_opt.optimize(&mut fmodel, None)?;

    println!("Optimization results:");
    print_results_table(
        &fmodel.flow_field().wind_directions,
        &fmodel.flow_field().wind_speeds,
        &result.farm_powers_opt,
        &result.farm_powers_baseline,
    );

    Ok(())
}
