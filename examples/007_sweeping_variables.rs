/// Example 7: Sweeping Variables
///
/// Demonstrate methods for sweeping across variables. Wind directions, wind speeds,
/// turbulence intensities, as well as control inputs are passed as arrays and can be
/// swept and run in one call to run().
///
/// This example demonstrates sweeping:
///     1) Wind speeds
///     2) Wind directions
///     3) Turbulence intensities
///     4) Yaw angles

use florus::{Array1, Array2, FlorisModel};
use florus::core::turbines::TurbineLibrary;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set to a 2 turbine layout
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    println!("Layout: 2 turbines at [0, {:.0}]m", 5.0 * d);
    println!("  n_turbines: {}", fmodel.n_turbines());

    //////////////////////////////////////////////////
    // Sweep wind speeds
    //////////////////////////////////////////////////

    println!("\n========== Sweep Wind Speeds ==========");
    let wind_speeds: Vec<f64> = (50..100).map(|x| x as f64 / 10.0).collect(); // 5.0 to 9.9 m/s
    let n_ws = wind_speeds.len();
    let wind_directions = vec![270.0; n_ws];
    let turbulence_intensities = vec![0.06; n_ws];

    fmodel.set_wind_conditions(
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(wind_directions),
        Array1::from_vec(turbulence_intensities),
    )?;

    fmodel.run()?;
    let turbine_powers = fmodel.get_turbine_powers();

    println!("Wind speed sweep results (first 5):");
    for i in 0..5.min(n_ws) {
        println!(
            "  WS={:.1} m/s: T1={:.0} kW, T2={:.0} kW",
            wind_speeds[i],
            turbine_powers[[i, 0]] / 1000.0,
            turbine_powers[[i, 1]] / 1000.0
        );
    }

    println!("\nExample 7 completed successfully!");

    //////////////////////////////////////////////////
    // Sweep wind directions (separate model)
    //////////////////////////////////////////////////

    println!("\n========== Sweep Wind Directions ==========");
    let mut fmodel_wd = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_wd.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    
    let wind_directions: Vec<f64> = (250..290).map(|x| x as f64).collect();
    let n_wd = wind_directions.len();
    let wind_speeds = vec![8.0; n_wd];
    let turbulence_intensities = vec![0.06; n_wd];

    fmodel_wd.set_wind_conditions(
        Array1::from_vec(wind_speeds),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(turbulence_intensities),
    )?;

    fmodel_wd.run()?;
    let turbine_powers = fmodel_wd.get_turbine_powers();

    println!("Wind direction sweep results (sample):");
    for &wd in &[250.0, 260.0, 270.0, 280.0, 289.0] {
        let idx = (wd - 250.0) as usize;
        if idx < n_wd {
            println!(
                "  WD={:.0}°: T1={:.0} kW, T2={:.0} kW",
                wd,
                turbine_powers[[idx, 0]] / 1000.0,
                turbine_powers[[idx, 1]] / 1000.0
            );
        }
    }

    //////////////////////////////////////////////////
    // Sweep turbulence intensities (separate model)
    //////////////////////////////////////////////////

    println!("\n========== Sweep Turbulence Intensities ==========");
    let mut fmodel_ti = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_ti.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    
    let turbulence_intensities: Vec<f64> = (3..20).map(|x| x as f64 / 100.0).collect(); // 0.03 to 0.19
    let n_ti = turbulence_intensities.len();
    let wind_speeds = vec![8.0; n_ti];
    let wind_directions = vec![270.0; n_ti];

    fmodel_ti.set_wind_conditions(
        Array1::from_vec(wind_speeds),
        Array1::from_vec(wind_directions),
        Array1::from_vec(turbulence_intensities.clone()),
    )?;

    fmodel_ti.run()?;
    let turbine_powers = fmodel_ti.get_turbine_powers();

    println!("Turbulence intensity sweep results (sample):");
    for &ti in &[0.03, 0.06, 0.10, 0.15, 0.19] {
        let idx = ((ti - 0.03) * 100.0) as usize;
        if idx < n_ti {
            println!(
                "  TI={:.2}: T1={:.0} kW, T2={:.0} kW",
                ti,
                turbine_powers[[idx, 0]] / 1000.0,
                turbine_powers[[idx, 1]] / 1000.0
            );
        }
    }

    //////////////////////////////////////////////////
    // Sweep the upstream yaw angle (separate model)
    //////////////////////////////////////////////////

    println!("\n========== Sweep Upstream Yaw Angle ==========");
    let mut fmodel_yaw = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_yaw.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    
    let n_yaw = 21;
    let wind_directions = vec![270.0; n_yaw];
    let wind_speeds = vec![8.0; n_yaw];
    let turbulence_intensities = vec![0.06; n_yaw];

    fmodel_yaw.set_wind_conditions(
        Array1::from_vec(wind_speeds),
        Array1::from_vec(wind_directions),
        Array1::from_vec(turbulence_intensities),
    )?;

    // Create yaw angles array: upstream turbine sweeps from -30 to 30 degrees
    let yaw_angles_upstream: Vec<f64> = (-30..=30).step_by(3).map(|x| x as f64).collect();
    let mut yaw_angles = Array2::zeros((n_yaw, 2));
    for (i, &yaw) in yaw_angles_upstream.iter().enumerate() {
        yaw_angles[[i, 0]] = yaw;
        yaw_angles[[i, 1]] = 0.0; // Downstream turbine has no yaw
    }

    fmodel_yaw.set_yaw_angles(yaw_angles)?;
    fmodel_yaw.run()?;
    let turbine_powers = fmodel_yaw.get_turbine_powers();

    println!("Yaw angle sweep results:");
    for (i, &yaw) in yaw_angles_upstream.iter().enumerate() {
        println!(
            "  Yaw={:.0}°: T1={:.0} kW, T2={:.0} kW",
            yaw,
            turbine_powers[[i, 0]] / 1000.0,
            turbine_powers[[i, 1]] / 1000.0
        );
    }

    Ok(())
}
