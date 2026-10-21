//! Example: The load optimization of a row of turbines
//!
//! This example optimizes the power setpoints of a row of three turbines so
//! that the net revenue (the farm revenue minus the variable operating cost)
//! is maximized, for a set of wind directions, ambient load turbulence
//! intensities and electricity prices.
//!
//! Scenario aligned 1:1 with
//! examples_load_optimization/002_row_opt_example.py.

use florus::core::turbines::{POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED};
use florus::floris_model::FlorisModel;
use florus::optimization::load_optimization::load_optimization::{
    compute_farm_revenue, compute_farm_voc, compute_net_revenue, optimize_power_setpoints,
};
use florus::types::{Array1, Array2};
use florus::utilities::linspace;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // The parameters
    let d = 126.0;
    let d_spacing = 7.0;
    let n_turbines = 3;
    let a = 4.0e-6; // Selected to demonstrate variation in derating selection

    let mut fmodel = FlorisModel::from_default()?;

    // A row of turbines
    fmodel.set_layout(
        &Array1::from_iter((0..n_turbines).map(|i| i as f64 * d * d_spacing)),
        &Array1::from_vec(vec![0.0; n_turbines]),
    )?;
    fmodel.set_operation_model("simple-derating")?;

    // The input conditions: two wind directions, two ambient turbulence
    // intensities and two electricity values
    let mut wind_directions: Vec<f64> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut ambient_lti: Vec<f64> = Vec::new();
    for w_i in [240.0, 270.0] {
        for t_i in [0.1, 0.25] {
            for v_i in [1.0e-5, 1.0e-6] {
                wind_directions.push(w_i);
                values.push(v_i);
                ambient_lti.push(t_i);
            }
        }
    }
    let n = wind_directions.len();

    let mut time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    time_series.values = Array1::from_vec(values.clone());

    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    // The initial power setpoints (no derating)
    let initial_power_setpoint = Array2::from_elem((n, n_turbines), 5.0e6);
    fmodel.set_operation(
        None,
        Some(initial_power_setpoint.clone()),
        None,
        None,
        None,
        None,
    )?;
    fmodel.run()?;

    let net_revenue_initial = compute_net_revenue(&fmodel, a, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;
    let farm_voc_initial = compute_farm_voc(&fmodel, a, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;
    let farm_revenue_initial = compute_farm_revenue(&fmodel)?;

    // The optimal derating levels (Python: np.linspace(POWER_SETPOINT_DEFAULT,
    // POWER_SETPOINT_DISABLED, 10))
    let power_setpoint_levels = linspace(POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED, 10);
    let (opt_power_setpoints, opt_net_revenue) = optimize_power_setpoints(
        &mut fmodel,
        a,
        &ambient_lti,
        0.3,
        10.0,
        1.0,
        1.0,
        Some(&initial_power_setpoint),
        &power_setpoint_levels,
    )?;

    // The final values
    fmodel.set_operation(
        None,
        Some(opt_power_setpoints.clone()),
        None,
        None,
        None,
        None,
    )?;
    fmodel.run()?;
    let net_revenue_opt = compute_net_revenue(&fmodel, a, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;
    let farm_voc_opt = compute_farm_voc(&fmodel, a, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;
    let farm_revenue_opt = compute_farm_revenue(&fmodel)?;

    println!("Wind direction [deg], ambient LTI, value [USD/MWh], initial farm revenue, initial farm VOC, initial net revenue");
    for i in 0..n {
        println!(
            "{:.0}, {:.4}, {:.8}, {:.3}, {:.3}, {:.3}",
            wind_directions[i],
            ambient_lti[i],
            values[i],
            farm_revenue_initial[i],
            farm_voc_initial[i],
            net_revenue_initial[i]
        );
    }

    println!("Optimized power setpoints [W]:");
    for i in 0..n {
        let row: Vec<String> = (0..n_turbines)
            .map(|t| format!("{:.3}", opt_power_setpoints[[i, t]]))
            .collect();
        println!("  condition {i}: {}", row.join(", "));
    }

    println!(
        "Optimized farm revenue [USD], optimized farm VOC [USD], optimized net revenue [USD]:"
    );
    for i in 0..n {
        println!(
            "  condition {i}: {:.3}, {:.3}, {:.3}",
            farm_revenue_opt[i], farm_voc_opt[i], net_revenue_opt[i]
        );
    }

    println!("Optimized net revenue per condition [USD]:");
    for i in 0..n {
        println!(
            "  condition {i}: {:.3} (initial {:.3})",
            opt_net_revenue[i], net_revenue_initial[i]
        );
    }

    Ok(())
}
