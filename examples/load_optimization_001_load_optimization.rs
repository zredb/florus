//! Example: The load turbulence intensity and the variable operating cost
//!
//! This example computes the load turbulence intensity (LTI) and the variable
//! operating cost (VOC) of a two-turbine farm for a sweep of the wind
//! direction, a sweep of the power setpoint of the upstream turbine and a sweep
//! of the ambient load turbulence intensity.
//!
//! Scenario aligned 1:1 with
//! examples_load_optimization/001_lti_and_voc.py.

use florus::floris_model::FlorisModel;
use florus::optimization::load_optimization::load_optimization::{
    compute_lti, compute_turbine_voc,
};
use florus::types::{Array1, Array2};
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_default()?;
    fmodel.set_operation_model("simple-derating")?;

    // The wind direction sweep
    let wind_directions: Vec<f64> = (0..360).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let ambient_lti = vec![0.1; n];
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;

    let d = fmodel.core().farm.rotor_diameters[0];
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, d * 7.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    let load_ti = compute_lti(&fmodel, &ambient_lti, 0.3, 10.0)?;
    let voc = compute_turbine_voc(&fmodel, 2.0e-5, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;

    println!("Wind direction [deg], LTI t0, LTI t1, VOC t0 [USD], VOC t1 [USD]");
    for i in 0..n {
        println!(
            "{:.0}, {:.6}, {:.6}, {:.3}, {:.3}",
            wind_directions[i],
            load_ti[[i, 0]],
            load_ti[[i, 1]],
            voc[[i, 0]],
            voc[[i, 1]]
        );
    }

    // The power setpoint sweep
    let n = 50;
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    let setpoints: Vec<f64> = (0..n)
        .map(|i| 5.0e6 + (0.001 - 5.0e6) * i as f64 / (n - 1) as f64)
        .collect();
    let power_setpoints = Array2::from_shape_vec(
        (n, 2),
        setpoints
            .iter()
            .flat_map(|s| vec![*s, 1.0e12])
            .collect::<Vec<f64>>(),
    )?;
    let ambient_lti = vec![0.1; n];

    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, d * 7.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.set_operation(None, Some(power_setpoints), None, None, None, None)?;
    fmodel.run()?;

    let load_ti = compute_lti(&fmodel, &ambient_lti, 0.3, 10.0)?;
    let voc = compute_turbine_voc(&fmodel, 2.0e-5, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;

    println!("Power setpoint [W], LTI t0, LTI t1, VOC t0 [USD], VOC t1 [USD]");
    for i in 0..n {
        println!(
            "{:.1}, {:.6}, {:.6}, {:.3}, {:.3}",
            setpoints[i],
            load_ti[[i, 0]],
            load_ti[[i, 1]],
            voc[[i, 0]],
            voc[[i, 1]]
        );
    }

    // The load ambient TI sweep
    let ambient_lti: Vec<f64> = (0..n)
        .map(|i| 0.05 + 0.2 * i as f64 / (n - 1) as f64)
        .collect();
    let power_setpoints = Array2::from_elem((n, 2), 1.0e12);
    fmodel.set_operation(None, Some(power_setpoints), None, None, None, None)?;
    fmodel.run()?;

    let load_ti = compute_lti(&fmodel, &ambient_lti, 0.3, 10.0)?;
    let voc = compute_turbine_voc(&fmodel, 2.0e-5, &ambient_lti, 0.3, 10.0, 1.0, 1.0)?;

    println!("Load ambient TI, LTI t0, LTI t1, VOC t0 [USD], VOC t1 [USD]");
    for i in 0..n {
        println!(
            "{:.4}, {:.6}, {:.6}, {:.3}, {:.3}",
            ambient_lti[i],
            load_ti[[i, 0]],
            load_ti[[i, 1]],
            voc[[i, 0]],
            voc[[i, 1]]
        );
    }

    Ok(())
}
