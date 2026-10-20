//! Example: Optimize yaw and compare AEP
//!
//! This example demonstrates how to perform a yaw optimization and evaluate
//! the performance over a full wind rose.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/004_optimize_yaw_aep.py.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::{YawOptimization, YawOptimizationConfig};
use florus::types::{Array1, Array2, Float};
use florus::wind_data::{TimeSeries, WindRose};
use ndarray::Array;
use std::collections::BTreeMap;
use std::io::BufRead;
use std::time::Instant;

fn main() -> anyhow::Result<()> {
    // Load the wind rose from csv (same as 006_get_farm_aep)
    let wind_rose = read_wind_rose_csv("examples/inputs/wind_rose.csv", 0.06)?;

    // Load FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // 2x2 grid of turbines, 5D spacing
    let d = fmodel.core().farm.rotor_diameters_sorted[[0, 0]];
    let spacing = 5.0 * d;
    let layout_x = Array1::from_vec(vec![0.0, 0.0, spacing, spacing]);
    let layout_y = Array1::from_vec(vec![0.0, spacing, 0.0, spacing]);
    fmodel.set_layout(&layout_x, &layout_y)?;

    let n_turbines = layout_x.len();

    // Optimize yaw at 8 m/s for all wind directions in the rose
    let time_series = TimeSeries::new(
        wind_rose.wind_directions.clone(),
        Array1::from_elem(wind_rose.wind_directions.len(), 8.0),
        Array1::from_elem(wind_rose.wind_directions.len(), 0.06),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    let start = Instant::now();
    let config = YawOptimizationConfig {
        minimum_yaw_angle: 0.0,
        maximum_yaw_angle: 20.0,
        ..Default::default()
    };
    let mut yaw_opt = YawOptimizationSR::new();
    let df_opt = yaw_opt.optimize(&mut fmodel, Some(config))?;
    let t_tot = start.elapsed().as_secs_f64();
    println!("Optimization finished in {t_tot:.2} seconds.");

    // Calculate the AEP in the baseline case
    fmodel.set_wind_data(&wind_rose)?;
    fmodel.run()?;
    let aep_baseline = fmodel.get_farm_aep(&wind_rose, 8760.0)?;

    // Interpolate the optimal yaw angles for the full wind rose using the
    // rule of thumb: full steering between 6 and 12 m/s, linear ramp
    // outside that range, no steering below 4 or above 14 m/s.
    let n_findex = fmodel.flow_field().n_findex;
    let wind_speeds = fmodel.flow_field().wind_speeds.clone();
    let wind_directions = fmodel.flow_field().wind_directions.clone();

    let mut yaw_angles_wind_rose = Array2::zeros((n_findex, n_turbines));
    let wind_speed_low_no_steer = 4.0;
    let wind_speed_high_no_steer = 14.0;
    let wind_speed_low_steer = 6.0;
    let wind_speed_high_steer = 12.0;

    for i in 0..n_findex {
        let ws = wind_speeds[i];
        let wd = wind_directions[i];

        // Find the optimal yaw at 8 m/s for this wind direction
        // (df_opt row with matching wind_direction)
        let opt_fi = (0..df_opt.yaw_angles.shape()[0])
            .find(|&k| (wind_rose.wind_directions[k] - wd).abs() < 1e-6)
            .unwrap_or(0);
        let yaw_opt_full: Vec<Float> = (0..n_turbines)
            .map(|t| df_opt.yaw_angles[[opt_fi, t]])
            .collect();

        let factor = if ws < wind_speed_low_no_steer || ws > wind_speed_high_no_steer {
            0.0
        } else if ws < wind_speed_low_steer {
            (wind_speed_low_steer - ws) / (wind_speed_low_steer - wind_speed_low_no_steer)
        } else if ws > wind_speed_high_steer {
            (wind_speed_high_no_steer - ws) / (wind_speed_high_no_steer - wind_speed_high_steer)
        } else {
            1.0
        };

        for t in 0..n_turbines {
            yaw_angles_wind_rose[[i, t]] = yaw_opt_full[t] * factor;
        }
    }

    // Apply the optimal yaw angles and get the AEP
    fmodel.set_yaw_angles(yaw_angles_wind_rose)?;
    fmodel.run()?;
    let aep_opt = fmodel.get_farm_aep(&wind_rose, 8760.0)?;
    let aep_uplift = 100.0 * (aep_opt / aep_baseline - 1.0);

    println!("Baseline AEP: {:.2} GWh.", aep_baseline / 1e9);
    println!("Optimal AEP: {:.2} GWh.", aep_opt / 1e9);
    println!("Relative AEP uplift by wake steering: {aep_uplift:.3} %.");

    Ok(())
}

/// Replicate Python `WindRose.read_csv_long(..., freq_col=..., ti_col_or_value=0.06)`.
fn read_wind_rose_csv(path: &str, ti_value: f64) -> anyhow::Result<WindRose> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);

    let mut rows: Vec<(f64, f64, f64)> = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue;
        }
        let mut cols = line.split(',');
        let ws: f64 = cols.next().unwrap_or_default().parse()?;
        let wd: f64 = cols.next().unwrap_or_default().parse()?;
        let freq: f64 = cols.next().unwrap_or_default().parse()?;
        rows.push((wd, ws, freq));
    }

    let freq_sum: f64 = rows.iter().map(|r| r.2).sum();
    for r in rows.iter_mut() {
        r.2 /= freq_sum;
    }

    let mut unique_wd: Vec<f64> = rows.iter().map(|r| r.0).collect();
    unique_wd.sort_by(|a, b| a.partial_cmp(b).unwrap());
    unique_wd.dedup();
    let mut unique_ws: Vec<f64> = rows.iter().map(|r| r.1).collect();
    unique_ws.sort_by(|a, b| a.partial_cmp(b).unwrap());
    unique_ws.dedup();

    let wd_index: BTreeMap<u64, usize> = unique_wd
        .iter()
        .enumerate()
        .map(|(i, &v)| (v.to_bits(), i))
        .collect();
    let ws_index: BTreeMap<u64, usize> = unique_ws
        .iter()
        .enumerate()
        .map(|(i, &v)| (v.to_bits(), i))
        .collect();

    let n_dir = unique_wd.len();
    let n_ws = unique_ws.len();
    let mut freq_table = Array2::zeros((n_dir, n_ws));
    for &(wd, ws, freq) in &rows {
        let i = *wd_index.get(&wd.to_bits()).unwrap();
        let j = *ws_index.get(&ws.to_bits()).unwrap();
        freq_table[[i, j]] += freq;
    }

    Ok(WindRose::new(
        Array1::from_vec(unique_wd),
        Array1::from_vec(unique_ws),
        Array2::from_elem((n_dir, n_ws), ti_value),
        Some(freq_table),
        None,
        false,
        None,
        None,
    )?)
}
