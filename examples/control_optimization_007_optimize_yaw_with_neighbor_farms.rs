//! Example: Optimize yaw with neighbor farm
//!
//! This example demonstrates how to optimize the yaw angles of a subset of
//! turbines in order to maximize the annual energy production (AEP) of a wind
//! farm. In this case, the wind farm is part of a larger collection of
//! turbines, some of which are part of a neighboring farm. The optimization
//! is performed in two ways: first by accounting for the wakes of the
//! neighboring farm (while not including those turbines in the optimization
//! objective), and second by removing the neighboring farm from the model
//! for the optimization.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/007_optimize_yaw_with_neighbor_farms.py.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::{YawOptimization, YawOptimizationConfig};
use florus::types::{Array1, Array2, Float};
use florus::wind_data::{TimeSeries, WindRose};
use ndarray::Array;
use std::collections::BTreeMap;
use std::io::BufRead;

fn main() -> anyhow::Result<()> {
    // Load the wind rose from csv
    let wind_rose = read_wind_rose_csv("examples/inputs/wind_rose.csv", 0.06)?;

    // Load FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Layout: the first 20 turbines are the farm of interest (4 rows x 5),
    // the remaining 9 belong to a neighboring farm.
    let x_raw: [f64; 29] = [
        0.0, 756.0, 1512.0, 2268.0, 3024.0, 0.0, 756.0, 1512.0, 2268.0, 3024.0, 0.0, 756.0, 1512.0,
        2268.0, 3024.0, 0.0, 756.0, 1512.0, 2268.0, 3024.0, 4500.0, 5264.0, 6028.0, 4878.0, 0.0,
        756.0, 1512.0, 2268.0, 3024.0,
    ];
    let y_raw: [f64; 29] = [
        0.0, 0.0, 0.0, 0.0, 0.0, 504.0, 504.0, 504.0, 504.0, 504.0, 1008.0, 1008.0, 1008.0, 1008.0,
        1008.0, 1512.0, 1512.0, 1512.0, 1512.0, 1512.0, 4500.0, 4059.0, 3618.0, 5155.0, -504.0,
        -504.0, -504.0, -504.0, -504.0,
    ];
    let layout_x = Array1::from_iter(x_raw.iter().map(|v| v / 1.5));
    let layout_y = Array1::from_iter(y_raw.iter().map(|v| v / 1.5));
    let n_turbines = layout_x.len();

    // Turbine weights: only optimize the first 10 turbines
    let mut turbine_weights = vec![0.0; n_turbines];
    turbine_weights[0..10].iter_mut().for_each(|w| *w = 1.0);

    // Reinitialize the FLORIS layout
    fmodel.set_layout(&layout_x, &layout_y)?;

    // Optimize the yaw angles at 8 m/s
    let time_series = TimeSeries::new(
        wind_rose.wind_directions.clone(),
        Array1::from_elem(wind_rose.wind_directions.len(), 8.0),
        Array1::from_elem(wind_rose.wind_directions.len(), 0.06),
    )?;
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    // Number of optimization conditions (the 8 m/s time series: one per wind
    // direction in the rose)
    let n_opt = fmodel.flow_field().n_findex;

    // CASE 1: account for the neighboring farm's wakes by using turbine
    // weights, and freeze the neighboring turbines' yaw (min = max = 0).
    let mut minimum_yaw_angle = Array2::zeros((n_opt, n_turbines));
    let mut maximum_yaw_angle = Array2::zeros((n_opt, n_turbines));
    for fi in 0..n_opt {
        for ti in 0..10 {
            maximum_yaw_angle[[fi, ti]] = 30.0;
        }
    }

    let config = YawOptimizationConfig {
        minimum_yaw_angle: 0.0,
        maximum_yaw_angle: 30.0,
        minimum_yaw_angle_array: Some(minimum_yaw_angle),
        maximum_yaw_angle_array: Some(maximum_yaw_angle),
        turbine_weights: Some(turbine_weights.clone()),
        ..Default::default()
    };
    let mut yaw_opt = YawOptimizationSR::new();
    let df_opt_with_neighbor = yaw_opt.optimize(&mut fmodel, Some(config))?;

    // CASE 2: ignore the neighboring farm's wakes by limiting the model to
    // the turbines of the farm being optimized.
    let mut fmodel_subset = fmodel.clone();
    fmodel_subset.set_layout(
        &Array1::from_iter(layout_x.iter().take(10).copied()),
        &Array1::from_iter(layout_y.iter().take(10).copied()),
    )?;
    let config = YawOptimizationConfig {
        minimum_yaw_angle: 0.0,
        maximum_yaw_angle: 30.0,
        ..Default::default()
    };
    let mut yaw_opt = YawOptimizationSR::new();
    let df_opt_without_neighbor = yaw_opt.optimize(&mut fmodel_subset, Some(config))?;
    // The subset model has only the 10 farm turbines; its yaw angles are
    // mapped back onto the full layout (the remaining turbines keep yaw 0).
    // The subset ran on the 8 m/s time series (one condition per wind
    // direction in the rose), so the yaw is looked up per direction and
    // then ramped over the wind speeds below.
    let mut yaw_opt_without_neighbor_full = Array2::zeros((n_opt, n_turbines));
    for fi in 0..n_opt {
        for t in 0..10 {
            yaw_opt_without_neighbor_full[[fi, t]] = df_opt_without_neighbor.yaw_angles[[fi, t]];
        }
    }

    // Baseline AEP (with the neighboring farm present, weighted)
    let mut fmodel_full = fmodel.clone();
    fmodel_full.set_wind_data(&wind_rose)?;
    fmodel_full.run()?;
    let aep_baseline =
        fmodel_full.get_farm_aep_weighted(&wind_rose, 8760.0, Some(&turbine_weights))?;

    // Apply the optimal yaw angles to the wind rose with the rule of thumb
    let n_findex = fmodel_full.flow_field().n_findex;
    let wind_speeds = fmodel_full.flow_field().wind_speeds.clone();
    let wind_directions = fmodel_full.flow_field().wind_directions.clone();
    let mut yaw_with_neighbor = Array2::zeros((n_findex, n_turbines));
    let mut yaw_without_neighbor = Array2::zeros((n_findex, n_turbines));
    for i in 0..n_findex {
        let ws = wind_speeds[i];
        let wd = wind_directions[i];

        // The optimization ran on the 8 m/s time series (one condition per
        // wind direction in the rose); look up the matching direction.
        let opt_fi = (0..df_opt_with_neighbor.yaw_angles.shape()[0])
            .find(|&k| (wind_rose.wind_directions[k] - wd).abs() < 1e-6)
            .unwrap_or(0);

        let factor = if ws < 4.0 || ws > 14.0 {
            0.0
        } else if ws < 6.0 {
            (6.0 - ws) / 2.0
        } else if ws > 12.0 {
            (14.0 - ws) / 2.0
        } else {
            1.0
        };

        for t in 0..n_turbines {
            yaw_with_neighbor[[i, t]] = df_opt_with_neighbor.yaw_angles[[opt_fi, t]] * factor;
            yaw_without_neighbor[[i, t]] = yaw_opt_without_neighbor_full[[opt_fi, t]] * factor;
        }
    }

    // AEP accounting for the neighboring farm
    fmodel_full.set_yaw_angles(yaw_with_neighbor)?;
    fmodel_full.run()?;
    let aep_opt_with_neighbor =
        fmodel_full.get_farm_aep_weighted(&wind_rose, 8760.0, Some(&turbine_weights))?;

    // AEP ignoring the neighboring farm
    fmodel_full.set_yaw_angles(yaw_without_neighbor)?;
    fmodel_full.run()?;
    let aep_opt_without_neighbor =
        fmodel_full.get_farm_aep_weighted(&wind_rose, 8760.0, Some(&turbine_weights))?;

    println!("Baseline AEP: {:.2} GWh.", aep_baseline / 1e9);
    println!(
        "Optimal AEP (Not accounting for neighboring farm): {:.2} GWh.",
        aep_opt_without_neighbor / 1e9
    );
    println!(
        "Optimal AEP (Accounting for neighboring farm): {:.2} GWh.",
        aep_opt_with_neighbor / 1e9
    );

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
