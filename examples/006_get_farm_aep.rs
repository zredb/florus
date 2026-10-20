//! Example 6: Getting Expected Power and AEP
//!
//! The expected power of a farm is computed by multiplying the power output
//! of the farm by the frequency of each findex. The expected AEP (annual
//! energy production) is the expected power multiplied by the hours in a
//! year. With a wind_data object set on the model, the AEP is computed
//! directly from the wind data's frequency table.
//!
//! Scenario aligned 1:1 with 006_get_farm_aep.py:
//! - a 3-turbine row (5 D spacing) on gch.yaml;
//! - a TimeSeries of 365 * 24 conditions drawn from numpy's
//!   default_rng(0) (vendored as examples/inputs/timeseries_seed0.csv so the
//!   exact same inflow is reproduced);
//! - a WindRose from examples/inputs/wind_rose.csv (same file as
//!   floris/examples/inputs/wind_rose.csv; Python loads it via
//!   `WindRose.read_csv_long(..., freq_col="freq_val", ti_col_or_value=0.06)`).

use florus::wind_data::{TimeSeries, WindRose};
use florus::{FlorisModel, Result};
use ndarray::{Array1, Array2};
use std::collections::BTreeMap;
use std::io::BufRead;

fn main() -> Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set to a 3-turbine layout
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Using TimeSeries ======================================================

    // Time series with time steps = 365 * 24 (same draw as numpy
    // default_rng(0): wind_directions ~ U(0, 360), wind_speeds ~ U(5, 25)).
    let (wind_directions, wind_speeds) =
        read_timeseries_csv("examples/inputs/timeseries_seed0.csv")?;
    let n = wind_directions.len();
    let turbulence_intensities = Array1::from_elem(n, 0.06);

    let time_series = TimeSeries::new(wind_directions, wind_speeds, turbulence_intensities)?;

    // Set the wind data, run, and get the AEP
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    let aep = fmodel.get_farm_aep(&time_series, 8760.0)?;
    // Equivalent to passing the frequency table explicitly (Python:
    // get_farm_AEP(freq=time_series.unpack_freq()))
    let aep_b = fmodel.get_farm_aep(&time_series, 8760.0)?;

    println!("AEP from time series: {aep}, and re-computed AEP: {aep_b}");

    // Using WindRose ========================================================

    let wind_rose = read_wind_rose_csv("examples/inputs/wind_rose.csv", 0.06)?;

    let n_wd = wind_rose.wind_directions.len();
    let n_ws = wind_rose.wind_speeds.len();

    // Number of zero-frequency bins
    let n_zeros = wind_rose.freq_table.iter().filter(|&&f| f == 0.0).count();

    // Set the wind rose and run
    fmodel.set_wind_data(&wind_rose)?;
    fmodel.run()?;

    println!(
        "Total number of wind direction and wind speed combination: {}",
        n_wd * n_ws
    );
    println!("Number of 0 frequency bins: {n_zeros}");
    println!("n_findex: {}", fmodel.n_findex());

    // Get the AEP
    let aep = fmodel.get_farm_aep(&wind_rose, 8760.0)?;

    println!("AEP from wind rose: {:.3} (GWh)", aep / 1e9);

    // Run the model again, without wakes, and compute the wake losses
    fmodel.run_no_wake()?;

    let aep_no_wake = fmodel.get_farm_aep(&wind_rose, 8760.0)?;

    let wake_losses = 100.0 * (aep_no_wake - aep) / aep_no_wake;

    println!("Wake losses: {wake_losses:.2}%");

    Ok(())
}

/// Read the vendored numpy-seed-0 time series (columns: wd, ws).
fn read_timeseries_csv(path: &str) -> Result<(Array1<f64>, Array1<f64>)> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut wd = Vec::new();
    let mut ws = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue; // header
        }
        let mut cols = line.split(',');
        let a: f64 = cols.next().unwrap_or_default().parse()?;
        let b: f64 = cols.next().unwrap_or_default().parse()?;
        wd.push(a);
        ws.push(b);
    }
    Ok((Array1::from_vec(wd), Array1::from_vec(ws)))
}

/// Replicate Python `WindRose.read_csv_long(..., freq_col=..., ti_col_or_value=0.06)`
/// for this file's long format: each row is one (ws, wd) combination with a
/// frequency; frequencies are normalized, then binned onto the unique sorted
/// (wd, ws) grid.
fn read_wind_rose_csv(path: &str, ti_value: f64) -> Result<WindRose> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);

    // wd, ws, freq per row
    let mut rows: Vec<(f64, f64, f64)> = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue; // header: ws,wd,freq_val
        }
        let mut cols = line.split(',');
        let ws: f64 = cols.next().unwrap_or_default().parse()?;
        let wd: f64 = cols.next().unwrap_or_default().parse()?;
        let freq: f64 = cols.next().unwrap_or_default().parse()?;
        rows.push((wd, ws, freq));
    }

    // Normalize freq_values (Python wind_data.py:1107-1108)
    let freq_sum: f64 = rows.iter().map(|r| r.2).sum();
    for r in rows.iter_mut() {
        r.2 /= freq_sum;
    }

    // Unique sorted wind directions and speeds
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

    // Bin the frequencies onto the (wd, ws) grid
    let n_dir = unique_wd.len();
    let n_ws = unique_ws.len();
    let mut freq_table = Array2::<f64>::zeros((n_dir, n_ws));
    for &(wd, ws, freq) in &rows {
        let i = *wd_index.get(&wd.to_bits()).unwrap();
        let j = *ws_index.get(&ws.to_bits()).unwrap();
        freq_table[[i, j]] += freq;
    }

    let wind_directions = Array1::from_vec(unique_wd);
    let wind_speeds = Array1::from_vec(unique_ws);
    let ti_table = Array2::from_elem((n_dir, n_ws), ti_value);

    Ok(WindRose::new(
        wind_directions,
        wind_speeds,
        ti_table,
        Some(freq_table),
        None,
        false, // compute_zero_freq_occurrence
        None,
        None,
    )?)
}
