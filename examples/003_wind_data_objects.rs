//! Example 003: Wind Data Objects
//!
//! This example demonstrates the use of wind data objects in FLORUS:
//! TimeSeries, WindRose, and WindTIRose.
//!
//! For each of the WindData objects, examples are shown of:
//!    1) Initializing the object
//!    2) Broadcasting values
//!    3) Converting between objects
//!    4) Setting TI and value
//!    5) Plotting (where applicable)
//!    6) Setting the FLORIS model using the object
//!
//! Corresponds to `003_wind_data_objects.py` in Python FLORIS.

use florus::{FlorisModel, Result};
use florus::wind_data::{TimeSeries, WindRose, WindTIRose};
use florus::core::InterpMethod;
use ndarray::{Array1, Array2};

fn main() -> Result<()> {
    println!("=== FLORUS Wind Data Objects Example ===\n");

    //==========================================================
    // Initializing
    //==========================================================
    println!("=== Section 1: Initializing Wind Data Objects ===\n");

    // FLORUS provides a set of wind data objects to hold the ambient wind conditions in a
    // convenient classes that include capabilities and methods to manipulate and visualize
    // the data.

    // The TimeSeries class is used to hold time series data, such as wind speed, wind direction,
    // and turbulence intensity.

    // There is also a "value" wind data variable, which represents the value of the power
    // generated at each time step or wind condition (e.g., the price of electricity). This can
    // then be used in later optimization methods to optimize for quantities besides AEP.

    // Generate wind speeds, directions, turbulence intensities, and values via random signals
    let n = 100;
    use rand::Rng;
    let mut rng = rand::thread_rng();
    
    let wind_speeds: Vec<f64> = (0..n).map(|_| 8.0 + 2.0 * rng.gen::<f64>()).collect();
    let wind_directions: Vec<f64> = (0..n).map(|_| 270.0 + 30.0 * rng.gen::<f64>()).collect();
    let turbulence_intensities: Vec<f64> = (0..n).map(|_| 0.06 + 0.02 * rng.gen::<f64>()).collect();
    let values: Vec<f64> = (0..n).map(|_| 25.0 + 10.0 * rng.gen::<f64>()).collect();

    println!("1. Creating TimeSeries object...");
    let time_series = TimeSeries::new(
        Array1::from(wind_directions.clone()),
        Array1::from(wind_speeds.clone()),
        Array1::from(turbulence_intensities.clone()),
    )?;
    println!("   ✓ TimeSeries created with {} time steps", time_series.wind_directions.len());

    // The WindRose class is used to hold wind rose data, such as wind speed, wind direction,
    // and frequency. TI and value are represented as bin averages per wind direction and
    // speed bin.
    println!("\n2. Creating WindRose object...");
    let wd_bins: Vec<f64> = (0..360).step_by(3).map(|x| x as f64).collect();
    let ws_bins: Vec<f64> = (4..20).step_by(2).map(|x| x as f64).collect();
    
    let n_dir = wd_bins.len();
    let n_ws = ws_bins.len();
    
    // Make TI table 6% TI for all wind directions and speeds
    let ti_table = Array2::from_elem((n_dir, n_ws), 0.06);
    
    // Make value table 25 for all wind directions and speeds
    let value_table = Array2::from_elem((n_dir, n_ws), 25.0);
    
    // Uniform frequency
    let mut freq_table = Array2::from_elem((n_dir, n_ws), 1.0);
    let total: f64 = freq_table.sum();
    freq_table.mapv_inplace(|x| x / total);

    let wind_rose = WindRose::new(
        Array1::from(wd_bins.clone()),
        Array1::from(ws_bins.clone()),
        ti_table,
        Some(freq_table),
        Some(value_table),
        false,  // compute_zero_freq_occurrence
        None,   // heterogeneous_map
        None,   // multidim_conditions
    )?;
    println!("   ✓ WindRose created with {} directions × {} speeds", n_dir, n_ws);

    // The WindTIRose class is similar to the WindRose table except that TI is also binned
    // making the frequency table a 3D array.
    println!("\n3. Creating WindTIRose object...");
    let ti_bins: Vec<f64> = (5..15).map(|x| x as f64 / 100.0).collect();
    let n_ti = ti_bins.len();
    
    // Uniform TI table (3D)
    let ti_table_3d = ndarray::Array3::from_elem((n_dir, n_ws, n_ti), 1.0);
    
    // Uniform frequency
    let freq_table_2d = Array2::from_elem((n_dir, n_ws), 1.0 / (n_dir * n_ws) as f64);
    
    // Uniform value
    let value_table_2d = Array2::from_elem((n_dir, n_ws), 25.0);

    let wind_ti_rose = WindTIRose::new(
        Array1::from(wd_bins.clone()),
        Array1::from(ws_bins.clone()),
        Array1::from(ti_bins),
        ti_table_3d,
        Some(freq_table_2d),
        Some(value_table_2d),
    )?;
    println!("   ✓ WindTIRose created with {} directions × {} speeds × {} TI bins", n_dir, n_ws, n_ti);

    //==========================================================
    // Broadcasting
    //==========================================================
    println!("\n\n=== Section 2: Broadcasting ===\n");

    // A convenience method of the wind data objects is that, unlike the lower-level
    // FlorisModel.set() method, the wind data objects can broadcast upward data provided
    // as a scalar to the full array. This is useful for setting the same wind conditions
    // for all turbines in a wind farm.

    // For TimeSeries, as long as one condition is given as an array, the other 2
    // conditions can be given as scalars. The TimeSeries object will broadcast the
    // scalars to the full array (uniform)
    println!("4. Creating TimeSeries with broadcasting...");
    let wind_dirs_varied: Vec<f64> = (0..n).map(|_| 270.0 + 30.0 * rng.gen::<f64>()).collect();
    let wind_speeds_uniform = vec![8.0; n];
    let tis_uniform = vec![0.06; n];
    let time_series_broadcast = TimeSeries::new(
        Array1::from(wind_dirs_varied),
        Array1::from(wind_speeds_uniform),
        Array1::from(tis_uniform),
    )?;
    println!("   ✓ TimeSeries created with varied wind directions, uniform speed and TI");
    println!("      Wind speeds: all {:.1} m/s", time_series_broadcast.wind_speeds[0]);
    println!("      TIs: all {:.2}", time_series_broadcast.turbulence_intensities[0]);

    // For WindRose, wind directions and wind speeds must be given as arrays, but the
    // ti_table can be supplied as a scalar which will apply uniformly to all wind
    // directions and speeds. Not supplying a freq table will similarly generate
    // a uniform frequency table.
    println!("\n5. Creating WindRose with simplified inputs...");
    let wd_bins_simple: Vec<f64> = (0..360).step_by(3).map(|x| x as f64).collect();
    let ws_bins_simple: Vec<f64> = (4..20).step_by(2).map(|x| x as f64).collect();
    let n_dir_simple = wd_bins_simple.len();
    let n_ws_simple = ws_bins_simple.len();
    let ti_table_uniform = Array2::from_elem((n_dir_simple, n_ws_simple), 0.06);
    
    let wind_rose_simple = WindRose::new(
        Array1::from(wd_bins_simple),
        Array1::from(ws_bins_simple),
        ti_table_uniform,
        None,  // No freq table (will use uniform)
        None,  // No value table
        false,
        None,
        None,
    )?;
    println!("   ✓ WindRose created with uniform TI = {:.2}", wind_rose_simple.ti_table[(0, 0)]);

    //==========================================================
    // Wind Rose from Time Series
    //==========================================================
    println!("\n\n=== Section 3: Wind Rose from Time Series ===\n");

    // The TimeSeries class has a method to generate a wind rose from a time series based on binning
    println!("6. Converting TimeSeries to WindRose...");
    let wind_rose_from_ts = time_series.to_wind_rose(3.0, 2.0);
    println!("   ✓ WindRose created from TimeSeries");
    println!("      Directions: {} bins", wind_rose_from_ts.wind_directions.len());
    println!("      Speeds: {} bins", wind_rose_from_ts.wind_speeds.len());

    //==========================================================
    // Aggregating and Resampling the Wind Rose
    //==========================================================
    println!("\n\n=== Section 4: Aggregating and Resampling ===\n");

    // The downsample function allows for aggregation of the wind rose data into
    // fewer wind direction and wind speed bins.
    println!("7. Downsampling WindRose...");
    let wind_rose_agg = wind_rose.downsample(Some(10.0), Some(2.0), None);
    println!("   ✓ WindRose downsampled");
    println!("      Original: {}×{} bins", wind_rose.wind_directions.len(), wind_rose.wind_speeds.len());
    println!("      Aggregated: {}×{} bins", wind_rose_agg.wind_directions.len(), wind_rose_agg.wind_speeds.len());

    // For upsampling, the upsample function can be used to interpolate
    // the wind rose data to a finer grid.
    println!("\n8. Upsampling WindRose...");
    let wind_rose_resample = wind_rose.upsample(0.5, 0.25, &InterpMethod::Linear);
    println!("   ✓ WindRose upsampled");
    println!("      Original: {}×{} bins", wind_rose.wind_directions.len(), wind_rose.wind_speeds.len());
    println!("      Resampled: {}×{} bins", wind_rose_resample.wind_directions.len(), wind_rose_resample.wind_speeds.len());

    //==========================================================
    // Setting turbulence intensity
    //==========================================================
    println!("\n\n=== Section 5: Setting Turbulence Intensity ===\n");

    // Each of the wind data objects also has the ability to set the turbulence intensity
    // according to a function of wind speed and direction. This can be done using a custom
    // function by using the assign_ti_using_wd_ws_function method. There is also a method
    // called assign_ti_using_IEC_method which assigns TI based on the IEC 61400-1 standard.
    println!("9. Assigning TI using IEC method...");
    let mut wind_rose_iec = wind_rose.clone();
    wind_rose_iec.assign_ti_using_iec_method(None);
    println!("   ✓ TI assigned using IEC 61400-1 standard");
    println!("      Sample TI at 270°, 8 m/s: {:.3}", 
             wind_rose_iec.ti_table[(
                 wind_rose_iec.wind_directions.iter().position(|&x| (x - 270.0).abs() < 0.1).unwrap_or(0),
                 wind_rose_iec.wind_speeds.iter().position(|&x| (x - 8.0).abs() < 0.1).unwrap_or(0)
             )]);

    //==========================================================
    // Setting value
    //==========================================================
    println!("\n\n=== Section 6: Setting Value ===\n");

    // Similarly, each of the wind data objects also has the ability to set the value according to
    // a function of wind speed and direction. This can be done using a custom function by using
    // the assign_value_using_wd_ws_function method. There is also a method called
    // assign_value_piecewise_linear which assigns value based on a linear piecewise function of
    // wind speed.
    println!("10. Assigning value using piecewise linear function...");
    let mut wind_rose_value = wind_rose.clone();
    // Parameters: value at zero ws, ws knee point, slope 1, slope 2, limit to zero, normalize
    wind_rose_value.assign_value_piecewise_linear(0.0, 5.0, 1.0, -0.5, false, false);
    println!("   ✓ Value assigned using piecewise linear function");
    println!("      (Approximates normalized mean electricity price vs. wind speed)");

    //==========================================================
    // Setting the FLORIS model via wind data
    //==========================================================
    println!("\n\n=== Section 7: Setting FLORIS Model via Wind Data ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set the wind conditions using the TimeSeries object
    println!("11. Setting FlorisModel with TimeSeries...");
    fmodel.set_wind_data(&time_series)?;
    println!("   ✓ Model set with TimeSeries ({} conditions)", fmodel.n_findex());

    // Set the wind conditions using the WindRose object
    println!("\n12. Setting FlorisModel with WindRose...");
    fmodel.set_wind_data(&wind_rose)?;
    println!("   ✓ Model set with WindRose ({} conditions)", fmodel.n_findex());

    // Note that in the case of the wind_rose, under the default settings, wind direction and wind speed
    // bins for which frequency is zero are not simulated. This can be changed by setting the
    // compute_zero_freq_occurrence parameter to True.
    println!("\n13. Testing compute_zero_freq_occurrence parameter...");
    let wind_directions_test = vec![200.0, 300.0];
    let wind_speeds_test = vec![5.0, 10.0];
    // Use non-zero frequencies for all bins
    let freq_table_test = ndarray::arr2(&[[0.25, 0.25], [0.25, 0.25]]);
    
    let wind_rose_zero_freq = WindRose::new(
        Array1::from(wind_directions_test.clone()),
        Array1::from(wind_speeds_test.clone()),
        Array2::from_elem((2, 2), 0.06),
        Some(freq_table_test.clone()),
        None,
        false,  // compute_zero_freq_occurrence = False
        None,
        None,
    )?;
    fmodel.set_wind_data(&wind_rose_zero_freq)?;
    println!("   With non-zero frequencies: {} conditions", fmodel.n_findex());

    // Set the wind conditions using the WindTIRose object
    println!("\n14. Setting FlorisModel with WindTIRose...");
    fmodel.set_wind_data(&wind_ti_rose)?;
    println!("   ✓ Model set with WindTIRose ({} conditions)", fmodel.n_findex());

    println!("\n=== Wind Data Objects Example Complete ===");
    println!("\nSummary:");
    println!("  - TimeSeries: For time-series wind data");
    println!("  - WindRose: For binned wind climate data");
    println!("  - WindTIRose: For binned data with TI distribution");
    println!("  - All support broadcasting, conversion, and model integration");

    Ok(())
}
