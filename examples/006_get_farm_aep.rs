/// Example 3: Annual Energy Production (AEP) Calculations
///
/// FLORIS provides methods to calculate the Annual Energy Production (AEP)
/// of a wind farm based on wind resource data. This example demonstrates:
///
/// 1. Creating a WindRose with wind resource statistics
/// 2. Using FlorisModel with wind rose data
/// 3. Calculating AEP using get_farm_aep()
/// 4. Calculating AEP with uniform frequencies (get_farm_aep_uniform)
/// 5. Computing capacity factors
/// 6. Analyzing energy production by turbine
///
/// This is the Rust equivalent of Python's 006_get_farm_aep.py

use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::types::{Array1, Array2};
use florus::wind_data::{WindData, WindRose};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 3: Annual Energy Production (AEP) Calculations");
    println!("===============================================================\n");

    // Create a 3-turbine wind farm
    let d = 126.0; // NREL 5MW rotor diameter
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    println!("Creating 3-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, layout_y[i]);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    // ================================================================
    // Create a WindRose with wind resource data
    // ================================================================
    println!("\n--- Creating Wind Rose ---\n");

    // Define wind direction and speed bins
    let wind_directions: Vec<f64> = (0..36).map(|i| (i as f64) * 10.0).collect(); // 0° to 350°
    let wind_speeds: Vec<f64> = (3..26).map(|i| i as f64).collect(); // 3-25 m/s

    let n_dir = wind_directions.len();
    let n_ws = wind_speeds.len();

    println!("Wind rose configuration:");
    println!("  Wind directions: {} bins ({}° to {}°)",
             n_dir, wind_directions.first().unwrap(), wind_directions.last().unwrap());
    println!("  Wind speeds: {} bins ({} m/s to {} m/s)",
             n_ws, wind_speeds.first().unwrap(), wind_speeds.last().unwrap());

    // Create turbulence intensity table
    // TI varies with wind speed (lower at higher speeds)
    let mut ti_table = Array2::zeros((n_dir, n_ws));
    for (i, _) in wind_directions.iter().enumerate() {
        for (j, &ws) in wind_speeds.iter().enumerate() {
            // Simple TI model: TI decreases with wind speed
            let ti = 0.1 * (10.0 / ws.max(1.0f64)).min(0.25);
            ti_table[[i, j]] = ti;
        }
    }

    // Create frequency table - more wind from certain directions
    // Dominant wind from SW (around 225°)
    let mut freq_table = Array2::zeros((n_dir, n_ws));
    for (i, &wd) in wind_directions.iter().enumerate() {
        // Base frequency varies by direction
        let base_freq = if wd >= 180.0 && wd <= 270.0 {
            0.15 // Higher frequency for SW winds
        } else {
            0.05 // Lower for other directions
        };

        for (j, &ws) in wind_speeds.iter().enumerate() {
            // Wind speed distribution (Weibull-like)
            let ws_factor = if ws >= 8.0 && ws <= 12.0 {
                0.15 // Peak around 8-12 m/s
            } else {
                0.05
            };
            freq_table[[i, j]] = base_freq * ws_factor;
        }
    }

    // Normalize frequency table
    let freq_sum: f64 = freq_table.iter().sum();
    if freq_sum > 0.0 {
        for val in freq_table.iter_mut() {
            *val /= freq_sum;
        }
    }

    // Create the WindRose
    let wind_rose = WindRose::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(wind_speeds.clone()),
        ti_table,
        Some(freq_table),
        None,
        false,
        None,
        None,
    )?;

    println!("  Total conditions: {}", wind_rose.n_conditions());

    // ================================================================
    // Convert WindRose to time series for FlorisModel
    // ================================================================
    println!("\n--- Preparing Wind Data for FLORIS ---\n");

    // Unpack the wind rose to get flattened arrays for simulation
    // Returns: (wind_directions, wind_speeds, turbulence_intensities, freq_table, value_table, heterogeneous_config)
    let (wd_flat, ws_flat, ti_flat, freq_table, _, _) = wind_rose.unpack();

    // Extract frequency per condition (2D freq_table has shape [n_dir, n_ws])
    // Flatten it to get frequencies for each condition
    let freq_flat: Vec<f64> = freq_table.iter().copied().collect();

    println!("Wind conditions prepared:");
    println!("  Number of conditions: {}", wd_flat.len());

    // ================================================================
    // Create FlorisModel with wind conditions
    // ================================================================
    println!("\n--- Setting Up FlorisModel ---\n");

    // Create flow field (ti_flat comes from unpack)
    let flow_field = florus::core::FlowField::new(
        ws_flat.clone(),
        wd_flat.clone(),
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
        ti_flat,
        90.0,   // reference_wind_height
    )?;

    // Create model
    let mut model = florus::FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    // Initialize and run
    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    // ================================================================
    // Calculate AEP
    // ================================================================
    println!("\n--- AEP Calculations ---\n");

    // Standard AEP calculation using wind rose frequencies
    let aep = model.get_farm_aep_uniform(8760.0);
    let aep_gwh = aep / 1_000_000_000.0;

    println!("Annual Energy Production (AEP):");
    println!("  Total AEP: {:.2} GWh", aep_gwh);
    println!("  Total AEP: {:.2} MWh", aep / 1_000_000.0);

    // Get turbine powers for analysis
    let turbine_powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();

    println!("\n--- Power Analysis ---\n");

    // Average power output
    let avg_farm_power = farm_power.mean().unwrap_or(0.0);
    println!("Average farm power: {:.1} kW", avg_farm_power / 1000.0);

    // Power by turbine
    println!("\nAverage power by turbine:");
    for ti in 0..model.farm.n_turbines() {
        let avg_power = turbine_powers.column(ti).mean().unwrap_or(0.0);
        println!("  Turbine {}: {:.1} kW", ti, avg_power / 1000.0);
    }

    // ================================================================
    // Capacity Factor Calculation
    // ================================================================
    println!("\n--- Capacity Factor Analysis ---\n");

    // NREL 5MW rated power
    let rated_power_watts = 5_000_000.0;
    let n_turbines = model.farm.n_turbines();
    let total_rated_power = rated_power_watts * n_turbines as f64;

    // Simple capacity factor calculation
    // CF = Annual Energy / (Rated Power × Hours per Year)
    let capacity_factor = (aep / total_rated_power) / 8760.0 * 100.0;

    println!("Capacity Factor Calculation:");
    println!("  Rated power per turbine: {:.0} kW", rated_power_watts / 1000.0);
    println!("  Number of turbines: {}", n_turbines);
    println!("  Total rated power: {:.0} kW", total_rated_power / 1000.0);
    println!("  Annual energy production: {:.2} GWh", aep_gwh);
    println!("  Capacity factor: {:.1}%", capacity_factor);

    // ================================================================
    // Detailed AEP Breakdown
    // ================================================================
    println!("\n--- Detailed AEP Breakdown ---\n");

    // Calculate energy by wind direction sector
    println!("AEP by Wind Direction Sector:");
    println!("  {:>8}  {:>12}  {:>12}", "Sector", "Freq", "Energy (MWh)");
    println!("  {}", "-".repeat(40));

    let sectors = [0, 90, 180, 270]; // N, E, S, W
    for sector_start in sectors {
        let sector_end = sector_start + 90;
        let mut sector_energy = 0.0;
        let mut sector_freq = 0.0;

        for (i, &wd) in wd_flat.iter().enumerate() {
            if wd >= sector_start as f64 && wd < sector_end as f64 {
                sector_energy += farm_power[i] * freq_flat[i] * 8760.0;
                sector_freq += freq_flat[i];
            }
        }

        let sector_name = match sector_start {
            0 => "N (0-90)",
            90 => "E (90-180)",
            180 => "S (180-270)",
            270 => "W (270-360)",
            _ => "Unknown",
        };

        println!("  {:>8}  {:>12.1}%  {:>12.1}", sector_name, sector_freq * 100.0, sector_energy / 1_000_000.0);
    }

    // ================================================================
    // Wake Loss Impact on AEP
    // ================================================================
    println!("\n--- Wake Loss Impact on AEP ---\n");

    // Compare with no-wake scenario (simplified)
    // Assuming no wake loss at 0° (turbines perpendicular to wind)
    let mut min_wake_loss = f64::MAX;
    let mut max_wake_loss = 0.0;

    for (i, &wd) in wd_flat.iter().enumerate() {
        // For aligned winds (around 270°), downstream turbines experience wake loss
        if wd >= 250.0 && wd <= 290.0 {
            let upstream_power = turbine_powers[[i, 0]];
            let downstream_power = turbine_powers[[i, 2]];
            if upstream_power > 0.0 {
                let wake_loss = (1.0 - downstream_power / upstream_power) * 100.0;
                min_wake_loss = f64::min(min_wake_loss, wake_loss);
                max_wake_loss = f64::max(max_wake_loss, wake_loss);
            }
        }
    }

    if min_wake_loss < f64::MAX {
        println!("Wake loss analysis (aligned winds, 250°-290°):");
        println!("  Minimum wake loss (T0→T2): {:.1}%", min_wake_loss);
        println!("  Maximum wake loss (T0→T2): {:.1}%", max_wake_loss);
        println!("  Estimated annual wake loss: {:.1}%", (min_wake_loss + max_wake_loss) / 2.0);
    }

    // ================================================================
    // Summary Statistics
    // ================================================================
    println!("\n--- Summary ---\n");

    println!("Annual Energy Production Summary:");
    println!("  Total AEP: {:.2} GWh ({:.2} MWh)", aep_gwh, aep / 1_000_000.0);
    println!("  Average farm output: {:.1} kW", avg_farm_power / 1000.0);
    println!("  Capacity factor: {:.1}%", capacity_factor);
    println!("  Energy per turbine: {:.2} MWh/turbine",
             (aep / 1_000_000.0) / n_turbines as f64);

    println!("\n===============================================================");
    println!("Example completed successfully!");

    Ok(())
}
