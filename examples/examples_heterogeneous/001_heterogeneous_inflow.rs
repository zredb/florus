/// Example 12: Heterogeneous Wind Inflow
///
/// This example demonstrates FLORIS-RS heterogeneous inflow capabilities:
/// 1. Creating a HeterogeneousMap with spatial wind speed variations
/// 2. Using heterogeneous maps with WindRose
/// 3. Comparing homogeneous vs heterogeneous wind conditions
/// 4. Modeling terrain effects and speed-ups
///
/// Heterogeneous inflow represents spatial variations in wind speed
/// across a wind farm, such as those caused by terrain or obstacles.

use florus::core::{Farm, FlowField};
use florus::heterogeneous_map::HeterogeneousMap;
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;
use florus::core::base::InterpMethod;
use florus::FlorisModel;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 12: Heterogeneous Wind Inflow");
    println!("================================================\n");

    // Create a 4-turbine wind farm
    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d, 15.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    println!("Creating 4-turbine wind farm:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, layout_y[i]);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    // ============================================================
    // HOMOGENEOUS BASELINE
    // ============================================================
    println!("\n--- Homogeneous Wind Baseline ---\n");

    // Standard wind rose without heterogeneous map
    let wind_directions = Array1::from_vec(vec![270.0, 280.0, 290.0]);
    let wind_speeds = Array1::from_vec(vec![8.0, 10.0, 12.0]);
    let ti_table = Array2::from_elem((3, 3), 0.06);

    let wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        ti_table,
        None,
        None,
        false,
        None,
        None,
    )?;

    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.12,
        1.225,
        Array1::from_vec(vec![0.06, 0.07, 0.08]),
        90.0,
    )?;

    let mut model = FlorisModel {
        farm: farm.clone(),
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let homogeneous_powers = model.get_turbine_powers();
    let homogeneous_farm_power: f64 = homogeneous_powers.iter().sum();

    println!("Homogeneous Results:");
    for fi in 0..3 {
        let findex_power: f64 = (0..4).map(|ti| homogeneous_powers[[fi, ti]]).sum();
        println!("  Wind {:.0} m/s, {:.0} deg: {:.2} MW", 
            wind_speeds[fi], wind_directions[fi], findex_power / 1_000_000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", homogeneous_farm_power / 1_000_000.0);

    // ============================================================
    // CREATE HETEROGENEOUS MAP
    // ============================================================
    println!("--- Creating Heterogeneous Map ---\n");

    println!("Heterogeneous inflow models spatial wind speed variations:");
    println!("  - Terrain effects (hills, ridges)");
    println!("  - Speed-ups near cliffs or escarpments");
    println!("  - Slowdowns behind obstacles");
    println!("  - Forest or vegetation effects\n");

    // Create a grid of points for the heterogeneous map
    // This models a speed-up effect on the left side of the farm
    let het_x = Array1::from_vec(vec![-500.0, 0.0, 500.0, 1000.0, 1500.0]);
    let het_y = Array1::from_vec(vec![-500.0, -500.0, -500.0, -500.0, -500.0]);

    // Speed multipliers: 1.0 = no change, >1.0 = speed-up, <1.0 = slow-down
    // This represents a speed-up near x=0 (e.g., a hill or ridge)
    let speed_multipliers = Array2::from_shape_vec(
        (1, 5),
        vec![
            1.15, 1.10, 1.00, 0.95, 0.90,  // For 270 deg wind
        ],
    )?;

    // Create heterogeneous map with wind direction-specific multipliers
    let het_x_directions = Array1::from_vec(vec![-500.0, 0.0, 500.0, 1000.0, 1500.0]);
    let het_y_directions = Array1::from_vec(vec![-500.0, -500.0, -500.0, -500.0, -500.0]);

    let speed_multipliers_directions = Array2::from_shape_vec(
        (3, 5),
        vec![
            1.15, 1.10, 1.00, 0.95, 0.90,  // For 270 deg
            1.12, 1.08, 1.02, 0.97, 0.92, // For 280 deg
            1.10, 1.05, 1.00, 0.98, 0.95, // For 290 deg
        ],
    )?;

    let heterogeneous_map = HeterogeneousMap::new(
        het_x_directions.clone(),
        het_y_directions.clone(),
        speed_multipliers_directions,
        None,
        Some(Array1::from_vec(vec![270.0, 280.0, 290.0])),  // Wind directions
        None,  // No wind speed dependence
        InterpMethod::Linear,
    )?;

    println!("Heterogeneous Map Configuration:");
    println!("  Grid points: {:?}", het_x_directions);
    println!("  Y positions: {:?}", het_y_directions);
    println!("  Wind directions: [270, 280, 290] deg\n");

    println!("Speed Multipliers by Wind Direction:");
    println!("  270 deg: [1.15, 1.10, 1.00, 0.95, 0.90]");
    println!("  280 deg: [1.12, 1.08, 1.02, 0.97, 0.92]");
    println!("  290 deg: [1.10, 1.05, 1.00, 0.98, 0.95]\n");

    // ============================================================
    // HETEROGENEOUS WIND ROSE
    // ============================================================
    println!("--- Wind Rose with Heterogeneous Map ---\n");

    let wind_rose_het = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        Array2::from_elem((3, 3), 0.06),
        None,
        None,
        false,
        Some(heterogeneous_map),
        None,
    )?;

    let flow_field_het = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,
        0.12,
        1.225,
        Array1::from_vec(vec![0.06, 0.07, 0.08]),
        90.0,
    )?;

    let mut model_het = FlorisModel {
        farm: farm.clone(),
        flow_field: flow_field_het,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model_het.initialize_grid()?;
    model_het.initialize_flow_field()?;
    model_het.run()?;

    let heterogeneous_powers = model_het.get_turbine_powers();
    let heterogeneous_farm_power: f64 = heterogeneous_powers.iter().sum();

    println!("Heterogeneous Results:");
    for fi in 0..3 {
        let findex_power: f64 = (0..4).map(|ti| heterogeneous_powers[[fi, ti]]).sum();
        println!("  Wind {:.0} m/s, {:.0} deg: {:.2} MW", 
            wind_speeds[fi], wind_directions[fi], findex_power / 1_000_000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", heterogeneous_farm_power / 1_000_000.0);

    // ============================================================
    // COMPARISON
    // ============================================================
    println!("--- Comparison: Homogeneous vs Heterogeneous ---\n");

    let power_change = heterogeneous_farm_power - homogeneous_farm_power;
    let percent_change = if homogeneous_farm_power > 0.0 {
        (power_change / homogeneous_farm_power) * 100.0
    } else {
        0.0
    };

    println!("{:>25} {:>15} {:>15}", "", "Homogeneous", "Heterogeneous");
    println!("{}", "-".repeat(55));
    println!("{:>25} {:>15.2} {:>15.2}", "Total Farm Power (MW):", 
        homogeneous_farm_power / 1_000_000.0, 
        heterogeneous_farm_power / 1_000_000.0);

    for fi in 0..3 {
        let homo_power: f64 = (0..4).map(|ti| homogeneous_powers[[fi, ti]]).sum();
        let het_power: f64 = (0..4).map(|ti| heterogeneous_powers[[fi, ti]]).sum();
        println!("{:>25} {:>15.2} {:>15.2}", 
            format!("Findex {} ({} m/s):", fi, wind_speeds[fi]),
            homo_power / 1_000_000.0,
            het_power / 1_000_000.0);
    }

    println!("\nPower Change: {:+.2} MW ({:+.1}%)\n", 
        power_change / 1_000_000.0, percent_change);

    // ============================================================
    // TURBINE-LEVEL COMPARISON
    // ============================================================
    println!("--- Turbine-Level Power Comparison ---\n");

    println!("Position effects on individual turbines:");
    println!("  Turbine 0 (x=0): Experiences speed-up from incoming wind");
    println!("  Turbine 2 (x=10D): Near neutral zone");
    println!("  Turbine 3 (x=15D): May experience slight slow-down\n");

    println!("{:>10} {:>12} {:>12} {:>12}", "Turbine", "Homo (kW)", "Het (kW)", "Change (%)");
    println!("{}", "-".repeat(50));

    for ti in 0..4 {
        let homo_power: f64 = (0..3).map(|fi| homogeneous_powers[[fi, ti]]).sum::<f64>() / 3.0;
        let het_power: f64 = (0..3).map(|fi| heterogeneous_powers[[fi, ti]]).sum::<f64>() / 3.0;
        let change_pct = if homo_power > 0.0 {
            ((het_power - homo_power) / homo_power) * 100.0
        } else {
            0.0
        };
        println!("{:>10} {:>12.1} {:>12.1} {:>+12.1}", 
            format!("T{}", ti), homo_power / 1000.0, het_power / 1000.0, change_pct);
    }

    // ============================================================
    // USE CASES FOR HETEROGENEOUS MAPS
    // ============================================================
    println!("\n--- Use Cases for Heterogeneous Inflow ---\n");

    println!("Common applications of heterogeneous wind maps:");
    println!("\n  1. TERRAIN EFFECTS:");
    println!("     - Hills and ridges cause speed-ups on windward slopes");
    println!("     - Valleys may have sheltering effects");
    println!("     - Cliffs and escarpments create acceleration zones\n");

    println!("  2. OFFSHORE APPLICATIONS:");
    println!("     - Land/sea breeze effects near coastlines");
    println!("     - Atmospheric stability variations");
    println!("     - Fetch effects based on wind direction\n");

    println!("  3. FOREST AND VEGETATION:");
    println!("     - Wind speed reduction in forested areas");
    println!("     - Edge effects at forest boundaries");
    println!("     - Canopy gaps and clearings\n");

    println!("  4. OBSTACLE EFFECTS:");
    println!("     - Buildings, structures, vegetation");
    println!("     - Wake effects from non-turbine obstacles");
    println!("     - Wind breaks and shelterbelts\n");

    println!("  5. COMPLEX TERRAIN:");
    println!("     - 3D wind field modeling");
    println!("     - Channeling effects in valleys");
    println!("     - Rotational effects around terrain features\n");

    // ============================================================
    // CREATING HETEROGENEOUS MAPS
    // ============================================================
    println!("--- Creating Heterogeneous Maps ---\n");

    println!("HeterogeneousMap Constructor Parameters:");
    println!("  - x, y: Grid coordinates [n_points]");
    println!("  - speed_multipliers: Array2 [n_conditions, n_points]");
    println!("  - z: Optional height coordinates for 3D grids");
    println!("  - wind_directions: Optional wind directions for condition selection");
    println!("  - wind_speeds: Optional wind speeds for condition selection");
    println!("  - interp_method: Linear or nearest-neighbor interpolation\n");

    // Example: Creating a constant speed multiplier map
    println!("Example: Simple uniform speed multiplier");
    let simple_x = Array1::from_vec(vec![0.0, 100.0, 200.0]);
    let simple_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let simple_multipliers = Array2::from_shape_vec((1, 3), vec![1.05, 1.05, 1.05])?;

    let simple_het_map = HeterogeneousMap::new(
        simple_x,
        simple_y,
        simple_multipliers,
        None,
        None,
        None,
        InterpMethod::Linear,
    )?;

    println!("  Created: 3-point grid with 5% uniform speed-up\n");

    // Example: Creating a 2D grid with varying multipliers
    println!("Example: 2D grid with spatial variation");
    let grid_x = Array1::from_vec(vec![0.0, 500.0, 1000.0, 0.0, 500.0, 1000.0]);
    let grid_y = Array1::from_vec(vec![0.0, 0.0, 0.0, 500.0, 500.0, 500.0]);
    let grid_multipliers = Array2::from_shape_vec(
        (1, 6),
        vec![1.20, 1.10, 1.00, 1.15, 1.08, 0.95],
    )?;

    let grid_het_map = HeterogeneousMap::new(
        grid_x,
        grid_y,
        grid_multipliers,
        None,
        None,
        None,
        InterpMethod::Linear,
    )?;

    println!("  Grid shape: 2 x 3 (6 points)");
    println!("  X range: 0-1000 m");
    println!("  Y range: 0-500 m");
    println!("  Speed range: 0.95 - 1.20\n");

    // ============================================================
    // SUMMARY
    // ============================================================
    println!("--- Summary ---\n");

    println!("HeterogeneousMap Key Points:");
    println!("  - Models spatial wind speed variations across a wind farm");
    println!("  - Can be defined per wind direction");
    println!("  - Speed multipliers: >1.0 = acceleration, <1.0 = deceleration");
    println!("  - Supports 2D (x, y) and 3D (x, y, z) grids");
    println!("  - Interpolates between grid points for smooth transitions\n");

    println!("Common Applications:");
    println!("  - Terrain effects (hills, ridges, valleys)");
    println!("  - Offshore coastal and stability effects");
    println!("  - Forest and vegetation effects");
    println!("  - Obstacle and sheltering effects\n");

    println!("Integration with FLORIS:");
    println!("  - Attach to WindRose via heterogeneous_map field");
    println!("  - Automatically applies multipliers during flow field calculation");
    println!("  - Affects both power and thrust calculations\n");

    println!("================================================");
    println!("Example completed successfully!");

    Ok(())
}
