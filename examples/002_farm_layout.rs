/// Example 8: Complex Wind Farm Layouts
///
/// This example demonstrates various wind farm layout configurations
/// and their impact on wake losses and power production.
///
/// Layouts covered:
/// 1. Linear array (baseline)
/// 2. Staggered/offset layout
/// 3. Grid pattern
/// 4. Cluster pattern
/// 5. Optimal spacing analysis

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 8: Complex Wind Farm Layouts");
    println!("==========================================\n");

    let d = 126.0; // NREL 5MW rotor diameter
    let rated_ws = 11.4; // Rated wind speed

    // ============================================================
    // Layout 1: Linear Array
    // ============================================================
    println!("--- Layout 1: Linear Array ---\n");

    let spacing = 5.0 * d;
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    println!("Linear array configuration:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, layout_y[i]);
    }
    println!("  Spacing: {:.0} m ({:.1}D)\n", spacing, spacing / d);

    let farm_linear = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;
    let flow_field = florus::core::FlowField::new(
        Array1::from_vec(vec![rated_ws]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm: farm_linear,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let powers = model.get_turbine_powers();
    let farm_power: f64 = powers.iter().sum();

    println!("Results:");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, powers[[0, ti]] / 1000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", farm_power / 1_000_000.0);

    // ============================================================
    // Layout 2: Staggered/Offset Layout
    // ============================================================
    println!("--- Layout 2: Staggered Layout ---\n");

    let stagger_offset = 3.0 * d; // Offset turbines
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0, stagger_offset, 0.0, stagger_offset]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    println!("Staggered configuration:");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, layout_y[i]);
    }
    println!("  Row spacing: {:.0} m", spacing);
    println!("  Lateral offset: {:.0} m\n", stagger_offset);

    let farm_staggered = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;
    let flow_field = florus::core::FlowField::new(
        Array1::from_vec(vec![rated_ws]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm: farm_staggered,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let powers = model.get_turbine_powers();
    let farm_power_staggered: f64 = powers.iter().sum();

    println!("Results:");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, powers[[0, ti]] / 1000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", farm_power_staggered / 1_000_000.0);

    // ============================================================
    // Layout 3: 3x3 Grid Pattern
    // ============================================================
    println!("--- Layout 3: 3x3 Grid Pattern ---\n");

    let grid_spacing = 7.0 * d;
    let mut layout_x = Vec::new();
    let mut layout_y = Vec::new();

    for row in 0..3 {
        for col in 0..3 {
            layout_x.push(col as f64 * grid_spacing);
            layout_y.push(row as f64 * grid_spacing);
        }
    }

    let turbine_types = vec!["nrel_5MW".to_string(); 9];

    println!("3x3 Grid configuration (9 turbines):");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, y);
    }
    println!("  Grid spacing: {:.0} m ({:.1}D)\n", grid_spacing, grid_spacing / d);

    let farm_grid = Farm::new(
        Array1::from_vec(layout_x.clone()),
        Array1::from_vec(layout_y.clone()),
        turbine_types.clone(),
    )?;
    let flow_field = florus::core::FlowField::new(
        Array1::from_vec(vec![rated_ws]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm: farm_grid,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let powers = model.get_turbine_powers();
    let farm_power_grid: f64 = powers.iter().sum();

    println!("Results:");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, powers[[0, ti]] / 1000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", farm_power_grid / 1_000_000.0);

    // ============================================================
    // Layout 4: Cluster Pattern
    // ============================================================
    println!("--- Layout 4: Cluster Pattern ---\n");

    let cluster_spacing = 3.0 * d; // Close within cluster
    let cluster_distance = 10.0 * d; // Distance between clusters
    let mut layout_x = Vec::new();
    let mut layout_y = Vec::new();

    // Cluster 1: turbines at origin
    layout_x.push(0.0);
    layout_y.push(0.0);
    layout_x.push(cluster_spacing);
    layout_y.push(0.0);
    layout_x.push(0.0);
    layout_y.push(cluster_spacing);

    // Cluster 2: downstream
    layout_x.push(cluster_distance);
    layout_y.push(0.0);
    layout_x.push(cluster_distance + cluster_spacing);
    layout_y.push(0.0);
    layout_x.push(cluster_distance);
    layout_y.push(cluster_spacing);

    let turbine_types = vec!["nrel_5MW".to_string(); 6];

    println!("Cluster configuration (2 clusters of 3):");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, y);
    }
    println!("  Within-cluster spacing: {:.0} m", cluster_spacing);
    println!("  Cluster spacing: {:.0} m\n", cluster_distance);

    let farm_cluster = Farm::new(
        Array1::from_vec(layout_x.clone()),
        Array1::from_vec(layout_y.clone()),
        turbine_types.clone(),
    )?;
    let flow_field = florus::core::FlowField::new(
        Array1::from_vec(vec![rated_ws]),
        Array1::from_vec(vec![270.0]),
        0.0, 0.14, 1.225,
        Array1::from_vec(vec![0.06]),
        90.0,
    )?;

    let mut model = florus::FlorisModel {
        farm: farm_cluster,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;

    let powers = model.get_turbine_powers();
    let farm_power_cluster: f64 = powers.iter().sum();

    println!("Results:");
    for ti in 0..model.farm.n_turbines() {
        println!("  Turbine {}: {:.1} kW", ti, powers[[0, ti]] / 1000.0);
    }
    println!("  Total Farm Power: {:.2} MW\n", farm_power_cluster / 1_000_000.0);

    // ============================================================
    // Spacing Analysis
    // ============================================================
    println!("--- Spacing Analysis ---\n");

    println!("Testing different spacing ratios:");
    println!("  {:>8}  {:>12}  {:>12}", "Spacing", "Farm Power", "Wake Loss");
    println!("  {}", "-".repeat(40));

    let spacing_ratios = vec![3.0, 5.0, 7.0, 10.0, 15.0];

    let mut best_spacing = 0.0;
    let mut best_power = 0.0;

    for &ratio in &spacing_ratios {
        let spacing_test = ratio * d;
        let layout_x_test = Array1::from_vec(vec![0.0, spacing_test, 2.0 * spacing_test]);
        let layout_y_test = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let turbine_types_test = vec!["nrel_5MW".to_string(); 3];

        let farm_test = Farm::new(layout_x_test, layout_y_test, turbine_types_test)?;
        let flow_field_test = florus::core::FlowField::new(
            Array1::from_vec(vec![rated_ws]),
            Array1::from_vec(vec![270.0]),
            0.0, 0.14, 1.225,
            Array1::from_vec(vec![0.06]),
            90.0,
        )?;

        let mut model_test = florus::FlorisModel {
            farm: farm_test,
            flow_field: flow_field_test,
            state: florus::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };

        model_test.initialize_grid()?;
        model_test.initialize_flow_field()?;
        model_test.run()?;

        let powers_test = model_test.get_turbine_powers();
        let farm_power_test: f64 = powers_test.iter().sum();

        // Calculate wake loss (compare to no-wake: 3 turbines × rated power)
        let no_wake_power = 3.0 * 5_000_000.0; // 3 turbines × 5 MW
        let wake_loss = (1.0 - farm_power_test / no_wake_power) * 100.0;

        println!("  {:>6.0}D  {:>10.2} MW  {:>10.1}%", ratio, farm_power_test / 1_000_000.0, wake_loss);

        if farm_power_test > best_power {
            best_power = farm_power_test;
            best_spacing = ratio;
        }
    }

    println!("\n  Best spacing ratio: {:.0}D\n", best_spacing);

    // ============================================================
    // Summary
    // ============================================================
    println!("--- Summary ---\n");

    println!("Layout Comparison (9 m/s wind, 270° direction):");
    println!("  {:>20}  {:>12}", "Layout", "Farm Power");
    println!("  {}", "-".repeat(35));
    println!("  {:>20}  {:>10.2} MW", "Linear (4 turbines)", farm_power / 1_000_000.0);
    println!("  {:>20}  {:>10.2} MW", "Staggered (4 turbines)", farm_power_staggered / 1_000_000.0);
    println!("  {:>20}  {:>10.2} MW", "3x3 Grid (9 turbines)", farm_power_grid / 1_000_000.0);
    println!("  {:>20}  {:>10.2} MW", "Cluster (6 turbines)", farm_power_cluster / 1_000_000.0);

    println!("\nLayout Design Considerations:");
    println!("  1. Staggered layouts reduce wake interactions");
    println!("  2. Optimal spacing balances land use vs wake losses");
    println!("  3. Grid patterns efficient for large farms");
    println!("  4. Cluster layouts useful for space-constrained sites");
    println!("  5. Consider predominant wind direction in layout design");

    println!("\n==========================================");
    println!("Example completed successfully!");

    Ok(())
}
