use florus::core::Farm;
use florus::floris_config::SolverConfig;
use florus::types::Array1;
use florus::optimization::layout_optimization::{Boundary, LayoutOptimizationBoundaryGrid, LayoutOptimizationConfig};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 15: Layout Optimization");
    println!("==================================\n");

    let d = 126.0;

    println!("Layout Optimization Setup:");
    println!("  Turbine type: nrel_5MW");
    println!("  Rotor diameter: {:.0} m", d);

    // Create a simple farm layout
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 3]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];
    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    // Create a flow field
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06]);
    let flow_field = florus::core::FlowField::new(
        wind_speeds,
        wind_directions,
        0.0,
        0.14,
        1.225,
        turbulence_intensities,
        90.0,
    )?;

    // Create the FlorisModel
    let model = florus::FlorisModel {
        farm: farm.clone(),
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver: SolverConfig::default(),
        model_manager: None,
    };

    let boundary = Boundary::Rectangle {
        min_x: 0.0,
        max_x: 2000.0,
        min_y: 0.0,
        max_y: 1000.0,
    };

    println!("\nOptimization Boundary:");
    println!("  X range: {:.0} - {:.0} m", 0.0, 2000.0);
    println!("  Y range: {:.0} - {:.0} m", 0.0, 1000.0);

    println!("\n--- Optimization Approaches ---\n");

    println!("1. BOUNDARY GRID OPTIMIZATION:");
    println!("   - Places turbines on a grid within boundary");
    println!("   - Good for initial layout exploration");
    println!("   - Adjustable grid resolution");

    // Create boundary grid optimizer with builder pattern
    let boundary_optimizer = LayoutOptimizationBoundaryGrid::new(
        &model,
        boundary.clone(),
        20, // grid_resolution
    )?.with_min_dist(5.0 * d);

    println!("   Optimizer created with boundary constraints");
    println!("   Grid resolution: 20x20 points");
    println!("   Minimum spacing: {:.0} m (5D)", 5.0 * d);

    println!("\n2. LAYOUT OPTIMIZATION STRATEGIES:");
    println!("   - Grid-based: Fast exploration on defined grid");
    println!("   - Mixed-integer: Combines grid and continuous optimization");
    println!("   - Supports min distance constraints");
    println!("   - Can use value-based optimization (AVP) or energy (AEP)");

    println!("\n--- Configuration ---\n");

    println!("LayoutOptimizationConfig:");
    println!("  min_dist: Minimum inter-turbine spacing (use Some(value))");
    println!("  boundaries: Boundary definition");
    println!("  use_value: true for AVP optimization, false for AEP");

    println!("\n--- Expected Workflow ---\n");

    println!("1. DEFINE BOUNDARY:");
    println!("   Boundary::Rectangle {{ min_x, max_x, min_y, max_y }}");
    println!("   Boundary::Polygon(vec![(x1,y1), (x2,y2), ...])");

    println!("\n2. CREATE FLORIS MODEL:");
    println!("   let model = FlorisModel {{ ... }};");

    println!("\n3. CREATE OPTIMIZER:");
    println!("   LayoutOptimizationBoundaryGrid::new(&model, boundary, resolution)?");
    println!("       .with_min_dist(5.0 * d)");

    println!("\n4. RUN OPTIMIZATION:");
    println!("   let result = optimizer.optimize()?;");

    println!("\n--- Optimization Tips ---\n");

    println!("1. Start with coarse grid, refine as needed");
    println!("2. Consider prevailing wind direction");
    println!("3. Balance spacing vs. energy production");
    println!("4. Include constraints: noise, wildlife, etc.");
    println!("5. Validate against local regulations");

    println!("\n==================================");
    println!("Example completed successfully!");
    Ok(())
}
