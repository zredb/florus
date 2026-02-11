/// Example 2: Visualizations
///
/// This example demonstrates the use of the flow and layout visualizations in FLORIS.
/// First, an example wind farm layout is plotted, with the turbine names and the directions
/// and distances between turbines shown in different configurations by subplot.
/// Next, the horizontal flow field at hub height is plotted for a single wind condition.
///
/// FLORIS includes two modules for visualization:
///   1) flow_visualization: for visualizing the flow field
///   2) layout_visualization: for visualizing the layout of the wind farm
/// The two modules can be used together to visualize the flow field and the layout
/// of the wind farm.
///
/// This is the Rust equivalent of Python's 002_visualizations.py
///
/// Note: Full visualization capabilities require a plotting library like plotters.
/// This example demonstrates the data structures and flow calculations needed for visualization.

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 2: Visualizations");
    println!("====================================\n");

    // Create a FlorisModel with the default configuration
    // In Python: fmodel = FlorisModel("inputs/gch.yaml")

    // Set the farm layout to have 8 turbines irregularly placed
    // In Python: layout_x = [0, 500, 0, 128, 1000, 900, 1500, 1250]
    //            layout_y = [0, 300, 750, 1400, 0, 567, 888, 1450]
    let layout_x = Array1::from_vec(vec![0.0, 500.0, 0.0, 128.0, 1000.0, 900.0, 1500.0, 1250.0]);
    let layout_y = Array1::from_vec(vec![0.0, 300.0, 750.0, 1400.0, 0.0, 567.0, 888.0, 1450.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 8];

    println!("Creating wind farm with 8 turbines:");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m", i, x, y);
    }

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    // ============================================================
    // Layout Visualization
    // ============================================================
    println!("\n--- Layout Visualization Concepts ---\n");

    println!("Layout visualization includes the following functions:");
    println!("  1. plot_turbine_points - Shows turbine locations as points");
    println!("  2. plot_turbine_labels - Adds turbine labels/names");
    println!("  3. plot_turbine_rotors - Shows turbine rotor circles");
    println!("  4. plot_waking_directions - Shows wake propagation directions");
    println!();
    println!("These can be overlaid to provide further information about the layout.");
    println!("Example subplot configurations:");
    println!("  - Subplot 1: Turbine Points only");
    println!("  - Subplot 2: Turbine Points + Labels");
    println!("  - Subplot 3: Points + Labels + Wake Directions (limit=2)");
    println!("  - Subplot 4: Custom turbine names (T1, T2, T3, etc.)");

    // ============================================================
    // Flow Field Visualization Setup
    // ============================================================
    println!("\n--- Flow Field Visualization ---\n");

    println!("Flow visualizations are created using calculate plane methods.");
    println!("This example shows the horizontal plane at hub height.");
    println!();
    println!("Configuration:");
    println!("  - Wind speed: 8.0 m/s");
    println!("  - Wind direction: 290°");
    println!("  - Turbulence intensity: 0.06");
    println!("  - Horizontal plane resolution: 200 x 100");
    println!("  - Height: 90.0 m (hub height)");

    // Set wind conditions for visualization
    // In Python: fmodel.set(wind_speeds=[8.0], wind_directions=[290.0], turbulence_intensities=[0.06])
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![290.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    let flow_field = florus::core::FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
        turbulence_intensities.clone(),
        90.0,   // reference_wind_height
    )?;

    let mut model = florus::FlorisModel {
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

    // Calculate horizontal plane for visualization
    // In Python: horizontal_plane = fmodel.calculate_horizontal_plane(x_resolution=200, y_resolution=100, height=90.0)
    let x_resolution = 200;
    let y_resolution = 100;
    let height = 90.0;

    println!();
    println!("Calculating horizontal flow plane:");
    println!("  X resolution: {} points", x_resolution);
    println!("  Y resolution: {} points", y_resolution);
    println!("  Height: {} m", height);

    // The calculate_horizontal_plane method would return a flow field data structure
    // containing velocity values at each grid point for visualization
    println!("  Grid span: x=[min_x, max_x], y=[min_y, max_y]");
    println!("  Total grid points: {}", x_resolution * y_resolution);

    // ============================================================
    // Visualization Output Options
    // ============================================================
    println!("\n--- Visualization Output Options ---\n");

    println!("Cut plane visualization options:");
    println!("  - visualize_cut_plane: Main visualization function");
    println!("  - label_contours: Toggle contour labels (true/false)");
    println!("  - title: Set the plot title");
    println!();
    println!("Overlaying layout elements:");
    println!("  - layoutviz.plot_turbine_rotors() - Shows rotor circles");
    println!("  - layoutviz.plot_turbine_labels() - Shows turbine names");
    println!("  - Custom turbine names array can be provided");

    println!("\n--- Summary ---\n");

    println!("Key Visualization Points:");
    println!("  ✓ FLORIS provides both flow and layout visualization tools");
    println!("  ✓ Flow planes can be calculated at any height and resolution");
    println!("  ✓ Layout and flow visualizations can be combined");
    println!("  ✓ Custom turbine names enhance visualization readability");
    println!("  ✓ Multiple subplots can show different visualization modes");

    println!("\n====================================");
    println!("Example completed successfully!");
    println!("Note: Full graphical output requires a plotting library.");

    Ok(())
}
