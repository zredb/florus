use florus::Array1;

/// Example 1: Opening FLORIS and Computing Power
///
/// This example demonstrates the key concepts in FLORIS-RS:
///
/// 1. Initializing a FLORIS model
/// 2. Changing the wind farm layout
/// 3. Changing the wind speed, wind direction and turbulence intensity
/// 4. Running the FLORIS simulation
/// 5. Getting the power output of the turbines
///
/// This is the Rust equivalent of Python's 001_opening_floris_computing_power.py

fn main() -> anyhow::Result<()> {
    // The FlorisModel class is the entry point for most usage.
    // Initialize using an input yaml file or create programmatically

    // In this example, we'll create the model programmatically
    // (equivalent to Python's fmodel = FlorisModel("inputs/gch.yaml"))

    println!("\n--- Setting Wind Conditions ---");

    // Create the FlorisModel
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();

    model.set_layout(
        &Array1::from_vec(vec![0.0, 500.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0, 8.0, 10.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06]),
    )?;

    model.run()?;

    // ============================================================
    // Getting the power output
    // ============================================================
    // In Python:
    // turbine_powers = fmodel.get_turbine_powers() / 1000.0
    // farm_power = fmodel.get_farm_power() / 1000.0

    let turbine_powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();

    println!("\n--- Results ---");

    // The turbine power matrix has dimensions (n_findex, n_turbines)
    println!(
        "Turbine power matrix shape: ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );
    println!("Turbine powers (kW):");
    for ti in 0..model.farm.n_turbines() {
        println!(
            "  Turbine {}: {:.1} kW",
            ti,
            turbine_powers[[0, ti]] / 1000.0
        );
    }

    // Farm power is a 1D array of length n_findex
    println!("\nFarm power: {:.1} kW", farm_power[[0]] / 1000.0);

    // ============================================================
    // Demonstrating multiple conditions (n_findex = 4)
    // ============================================================
    println!("\n--- Multiple Conditions (n_findex = 4) ---");

    // In Python:
    // fmodel.set(
    //     wind_directions=np.array([270.0, 270.0, 270.0, 270.0]),
    //     wind_speeds=[8.0, 8.0, 10.0, 10.0],
    //     turbulence_intensities=np.array([0.06, 0.06, 0.06, 0.06])
    // )

    // Set up 4 conditions
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0, 8.0, 10.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06]),
    )?;

    println!("Simulating 4 conditions:");
    for i in 0..4 {
        let ws = if i < 2 { 8.0 } else { 10.0 };
        println!("  Condition {}: {} m/s, 270°", i + 1, ws);
    }

    // Run the simulation
    model.run()?;

    // Get results
    let turbine_powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();

    println!("\nResults for 4 conditions:");
    println!(
        "Turbine power matrix shape: ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );

    println!("\nTurbine powers (kW):");
    for fi in 0..4 {
        print!("  Condition {}: ", fi + 1);
        for ti in 0..model.farm.n_turbines() {
            print!("T{}: {:.0}  ", ti, turbine_powers[[fi, ti]] / 1000.0);
        }
        println!();
    }

    println!("\nFarm power for each condition (kW):");
    for fi in 0..4 {
        println!(
            "  Condition {}: {:.1} kW",
            fi + 1,
            farm_power[[fi]] / 1000.0
        );
    }

    println!("\n===========================================================");
    println!("Example completed successfully!");

    Ok(())
}
