use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 17: Offshore Wind Considerations");
    println!("====================================\n");

    let d = 126.0;
    let turbine_types = vec!["nrel_5MW".to_string()];

    println!("Offshore Wind Farm Considerations:");
    println!("  - Different atmospheric stability");
    println!("  - Wake recovery patterns");
    println!("  - Layout optimization strategies\n");

    println!("--- Key Offshore Differences ---\n");

    println!("1. ATMOSPHERIC STABILITY:");
    println!("   - Offshore: More stable (lower TI at high wind speeds)");
    println!("   - Onshore: More unstable (higher turbulence)");
    println!("   - Affects wake recovery and mixing\n");

    println!("2. WIND CONDITIONS:");
    println!("   - Higher average wind speeds offshore");
    println!("   - Lower turbulence intensity");
    println!("   - More consistent wind direction\n");

    println!("3. LAYOUT CONSIDERATIONS:");
    println!("   - Typical spacing: 7-10D for offshore");
    println!("   - Greater spacing due to larger wakes");
    println!("   - Wake steering more effective\n");

    println!("--- Offshore Layout Example ---\n");

    let spacing = 8.0 * d;
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 4]);

    println!("4-turbine offshore layout:");
    println!("  Spacing: {:.0} m ({:.1}D)", spacing, spacing / d);

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    println!("  Created successfully!\n");

    println!("--- Offshore TI Model ---\n");

    println!("Turbulence Intensity Models:");
    println!("  1. IEC: International Electrotechnical Commission");
    println!("  2. KTI: Kaimal turbulence model");
    println!("  3. NTI: Non-turbulent intensity\n");

    println!("Offshore typically uses:");
    println!("  Lower TI (0.06-0.08) vs onshore (0.08-0.12)");
    println!("  Less variation with wind speed");
    println!("  More stable conditions\n");

    println!("--- Wake Recovery ---\n");

    println!("Offshore Wake Behavior:");
    println!("  - Slower wake recovery due to stable conditions");
    println!("  - Wakes travel farther");
    println!("  - Greater wake losses\n");

    println!("Recovery distance factors:");
    println!("  Offshore: 15-30 rotor diameters");
    println!("  Onshore: 10-20 rotor diameters\n");

    println!("====================================");
    println!("Example completed successfully!");
    Ok(())
}
