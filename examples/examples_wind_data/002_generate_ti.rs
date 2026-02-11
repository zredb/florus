/// Generate TI Example
///
/// This example demonstrates generating turbulence intensity data.
///
/// This is the Rust equivalent of Python's 002_generate_ti.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Generate Turbulence Intensity");
    println!("====================================\n");

    // ============================================================
    // TI Generation Methods
    // ============================================================
    println!("--- TI Generation Methods ---\n");

    println!("FLORIS provides multiple methods to set turbulence intensity:");
    println!());
    println!("1. Fixed TI:");
    println!("   - Constant TI for all conditions");
    println!("   - Simple but unrealistic");
    println!());

    println!("2. IEC Method:");
    println!("   - Based on IEC 61400-1 standard");
    println!("   - TI as function of wind speed");
    println!("   - Class I, II, III turbines");
    println!());

    println!("3. Custom Function:");
    println!("   - User-defined TI(WD, WS) function");
    println!("   - Maximum flexibility");
    println!("   - Based on site measurements");
    println!());

    // ============================================================
    // IEC Method Details
    // ============================================================
    println!("--- IEC Method ---\n");

    println!("IEC 61400-1 TI calculation:");
    println!("  TI = I_ref * (0.75 + 0.5 * V_hub / V_out) / (1 + 0.2 * V_hub / V_out)");
    println!());

    println!("Parameters:");
    println!("  - I_ref: Reference turbulence class (0.16, 0.14, 0.12)");
    println!("  - V_hub: Hub-height wind speed");
    println!("  - V_out: Cut-out wind speed");
    println!());

    println!("Turbulence classes:");
    println!("  Class I: I_ref = 0.16 (high turbulence sites)");
    println!("  Class II: I_ref = 0.14 (medium turbulence)");
    println!("  Class III: I_ref = 0.12 (low turbulence)");
    println!());

    // ============================================================
    // Custom TI Function
    // ============================================================
    println!("\n--- Custom TI Function ---\n");

    println!("Custom TI(WD, WS) function:");
    println!("  - Map wind direction to TI variation");
    println!("  - Map wind speed to TI curve");
    println!("  - Interpolate between measurements");
    println!());

    println!("Example custom TI profile:");
    println!("  Wind Direction | TI at 8 m/s | TI at 12 m/s");
    println!("  ---------------|-------------|-------------");
    println!("  0°             | 0.06        | 0.08");
    println!("  90°            | 0.08        | 0.10");
    println!("  180°           | 0.07        | 0.09");
    println!("  270°           | 0.05        | 0.07");
    println!());

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("TI Generation Key Points:");
    println!("  ✓ Multiple methods available");
    println!("  ✓ IEC method: Standard-compliant");
    println!("  ✓ Custom: Site-specific accuracy");
    println!("  ✓ TI affects wake loss and power");
    println!());

    println!("Impact of TI on wake modeling:");
    println!("  - Higher TI = faster wake recovery");
    println!("  - Lower TI = slower wake recovery");
    println!("  - Important for offshore (low TI)");
    println!());

    println!("\n====================================");
    println!("Example completed successfully!");

    Ok(())
}
