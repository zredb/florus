use florus::core::Farm;
use florus::types::{Array1};
use florus::wind_data::WindTIRose;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 18: Turbulence Intensity Rose");
    println!("==================================\n");

    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 3]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    println!("Turbulence Intensity Rose Analysis:");
    println!("  Farm: 3 turbines at 5D spacing\n");

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    println!("--- TI Rose Configuration ---\n");

    let wind_directions = Array1::from_vec(vec![
        0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0,
    ]);

    let wind_speeds = Array1::from_vec(vec![5.0, 7.5, 10.0, 12.5, 15.0]);

    // Define turbulence intensities bins (3 bins for this example)
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.08, 0.10]);

    // Create 3D TI table: [n_directions, n_speeds, n_tis]
    let ti_table_offshore = Array3::from_shape_fn((8, 5, 3), |(i, j, k)| {
        let base_ti = match i {
            0..=1 => 0.06,
            2..=3 => 0.07,
            4..=5 => 0.09,
            _ => 0.07,
        };
        base_ti + (j as f64 * 0.01) + (k as f64 * 0.02)
    });

    // Print wind directions with average TI
    println!("Wind Directions:");
    for i in 0..wind_directions.len() {
        let avg_ti: f64 = (0..wind_speeds.len()).map(|j| {
            // Average over TI bins
            let ti_sum: f64 = (0..turbulence_intensities.len())
                .map(|k| ti_table_offshore[[i, j, k]])
                .sum();
            ti_sum / turbulence_intensities.len() as f64
        }).sum::<f64>() / wind_speeds.len() as f64;
        println!("  {:>4.0} deg: avg TI = {:.2}", wind_directions[i], avg_ti);
    }

    println!("\n--- WindTIRose ---\n");

    let ti_rose = WindTIRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        turbulence_intensities,
        ti_table_offshore,
        None,
        None,
    )?;

    println!("WindTIRose Configuration:");
    println!("  Directions: {}", ti_rose.wind_directions.len());
    println!("  Speed bins: {}", ti_rose.wind_speeds.len());
    println!("  TI bins: {}", ti_rose.turbulence_intensities.len());

    println!("\n--- TI Effects on Wake ---\n");

    println!("Turbulence Intensity Impact:");
    println!("  High TI (>0.10):");
    println!("    - Faster wake recovery");
    println!("    - Less wake loss");
    println!("    - Higher turbine loads\n");

    println!("  Low TI (<0.06):");
    println!("    - Slower wake recovery");
    println!("    - Greater wake losses");
    println!("    - Lower turbine loads\n");

    println!("--- TI by Wind Sector ---\n");

    println!("{:>8} {:>10} {:>15}", "Dir", "TI @ 8m/s", "TI @ 12m/s");
    println!("{}", "-".repeat(35));

    for i in 0..wind_directions.len() {
        // Get TI at index 1 (7.5 m/s) and index 2 (10 m/s) for middle TI bin
        let ti_8 = ti_rose.ti_table[[i, 1, 1]];  // 7.5 m/s, middle TI bin
        let ti_12 = ti_rose.ti_table[[i, 2, 1]]; // 10 m/s, middle TI bin
        println!("{:>8.0} deg {:>10.2} {:>15.2}", wind_directions[i], ti_8, ti_12);
    }

    println!("\n--- TI Classification ---\n");

    println!("IEC TI Categories:");
    println!("  Class I: High wind, low TI");
    println!("  Class II: Medium conditions");
    println!("  Class III: Low wind, high TI\n");

    println!("Site Assessment:");
    println!("  Offshore: TI typically 0.06-0.08");
    println!("  Onshore: TI typically 0.08-0.12");
    println!("  Complex terrain: TI > 0.15 possible\n");

    println!("==================================");
    println!("Example completed successfully!");
    Ok(())
}
