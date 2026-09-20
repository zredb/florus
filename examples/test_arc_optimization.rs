/// Test to verify TurbineType Arc optimization
use florus::core::Farm;
use std::sync::Arc;

fn main() {
    println!("=== Testing TurbineType Arc Optimization ===\n");

    // Create a farm with multiple turbines of the same type
    let layout_x = vec![0.0, 500.0, 1000.0];
    let layout_y = vec![0.0, 0.0, 0.0];

    let farm = Farm::new(
        ndarray::Array1::from_vec(layout_x),
        ndarray::Array1::from_vec(layout_y),
        &vec!["nrel_5MW".to_string()], // All turbines use the same type
    )
    .expect("Failed to create farm");

    println!("Created farm with {} turbines", farm.n_turbines());

    // Verify that all turbines share the same TurbineType via Arc
    if farm.turbines.len() >= 2 {
        let t1 = &farm.turbines[0];
        let t2 = &farm.turbines[1];

        // Check if they share the same Arc instance
        if Arc::ptr_eq(&t1.turbine_type, &t2.turbine_type) {
            println!("✅ SUCCESS: Turbines 0 and 1 share the same TurbineType (Arc)");
        } else {
            println!("❌ FAIL: Turbines have different TurbineType instances");
        }

        // Verify data is accessible
        println!("\nTurbine 0 configuration:");
        println!("  - Rotor diameter: {} m", t1.turbine_type.rotor_diameter);
        println!("  - Hub height: {} m", t1.turbine_type.hub_height);
        println!("  - Name: {}", t1.turbine_type.name);

        println!("\nTurbine 1 configuration:");
        println!("  - Rotor diameter: {} m", t2.turbine_type.rotor_diameter);
        println!("  - Hub height: {} m", t2.turbine_type.hub_height);
        println!("  - Name: {}", t2.turbine_type.name);

        // Verify Arc reference count
        let ref_count = Arc::strong_count(&t1.turbine_type);
        println!("\nArc strong count for TurbineType: {}", ref_count);
        println!(
            "  (Expected: {} for {} turbines)",
            farm.turbines.len(),
            farm.turbines.len()
        );

        if ref_count == farm.turbines.len() as usize {
            println!("✅ Reference count matches number of turbines!");
        } else {
            println!("⚠️  Warning: Reference count mismatch");
        }
    }

    println!("\n=== Test Complete ===");
}
