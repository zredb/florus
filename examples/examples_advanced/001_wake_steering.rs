use florus::core::Farm;
use florus::types::Array1;
use florus::optimization::{yaw_cosine_loss, estimate_wake_deflection_angle, YawAngleBounds};

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 16: Wake Steering Deep Dive");
    println!("======================================\n");

    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d, 15.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 4]);
    let turbine_types = vec!["nrel_5MW".to_string(); 4];

    println!("Wake Steering Demonstration:");
    println!("  Farm: 4 turbines at 5D spacing");
    println!("  Wind: 8 m/s from West (270 deg)\n");

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("--- Wake Deflection Physics ---\n");

    println!("Wake Steering Fundamentals:");
    println!("  When a turbine yaws, its wake is deflected downwind");
    println!("  The deflection depends on yaw angle and thrust coefficient\n");

    let ct_values = [0.6, 0.7, 0.8, 0.9];

    println!("Deflection vs Thrust Coefficient:");
    println!("{:>8} {:>10} {:>10}", "Ct", "20 deg", "30 deg");
    println!("{}", "-".repeat(30));

    for &ct in &ct_values {
        let defl_20 = estimate_wake_deflection_angle(20.0_f64.to_radians(), ct, d, 5.0 * d, 0.1, 1.0);
        let defl_30 = estimate_wake_deflection_angle(30.0_f64.to_radians(), ct, d, 5.0 * d, 0.1, 1.0);
        println!("{:>8.1} {:>10.2} m {:>10.2} m", ct, defl_20, defl_30);
    }

    println!("\n--- Cosine Loss ---\n");

    println!("Yaw Angle vs Power Factor:");
    println!("{:>8} {:>12} {:>12}", "Yaw", "cos(yaw)", "cos^3(yaw)");
    println!("{}", "-".repeat(34));

    for yaw in [0, 5, 10, 15, 20, 25, 30] {
        let yaw_rad = yaw as f64;
        let factor = yaw_cosine_loss(yaw_rad.to_radians(), 3.0);
        println!("{:>8} {:>12.4} {:>12.4}", yaw, yaw_rad.to_radians().cos(), factor);
    }

    println!("\n--- Yaw Optimization ---\n");

    let bounds = YawAngleBounds::new(-30.0, 30.0);
    println!("Yaw Angle Bounds:");
    println!("  Minimum: {} deg", bounds.min_yaw);
    println!("  Maximum: {} deg", bounds.max_yaw);

    println!("\nOptimization Strategies:");
    println!("  1. SERIAL REFINE: Turbine-by-turbine optimization");
    println!("  2. GEOMETRIC: Fast approximation based on layout");
    println!("  3. SCIPY: Gradient-based optimization");

    println!("\n--- Trade-off Analysis ---\n");

    println!("Wake Steering Trade-off:");
    println!("  Upstream turbine loses power due to yaw misalignment");
    println!("  Downstream turbines gain from reduced wake");
    println!("  Net gain when deflection benefit > cosine loss\n");

    let upstream_loss = 1.0 - yaw_cosine_loss(20.0_f64.to_radians(), 3.0);
    let wake_recovery = 0.15;
    println!("Example at 20 deg yaw:");
    println!("  Upstream power loss: {:.1}%", upstream_loss * 100.0);
    println!("  Typical wake recovery: {:.0}%", wake_recovery * 100.0);
    println!("  Net effect: Depends on turbine spacing and layout");

    println!("\n--- Practical Guidelines ---\n");

    println!("Wake Steering Best Practices:");
    println!("  1. Yaw angles of 15-25 deg typically optimal");
    println!("  2. More effective in offshore (stable atmospheres)");
    println!("  3. Limited benefit at high turbulence sites");
    println!("  4. Consider fatigue loads from frequent yaw changes");
    println!("  5. Optimal yaw varies with wind direction");

    println!("\n======================================");
    println!("Example completed successfully!");
    Ok(())
}
