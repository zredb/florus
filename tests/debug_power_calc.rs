mod common;

use florus::floris_config::{
    CombinationModelConfig, DeflectionModelConfig, FarmConfig, FlorisConfig, FlowFieldConfig,
    LoggingConfig, SolverConfig, TurbulenceModelConfig, VelocityModelConfig, WakeConfig,
    WakeModelStringsConfig,
};
use florus::floris_model::FlorisModel;
use florus::types::Array2;
use ndarray::Array;

const WIND_DIRECTIONS: [f64; 16] = [
    270.0, 270.0, 270.0, 270.0, 360.0, 360.0, 360.0, 360.0, 285.0, 285.0, 285.0, 285.0, 315.0,
    315.0, 315.0, 315.0,
];

const WIND_SPEEDS: [f64; 16] = [
    8.0, 9.0, 10.0, 11.0, 8.0, 9.0, 10.0, 11.0, 8.0, 9.0, 10.0, 11.0, 8.0, 9.0, 10.0, 11.0,
];

const TURBULENCE_INTENSITIES: [f64; 16] = [
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
];

const ROTOR_DIAMETER: f64 = 126.0;
const X_COORDS: [f64; 3] = [0.0, 5.0 * ROTOR_DIAMETER, 10.0 * ROTOR_DIAMETER];
const Y_COORDS: [f64; 3] = [0.0, 0.0, 0.0];

fn create_model(velocity_model: VelocityModelConfig) -> FlorisModel {
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    let flow_field = FlowFieldConfig {
        wind_speeds: WIND_SPEEDS.to_vec(),
        wind_directions: WIND_DIRECTIONS.to_vec(),
        turbulence_intensities: TURBULENCE_INTENSITIES.to_vec(),
        wind_shear: 0.12,
        wind_veer: 0.0,
        air_density: 1.225,
        reference_wind_height: 90.0,
        multidim_conditions: None,
    };

    let config = FlorisConfig {
        name: "debug_test".to_string(),
        description: Some("Debug power calculation".to_string()),
        floris_version: "v4".to_string(),
        logging: LoggingConfig::default(),
        solver: SolverConfig::default(),
        farm: FarmConfig {
            layout_x: X_COORDS.to_vec(),
            layout_y: Y_COORDS.to_vec(),
            turbine_type: turbine_types,
        },
        flow_field,
        wake: WakeConfig {
            model_strings: WakeModelStringsConfig {
                velocity_model,
                deflection_model: DeflectionModelConfig::Jimenez {
                    ad: 0.0,
                    bd: 0.0,
                    kd: 0.05,
                },
                combination_model: CombinationModelConfig::SOSFS,
                turbulence_model: TurbulenceModelConfig::None,
            },
            enable_secondary_steering: true,
            enable_yaw_added_recovery: true,
            enable_transverse_velocities: true,
            enable_active_wake_mixing: false,
            enable_wake_mixing: false,
            use_parallel_calc: false,
        },
        turbine_library: "turbine_library".to_string(),
    };

    FlorisModel::from_config(config).expect("Failed to create FlorisModel")
}

#[test]
fn debug_power_calculation() {
    let mut model = create_model(VelocityModelConfig::Jensen { we: 0.05 });

    println!("=== BEFORE RUN ===");
    let ff = model.flow_field();
    println!("u shape: {:?}", ff.u.shape());
    println!("u_sorted shape: {:?}", ff.u_sorted.shape());
    assert!(ff.u.is_empty(), "u should be empty before run()");

    println!("\n=== RUNNING MODEL ===");
    model.run().expect("Model run failed");

    println!("\n=== AFTER RUN ===");
    let ff = model.flow_field();
    println!("u shape: {:?}", ff.u.shape());
    println!("u_sorted shape: {:?}", ff.u_sorted.shape());

    // run() must populate the flow field: get_turbine_powers() reads u
    assert!(!ff.u.is_empty(), "u should be populated after run()");

    // Print u values for first turbine, first findex
    println!("\nu[fi=0, ti=0, :, :]:");
    for iy in 0..ff.u.shape()[2] {
        for iz in 0..ff.u.shape()[3] {
            println!("  u[0,0,{},{}] = {:.6}", iy, iz, ff.u[[0, 0, iy, iz]]);
        }
    }

    println!("\nu_sorted[fi=0, ti=0, :, :]:");
    for iy in 0..ff.u_sorted.shape()[2] {
        for iz in 0..ff.u_sorted.shape()[3] {
            println!(
                "  u_sorted[0,0,{},{}] = {:.6}",
                iy,
                iz,
                ff.u_sorted[[0, 0, iy, iz]]
            );
        }
    }

    // Print u values for all turbines, first findex
    println!("\nu[fi=0, all turbines, 0, 0]:");
    for ti in 0..ff.u.shape()[1] {
        println!("  u[0,{},0,0] = {:.6}", ti, ff.u[[0, ti, 0, 0]]);
    }

    println!("\nu_sorted[fi=0, all turbines, 0, 0]:");
    for ti in 0..ff.u_sorted.shape()[1] {
        println!(
            "  u_sorted[0,{},0,0] = {:.6}",
            ti,
            ff.u_sorted[[0, ti, 0, 0]]
        );
    }

    println!("\n=== POWER CALCULATION ===");
    let powers = model.get_turbine_powers();
    println!("powers shape: {:?}", powers.shape());
    println!("powers[fi=0]: {:?}", powers.row(0).to_vec());
    println!("powers[fi=1]: {:?}", powers.row(1).to_vec());

    // Also check the flow field u is actually getting values
    println!("\n=== U FIELD SUMMARY ===");
    let u_min = ff.u.iter().cloned().fold(f64::INFINITY, f64::min);
    let u_max = ff.u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let u_mean: f64 = ff.u.iter().sum::<f64>() / ff.u.len() as f64;
    println!("u: min={:.6}, max={:.6}, mean={:.6}", u_min, u_max, u_mean);

    let u_sorted_min = ff.u_sorted.iter().cloned().fold(f64::INFINITY, f64::min);
    let u_sorted_max = ff
        .u_sorted
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let u_sorted_mean: f64 = ff.u_sorted.iter().sum::<f64>() / ff.u_sorted.len() as f64;
    println!(
        "u_sorted: min={:.6}, max={:.6}, mean={:.6}",
        u_sorted_min, u_sorted_max, u_sorted_mean
    );
}

#[test]
fn debug_direct_solver() {
    // Compare with the manual solver approach that works in regression tests
    use florus::core::wake::WakeModelManager;
    use florus::core::TurbineGrid;
    use florus::floris_config::{SolverType, WakeModelStringsConfig};
    use florus::types::Array1;

    let turbine_types = vec!["nrel_5MW".to_string(); 3];
    let flow_field_config = FlowFieldConfig {
        wind_speeds: WIND_SPEEDS.to_vec(),
        wind_directions: WIND_DIRECTIONS.to_vec(),
        turbulence_intensities: TURBULENCE_INTENSITIES.to_vec(),
        wind_shear: 0.12,
        wind_veer: 0.0,
        air_density: 1.225,
        reference_wind_height: 90.0,
        multidim_conditions: None,
    };

    let config = FlorisConfig {
        name: "debug_direct".to_string(),
        description: Some("Debug direct solver".to_string()),
        floris_version: "v4".to_string(),
        logging: LoggingConfig::default(),
        solver: SolverConfig::default(),
        farm: FarmConfig {
            layout_x: X_COORDS.to_vec(),
            layout_y: Y_COORDS.to_vec(),
            turbine_type: turbine_types,
        },
        flow_field: flow_field_config,
        wake: WakeConfig {
            model_strings: WakeModelStringsConfig {
                velocity_model: VelocityModelConfig::Jensen { we: 0.05 },
                deflection_model: DeflectionModelConfig::Jimenez {
                    ad: 0.0,
                    bd: 0.0,
                    kd: 0.05,
                },
                combination_model: CombinationModelConfig::SOSFS,
                turbulence_model: TurbulenceModelConfig::None,
            },
            enable_secondary_steering: true,
            enable_yaw_added_recovery: true,
            enable_transverse_velocities: true,
            enable_active_wake_mixing: false,
            enable_wake_mixing: false,
            use_parallel_calc: false,
        },
        turbine_library: "turbine_library".to_string(),
    };

    let mut model = FlorisModel::from_config(config).expect("Failed to create FlorisModel");

    // Use the manual approach
    model.initialize_grid().expect("init grid failed");
    model
        .initialize_flow_field()
        .expect("init flow field failed");

    let grid = model.core().grid.as_ref().unwrap();
    let turbine_grid = grid.as_any().downcast_ref::<TurbineGrid>().unwrap();

    let farm = model.core().farm.clone();
    let mut flow_field = model.core().flow_field.clone();
    let wake = model.core().wake.clone();

    florus::core::solver::sequential_solver(&farm, &mut flow_field, turbine_grid, &wake)
        .expect("solver failed");

    println!("=== DIRECT SOLVER RESULTS ===");
    println!(
        "u_sorted[fi=0, ti=0, 0, 0]: {:.6}",
        flow_field.u_sorted[[0, 0, 0, 0]]
    );
    println!(
        "u_sorted[fi=0, ti=1, 0, 0]: {:.6}",
        flow_field.u_sorted[[0, 1, 0, 0]]
    );

    // Compute power from u_sorted manually (like regression tests)
    use florus::core::turbines::cp_ct_table::TableConditions;
    let turbines = &farm.turbines;
    for ti in 0..3 {
        let mut v_cubed_sum = 0.0;
        let grid_points = flow_field.u_sorted.shape()[2] * flow_field.u_sorted.shape()[3];
        for iy in 0..flow_field.u_sorted.shape()[2] {
            for iz in 0..flow_field.u_sorted.shape()[3] {
                v_cubed_sum += flow_field.u_sorted[[0, ti, iy, iz]].powi(3);
            }
        }
        let rotor_avg_velocity = (v_cubed_sum / grid_points as f64).powf(1.0 / 3.0);
        println!(
            "turbine {}: rotor_avg_velocity={:.6}",
            ti, rotor_avg_velocity
        );

        let mut conditions = TableConditions::default();
        conditions.wind_speed = rotor_avg_velocity;
        let cp = turbines[ti]
            .turbine_type
            .power_thrust_table
            .cp_ct_table
            .get_cp(&conditions)
            .unwrap_or(0.0);
        println!("  cp = {:.6}", cp);
        println!("  cp * 1000 = {:.6}", cp * 1000.0);
    }
}
