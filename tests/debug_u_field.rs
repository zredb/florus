use florus::floris_model::FlorisModel;
use florus::types::{Array1, Float};
use ndarray::Array;

#[test]
fn debug_u_field_after_run() {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();

    fmodel
        .set(
            Some(Array1::from_vec(vec![8.0, 8.0, 8.0])),
            Some(Array1::from_vec(vec![270.0, 270.0, 270.0])),
            None, // wind_shear
            None, // wind_veer
            None, // reference_wind_height
            Some(Array1::from_vec(vec![0.1, 0.1, 0.1])),
            None, // air_density
            None, // layout_x
            None, // layout_y
            None, // yaw_angles
            None, // power_setpoints
            None, // awc_modes
            None, // awc_amplitudes
            None, // awc_frequencies
            None, // disable_turbines
        )
        .unwrap();

    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers();
    println!("Powers: {:?}", powers);

    let flow_field = fmodel.flow_field();
    println!("u shape: {:?}", flow_field.u.shape());
    println!("u_sorted shape: {:?}", flow_field.u_sorted.shape());

    // Check if u has meaningful values
    for ti in 0..3 {
        let u_vals: Vec<Float> = (0..flow_field.u.shape()[2])
            .flat_map(|iy| {
                (0..flow_field.u.shape()[3]).map(move |iz| flow_field.u[[0, ti, iy, iz]])
            })
            .collect();
        println!(
            "Turbine {} u values (first 5): {:?}",
            ti,
            &u_vals[..5.min(u_vals.len())]
        );
    }
}
