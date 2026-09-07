//! Temporary comparison dump: Rust horizontal plane vs Python reference.
use florus::{FlorisModel, Result};
use std::fs::File;
use std::io::Write;

fn main() -> Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    let layout_x = vec![0.0, 500.0, 0.0, 128.0, 1000.0, 900.0, 1500.0, 1250.0];
    let layout_y = vec![0.0, 300.0, 750.0, 1400.0, 0.0, 567.0, 888.0, 1450.0];
    fmodel.set_layout(
        &ndarray::arr1(&layout_x),
        &ndarray::arr1(&layout_y),
    )?;

    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[290.0]),
        ndarray::arr1(&[0.06]),
    )?;

    let plane = fmodel.calculate_horizontal_plane(90.0, 200, 100)?;

    let x1 = &plane.data.x1;
    let x2 = &plane.data.x2;
    let u = &plane.data.u;
    let v = &plane.data.v;
    let w = &plane.data.w;

    println!("n_points = {}", u.len());
    let umin = u.iter().cloned().fold(f64::INFINITY, f64::min);
    let umax = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let umean = u.iter().sum::<f64>() / u.len() as f64;
    println!("u: min={:.8} max={:.8} mean={:.8}", umin, umax, umean);

    let mut f = File::create("output/rs_plane.csv")?;
    writeln!(f, "x1,x2,u,v,w")?;
    for i in 0..u.len() {
        writeln!(f, "{},{},{},{},{}", x1[i], x2[i], u[i], v[i], w[i])?;
    }
    println!("WROTE output/rs_plane.csv");

    Ok(())
}
