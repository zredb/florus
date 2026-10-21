//! Example: Approximate Model Parameters
//!
//! This example demonstrates how to use the UncertainFlorisModel class to
//! analyze the impact of uncertain wind direction on power results.
//!
//! The example generates a random walk on the wind speed, runs it with the
//! full `FlorisModel` and with the approximate model (the wind speeds are
//! rounded to 0.5 m/s, so that the conditions can be reused) and compares
//! the expected power of the two models.
//!
//! Scenario aligned 1:1 with examples_uncertain/002_approx_floris_model.py.
//!
//! Note: the Python example draws the random walk with
//! `numpy.random.default_rng(0)`, whose stream cannot be reproduced without
//! numpy's PCG64 implementation and its ziggurat normal distribution. The
//! Rust example uses its own (documented) random number generator, so the
//! wind speed series and therefore the expected power differ slightly from
//! Python's (the two agree to within a few tenths of a percent).

use florus::types::Array1;
use florus::uncertain::ApproxFlorisModel;
use florus::{FlorisModel, Result, TimeSeries};

/// A small deterministic random number generator (xorshift64*).
///
/// It is used to draw the random walk on the wind speed; Python uses
/// `numpy.random.default_rng(0)` instead.
struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next uniformly distributed number in [0, 1)
    fn next_f64(&mut self) -> f64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        let value = x.wrapping_mul(0x2545F4914F6CDD1D);
        ((value >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// The next normally distributed number (Box-Muller)
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

fn main() -> Result<()> {
    // Generate time series data using a random walk on wind speeds with
    // constant wind direction
    let n = 5000;
    let n_turbines = 25;

    // Random walk on wind speed with values between 5 and 20 m/s
    let ws_ll = 5.0;
    let ws_ul = 20.0;
    let mut rng = XorShift64::new(0x9E37_79B9_7F4A_7C15);
    let mut wind_speeds = vec![10.0; n];
    for i in 1..n {
        wind_speeds[i] = wind_speeds[i - 1] + 0.25 * rng.next_normal();
        if wind_speeds[i] < ws_ll {
            wind_speeds[i] = ws_ll;
        }
        if wind_speeds[i] > ws_ul {
            wind_speeds[i] = ws_ul;
        }
    }

    let time_series = TimeSeries::new(
        Array1::from_elem(n, 270.0),
        Array1::from_vec(wind_speeds),
        Array1::from_elem(n, 0.06),
    )?;

    // Instantiate a FlorisModel and an ApproxFlorisModel
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut afmodel =
        ApproxFlorisModel::new(FlorisModel::from_file("examples/inputs/gch.yaml")?, 0.5)?;

    // Set both models to an n_turbine layout and use the above time series
    let layout_x = Array1::from_iter((0..n_turbines).map(|i| i as f64 * 500.0));
    let layout_y = Array1::from_elem(n_turbines, 0.0);
    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.set_wind_data(&time_series)?;
    afmodel.set_layout(&layout_x, &layout_y)?;
    afmodel.set_wind_data(&time_series)?;

    // Now time both runs to show the speedup from approximating the wind speed
    let start = std::time::Instant::now();
    fmodel.run()?;
    println!(
        "FlorisModel run time: {:.3} s",
        start.elapsed().as_secs_f64()
    );

    let start = std::time::Instant::now();
    afmodel.run()?;
    println!(
        "ApproxFlorisModel run time: {:.3} s",
        start.elapsed().as_secs_f64()
    );

    // The number of conditions that have to be solved
    println!(
        "FlorisModel conditions: {}, ApproxFlorisModel conditions: {}",
        fmodel.flow_field().n_findex,
        afmodel.model.n_unique()
    );

    // Compare the expected power results
    println!(
        "Expected power from FlorisModel: {:.2} MW",
        fmodel.get_expected_farm_power(None, None)? / 1e6
    );
    println!(
        "Expected power from ApproxFlorisModel: {:.2} MW",
        afmodel.get_expected_farm_power(None)? / 1e6
    );

    Ok(())
}
