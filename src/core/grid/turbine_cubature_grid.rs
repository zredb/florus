/// Turbine grid with cubature integration points
use super::*;
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::{reverse_rotate_coordinates_rel_west, rotate_coordinates_rel_west};
use ndarray::Array;
use std::f64::consts::FRAC_1_SQRT_2;
///与 TurbineGrid 并列的一种转子采样网格，其核心区别在于：它采用数值积分（cubature） 的方式，在转子圆面上布置非均匀的点并赋予不同权重，旨在更精确地估算转子上的平均入流风速，而非简单的矩形网格采样。
/// 它根据 cubature 积分规则，在圆面上布置特定位置的点，每个点有对应的积分权重，从而更好地逼近真实的转子平均风速。
#[derive(Debug, Clone)]
pub struct TurbineCubatureGrid {
    pub turbine_coordinates: Array2,
    pub turbine_diameters: Array1,
    pub wind_directions: Array1,
    pub grid_resolution: usize,
    pub n_turbines: usize,
    pub n_findex: usize,
    pub x_sorted: Array4,
    pub y_sorted: Array4,
    pub z_sorted: Array4,
    pub x_sorted_inertial_frame: Array4,
    pub y_sorted_inertial_frame: Array4,
    pub z_sorted_inertial_frame: Array4,
    pub cubature_weights: Array2,
    pub sorted_indices: Array2,
    pub sorted_coord_indices: Array2,
    pub unsorted_indices: Array2,
    pub x_center_of_rotation: Array2,
    pub y_center_of_rotation: Array2,
}

impl TurbineCubatureGrid {
    pub fn new(
        turbine_coordinates: Array2,
        turbine_diameters: Array1,
        wind_directions: Array1,
        grid_resolution: usize,
    ) -> crate::Result<Self> {
        // `grid_resolution` is the cubature ORDER N (Python grid.py:265-272,
        // 373-376): the rotor disk is sampled at N^2 points and N must be an
        // integer in [1, 10].
        if !(1..=10).contains(&grid_resolution) {
            return Err(anyhow::anyhow!(
                "Order of cubature integration must be between '1' and '10', given {}.",
                grid_resolution
            ));
        }

        let n_turbines = turbine_coordinates.shape()[0];
        let n_findex = wind_directions.len();

        // Polar cubature coefficients (Python get_cubature_coefficients,
        // grid.py:358-435): N radii |r_i| crossed with N unit direction
        // vectors (q_j, t_j) place N^2 points on the unit disk.
        let (r, t, q, a) = Self::get_cubature_coefficients(grid_resolution);
        let n = grid_resolution;
        let n_points = n * n;

        // Rotate coordinates
        let (
            rotated_coords_x,
            rotated_coords_y,
            _rotated_coords_z,
            x_center_of_rotation,
            y_center_of_rotation,
        ) = rotate_coordinates_rel_west(&wind_directions, &turbine_coordinates)?;

        // Get hub heights from turbine coordinates (z-coordinate)
        let hub_heights = turbine_coordinates.column(2).to_owned();

        // Grid points (Python grid.py:296-329):
        //   yv = np.kron(r, q); zv = np.kron(r, t)
        //   _y[:, ti] += yv * D[ti] / 2; _z[:, ti] += zv * D[ti] / 2
        // i.e. point k = (i * N + j) sits at (r_i * q_j, r_i * t_j) scaled by
        // the rotor radius, offset from the rotated turbine center.
        // Shape: (n_findex, n_turbines, n_points, 1)
        let mut x_sorted = Array::zeros((n_findex, n_turbines, n_points, 1));
        let mut y_sorted = Array::zeros((n_findex, n_turbines, n_points, 1));
        let mut z_sorted = Array::zeros((n_findex, n_turbines, n_points, 1));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let rotor_radius = turbine_diameters[ti] / 2.0;
                for pi in 0..n_points {
                    let (i, j) = (pi / n, pi % n);
                    x_sorted[[fi, ti, pi, 0]] = rotated_coords_x[[fi, ti]];
                    y_sorted[[fi, ti, pi, 0]] =
                        rotated_coords_y[[fi, ti]] + r[i] * q[j] * rotor_radius;
                    z_sorted[[fi, ti, pi, 0]] = hub_heights[ti] + r[i] * t[j] * rotor_radius;
                }
            }
        }

        // Integration weights (Python grid.py:301-304):
        //   np.kron(A, np.ones((1, N))) * B / np.pi,  B = np.pi / N
        // Row i holds N copies of A_i / N, so the (N, N) matrix sums to 1
        // (the A_i are normalized ring-area fractions). The weights are
        // turbine-independent; `simple_cubature` re-normalizes on use
        // (rotor_velocity.py:60-64), which is scale-invariant.
        let mut cubature_weights = Array::zeros((n, n));
        for i in 0..n {
            for j in 0..n {
                cubature_weights[[i, j]] = a[i] / n as Float;
            }
        }

        // Sort turbines by x coordinate (upstream to downstream)
        let mut sorted_indices = Array::zeros((n_findex, n_turbines));
        let mut sorted_coord_indices = Array::zeros((n_findex, n_turbines));
        let mut unsorted_indices = Array::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            let mut indices: Vec<(usize, Float)> = (0..n_turbines)
                .map(|i| (i, rotated_coords_x[[fi, i]]))
                .collect();
            indices.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            for (new_i, (old_i, _)) in indices.iter().enumerate() {
                sorted_indices[[fi, new_i]] = *old_i as Float;
                sorted_coord_indices[[fi, new_i]] = *old_i as Float;
                unsorted_indices[[fi, *old_i]] = new_i as Float;
            }
        }

        // Reverse rotate to get inertial frame coordinates
        let (x_sorted_inertial_frame, y_sorted_inertial_frame, z_sorted_inertial_frame) =
            reverse_rotate_coordinates_rel_west(
                &wind_directions,
                &x_sorted,
                &y_sorted,
                &z_sorted,
                &x_center_of_rotation,
                &y_center_of_rotation,
            )?;

        Ok(Self {
            turbine_coordinates,
            turbine_diameters,
            wind_directions,
            grid_resolution,
            n_turbines,
            n_findex,
            x_sorted,
            y_sorted,
            z_sorted,
            x_sorted_inertial_frame,
            y_sorted_inertial_frame,
            z_sorted_inertial_frame,
            cubature_weights,
            sorted_indices,
            sorted_coord_indices,
            unsorted_indices,
            x_center_of_rotation,
            y_center_of_rotation,
        })
    }

    /// Python `TurbineCubatureGrid.get_cubature_coefficients` (grid.py:358-435).
    ///
    /// Returns `(r, t, q, A)`: the radii, the z-direction cosines, the
    /// y-direction cosines, and the normalized ring-area fractions for the
    /// N-th order unit-disk cubature. `B = pi / N` is folded into the weight
    /// construction at the call site.
    ///
    /// The 26-digit literals are transcribed verbatim from the Python source;
    /// clippy::excessive_precision is allowed here so the tables stay diffable
    /// against the reference implementation (f64 keeps only 17 significant
    /// digits — the extra digits are inert). `FRAC_1_SQRT_2` is used where the
    /// literal is bit-identical to the std constant.
    #[allow(clippy::excessive_precision)]
    fn get_cubature_coefficients(n: usize) -> (Vec<Float>, Vec<Float>, Vec<Float>, Vec<Float>) {
        match n {
            1 => (vec![0.0], vec![0.0], vec![1.0], vec![1.0]),
            2 => (
                vec![-FRAC_1_SQRT_2, FRAC_1_SQRT_2],
                vec![-FRAC_1_SQRT_2, FRAC_1_SQRT_2],
                vec![FRAC_1_SQRT_2, FRAC_1_SQRT_2],
                vec![0.5, 0.5],
            ),
            3 => (
                vec![
                    -0.8164965809277260327324280,
                    0.0,
                    0.8164965809277260327324280,
                ],
                vec![
                    -0.8660254037844386467637232,
                    0.0,
                    0.8660254037844386467637232,
                ],
                vec![0.5, 1.0, 0.5],
                vec![0.375, 0.25, 0.375],
            ),
            4 => (
                vec![
                    -0.8880738339771152621607646,
                    -0.4597008433809830609776340,
                    0.4597008433809830609776340,
                    0.8880738339771152621607646,
                ],
                vec![
                    -0.9238795325112867561281832,
                    -0.3826834323650897717284600,
                    0.3826834323650897717284600,
                    0.9238795325112867561281832,
                ],
                vec![
                    0.3826834323650897717284600,
                    0.9238795325112867561281832,
                    0.9238795325112867561281832,
                    0.3826834323650897717284600,
                ],
                vec![0.25, 0.25, 0.25, 0.25],
            ),
            5 => (
                vec![
                    -0.9192110607898045793726291,
                    -0.5958615826865180525340234,
                    0.0,
                    0.5958615826865180525340234,
                    0.9192110607898045793726291,
                ],
                vec![
                    -0.9510565162951535721164393,
                    -0.5877852522924731291687060,
                    0.0,
                    0.5877852522924731291687060,
                    0.9510565162951535721164393,
                ],
                vec![
                    0.3090169943749474241022934,
                    0.8090169943749474241022934,
                    1.0,
                    0.8090169943749474241022934,
                    0.3090169943749474241022934,
                ],
                vec![
                    0.1882015313502336375250377,
                    0.2562429130942108069194067,
                    0.1111111111111111111111111,
                    0.2562429130942108069194067,
                    0.1882015313502336375250377,
                ],
            ),
            6 => (
                vec![
                    -0.9419651451198933233901941,
                    -FRAC_1_SQRT_2,
                    -0.3357106870197288066698994,
                    0.3357106870197288066698994,
                    FRAC_1_SQRT_2,
                    0.9419651451198933233901941,
                ],
                vec![
                    -0.9659258262890682867497432,
                    -FRAC_1_SQRT_2,
                    -0.2588190451025207623488988,
                    0.2588190451025207623488988,
                    FRAC_1_SQRT_2,
                    0.9659258262890682867497432,
                ],
                vec![
                    0.2588190451025207623488988,
                    FRAC_1_SQRT_2,
                    0.9659258262890682867497432,
                    0.9659258262890682867497432,
                    FRAC_1_SQRT_2,
                    0.2588190451025207623488988,
                ],
                vec![
                    0.1388888888888888888888889,
                    0.2222222222222222222222222,
                    0.1388888888888888888888889,
                    0.1388888888888888888888889,
                    0.2222222222222222222222222,
                    0.1388888888888888888888889,
                ],
            ),
            7 => (
                vec![
                    -0.9546790248493448767148503,
                    -0.7684615381131740734708478,
                    -0.4608042298407784190147371,
                    0.0,
                    0.4608042298407784190147371,
                    0.7684615381131740734708478,
                    0.9546790248493448767148503,
                ],
                vec![
                    -0.9749279121818236070181317,
                    -0.7818314824680298087084445,
                    -0.4338837391175581204757683,
                    0.0,
                    0.4338837391175581204757683,
                    0.7818314824680298087084445,
                    0.9749279121818236070181317,
                ],
                vec![
                    0.2225209339563144042889026,
                    0.6234898018587335305250049,
                    0.9009688679024191262361023,
                    1.0,
                    0.9009688679024191262361023,
                    0.6234898018587335305250049,
                    0.2225209339563144042889026,
                ],
                vec![
                    0.1102311055883841876377392,
                    0.1940967344215859403901162,
                    0.1644221599900298719721446,
                    0.0625,
                    0.1644221599900298719721446,
                    0.1940967344215859403901162,
                    0.1102311055883841876377392,
                ],
            ),
            8 => (
                vec![
                    -0.9646596061808674528345806,
                    -0.8185294874300058668603761,
                    -0.5744645143153507855310459,
                    -0.2634992299855422962484895,
                    0.2634992299855422962484895,
                    0.5744645143153507855310459,
                    0.8185294874300058668603761,
                    0.9646596061808674528345806,
                ],
                vec![
                    -0.9807852804032304491261822,
                    -0.8314696123025452370787884,
                    -0.5555702330196022247428308,
                    -0.1950903220161282678482849,
                    0.1950903220161282678482849,
                    0.5555702330196022247428308,
                    0.8314696123025452370787884,
                    0.9807852804032304491261822,
                ],
                vec![
                    0.1950903220161282678482849,
                    0.5555702330196022247428308,
                    0.8314696123025452370787884,
                    0.9807852804032304491261822,
                    0.9807852804032304491261822,
                    0.8314696123025452370787884,
                    0.5555702330196022247428308,
                    0.1950903220161282678482849,
                ],
                vec![
                    0.0869637112843634643432660,
                    0.1630362887156365356567340,
                    0.1630362887156365356567340,
                    0.0869637112843634643432660,
                    0.0869637112843634643432660,
                    0.1630362887156365356567340,
                    0.1630362887156365356567340,
                    0.0869637112843634643432660,
                ],
            ),
            9 => (
                vec![
                    -0.9710282199223060261836893,
                    -0.8503863747508400503582112,
                    -0.6452980455813291706201889,
                    -0.3738447061866471744516959,
                    0.0,
                    0.3738447061866471744516959,
                    0.6452980455813291706201889,
                    0.8503863747508400503582112,
                    0.9710282199223060261836893,
                ],
                vec![
                    -0.9848077530122080593667430,
                    -0.8660254037844386467637232,
                    -0.6427876096865393263226434,
                    -0.3420201433256687330440996,
                    0.0,
                    0.3420201433256687330440996,
                    0.6427876096865393263226434,
                    0.8660254037844386467637232,
                    0.9848077530122080593667430,
                ],
                vec![
                    0.1736481776669303488517166,
                    0.5,
                    0.7660444431189780352023927,
                    0.9396926207859083840541093,
                    1.0,
                    0.9396926207859083840541093,
                    0.7660444431189780352023927,
                    0.5,
                    0.1736481776669303488517166,
                ],
                vec![
                    0.0718567803956129706617061,
                    0.1406780075747310300960863,
                    0.1559132614878706270409275,
                    0.1115519505417853722012801,
                    0.04,
                    0.1115519505417853722012801,
                    0.1559132614878706270409275,
                    0.1406780075747310300960863,
                    0.0718567803956129706617061,
                ],
            ),
            10 => (
                vec![
                    -0.9762632447087885713212574,
                    -0.8770602345636481685478274,
                    -FRAC_1_SQRT_2,
                    -0.4803804169063914437972190,
                    -0.2165873427295972057980989,
                    0.2165873427295972057980989,
                    0.4803804169063914437972190,
                    FRAC_1_SQRT_2,
                    0.8770602345636481685478274,
                    0.9762632447087885713212574,
                ],
                vec![
                    -0.9876883405951377261900402,
                    -0.8910065241883678623597096,
                    -FRAC_1_SQRT_2,
                    -0.4539904997395467915604084,
                    -0.1564344650402308690101053,
                    0.1564344650402308690101053,
                    0.4539904997395467915604084,
                    FRAC_1_SQRT_2,
                    0.8910065241883678623597096,
                    0.9876883405951377261900402,
                ],
                vec![
                    0.1564344650402308690101053,
                    0.4539904997395467915604084,
                    FRAC_1_SQRT_2,
                    0.8910065241883678623597096,
                    0.9876883405951377261900402,
                    0.9876883405951377261900402,
                    0.8910065241883678623597096,
                    FRAC_1_SQRT_2,
                    0.4539904997395467915604084,
                    0.1564344650402308690101053,
                ],
                vec![
                    0.0592317212640472718785660,
                    0.1196571676248416170103229,
                    0.1422222222222222222222222,
                    0.1196571676248416170103229,
                    0.0592317212640472718785660,
                    0.0592317212640472718785660,
                    0.1196571676248416170103229,
                    0.1422222222222222222222222,
                    0.1196571676248416170103229,
                    0.0592317212640472718785660,
                ],
            ),
            _ => unreachable!("validated to [1, 10] at the call site"),
        }
    }
}

impl Grid for TurbineCubatureGrid {
    fn n_turbines(&self) -> usize {
        self.n_turbines
    }
    fn n_findex(&self) -> usize {
        self.n_findex
    }
    fn x_sorted(&self) -> &Array4 {
        &self.x_sorted
    }
    fn y_sorted(&self) -> &Array4 {
        &self.y_sorted
    }
    fn z_sorted(&self) -> &Array4 {
        &self.z_sorted
    }
    fn x_sorted_inertial_frame(&self) -> &Array4 {
        &self.x_sorted_inertial_frame
    }
    fn y_sorted_inertial_frame(&self) -> &Array4 {
        &self.y_sorted_inertial_frame
    }
    fn z_sorted_inertial_frame(&self) -> &Array4 {
        &self.z_sorted_inertial_frame
    }
    fn cubature_weights(&self) -> Option<&Array2> {
        Some(&self.cubature_weights)
    }
    fn average_method(&self) -> AveragingMethod {
        AveragingMethod::SimpleCubature
    }
    // The default impl indexes z_sorted[[0, ti, 0, 1]], which assumes the
    // TurbineGrid's (n_grid, n_grid) rotor layout; the cubature grid's z
    // axis has length 1 and its points sit off-hub. The true hub heights
    // are the turbine coordinate z values.
    fn hub_heights(&self) -> Array1 {
        self.turbine_coordinates.column(2).to_owned()
    }
    fn sorted_indices(&self) -> &Array2 {
        &self.sorted_indices
    }
    fn sorted_coord_indices(&self) -> &Array2 {
        &self.sorted_coord_indices
    }
    fn resolution(&self) -> usize {
        self.grid_resolution
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl Default for TurbineCubatureGrid {
    fn default() -> Self {
        let empty_1d = Array::zeros((0,));
        let empty_2d = Array::zeros((0, 0));
        let empty_4d = Array::zeros((0, 0, 0, 0));
        Self {
            turbine_coordinates: empty_2d.clone(),
            turbine_diameters: empty_1d.clone(),
            wind_directions: empty_1d.clone(),
            grid_resolution: 0,
            n_turbines: 0,
            n_findex: 0,
            x_sorted: empty_4d.clone(),
            y_sorted: empty_4d.clone(),
            z_sorted: empty_4d.clone(),
            x_sorted_inertial_frame: empty_4d.clone(),
            y_sorted_inertial_frame: empty_4d.clone(),
            z_sorted_inertial_frame: empty_4d.clone(),
            cubature_weights: empty_2d.clone(),
            sorted_indices: empty_2d.clone(),
            sorted_coord_indices: empty_2d.clone(),
            unsorted_indices: empty_2d.clone(),
            x_center_of_rotation: empty_2d.clone(),
            y_center_of_rotation: empty_2d.clone(),
        }
    }
}
