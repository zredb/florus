//! Controller-dependent (TUM) turbine operation model
//!
//! Port of `floris.core.turbine.controller_dependent_operation_model.py`
//! (Tamaro, Campagnolo, Bottasso — TU München, wes-2023-133).
//!
//! The model predicts power/thrust from blade-element-momentum closed-form
//! expressions (Eq. 25 of Tamaro et al.) evaluated at the tip-speed-ratio and
//! pitch angle that the turbine's standard controller would command
//! (`control_trajectory`, region 2 / 2.5 / 3 logic). The yawed/unyawed power
//! or thrust ratio from the BEM expressions is applied to the Cp/Ct surface
//! from the turbine's `.npz` lookup tables, interpolated with a
//! not-a-knot tensor-product cubic spline (scipy
//! `RegularGridInterpolator(..., method="cubic")` = `make_ndbspl`, which
//! documents not-a-knot boundary conditions in all dimensions).
//!
//! Numerical root finding uses a damped-Newton stand-in for
//! `scipy.optimize.fsolve` (MINPACK hybrd, n=1); all three solve levels find
//! the same smooth, monotone roots.

use crate::core::turbines::operation_models::{
    ControllerDependentParams, OperationModel, TurbineContext, TurbineParameters,
};
use crate::types::{Array2, Array4, Float};
use ndarray::s;
///由慕尼黑工业大学（TUM）开发的高级操作模型。它不再依赖简单的余弦损失公式，而是试图模拟风机控制器在偏航、倾斜及降额工况下的实际调节行为，从而更物理地计算功率与推力.
/// 核心逻辑：模拟控制轨迹
///该模型的关键在于 control_trajectory 方法。它会根据来流风速和操作指令，预测风机实际会运行在哪个叶尖速比（TSR）和桨距角上。其控制逻辑假设了标准的 Region 2（最优 Cp 跟踪）和 Region 3（额定功率保持），并考虑了 Region 2.5（过渡区）的存在。得到 TSR 和桨距角后，模型再查询 Cp / Ct 查找表来计算最终性能。

///区别于余弦损失的物理考量
///与 cosine-loss 的固定指数衰减不同，该模型引入了两个重要的物理效应：

///垂直风切变的局部计算：模型会计算每个转子位置处的局部垂直风切变。这意味着，当风机朝不同方向偏航时，由于切变的存在，其性能表现会不对称（例如，向左侧偏航和向右侧偏航的功率损失可能不同）。风切变会受到上游尾流的影响而改变。

///同时处理偏航与倾斜：模型明确包含偏航角（yaw）和倾斜角（tilt）对性能的影响，而不仅仅处理偏航。

///所需的关键参数
///由于需要模拟控制器和气动细节，该模型要求风机定义中包含比简单模型多得多的参数：

///控制器参数：额定转速（rated_rpm）、发电机效率（generator_efficiency）、额定功率（rated_power）。

///叶片气动参数：叶片扭角（beta）、阻力系数（cd）、升力线斜率（cl_alfa）、转子实度（rotor_solidity）。

///查找表：Cp 和 Ct 作为叶尖速比和桨距角函数的二维查找表（cp_ct_data）。

///网格要求：为了捕捉局部风切变效应，必须确保 turbine_grid_points > 1。

///已知局限与注意事项
///计算速度较慢：因为需要求解平衡方程和进行插值，它比 cosine-loss 等简单模型慢。

///演示性数据：FLORIS 自带的参考风机（如 IEA 15MW）的 Cp/Ct 表是演示性质的，由这些表反推得到的功率/推力曲线与官方参考曲///线并非完全一致。因此，使用该模型得到的结果应被视为“演示性”的，而非绝对精确的基准。

///参数获取：对于非内置风机，需要自行估算或从文献中获取上述叶片与控制参数，门槛较高。
#[derive(Debug, Clone, Default)]
pub struct ControllerDependentTurbine;

/// Which Python entry point's conventions `solve_pair` applies to the BEM
/// solve results (see `ControllerDependentTurbine::solve_pair` docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SolveQuantity {
    /// `power()`: strict `ier == 1` gate, Cp on success, `-1e3` otherwise.
    Power,
    /// `thrust_coefficient()`: convergence ignored; yawed clipped, unyawed raw.
    Thrust,
}

// ===================== numerics: 1D not-a-knot cubic spline =====================

/// 1D cubic spline with not-a-knot end conditions, matching scipy's
/// `make_interp_spline` default bc (the unique not-a-knot spline).
struct CubicSpline {
    x: Vec<Float>,
    y: Vec<Float>,
    /// Second derivatives at the knots (moment form).
    m: Vec<Float>,
}

impl CubicSpline {
    fn fit(x: &[Float], y: &[Float]) -> CubicSpline {
        let n = x.len();
        debug_assert!(n >= 2);
        if n < 4 {
            // Degenerate for not-a-knot: fall back to the natural spline.
            return Self::fit_natural(x, y);
        }

        let h: Vec<Float> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();

        // Moment system M[0..n-1] (second derivatives), not-a-knot ends:
        //   row 0 (left):   M0 - (1 + h0/h1) M1 + (h0/h1) M2 = 0
        //   row i (1..=n-2): h[i-1] M[i-1] + 2(h[i-1]+h[i]) M[i] + h[i] M[i+1]
        //                    = 6 * (slope[i] - slope[i-1])
        //   row n-1 (right): (h[n-2]/h[n-3]) M[n-3] - (1+h[n-2]/h[n-3]) M[n-2] + M[n-1] = 0
        // Rows 0 and n-1 reach two entries to the right/left of the diagonal,
        // so solve the (n x n) system densely (n <= ~40 for the surfaces).
        let mut rows: Vec<(usize, usize, Float)> = Vec::with_capacity(3 * n);
        let mut rhs = vec![0.0; n];
        // Row 0
        rows.push((0, 0, 1.0));
        rows.push((0, 1, -(1.0 + h[0] / h[1])));
        rows.push((0, 2, h[0] / h[1]));
        // Interior rows
        for i in 1..=n - 2 {
            let d = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
            rows.push((i, i - 1, h[i - 1]));
            rows.push((i, i, 2.0 * (h[i - 1] + h[i])));
            rows.push((i, i + 1, h[i]));
            rhs[i] = d;
        }
        // Row n-1
        let (hm2, hm3) = (h[n - 2], h[n - 3]);
        rows.push((n - 1, n - 3, hm2 / hm3));
        rows.push((n - 1, n - 2, -(1.0 + hm2 / hm3)));
        rows.push((n - 1, n - 1, 1.0));

        // Gaussian elimination on the pentadiagonal-in-two-rows system.
        // Rows 0 and n-1 break the tridiagonal pattern; solve with a small
        // dense LU (n <= 40 for the turbine surfaces).
        let mut mat = vec![vec![0.0; n]; n];
        for (i, j, v) in &rows {
            mat[*i][*j] = *v;
        }
        let m = dense_solve(&mut mat, &mut rhs);

        CubicSpline {
            x: x.to_vec(),
            y: y.to_vec(),
            m,
        }
    }

    fn fit_natural(x: &[Float], y: &[Float]) -> CubicSpline {
        let n = x.len();
        if n == 2 {
            // straight line (zero second derivatives)
            return CubicSpline {
                x: x.to_vec(),
                y: y.to_vec(),
                m: vec![0.0; 2],
            };
        }
        let h: Vec<Float> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
        let mut mat = vec![vec![0.0; n]; n];
        let mut rhs = vec![0.0; n];
        mat[0][0] = 1.0;
        mat[n - 1][n - 1] = 1.0;
        for i in 1..n - 1 {
            let d = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
            mat[i][i - 1] = h[i - 1];
            mat[i][i] = 2.0 * (h[i - 1] + h[i]);
            mat[i][i + 1] = h[i];
            rhs[i] = d;
        }
        let m = dense_solve(&mut mat, &mut rhs);
        CubicSpline {
            x: x.to_vec(),
            y: y.to_vec(),
            m,
        }
    }

    /// Evaluate at `t`; outside the knot range the boundary interval's cubic
    /// extrapolates (scipy RGI with fill_value=None extrapolates too).
    fn eval(&self, t: Float) -> Float {
        let n = self.x.len();
        if n < 2 {
            return self.y.first().copied().unwrap_or(0.0);
        }
        // interval selection (clamped; extrapolates via end intervals)
        let i = if t <= self.x[0] {
            0
        } else if t >= self.x[n - 1] {
            n - 2
        } else {
            // binary search: largest i with x[i] <= t
            let (mut lo, mut hi) = (0usize, n - 1);
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if self.x[mid] <= t {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        let (x0, x1) = (self.x[i], self.x[i + 1]);
        let h = x1 - x0;
        let (m0, m1) = (self.m[i], self.m[i + 1]);
        let (y0, y1) = (self.y[i], self.y[i + 1]);
        let a = x1 - t;
        let b = t - x0;
        m0 * a * a * a / (6.0 * h)
            + m1 * b * b * b / (6.0 * h)
            + (y0 / h - m0 * h / 6.0) * a
            + (y1 / h - m1 * h / 6.0) * b
    }
}

/// Small dense Gaussian elimination with partial pivoting (n <= ~40).
fn dense_solve(a: &mut [Vec<Float>], b: &mut [Float]) -> Vec<Float> {
    let n = b.len();
    for k in 0..n {
        // pivot
        let mut piv = k;
        for i in k + 1..n {
            if a[i][k].abs() > a[piv][k].abs() {
                piv = i;
            }
        }
        if a[piv][k].abs() < 1e-300 {
            continue;
        }
        if piv != k {
            a.swap(k, piv);
            b.swap(k, piv);
        }
        for i in k + 1..n {
            let f = a[i][k] / a[k][k];
            if f == 0.0 {
                continue;
            }
            for j in k..n {
                a[i][j] -= f * a[k][j];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s -= a[i][j] * x[j];
        }
        if a[i][i].abs() > 1e-300 {
            x[i] = s / a[i][i];
        }
    }
    x
}

// ===================== numerics: misc =====================

/// `np.interp`: 1D linear interpolation clamped to the endpoints.
fn np_interp(x: Float, xp: &[Float], fp: &[Float]) -> Float {
    let n = xp.len();
    debug_assert!(n == fp.len() && n >= 2);
    if x <= xp[0] {
        return fp[0];
    }
    if x >= xp[n - 1] {
        return fp[n - 1];
    }
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xp[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (x - xp[lo]) / (xp[hi] - xp[lo]);
    fp[lo] + t * (fp[hi] - fp[lo])
}

/// `scipy.optimize.fsolve` equivalent for smooth 1D residuals: damped Newton
/// with a central-difference Jacobian and step halving. Returns
/// `(x, |f(x)|, converged)`.
fn fsolve1d(
    f: &dyn Fn(Float) -> Float,
    x0: Float,
    xtol: Float,
    maxfev: usize,
) -> (Float, Float, bool) {
    let max_iter = (maxfev / 3).max(20);
    let mut x = x0;
    let mut fx = f(x);
    if fx.abs() <= 1e-13 {
        return (x, fx.abs(), true);
    }
    for _ in 0..max_iter {
        let h = 1e-7 * x.abs().max(1.0);
        let df = (f(x + h) - f(x - h)) / (2.0 * h);
        if !df.is_finite() || df.abs() < 1e-300 {
            return (x, fx.abs(), false);
        }
        let mut step = -fx / df;
        // hybrd-like trust clamp to avoid divergence
        let max_step = 1e3 * x.abs().max(1.0);
        if step.abs() > max_step {
            step = step.signum() * max_step;
        }
        let step0 = step;
        let mut advanced = false;
        for _ in 0..40 {
            let xn = x + step;
            let fnv = f(xn);
            if fnv.abs() < fx.abs() {
                x = xn;
                fx = fnv;
                advanced = true;
                break;
            }
            step *= 0.5;
        }
        if !advanced {
            return (x, fx.abs(), false);
        }
        if step0.abs() <= xtol * x.abs().max(1.0) || fx.abs() <= 1e-13 {
            return (x, fx.abs(), true);
        }
    }
    (x, fx.abs(), fx.abs() <= xtol * x.abs().max(1.0))
}

/// `np.polyfit(x, y, 1)[0]` — closed-form least-squares slope.
fn polyfit_slope(x: &[Float], y: &[Float]) -> Float {
    let n = x.len() as Float;
    let xm: Float = x.iter().sum::<Float>() / n;
    let ym: Float = y.iter().sum::<Float>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..x.len() {
        num += (x[i] - xm) * (y[i] - ym);
        den += (x[i] - xm) * (x[i] - xm);
    }
    if den.abs() < 1e-300 {
        return 0.0;
    }
    num / den
}

/// `np.linspace(start, stop, num)`.
fn linspace(start: Float, stop: Float, num: usize) -> Vec<Float> {
    if num == 1 {
        return vec![start];
    }
    (0..num)
        .map(|i| start + (stop - start) * i as Float / (num - 1) as Float)
        .collect()
}

// ===================== tensor-product cubic surface =====================

/// Cubic (not-a-knot) tensor-product interpolation over a regular
/// (tsr, pitch) grid — scipy `RegularGridInterpolator(method="cubic")`.
struct CubicRgi<'a> {
    tsr_lut: Vec<Float>,
    pitch_lut: Vec<Float>,
    /// Moments of the tsr-axis splines, one per pitch column.
    tsr_moments: Vec<Vec<Float>>,
    /// Column values (values[:, j]) reference, shape (n_tsr, n_pitch).
    values: &'a Array2,
}

impl<'a> CubicRgi<'a> {
    fn new(tsr_lut: &[Float], pitch_lut: &[Float], values: &'a Array2) -> CubicRgi<'a> {
        let n_pitch = pitch_lut.len();
        let mut tsr_moments = Vec::with_capacity(n_pitch);
        for j in 0..n_pitch {
            let col: Vec<Float> = values.slice(s![.., j]).to_vec();
            let sp = CubicSpline::fit(tsr_lut, &col);
            tsr_moments.push(sp.m);
        }
        CubicRgi {
            tsr_lut: tsr_lut.to_vec(),
            pitch_lut: pitch_lut.to_vec(),
            tsr_moments,
            values,
        }
    }

    /// Evaluate the tensor-product spline at (tsr, pitch). Separability:
    /// interpolate along tsr for each pitch column, then along pitch.
    fn eval(&self, tsr: Float, pitch: Float) -> Float {
        let n_pitch = self.pitch_lut.len();
        let n_tsr = self.tsr_lut.len();
        let mut v = vec![0.0; n_pitch];
        for j in 0..n_pitch {
            // spline over tsr with column j values/moments, evaluated at tsr
            v[j] = eval_with_moments(
                &self.tsr_lut,
                &self.values.slice(s![.., j]).to_vec(),
                &self.tsr_moments[j],
                tsr,
            );
            let _ = n_tsr;
        }
        // fit a spline over pitch with v, evaluate at pitch
        let sp = CubicSpline::fit(&self.pitch_lut, &v);
        sp.eval(pitch)
    }
}

/// Interval-clamped moment-form cubic evaluation (mirrors CubicSpline::eval
/// but with externally supplied knots/values/moments to avoid reallocating).
fn eval_with_moments(x: &[Float], y: &[Float], m: &[Float], t: Float) -> Float {
    let n = x.len();
    if n < 2 {
        return y.first().copied().unwrap_or(0.0);
    }
    let i = if t <= x[0] {
        0
    } else if t >= x[n - 1] {
        n - 2
    } else {
        let (mut lo, mut hi) = (0usize, n - 1);
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if x[mid] <= t {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    };
    let (x0, x1) = (x[i], x[i + 1]);
    let h = x1 - x0;
    let (m0, m1) = (m[i], m[i + 1]);
    let (y0, y1) = (y[i], y[i + 1]);
    let a = x1 - t;
    let b = t - x0;
    m0 * a * a * a / (6.0 * h)
        + m1 * b * b * b / (6.0 * h)
        + (y0 / h - m0 * h / 6.0) * a
        + (y1 / h - m1 * h / 6.0) * b
}

// ===================== BEM closed forms (Tamaro et al., Eq. 25) =====================

use crate::utilities::{cosd, sind};

/// Arguments shared by `find_cp` and `get_ct` (the Python `data` tuple).
#[derive(Debug, Clone, Copy)]
struct BemArgs {
    sigma: Float,
    cd: Float,
    cl_alfa: Float,
    /// yaw angle in degrees
    gamma: Float,
    /// tilt angle in degrees (already sign-flipped by the caller)
    delta: Float,
    /// local vertical shear exponent
    k: Float,
    cos_mu: Float,
    sin_mu: Float,
    tsr: Float,
    /// theta = rad(pitch + beta), radians
    theta: Float,
    /// overall misalignment in radians
    mu: Float,
}

/// Python ControllerDependentTurbine.find_cp (lines 868-934).
fn find_cp(a_args: &BemArgs, ct: Float) -> Float {
    let (mut mu, mut sin_mu, mut cos_mu) = (a_args.mu, a_args.sin_mu, a_args.cos_mu);
    if mu == 0.0 {
        mu = 1e-6;
        sin_mu = mu.sin();
        cos_mu = mu.cos();
    }
    // Axial induction (Eq. 25a)
    let a = 1.0
        - (1.0 + (1.0 - ct - (1.0 / 16.0) * sin_mu * sin_mu * ct * ct).sqrt())
            / (2.0 * (1.0 + (1.0 / 16.0) * ct * sin_mu * sin_mu));

    let sg = sind(a_args.gamma);
    let cg = cosd(a_args.gamma);
    let sd = sind(a_args.delta);
    let cd = cosd(a_args.delta);
    let k_1s =
        -1.0 * (15.0 * std::f64::consts::PI / 32.0) * ((mu + sin_mu * (ct / 2.0)) / 2.0).tan();
    let (tsr, theta, cd_v, cl_alfa, k) = (
        a_args.tsr,
        a_args.theta,
        a_args.cd,
        a_args.cl_alfa,
        a_args.k,
    );
    let pi = std::f64::consts::PI;

    let bracket = pi * cos_mu * cos_mu * tsr * cl_alfa * (a - 1.0) * (a - 1.0)
        - (tsr
            * cd_v
            * pi
            * (cd * cd * cg * cg * sd * sd * k * k + 3.0 * cd * cd * sg * sg * k * k
                - 8.0 * cd * tsr * sg * k
                + 8.0 * tsr * tsr))
            / 16.0
        - (pi * tsr * sin_mu * sin_mu * cd_v) / 2.0
        - (2.0 * pi * cos_mu * tsr * tsr * cl_alfa * theta) / 3.0
        + (pi * cos_mu * cos_mu * k_1s * k_1s * tsr * a * a * cl_alfa) / 4.0
        + (2.0 * pi * cos_mu * tsr * tsr * a * cl_alfa * theta) / 3.0
        + (2.0 * pi * cd * cos_mu * tsr * sg * cl_alfa * k * theta) / 3.0
        + (cd
            * cd
            * cos_mu
            * cos_mu
            * tsr
            * cl_alfa
            * k
            * k
            * pi
            * (a - 1.0)
            * (a - 1.0)
            * (cg * cg * sd * sd + sg * sg))
            / (4.0 * sin_mu * sin_mu)
        - (2.0 * pi * cd * cos_mu * tsr * sg * a * cl_alfa * k * theta) / 3.0
        + (cd
            * cd
            * cos_mu
            * cos_mu
            * k_1s
            * k_1s
            * tsr
            * a
            * a
            * cl_alfa
            * k
            * k
            * pi
            * (3.0 * cg * cg * sd * sd + sg * sg))
            / (24.0 * sin_mu * sin_mu)
        - (pi * cd * cg * cos_mu * cos_mu * k_1s * tsr * sd * a * cl_alfa * k) / sin_mu
        + (pi * cd * cg * cos_mu * cos_mu * k_1s * tsr * sd * a * a * cl_alfa * k) / sin_mu
        + (pi * cd * cg * cos_mu * k_1s * tsr * tsr * sd * a * cl_alfa * k * theta)
            / (5.0 * sin_mu)
        - (pi * cd * cd * cg * cos_mu * k_1s * tsr * sd * sg * a * cl_alfa * k * k * theta)
            / (10.0 * sin_mu);

    a_args.sigma * bracket / (2.0 * pi)
}

/// Python ControllerDependentTurbine.get_ct (lines 936-986) — residual in Ct.
fn get_ct_resid(a_args: &BemArgs, x: Float) -> Float {
    let (mut sin_mu, mut cos_mu) = (a_args.sin_mu, a_args.cos_mu);
    if a_args.mu == 0.0 {
        sin_mu = 1e-6_f64.sin();
        cos_mu = 1e-6_f64.cos();
    }
    let cd = cosd(a_args.delta);
    let cg = cosd(a_args.gamma);
    let sd = sind(a_args.delta);
    let sg = sind(a_args.gamma);

    // Axial induction with x as Ct
    let a = 1.0
        - (1.0 + (1.0 - x - (1.0 / 16.0) * x * x * sin_mu * sin_mu).sqrt())
            / (2.0 * (1.0 + (1.0 / 16.0) * x * sin_mu * sin_mu));

    let k_1s = -1.0
        * (15.0 * std::f64::consts::PI / 32.0)
        * ((a_args.mu + sin_mu * (x / 2.0)) / 2.0).tan();
    let pi = std::f64::consts::PI;

    let i1 = -(pi * cos_mu * (a_args.tsr - cd * sg * a_args.k) * (a - 1.0)
        + (cd
            * cg
            * cos_mu
            * k_1s
            * sd
            * a
            * a_args.k
            * pi
            * (2.0 * a_args.tsr - cd * sg * a_args.k))
            / (8.0 * sin_mu))
        / (2.0 * pi);

    let i2 = (pi * sin_mu * sin_mu
        + (pi
            * (cd * cd * cg * cg * sd * sd * a_args.k * a_args.k
                + 3.0 * cd * cd * sg * sg * a_args.k * a_args.k
                - 8.0 * cd * a_args.tsr * sg * a_args.k
                + 8.0 * a_args.tsr * a_args.tsr))
            / 12.0)
        / (2.0 * pi);

    (a_args.sigma * (a_args.cd + a_args.cl_alfa) * i1
        - a_args.sigma * a_args.cl_alfa * a_args.theta * i2)
        - x
}

/// Solve the BEM Ct root for the given misalignment configuration.
fn solve_ct(args: &BemArgs, x0: Float) -> (Float, bool) {
    let (ct, _resid, ok) = fsolve1d(&|x| get_ct_resid(args, x), x0, 1.49e-8, 400);
    (ct, ok)
}

// ===================== local vertical shear =====================

/// Python compute_local_vertical_shear (lines 418-451): per-turbine linear
/// vertical shear exponent from the rotor velocity profile. The even-length
/// normalization is reproduced verbatim from the Python source
/// (`mean / (m[mid] + m[mid-1]) / 2`).
fn compute_local_vertical_shear(velocities: &Array4) -> crate::Result<Array2> {
    let shape = velocities.shape();
    let (n_findex, n_turbines, n_y, n_z) = (shape[0], shape[1], shape[2], shape[3]);
    if n_z == 1 {
        anyhow::bail!(
            "The ControllerDependentTurbine computes a local shear based on inflow wind \
             speeds across the rotor. The provided velocities do not contain a vertical \
             profile. This can occur if n_grid is set to 1 in the FLORIS input yaml."
        );
    }
    let mut shear = Array2::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            // mean over the y axis -> profile over z
            let mut mean_speed = vec![0.0; n_z];
            for (iz, m) in mean_speed.iter_mut().enumerate() {
                let mut sum = 0.0;
                for iy in 0..n_y {
                    sum += velocities[[fi, ti, iy, iz]];
                }
                *m = sum / n_y as Float;
            }
            // normalized by the hub-height (middle) speed
            let u_u_hh: Vec<Float> = if mean_speed.len() % 2 != 0 {
                let mid = mean_speed.len() / 2;
                let hh = mean_speed[mid];
                mean_speed.iter().map(|v| v / hh).collect()
            } else {
                let mid = mean_speed.len() / 2;
                let denom = mean_speed[mid] + mean_speed[mid - 1];
                mean_speed.iter().map(|v| v / denom / 2.0).collect()
            };
            // zg_R = linspace(-1, 1, n+2), interior points
            let zg = linspace(-1.0, 1.0, mean_speed.len() + 2);
            let xs: Vec<Float> = zg[1..zg.len() - 1].to_vec();
            let ys: Vec<Float> = u_u_hh.iter().map(|v| 1.0 - v).collect();
            shear[[fi, ti]] = -polyfit_slope(&xs, &ys);
        }
    }
    Ok(shear)
}

// ===================== control trajectory =====================

struct TrajectoryInput<'a> {
    params: &'a ControllerDependentParams,
    rotor_effective_velocities: &'a Array2,
    yaw_angles: &'a Array2,
    tilt_angles: &'a Array2,
    shear: &'a Array2,
    power_demanded: &'a Array2,
    /// Cubic interpolator over the Cp surface (borrowed luts).
    cp_surface: &'a dyn Fn(Float, Float) -> Float,
}

/// Python control_trajectory (lines 453-866). Returns (pitch_out, tsr_out).
#[allow(clippy::too_many_arguments)]
fn control_trajectory(inp: &TrajectoryInput) -> crate::Result<(Array2, Array2)> {
    let p = inp.params;
    let (n_findex, n_turbines) = inp.rotor_effective_velocities.dim();

    // Optimum point of the Cp surface (first row-major occurrence of the max)
    let mut max_cp = Float::NEG_INFINITY;
    let (mut tsr_opt, mut pitch_opt) = (0.0, 0.0);
    for i in 0..p.cp_lut.shape()[0] {
        for j in 0..p.cp_lut.shape()[1] {
            let v = p.cp_lut[[i, j]];
            if v > max_cp {
                max_cp = v;
                tsr_opt = p.tsr_lut[i];
                pitch_opt = p.pitch_lut[j];
            }
        }
    }

    let air_density = p.ref_air_density; // Python uses ref density here (line 83)
    let omega_cut_in = 0.0; // RPM
    let omega_max = p.rated_rpm;
    let rated_power_aero_full = p.rated_power_kw * 1000.0 / p.generator_efficiency;

    // Torque-rpm curve and region 2.5 detection
    let omega_array = linspace(omega_cut_in, omega_max, 161);
    let omega_rad: Vec<Float> = omega_array
        .iter()
        .map(|w| w * std::f64::consts::PI / 30.0)
        .collect();
    let q: Vec<Float> = omega_rad
        .iter()
        .map(|w| {
            0.5 * air_density * w * w * p.rotor_radius.powi(5) * std::f64::consts::PI * max_cp
                / tsr_opt.powi(3)
        })
        .collect();
    let paero: Vec<Float> = omega_rad
        .iter()
        .zip(q.iter())
        .map(|(w, tq)| tq * w)
        .collect();

    let mut region_2_and_a_half = false;
    let mut rated_power_aero = rated_power_aero_full;
    // np.append copies
    let mut q_lut: Vec<Float> = q;
    let mut omega_array_lut: Vec<Float> = omega_rad;
    let mut paero_lut: Vec<Float> = paero;
    if *paero_lut.last().unwrap() < rated_power_aero {
        region_2_and_a_half = true;
        let q_extra = rated_power_aero / (omega_max * std::f64::consts::PI / 30.0);
        q_lut.push(q_extra);
        omega_array_lut.push(*omega_array_lut.last().unwrap());
        paero_lut.push(rated_power_aero);
    } else {
        rated_power_aero = *paero_lut.last().unwrap();
    }
    let _ = rated_power_aero;

    // omega_lut_pow (RPM) aliases omega_lut_torque in Python; the
    // Region2andAhalf fix increments BOTH lines, i.e. the shared last entry
    // moves by 2e-2 per turbine loop iteration.
    let mut omega_lut_pow: Vec<Float> = omega_array_lut
        .iter()
        .map(|w| w * 30.0 / std::f64::consts::PI)
        .collect();
    // omega_rated / u_rated per (findex, turbine)
    let mut omega_rated = Array2::zeros((n_findex, n_turbines));
    let mut u_rated = Array2::zeros((n_findex, n_turbines));
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let demanded = inp.power_demanded[[fi, ti]];
            let omega_rpm = np_interp(demanded, &paero_lut, &omega_lut_pow);
            omega_rated[[fi, ti]] = omega_rpm * std::f64::consts::PI / 30.0;
            u_rated[[fi, ti]] = (demanded
                / (0.5
                    * air_density
                    * std::f64::consts::PI
                    * p.rotor_radius
                    * p.rotor_radius
                    * max_cp))
                .powf(1.0 / 3.0);
        }
    }

    let mut pitch_out = Array2::zeros((n_findex, n_turbines));
    let mut tsr_out = Array2::zeros((n_findex, n_turbines));

    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let u_v = inp.rotor_effective_velocities[[fi, ti]];
            let yaw = inp.yaw_angles[[fi, ti]];
            let tilt = inp.tilt_angles[[fi, ti]];
            let shear_ij = inp.shear[[fi, ti]];

            let mut tsr_v = if u_v > u_rated[[fi, ti]] {
                omega_rated[[fi, ti]] * p.rotor_radius / u_v * cosd(yaw).powf(0.5)
            } else {
                tsr_opt * cosd(yaw)
            };

            if region_2_and_a_half {
                // Python mutates the aliased arrays twice per iteration
                let last = omega_lut_pow.len() - 1;
                omega_lut_pow[last] += 1e-2; // omega_lut_torque[-1] (alias)
                omega_lut_pow[last] += 1e-2; // omega_lut_pow[-1]
            }

            // --- solve the aero/electric power balance for tsr (get_tsr) ---
            let (tsr_sol, resid, _ok) = fsolve1d(
                &|x| {
                    get_tsr_residual(
                        x,
                        p,
                        air_density,
                        shear_ij,
                        yaw,
                        tilt,
                        u_v,
                        pitch_opt,
                        &omega_lut_pow,
                        &q_lut,
                        inp.cp_surface,
                    )
                },
                tsr_v,
                1.49e-8,
                400,
            );
            let mut tsr_solution = tsr_sol;
            if resid > 10.0 || tsr_solution < 4.0 {
                tsr_solution = 1000.0;
            }
            let mut tsr_o = tsr_solution;
            let omega = tsr_o * u_v / p.rotor_radius;

            let mut pitch0 = pitch_opt;
            if omega >= omega_rated[[fi, ti]] {
                // Region 3: solve for the pitch at the rated omega
                tsr_o = omega_rated[[fi, ti]] * p.rotor_radius / u_v;
                let (pitch_sol, _r, _ok) = fsolve1d(
                    &|x| {
                        get_pitch_residual(
                            x,
                            p,
                            air_density,
                            shear_ij,
                            yaw,
                            tilt,
                            u_v,
                            omega_rated[[fi, ti]],
                            &omega_array_lut,
                            &q_lut,
                            inp.cp_surface,
                        )
                    },
                    u_v,
                    1e-10,
                    2000,
                );
                pitch0 = if pitch_sol < pitch_opt {
                    pitch_opt
                } else {
                    pitch_sol
                };
            }

            pitch_out[[fi, ti]] = pitch0;
            tsr_out[[fi, ti]] = tsr_o;
            let _ = &mut tsr_v;
        }
    }

    Ok((pitch_out, tsr_out))
}

/// Python get_tsr (lines 492-603) — residual in tip-speed ratio.
#[allow(clippy::too_many_arguments)]
fn get_tsr_residual(
    x: Float,
    p: &ControllerDependentParams,
    air_density: Float,
    shear: Float,
    gamma: Float,
    tilt: Float,
    u: Float,
    pitch_in: Float,
    omega_lut_pow: &[Float],
    torque_lut: &[Float],
    cp_surface: &dyn Fn(Float, Float) -> Float,
) -> Float {
    let r = p.rotor_radius;
    let omega = x * u / r;
    let omega_rpm = omega * 30.0 / std::f64::consts::PI;
    // omega_lut_torque = omega_lut_pow * pi/30 (rad/s)
    let omega_lut_torque: Vec<Float> = omega_lut_pow
        .iter()
        .map(|w| w * std::f64::consts::PI / 30.0)
        .collect();
    let torque_nm = np_interp(omega, &omega_lut_torque, torque_lut);

    let theta = pitch_in.to_radians() + p.beta.to_radians();

    // Yawed case
    let mu = (cosd(gamma) * cosd(tilt)).acos();
    let args = BemArgs {
        sigma: p.sigma,
        cd: p.cd,
        cl_alfa: p.cl_alfa,
        gamma,
        delta: tilt,
        k: shear,
        cos_mu: mu.cos(),
        sin_mu: mu.sin(),
        tsr: x,
        theta,
        mu,
    };
    let (ct, _ok) = solve_ct(&args, 0.1);
    let cp = find_cp(&args, ct);

    // Unyawed case
    let mu0 = cosd(tilt).acos();
    let args0 = BemArgs {
        sigma: p.sigma,
        cd: p.cd,
        cl_alfa: p.cl_alfa,
        gamma: 0.0,
        delta: tilt,
        k: shear,
        cos_mu: mu0.cos(),
        sin_mu: mu0.sin(),
        tsr: x,
        theta,
        mu: mu0,
    };
    let (ct0, _ok0) = solve_ct(&args0, 0.1);
    let cp0 = find_cp(&args0, ct0);

    let eta_p = cp / cp0;
    let cp_now = cp_surface(x, pitch_in);
    let cp_g1 = cp_now * eta_p;
    let aero_pow = 0.5 * air_density * (std::f64::consts::PI * r * r) * u.powi(3) * cp_g1;
    let electric_pow = torque_nm * (omega_rpm * std::f64::consts::PI / 30.0);
    aero_pow - electric_pow
}

/// Python get_pitch (lines 606-718) — residual in pitch angle (degrees).
#[allow(clippy::too_many_arguments)]
fn get_pitch_residual(
    x: Float,
    p: &ControllerDependentParams,
    air_density: Float,
    shear: Float,
    gamma: Float,
    tilt: Float,
    u: Float,
    omega_rated: Float,
    omega_array_lut: &[Float], // rad/s
    torque_lut: &[Float],
    cp_surface: &dyn Fn(Float, Float) -> Float,
) -> Float {
    let r = p.rotor_radius;
    let omega_rpm = omega_rated * 30.0 / std::f64::consts::PI;
    let tsr = omega_rated * r / u;

    let theta = x.to_radians() + p.beta.to_radians();

    // torque from the region-2.5 curve (omega_lut_torque * 30/pi converts the
    // rad/s LUT back to rpm, Python line 630-632)
    let omega_lut_rpm: Vec<Float> = omega_array_lut
        .iter()
        .map(|w| w * 30.0 / std::f64::consts::PI)
        .collect();
    let torque_nm = np_interp(omega_rpm, &omega_lut_rpm, torque_lut);

    // Yawed case
    let mu = (cosd(gamma) * cosd(tilt)).acos();
    let args = BemArgs {
        sigma: p.sigma,
        cd: p.cd,
        cl_alfa: p.cl_alfa,
        gamma,
        delta: tilt,
        k: shear,
        cos_mu: mu.cos(),
        sin_mu: mu.sin(),
        tsr,
        theta,
        mu,
    };
    let (ct, _ok) = solve_ct(&args, 0.1);
    let cp = find_cp(&args, ct);

    // Unyawed case
    let mu0 = cosd(tilt).acos();
    let args0 = BemArgs {
        sigma: p.sigma,
        cd: p.cd,
        cl_alfa: p.cl_alfa,
        gamma: 0.0,
        delta: tilt,
        k: shear,
        cos_mu: mu0.cos(),
        sin_mu: mu0.sin(),
        tsr,
        theta,
        mu: mu0,
    };
    let (ct0, _ok0) = solve_ct(&args0, 0.1);
    let cp0 = find_cp(&args0, ct0);

    let eta_p = cp / cp0;
    let cp_now = cp_surface(tsr, x);
    let cp_g1 = cp_now * eta_p;
    let aero_pow = 0.5 * air_density * (std::f64::consts::PI * r * r) * u.powi(3) * cp_g1;
    let electric_pow = torque_nm * (omega_rpm * std::f64::consts::PI / 30.0);
    aero_pow - electric_pow
}

// ===================== the operation model =====================

impl ControllerDependentTurbine {
    /// Shared solve: (tilt_negated, rotor_effective_velocities, shear,
    /// pitch_out, tsr_out).
    ///
    /// `apply_density_correction` mirrors the Python asymmetry between the
    /// two entry points: `power()` corrects the rotor-average velocity with
    /// the ambient air density (`rotor_velocity_air_density_correction`,
    /// line 65-69), while `thrust_coefficient()` reassigns
    /// `air_density = power_thrust_table["ref_air_density"]` BEFORE the
    /// correction (line 280-286), making it a no-op — the control trajectory
    /// there sees the UNCORRECTED rotor-average velocity.
    fn prepare(
        params: &TurbineParameters,
        ctx: &TurbineContext,
        apply_density_correction: bool,
        apply_tilt_interp: bool,
    ) -> crate::Result<(Array2, Array2, Array2, Array2, Array2)> {
        let cd = params
            .controller_dependent
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("controller_dependent_turbine_parameters missing"))?;
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for ControllerDependentTurbine"))?;
        let tilt_angles = ctx.tilt_angles.ok_or_else(|| {
            anyhow::anyhow!("Tilt angles required for ControllerDependentTurbine")
        })?;
        if cd.tsr_lut.len() < 4 || cd.pitch_lut.len() < 4 {
            anyhow::bail!("cp/ct surface too small for cubic interpolation");
        }

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.n_turbines();

        // Python sign convention: negative tilt creates tower clearance
        let mut tilt_negated = -tilt_angles.clone();

        // rotor average + air-density correction (power path only; see
        // `prepare` docs — the Python thrust path skips it)
        let velocities_dyn = ctx.velocities.view().into_dyn().to_owned();
        let cubature_weights = ctx.cubature_weights_1d();
        let rotor_avg = crate::core::rotor_velocity::average_velocity(
            &velocities_dyn,
            None,
            &ctx.average_method,
            cubature_weights.as_ref(),
        )?;
        let mut rotor_eff = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let density_factor = if apply_density_correction {
                (ctx.air_density[fi] / params.ref_air_density).powf(1.0 / 3.0)
            } else {
                1.0
            };
            for ti in 0..n_turbines {
                rotor_eff[[fi, ti]] = rotor_avg[[fi, ti]] * density_factor;
            }
        }

        // Python thrust_coefficient (controller_dependent_operation_model.py:
        // 250-269): after the sign flip, the tilt of turbines requesting the
        // correction is REPLACED by the floating-tilt interpolation at the
        // UNCORRECTED rotor-average velocity (the interp output overrides the
        // negated static tilt). The power path does NOT apply this — its
        // kwargs swallow `tilt_interp` (the Python power() signature has no
        // tilt handling), so `apply_tilt_interp` is false there.
        if apply_tilt_interp {
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    let static_tilt = tilt_negated[[fi, ti]];
                    tilt_negated[[fi, ti]] =
                        ctx.effective_tilt(params, fi, ti, rotor_avg[[fi, ti]], static_tilt);
                }
            }
        }

        let shear = compute_local_vertical_shear(ctx.velocities)?;

        // power_demanded (Python lines 476-489)
        let mut power_demanded = Array2::zeros((n_findex, n_turbines));
        match ctx.power_setpoints {
            None => {
                let base = cd.rated_power_kw * 1000.0 / cd.generator_efficiency;
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        power_demanded[[fi, ti]] = base;
                    }
                }
            }
            Some(sp) => {
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        power_demanded[[fi, ti]] = sp[[fi, ti]] / cd.generator_efficiency;
                    }
                }
            }
        }

        // cubic Cp interpolator
        let cp_lut = &cd.cp_lut;
        let tsr_lut: Vec<Float> = cd.tsr_lut.to_vec();
        let pitch_lut: Vec<Float> = cd.pitch_lut.to_vec();
        let cp_rgi = CubicRgi::new(&tsr_lut, &pitch_lut, cp_lut);

        let inp = TrajectoryInput {
            params: cd,
            rotor_effective_velocities: &rotor_eff,
            yaw_angles,
            tilt_angles: &tilt_negated,
            shear: &shear,
            power_demanded: &power_demanded,
            cp_surface: &|t, p| cp_rgi.eval(t, p),
        };
        let (pitch_out, tsr_out) = control_trajectory(&inp)?;

        Ok((tilt_negated, rotor_eff, shear, pitch_out, tsr_out))
    }

    /// Yawed/unyawed BEM solve pair for each (findex, turbine). Mirrors the
    /// per-quantity Python conventions:
    ///
    /// - `Power` (lines 107-197): both loops gate on `fsolve(..., ier == 1)`
    ///   and fall back to `-1e3` on non-convergence.
    /// - `Thrust` (lines 307-356): `fsolve` return is used REGARDLESS of the
    ///   convergence flag; the yawed Ct is clipped to `[0.0001, 0.9999]`
    ///   (line 330) while the unyawed Ct is used raw (line 356 — the clip
    ///   is commented out in the Python source).
    fn solve_pair(
        params: &TurbineParameters,
        prepared: &(Array2, Array2, Array2, Array2, Array2),
        yaw_angles: &Array2,
        quantity: SolveQuantity,
        x0: Float,
    ) -> crate::Result<(Array2, Array2)> {
        let cd = params.controller_dependent.as_ref().unwrap();
        let (tilt_negated, _rotor_eff, shear, pitch_out, tsr_out) = prepared;
        let (n_findex, n_turbines) = pitch_out.dim();
        let mut yawed = Array2::zeros((n_findex, n_turbines));
        let mut unyawed = Array2::zeros((n_findex, n_turbines));

        let theta_array: Array2 = pitch_out.mapv(|p| (p + cd.beta).to_radians());

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let gamma = yaw_angles[[fi, ti]];
                let delta = tilt_negated[[fi, ti]];
                let k = shear[[fi, ti]];
                let tsr = tsr_out[[fi, ti]];
                let theta = theta_array[[fi, ti]];

                // yawed
                let mu = (cosd(gamma) * cosd(delta)).acos();
                let args = BemArgs {
                    sigma: cd.sigma,
                    cd: cd.cd,
                    cl_alfa: cd.cl_alfa,
                    gamma,
                    delta,
                    k,
                    cos_mu: mu.cos(),
                    sin_mu: mu.sin(),
                    tsr,
                    theta,
                    mu,
                };
                let (ct, ok) = solve_ct(&args, x0);
                yawed[[fi, ti]] = match quantity {
                    SolveQuantity::Power => {
                        if ok {
                            find_cp(&args, ct)
                        } else {
                            -1e3
                        }
                    }
                    SolveQuantity::Thrust => ct.clamp(0.0001, 0.9999),
                };

                // unyawed
                let mu0 = cosd(delta).acos();
                let args0 = BemArgs {
                    sigma: cd.sigma,
                    cd: cd.cd,
                    cl_alfa: cd.cl_alfa,
                    gamma: 0.0,
                    delta,
                    k,
                    cos_mu: mu0.cos(),
                    sin_mu: mu0.sin(),
                    tsr,
                    theta,
                    mu: mu0,
                };
                let (ct0, ok0) = solve_ct(&args0, x0);
                unyawed[[fi, ti]] = match quantity {
                    SolveQuantity::Power => {
                        if ok0 {
                            find_cp(&args0, ct0)
                        } else {
                            -1e3
                        }
                    }
                    // Python uses the raw unyawed Ct (no clip, no ier gate)
                    SolveQuantity::Thrust => ct0,
                };
            }
        }
        Ok((yawed, unyawed))
    }
}

impl OperationModel for ControllerDependentTurbine {
    fn model_name(&self) -> &'static str {
        "controller-dependent"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let cd = params
            .controller_dependent
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("controller_dependent_turbine_parameters missing"))?;
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for ControllerDependentTurbine"))?;

        let prepared = ControllerDependentTurbine::prepare(params, ctx, true, false)?;
        let (p, p0) = ControllerDependentTurbine::solve_pair(
            params,
            &prepared,
            yaw_angles,
            SolveQuantity::Power,
            0.2,
        )?;

        let ratio = &p / &p0;

        // cubic Cp interpolation at (tsr, pitch) * ratio
        let (_tilt, rotor_eff, _shear, pitch_out, tsr_out) = &prepared;
        let (n_findex, n_turbines) = ratio.dim();
        let mut power_coefficient = Array2::zeros((n_findex, n_turbines));
        {
            let tsr_lut: Vec<Float> = cd.tsr_lut.to_vec();
            let pitch_lut: Vec<Float> = cd.pitch_lut.to_vec();
            let rgi = CubicRgi::new(&tsr_lut, &pitch_lut, &cd.cp_lut);
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    let cp_interp = rgi.eval(tsr_out[[fi, ti]], pitch_out[[fi, ti]]);
                    power_coefficient[[fi, ti]] = cp_interp * ratio[[fi, ti]];
                }
            }
        }

        // power = 0.5 * rho_ref * v^3 * pi * R^2 * cp * gen_eff
        let mut power = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                power[[fi, ti]] = 0.5
                    * params.ref_air_density
                    * rotor_eff[[fi, ti]].powi(3)
                    * std::f64::consts::PI
                    * cd.rotor_radius
                    * cd.rotor_radius
                    * power_coefficient[[fi, ti]]
                    * cd.generator_efficiency;
            }
        }

        let rated = cd.rated_power_kw * 1000.0;
        if let Some(mx) = power
            .iter()
            .cloned()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
        {
            if mx > rated * 1.01 {
                log::warn!(
                    "Powers more than 1% above rated detected. Consider checking Cp-Ct data."
                );
            }
        }
        let mut clipped = power;
        for v in clipped.iter_mut() {
            *v = v.clamp(0.0, rated);
        }
        Ok(clipped)
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let cd = params
            .controller_dependent
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("controller_dependent_turbine_parameters missing"))?;
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for ControllerDependentTurbine"))?;

        // Python thrust_coefficient: air_density := ref BEFORE the rotor
        // density correction (a no-op), so the trajectory sees the raw
        // rotor-average velocity (lines 280-298).
        let prepared = ControllerDependentTurbine::prepare(params, ctx, false, true)?;

        // Python clips the yawed Ct (line 330) but not the unyawed one
        // (line 356), and uses the fsolve return regardless of convergence.
        let (ct1, ct0) = ControllerDependentTurbine::solve_pair(
            params,
            &prepared,
            yaw_angles,
            SolveQuantity::Thrust,
            0.2,
        )?;
        let ratio = &ct1 / &ct0;

        let (_tilt, _rotor_eff, _shear, pitch_out, tsr_out) = &prepared;
        let (n_findex, n_turbines) = ratio.dim();
        let mut thrust_coefficient = Array2::zeros((n_findex, n_turbines));
        {
            let tsr_lut: Vec<Float> = cd.tsr_lut.to_vec();
            let pitch_lut: Vec<Float> = cd.pitch_lut.to_vec();
            let rgi = CubicRgi::new(&tsr_lut, &pitch_lut, &cd.ct_lut);
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    let ct_interp = rgi.eval(tsr_out[[fi, ti]], pitch_out[[fi, ti]]);
                    thrust_coefficient[[fi, ti]] = ct_interp * ratio[[fi, ti]];
                }
            }
        }
        Ok(thrust_coefficient)
    }

    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let thrust_coefficients = self.thrust_coefficient(params, ctx)?;
        let tilt_angles = ctx.tilt_angles.ok_or_else(|| {
            anyhow::anyhow!("Tilt angles required for ControllerDependentTurbine")
        })?;

        // Eq. (25a) with yaw zeroed; MU from the ORIGINAL (un-negated) tilt
        let (n_findex, n_turbines) = thrust_coefficients.dim();
        let mut axial_induction = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let ct = thrust_coefficients[[fi, ti]];
                let mu = cosd(tilt_angles[[fi, ti]]).acos();
                let sin_mu = mu.sin();
                let a = 1.0
                    - (1.0 + (1.0 - ct - (1.0 / 16.0) * ct * ct * sin_mu * sin_mu).sqrt())
                        / (2.0 * (1.0 + (1.0 / 16.0) * ct * sin_mu * sin_mu));
                axial_induction[[fi, ti]] = a.clamp(0.0001, 0.9999);
            }
        }
        Ok(axial_induction)
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
