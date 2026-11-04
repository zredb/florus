//! A port of scipy's SLSQP solver.
//!
//! The solver is a translation of scipy's `scipy/optimize/__slsqp.c` (which
//! is itself a translation of the Fortran code of Dieter Kraft, "A software
//! package for sequential quadratic programming", DFVLR-FB 88-28), so that
//! the layout optimization of FLORIS gives the same result as Python's
//! `scipy.optimize.minimize(method="SLSQP")`.
//!
//! The port covers the inequality constrained problems (the equality
//! constraints are not used by the layout optimization and are not
//! supported).
//!
//! Reference: scipy `optimize/__slsqp.c`, `optimize/__nnls.c` and
//! `optimize/_slsqp_py.py`.

use crate::types::Float;

/// The exit mode of the solver (scipy's `exit_modes`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlsqpMode {
    /// Gradient evaluation required (g & a).
    GradientRequired,
    /// Optimization terminated successfully.
    Success,
    /// More equality constraints than independent variables.
    TooManyEqualityConstraints,
    /// More than 3*n iterations in the LSQ subproblem.
    LsqIterationLimit,
    /// Inequality constraints incompatible.
    IncompatibleInequalityConstraints,
    /// Singular matrix E in the LSQ subproblem.
    SingularE,
    /// Singular matrix C in the LSQ subproblem.
    SingularC,
    /// Rank-deficient equality constraint subproblem.
    RankDeficient,
    /// Positive directional derivative for the line search.
    PositiveDirectionalDerivative,
    /// Iteration limit reached.
    IterationLimit,
}

impl SlsqpMode {
    /// The exit code of scipy's C implementation.
    pub fn code(self) -> i32 {
        match self {
            SlsqpMode::GradientRequired => -1,
            SlsqpMode::Success => 0,
            SlsqpMode::TooManyEqualityConstraints => 2,
            SlsqpMode::LsqIterationLimit => 3,
            SlsqpMode::IncompatibleInequalityConstraints => 4,
            SlsqpMode::SingularE => 5,
            SlsqpMode::SingularC => 6,
            SlsqpMode::RankDeficient => 7,
            SlsqpMode::PositiveDirectionalDerivative => 8,
            SlsqpMode::IterationLimit => 9,
        }
    }

    /// The message of scipy's exit modes.
    pub fn message(self) -> &'static str {
        match self {
            SlsqpMode::GradientRequired => "Gradient evaluation required (g & a)",
            SlsqpMode::Success => "Optimization terminated successfully",
            SlsqpMode::TooManyEqualityConstraints => {
                "More equality constraints than independent variables"
            }
            SlsqpMode::LsqIterationLimit => "More than 3*n iterations in LSQ subproblem",
            SlsqpMode::IncompatibleInequalityConstraints => "Inequality constraints incompatible",
            SlsqpMode::SingularE => "Singular matrix E in LSQ subproblem",
            SlsqpMode::SingularC => "Singular matrix C in LSQ subproblem",
            SlsqpMode::RankDeficient => "Rank-deficient equality constraint subproblem",
            SlsqpMode::PositiveDirectionalDerivative => {
                "Positive directional derivative for linesearch"
            }
            SlsqpMode::IterationLimit => "Iteration limit reached",
        }
    }
}

/// The result of the optimization.
#[derive(Debug, Clone)]
pub struct SlsqpResult {
    /// The solution.
    pub x: Vec<Float>,
    /// The objective function at the solution.
    pub fun: Float,
    /// The number of the iterations.
    pub nit: usize,
    /// The number of the function evaluations.
    pub nfev: usize,
    /// The number of the gradient evaluations.
    pub njev: usize,
    /// The exit mode.
    pub mode: SlsqpMode,
    /// The message of the exit mode.
    pub message: String,
}

/// The internal state of the solver (scipy's `SLSQP_vars`).
struct SlsqpVars {
    acc: Float,
    alpha: Float,
    f0: Float,
    gs: Float,
    h1: Float,
    h2: Float,
    h3: Float,
    h4: Float,
    t: Float,
    t0: Float,
    tol: Float,
    iter: i32,
    itermax: i32,
    line: i32,
    mode: i32,
    meq: usize,
    m: usize,
    reset: i32,
}

/// Minimize a function subject to the inequality constraints `c(x) >= 0`.
///
/// * `x0` - the initial guess.
/// * `fun` - the objective function (it is minimized).
/// * `grad` - the gradient of the objective function.
/// * `con` - the constraint functions (the values must be greater than or
///   equal to zero).
/// * `con_jac` - the Jacobian of the constraint functions (m rows, n
///   columns).
/// * `xl`, `xu` - the bounds of the variables.
/// * `maxiter` - the maximum number of the iterations.
/// * `acc` - the requested accuracy (the convergence criterion).
pub fn minimize_slsqp(
    x0: &[Float],
    fun: &mut dyn FnMut(&[Float]) -> Float,
    grad: &mut dyn FnMut(&[Float]) -> Vec<Float>,
    con: &mut dyn FnMut(&[Float]) -> Vec<Float>,
    con_jac: &mut dyn FnMut(&[Float]) -> Vec<Vec<Float>>,
    xl: &[Float],
    xu: &[Float],
    maxiter: usize,
    acc: Float,
) -> SlsqpResult {
    let n = x0.len();
    let d0 = con(x0);
    let m = d0.len();
    let meq = 0usize;

    let mut vars = SlsqpVars {
        acc,
        alpha: 0.0,
        f0: 0.0,
        gs: 0.0,
        h1: 0.0,
        h2: 0.0,
        h3: 0.0,
        h4: 0.0,
        t: 0.0,
        t0: 0.0,
        tol: 0.0,
        iter: 0,
        itermax: maxiter as i32,
        line: 0,
        mode: 0,
        meq,
        m,
        reset: 0,
    };

    let mut sol = x0.to_vec();
    let mut mult = vec![0.0; m + 2 * n + 2];
    let mut indices = vec![0i32; (m + 2 * n + 2).max(1)];

    // The function, the gradient and the constraints
    let mut funx = fun(x0);
    let mut gradx = grad(x0);
    let mut d = d0;
    let mut cmat = {
        let jac = con_jac(x0);
        let mut c = vec![0.0; m.max(1) * n];
        for i in 0..m {
            for j in 0..n {
                c[i + j * m] = jac[i][j];
            }
        }
        c
    };

    let mut nfev = 1;
    let mut njev = 1;

    // The work arrays
    let n2 = n * (n + 1) / 2;
    let mut bfgs = vec![0.0; n2];
    let mut x0v = vec![0.0; n];
    let mut mu = vec![0.0; m];
    let mut s = vec![0.0; n + 1];
    let mut u = vec![0.0; n + 1];
    let mut v = vec![0.0; n + 1];

    // The state machine of the Fortran code
    let mut state = 0u8;
    let mut badlin = false;
    let mut reset = 0;

    loop {
        match state {
            // MODE0: the initialization
            0 => {
                vars.acc = vars.acc.abs();
                vars.tol = 10.0 * vars.acc;
                vars.iter = 0;
                reset = 0;
                state = 10; // RESET_BFGS
            }

            // RESET_BFGS: the reset of the BFGS matrix
            10 => {
                reset += 1;
                if reset > 5 {
                    vars.reset = reset;
                    state = 30; // LABEL255
                    continue;
                }
                vars.reset = reset;
                reset_bfgs(&mut bfgs, n);
                state = 20; // ITER_START
            }

            // ITER_START: the main iteration
            20 => {
                vars.mode = 9;
                if vars.iter >= vars.itermax {
                    break;
                }
                vars.iter += 1;

                // The search direction as the solution of the QP problem
                for i in 0..n {
                    u[i] = -sol[i] + xl[i];
                    v[i] = -sol[i] + xu[i];
                }
                vars.h4 = 1.0;
                let mut mode = lsq(
                    vars.m,
                    vars.meq,
                    n,
                    false,
                    0.0,
                    &bfgs,
                    &gradx,
                    &cmat,
                    &d,
                    &u,
                    &v,
                    &mut s,
                    &mut mult,
                    &mut indices,
                );

                // If the original SQP problem is inconsistent, the problem is
                // augmented with a dummy variable
                badlin = false;
                if mode == 6 && n == vars.meq {
                    mode = 4;
                }
                if mode == 4 {
                    badlin = true;
                    s.iter_mut().for_each(|value| *value = 0.0);
                    vars.h3 = 0.0;
                    let mut rho = 100.0;
                    let mut inconsistent = 0;
                    loop {
                        mode = lsq(
                            vars.m,
                            vars.meq,
                            n,
                            true,
                            rho,
                            &bfgs,
                            &gradx,
                            &cmat,
                            &d,
                            &u,
                            &v,
                            &mut s,
                            &mut mult,
                            &mut indices,
                        );
                        vars.h4 = 1.0 - s[n];
                        if mode == 4 {
                            rho *= 10.0;
                            inconsistent += 1;
                            if inconsistent > 5 {
                                break;
                            }
                            continue;
                        } else if mode != 1 {
                            break;
                        }
                        break;
                    }
                } else if mode != 1 {
                    break;
                }

                // The multipliers for the L1 test
                for i in 0..n {
                    v[i] = gradx[i];
                }
                for i in 0..n {
                    let mut sum = 0.0;
                    for j in 0..vars.m {
                        sum += cmat[j + i * vars.m] * mult[j];
                    }
                    v[i] -= sum;
                }

                vars.f0 = funx;
                x0v.copy_from_slice(&sol[..n]);
                vars.gs = dot(&gradx[..n], &s[..n]);
                vars.h1 = vars.gs.abs();
                vars.h2 = 0.0;
                for j in 0..vars.m {
                    let h3 = if j < vars.meq { d[j] } else { 0.0 };
                    vars.h2 += (-d[j]).max(h3);
                    let h3 = mult[j].abs();
                    mu[j] = h3.max((mu[j] + h3) / 2.0);
                    vars.h1 += h3 * d[j].abs();
                }

                // The convergence check
                vars.mode = 0;
                if vars.h1 < vars.acc && vars.h2 < vars.acc && !badlin && !funx.is_nan() {
                    break;
                }
                vars.h1 = 0.0;
                for j in 0..vars.m {
                    let h3 = if j < vars.meq { d[j] } else { 0.0 };
                    vars.h1 += mu[j] * (-d[j]).max(h3);
                }
                vars.t0 = funx + vars.h1;
                vars.h3 = vars.gs - vars.h1 * vars.h4;
                vars.mode = 8;
                if vars.h3 >= 0.0 {
                    state = 10; // RESET_BFGS
                    continue;
                }

                // The line search with the L1 test function
                vars.line = 0;
                vars.alpha = 1.0;
                state = 11; // LINE_SEARCH
            }

            // LINE_SEARCH
            11 => {
                vars.line += 1;
                vars.h3 *= vars.alpha;
                for i in 0..n {
                    s[i] *= vars.alpha;
                }
                for i in 0..n {
                    sol[i] = x0v[i] + s[i];
                }

                // The function and the constraint evaluation
                funx = fun(&sol);
                nfev += 1;
                d = con(&sol);
                state = 12; // MODE1
            }

            // MODE1: the check of the line search
            12 => {
                vars.t = funx;
                for j in 0..vars.m {
                    let h1 = if j < vars.meq { d[j] } else { 0.0 };
                    vars.t += mu[j] * (-d[j]).max(h1);
                }
                vars.h1 = vars.t - vars.t0;

                if !(vars.h1 <= vars.h3 / 10.0 || vars.line > 10) {
                    vars.alpha = (vars.h3 / (2.0 * (vars.h3 - vars.h1))).max(0.1);
                    state = 11; // LINE_SEARCH
                    continue;
                }

                // The convergence check
                vars.h3 = 0.0;
                for j in 0..vars.m {
                    let h1 = if j < vars.meq { d[j] } else { 0.0 };
                    vars.h3 += (-d[j]).max(h1);
                }
                if ((funx - vars.f0).abs() < vars.acc || nrm2(&s[..n]) < vars.acc)
                    && vars.h3 < vars.acc
                    && !badlin
                    && !funx.is_nan()
                {
                    vars.mode = 0;
                    break;
                } else {
                    vars.mode = -1;
                }

                // The gradient and the constraint Jacobian
                gradx = grad(&sol);
                njev += 1;
                let jac = con_jac(&sol);
                for i in 0..vars.m {
                    for j in 0..n {
                        cmat[i + j * vars.m] = jac[i][j];
                    }
                }
                state = 21; // MODEM1
            }

            // MODEM1: the update of the Cholesky factors of the Hessian
            21 => {
                // u = gradx - C^T mult - v
                for i in 0..n {
                    u[i] = gradx[i];
                }
                for i in 0..n {
                    let mut sum = 0.0;
                    for j in 0..vars.m {
                        sum += cmat[j + i * vars.m] * mult[j];
                    }
                    u[i] -= sum;
                }
                for i in 0..n {
                    u[i] -= v[i];
                }

                // L'*S
                for i in 0..n {
                    v[i] = s[i];
                }
                dtpmv_lt_unit(&bfgs, &mut v[..n], n);

                // D*L'*S
                let mut j = 0usize;
                for i in 0..n {
                    v[i] = bfgs[j] * v[i];
                    j += n - i;
                }

                // L*D*L'*S
                dtpmv_l_unit(&bfgs, &mut v[..n], n);

                vars.h1 = dot(&s[..n], &u[..n]);
                vars.h2 = dot(&s[..n], &v[..n]);
                vars.h3 = 0.2 * vars.h2;
                if vars.h1 < vars.h3 {
                    vars.h4 = (vars.h2 - vars.h3) / (vars.h2 - vars.h1);
                    vars.h1 = vars.h3;
                    let tmp = 1.0 - vars.h4;
                    for i in 0..n {
                        u[i] = vars.h4 * u[i] + tmp * v[i];
                    }
                }

                // The test for a singular update
                if vars.h1 == 0.0 || vars.h2 == 0.0 {
                    state = 10; // RESET_BFGS
                    continue;
                }

                ldl_update(n, &mut bfgs, &mut u[..n], 1.0 / vars.h1, &mut v[..n]);
                ldl_update(n, &mut bfgs, &mut v[..n], -1.0 / vars.h2, &mut u[..n]);

                state = 20; // ITER_START
            }

            // LABEL255: the relaxed convergence check
            30 => {
                vars.h3 = 0.0;
                for j in 0..vars.m {
                    let h1 = if j < vars.meq { d[j] } else { 0.0 };
                    vars.h3 += (-d[j]).max(h1);
                }
                if ((funx - vars.f0).abs() < vars.tol || nrm2(&s[..n]) < vars.tol)
                    && vars.h3 < vars.tol
                    && !badlin
                    && !funx.is_nan()
                {
                    vars.mode = 0;
                } else {
                    vars.mode = 8;
                }
                break;
            }

            _ => break,
        }
    }

    let mode = match vars.mode {
        0 => SlsqpMode::Success,
        -1 => SlsqpMode::GradientRequired,
        2 => SlsqpMode::TooManyEqualityConstraints,
        3 => SlsqpMode::LsqIterationLimit,
        4 => SlsqpMode::IncompatibleInequalityConstraints,
        5 => SlsqpMode::SingularE,
        6 => SlsqpMode::SingularC,
        7 => SlsqpMode::RankDeficient,
        8 => SlsqpMode::PositiveDirectionalDerivative,
        _ => SlsqpMode::IterationLimit,
    };

    SlsqpResult {
        x: sol,
        fun: funx,
        nit: vars.iter.max(0) as usize,
        nfev,
        njev,
        mode,
        message: mode.message().to_string(),
    }
}

/// The reset of the BFGS matrix to the identity (the packed format).
fn reset_bfgs(bfgs: &mut [Float], n: usize) {
    bfgs.iter_mut().for_each(|value| *value = 0.0);
    let mut j = 0usize;
    for i in 0..n {
        bfgs[j] = 1.0;
        j += n - i;
    }
}

/// The solution of the quadratic subproblem (the port of scipy's `lsq`).
///
/// The equality constraints are not supported (`meq` has to be zero).
#[allow(clippy::too_many_arguments)]
fn lsq(
    m: usize,
    meq: usize,
    n: usize,
    augment: bool,
    aug_weight: Float,
    lf: &[Float],
    gradx: &[Float],
    cmat: &[Float],
    d: &[Float],
    xl: &[Float],
    xu: &[Float],
    x: &mut [Float],
    y: &mut [Float],
    jw: &mut [i32],
) -> i32 {
    debug_assert_eq!(meq, 0);
    let orign = n;
    let mineq = m - meq;
    let ld = if augment { n + 1 } else { n };

    // The A matrix and the b vector of the least squares problem
    let mut wa = vec![0.0; ld * ld];
    let mut wb = vec![0.0; ld];

    let mut cursor = 0usize;
    for j in 0..orign {
        let diag = lf[cursor].sqrt();
        cursor += 1;
        wa[j + j * ld] = diag;
        for i in (j + 1)..orign {
            wa[j + i * ld] = lf[cursor] * diag;
            cursor += 1;
        }
    }

    // b = -1/sqrt(d) * inv(Lf) * gradx
    wb[..orign].copy_from_slice(&gradx[..orign]);
    dtpsv_l_unit(lf, &mut wb[..orign], orign);
    cursor = 0;
    for i in 0..orign {
        wb[i] /= -lf[cursor].sqrt();
        cursor += orign - i;
    }

    if augment {
        wa[ld * ld - 1] = aug_weight;
        // The dummy variable is added to the problem
        x[orign] = 1.0;
    }

    // The number of the variables of the least squares problem (the dummy
    // variable is included when the problem is augmented)
    let n = if augment { orign + 1 } else { orign };

    // The bounds of the variables (the dummy variable is bounded by the
    // unit interval)
    let lower = |i: usize| if i < orign { xl[i] } else { 0.0 };
    let upper = |i: usize| if i < orign { xu[i] } else { 1.0 };

    // The inequality constraints: the bounds of the step become the rows of
    // the G matrix
    let mut wg = vec![0.0; (mineq + 2 * n) * ld];
    let mut wh = vec![0.0; mineq + 2 * n];

    let mut nan_count = 0usize;
    let mut nrow = mineq;
    for i in 0..mineq {
        wh[i] = -d[meq + i];
    }
    for i in 0..n {
        if lower(i).is_nan() {
            nan_count += 1;
        } else {
            wh[nrow] = lower(i);
            nrow += 1;
        }
    }
    for i in 0..n {
        if upper(i).is_nan() {
            nan_count += 1;
        } else {
            wh[nrow] = -upper(i);
            nrow += 1;
        }
    }
    let n_wg_rows = mineq + 2 * n - nan_count;

    // The top part of the G matrix (the constraint normals)
    for j in 0..orign {
        for i in 0..mineq {
            wg[i + j * n_wg_rows] = cmat[meq + i + j * m];
        }
    }
    if augment {
        for i in 0..mineq {
            wg[i + orign * n_wg_rows] = (-d[meq + i]).max(0.0);
        }
    }

    // The rows of the bounds
    let mut nrow = mineq;
    for i in 0..n {
        if !lower(i).is_nan() {
            wg[nrow + i * n_wg_rows] = 1.0;
            nrow += 1;
        }
    }
    for i in 0..n {
        if !upper(i).is_nan() {
            wg[nrow + i * n_wg_rows] = -1.0;
            nrow += 1;
        }
    }

    // The solution of the least squares problem
    let mut multipliers = vec![0.0; n_wg_rows];
    let mode = lsei(
        ld,
        meq,
        n_wg_rows,
        n,
        &wa,
        &wb,
        &mut wg,
        &mut wh,
        x,
        jw,
        &mut multipliers,
    );

    if mode == 1 {
        // The Lagrange multipliers of the inequality constraints
        for i in 0..meq {
            y[i] = 0.0;
        }
        for i in 0..mineq {
            y[meq + i] = multipliers[i];
        }
        for i in 0..2 * n {
            y[m + i] = Float::NAN;
        }
    }

    // The solution is clamped to the bounds
    for i in 0..n {
        if !lower(i).is_nan() && x[i] < lower(i) {
            x[i] = lower(i);
        } else if !upper(i).is_nan() && x[i] > upper(i) {
            x[i] = upper(i);
        }
    }

    mode
}

/// The solution of the equality and inequality constrained least squares
/// problem (the port of scipy's `lsei`).
///
/// The equality constraints are not supported (`me` has to be zero).
#[allow(clippy::too_many_arguments)]
fn lsei(
    ma: usize,
    me: usize,
    mg: usize,
    n: usize,
    a: &[Float],
    b: &[Float],
    g: &mut [Float],
    h: &mut [Float],
    x: &mut [Float],
    jw: &mut [i32],
    multipliers: &mut [Float],
) -> i32 {
    debug_assert_eq!(me, 0);
    let nvars = n - me;

    for i in 0..n {
        x[i] = 0.0;
    }
    if me > n {
        return 2;
    }

    // The modified right hand side (the equality part is empty)
    let wb = b.to_vec();

    // The transformed A2 and G2
    let mut a2 = vec![0.0; ma * nvars];
    let mut g2 = vec![0.0; mg * nvars];
    for j in 0..nvars {
        for i in 0..ma {
            a2[i + j * ma] = a[i + j * ma];
        }
        for i in 0..mg {
            g2[i + j * mg] = g[i + j * mg];
        }
    }

    if mg == 0 {
        // The least squares problem without the inequality constraints
        if ma < nvars {
            // The rank of A is smaller than the number of the variables
            return 7;
        }
        let mut a_qr = a2.clone();
        let (qr, tau) = geqr2(ma, nvars, &mut a_qr);
        let mut wb_qr = wb.clone();
        apply_qt(ma, 1, ma.min(nvars), &qr, &tau, &mut wb_qr, ma);
        let mut x_qr = vec![0.0; nvars];
        for i in 0..nvars {
            x_qr[i] = wb_qr[i];
        }
        trsv_upper(&qr, &mut x_qr, ma, nvars);
        x[..nvars].copy_from_slice(&x_qr[..nvars]);
        return 1;
    }

    // The inequality constrained least squares problem
    lsi(ma, mg, nvars, &a2, &wb, &g2, h, x, jw, multipliers)
}

/// The solution of the inequality constrained least squares problem (the port
/// of scipy's `lsi`).
#[allow(clippy::too_many_arguments)]
fn lsi(
    ma: usize,
    mg: usize,
    n: usize,
    a: &[Float],
    b: &[Float],
    g: &[Float],
    h: &mut [Float],
    x: &mut [Float],
    jw: &mut [i32],
    multipliers: &mut [Float],
) -> i32 {
    const EPSMACH: Float = 2.220446049250313e-16;

    // The QR decomposition of A and the application of Q^T to b
    let mut a_qr = a.to_vec();
    let (qr, tau) = geqr2(ma, n, &mut a_qr);
    let mut b_work = b.to_vec();
    apply_qt(ma, 1, ma.min(n), &qr, &tau, &mut b_work, ma);

    // The check for a rank deficiency
    let tmp_int = ma.min(n);
    for i in 0..tmp_int {
        if !(qr[i + i * ma].abs() >= EPSMACH) {
            return 5;
        }
    }

    // Solve X*R = G
    let mut g_work = g.to_vec();
    trsm_right_upper(&qr, &mut g_work, ma, mg, n);

    // h = h - X*b
    for i in 0..mg {
        let mut sum = 0.0;
        for j in 0..n {
            sum += g_work[i + j * mg] * b_work[j];
        }
        h[i] -= sum;
    }

    // The least distance problem
    let mode = ldp(mg, n, &g_work, h, x, jw, multipliers);
    if mode != 1 {
        return mode;
    }

    // The solution of the original problem
    for i in 0..n {
        x[i] += b_work[i];
    }
    trsv_upper(&qr, x, ma, n);

    mode
}

/// The solution of the least distance problem (the port of scipy's `ldp`).
fn ldp(
    m: usize,
    n: usize,
    g: &[Float],
    h: &[Float],
    x: &mut [Float],
    _indices: &mut [i32],
    multipliers: &mut [Float],
) -> i32 {
    if n == 0 {
        return 2;
    }
    for i in 0..n {
        x[i] = 0.0;
    }
    if m == 0 {
        return 1;
    }

    // The dual problem: the transposed and stacked G and h
    let mut a = vec![0.0; (n + 1) * m];
    for j in 0..m {
        for i in 0..n {
            a[i + j * (n + 1)] = g[j + i * m];
        }
        a[n + j * (n + 1)] = h[j];
    }
    let mut b = vec![0.0; n + 1];
    b[n] = 1.0;

    let mut y = vec![0.0; m];
    let mut w = vec![0.0; m];
    let mut zz = vec![0.0; n + 1];
    let mut idx = vec![0i32; m];
    let mut rnorm = 0.0;
    let mode = nnls(
        n + 1,
        m,
        &mut a,
        &mut b,
        &mut y,
        &mut w,
        &mut zz,
        &mut idx,
        3 * m,
        &mut rnorm,
    );
    if mode != 1 {
        return mode;
    }
    if rnorm <= 0.0 {
        return 4;
    }

    // The primal problem
    let fac = 1.0 - dot(&h[..m], &y[..m]);
    if !((1.0 + fac) - 1.0 > 0.0) {
        return 4;
    }
    let fac = 1.0 / fac;
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..m {
            sum += g[j + i * m] * y[j];
        }
        x[i] = fac * sum;
    }

    // The Lagrange multipliers of the primal problem
    for i in 0..m {
        multipliers[i] = fac * y[i];
    }

    mode
}

/// The non-negative least squares solution (the port of scipy's `__nnls`,
/// the algorithm of Lawson and Hanson).
#[allow(clippy::too_many_arguments)]
fn nnls(
    m: usize,
    n: usize,
    a: &mut [Float],
    b: &mut [Float],
    x: &mut [Float],
    w: &mut [Float],
    zz: &mut [Float],
    indices: &mut [i32],
    maxiter: usize,
    rnorm: &mut Float,
) -> i32 {
    if m == 0 || n == 0 {
        return 2;
    }

    for i in 0..n {
        indices[i] = i as i32;
        x[i] = 0.0;
    }

    let mut indz = 0usize;
    let mut iteration = 0usize;
    let mut iz;

    while indz < m.min(n) {
        // The dual vector components of the set Z
        for i in indz..n {
            let j = indices[i] as usize;
            let mut sum = 0.0;
            for k in indz..m {
                sum += a[k + j * m] * b[k];
            }
            w[j] = sum;
        }

        // The next linearly independent column
        let mut pivot;
        let mut tau;
        loop {
            let mut wmax = 0.0;
            let mut izmax = indz;
            for k in indz..n {
                let j = indices[k] as usize;
                if w[j] > wmax {
                    wmax = w[j];
                    izmax = k;
                }
            }
            if wmax <= 0.0 {
                // The KKT certificate: the solution is optimal
                *rnorm = if indz < m { nrm2(&b[indz..m]) } else { 0.0 };
                return 1;
            }
            iz = izmax;
            let j = indices[iz] as usize;

            // The Householder reflector of the column j (the non-negative
            // beta variant of LAPACK's `dlarfgp`)
            let mut column = vec![0.0; m - indz - 1];
            for i in 0..(m - indz - 1) {
                column[i] = a[indz + 1 + i + j * m];
            }
            let (alpha, t) = larfgp(m - indz, a[indz + j * m], &mut column);
            pivot = alpha;
            tau = t;
            for i in 0..(m - indz - 1) {
                a[indz + 1 + i + j * m] = column[i];
            }

            // The norm of the column above the pivot
            let unorm = if indz > 0 {
                nrm2(&a[j * m..j * m + indz])
            } else {
                0.0
            };
            let spacing = if unorm > 0.0 {
                next_after(unorm, 2.0 * unorm) - unorm
            } else {
                0.0
            };

            if pivot.abs() > 100.0 * spacing {
                // The column is sufficiently independent
                for i in 0..m {
                    zz[i] = b[i];
                }
                let pivot2 = a[indz + j * m];
                a[indz + j * m] = 1.0;
                let mut v = vec![0.0; m - indz];
                v[0] = 1.0;
                for i in 1..(m - indz) {
                    v[i] = a[indz + i + j * m];
                }
                apply_reflector_left(&v, tau, &mut zz[indz..m]);
                let ztest = zz[indz] / pivot;
                if ztest > 0.0 {
                    break;
                } else {
                    a[indz + j * m] = pivot2;
                }
            }
            // The column is rejected
            w[j] = 0.0;
        }

        // The index j is moved from the set Z to the set P
        let j = indices[iz] as usize;
        for i in 0..m {
            b[i] = zz[i];
        }
        indices[iz] = indices[indz];
        indices[indz] = j as i32;
        indz += 1;

        // The Householder transformation of the remaining columns
        if indz < n {
            let mut v = vec![0.0; m - indz + 1];
            v[0] = 1.0;
            for i in 1..(m - indz + 1) {
                v[i] = a[indz - 1 + i + j * m];
            }
            for k in indz..n {
                let jk = indices[k] as usize;
                let mut col = vec![0.0; m - indz + 1];
                for i in 0..(m - indz + 1) {
                    col[i] = a[indz - 1 + i + jk * m];
                }
                apply_reflector_left(&v, tau, &mut col);
                for i in 0..(m - indz + 1) {
                    a[indz - 1 + i + jk * m] = col[i];
                }
            }
        }
        // The pivot element is restored and the subdiagonal is zeroed
        a[indz - 1 + j * m] = pivot;
        if indz < m {
            for i in indz..m {
                a[j * m + i] = 0.0;
            }
        }
        w[j] = 0.0;

        // The permuted triangular system is solved
        for k in 0..indz {
            let ip = indz - 1 - k;
            if k != 0 {
                let jk = indices[ip + 1] as usize;
                for i in 0..=ip {
                    zz[i] -= a[i + jk * m] * zz[ip + 1];
                }
            }
            let jk = indices[ip] as usize;
            zz[ip] /= a[ip + jk * m];
        }

        // The inner loop
        let mut jj = 0usize;
        loop {
            iteration += 1;
            if iteration >= maxiter {
                return 3;
            }

            // The feasibility of the new coefficients
            let mut alpha = 2.0;
            for ip in 0..indz {
                let k = indices[ip] as usize;
                if zz[ip] <= 0.0 {
                    let t = -x[k] / (zz[ip] - x[k]);
                    if alpha > t {
                        alpha = t;
                        jj = ip;
                    }
                }
            }
            if alpha == 2.0 {
                break;
            }

            // The interpolation between the old x and zz
            for ip in 0..indz {
                let k = indices[ip] as usize;
                x[k] += alpha * (zz[ip] - x[k]);
            }

            // The coefficient i is moved from the set P to the set Z
            let mut i = indices[jj] as usize;
            loop {
                x[i] = 0.0;

                if jj != indz - 1 {
                    jj += 1;
                    for j in jj..indz {
                        let ii = indices[j] as usize;
                        indices[j - 1] = ii as i32;
                        let (cc, ss, r) = dartg(a[j - 1 + ii * m], a[j + ii * m]);
                        a[j - 1 + ii * m] = r;
                        a[j + ii * m] = 0.0;
                        for k in 0..n {
                            if k != ii {
                                let tmp = a[j - 1 + k * m];
                                a[j - 1 + k * m] = cc * tmp + ss * a[j + k * m];
                                a[j + k * m] = -ss * tmp + cc * a[j + k * m];
                            }
                        }
                        let tmp = b[j - 1];
                        b[j - 1] = cc * tmp + ss * b[j];
                        b[j] = -ss * tmp + cc * b[j];
                    }
                }
                indz -= 1;
                indices[indz] = i as i32;

                // The remaining coefficients of the set P
                let mut nobreak = false;
                jj = 0;
                while jj < indz {
                    i = indices[jj] as usize;
                    if x[i] <= 0.0 {
                        break;
                    }
                    if jj == indz - 1 {
                        nobreak = true;
                    }
                    jj += 1;
                }
                if nobreak {
                    break;
                }
            }

            for i in 0..m {
                zz[i] = b[i];
            }
            for k in 0..indz {
                let ip = indz - 1 - k;
                if k != 0 {
                    let jk = indices[ip + 1] as usize;
                    for i in 0..=ip {
                        zz[i] -= a[i + jk * m] * zz[ip + 1];
                    }
                }
                let jk = indices[ip] as usize;
                zz[ip] /= a[ip + jk * m];
            }
        }

        // The outer loop
        for k in 0..indz {
            let i = indices[k] as usize;
            x[i] = zz[k];
        }
    }

    *rnorm = if indz < m { nrm2(&b[indz..m]) } else { 0.0 };
    1
}

/// The update of the LDL' factors of a matrix by a rank-one matrix
/// `sigma * z * z'` (the port of scipy's `ldl_update`).
fn ldl_update(n: usize, a: &mut [Float], z: &mut [Float], sigma: Float, w: &mut [Float]) {
    const EPSMACH: Float = 2.220446049250313e-16;
    if sigma == 0.0 {
        return;
    }
    let mut t = 1.0 / sigma;

    if sigma <= 0.0 {
        // The negative update
        w[..n].copy_from_slice(&z[..n]);
        let mut ij = 0usize;
        for i in 0..n {
            let v = w[i];
            t += v * v / a[ij];
            for j in (i + 1)..n {
                ij += 1;
                w[j] -= v * a[ij];
            }
            ij += 1;
        }
        if t >= 0.0 {
            t = EPSMACH / sigma;
        }
        for i in 0..n {
            let j = n - i - 1;
            ij -= i + 1;
            let u = w[j];
            w[j] = t;
            t -= u * u / a[ij];
        }
    }

    // The positive update
    let mut ij = 0usize;
    for i in 0..n {
        let v = z[i];
        let delta = v / a[ij];
        let tp = if sigma < 0.0 { w[i] } else { t + delta * v };
        let alpha = tp / t;
        a[ij] *= alpha;
        if i == n - 1 {
            return;
        }
        let beta = delta / tp;
        if alpha <= 4.0 {
            for j in (i + 1)..n {
                ij += 1;
                z[j] -= v * a[ij];
                a[ij] += beta * z[j];
            }
        } else {
            let gamma = t / tp;
            for j in (i + 1)..n {
                ij += 1;
                let u = a[ij];
                a[ij] = gamma * u + beta * z[j];
                z[j] -= v * u;
            }
        }
        ij += 1;
        t = tp;
    }
}

// ---------------------------------------------------------------------------
// The linear algebra helpers
// ---------------------------------------------------------------------------

/// The dot product of two vectors.
fn dot(x: &[Float], y: &[Float]) -> Float {
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum()
}

/// The Euclidean norm of a vector (the BLAS `dnrm2`).
fn nrm2(x: &[Float]) -> Float {
    let mut scale = 0.0;
    let mut ssq = 1.0;
    for &v in x {
        if v != 0.0 {
            let av = v.abs();
            if scale < av {
                ssq = 1.0 + ssq * (scale / av) * (scale / av);
                scale = av;
            } else {
                ssq += (av / scale) * (av / scale);
            }
        }
    }
    scale * ssq.sqrt()
}

/// The next representable number after `x` (towards `y`).
fn next_after(x: Float, y: Float) -> Float {
    if x == y {
        return y;
    }
    if x.is_nan() || y.is_nan() {
        return Float::NAN;
    }
    let bits = x.to_bits();
    let next = if (x < y) == (x >= 0.0) {
        bits + 1
    } else {
        bits - 1
    };
    Float::from_bits(next)
}

/// The index of the element (i, j) of a lower triangular matrix in the
/// LAPACK packed format (the columns are stored one after the other).
fn packed_index(i: usize, j: usize, n: usize) -> usize {
    j * n - j * (j + 1) / 2 + i
}

/// `v = L^T v`, where L is the unit lower triangular matrix whose strictly
/// lower part is stored in the packed array `ap` (LAPACK's `dtpmv` with
/// DIAG='U').
fn dtpmv_lt_unit(ap: &[Float], v: &mut [Float], n: usize) {
    for j in 0..n {
        let mut sum = v[j];
        for i in (j + 1)..n {
            sum += ap[packed_index(i, j, n)] * v[i];
        }
        v[j] = sum;
    }
}

/// `v = L v`, where L is the unit lower triangular matrix whose strictly
/// lower part is stored in the packed array `ap` (LAPACK's `dtpmv` with
/// DIAG='U').
fn dtpmv_l_unit(ap: &[Float], v: &mut [Float], n: usize) {
    for i in 0..n {
        let mut sum = v[i];
        for j in 0..i {
            sum += ap[packed_index(i, j, n)] * v[j];
        }
        v[i] = sum;
    }
}

/// `x = L^-1 b`, where L is the unit lower triangular matrix whose strictly
/// lower part is stored in the packed array `ap` (LAPACK's `dtpsv` with
/// DIAG='U').
fn dtpsv_l_unit(ap: &[Float], b: &mut [Float], n: usize) {
    for i in 0..n {
        let mut sum = b[i];
        for j in 0..i {
            sum -= ap[packed_index(i, j, n)] * b[j];
        }
        b[i] = sum;
    }
}

/// The Householder reflector of a vector (LAPACK's `dlarfg`).
///
/// It returns the new value of the first element of the vector (the
/// reflector `v` has the implicit first element 1) and `tau`.
fn larfg(n: usize, alpha: Float, x: &mut [Float]) -> (Float, Float) {
    if n <= 1 {
        return (alpha, 0.0);
    }
    let xnorm = nrm2(x);
    if xnorm == 0.0 {
        return (alpha, 0.0);
    }
    let beta = -alpha.hypot(xnorm).copysign(alpha);
    let tau = (beta - alpha) / beta;
    let scale = 1.0 / (alpha - beta);
    for value in x.iter_mut() {
        *value *= scale;
    }
    (beta, tau)
}

/// The Givens rotation (LAPACK's `dlartgp`).
///
/// It returns the cosine and the sine of the rotation and the resulting `r`
/// (which is never negative).
fn dartg(f: Float, g: Float) -> (Float, Float, Float) {
    if g == 0.0 {
        return (if f >= 0.0 { 1.0 } else { -1.0 }, 0.0, f.abs());
    }
    if f == 0.0 {
        return (0.0, if g >= 0.0 { 1.0 } else { -1.0 }, g.abs());
    }
    let r = (f * f + g * g).sqrt();
    (f / r, g / r, r)
}

/// The Householder reflector of a vector with a non-negative beta (LAPACK's
/// `dlarfgp`, which the non-negative least squares solver uses).
///
/// It returns the new value of the first element of the vector (the reflector
/// `v` has the implicit first element 1) and `tau`.
fn larfgp(n: usize, alpha: Float, x: &mut [Float]) -> (Float, Float) {
    if n == 0 {
        return (alpha, 0.0);
    }
    let eps = Float::EPSILON;
    let xnorm = nrm2(x);
    if xnorm <= eps * alpha.abs() {
        if alpha >= 0.0 {
            return (alpha, 0.0);
        }
        // The reflector is the identity with the sign of the first element
        // flipped
        for value in x.iter_mut() {
            *value = 0.0;
        }
        return (-alpha, 2.0);
    }

    let mut beta = alpha.hypot(xnorm).copysign(alpha);
    let mut alpha = alpha + beta;
    let tau;
    if beta < 0.0 {
        beta = -beta;
        tau = -alpha / beta;
    } else {
        alpha = xnorm * (xnorm / alpha);
        tau = alpha / beta;
        alpha = -alpha;
    }
    for value in x.iter_mut() {
        *value /= alpha;
    }
    (beta, tau)
}

/// The application of a Householder reflector to the columns of a matrix
/// (LAPACK's `dlarf` with SIDE='L').
fn apply_reflector_left(v: &[Float], tau: Float, c: &mut [Float]) {
    if tau == 0.0 {
        return;
    }
    let m = v.len();
    let n = c.len() / m;
    let mut work = vec![0.0; n];
    for j in 0..n {
        let mut sum = 0.0;
        for i in 0..m {
            sum += c[i + j * m] * v[i];
        }
        work[j] = sum * tau;
    }
    for j in 0..n {
        for i in 0..m {
            c[i + j * m] -= work[j] * v[i];
        }
    }
}

/// The QR decomposition of a matrix (LAPACK's `dgeqr2`).
///
/// The matrix is overwritten with the upper triangular factor R and the
/// Householder reflectors; the reflectors are returned in `tau`.
fn geqr2(m: usize, n: usize, a: &mut [Float]) -> (Vec<Float>, Vec<Float>) {
    let k = m.min(n);
    let mut tau = vec![0.0; k];
    for i in 0..k {
        // The reflector of the column i
        let mut column = vec![0.0; m - i - 1];
        for j in 0..(m - i - 1) {
            column[j] = a[i + 1 + j + i * m];
        }
        let (alpha, t) = larfg(m - i, a[i + i * m], &mut column);
        a[i + i * m] = alpha;
        for j in 0..(m - i - 1) {
            a[i + 1 + j + i * m] = column[j];
        }
        tau[i] = t;

        // The reflector is applied to the remaining columns
        if i + 1 < n {
            let mut v = vec![0.0; m - i];
            v[0] = 1.0;
            for j in 1..(m - i) {
                v[j] = a[i + j + i * m];
            }
            let mut c = vec![0.0; (m - i) * (n - i - 1)];
            for j in 0..(n - i - 1) {
                for l in 0..(m - i) {
                    c[l + j * (m - i)] = a[i + l + (i + 1 + j) * m];
                }
            }
            apply_reflector_left(&v, t, &mut c);
            for j in 0..(n - i - 1) {
                for l in 0..(m - i) {
                    a[i + l + (i + 1 + j) * m] = c[l + j * (m - i)];
                }
            }
        }
    }
    (a.to_vec(), tau)
}

/// The application of Q^T to a matrix (LAPACK's `dorm2r` with SIDE='L' and
/// TRANS='T').
///
/// `m` and `k` are the dimensions of A, `ncols` is the number of the columns
/// of C and `ldc` its leading dimension.
fn apply_qt(
    m: usize,
    ncols: usize,
    k: usize,
    a: &[Float],
    tau: &[Float],
    c: &mut [Float],
    ldc: usize,
) {
    for i in (0..k).rev() {
        let mut v = vec![0.0; m - i];
        v[0] = 1.0;
        for j in 1..(m - i) {
            v[j] = a[i + j + i * ldc];
        }
        // The reflector is applied to the rows i..m of C
        let mut col = vec![0.0; (m - i) * ncols];
        for j in 0..ncols {
            for l in 0..(m - i) {
                col[l + j * (m - i)] = c[i + l + j * ldc];
            }
        }
        apply_reflector_left(&v, tau[i], &mut col);
        for j in 0..ncols {
            for l in 0..(m - i) {
                c[i + l + j * ldc] = col[l + j * (m - i)];
            }
        }
    }
}

/// The solution of the upper triangular system `A x = b` (LAPACK's `dtrsv`
/// with UPLO='U', TRANS='N').
fn trsv_upper(a: &[Float], b: &mut [Float], lda: usize, n: usize) {
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= a[i + j * lda] * b[j];
        }
        b[i] = sum / a[i + i * lda];
    }
}

/// The solution of the right sided upper triangular system `X A = B` (LAPACK's
/// `dtrsm` with SIDE='R', UPLO='U', TRANSA='N').
///
/// `B` is m x n and `A` is the n x n upper triangular matrix. The columns of
/// `X` are computed from the previous ones (the forward substitution over the
/// columns of `B`).
fn trsm_right_upper(a: &[Float], b: &mut [Float], lda: usize, m: usize, n: usize) {
    for j in 0..n {
        for k in 0..m {
            let mut sum = b[k + j * m];
            for l in 0..j {
                sum -= b[k + l * m] * a[l + j * lda];
            }
            b[k + j * m] = sum / a[j + j * lda];
        }
    }
}

#[cfg(test)]
#[allow(clippy::approx_constant)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// The reference values are from scipy:
    ///   minimize(f, [0, 0], method="SLSQP", bounds=[(0, 2), (0, 2)],
    ///            constraints=[{"type": "ineq", "fun": c}], options={...})
    #[test]
    fn test_slsqp_matches_scipy() {
        let fun = |x: &[Float]| (x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2);
        let grad = |x: &[Float]| vec![2.0 * (x[0] - 1.0), 2.0 * (x[1] - 2.0)];
        let con = |x: &[Float]| vec![1.5 - x[0] - x[1], x[0], 2.0 - x[0], x[1], 2.0 - x[1]];
        let con_jac = |_x: &[Float]| {
            vec![
                vec![-1.0, -1.0],
                vec![1.0, 0.0],
                vec![-1.0, 0.0],
                vec![0.0, 1.0],
                vec![0.0, -1.0],
            ]
        };

        let mut fun = fun;
        let mut grad = grad;
        let mut con = con;
        let mut con_jac = con_jac;
        let result = minimize_slsqp(
            &[0.0, 0.0],
            &mut fun,
            &mut grad,
            &mut con,
            &mut con_jac,
            &[0.0, 0.0],
            &[2.0, 2.0],
            100,
            1e-9,
        );

        assert_eq!(result.mode, SlsqpMode::Success);
        assert_relative_eq!(result.x[0], 0.25, epsilon = 1e-9);
        assert_relative_eq!(result.x[1], 1.25, epsilon = 1e-9);
        assert_relative_eq!(result.fun, 1.125, epsilon = 1e-9);
        assert_eq!(result.nit, 4);
        assert_eq!(result.nfev, 4);
        assert_eq!(result.njev, 4);
    }

    /// The reference values are from scipy's `scipy.optimize.nnls` (the
    /// algorithm of Lawson and Hanson).
    #[test]
    fn test_nnls_matches_scipy() {
        let mut a = vec![
            1.78862847,
            -0.27738820,
            -0.04381817,
            0.88131804,
            -0.54535995,
            0.43650985,
            -0.35475898,
            -0.47721803,
            1.70957306,
            -1.54647732,
            0.09649747,
            -0.08274148,
            -1.31386475,
            0.05003364,
            0.98236743,
            -1.86349270,
            -0.62700068,
            0.88462238,
            -0.40467741,
            -1.10106763,
        ];
        let mut b = vec![
            -1.18504653,
            -0.20564990,
            1.48614836,
            0.23671627,
            -1.02378514,
        ];
        let mut x = vec![0.0; 4];
        let mut w = vec![0.0; 4];
        let mut zz = vec![0.0; 5];
        let mut indices = vec![0i32; 4];
        let mut rnorm = 0.0;

        let mode = nnls(
            5,
            4,
            &mut a,
            &mut b,
            &mut x,
            &mut w,
            &mut zz,
            &mut indices,
            15,
            &mut rnorm,
        );
        assert_eq!(mode, 1);
        assert_relative_eq!(x[0], 0.3222906, epsilon = 1e-6);
        assert_relative_eq!(x[1], 0.0, epsilon = 1e-9);
        assert_relative_eq!(x[2], 0.0, epsilon = 1e-9);
        assert_relative_eq!(x[3], 0.93544373, epsilon = 1e-6);
        assert_relative_eq!(rnorm, 0.9038348271765462, epsilon = 1e-7);
    }

    /// The reconstruction of a matrix from its LDL' factorization (the packed
    /// format of the BFGS matrix).
    fn unpack_ldl(ap: &[Float], n: usize) -> Vec<Float> {
        // B = L D L' with L unit lower triangular
        let mut b = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    let l_ik = if k > i {
                        0.0
                    } else if k == i {
                        1.0
                    } else {
                        ap[packed_index(i, k, n)]
                    };
                    let l_jk = if k > j {
                        0.0
                    } else if k == j {
                        1.0
                    } else {
                        ap[packed_index(j, k, n)]
                    };
                    sum += l_ik * ap[packed_index(k, k, n)] * l_jk;
                }
                b[i + j * n] = sum;
            }
        }
        b
    }

    /// The update of the LDL' factors has to add the rank-one matrix
    /// `sigma * z * z'` to the factorization.
    #[test]
    fn test_ldl_update_adds_the_rank_one_matrix() {
        let n = 3;
        // The LDL' factorization of a symmetric positive definite matrix
        // (the packed format: the strictly lower part of L and the diagonal D)
        let mut ap = vec![0.0; n * (n + 1) / 2];
        ap[0] = 2.0; // D(0, 0)
        ap[1] = 0.3; // L(1, 0)
        ap[2] = -0.2; // L(2, 0)
        ap[3] = 1.5; // D(1, 1)
        ap[4] = 0.4; // L(2, 1)
        ap[5] = 1.1; // D(2, 2)

        let before = unpack_ldl(&ap, n);
        let z = vec![0.5, -1.0, 2.0];
        let mut z_work = z.clone();
        let mut w = vec![0.0; n];
        let sigma = 0.7;
        ldl_update(n, &mut ap, &mut z_work, sigma, &mut w);
        let after = unpack_ldl(&ap, n);

        for i in 0..n {
            for j in 0..=i {
                let expected = before[i + j * n] + sigma * z[i] * z[j];
                assert_relative_eq!(after[i + j * n], expected, epsilon = 1e-12);
            }
        }
    }

    /// The solution of the quadratic subproblem (the port of scipy's `lsq`).
    ///
    /// The reference value is from scipy:
    ///   minimize(0.5 s' B s + g' s, [0, 0], method="SLSQP", ...)
    /// with B = [[2, 0.5], [0.5, 1]] and g = [-2, -1].
    #[test]
    fn test_lsq_matches_scipy() {
        let n = 2;
        let m = 5;

        // The LDL' factorization of B = [[2, 0.5], [0.5, 1]] (the packed
        // format: the strictly lower part of L and the diagonal of D)
        let lf = vec![2.0, 0.25, 0.875];
        let gradx = vec![-2.0, -1.0];
        let cmat = vec![
            // The constraint normals (m rows, n columns)
            -1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0,
        ];
        let d = vec![1.5, 0.0, 2.0, 0.0, 2.0];
        let xl = vec![0.0, 0.0];
        let xu = vec![2.0, 2.0];

        let mut x = vec![0.0; n];
        let mut y = vec![0.0; m + 2 * n + 2];
        let mut jw = vec![0i32; m + 2 * n + 2];
        let mode = lsq(
            m, 0, n, false, 0.0, &lf, &gradx, &cmat, &d, &xl, &xu, &mut x, &mut y, &mut jw,
        );

        assert_eq!(mode, 1);
        assert_relative_eq!(x[0], 0.85714284, epsilon = 1e-6);
        assert_relative_eq!(x[1], 0.57142861, epsilon = 1e-6);
    }

    /// The solution of the quadratic subproblem with an active constraint.
    ///
    /// The QP is `min 0.5 s' B s + g' s` with B = [[1, 0], [0, 2]] and
    /// g = [-2, -1], subject to `s0 + s1 <= 0` and `0 <= s0 <= 2`,
    /// `-1.5 <= s1 <= 0.5`. The solution is s = [1/3, -1/3].
    #[test]
    fn test_lsq_with_an_active_constraint() {
        let n = 2;
        let m = 5;

        // The LDL' factorization of B = [[1, 0], [0, 2]]
        let lf = vec![1.0, 0.0, 2.0];
        let gradx = vec![-2.0, -1.0];
        let cmat = vec![
            // The constraint normals (m rows, n columns)
            -1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0,
        ];
        // The constraint values at the current point
        let d = vec![0.0, 0.0, 2.0, 1.5, 0.5];
        // The bounds of the step
        let xl = vec![0.0, -1.5];
        let xu = vec![2.0, 0.5];

        let mut x = vec![0.0; n];
        let mut y = vec![0.0; m + 2 * n + 2];
        let mut jw = vec![0i32; m + 2 * n + 2];
        let mode = lsq(
            m, 0, n, false, 0.0, &lf, &gradx, &cmat, &d, &xl, &xu, &mut x, &mut y, &mut jw,
        );

        assert_eq!(mode, 1);
        assert_relative_eq!(x[0], 1.0 / 3.0, epsilon = 1e-6);
        assert_relative_eq!(x[1], -1.0 / 3.0, epsilon = 1e-6);
    }

    /// The solution of the least distance problem (the port of scipy's `ldp`):
    /// `min 0.5 |x|^2` subject to `G x >= h`.
    #[test]
    fn test_ldp_matches_the_least_distance_problem() {
        // G = I and h = (1, 1): the solution is x = (1, 1)
        let g = vec![1.0, 0.0, 0.0, 1.0];
        let h = vec![1.0, 1.0];
        let mut x = vec![0.0; 2];
        let mut indices = vec![0i32; 2];
        let mut multipliers = vec![0.0; 2];
        let mode = ldp(2, 2, &g, &h, &mut x, &mut indices, &mut multipliers);
        assert_eq!(mode, 1);
        assert_relative_eq!(x[0], 1.0, epsilon = 1e-9);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-9);
    }

    /// The solution of the inequality constrained least squares problem (the
    /// port of scipy's `lsi`): `min |A x - b|` subject to `G x >= h`.
    #[test]
    fn test_lsi_matches_the_constrained_least_squares() {
        // A = I and b = 0: the solution is the projection of 0 onto
        // x >= (1, 1), which is x = (1, 1)
        let a = vec![1.0, 0.0, 0.0, 1.0];
        let b = vec![0.0, 0.0];
        let g = vec![1.0, 0.0, 0.0, 1.0];
        let mut h = vec![1.0, 1.0];
        let mut x = vec![0.0; 2];
        let mut jw = vec![0i32; 2];
        let mut multipliers = vec![0.0; 2];

        let mode = lsi(
            2,
            2,
            2,
            &a,
            &b,
            &g,
            &mut h,
            &mut x,
            &mut jw,
            &mut multipliers,
        );
        assert_eq!(mode, 1);
        assert_relative_eq!(x[0], 1.0, epsilon = 1e-9);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-9);
    }

    /// The non-negative least squares solution of the dual problem of the
    /// least distance problem (the reference values are from scipy's
    /// `scipy.optimize.nnls`).
    #[test]
    fn test_nnls_of_the_dual_problem_matches_scipy() {
        // The dual problem: A = [[G'],[h']] and b = (0, 0, 1) (the column
        // major layout: the element (i, j) is at i + j * 3)
        let mut a = vec![
            -1.0,
            0.7071067811865475,
            2.7071067811865475,
            -1.0,
            0.0,
            -1.414213562373095,
            0.0,
            1.0,
            -0.7071067811865475,
            0.0,
            -1.0,
            -1.5,
            0.0,
            -1.0,
            0.20710678118654746,
            0.0,
            0.0,
            -2.0,
            1.0,
            -1.0,
            -2.2071067811865475,
            0.0,
            1.0,
            0.0,
            0.0,
            -1.0,
            0.20710678118654746,
        ];
        let mut b = vec![0.0, 0.0, 1.0];
        let mut x = vec![0.0; 9];
        let mut w = vec![0.0; 9];
        let mut zz = vec![0.0; 3];
        let mut indices = vec![0i32; 9];
        let mut rnorm = 0.0;

        let mode = nnls(
            3,
            9,
            &mut a,
            &mut b,
            &mut x,
            &mut w,
            &mut zz,
            &mut indices,
            27,
            &mut rnorm,
        );
        assert_eq!(mode, 1);
        // scipy: y = [2, 0, 0, 0, 0, 0, 2, 0.58578644, 0], the residual norm 0
        assert_relative_eq!(x[0], 2.0, epsilon = 1e-6);
        assert_relative_eq!(x[6], 2.0, epsilon = 1e-6);
        assert_relative_eq!(x[7], 0.58578644, epsilon = 1e-6);
        assert_relative_eq!(rnorm, 0.0, epsilon = 1e-9);
    }
}
