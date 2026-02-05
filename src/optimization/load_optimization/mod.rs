//! Load optimization module for FLORIS-RS
//!
//! This module implements load-aware turbine optimization that balances power production
//! with variable operating costs (VOC). It extends FLORIS beyond power maximization to:
//!
//! 1. Compute Load Turbulence Intensity (LTI) - IEC 61400-1 Ed.4 standard-based turbulence calculation
//! 2. Calculate Variable Operating Costs (VOC) - Cost model based on wind speed variation and turbine thrust
//! 3. Optimize Power Setpoints - Sequential derating optimization to maximize net revenue
//!
//! This module corresponds to floris/optimization/load_optimization/ in Python FLORIS v4.6.

pub mod load_optimization;

// Re-export all public functions and constants
pub use load_optimization::{
    compute_lti,
    compute_turbine_voc,
    compute_farm_voc,
    compute_farm_revenue,
    compute_net_revenue,
    find_a_to_satisfy_rev_voc_ratio,
    find_a_to_satisfy_target_voc_per_mw,
    optimize_power_setpoints,
    POWER_SETPOINT_DEFAULT,
    POWER_SETPOINT_DISABLED,
};
