//! Active Wake Control turbine model
//!
//! Placeholder implementation for AWC/helix wake mixing control
//! Current implementation delegates to CosineLossTurbine
//!
//! Full AWC model requires:
//! - Helix tuning parameters from turbine YAML
//! - Helix power/thrust coefficient tables
//! - Implementation of helix wake mixing model

use crate::types::Float;
use crate::types::Array2;
use crate::core::turbine::operation_models::base::*;
use super::CosineLossTurbine;

#[derive(Debug, Clone, Default)]
pub struct AWCTurbine;

impl OperationModel for AWCTurbine {
    fn model_name(&self) -> &'static str {
        "awc"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        // AWC model requires specific parameters from turbine type
        // Currently delegated to CosineLossTurbine as placeholder
        let cosine = CosineLossTurbine;
        cosine.power(params, ctx)
    }

    fn thrust_coefficient(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let cosine = CosineLossTurbine;
        cosine.thrust_coefficient(params, ctx)
    }

    fn axial_induction(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let cosine = CosineLossTurbine;
        cosine.axial_induction(params, ctx)
    }
}
