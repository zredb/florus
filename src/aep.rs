/// AEP (Annual Energy Production) calculation module
///
/// Provides functions to calculate farm energy production based on wind data
use crate::types::{Array1, Float};
use crate::wind_data::{ WindData};
use crate::FlorisModel;
use crate::core::Farm;
use crate::core::FlowField;

/// Calculate annual energy production from time series wind data
pub fn calculate_aep_from_time_series(
    farm: &Farm,
    time_series: &dyn WindData,
    _frequency_cutoff: Option<Float>,
) -> AEPResult {
    let wind_speeds = time_series.wind_speeds();
    let wind_directions = time_series.wind_directions();
    let turbulence_intensities = time_series.turbulence_intensities();

    let mut total_energy = 0.0; // Wh
    let mut energy_by_turbine = vec![0.0; farm.n_turbines()];

    let n_conditions = time_series.n_conditions();

    for fi in 0..n_conditions {
        // For time series, each condition represents 1 hour
        // (TimeSeries stores individual time steps, not aggregated frequencies)
        let hours = 1.0;

        // Extract conditions
        let ws = wind_speeds[fi];
        let wd = wind_directions[fi];
        let ti = turbulence_intensities[fi];

        // Create flow field for this condition
        let flow_field = FlowField::new(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![wd]),
            0.0,
            0.12,
            1.225,
            Array1::from_vec(vec![ti]),
            90.0,
        ).unwrap();

        // Create and run model
        let mut model = FlorisModel {
            farm: farm.clone(),
            flow_field,
            state: crate::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            turbine_grid_points: 3,
            model_manager: None,
        };

        // Initialize and run
        if let Err(e) = model.initialize_grid() {
            eprintln!("Grid initialization failed: {:?}", e);
            continue;
        }
        if let Err(e) = model.initialize_flow_field() {
            eprintln!("Flow field initialization failed: {:?}", e);
            continue;
        }
        if let Err(e) = model.run() {
            eprintln!("Model run failed: {:?}", e);
            continue;
        }

        // Get power for each turbine
        let powers = model.get_turbine_powers();
        let shape = powers.shape();
        let n_turbines = shape[1];

        // Energy = Power * Time (hours)
        // Note: powers has shape [1, n_turbines] since each model run is for 1 findex
        let energy_kwh = powers.mapv(|p| p / 1000.0 * hours);

        // Sum energy for this condition
        for ti_idx in 0..n_turbines {
            energy_by_turbine[ti_idx] += energy_kwh[[0, ti_idx]];
            total_energy += energy_kwh[[0, ti_idx]];
        }
    }
    
    let energy_by_turbine_kwh: Vec<Float> = energy_by_turbine.iter().map(|e| e / 1000.0).collect();
    let energy_by_turbine_mwh: Vec<Float> = energy_by_turbine_kwh.iter().map(|e| e / 1000.0).collect();

    AEPResult {
        total_energy_wh: total_energy,
        energy_by_turbine_wh: energy_by_turbine.clone(),
        total_energy_kwh: total_energy / 1000.0,
        energy_by_turbine_kwh,
        total_energy_mwh: total_energy / 1000.0,
        energy_by_turbine_mwh,
        conditions_processed: n_conditions,
    }
}

/// Result of AEP calculation
#[derive(Debug, Clone)]
pub struct AEPResult {
    pub total_energy_wh: Float,
    pub energy_by_turbine_wh: Vec<Float>,
    pub total_energy_kwh: Float,
    pub energy_by_turbine_kwh: Vec<Float>,
    pub total_energy_mwh: Float,
    pub energy_by_turbine_mwh: Vec<Float>,
    pub conditions_processed: usize,
}

impl AEPResult {
    /// Get total annual energy in MWh
    pub fn total_mwh(&self) -> Float {
        self.total_energy_mwh
    }
    
    /// Get energy by turbine in MWh
    pub fn by_turbine_mwh(&self) -> &[Float] {
        &self.energy_by_turbine_mwh
    }
    
    /// Get capacity factor
    pub fn capacity_factor(&self, rated_power_watts: Float, n_turbines: usize) -> Float {
        let annual_production_wh = self.total_energy_wh;
        let max_production_wh = rated_power_watts * n_turbines as Float * 8760.0;
        if max_production_wh > 0.0 {
            annual_production_wh / max_production_wh * 100.0
        } else {
            0.0
        }
    }
}

/// Calculate power at specific wind conditions
pub fn calculate_power_at_conditions(
    farm: &Farm,
    wind_speed: Float,
    wind_direction: Float,
    turbulence_intensity: Float,
) -> Vec<Float> {
    let flow_field = FlowField::new(
        Array1::from_vec(vec![wind_speed]),
        Array1::from_vec(vec![wind_direction]),
        0.0,
        0.12,
        1.225,
        Array1::from_vec(vec![turbulence_intensity]),
        90.0,
    ).unwrap();
    
    let mut model = FlorisModel {
        farm: farm.clone(),
        flow_field,
        state: crate::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        turbine_grid_points: 3,
        model_manager: None,
    };
    
    let _ = model.initialize_grid();
    let _ = model.initialize_flow_field();
    let _ = model.run();
    
    let powers = model.get_turbine_powers();
    let shape = powers.shape();
    let n_turbines = shape[1];
    
    let mut result = Vec::with_capacity(n_turbines);
    for ti in 0..n_turbines {
        result.push(powers[[0, ti]]);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;
    use crate::wind_data::TimeSeries;
    
    #[test]
    fn test_aep_result_format() {
        let result = AEPResult {
            total_energy_wh: 1000000.0,
            energy_by_turbine_wh: vec![500000.0, 500000.0],
            total_energy_kwh: 1000.0,
            energy_by_turbine_kwh: vec![500.0, 500.0],
            total_energy_mwh: 1.0,
            energy_by_turbine_mwh: vec![0.5, 0.5],
            conditions_processed: 10,
        };
        
        assert_eq!(result.total_mwh(), 1.0);
        assert_eq!(result.by_turbine_mwh().len(), 2);
    }
    
    #[test]
    fn test_capacity_factor() {
        // 50% CF for 5MW turbine: 5MW * 8760h * 0.5 = 21,900 MWh = 21,900,000,000 Wh
        let result = AEPResult {
            total_energy_wh: 21900000000.0,
            energy_by_turbine_wh: vec![21900000000.0],
            total_energy_kwh: 21900000.0,
            energy_by_turbine_kwh: vec![21900000.0],
            total_energy_mwh: 21900.0,
            energy_by_turbine_mwh: vec![21900.0],
            conditions_processed: 10,
        };
        
        let cf = result.capacity_factor(5_000_000.0, 1);
        assert!((cf - 50.0).abs() < 0.1);
    }
    
    #[test]
    fn test_time_series_aep() {
        let layout_x = Array1::from_vec(vec![0.0, 630.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];

        let farm = crate::core::Farm::new(layout_x, layout_y, turbine_types).unwrap();

        let ws = Array1::from_vec(vec![8.0, 8.0, 10.0, 10.0]);
        let wd = Array1::from_vec(vec![270.0, 270.0, 270.0, 270.0]);
        let ti = Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06]);
      //  let freq = Array1::from_vec(vec![2190.0, 2190.0, 2190.0, 2190.0]); // 4 seasons, 2190 hours each

        let time_series = TimeSeries::new(wd, ws, ti).unwrap();

        // Create a simple AEP calculation test
        let result = calculate_aep_from_time_series(&farm, &time_series, None);

        assert!(result.conditions_processed > 0, "No conditions processed");
        assert!(result.total_energy_mwh > 0.0, "Total energy is 0: {:?}", result);
    }
}
