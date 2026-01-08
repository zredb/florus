/// State management for solver

use serde::{Deserialize, Serialize};

/// Solver state for tracking computation progress
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct State {
    /// Whether the flow field has been initialized
    pub initialized: bool,
    
    /// Current iteration count
    pub iteration: usize,
    
    /// Maximum iterations allowed
    pub max_iterations: usize,
    
    /// Convergence tolerance
    pub tolerance: f64,
    
    /// Whether the solution has converged
    pub converged: bool,
}

impl State {
    /// Create a new state
    pub fn new() -> Self {
        Self {
            initialized: false,
            iteration: 0,
            max_iterations: 100,
            tolerance: 1e-6,
            converged: false,
        }
    }
    
    /// Reset state
    pub fn reset(&mut self) {
        self.initialized = false;
        self.iteration = 0;
        self.converged = false;
    }
}

impl Default for State {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_state_new() {
        let state = State::new();
        
        assert!(!state.initialized);
        assert_eq!(state.iteration, 0);
        assert_eq!(state.max_iterations, 100);
        assert_relative_eq!(state.tolerance, 1e-6);
        assert!(!state.converged);
    }

    #[test]
    fn test_state_reset() {
        let mut state = State::new();
        state.initialized = true;
        state.iteration = 50;
        state.converged = true;
        
        state.reset();
        
        assert!(!state.initialized);
        assert_eq!(state.iteration, 0);
        assert!(!state.converged);
    }

    #[test]
    fn test_state_default() {
        let state = State::default();
        
        assert!(!state.initialized);
        assert_eq!(state.iteration, 0);
        assert!(!state.converged);
    }

    #[test]
    fn test_state_clone() {
        let state = State::new();
        let cloned = state.clone();
        
        assert_eq!(state.initialized, cloned.initialized);
        assert_eq!(state.iteration, cloned.iteration);
        assert_eq!(state.tolerance, cloned.tolerance);
        assert_eq!(state.converged, cloned.converged);
    }
}
