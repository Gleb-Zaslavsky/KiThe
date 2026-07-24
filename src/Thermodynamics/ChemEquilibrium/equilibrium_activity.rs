//! Activity models for chemical species in the canonical equilibrium solver.
//!
//! # Purpose
//!
//! This module defines the **activity model** for each thermodynamic phase in
//! the equilibrium calculation. The activity `a_i` relates the chemical potential
//! to the standard state:
//!
//! ```text
//! μ_i(T, P, x) = μ_i°(T) + R·T·ln(a_i)
//! ```
//!
//! The module provides two models: ideal gas and ideal solution. Both are used
//! by the residual and Jacobian generators in [`equilibrium_log_moles`](super::equilibrium_log_moles)
//! and by the symbolic engine in [`equilibrium_rst_backend`](super::equilibrium_rst_backend).
//!
//! # Physical and Mathematical Background
//!
//! ## Ideal Gas
//!
//! For an ideal gas mixture:
//!
//! ```text
//! a_i = x_i · P / P0 = (n_i / N_phase) · (P / P0)
//! ln(a_i) = ln(n_i) - ln(N_phase) + ln(P / P0)
//! ```
//!
//! where `x_i` is the mole fraction, `P` is the system pressure, and `P0` is
//! the reference pressure (typically 1 atm = 101325 Pa).
//!
//! ## Ideal Solution
//!
//! For an ideal condensed solution:
//!
//! ```text
//! a_i = x_i = n_i / N_phase
//! ln(a_i) = ln(n_i) - ln(N_phase)
//! ```
//!
//! For a pure condensed phase, `x_i = 1` and `ln(a_i) = 0`.
//!
//! ## Phase Offset
//!
//! The `log_phase_offset` method returns the pressure-dependent part:
//!
//! - Ideal gas: `ln(P / P0)`
//! - Ideal solution: `0.0`
//!
//! This offset is **constant with respect to log-mole coordinates**, so the
//! analytical Jacobian differentiates only the `ln(n_i) - ln(N_phase)` terms.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`PhaseActivityModel`] | Enum: IdealGas or IdealSolution |
//!
//! # Key Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`PhaseActivityModel::log_activity`] | Full log-activity from physical moles |
//! | [`PhaseActivityModel::log_phase_offset`] | Pressure-dependent offset only |
//! | [`phase_activity_models`] | Maps phase descriptors to activity models |
//!
//! # Dataflow
//!
//! ```text
//!   Phase descriptor (from EquilibriumLogMoles or EquilibriumProblem)
//!     │
//!     v
//!   PhaseActivityModel::log_phase_offset(P, P0)
//!     ├── IdealGas → ln(P / P0)
//!     └── IdealSolution → 0.0
//!     │
//!     v
//!   PhaseActivityModel::log_activity(n_i, N_phase, P, P0)
//!     ├── Validates inputs (finite, positive)
//!     ├── Computes ln(n_i / N_phase) + phase_offset
//!     └── Returns ln(a_i)
//!     │
//!     v
//!   Used by: evaluate_equilibrium_logmole_residual()
//!   Used by: equilibrium_logmole_residual2() (deprecated)
//!   Used by: RST symbolic engine
//! ```
//!
//! # Examples
//!
//! ```rust
//! use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
//!
//! let model = PhaseActivityModel::IdealGas;
//! let offset = model.log_phase_offset(101325.0, 101325.0).unwrap();
//! assert_eq!(offset, 0.0); // P == P0
//!
//! let activity = model.log_activity(0.5, 1.0, 101325.0, 101325.0).unwrap();
//! assert!((activity - (-0.693)).abs() < 1e-3); // ln(0.5)
//! ```
//!
//! # Non-obvious Details
//!
//! - Both `log_activity` and `log_phase_offset` validate their inputs and return
//!   [`ReactionExtentError::InvalidConditions`] for non-finite or non-positive values.
//! - The activity model is **per-phase**, not per-species. All species in the same
//!   phase share the same activity model.
//! - For pure condensed phases (IdealSolution with one species), `log_activity`
//!   returns `0.0` because `n_i = N_phase` and the phase offset is zero.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — residual/Jacobian evaluation
//! - [`equilibrium_rst_backend`](super::equilibrium_rst_backend) — symbolic backend
//! - [`equilibrium_problem`](super::equilibrium_problem) — problem definition
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Phase;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;

/// Supported activity correction for one declared phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseActivityModel {
    /// `a_i = x_i * P / P0`.
    IdealGas,
    /// `a_i = x_i`; for a pure condensed phase this is exactly one.
    IdealSolution,
}

impl PhaseActivityModel {
    /// Constant contribution to `ln(a_i)` at the requested pressure.
    pub fn log_phase_offset(
        self,
        pressure: f64,
        reference_pressure: f64,
    ) -> Result<f64, ReactionExtentError> {
        for (parameter, value) in [
            ("pressure", pressure),
            ("reference_pressure", reference_pressure),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidConditions { parameter, value });
            }
        }
        Ok(match self {
            Self::IdealGas => (pressure / reference_pressure).ln(),
            Self::IdealSolution => 0.0,
        })
    }

    /// Evaluates the full log activity of one component from physical moles.
    pub fn log_activity(
        self,
        species_moles: f64,
        phase_moles: f64,
        pressure: f64,
        reference_pressure: f64,
    ) -> Result<f64, ReactionExtentError> {
        if !species_moles.is_finite() || species_moles <= 0.0 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "activity_species_moles",
                message: format!(
                    "species mole number must be finite and positive, got {species_moles}"
                ),
            });
        }
        if !phase_moles.is_finite() || phase_moles <= 0.0 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "activity_phase_moles",
                message: format!(
                    "phase mole number must be finite and positive, got {phase_moles}"
                ),
            });
        }
        Ok((species_moles / phase_moles).ln()
            + self.log_phase_offset(pressure, reference_pressure)?)
    }
}

/// Builds models in declared phase order.
pub fn phase_activity_models(phases: &[Phase]) -> Vec<PhaseActivityModel> {
    phases.iter().map(|phase| phase.kind).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gas_and_solution_activity_offsets_follow_their_contracts() {
        assert!(
            (PhaseActivityModel::IdealGas
                .log_phase_offset(202_650.0, 101_325.0)
                .unwrap()
                - 2.0_f64.ln())
            .abs()
                < 1e-12
        );
        assert_eq!(
            PhaseActivityModel::IdealSolution
                .log_phase_offset(202_650.0, 101_325.0)
                .unwrap(),
            0.0
        );
    }

    #[test]
    fn pure_condensed_phase_has_unit_activity() {
        assert_eq!(
            PhaseActivityModel::IdealSolution
                .log_activity(3.0, 3.0, 101_325.0, 101_325.0)
                .unwrap(),
            0.0
        );
    }
}
