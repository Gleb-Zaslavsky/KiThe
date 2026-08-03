//! Typed route-specific controls for the fixed-pressure, fixed-enthalpy solve.
//!
//! The public workflow keeps a compatibility options struct because it is a
//! convenient builder for callers. These smaller value objects describe what
//! each numerical route is actually allowed to consume: acceptance belongs to
//! the coupled formulation, scalar/bracket controls belong to the nested
//! reference route, and the monolithic route receives no bracket policy.

use std::time::Duration;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::EnthalpyScale;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;

/// Policy for sampled reversals of the nested `H_eq(P,T)` branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhMonotonicityPolicy {
    /// Stop rather than silently selecting an arbitrary root of a branch with
    /// multiple crossings or a phase-transition-induced reversal.
    RejectObservedNonMonotonicity,
    /// Keep the safeguarded sign-bracket search for compatibility. The
    /// resulting root is bracketed, but no uniqueness is claimed.
    AllowBracketedSignSearch,
}

impl Default for PhMonotonicityPolicy {
    fn default() -> Self {
        Self::RejectObservedNonMonotonicity
    }
}

/// Shared physical acceptance contract for one P,H candidate.
///
/// This is deliberately independent of either outer route. Both monolithic
/// and nested solving must use the same enthalpy acceptance semantics even
/// though only the nested route has scalar bracket controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhAcceptanceOptions {
    scaled_enthalpy_tolerance: f64,
    absolute_enthalpy_tolerance_joules: f64,
}

impl PhAcceptanceOptions {
    /// Creates a validated acceptance contract.
    pub fn new(
        scaled_enthalpy_tolerance: f64,
        absolute_enthalpy_tolerance_joules: f64,
    ) -> Result<Self, ReactionExtentError> {
        if !scaled_enthalpy_tolerance.is_finite() || scaled_enthalpy_tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "scaled_enthalpy_tolerance",
                message: "tolerance must be finite and positive".to_string(),
            });
        }
        if !absolute_enthalpy_tolerance_joules.is_finite()
            || absolute_enthalpy_tolerance_joules <= 0.0
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "absolute_enthalpy_tolerance_joules",
                message: "tolerance must be finite and positive".to_string(),
            });
        }
        Ok(Self {
            scaled_enthalpy_tolerance,
            absolute_enthalpy_tolerance_joules,
        })
    }

    /// Default acceptance contract used by [`PhTemperatureSolveOptions`].
    pub const fn default_values() -> Self {
        Self {
            scaled_enthalpy_tolerance: 1.0e-8,
            absolute_enthalpy_tolerance_joules: 1.0e-6,
        }
    }

    /// Dimensionless enthalpy residual tolerance.
    pub fn scaled_enthalpy_tolerance(self) -> f64 {
        self.scaled_enthalpy_tolerance
    }

    /// Absolute enthalpy tolerance in joules.
    pub fn absolute_enthalpy_tolerance_joules(self) -> f64 {
        self.absolute_enthalpy_tolerance_joules
    }

    /// Absolute joule limit corresponding to a validated enthalpy scale.
    pub fn accepted_error_limit_joules(self, scale: EnthalpyScale) -> f64 {
        self.absolute_enthalpy_tolerance_joules
            .max(self.scaled_enthalpy_tolerance * scale.joules())
    }

    /// Tests the complete absolute-plus-relative acceptance contract.
    pub fn accepts_error(self, error_joules: f64, scale: EnthalpyScale) -> bool {
        error_joules.is_finite()
            && error_joules.abs() <= self.accepted_error_limit_joules(scale)
    }

    /// Crate-internal spelling used by the monolithic runner.
    pub(crate) fn accepts_enthalpy_error(self, error_joules: f64, scale: EnthalpyScale) -> bool {
        self.accepts_error(error_joules, scale)
    }
}

impl Default for PhAcceptanceOptions {
    fn default() -> Self {
        Self::default_values()
    }
}

/// Controls consumed only by the coupled monolithic P,H formulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhMonolithicOptions {
    acceptance: PhAcceptanceOptions,
}

impl PhMonolithicOptions {
    /// Creates monolithic controls from the common physical acceptance gate.
    pub const fn new(acceptance: PhAcceptanceOptions) -> Self {
        Self { acceptance }
    }

    /// Acceptance gate used by the coupled runner.
    pub const fn acceptance(self) -> PhAcceptanceOptions {
        self.acceptance
    }
}

/// Controls consumed only by the nested scalar/bracket reference route.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhNestedOptions {
    acceptance: PhAcceptanceOptions,
    temperature_tolerance: f64,
    max_iterations: usize,
    max_temperature_evaluations: usize,
    monotonicity_policy: PhMonotonicityPolicy,
    max_wall_time: Option<Duration>,
}

impl PhNestedOptions {
    /// Creates validated nested scalar controls.
    pub fn new(
        acceptance: PhAcceptanceOptions,
        temperature_tolerance: f64,
        max_iterations: usize,
        max_temperature_evaluations: usize,
        monotonicity_policy: PhMonotonicityPolicy,
        max_wall_time: Option<Duration>,
    ) -> Result<Self, ReactionExtentError> {
        if !temperature_tolerance.is_finite() || temperature_tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_tolerance",
                message: "tolerance must be finite and positive".to_string(),
            });
        }
        if max_iterations == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_iterations",
                message: "maximum iterations must be greater than zero".to_string(),
            });
        }
        if max_temperature_evaluations < 2 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_temperature_evaluations",
                message: "at least two evaluations are required to test a bracket".to_string(),
            });
        }
        if max_wall_time.is_some_and(|limit| limit.is_zero()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_wall_time",
                message: "wall-time budget must be positive when provided".to_string(),
            });
        }
        Ok(Self {
            acceptance,
            temperature_tolerance,
            max_iterations,
            max_temperature_evaluations,
            monotonicity_policy,
            max_wall_time,
        })
    }

    /// Shared enthalpy acceptance gate used by nested trials.
    pub const fn acceptance(self) -> PhAcceptanceOptions {
        self.acceptance
    }

    /// Temperature convergence tolerance in K.
    pub const fn temperature_tolerance(self) -> f64 {
        self.temperature_tolerance
    }

    /// Maximum number of interior scalar iterations.
    pub const fn max_iterations(self) -> usize {
        self.max_iterations
    }

    /// Global number of scalar temperature evaluations.
    pub const fn max_temperature_evaluations(self) -> usize {
        self.max_temperature_evaluations
    }

    /// Sampled branch monotonicity policy.
    pub const fn monotonicity_policy(self) -> PhMonotonicityPolicy {
        self.monotonicity_policy
    }

    /// Optional wall-time budget for the nested transaction.
    pub const fn max_wall_time(self) -> Option<Duration> {
        self.max_wall_time
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn route_options_keep_acceptance_shared_but_do_not_mix_controls() {
        let acceptance = PhAcceptanceOptions::new(1.0e-7, 2.0e-6).unwrap();
        let monolithic = PhMonolithicOptions::new(acceptance);
        let nested = PhNestedOptions::new(
            acceptance,
            1.0e-5,
            12,
            14,
            PhMonotonicityPolicy::AllowBracketedSignSearch,
            None,
        )
        .unwrap();

        assert_eq!(monolithic.acceptance(), nested.acceptance());
        assert_eq!(nested.temperature_tolerance(), 1.0e-5);
        assert_eq!(nested.max_iterations(), 12);
        assert_eq!(nested.max_temperature_evaluations(), 14);
        assert_eq!(
            nested.monotonicity_policy(),
            PhMonotonicityPolicy::AllowBracketedSignSearch
        );
    }

    #[test]
    fn invalid_route_specific_controls_fail_at_construction() {
        let acceptance = PhAcceptanceOptions::default_values();
        assert!(PhNestedOptions::new(acceptance, 0.0, 1, 2, Default::default(), None).is_err());
        assert!(PhNestedOptions::new(acceptance, 1.0, 0, 2, Default::default(), None).is_err());
        assert!(PhNestedOptions::new(acceptance, 1.0, 1, 1, Default::default(), None).is_err());
    }
}
