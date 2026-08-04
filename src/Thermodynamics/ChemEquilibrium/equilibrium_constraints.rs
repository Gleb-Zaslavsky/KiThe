//! Typed thermodynamic constraints and small, pure energy helpers.
//!
//! This module deliberately does not solve an equilibrium problem. It defines
//! the boundary that a future `P,H` workflow will consume while the existing
//! fixed-`P,T` workflow remains the numerical source of truth.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;

/// A validated, finite temperature interval shared by thermochemical records.
///
/// The interval is inclusive. It represents the intersection of the selected
/// records' domains, not an extrapolation permission for NASA/NIST formulas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TemperatureBounds {
    lower: f64,
    upper: f64,
}

/// Extensive total enthalpy expressed in joules.
///
/// This newtype is the unit boundary for the `P,H` workflow. GUI or service
/// layers may convert molar, mass-specific, or normalized inputs before they
/// construct this value, but the equilibrium engine never receives those
/// ambiguous representations directly.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct TotalEnthalpyJoules(f64);

impl TotalEnthalpyJoules {
    /// Creates a finite extensive enthalpy value in joules.
    pub fn new(value: f64) -> Result<Self, ReactionExtentError> {
        if !value.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "target_enthalpy",
                message: "total enthalpy in joules must be finite".to_string(),
            });
        }
        Ok(Self(value))
    }

    /// Returns the numeric value in joules.
    pub fn joules(self) -> f64 {
        self.0
    }
}

impl TemperatureBounds {
    /// Creates bounds in kelvin.
    pub fn new(lower: f64, upper: f64) -> Result<Self, ReactionExtentError> {
        if !lower.is_finite() || !upper.is_finite() || lower <= 0.0 || upper <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_bounds",
                message: "temperature bounds must be finite and strictly positive".to_string(),
            });
        }
        if lower > upper {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_bounds",
                message: "temperature lower bound must not exceed upper bound".to_string(),
            });
        }
        Ok(Self { lower, upper })
    }

    /// Lower bound in kelvin.
    pub fn lower(self) -> f64 {
        self.lower
    }

    /// Upper bound in kelvin.
    pub fn upper(self) -> f64 {
        self.upper
    }

    /// Returns whether a temperature belongs to this interval.
    pub fn contains(self, temperature: f64) -> bool {
        temperature.is_finite() && (self.lower..=self.upper).contains(&temperature)
    }

    /// Intersects two record domains without permitting extrapolation.
    pub fn intersect(self, other: Self) -> Result<Self, ReactionExtentError> {
        Self::new(self.lower.max(other.lower), self.upper.min(other.upper))
    }
}

/// The thermodynamic constraint imposed on an equilibrium calculation.
///
/// `PT` is the existing fixed-temperature problem. `PH` stores total
/// extensive enthalpy in joules and a temperature seed for the future outer
/// temperature solve. The seed is not the final equilibrium temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EquilibriumConstraint {
    /// Fixed temperature and pressure.
    PT { conditions: EquilibriumConditions },
    /// Fixed pressure and total enthalpy.
    PH {
        /// System pressure in pascals.
        pressure: f64,
        /// Standard-state reference pressure in pascals.
        reference_pressure: f64,
        /// Target total enthalpy in joules.
        target_enthalpy: TotalEnthalpyJoules,
        /// Initial temperature estimate in kelvin; not the solved result.
        initial_temperature: f64,
    },
}

impl EquilibriumConstraint {
    /// Constructs a fixed-`P,T` constraint from already validated conditions.
    pub fn pt(conditions: EquilibriumConditions) -> Self {
        Self::PT { conditions }
    }

    /// Constructs and validates a fixed-`P,H` constraint.
    pub fn ph(
        pressure: f64,
        reference_pressure: f64,
        target_enthalpy: f64,
        initial_temperature: f64,
    ) -> Result<Self, ReactionExtentError> {
        let target_enthalpy = TotalEnthalpyJoules::new(target_enthalpy)?;
        Self::ph_joules(
            pressure,
            reference_pressure,
            target_enthalpy,
            initial_temperature,
        )
    }

    /// Constructs a fixed-`P,H` constraint from an explicitly typed joule
    /// target. This is the canonical constructor for new production callers.
    pub fn ph_joules(
        pressure: f64,
        reference_pressure: f64,
        target_enthalpy: TotalEnthalpyJoules,
        initial_temperature: f64,
    ) -> Result<Self, ReactionExtentError> {
        for (field, value) in [
            ("pressure", pressure),
            ("reference_pressure", reference_pressure),
            ("initial_temperature", initial_temperature),
        ] {
            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field,
                    message: "value must be finite".to_string(),
                });
            }
        }
        if pressure <= 0.0 || reference_pressure <= 0.0 || initial_temperature <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_constraint",
                message: "pressure, reference pressure, and initial temperature must be positive"
                    .to_string(),
            });
        }

        Ok(Self::PH {
            pressure,
            reference_pressure,
            target_enthalpy,
            initial_temperature,
        })
    }

    /// Returns the pressure in pascals for either constraint.
    pub fn pressure(self) -> f64 {
        match self {
            Self::PT { conditions } => conditions.pressure(),
            Self::PH { pressure, .. } => pressure,
        }
    }

    /// Returns the standard-state reference pressure in pascals.
    pub fn reference_pressure(self) -> f64 {
        match self {
            Self::PT { conditions } => conditions.reference_pressure(),
            Self::PH {
                reference_pressure, ..
            } => reference_pressure,
        }
    }

    /// Returns the fixed temperature, if this is a `PT` constraint.
    pub fn fixed_temperature(self) -> Option<f64> {
        match self {
            Self::PT { conditions } => Some(conditions.temperature()),
            Self::PH { .. } => None,
        }
    }

    /// Returns the target extensive enthalpy, if this is a `PH` constraint.
    pub fn target_enthalpy(self) -> Option<f64> {
        match self {
            Self::PT { .. } => None,
            Self::PH {
                target_enthalpy, ..
            } => Some(target_enthalpy.joules()),
        }
    }

    /// Returns the typed total enthalpy target in joules.
    pub fn target_enthalpy_joules(self) -> Option<TotalEnthalpyJoules> {
        match self {
            Self::PT { .. } => None,
            Self::PH {
                target_enthalpy, ..
            } => Some(target_enthalpy),
        }
    }

    /// Returns the user-provided `PH` temperature seed, if present.
    pub fn initial_temperature(self) -> Option<f64> {
        match self {
            Self::PT { .. } => None,
            Self::PH {
                initial_temperature,
                ..
            } => Some(initial_temperature),
        }
    }

    /// Materializes conditions at a trial temperature.
    ///
    /// This is the bridge used by the future outer scalar solver. Calling it
    /// for `PT` is also useful when one generic workflow handles both modes.
    pub fn conditions_at(
        self,
        temperature: f64,
    ) -> Result<EquilibriumConditions, ReactionExtentError> {
        let temperature = match self {
            Self::PT { conditions } => conditions.temperature(),
            Self::PH { .. } => temperature,
        };
        EquilibriumConditions::new(temperature, self.pressure(), self.reference_pressure())
    }

    /// Validates that a fixed or trial temperature lies in a record domain.
    pub fn validate_temperature(
        self,
        temperature: f64,
        bounds: TemperatureBounds,
    ) -> Result<(), ReactionExtentError> {
        if !bounds.contains(temperature) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature",
                message: format!(
                    "temperature {temperature} K is outside [{}, {}] K",
                    bounds.lower(),
                    bounds.upper()
                ),
            });
        }
        Ok(())
    }
}

/// Computes additive total enthalpy in joules from moles and J/mol values.
///
/// This helper intentionally models only the currently supported ideal/pure
/// additive contract. Unsupported excess or mixing enthalpy belongs in a
/// separate phase-model capability rather than being silently omitted here.
pub fn additive_total_enthalpy(
    moles: &[f64],
    molar_enthalpies: &[f64],
) -> Result<f64, ReactionExtentError> {
    if moles.len() != molar_enthalpies.len() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "enthalpy_dimensions",
            message: format!(
                "moles has length {}, molar enthalpies has length {}",
                moles.len(),
                molar_enthalpies.len()
            ),
        });
    }

    let mut total = 0.0;
    for (index, (&moles_i, &enthalpy_i)) in moles.iter().zip(molar_enthalpies).enumerate() {
        if !moles_i.is_finite() || moles_i < 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "moles",
                message: format!("moles[{index}] must be finite and non-negative"),
            });
        }
        if !enthalpy_i.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "molar_enthalpy",
                message: format!("molar_enthalpies[{index}] must be finite"),
            });
        }
        total += moles_i * enthalpy_i;
    }

    if !total.is_finite() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "total_enthalpy",
            message: "total enthalpy is not finite".to_string(),
        });
    }
    Ok(total)
}

/// Validated positive scale for the dimensionless enthalpy residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnthalpyScale(f64);

impl EnthalpyScale {
    /// Smallest allowed scale in joules. It prevents division by zero while
    /// preserving the actual root of `H - H_target = 0`.
    pub const MINIMUM_JOULES: f64 = 1.0e-12;

    /// Builds a scale from target and initial-state magnitudes in joules.
    pub fn from_magnitudes(
        target_enthalpy: f64,
        initial_moles: &[f64],
        initial_molar_enthalpies: &[f64],
    ) -> Result<Self, ReactionExtentError> {
        let initial_enthalpy_magnitude = initial_moles
            .iter()
            .zip(initial_molar_enthalpies)
            .map(|(&moles_i, &enthalpy_i)| {
                if !moles_i.is_finite() || moles_i < 0.0 || !enthalpy_i.is_finite() {
                    None
                } else {
                    Some(moles_i * enthalpy_i.abs())
                }
            })
            .try_fold(0.0, |sum, term| term.map(|term| sum + term))
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "enthalpy_scale",
                message: "initial enthalpy scale inputs must be finite and non-negative"
                    .to_string(),
            })?;

        if !target_enthalpy.is_finite() || !initial_enthalpy_magnitude.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_scale",
                message: "enthalpy scale inputs must be finite".to_string(),
            });
        }

        Ok(Self(
            target_enthalpy
                .abs()
                .max(initial_enthalpy_magnitude)
                .max(Self::MINIMUM_JOULES),
        ))
    }

    /// Creates a scale from one explicit positive finite value in joules.
    pub fn new(value: f64) -> Result<Self, ReactionExtentError> {
        if !value.is_finite() || value <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_scale",
                message: "enthalpy scale must be finite and positive".to_string(),
            });
        }
        Ok(Self(value.max(Self::MINIMUM_JOULES)))
    }

    /// Scale in joules.
    pub fn joules(self) -> f64 {
        self.0
    }

    /// Converts a raw enthalpy error in joules into a dimensionless residual.
    pub fn scale_error(self, error_joules: f64) -> Result<f64, ReactionExtentError> {
        if !error_joules.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_error",
                message: "enthalpy error must be finite".to_string(),
            });
        }
        Ok(error_joules / self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ph_constraint_keeps_seed_distinct_from_solved_temperature() {
        let constraint = EquilibriumConstraint::ph(101_325.0, 101_325.0, -2_000.0, 900.0).unwrap();

        assert_eq!(constraint.fixed_temperature(), None);
        assert_eq!(constraint.initial_temperature(), Some(900.0));
        assert_eq!(constraint.target_enthalpy(), Some(-2_000.0));
        assert_eq!(
            constraint.target_enthalpy_joules().unwrap().joules(),
            -2_000.0
        );
        assert_eq!(
            constraint.conditions_at(1_100.0).unwrap().temperature(),
            1_100.0
        );
    }

    #[test]
    fn temperature_bounds_intersect_without_extrapolation() {
        let first = TemperatureBounds::new(200.0, 2_000.0).unwrap();
        let second = TemperatureBounds::new(300.0, 1_500.0).unwrap();
        let intersection = first.intersect(second).unwrap();

        assert_eq!(
            intersection,
            TemperatureBounds::new(300.0, 1_500.0).unwrap()
        );
        assert!(intersection.contains(300.0));
        assert!(!intersection.contains(1_500.1));
        assert!(first
            .intersect(TemperatureBounds::new(2_100.0, 3_000.0).unwrap())
            .is_err());
    }

    #[test]
    fn additive_enthalpy_and_scale_have_explicit_units_contract() {
        let total = additive_total_enthalpy(&[2.0, 0.5], &[10.0, -20.0]).unwrap();
        assert_eq!(total, 10.0);

        let scale = EnthalpyScale::from_magnitudes(100.0, &[2.0, 0.5], &[10.0, -20.0]).unwrap();
        assert_eq!(scale.joules(), 100.0);
        assert_eq!(scale.scale_error(-25.0).unwrap(), -0.25);
    }

    #[test]
    fn enthalpy_helpers_reject_dimension_and_domain_errors() {
        assert!(additive_total_enthalpy(&[1.0], &[]).is_err());
        assert!(additive_total_enthalpy(&[-1.0], &[10.0]).is_err());
        assert!(additive_total_enthalpy(&[1.0], &[f64::NAN]).is_err());
        assert!(EquilibriumConstraint::ph(101_325.0, 101_325.0, 0.0, 0.0).is_err());
        assert!(EnthalpyScale::new(0.0).is_err());
    }

    #[test]
    fn typed_joule_constructor_is_the_canonical_ph_boundary() {
        let target = TotalEnthalpyJoules::new(-2_000.0).unwrap();
        let constraint =
            EquilibriumConstraint::ph_joules(101_325.0, 101_325.0, target, 900.0).unwrap();

        assert_eq!(constraint.target_enthalpy_joules(), Some(target));
        assert!(TotalEnthalpyJoules::new(f64::INFINITY).is_err());
    }
}
