//! Explicit boundary between physical and internally normalized extensive data.
//!
//! | Quantity | Internal normalized representation | Public physical meaning |
//! | --- | --- | --- |
//! | Component/phase/element moles | `n / s` | `n` mol |
//! | Total enthalpy and absolute energy error | `H / s` | `H` J |
//! | Absolute mole thresholds (`phase_eps`, absolute trace floor) | `eps / s` | `eps` mol |
//! | Relative trace fraction | unchanged | dimensionless fraction |
//! | Temperature, pressure, `p0`, mole fractions, activities | unchanged | intensive |
//! | `G0`, `H0`, `S0`, reaction affinity, chemical potential, TPD | unchanged | intensive |
//! | TPD create/keep hysteresis and phase topology | unchanged | physical decision semantics |
//!
//! This module deliberately performs no solve and chooses no retry policy.
//! It only makes exact extensive transformations explicit so a future solver
//! route cannot silently mix normalized internal units with physical reports.
//!
//! Continuation remains a physical-state concern. A normalized fresh solve is
//! still fresh; this type never promotes a normalized coordinate vector into a
//! physical continuation seed. Any future normalized continuation route must
//! reconstruct physical component amounts and retain accepted phase history
//! before it crosses the public workflow boundary.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TotalEnthalpyJoules;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::LogMolesInitialGuess;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::TraceSpeciesSeedPolicy;

/// Answer-independent mapping between physical and normalized extensive space.
///
/// The scale is the sum of positive physical input component moles. It must
/// never be inferred from an accepted solution, expected topology, external
/// reference data, or continuation state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtensiveNormalization {
    physical_inventory_scale: f64,
}

impl ExtensiveNormalization {
    /// Builds the mapping from one physical input mole vector.
    pub fn from_physical_moles(moles: &[f64]) -> Result<Self, ReactionExtentError> {
        if moles
            .iter()
            .any(|amount| !amount.is_finite() || *amount < 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "extensive_normalization_moles",
                message: "physical input moles must be finite and non-negative".to_owned(),
            });
        }
        let scale = moles.iter().copied().filter(|amount| *amount > 0.0).sum();
        Self::from_physical_inventory_scale(scale)
    }

    /// Builds the mapping from an already audited physical inventory scale.
    pub fn from_physical_inventory_scale(scale: f64) -> Result<Self, ReactionExtentError> {
        if !scale.is_finite() || scale <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "extensive_normalization_scale",
                message: "physical inventory scale must be finite and strictly positive".to_owned(),
            });
        }
        Ok(Self {
            physical_inventory_scale: scale,
        })
    }

    /// Physical total inventory represented by one normalized mole.
    pub fn physical_inventory_scale(self) -> f64 {
        self.physical_inventory_scale
    }

    /// Converts physical component moles to internal normalized moles.
    pub fn normalize_moles(self, physical_moles: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        self.transform_moles(physical_moles, 1.0 / self.physical_inventory_scale)
    }

    /// Reconstructs physical component moles from internal normalized moles.
    pub fn denormalize_moles(
        self,
        normalized_moles: &[f64],
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.transform_moles(normalized_moles, self.physical_inventory_scale)
    }

    /// Converts one physical log-mole seed to normalized solver coordinates.
    ///
    /// Uniform extensive scaling is translation in log space:
    /// `ln(n / s) = ln(n) - ln(s)`. No exponentiation is needed, so very small
    /// trace coordinates retain their finite representation.
    pub fn normalize_log_mole_guess(
        self,
        physical: &LogMolesInitialGuess,
    ) -> Result<LogMolesInitialGuess, ReactionExtentError> {
        self.transform_log_moles(physical.as_slice(), -self.physical_inventory_scale.ln())
    }

    /// Reconstructs one physical log-mole seed from normalized coordinates.
    pub fn denormalize_log_moles(
        self,
        normalized: &[f64],
    ) -> Result<LogMolesInitialGuess, ReactionExtentError> {
        self.transform_log_moles(normalized, self.physical_inventory_scale.ln())
    }

    /// Converts phase totals expressed in physical moles to normalized moles.
    ///
    /// Phase totals share the same non-negative extensive-unit contract as
    /// component moles; this named operation makes report conversion explicit.
    pub fn normalize_phase_totals(
        self,
        physical_phase_totals: &[f64],
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.normalize_moles(physical_phase_totals)
    }

    /// Reconstructs physical phase totals from normalized moles.
    pub fn denormalize_phase_totals(
        self,
        normalized_phase_totals: &[f64],
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.denormalize_moles(normalized_phase_totals)
    }

    /// Converts non-negative elemental inventories to normalized moles.
    ///
    /// The element basis is not changed: every element total is divided by the
    /// same inventory scale as the component vector.
    pub fn normalize_element_totals(
        self,
        physical_element_totals: &[f64],
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.normalize_moles(physical_element_totals)
    }

    /// Reconstructs physical elemental inventories from normalized moles.
    pub fn denormalize_element_totals(
        self,
        normalized_element_totals: &[f64],
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.denormalize_moles(normalized_element_totals)
    }

    /// Converts a typed physical composition while retaining its layout proof.
    pub fn normalize_initial_composition(
        self,
        layout: &MultiphaseEquilibriumLayout,
        physical: &MultiphaseInitialComposition,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        physical.validate_for(layout)?;
        MultiphaseInitialComposition::from_dense(layout, self.normalize_moles(physical.moles())?)
    }

    /// Converts physical total enthalpy to normalized total enthalpy.
    pub fn normalize_total_enthalpy(
        self,
        physical: TotalEnthalpyJoules,
    ) -> Result<TotalEnthalpyJoules, ReactionExtentError> {
        TotalEnthalpyJoules::new(physical.joules() / self.physical_inventory_scale)
    }

    /// Reconstructs physical total enthalpy from normalized total enthalpy.
    pub fn denormalize_total_enthalpy(
        self,
        normalized: TotalEnthalpyJoules,
    ) -> Result<TotalEnthalpyJoules, ReactionExtentError> {
        TotalEnthalpyJoules::new(normalized.joules() * self.physical_inventory_scale)
    }

    /// Converts an internally normalized extensive residual to physical units.
    pub fn denormalize_extensive_error(
        self,
        normalized_error: f64,
    ) -> Result<f64, ReactionExtentError> {
        self.transform_finite_value(
            normalized_error,
            self.physical_inventory_scale,
            "normalized_extensive_error",
        )
    }

    /// Converts a physical absolute-mole threshold for normalized solver space.
    ///
    /// This applies to policies such as an absolute trace floor or `phase_eps`.
    /// It intentionally does not apply to TPD/hysteresis thresholds because
    /// those are intensive thermodynamic quantities.
    pub fn normalize_physical_mole_threshold(
        self,
        physical_threshold_moles: f64,
    ) -> Result<f64, ReactionExtentError> {
        self.transform_positive_threshold(
            physical_threshold_moles,
            1.0 / self.physical_inventory_scale,
            "physical_mole_threshold",
        )
    }

    /// Reconstructs a physical absolute-mole threshold from normalized space.
    pub fn denormalize_mole_threshold(
        self,
        normalized_threshold_moles: f64,
    ) -> Result<f64, ReactionExtentError> {
        self.transform_positive_threshold(
            normalized_threshold_moles,
            self.physical_inventory_scale,
            "normalized_mole_threshold",
        )
    }

    /// Converts a physical absolute energy tolerance into normalized joules.
    pub fn normalize_physical_energy_tolerance(
        self,
        physical_tolerance_joules: f64,
    ) -> Result<f64, ReactionExtentError> {
        self.transform_positive_threshold(
            physical_tolerance_joules,
            1.0 / self.physical_inventory_scale,
            "physical_energy_tolerance",
        )
    }

    /// Maps a trace policy while preserving its public physical-mole meaning.
    ///
    /// Relative fractions are dimensionless and therefore unchanged. Their
    /// absolute fallback cap is physical moles and must be normalized.
    pub fn normalize_trace_seed_policy(
        self,
        physical_policy: TraceSpeciesSeedPolicy,
    ) -> Result<TraceSpeciesSeedPolicy, ReactionExtentError> {
        match physical_policy {
            TraceSpeciesSeedPolicy::Absolute { floor } => Ok(TraceSpeciesSeedPolicy::Absolute {
                floor: self.normalize_physical_mole_threshold(floor)?,
            }),
            TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                fraction,
                minimum_floor,
            } => Ok(TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                fraction,
                minimum_floor: self.normalize_physical_mole_threshold(minimum_floor)?,
            }),
        }
    }

    fn transform_moles(
        self,
        values: &[f64],
        multiplier: f64,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        values
            .iter()
            .enumerate()
            .map(|(index, &value)| {
                if !value.is_finite() || value < 0.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "extensive_normalization_moles",
                        message: format!("moles[{index}] must be finite and non-negative"),
                    });
                }
                let transformed = value * multiplier;
                if !transformed.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "extensive_normalization_moles",
                        message: format!(
                            "moles[{index}] overflowed during extensive normalization"
                        ),
                    });
                }
                Ok(transformed)
            })
            .collect()
    }

    fn transform_positive_threshold(
        self,
        value: f64,
        multiplier: f64,
        field: &'static str,
    ) -> Result<f64, ReactionExtentError> {
        if !value.is_finite() || value <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field,
                message: "threshold must be finite and strictly positive".to_owned(),
            });
        }
        self.transform_finite_value(value, multiplier, field)
    }

    fn transform_log_moles(
        self,
        values: &[f64],
        offset: f64,
    ) -> Result<LogMolesInitialGuess, ReactionExtentError> {
        let transformed = values
            .iter()
            .enumerate()
            .map(|(index, value)| {
                let shifted = *value + offset;
                if !shifted.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "extensive_normalization_log_moles",
                        message: format!(
                            "log-moles[{index}] became non-finite during extensive normalization"
                        ),
                    });
                }
                Ok(shifted)
            })
            .collect::<Result<Vec<_>, _>>()?;
        LogMolesInitialGuess::new(transformed)
    }

    fn transform_finite_value(
        self,
        value: f64,
        multiplier: f64,
        field: &'static str,
    ) -> Result<f64, ReactionExtentError> {
        if !value.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field,
                message: "value must be finite".to_owned(),
            });
        }
        let transformed = value * multiplier;
        if !transformed.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field,
                message: "value overflowed during extensive normalization".to_owned(),
            });
        }
        Ok(transformed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(left: f64, right: f64) {
        let relative = (left - right).abs() / left.abs().max(right.abs()).max(1.0e-30);
        assert!(
            relative <= 1.0e-12,
            "left={left:e}, right={right:e}, relative={relative:e}"
        );
    }

    #[test]
    fn moles_round_trip_and_uniform_scale_invariance_are_explicit() {
        let base = vec![2.0, 0.5, 0.0, 1.5];
        let factor = 1.0e4;
        let scaled = base
            .iter()
            .map(|amount| amount * factor)
            .collect::<Vec<_>>();
        let base_normalization = ExtensiveNormalization::from_physical_moles(&base).unwrap();
        let scaled_normalization = ExtensiveNormalization::from_physical_moles(&scaled).unwrap();

        for ((&base_normalized, &scaled_normalized), &original) in base_normalization
            .normalize_moles(&base)
            .unwrap()
            .iter()
            .zip(
                scaled_normalization
                    .normalize_moles(&scaled)
                    .unwrap()
                    .iter(),
            )
            .zip(base.iter())
        {
            assert_close(base_normalized, scaled_normalized);
            assert_close(
                base_normalization
                    .denormalize_moles(&[base_normalized])
                    .unwrap()[0],
                original,
            );
        }
    }

    #[test]
    fn extensive_round_trips_cover_both_directions_and_reviewed_scales() {
        let physical_moles = [2.0e-6, 0.5e-6, 0.0, 1.5e-6];
        let physical_elements = [3.0e-6, 1.0e-6];
        let physical_phases = [2.5e-6, 1.5e-6];
        let physical_h = TotalEnthalpyJoules::new(-5.0e-6).unwrap();

        for scale in [1.0e-6, 1.0e-2, 1.0, 1.0e2, 1.0e6] {
            let normalization =
                ExtensiveNormalization::from_physical_inventory_scale(scale).unwrap();
            let normalized_moles = normalization.normalize_moles(&physical_moles).unwrap();
            let reconstructed_moles = normalization.denormalize_moles(&normalized_moles).unwrap();
            for (&actual, &expected) in reconstructed_moles.iter().zip(physical_moles.iter()) {
                assert_close(actual, expected);
            }

            let normalized_elements = normalization
                .normalize_element_totals(&physical_elements)
                .unwrap();
            assert_close(
                normalization
                    .denormalize_element_totals(&normalized_elements)
                    .unwrap()[0],
                physical_elements[0],
            );

            let normalized_phases = normalization
                .normalize_phase_totals(&physical_phases)
                .unwrap();
            assert_close(
                normalization
                    .denormalize_phase_totals(&normalized_phases)
                    .unwrap()[1],
                physical_phases[1],
            );

            let normalized_h = normalization.normalize_total_enthalpy(physical_h).unwrap();
            assert_close(
                normalization
                    .denormalize_total_enthalpy(normalized_h)
                    .unwrap()
                    .joules(),
                physical_h.joules(),
            );

            let physical_threshold = 1.0e-30;
            let normalized_threshold = normalization
                .normalize_physical_mole_threshold(physical_threshold)
                .unwrap();
            assert_close(
                normalization
                    .denormalize_mole_threshold(normalized_threshold)
                    .unwrap(),
                physical_threshold,
            );
        }
    }

    #[test]
    fn enthalpy_and_physical_thresholds_transform_with_the_same_scale() {
        let normalization = ExtensiveNormalization::from_physical_inventory_scale(2.5e4).unwrap();
        let normalized_h = normalization
            .normalize_total_enthalpy(TotalEnthalpyJoules::new(-5.0e6).unwrap())
            .unwrap();
        assert_close(normalized_h.joules(), -200.0);
        assert_close(
            normalization
                .denormalize_total_enthalpy(normalized_h)
                .unwrap()
                .joules(),
            -5.0e6,
        );

        let physical_floor = 1.0e-30;
        let normalized_floor = normalization
            .normalize_physical_mole_threshold(physical_floor)
            .unwrap();
        assert_close(
            normalization
                .denormalize_mole_threshold(normalized_floor)
                .unwrap(),
            physical_floor,
        );
        assert_close(
            normalization.denormalize_extensive_error(4.0e-5).unwrap(),
            1.0,
        );
    }

    #[test]
    fn log_mole_seed_round_trip_is_an_exact_coordinate_translation() {
        let normalization = ExtensiveNormalization::from_physical_inventory_scale(1.0e4).unwrap();
        let physical = LogMolesInitialGuess::new(vec![1.0e-30_f64.ln(), 2.5_f64.ln()]).unwrap();
        let normalized = normalization.normalize_log_mole_guess(&physical).unwrap();
        let reconstructed = normalization
            .denormalize_log_moles(normalized.as_slice())
            .unwrap();

        for (actual, expected) in reconstructed.as_slice().iter().zip(physical.as_slice()) {
            assert_close(*actual, *expected);
        }
    }

    #[test]
    fn trace_policy_scales_only_its_physical_mole_parts() {
        let normalization = ExtensiveNormalization::from_physical_inventory_scale(100.0).unwrap();
        match normalization
            .normalize_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1.0e-20 })
            .unwrap()
        {
            TraceSpeciesSeedPolicy::Absolute { floor } => assert_close(floor, 1.0e-22),
            policy => panic!("expected absolute trace policy, got {policy:?}"),
        }
        match normalization
            .normalize_trace_seed_policy(TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                fraction: 1.0e-8,
                minimum_floor: 1.0e-20,
            })
            .unwrap()
        {
            TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                fraction,
                minimum_floor,
            } => {
                assert_close(fraction, 1.0e-8);
                assert_close(minimum_floor, 1.0e-22);
            }
            policy => panic!("expected relative trace policy, got {policy:?}"),
        }
    }

    #[test]
    fn phase_and_element_totals_share_the_component_mole_mapping() {
        let normalization = ExtensiveNormalization::from_physical_inventory_scale(10.0).unwrap();
        for (&actual, &expected) in normalization
            .normalize_phase_totals(&[5.0, 2.0])
            .unwrap()
            .iter()
            .zip([0.5, 0.2].iter())
        {
            assert_close(actual, expected);
        }
        for (&actual, &expected) in normalization
            .normalize_element_totals(&[6.0, 3.0])
            .unwrap()
            .iter()
            .zip([0.6, 0.3].iter())
        {
            assert_close(actual, expected);
        }
        for (&actual, &expected) in normalization
            .denormalize_phase_totals(&[0.5, 0.2])
            .unwrap()
            .iter()
            .zip([5.0, 2.0].iter())
        {
            assert_close(actual, expected);
        }
        for (&actual, &expected) in normalization
            .denormalize_element_totals(&[0.6, 0.3])
            .unwrap()
            .iter()
            .zip([6.0, 3.0].iter())
        {
            assert_close(actual, expected);
        }
    }

    #[test]
    fn construction_rejects_nonphysical_scales_and_moles() {
        for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ExtensiveNormalization::from_physical_inventory_scale(scale).is_err());
        }
        assert!(ExtensiveNormalization::from_physical_moles(&[1.0, -1.0]).is_err());
        assert!(ExtensiveNormalization::from_physical_moles(&[0.0, 0.0]).is_err());
        assert!(ExtensiveNormalization::from_physical_moles(&[f64::NAN]).is_err());
    }

    #[test]
    fn extreme_scales_transform_without_silent_overflow() {
        let tiny = ExtensiveNormalization::from_physical_inventory_scale(1.0e-250).unwrap();
        let tiny_moles = [1.0e-250, 2.0e-250];
        assert_eq!(tiny.normalize_moles(&tiny_moles).unwrap(), vec![1.0, 2.0]);
        for (&actual, &expected) in tiny
            .denormalize_moles(&[1.0, 2.0])
            .unwrap()
            .iter()
            .zip(tiny_moles.iter())
        {
            assert_close(actual, expected);
        }

        let huge = ExtensiveNormalization::from_physical_inventory_scale(1.0e250).unwrap();
        assert_eq!(
            huge.normalize_moles(&[1.0e250, 2.0e250]).unwrap(),
            vec![1.0, 2.0]
        );
        for (&actual, &expected) in huge
            .denormalize_moles(&[1.0, 2.0])
            .unwrap()
            .iter()
            .zip([1.0e250, 2.0e250].iter())
        {
            assert_close(actual, expected);
        }

        assert!(huge.denormalize_moles(&[f64::MAX]).is_err());
        assert!(
            tiny.normalize_total_enthalpy(TotalEnthalpyJoules::new(1.0e308).unwrap())
                .is_err()
        );
    }
}
