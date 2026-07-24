//! Cross-validation reports comparing the canonical equilibrium solution with
//! the independent K_eq solution.
//!
//! This module is intentionally narrow. It does not solve either problem; it
//! only compares already computed accepted results and records the largest
//! species, fraction, and thermodynamic disagreements.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::{
    EquilibriumConstantActivityModel, EquilibriumConstantProblem, MOLAR_GAS_CONSTANT,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolveResult;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumSolution,
};
use std::fmt;

/// Numerical tolerances for comparing canonical and K_eq solutions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumConstantCrossValidationTolerances {
    /// Largest allowed absolute mole-number mismatch.
    pub max_abs_species_mole_delta: f64,
    /// Largest allowed absolute mole-fraction mismatch.
    pub max_abs_species_fraction_delta: f64,
    /// Largest allowed absolute total-Gibbs mismatch.
    pub max_abs_total_gibbs_delta: f64,
}

impl Default for EquilibriumConstantCrossValidationTolerances {
    fn default() -> Self {
        Self {
            max_abs_species_mole_delta: 1e-5,
            max_abs_species_fraction_delta: 1e-5,
            max_abs_total_gibbs_delta: 1e-6,
        }
    }
}

impl EquilibriumConstantCrossValidationTolerances {
    /// Validates the tolerance bundle before report construction.
    pub fn validate(self) -> Result<Self, ReactionExtentError> {
        for (field, value) in [
            (
                "max_abs_species_mole_delta",
                self.max_abs_species_mole_delta,
            ),
            (
                "max_abs_species_fraction_delta",
                self.max_abs_species_fraction_delta,
            ),
            ("max_abs_total_gibbs_delta", self.max_abs_total_gibbs_delta),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "equilibrium_constant_cross_validation_tolerances",
                    message: format!("{field} must be finite and strictly positive"),
                });
            }
        }
        Ok(self)
    }
}

/// Per-species comparison between canonical and independent solutions.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantSpeciesComparison {
    /// Deterministic species index.
    pub species_index: usize,
    /// Species name in report order.
    pub species_name: String,
    /// Canonical equilibrium mole number.
    pub canonical_moles: f64,
    /// Independent K_eq equilibrium mole number.
    pub keq_moles: f64,
    /// Absolute mole-number mismatch.
    pub abs_mole_delta: f64,
    /// Canonical mole fraction.
    pub canonical_mole_fraction: f64,
    /// Independent K_eq mole fraction.
    pub keq_mole_fraction: f64,
    /// Absolute mole-fraction mismatch.
    pub abs_mole_fraction_delta: f64,
}

impl fmt::Display for EquilibriumConstantSpeciesComparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: canonical={:.6e}, keq={:.6e}, |Δn|={:.6e}, |Δx|={:.6e}",
            self.species_name,
            self.canonical_moles,
            self.keq_moles,
            self.abs_mole_delta,
            self.abs_mole_fraction_delta,
        )
    }
}

/// Cross-validation report comparing the canonical and independent solutions.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantCrossValidationReport {
    /// Thermodynamic conditions under which both solutions were evaluated.
    pub conditions: EquilibriumConditions,
    /// One row per species.
    pub species: Vec<EquilibriumConstantSpeciesComparison>,
    /// Canonical total Gibbs energy under the current ideal-gas model.
    pub canonical_total_gibbs: f64,
    /// Independent K_eq total Gibbs energy under the current ideal-gas model.
    pub keq_total_gibbs: f64,
    /// Absolute Gibbs-energy difference.
    pub abs_total_gibbs_delta: f64,
    /// Largest absolute mole mismatch among species.
    pub max_abs_species_mole_delta: f64,
    /// Largest absolute mole-fraction mismatch among species.
    pub max_abs_species_fraction_delta: f64,
    /// Largest absolute independent `ln(Q) - ln(K)` mismatch.
    pub max_abs_log_residual: f64,
    /// Whether the comparison satisfies the configured tolerances.
    pub accepted: bool,
}

/// One stable summary row for the cross-validation report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumConstantCrossValidationRow {
    /// Logical section name.
    pub section: &'static str,
    /// Stable row label.
    pub label: String,
    /// Human-readable value.
    pub value: String,
}

impl fmt::Display for EquilibriumConstantCrossValidationRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

impl EquilibriumConstantCrossValidationReport {
    /// Returns stable summary rows for CLI output and snapshot tests.
    pub fn summary_rows(&self) -> Vec<EquilibriumConstantCrossValidationRow> {
        let mut rows = vec![
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "temperature".to_string(),
                value: format!("{:.6}", self.conditions.temperature()),
            },
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "pressure".to_string(),
                value: format!("{:.6}", self.conditions.pressure()),
            },
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "accepted".to_string(),
                value: self.accepted.to_string(),
            },
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "max_abs_species_mole_delta".to_string(),
                value: format!("{:.6e}", self.max_abs_species_mole_delta),
            },
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "max_abs_species_fraction_delta".to_string(),
                value: format!("{:.6e}", self.max_abs_species_fraction_delta),
            },
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "max_abs_log_residual".to_string(),
                value: format!("{:.6e}", self.max_abs_log_residual),
            },
            EquilibriumConstantCrossValidationRow {
                section: "comparison",
                label: "abs_total_gibbs_delta".to_string(),
                value: format!("{:.6e}", self.abs_total_gibbs_delta),
            },
        ];

        for species in &self.species {
            rows.push(EquilibriumConstantCrossValidationRow {
                section: "species",
                label: species.species_name.clone(),
                value: species.to_string(),
            });
        }

        rows
    }
}

impl fmt::Display for EquilibriumConstantCrossValidationReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

/// High-level status for a cross-validation request.
#[derive(Debug, Clone, PartialEq)]
pub enum EquilibriumConstantCrossValidationStatus {
    /// Both solvers completed and the resulting compositions were compared.
    Compared(EquilibriumConstantCrossValidationReport),
    /// The canonical solver failed before a comparison could be made.
    CanonicalFailed {
        /// Canonical solver error summary.
        summary: String,
    },
    /// The independent validator failed numerically or physically.
    ValidatorFailed {
        /// Validator error summary.
        summary: String,
    },
    /// The independent validator was not applicable for the current request.
    ValidatorNotApplicable {
        /// Human-readable reason.
        message: String,
    },
}

impl fmt::Display for EquilibriumConstantCrossValidationStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Compared(report) => write!(f, "compared\n{report}"),
            Self::CanonicalFailed { summary } => {
                write!(f, "canonical solver failed: {summary}")
            }
            Self::ValidatorFailed { summary } => {
                write!(f, "validator failed: {summary}")
            }
            Self::ValidatorNotApplicable { message } => {
                write!(f, "validator not applicable: {message}")
            }
        }
    }
}

/// Compares one canonical equilibrium solution against one independent K_eq result.
pub fn compare_equilibrium_constant_solutions(
    problem: &EquilibriumConstantProblem,
    canonical: &EquilibriumSolution,
    keq: &EquilibriumConstantSolveResult,
    tolerances: EquilibriumConstantCrossValidationTolerances,
) -> Result<EquilibriumConstantCrossValidationReport, ReactionExtentError> {
    let tolerances = tolerances.validate()?;
    let species = problem.basis().species();
    if canonical.moles().len() != species.len() || keq.moles.len() != species.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "comparison has {} species, canonical moles have {}, and K_eq moles have {}",
            species.len(),
            canonical.moles().len(),
            keq.moles.len()
        )));
    }
    if !matches!(
        problem.activity_model(),
        EquilibriumConstantActivityModel::IdealGas
    ) {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "equilibrium_constant_cross_validation",
            message: "cross-validation currently supports only the ideal-gas activity model"
                .to_string(),
        });
    }

    let canonical_total_gibbs = total_gibbs(problem, canonical.conditions(), canonical.moles())?;
    let keq_total_gibbs = total_gibbs(problem, problem.conditions(), &keq.moles)?;
    let abs_total_gibbs_delta = (canonical_total_gibbs - keq_total_gibbs).abs();

    let canonical_total: f64 = canonical.moles().iter().sum();
    let keq_total: f64 = keq.moles.iter().sum();
    if canonical_total <= 0.0 || keq_total <= 0.0 {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "cross_validation_solution",
            message: "comparison requires positive total mole numbers".to_string(),
        });
    }

    let mut comparisons = Vec::with_capacity(species.len());
    let mut max_abs_species_mole_delta = 0.0_f64;
    let mut max_abs_species_fraction_delta = 0.0_f64;
    for (index, name) in species.iter().enumerate() {
        let canonical_moles = canonical.moles()[index];
        let keq_moles = keq.moles[index];
        let abs_mole_delta = (canonical_moles - keq_moles).abs();
        let canonical_mole_fraction = canonical_moles / canonical_total;
        let keq_mole_fraction = keq_moles / keq_total;
        let abs_mole_fraction_delta = (canonical_mole_fraction - keq_mole_fraction).abs();
        max_abs_species_mole_delta = max_abs_species_mole_delta.max(abs_mole_delta);
        max_abs_species_fraction_delta =
            max_abs_species_fraction_delta.max(abs_mole_fraction_delta);
        comparisons.push(EquilibriumConstantSpeciesComparison {
            species_index: index,
            species_name: name.clone(),
            canonical_moles,
            keq_moles,
            abs_mole_delta,
            canonical_mole_fraction,
            keq_mole_fraction,
            abs_mole_fraction_delta,
        });
    }

    let max_abs_log_residual = keq.validation.max_abs_log_residual;
    Ok(EquilibriumConstantCrossValidationReport {
        conditions: problem.conditions(),
        species: comparisons,
        canonical_total_gibbs,
        keq_total_gibbs,
        abs_total_gibbs_delta,
        max_abs_species_mole_delta,
        max_abs_species_fraction_delta,
        max_abs_log_residual,
        accepted: max_abs_species_mole_delta <= tolerances.max_abs_species_mole_delta
            && max_abs_species_fraction_delta <= tolerances.max_abs_species_fraction_delta
            && abs_total_gibbs_delta <= tolerances.max_abs_total_gibbs_delta
            && canonical.validation().max_abs_element_balance_error <= 1e-8
            && keq.validation.accepted,
    })
}

/// Wraps the comparison step and classifies which side failed when a report
/// cannot be produced.
pub fn classify_equilibrium_constant_cross_validation(
    problem: &EquilibriumConstantProblem,
    canonical: Result<EquilibriumSolution, ReactionExtentError>,
    validator: Result<Option<EquilibriumConstantSolveResult>, ReactionExtentError>,
    tolerances: EquilibriumConstantCrossValidationTolerances,
) -> Result<EquilibriumConstantCrossValidationStatus, ReactionExtentError> {
    match canonical {
        Err(error) => Ok(EquilibriumConstantCrossValidationStatus::CanonicalFailed {
            summary: error.to_string(),
        }),
        Ok(canonical) => match validator {
            Ok(None) => Ok(
                EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable {
                    message: "independent validator was disabled or skipped".to_string(),
                },
            ),
            Err(ReactionExtentError::ValidationNotApplicable { message, .. }) => {
                Ok(EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable { message })
            }
            Err(error) => Ok(EquilibriumConstantCrossValidationStatus::ValidatorFailed {
                summary: error.to_string(),
            }),
            Ok(Some(keq)) => Ok(EquilibriumConstantCrossValidationStatus::Compared(
                compare_equilibrium_constant_solutions(problem, &canonical, &keq, tolerances)?,
            )),
        },
    }
}

fn total_gibbs(
    problem: &EquilibriumConstantProblem,
    conditions: EquilibriumConditions,
    moles: &[f64],
) -> Result<f64, ReactionExtentError> {
    if moles.len() != problem.standard_gibbs().len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "total Gibbs evaluation has {} mole values for {} species",
            moles.len(),
            problem.standard_gibbs().len()
        )));
    }
    let temperature = conditions.temperature();
    let pressure_ratio = conditions.pressure() / conditions.reference_pressure();
    let total: f64 = moles.iter().sum();
    if total <= 0.0 || !total.is_finite() {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "cross_validation_moles",
            message: "total Gibbs evaluation requires positive total moles".to_string(),
        });
    }
    let mut total_gibbs = 0.0;
    for (index, (&mole, gibbs)) in moles.iter().zip(problem.standard_gibbs()).enumerate() {
        if mole <= 0.0 || !mole.is_finite() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "cross_validation_moles",
                message: format!("species {index} has an invalid mole number {mole}"),
            });
        }
        let x_i = mole / total;
        let mu_i =
            gibbs(temperature) + MOLAR_GAS_CONSTANT * temperature * (x_i * pressure_ratio).ln();
        total_gibbs += mole * mu_i;
    }
    Ok(total_gibbs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::{
        EquilibriumConstantActivityModel, EquilibriumConstantProblem, MOLAR_GAS_CONSTANT,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolver;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
        EquilibriumConstantValidationTolerances, validate_equilibrium_constants,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Phase;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::PhaseKind;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_reaction_basis::{
        ReactionBasisTolerances, ValidatedReactionBasis,
    };
    use nalgebra::DMatrix;
    use std::rc::Rc;

    fn constant_gibbs(value: f64) -> GibbsFn {
        Rc::new(move |_| value)
    }

    fn comparison_fixture() -> (
        EquilibriumConstantProblem,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolveResult,
    ){
        let target_k = 0.75_f64;
        let temperature = 1_500.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let basis = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 1, &[-1.0, 2.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0],
            vec![constant_gibbs(-delta_g), constant_gibbs(0.0)],
            EquilibriumConditions::new(temperature, 2.0 * 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();
        let canonical_problem = EquilibriumProblem::new(
            vec!["A2".to_string(), "A".to_string()],
            vec![1.0, 0.0],
            LogMolesInitialGuess::from_initial_moles(&[1.0, 0.0]).unwrap(),
            DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            vec![constant_gibbs(-delta_g), constant_gibbs(0.0)],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            EquilibriumConditions::new(temperature, 2.0 * 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let canonical = crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles::solve_problem(
            canonical_problem,
        )
        .unwrap();
        let keq = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap();
        (problem, canonical, keq)
    }

    fn water_gas_shift_fixture() -> (
        EquilibriumConstantProblem,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolveResult,
    ){
        let target_k = 4.0_f64;
        let temperature = 1_400.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let basis = ValidatedReactionBasis::new(
            vec![
                "CO".to_string(),
                "H2O".to_string(),
                "CO2".to_string(),
                "H2".to_string(),
            ],
            &DMatrix::from_row_slice(
                4,
                3,
                &[
                    1.0, 0.0, 1.0, // CO
                    0.0, 2.0, 1.0, // H2O
                    1.0, 0.0, 2.0, // CO2
                    0.0, 2.0, 0.0, // H2
                ],
            ),
            DMatrix::from_column_slice(4, 1, &[-1.0, -1.0, 1.0, 1.0]),
            3,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 1.0, 0.0, 0.0],
            vec![
                constant_gibbs(0.0),
                constant_gibbs(-delta_g),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();
        let canonical_problem = EquilibriumProblem::new(
            vec![
                "CO".to_string(),
                "H2O".to_string(),
                "CO2".to_string(),
                "H2".to_string(),
            ],
            vec![1.0, 1.0, 0.0, 0.0],
            LogMolesInitialGuess::from_initial_moles(&[1.0, 1.0, 0.0, 0.0]).unwrap(),
            DMatrix::from_row_slice(
                4,
                3,
                &[
                    1.0, 0.0, 1.0, // CO
                    0.0, 2.0, 1.0, // H2O
                    1.0, 0.0, 2.0, // CO2
                    0.0, 2.0, 0.0, // H2
                ],
            ),
            vec![
                constant_gibbs(0.0),
                constant_gibbs(-delta_g),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1, 2, 3],
            }],
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let canonical = crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles::solve_problem(
            canonical_problem,
        )
        .unwrap();
        let keq = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap();
        (problem, canonical, keq)
    }

    fn inert_dilution_fixture() -> (
        EquilibriumConstantProblem,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolveResult,
    ){
        // Reaction: 2H2 + O2 <=> 2H2O, with N2 as an inert diluent.
        let target_k = 5.0_f64;
        let temperature = 1_800.0;
        let pressure = 2.0 * 101_325.0;
        let reference_pressure = 101_325.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let basis = ValidatedReactionBasis::new(
            vec![
                "H2".to_string(),
                "O2".to_string(),
                "H2O".to_string(),
                "N2".to_string(),
            ],
            &DMatrix::from_row_slice(
                4,
                3,
                &[
                    2.0, 0.0, 0.0, // H2
                    0.0, 2.0, 0.0, // O2
                    2.0, 1.0, 0.0, // H2O
                    0.0, 0.0, 2.0, // N2
                ],
            ),
            DMatrix::from_column_slice(4, 1, &[-2.0, -1.0, 2.0, 0.0]),
            3,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![2.0, 1.0, 0.0, 5.0],
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(-delta_g / 2.0),
                constant_gibbs(0.0),
            ],
            EquilibriumConditions::new(temperature, pressure, reference_pressure).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();
        let canonical_problem = EquilibriumProblem::new(
            vec![
                "H2".to_string(),
                "O2".to_string(),
                "H2O".to_string(),
                "N2".to_string(),
            ],
            vec![2.0, 1.0, 0.0, 5.0],
            LogMolesInitialGuess::from_initial_moles(&[2.0, 1.0, 0.0, 5.0]).unwrap(),
            DMatrix::from_row_slice(
                4,
                3,
                &[
                    2.0, 0.0, 0.0, // H2
                    0.0, 2.0, 0.0, // O2
                    2.0, 1.0, 0.0, // H2O
                    0.0, 0.0, 2.0, // N2
                ],
            ),
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(-delta_g / 2.0),
                constant_gibbs(0.0),
            ],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1, 2, 3],
            }],
            EquilibriumConditions::new(temperature, pressure, reference_pressure).unwrap(),
        )
        .unwrap();
        let canonical = crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles::solve_problem(
            canonical_problem,
        )
        .unwrap();
        let keq = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap();
        (problem, canonical, keq)
    }

    #[test]
    fn cross_validation_report_exposes_stable_rows_and_display() {
        let (problem, canonical, keq) = comparison_fixture();
        let report = compare_equilibrium_constant_solutions(
            &problem,
            &canonical,
            &keq,
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        let rows = report.summary_rows();
        assert!(
            rows.iter()
                .any(|row| row.section == "comparison" && row.label == "accepted")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "species" && row.label == "A2")
        );
        let rendered = format!("{report}");
        assert!(rendered.contains("[comparison] accepted = true"));
        assert!(rendered.contains("[species] A2 = A2:"));
    }

    #[test]
    fn cross_validation_status_display_is_human_readable() {
        let status = EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable {
            message: "no independent validator for this case".to_string(),
        };
        assert!(format!("{status}").contains("validator not applicable"));
    }

    fn methane_combustion_fixture_with_initial_moles(
        initial_moles: [f64; 4],
        temperature: f64,
        target_k: f64,
    ) -> (
        EquilibriumConstantProblem,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolveResult,
    ){
        // Species: CH4, O2, CO2, H2O
        // Reaction: CH4 + 2O2 <=> CO2 + 2H2O
        // Delta-nu = 0, so the quotient is especially stable to validate.
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let basis = ValidatedReactionBasis::new(
            vec![
                "CH4".to_string(),
                "O2".to_string(),
                "CO2".to_string(),
                "H2O".to_string(),
            ],
            &DMatrix::from_row_slice(
                4,
                3,
                &[
                    1.0, 4.0, 0.0, // CH4
                    0.0, 0.0, 2.0, // O2
                    1.0, 0.0, 2.0, // CO2
                    0.0, 2.0, 1.0, // H2O
                ],
            ),
            DMatrix::from_column_slice(4, 1, &[-1.0, -2.0, 1.0, 2.0]),
            3,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            initial_moles.to_vec(),
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(-delta_g),
                constant_gibbs(0.0),
            ],
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();
        let canonical_problem = EquilibriumProblem::new(
            vec![
                "CH4".to_string(),
                "O2".to_string(),
                "CO2".to_string(),
                "H2O".to_string(),
            ],
            initial_moles.to_vec(),
            LogMolesInitialGuess::from_initial_moles(&initial_moles).unwrap(),
            DMatrix::from_row_slice(
                4,
                3,
                &[
                    1.0, 4.0, 0.0, // CH4
                    0.0, 0.0, 2.0, // O2
                    1.0, 0.0, 2.0, // CO2
                    0.0, 2.0, 1.0, // H2O
                ],
            ),
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(-delta_g),
                constant_gibbs(0.0),
            ],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1, 2, 3],
            }],
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let canonical = crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles::solve_problem(
            canonical_problem,
        )
        .unwrap();
        let keq = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap();
        (problem, canonical, keq)
    }

    fn methane_combustion_fixture() -> (
        EquilibriumConstantProblem,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolveResult,
    ){
        methane_combustion_fixture_with_initial_moles([1.0, 2.0, 0.0, 0.0], 1_600.0, 8.0)
    }

    #[test]
    fn canonical_and_independent_solutions_compare_as_expected() {
        let (problem, canonical, keq) = comparison_fixture();
        let report = compare_equilibrium_constant_solutions(
            &problem,
            &canonical,
            &keq,
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        assert!(report.accepted);
        assert!(report.max_abs_species_mole_delta < 1e-5);
        assert!(report.max_abs_species_fraction_delta < 1e-5);
        assert!(report.abs_total_gibbs_delta < 1e-6);
        assert!(report.max_abs_log_residual < 1e-8);
        assert_eq!(report.species.len(), 2);
    }

    #[test]
    fn cross_validation_reports_a_large_disagreement() {
        let (problem, canonical, mut keq) = comparison_fixture();
        keq.moles[0] += 0.1;
        keq.moles[1] -= 0.1;
        keq.validation = validate_equilibrium_constants(
            &problem,
            &keq.moles,
            EquilibriumConstantValidationTolerances::default(),
        )
        .unwrap();
        keq.validation.accepted = false;

        let report = compare_equilibrium_constant_solutions(
            &problem,
            &canonical,
            &keq,
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        assert!(!report.accepted);
        assert!(report.max_abs_species_mole_delta > 0.0);
        assert!(report.abs_total_gibbs_delta > 0.0);
    }

    #[test]
    fn water_gas_shift_fixture_cross_validates_independently() {
        let (problem, canonical, keq) = water_gas_shift_fixture();
        let tolerances = EquilibriumConstantCrossValidationTolerances {
            max_abs_species_mole_delta: 1e-5,
            max_abs_species_fraction_delta: 1e-5,
            max_abs_total_gibbs_delta: 1e-5,
        };
        let report =
            compare_equilibrium_constant_solutions(&problem, &canonical, &keq, tolerances).unwrap();

        assert!(report.accepted);
        assert_eq!(report.species.len(), 4);
        assert!(report.max_abs_species_mole_delta < 1e-5);
        assert!(report.max_abs_species_fraction_delta < 1e-5);
        assert!(report.abs_total_gibbs_delta < 1e-5);
    }

    #[test]
    fn solve_if_applicable_returns_none_for_unsupported_multi_reaction_systems() {
        let basis = ValidatedReactionBasis::new(
            vec![
                "NO".to_string(),
                "N2".to_string(),
                "O2".to_string(),
                "NO2".to_string(),
            ],
            &DMatrix::from_row_slice(
                4,
                2,
                &[
                    1.0, 1.0, // NO
                    2.0, 0.0, // N2
                    0.0, 2.0, // O2
                    1.0, 2.0, // NO2
                ],
            ),
            DMatrix::from_column_slice(
                4,
                2,
                &[
                    -2.0, 1.0, 1.0, 0.0, // 2NO <=> N2 + O2
                    -4.0, 1.0, 0.0, 2.0, // 4NO <=> N2 + 2NO2
                ],
            ),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0, 0.0, 0.0],
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            EquilibriumConditions::new(1_500.0, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let solver = EquilibriumConstantSolver::default();
        assert!(solver.solve_if_applicable(&problem).unwrap().is_none());
    }

    #[test]
    fn solve_if_applicable_returns_none_when_disabled_even_for_supported_systems() {
        let (problem, _, _) = comparison_fixture();
        let solver = EquilibriumConstantSolver {
            mode: super::super::equilibrium_constant_solver::EquilibriumConstantSolverMode::Off,
            ..EquilibriumConstantSolver::default()
        };

        assert!(solver.solve_if_applicable(&problem).unwrap().is_none());
    }

    #[test]
    fn cross_validation_classification_reports_canonical_failure() {
        let (problem, _, keq) = comparison_fixture();
        let status = classify_equilibrium_constant_cross_validation(
            &problem,
            Err(ReactionExtentError::InvalidProblem {
                field: "canonical_solver",
                message: "failed before comparison".to_string(),
            }),
            Ok(Some(keq)),
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        assert!(matches!(
            status,
            EquilibriumConstantCrossValidationStatus::CanonicalFailed { .. }
        ));
    }

    #[test]
    fn cross_validation_classification_reports_validator_not_applicable() {
        let (problem, canonical, _) = comparison_fixture();
        let status = classify_equilibrium_constant_cross_validation(
            &problem,
            Ok(canonical),
            Ok(None),
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        assert!(matches!(
            status,
            EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable { .. }
        ));
    }

    #[test]
    fn cross_validation_classification_reports_validator_failure() {
        let (problem, canonical, _) = comparison_fixture();
        let status = classify_equilibrium_constant_cross_validation(
            &problem,
            Ok(canonical),
            Err(ReactionExtentError::InvalidProblem {
                field: "validator_solver",
                message: "failed before comparison".to_string(),
            }),
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        assert!(matches!(
            status,
            EquilibriumConstantCrossValidationStatus::ValidatorFailed { .. }
        ));
    }

    #[test]
    fn cross_validation_classification_reports_compared_result() {
        let (problem, canonical, keq) = comparison_fixture();
        let status = classify_equilibrium_constant_cross_validation(
            &problem,
            Ok(canonical),
            Ok(Some(keq)),
            EquilibriumConstantCrossValidationTolerances::default(),
        )
        .unwrap();

        assert!(matches!(
            status,
            EquilibriumConstantCrossValidationStatus::Compared(_)
        ));
    }

    #[test]
    fn inert_dilution_fixture_cross_validates_independently() {
        let (problem, canonical, keq) = inert_dilution_fixture();
        let tolerances = EquilibriumConstantCrossValidationTolerances {
            max_abs_species_mole_delta: 2e-5,
            max_abs_species_fraction_delta: 1e-5,
            max_abs_total_gibbs_delta: 1e-5,
        };
        let report =
            compare_equilibrium_constant_solutions(&problem, &canonical, &keq, tolerances).unwrap();

        assert!(report.accepted);
        assert_eq!(report.species.len(), 4);
        assert!(report.max_abs_species_mole_delta < 2e-5);
        assert!(report.max_abs_species_fraction_delta < 1e-5);
        assert!(report.abs_total_gibbs_delta < 1e-5);
    }

    #[test]
    fn methane_combustion_fixture_cross_validates_independently() {
        let (problem, canonical, keq) = methane_combustion_fixture();
        let tolerances = EquilibriumConstantCrossValidationTolerances {
            max_abs_species_mole_delta: 2e-5,
            max_abs_species_fraction_delta: 1e-5,
            max_abs_total_gibbs_delta: 1e-5,
        };
        let report =
            compare_equilibrium_constant_solutions(&problem, &canonical, &keq, tolerances).unwrap();

        assert!(report.accepted);
        assert_eq!(report.species.len(), 4);
        assert!(report.max_abs_species_mole_delta < 2e-5);
        assert!(report.max_abs_species_fraction_delta < 1e-5);
        assert!(report.abs_total_gibbs_delta < 1e-5);
    }

    #[test]
    fn methane_lean_and_rich_fixtures_cross_validate_independently() {
        let cases = [
            ("lean", [1.0, 4.0, 0.0, 0.0]),
            ("rich", [2.0, 3.0, 0.0, 0.0]),
        ];

        for (label, initial_moles) in cases {
            let (problem, canonical, keq) =
                methane_combustion_fixture_with_initial_moles(initial_moles, 1_600.0, 8.0);
            let report = compare_equilibrium_constant_solutions(
                &problem,
                &canonical,
                &keq,
                EquilibriumConstantCrossValidationTolerances {
                    max_abs_species_mole_delta: 3e-5,
                    max_abs_species_fraction_delta: 1e-5,
                    max_abs_total_gibbs_delta: 2e-5,
                },
            )
            .unwrap();

            assert!(
                report.accepted,
                "{label} methane case should cross-validate"
            );
            assert_eq!(report.species.len(), 4);
        }
    }

    #[test]
    fn methane_combustion_temperature_sweep_cross_validates_independently() {
        for temperature in [1_200.0, 1_600.0, 2_000.0] {
            let (problem, canonical, keq) = methane_combustion_fixture_with_initial_moles(
                [1.0, 2.0, 0.0, 0.0],
                temperature,
                8.0,
            );
            let report = compare_equilibrium_constant_solutions(
                &problem,
                &canonical,
                &keq,
                EquilibriumConstantCrossValidationTolerances {
                    max_abs_species_mole_delta: 3e-5,
                    max_abs_species_fraction_delta: 1e-5,
                    max_abs_total_gibbs_delta: 2e-5,
                },
            )
            .unwrap();

            assert!(
                report.accepted,
                "temperature {temperature} K should cross-validate"
            );
            assert_eq!(report.species.len(), 4);
        }
    }
}
