//! Immutable, phase-aware equilibrium results.
//!
//! A successful nonlinear solve first produces a canonical
//! [`EquilibriumSolution`]. This module attaches that numeric snapshot to the
//! phase-qualified layout, thermochemical provenance, and backend trace that
//! created it. The result is therefore safe to query without reconstructing
//! parallel vectors or guessing whether a bare substance name is ambiguous.

use std::collections::BTreeMap;
use std::fmt;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::EquilibriumSolveReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    EquilibriumPhaseDescriptor, PhaseEquilibriumBuildReport, PhaseEquilibriumMetadata,
    PhaseEquilibriumSolutionBundle,
};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

/// One stable row in a multiphase result summary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiphaseEquilibriumSummaryRow {
    /// Logical section, for example `conditions`, `phase`, or `backend`.
    pub section: &'static str,
    /// Stable row key.
    pub label: String,
    /// Human-readable value.
    pub value: String,
}

impl fmt::Display for MultiphaseEquilibriumSummaryRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

/// Immutable result of a fixed-active-set phase equilibrium calculation.
///
/// The current bridge solves every declared phase as active. The status vector
/// explicitly records that fact instead of inferring activity from a tiny
/// numerical amount. A later phase-control adapter can construct the same
/// type with an accepted `PhaseSet` without changing lookup semantics.
/// Immutable result of a fixed-active-set phase equilibrium calculation.
///
/// The current bridge solves every declared phase as active. The status vector
/// explicitly records that fact instead of inferring activity from a tiny
/// numerical amount. A later phase-control adapter can construct the same
/// type with an accepted `PhaseSet` without changing lookup semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiphaseEquilibriumSolution {
    /// Phase-qualified layout and lookup provenance.
    metadata: PhaseEquilibriumMetadata,
    /// Immutable thermochemical preparation evidence.
    build_report: PhaseEquilibriumBuildReport,
    /// Accepted canonical solution (log-moles + physical moles).
    accepted_solution: EquilibriumSolution,
    /// Total moles per phase in declared phase order.
    phase_totals: Vec<f64>,
    /// Mole fractions per component in `SystemLayout` component order.
    mole_fractions: Vec<f64>,
    /// Lifecycle status for each declared phase (Active, Inactive, Excluded, etc.).
    phase_statuses: Vec<PhaseStatus>,
    /// Ordered backend cascade evidence for the accepted result.
    solve_report: EquilibriumSolveReport,
    /// Optional independent equilibrium-constant cross-validation status.
    keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
}

impl MultiphaseEquilibriumSolution {
    /// Converts one accepted fixed-active bridge bundle into a queryable
    /// phase-aware result.
    ///
    /// Construction repeats the cheap boundary invariants deliberately. A
    /// future caller cannot accidentally combine the accepted numerical vector
    /// with provenance or a layout belonging to another resolved system.
    pub fn from_fixed_active_bundle(
        bundle: PhaseEquilibriumSolutionBundle,
    ) -> Result<Self, ReactionExtentError> {
        let metadata = bundle.metadata().clone();
        let build_report = bundle.build_report().clone();
        let accepted_solution = bundle.solution().clone();
        let solve_report = bundle.solve_report().clone();
        let keq_validation_status = bundle.keq_validation_status().cloned();

        if metadata.layout_fingerprint() != build_report.layout_fingerprint() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "multiphase_solution_layout",
                message: "accepted result and build report have different layout fingerprints"
                    .to_string(),
            });
        }
        if accepted_solution.conditions() != build_report.conditions() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "multiphase_solution_conditions",
                message: "accepted result and build report have different thermodynamic conditions"
                    .to_string(),
            });
        }
        if metadata.components().len() != accepted_solution.moles().len()
            || metadata.components().len() != build_report.components().len()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "multiphase result has {} components, solution has {} moles, and build report has {} rows",
                metadata.components().len(),
                accepted_solution.moles().len(),
                build_report.components().len(),
            )));
        }

        let mut phase_totals = Vec::with_capacity(metadata.phases().len());
        let mut mole_fractions = vec![0.0; metadata.components().len()];
        for phase in metadata.phases() {
            let range = phase.component_range();
            let total = accepted_solution.moles()[range.clone()].iter().sum::<f64>();
            if !total.is_finite() || total <= 0.0 {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "multiphase_phase_total",
                    message: format!(
                        "phase {:?} has invalid accepted total {total:e}",
                        phase.id().as_option()
                    ),
                });
            }
            for component_index in range {
                mole_fractions[component_index] =
                    accepted_solution.moles()[component_index] / total;
            }
            phase_totals.push(total);
        }

        Ok(Self {
            phase_statuses: vec![PhaseStatus::Active; metadata.phases().len()],
            metadata,
            build_report,
            accepted_solution,
            phase_totals,
            mole_fractions,
            solve_report,
            keq_validation_status,
        })
    }

    /// Fixed pressure-temperature conditions retained by the accepted snapshot.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.accepted_solution.conditions()
    }

    /// Canonical phase-qualified layout and provenance identity.
    pub fn metadata(&self) -> &PhaseEquilibriumMetadata {
        &self.metadata
    }

    /// Layout fingerprint that must match the originating resolved system.
    pub fn layout_fingerprint(&self) -> u64 {
        self.metadata.layout_fingerprint()
    }

    /// Original immutable standard-state and lookup evidence.
    pub fn build_report(&self) -> &PhaseEquilibriumBuildReport {
        &self.build_report
    }

    /// Canonical accepted log-mole/physical-mole snapshot.
    pub fn accepted_solution(&self) -> &EquilibriumSolution {
        &self.accepted_solution
    }

    /// All accepted component mole values in exact `SystemLayout` order.
    pub fn component_moles(&self) -> &[f64] {
        self.accepted_solution.moles()
    }

    /// Finds the accepted amount of one fully-qualified component.
    pub fn moles_for(&self, component: &PhaseComponentId) -> Option<f64> {
        self.metadata
            .component_index(component)
            .map(|index| self.accepted_solution.moles()[index])
    }

    /// Finds the local mole fraction of one fully-qualified component.
    pub fn mole_fraction_for(&self, component: &PhaseComponentId) -> Option<f64> {
        self.metadata
            .component_index(component)
            .map(|index| self.mole_fractions[index])
    }

    /// Ordered phase descriptors used by the accepted snapshot.
    pub fn phases(&self) -> &[EquilibriumPhaseDescriptor] {
        self.metadata.phases()
    }

    /// Accepted total mole amount in one semantic phase.
    pub fn phase_total(&self, phase: &PhaseId) -> Option<f64> {
        self.metadata
            .phase_index(phase)
            .map(|index| self.phase_totals[index.index()])
    }

    /// Explicit lifecycle state for one semantic phase.
    pub fn phase_status(&self, phase: &PhaseId) -> Option<PhaseStatus> {
        self.metadata
            .phase_index(phase)
            .map(|index| self.phase_statuses[index.index()])
    }

    /// Aggregates accepted amounts by bare substance as an explicit derived
    /// view. It is intentionally not used for solver identity because a name
    /// can occur in several physical phases.
    pub fn aggregate_moles_by_substance(&self) -> BTreeMap<String, f64> {
        let mut totals = BTreeMap::new();
        for (descriptor, &moles) in self
            .metadata
            .components()
            .iter()
            .zip(self.accepted_solution.moles())
        {
            *totals
                .entry(descriptor.substance().to_string())
                .or_insert(0.0) += moles;
        }
        totals
    }

    /// Complete backend cascade evidence for the accepted solve.
    pub fn solve_report(&self) -> &EquilibriumSolveReport {
        &self.solve_report
    }

    /// Optional independent equilibrium-constant validation evidence.
    pub fn keq_validation_status(&self) -> Option<&EquilibriumConstantCrossValidationStatus> {
        self.keq_validation_status.as_ref()
    }

    /// Stable summary rows for CLI, snapshots, and a future GUI.
    pub fn summary_rows(&self) -> Vec<MultiphaseEquilibriumSummaryRow> {
        let mut rows = vec![
            MultiphaseEquilibriumSummaryRow {
                section: "conditions",
                label: "temperature_k".to_string(),
                value: format!("{:.6}", self.conditions().temperature()),
            },
            MultiphaseEquilibriumSummaryRow {
                section: "conditions",
                label: "pressure_pa".to_string(),
                value: format!("{:.6}", self.conditions().pressure()),
            },
            MultiphaseEquilibriumSummaryRow {
                section: "layout",
                label: "fingerprint".to_string(),
                value: self.layout_fingerprint().to_string(),
            },
            MultiphaseEquilibriumSummaryRow {
                section: "backend",
                label: "accepted".to_string(),
                value: format!("{:?}", self.solve_report.accepted_backend),
            },
            MultiphaseEquilibriumSummaryRow {
                section: "validation",
                label: "residual_l2_norm".to_string(),
                value: format!(
                    "{:.6e}",
                    self.accepted_solution.validation().residual_l2_norm
                ),
            },
        ];

        for (index, phase) in self.metadata.phases().iter().enumerate() {
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "phase",
                label: phase
                    .id()
                    .as_option()
                    .clone()
                    .unwrap_or_else(|| "single".to_string()),
                value: format!(
                    "total={:.6e}, status={:?}",
                    self.phase_totals[index], self.phase_statuses[index]
                ),
            });
        }
        for (index, component) in self.metadata.components().iter().enumerate() {
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "component",
                label: component.label(),
                value: format!(
                    "moles={:.6e}, x={:.6e}",
                    self.accepted_solution.moles()[index],
                    self.mole_fractions[index]
                ),
            });
        }
        rows
    }
}

impl fmt::Display for MultiphaseEquilibriumSolution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}
