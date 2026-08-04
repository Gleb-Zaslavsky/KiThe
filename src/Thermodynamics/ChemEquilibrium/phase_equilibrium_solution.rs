//! Immutable, phase-aware equilibrium results.
//!
//! A successful nonlinear solve first produces a canonical
//! [`EquilibriumSolution`]. This module attaches that numeric snapshot to the
//! phase-qualified layout, thermochemical provenance, and backend trace that
//! created it. The result is therefore safe to query without reconstructing
//! parallel vectors or guessing whether a bare substance name is ambiguous.

use std::collections::BTreeMap;
use std::fmt;

use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, MultiStartSolveReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::{
    EquilibriumTimingCollector, EquilibriumTimingReport, EquilibriumTimingStage,
};
#[cfg(test)]
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    MultiphaseAcceptanceReport, PhaseControlledSolveReport, PhaseStatus,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    EquilibriumPhaseDescriptor, PhaseEquilibriumBuildReport, PhaseEquilibriumMetadata,
    PhaseEquilibriumSolutionBundle,
};

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
/// Every published value is tied to the bridge metadata that produced it.
/// Fixed active-set and bounded phase-control workflows differ only in their
/// phase-status and acceptance evidence, not in component lookup semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiphaseEquilibriumSolution {
    /// Phase-qualified layout and lookup provenance.
    metadata: PhaseEquilibriumMetadata,
    /// Immutable thermochemical preparation evidence.
    build_report: PhaseEquilibriumBuildReport,
    /// Accepted canonical numerical solution (log-moles + positive coordinates).
    accepted_solution: EquilibriumSolution,
    /// Published physical component moles in `SystemLayout` component order.
    /// Inactive, excluded, and disappeared phases are represented by zeroes;
    /// their positive solver trace coordinates never leak into this view.
    physical_component_moles: Vec<f64>,
    /// Published physical total moles per phase in declared phase order.
    phase_totals: Vec<f64>,
    /// Published physical mole fractions per component in `SystemLayout`
    /// component order. Components of inactive phases have fraction zero.
    mole_fractions: Vec<f64>,
    /// Numerical phase totals reconstructed from the accepted log-mole vector.
    /// This is diagnostic evidence for the nonlinear layer, not physical output.
    numerical_phase_totals: Vec<f64>,
    /// Lifecycle status for each declared phase (Active, Inactive, Excluded, etc.).
    phase_statuses: Vec<PhaseStatus>,
    /// Ordered backend cascade evidence for the accepted result.
    solve_report: EquilibriumSolveReport,
    /// Explicit multi-start seed comparison evidence, when requested.
    multi_start_report: Option<MultiStartSolveReport>,
    /// Optional independent equilibrium-constant cross-validation status.
    keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
    /// Bounded active-set transition evidence when phase control was used.
    phase_control_report: Option<PhaseControlledSolveReport>,
    /// Combined numerical and complementarity evidence when phase control was
    /// used to reach a stable phase set.
    acceptance_report: Option<MultiphaseAcceptanceReport>,
    /// Optional stage timing collected while this immutable result was built.
    timing: EquilibriumTimingReport,
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
        let multi_start_report = bundle.multi_start_report().cloned();
        let keq_validation_status = bundle.keq_validation_status().cloned();
        let timing = *bundle.timing_report();

        Self::from_parts(
            metadata,
            build_report,
            accepted_solution,
            solve_report,
            multi_start_report,
            keq_validation_status,
            None,
            None,
            None,
            timing,
        )
    }

    /// Builds the public fixed-active view from an internally prepared
    /// formulation whose accepted temperature differs from its build seed.
    ///
    /// This is crate-private so only the phase bridge can bind numerical
    /// evidence to matching lookup provenance and layout metadata.
    pub(crate) fn from_fixed_active_parts(
        metadata: PhaseEquilibriumMetadata,
        build_report: PhaseEquilibriumBuildReport,
        accepted_solution: EquilibriumSolution,
        solve_report: EquilibriumSolveReport,
        timing: EquilibriumTimingReport,
    ) -> Result<Self, ReactionExtentError> {
        Self::from_parts(
            metadata,
            build_report,
            accepted_solution,
            solve_report,
            None,
            None,
            None,
            None,
            None,
            timing,
        )
    }

    /// Converts accepted bounded phase-control evidence into the same public
    /// immutable result model used by fixed active-set solves.
    ///
    /// `phase_statuses` must remain aligned to the canonical bridge phase
    /// descriptors. This prevents an outer-loop report for one phase layout
    /// from being attached to moles or provenance belonging to another.
    pub(crate) fn from_phase_control_parts(
        metadata: PhaseEquilibriumMetadata,
        build_report: PhaseEquilibriumBuildReport,
        accepted_solution: EquilibriumSolution,
        solve_report: EquilibriumSolveReport,
        keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
        phase_control_report: PhaseControlledSolveReport,
        acceptance_report: MultiphaseAcceptanceReport,
        phase_statuses: Vec<PhaseStatus>,
        timing: EquilibriumTimingReport,
    ) -> Result<Self, ReactionExtentError> {
        if acceptance_report.phase_control != phase_control_report {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "multiphase_acceptance",
                message: "acceptance report does not retain the published phase-control report"
                    .to_string(),
            });
        }
        Self::from_parts(
            metadata,
            build_report,
            accepted_solution,
            solve_report,
            None,
            keq_validation_status,
            Some(phase_control_report),
            Some(acceptance_report),
            Some(phase_statuses),
            timing,
        )
    }

    fn from_parts(
        metadata: PhaseEquilibriumMetadata,
        build_report: PhaseEquilibriumBuildReport,
        accepted_solution: EquilibriumSolution,
        solve_report: EquilibriumSolveReport,
        multi_start_report: Option<MultiStartSolveReport>,
        keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
        phase_control_report: Option<PhaseControlledSolveReport>,
        acceptance_report: Option<MultiphaseAcceptanceReport>,
        phase_statuses: Option<Vec<PhaseStatus>>,
        timing: EquilibriumTimingReport,
    ) -> Result<Self, ReactionExtentError> {
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
        let phase_statuses =
            phase_statuses.unwrap_or_else(|| vec![PhaseStatus::Active; metadata.phases().len()]);
        if phase_statuses.len() != metadata.phases().len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "multiphase result has {} phase descriptors but {} statuses",
                metadata.phases().len(),
                phase_statuses.len(),
            )));
        }

        let mut phase_totals = Vec::with_capacity(metadata.phases().len());
        let mut numerical_phase_totals = Vec::with_capacity(metadata.phases().len());
        let mut physical_component_moles = vec![0.0; metadata.components().len()];
        let mut mole_fractions = vec![0.0; metadata.components().len()];
        for phase in metadata.phases() {
            let range = phase.component_range();
            let numerical_total = accepted_solution.moles()[range.clone()].iter().sum::<f64>();
            if !numerical_total.is_finite() || numerical_total <= 0.0 {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "multiphase_phase_total",
                    message: format!(
                        "phase {:?} has invalid accepted numerical total {numerical_total:e}",
                        phase.id().as_option()
                    ),
                });
            }
            numerical_phase_totals.push(numerical_total);

            let phase_index = phase.index().index();
            if phase_statuses[phase_index].is_active() {
                for component_index in range.clone() {
                    physical_component_moles[component_index] =
                        accepted_solution.moles()[component_index];
                    mole_fractions[component_index] =
                        accepted_solution.moles()[component_index] / numerical_total;
                }
                phase_totals.push(numerical_total);
            } else {
                // The solver must keep a positive coordinate for log-space
                // evaluation, but an absent physical phase has zero inventory.
                phase_totals.push(0.0);
            }
        }

        Ok(Self {
            metadata,
            build_report,
            accepted_solution,
            physical_component_moles,
            phase_totals,
            mole_fractions,
            numerical_phase_totals,
            phase_statuses,
            solve_report,
            multi_start_report,
            keq_validation_status,
            phase_control_report,
            acceptance_report,
            timing,
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

    /// Optional stage timing collected while this immutable result was built.
    pub fn timing_report(&self) -> &EquilibriumTimingReport {
        &self.timing
    }

    /// Updates only the wall-clock total after an enclosing public operation
    /// has completed. Stage measurements remain unchanged.
    pub(crate) fn with_timing_total(mut self, total: std::time::Duration) -> Self {
        if self.timing.enabled() {
            let mut collector = EquilibriumTimingCollector::from_report(self.timing);
            collector.set_total(total);
            self.timing = collector.finish();
        }
        self
    }

    /// Adds one enclosing-operation stage without rebuilding the result.
    pub(crate) fn with_timing_stage(
        mut self,
        stage: EquilibriumTimingStage,
        duration: std::time::Duration,
    ) -> Self {
        if self.timing.enabled() {
            let mut collector = EquilibriumTimingCollector::from_report(self.timing);
            collector.record(stage, duration);
            self.timing = collector.finish();
        }
        self
    }

    /// Canonical accepted numerical log-mole/positive-mole snapshot.
    ///
    /// This is the solver coordinate record. It may contain trace-floor
    /// amounts for absent phases and is therefore not the physical result view.
    pub fn accepted_solution(&self) -> &EquilibriumSolution {
        &self.accepted_solution
    }

    /// Replaces only the validation payload for a test-built malformed
    /// worker result.
    ///
    /// The production constructors never expose this operation. GUI tests use
    /// it to model a particularly important publication boundary: a worker
    /// may hand the UI an object labelled as accepted while its copied
    /// validation evidence no longer describes the physical mole vector.
    #[cfg(test)]
    pub(crate) fn with_validation_for_test(
        mut self,
        validation: EquilibriumCandidateReport,
    ) -> Self {
        self.accepted_solution = self.accepted_solution.with_validation_for_test(validation);
        self
    }

    /// Published physical component amounts in exact `SystemLayout` order.
    ///
    /// Inactive, excluded, and disappeared phases are exactly zero here even
    /// though their numerical log-mole coordinates remain positive internally.
    pub fn component_moles(&self) -> &[f64] {
        &self.physical_component_moles
    }

    /// Numerical positive-mole coordinates retained for diagnostics.
    pub fn numerical_component_moles(&self) -> &[f64] {
        self.accepted_solution.moles()
    }

    /// Finds the accepted amount of one fully-qualified component.
    pub fn moles_for(&self, component: &PhaseComponentId) -> Option<f64> {
        self.metadata
            .component_index(component)
            .map(|index| self.physical_component_moles[index])
    }

    /// Finds the numerical positive-mole coordinate of one component.
    pub fn numerical_moles_for(&self, component: &PhaseComponentId) -> Option<f64> {
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

    /// Numerical positive-mole total retained for one phase.
    pub fn numerical_phase_total(&self, phase: &PhaseId) -> Option<f64> {
        self.metadata
            .phase_index(phase)
            .map(|index| self.numerical_phase_totals[index.index()])
    }

    /// Explicit lifecycle state for one semantic phase.
    pub fn phase_status(&self, phase: &PhaseId) -> Option<PhaseStatus> {
        self.metadata
            .phase_index(phase)
            .map(|index| self.phase_statuses[index.index()])
    }

    /// Aggregates published physical amounts by bare substance as an explicit derived
    /// view. It is intentionally not used for solver identity because a name
    /// can occur in several physical phases.
    pub fn aggregate_moles_by_substance(&self) -> BTreeMap<String, f64> {
        let mut totals = BTreeMap::new();
        for (descriptor, moles) in self
            .metadata
            .components()
            .iter()
            .zip(self.physical_component_moles.iter().copied())
        {
            *totals
                .entry(descriptor.substance().to_string())
                .or_insert(0.0) += moles;
        }
        totals
    }

    /// Aggregates numerical positive-mole coordinates by bare substance.
    ///
    /// This diagnostic view is useful when explaining a trace-floor or a
    /// backend acceptance decision; it must not be used as physical inventory.
    pub fn aggregate_numerical_moles_by_substance(&self) -> BTreeMap<String, f64> {
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

    /// Total number of started nonlinear backend attempts represented by this
    /// accepted solution.
    ///
    /// A multi-start solve owns one backend trace per seed. In that case the
    /// seed-level aggregate is authoritative; otherwise the ordinary cascade
    /// report is sufficient.
    pub fn started_backend_attempts(&self) -> usize {
        self.multi_start_report.as_ref().map_or_else(
            || self.solve_report.started_attempt_count(),
            MultiStartSolveReport::started_backend_attempts,
        )
    }

    /// Total nonlinear iterations represented by this accepted solve,
    /// including every continuation multi-start seed when present.
    pub fn nonlinear_iterations(&self) -> usize {
        self.multi_start_report.as_ref().map_or_else(
            || self.solve_report.nonlinear_iterations(),
            MultiStartSolveReport::nonlinear_iterations,
        )
    }

    /// Number of accepted phase transitions represented by this solution.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_report
            .as_ref()
            .map_or(0, |report| report.transitions.len())
    }

    /// Explicit multi-start seed comparison evidence, when requested.
    pub fn multi_start_report(&self) -> Option<&MultiStartSolveReport> {
        self.multi_start_report.as_ref()
    }

    /// Optional independent equilibrium-constant validation evidence.
    pub fn keq_validation_status(&self) -> Option<&EquilibriumConstantCrossValidationStatus> {
        self.keq_validation_status.as_ref()
    }

    /// Bounded phase-control evidence when this result was solved through the
    /// active-set outer loop.
    pub fn phase_control_report(&self) -> Option<&PhaseControlledSolveReport> {
        self.phase_control_report.as_ref()
    }

    /// Final complementarity/validation evidence when phase control was used.
    pub fn acceptance_report(&self) -> Option<&MultiphaseAcceptanceReport> {
        self.acceptance_report.as_ref()
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
        if let Some(report) = &self.phase_control_report {
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "phase_control",
                label: "iterations".to_string(),
                value: report.iterations.to_string(),
            });
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "phase_control",
                label: "transitions".to_string(),
                value: report.transitions.len().to_string(),
            });
        }
        if let Some(report) = &self.acceptance_report {
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "acceptance",
                label: "complementarity_satisfied".to_string(),
                value: report.complementarity.satisfied.to_string(),
            });
        }
        if let Some(status) = &self.keq_validation_status {
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "keq_validation",
                label: "status".to_string(),
                value: match status {
                    EquilibriumConstantCrossValidationStatus::Compared(report) => {
                        format!("compared accepted={}", report.accepted)
                    }
                    EquilibriumConstantCrossValidationStatus::CanonicalFailed { .. } => {
                        "canonical_failed".to_string()
                    }
                    EquilibriumConstantCrossValidationStatus::ValidatorFailed { .. } => {
                        "validator_failed".to_string()
                    }
                    EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable { .. } => {
                        "not_applicable".to_string()
                    }
                },
            });
        }
        for (index, component) in self.metadata.components().iter().enumerate() {
            rows.push(MultiphaseEquilibriumSummaryRow {
                section: "component",
                label: component.label(),
                value: format!(
                    "moles={:.6e}, x={:.6e}",
                    self.physical_component_moles[index], self.mole_fractions[index]
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

#[cfg(test)]
mod tests {
    //! Boundary regressions for immutable phase-aware result publication.
    //!
    //! The fixture uses two independently valid local NASA gas requests.  The
    //! test then attempts to combine metadata from one request with the build
    //! evidence from the other, which must fail before a public result exists.

    use std::collections::HashMap;

    use crate::Thermodynamics::phase_layout::PhaseId;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, TraceSpeciesSeedPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
        build_phase_equilibrium_problem, PhaseEquilibriumBuildRequest,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{PhaseSpec, ResolvedPhaseSystem};
    use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};

    use super::MultiphaseEquilibriumSolution;

    fn accepted_local_nasa_gas(phase_name: &str) -> crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumSolutionBundle{
        let phase_id = PhaseId::new(Some(phase_name.to_string()));
        let phase = PhaseSpec::ideal_gas(
            phase_id.clone(),
            vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()],
        )
        .unwrap();
        let mut data = SubsData::new();
        data.substances = vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()];
        data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        data.search_substances().unwrap();
        data.parse_all_thermal_coeffs().unwrap();
        let resolved = ResolvedPhaseSystem::new(
            vec![phase],
            HashMap::from([(Some(phase_name.to_string()), data)]),
        )
        .unwrap();
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
        let composition =
            MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();

        build_phase_equilibrium_problem(
            PhaseEquilibriumBuildRequest::new(
                &resolved,
                EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap(),
                composition,
                TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
                Default::default(),
            )
            .unwrap(),
        )
        .unwrap()
        .solve()
        .unwrap()
    }

    #[test]
    fn reconstruction_rejects_metadata_and_build_report_from_different_layouts() {
        let gas = accepted_local_nasa_gas("gas");
        let other = accepted_local_nasa_gas("other_gas");
        assert_ne!(
            gas.metadata().layout_fingerprint(),
            other.build_report().layout_fingerprint()
        );

        let error = MultiphaseEquilibriumSolution::from_parts(
            gas.metadata().clone(),
            other.build_report().clone(),
            gas.solution().clone(),
            gas.solve_report().clone(),
            None,
            gas.keq_validation_status().cloned(),
            None,
            None,
            None,
            crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::
                EquilibriumTimingReport::default(),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError::InvalidCandidate {
                field: "multiphase_solution_layout",
                ..
            }
        ));
    }
}
