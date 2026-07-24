//! Typed boundary between resolved thermodynamic data and equilibrium algorithms.
//!
//! The phase subsystem owns lookup, physical-state selection, provenance, and
//! the canonical [`SystemLayout`]. The equilibrium subsystem owns activity
//! models, numerical coordinates, residuals, and nonlinear solver policy. This
//! module is the one-way bridge between those responsibilities:
//!
//! ```text
//! ResolvedPhaseSystem
//!     |  PhaseSpec + SubsData + provenance + SystemLayout
//!     v
//! PhaseEquilibriumBuildRequest
//!     |
//!     v
//! PhaseEquilibriumMetadata
//!     |  qualified components + ordered phases + stable indices
//!     v
//! EquilibriumProblem                 (built in the next adapter layer)
//! ```
//!
//! Structural validation precedes thermochemical extraction. The final builder
//! evaluates `G0(T)` on phase-local working copies, so missing or invalid data
//! cannot mutate the resolved input or publish a partial numerical problem.

use std::collections::{BTreeSet, HashMap};
use std::ops::Range;
use std::rc::Rc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_component::EquilibriumComponentDescriptor;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumLogMoles, EquilibriumSolverSettings, GibbsFn, Phase,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, EquilibriumSolution, LogMolesInitialGuess,
    TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::EquilibriumSolveReport;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::User_PhaseOrSolution::{
    PhaseModel, ResolvedPhaseSystem, ResolvedPhaseSystemReport,
};
use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances2::SearchSummaryRow;
use crate::Thermodynamics::phase_layout::{
    PhaseComponentId, PhaseId as SemanticPhaseId, SystemLayout,
};
use crate::Thermodynamics::physical_state::PhysicalState;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

/// Explicit version of the phase-model contract accepted by the bridge.
///
/// A policy value is carried by every build request so extending the phase
/// subsystem cannot silently make the equilibrium solver accept a model whose
/// chemical-potential equation has not been implemented and tested.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SupportedPhaseModelPolicy {
    /// One ideal-gas phase and any number of one-component pure condensed
    /// phases at fixed pressure and temperature.
    #[default]
    FixedPressureTemperatureV1,
}

/// Ordered solver projection of one semantic thermodynamic phase.
///
/// Each phase in the resolved system is projected into an `EquilibriumPhaseDescriptor`
/// that carries both the semantic identity (from the data subsystem) and the
/// dense index used inside ordered equilibrium arrays.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumPhaseDescriptor {
    /// Semantic phase identity retained from the data subsystem.
    id: SemanticPhaseId,
    /// Dense phase index used only inside ordered equilibrium arrays.
    index: PhaseIndex,
    /// Physical state (Gas, Liquid, Solid) used to resolve phase records.
    physical_state: PhysicalState,
    /// Declared thermodynamic model (IdealGas, PureCondensed, etc.).
    phase_model: PhaseModel,
    /// Activity law used by the equilibrium formulation.
    activity_model: PhaseActivityModel,
    /// Contiguous canonical component range owned by this phase.
    component_range: Range<usize>,
}

impl EquilibriumPhaseDescriptor {
    /// Semantic phase identity retained from the data subsystem.
    pub fn id(&self) -> &SemanticPhaseId {
        &self.id
    }

    /// Dense phase index used only inside ordered equilibrium arrays.
    pub fn index(&self) -> PhaseIndex {
        self.index
    }

    /// Physical state used to resolve the phase records.
    pub fn physical_state(&self) -> PhysicalState {
        self.physical_state
    }

    /// Declared thermodynamic model for this phase.
    pub fn phase_model(&self) -> PhaseModel {
        self.phase_model
    }

    /// Activity law used by the equilibrium formulation.
    pub fn activity_model(&self) -> PhaseActivityModel {
        self.activity_model
    }

    /// Contiguous canonical component range owned by this phase.
    pub fn component_range(&self) -> Range<usize> {
        self.component_range.clone()
    }
}

/// Immutable structural projection shared by bridge construction and results.
///
/// The metadata owns the canonical layout and all derived index maps. Callers
/// therefore never have to keep a `Vec<String>`, phase ranges, and provenance
/// aligned by convention.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhaseEquilibriumMetadata {
    /// Canonical ordered phase/component layout from the resolved system.
    layout: SystemLayout,
    /// Stable process-independent fingerprint of the declared layout.
    layout_fingerprint: u64,
    /// All phase-qualified components in canonical solver order.
    components: Vec<EquilibriumComponentDescriptor>,
    /// All phases in canonical solver order with their component ranges.
    phases: Vec<EquilibriumPhaseDescriptor>,
    /// Maps `PhaseComponentId` → dense component index for O(1) lookup.
    component_indices: HashMap<PhaseComponentId, usize>,
    /// Maps semantic `PhaseId` → dense `PhaseIndex` for O(1) lookup.
    phase_indices: HashMap<SemanticPhaseId, PhaseIndex>,
    /// Provenance evidence from the resolved phase system.
    provenance: ResolvedPhaseSystemReport,
}

impl PhaseEquilibriumMetadata {
    /// Projects one resolved phase system into equilibrium-owned descriptors.
    ///
    /// This operation is pure and transactional: all validation and allocation
    /// completes before a metadata value is returned.
    pub fn from_resolved(
        resolved: &ResolvedPhaseSystem,
        policy: SupportedPhaseModelPolicy,
    ) -> Result<Self, ReactionExtentError> {
        match policy {
            SupportedPhaseModelPolicy::FixedPressureTemperatureV1 => {}
        }

        let multiphase_layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
        if multiphase_layout.system_layout() != resolved.layout() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "resolved_phase_layout",
                message:
                    "resolved phase data and phase specifications use different component order"
                        .to_string(),
            });
        }

        let layout = resolved.layout().clone();
        let phase_count = resolved.phase_specs().len();
        let mut phases = Vec::with_capacity(phase_count);
        let mut components = Vec::with_capacity(layout.component_count());
        let mut phase_indices = HashMap::with_capacity(phase_count);
        let mut component_indices = HashMap::with_capacity(layout.component_count());

        for (phase_position, spec) in resolved.phase_specs().iter().enumerate() {
            let id = spec.id().clone();
            let index = PhaseIndex::new(phase_position, phase_count)?;
            let component_range = layout.phase_component_range(&id).cloned().ok_or_else(|| {
                ReactionExtentError::InvalidProblem {
                    field: "resolved_phase_layout",
                    message: format!("phase {:?} has no component range", id.as_option()),
                }
            })?;
            let activity_model = activity_model_for(spec.model());

            if phase_indices.insert(id.clone(), index).is_some() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "resolved_phase_layout",
                    message: format!("duplicate phase identity {:?}", id.as_option()),
                });
            }

            for component_position in component_range.clone() {
                let id = layout.components()[component_position].clone();
                if component_indices
                    .insert(id.clone(), component_position)
                    .is_some()
                {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "resolved_phase_layout",
                        message: format!("duplicate component identity '{}'", id.label()),
                    });
                }
                components.push(EquilibriumComponentDescriptor::new(
                    id,
                    spec.physical_state(),
                    spec.model(),
                    activity_model,
                ));
            }

            phases.push(EquilibriumPhaseDescriptor {
                id,
                index,
                physical_state: spec.physical_state(),
                phase_model: spec.model(),
                activity_model,
                component_range,
            });
        }

        if components.len() != layout.component_count() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "resolved_phase_layout",
                message: "not every resolved component belongs to a declared phase".to_string(),
            });
        }

        Ok(Self {
            layout,
            layout_fingerprint: multiphase_layout.fingerprint(),
            components,
            phases,
            component_indices,
            phase_indices,
            provenance: resolved.report().clone(),
        })
    }

    /// Canonical ordered phase/component layout.
    pub fn layout(&self) -> &SystemLayout {
        &self.layout
    }

    /// Stable fingerprint used to reject stale compositions and solutions.
    pub fn layout_fingerprint(&self) -> u64 {
        self.layout_fingerprint
    }

    /// Component descriptors in exact solver-vector order.
    pub fn components(&self) -> &[EquilibriumComponentDescriptor] {
        &self.components
    }

    /// Phase descriptors in canonical semantic phase order.
    pub fn phases(&self) -> &[EquilibriumPhaseDescriptor] {
        &self.phases
    }

    /// Finds the solver position of one qualified component.
    pub fn component_index(&self, id: &PhaseComponentId) -> Option<usize> {
        self.component_indices.get(id).copied()
    }

    /// Finds the dense solver index of one semantic phase.
    pub fn phase_index(&self, id: &SemanticPhaseId) -> Option<PhaseIndex> {
        self.phase_indices.get(id).copied()
    }

    /// Immutable lookup report retained from phase-system resolution.
    pub fn provenance(&self) -> &ResolvedPhaseSystemReport {
        &self.provenance
    }
}

/// Complete structural request for building one fixed-`P,T` equilibrium problem.
///
/// Construction validates that the physical composition and resolved records
/// share exactly the same phase-qualified layout. Standard-state Gibbs models
/// are intentionally extracted only by the next adapter pass.
#[derive(Debug)]
pub struct PhaseEquilibriumBuildRequest<'a> {
    resolved: &'a ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    initial_composition: MultiphaseInitialComposition,
    trace_seed_policy: TraceSpeciesSeedPolicy,
    model_policy: SupportedPhaseModelPolicy,
    metadata: PhaseEquilibriumMetadata,
}

impl<'a> PhaseEquilibriumBuildRequest<'a> {
    /// Creates a validated request without evaluating thermochemistry.
    pub fn new(
        resolved: &'a ResolvedPhaseSystem,
        conditions: EquilibriumConditions,
        initial_composition: MultiphaseInitialComposition,
        trace_seed_policy: TraceSpeciesSeedPolicy,
        model_policy: SupportedPhaseModelPolicy,
    ) -> Result<Self, ReactionExtentError> {
        let multiphase_layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
        initial_composition.validate_for(&multiphase_layout)?;
        let metadata = PhaseEquilibriumMetadata::from_resolved(resolved, model_policy)?;

        Ok(Self {
            resolved,
            conditions,
            initial_composition,
            trace_seed_policy,
            model_policy,
            metadata,
        })
    }

    /// Resolved phase-local thermochemical records.
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        self.resolved
    }

    /// Fixed thermodynamic conditions for the requested solve.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Validated physical mole numbers before trace seeding.
    pub fn initial_composition(&self) -> &MultiphaseInitialComposition {
        &self.initial_composition
    }

    /// Explicit numerical policy for zero-mole log-coordinate seeds.
    pub fn trace_seed_policy(&self) -> TraceSpeciesSeedPolicy {
        self.trace_seed_policy
    }

    /// Supported phase-model contract selected for this request.
    pub fn model_policy(&self) -> SupportedPhaseModelPolicy {
        self.model_policy
    }

    /// Immutable structural projection prepared transactionally at construction.
    pub fn metadata(&self) -> &PhaseEquilibriumMetadata {
        &self.metadata
    }
}

/// Immutable provenance and standard-state evidence for one bridge component.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumBridgeComponentReport {
    /// Phase-qualified component described by this report row.
    component: EquilibriumComponentDescriptor,
    /// Physical initial amount before trace seeding, in moles.
    initial_moles: f64,
    /// Standard-state Gibbs energy `G0(T)` checked during preparation, in J/mol.
    standard_gibbs_at_conditions: f64,
    /// Local thermochemical lookup row that supplied this component.
    thermo_source: SearchSummaryRow,
}

impl EquilibriumBridgeComponentReport {
    /// Phase-qualified component described by this report row.
    pub fn component(&self) -> &EquilibriumComponentDescriptor {
        &self.component
    }

    /// Physical initial amount before trace seeding.
    pub fn initial_moles(&self) -> f64 {
        self.initial_moles
    }

    /// Standard-state Gibbs energy `G0(T)` checked during preparation, in J/mol.
    pub fn standard_gibbs_at_conditions(&self) -> f64 {
        self.standard_gibbs_at_conditions
    }

    /// Local thermochemical lookup row that supplied this component.
    pub fn thermo_source(&self) -> &SearchSummaryRow {
        &self.thermo_source
    }
}

/// Read-only evidence emitted when a phase system becomes a solver problem.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseEquilibriumBuildReport {
    /// Fixed pressure-temperature conditions used for standard-state checks.
    conditions: EquilibriumConditions,
    /// Layout fingerprint shared by the input composition and bridge metadata.
    layout_fingerprint: u64,
    /// Element names in exact column order of the solver composition matrix.
    element_labels: Vec<String>,
    /// Conserved physical element totals before numerical trace seeding.
    element_totals: Vec<f64>,
    /// Component reports in exact `SystemLayout` order.
    components: Vec<EquilibriumBridgeComponentReport>,
}

impl PhaseEquilibriumBuildReport {
    /// Fixed pressure-temperature conditions used for standard-state checks.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Layout fingerprint shared by the input composition and bridge metadata.
    pub fn layout_fingerprint(&self) -> u64 {
        self.layout_fingerprint
    }

    /// Element names in exact column order of the solver composition matrix.
    pub fn element_labels(&self) -> &[String] {
        &self.element_labels
    }

    /// Conserved physical element totals before numerical trace seeding.
    pub fn element_totals(&self) -> &[f64] {
        &self.element_totals
    }

    /// Component reports in exact `SystemLayout` order.
    pub fn components(&self) -> &[EquilibriumBridgeComponentReport] {
        &self.components
    }
}

/// Complete immutable product of `ResolvedPhaseSystem -> EquilibriumProblem` preparation.
///
/// The bundle prevents callers from accidentally pairing a numerical problem
/// with a layout, source report, or initial inventory from another resolution.
pub struct PhaseEquilibriumProblemBundle {
    /// Canonical numerical problem ready for the solver.
    problem: EquilibriumProblem,
    /// Phase-qualified layout and lookup metadata.
    metadata: PhaseEquilibriumMetadata,
    /// Immutable thermochemical preparation evidence.
    report: PhaseEquilibriumBuildReport,
    /// RST builds its own Jacobian from these expressions. Keeping them beside
    /// the numeric G0 closures avoids reopening mutable SubsData during solve.
    symbolic_standard_gibbs: Vec<Expr>,
}

impl PhaseEquilibriumProblemBundle {
    /// Numerical problem ready for canonical preparation and solver selection.
    pub fn problem(&self) -> &EquilibriumProblem {
        &self.problem
    }

    /// Phase-qualified ordering and immutable lookup provenance.
    pub fn metadata(&self) -> &PhaseEquilibriumMetadata {
        &self.metadata
    }

    /// Inspectable preparation evidence without invoking a solver.
    pub fn report(&self) -> &PhaseEquilibriumBuildReport {
        &self.report
    }

    /// Consumes the bundle when a solver must own the numerical problem.
    pub fn into_problem(self) -> EquilibriumProblem {
        self.problem
    }

    /// Solves the complete declared phase set through the canonical backend
    /// cascade using its default numerical settings.
    ///
    /// This is deliberately a consuming operation. A failed numerical attempt
    /// cannot leave a mutable partial solution attached to the prepared
    /// thermochemical bundle, while a successful attempt returns all lookup,
    /// preparation, validation, and backend evidence as one immutable value.
    pub fn solve(self) -> Result<PhaseEquilibriumSolutionBundle, ReactionExtentError> {
        self.solve_with(|_| {})
    }

    /// Solves the complete declared phase set after configuring only numerical
    /// controls.
    ///
    /// The callback cannot alter the resolved phase layout, initial physical
    /// inventory, standard-state Gibbs functions, or provenance. It receives
    /// only the canonical solver settings, matching
    /// [`EquilibriumLogMoles::solve_problem_with`].
    pub fn solve_with<F>(
        self,
        configure: F,
    ) -> Result<PhaseEquilibriumSolutionBundle, ReactionExtentError>
    where
        F: FnOnce(&mut EquilibriumSolverSettings),
    {
        let mut solver = EquilibriumLogMoles::from_problem(self.problem)?;
        solver.gibbs_sym = self.symbolic_standard_gibbs;
        configure(&mut solver.solver_settings);
        solver.solve()?;

        let solution = solver.accepted_solution()?;
        let solve_report = solver.last_solve_report.clone().ok_or_else(|| {
            ReactionExtentError::InvalidCandidate {
                field: "phase_equilibrium_solution",
                message: "accepted bridge solve did not publish a backend report".to_string(),
            }
        })?;

        Ok(PhaseEquilibriumSolutionBundle {
            metadata: self.metadata,
            build_report: self.report,
            solution,
            solve_report,
            keq_validation_status: solver.last_keq_validation_status.clone(),
        })
    }
}

/// Immutable accepted result of a phase-system equilibrium solve.
///
/// The result owns the same metadata and preparation report that produced the
/// numerical problem. This keeps phase-qualified identities, source
/// provenance, elemental totals, the accepted solution, and the backend trace
/// in one transactional publication unit.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseEquilibriumSolutionBundle {
    /// Phase-qualified layout and lookup provenance for this accepted result.
    metadata: PhaseEquilibriumMetadata,
    /// Immutable thermochemical preparation evidence for this accepted result.
    build_report: PhaseEquilibriumBuildReport,
    /// Accepted physical/log-mole solution in `SystemLayout` component order.
    solution: EquilibriumSolution,
    /// Ordered backend cascade evidence for the accepted result.
    solve_report: EquilibriumSolveReport,
    /// Optional independent equilibrium-constant validation evidence.
    keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
}

impl PhaseEquilibriumSolutionBundle {
    /// Phase-qualified layout and lookup provenance for this accepted result.
    pub fn metadata(&self) -> &PhaseEquilibriumMetadata {
        &self.metadata
    }

    /// Immutable thermochemical preparation evidence for this accepted result.
    pub fn build_report(&self) -> &PhaseEquilibriumBuildReport {
        &self.build_report
    }

    /// Accepted physical/log-mole solution in `SystemLayout` component order.
    pub fn solution(&self) -> &EquilibriumSolution {
        &self.solution
    }

    /// Ordered backend cascade evidence for the accepted result.
    pub fn solve_report(&self) -> &EquilibriumSolveReport {
        &self.solve_report
    }

    /// Optional independent equilibrium-constant validation evidence.
    pub fn keq_validation_status(&self) -> Option<&EquilibriumConstantCrossValidationStatus> {
        self.keq_validation_status.as_ref()
    }

    /// Converts this accepted fixed-active-set result into the public
    /// phase-aware query model.
    pub fn into_multiphase_solution(
        self,
    ) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
        MultiphaseEquilibriumSolution::from_fixed_active_bundle(self)
    }
}

/// Builds the canonical numerical equilibrium input from resolved phase data.
///
/// Every phase payload is cloned into a private working copy while extracting
/// `G0(T)` closures and elemental composition. The resolved phase system stays
/// immutable, and no bundle is returned until all components, provenance rows,
/// element coefficients, and fixed-temperature standard-state values validate.
pub fn build_phase_equilibrium_problem(
    request: PhaseEquilibriumBuildRequest<'_>,
) -> Result<PhaseEquilibriumProblemBundle, ReactionExtentError> {
    let metadata = request.metadata.clone();
    let conditions = request.conditions;
    let mut phase_data = HashMap::with_capacity(metadata.phases().len());
    let mut all_elements = BTreeSet::new();

    for phase in metadata.phases() {
        let payload = request
            .resolved
            .phase_data()
            .get(phase.id().as_option())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "resolved_phase_data",
                message: format!(
                    "missing data payload for phase {:?}",
                    phase.id().as_option()
                ),
            })?;
        let prepared = prepare_phase_thermochemistry(payload)?;
        all_elements.extend(
            prepared
                .element_compositions
                .values()
                .flat_map(|composition| composition.keys().cloned()),
        );
        phase_data.insert(phase.id().clone(), prepared);
    }

    if all_elements.is_empty() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "element_composition",
            message: "resolved phase system contains no elemental composition".to_string(),
        });
    }
    let element_labels = all_elements.into_iter().collect::<Vec<_>>();
    let component_count = metadata.components().len();
    let mut element_composition = DMatrix::zeros(component_count, element_labels.len());
    let mut gibbs = Vec::with_capacity(component_count);
    let mut symbolic_standard_gibbs = Vec::with_capacity(component_count);
    let mut component_reports = Vec::with_capacity(component_count);
    let mut composition_by_substance = HashMap::<String, HashMap<String, f64>>::new();

    for (index, component) in metadata.components().iter().enumerate() {
        let phase = phase_data.get(&component.id().phase).ok_or_else(|| {
            ReactionExtentError::InvalidProblem {
                field: "resolved_phase_data",
                message: format!(
                    "missing prepared phase data for component '{}'",
                    component.label()
                ),
            }
        })?;
        let composition = phase
            .element_compositions
            .get(component.substance())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "element_composition",
                message: format!(
                    "missing elemental composition for component '{}'",
                    component.label()
                ),
            })?;

        if let Some(reference) = composition_by_substance.get(component.substance()) {
            if reference != composition {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "element_composition",
                    message: format!(
                        "component '{}' has a molecular composition inconsistent with another phase record",
                        component.substance()
                    ),
                });
            }
        } else {
            composition_by_substance.insert(component.substance().to_string(), composition.clone());
        }

        for (element_index, element) in element_labels.iter().enumerate() {
            element_composition[(index, element_index)] =
                composition.get(element).copied().unwrap_or(0.0);
        }

        let thermo_source = thermo_source_for(&metadata, component)?;
        let function = phase
            .gibbs_functions
            .get(component.substance())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "standard_gibbs",
                message: format!(
                    "missing standard Gibbs function for component '{}' from {}:{}",
                    component.label(),
                    thermo_source.library(),
                    thermo_source.record_key()
                ),
            })?;
        let standard_gibbs_at_conditions = function(conditions.temperature());
        if !standard_gibbs_at_conditions.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "standard_gibbs",
                message: format!(
                    "component '{}' from {}:{} returned non-finite G0 at {} K",
                    component.label(),
                    thermo_source.library(),
                    thermo_source.record_key(),
                    conditions.temperature()
                ),
            });
        }

        let function = std::sync::Arc::clone(function);
        gibbs.push(Rc::new(move |temperature: f64| function(temperature)) as GibbsFn);
        let symbolic = phase
            .symbolic_gibbs_functions
            .get(component.substance())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "symbolic_standard_gibbs",
                message: format!(
                    "missing symbolic standard Gibbs expression for component '{}' from {}:{}",
                    component.label(),
                    thermo_source.library(),
                    thermo_source.record_key()
                ),
            })?;
        symbolic_standard_gibbs.push(symbolic.clone());
        component_reports.push(EquilibriumBridgeComponentReport {
            component: component.clone(),
            initial_moles: request.initial_composition.moles()[index],
            standard_gibbs_at_conditions,
            thermo_source,
        });
    }

    let phases = metadata
        .phases()
        .iter()
        .map(|phase| Phase {
            kind: phase.activity_model(),
            species: phase.component_range().collect(),
        })
        .collect::<Vec<_>>();
    let initial_moles = request.initial_composition.moles().to_vec();
    let initial_log_moles =
        LogMolesInitialGuess::from_moles_with_policy(&initial_moles, request.trace_seed_policy)?;
    let element_totals = element_totals(&initial_moles, &element_composition);
    let problem = EquilibriumProblem::new(
        metadata.components().to_vec(),
        initial_moles,
        initial_log_moles,
        element_composition,
        gibbs,
        phases,
        conditions,
    )?;
    let report = PhaseEquilibriumBuildReport {
        conditions,
        layout_fingerprint: metadata.layout_fingerprint(),
        element_labels,
        element_totals,
        components: component_reports,
    };

    Ok(PhaseEquilibriumProblemBundle {
        problem,
        metadata,
        report,
        symbolic_standard_gibbs,
    })
}

#[derive(Clone)]
struct PreparedPhaseThermochemistry {
    gibbs_functions: HashMap<String, std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>>,
    symbolic_gibbs_functions: HashMap<String, Expr>,
    element_compositions: HashMap<String, HashMap<String, f64>>,
}

fn prepare_phase_thermochemistry(
    resolved_payload: &SubsData,
) -> Result<PreparedPhaseThermochemistry, ReactionExtentError> {
    let mut working = resolved_payload.clone();
    let gibbs_functions = working.calculate_dG0_fun_one_phase()?;
    let symbolic_gibbs_functions = working.calculate_dG0_sym_one_phase()?;
    let (_, compositions, _) =
        SubsData::calculate_elem_composition_and_molar_mass_local(&mut working, None)?;
    if compositions.len() != working.substances().len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "phase thermochemistry produced {} compositions for {} substances",
            compositions.len(),
            working.substances().len()
        )));
    }

    let mut element_compositions = HashMap::with_capacity(compositions.len());
    for (substance, composition) in working.substances().iter().cloned().zip(compositions) {
        if element_compositions
            .insert(substance.clone(), composition)
            .is_some()
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_thermochemistry",
                message: format!("duplicate local substance '{substance}'"),
            });
        }
    }

    Ok(PreparedPhaseThermochemistry {
        gibbs_functions: gibbs_functions
            .into_iter()
            .map(|(substance, function)| {
                (
                    substance,
                    std::sync::Arc::from(function)
                        as std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>,
                )
            })
            .collect(),
        symbolic_gibbs_functions,
        element_compositions,
    })
}

fn thermo_source_for(
    metadata: &PhaseEquilibriumMetadata,
    component: &EquilibriumComponentDescriptor,
) -> Result<SearchSummaryRow, ReactionExtentError> {
    metadata
        .provenance()
        .phase(&component.id().phase)
        .and_then(|phase| {
            phase.search().rows().iter().find(|row| {
                row.substance() == component.substance()
                    && row.property() == "Thermo"
                    && row.state() == "Found"
            })
        })
        .cloned()
        .ok_or_else(|| ReactionExtentError::InvalidProblem {
            field: "thermochemical_provenance",
            message: format!(
                "component '{}' has no resolved thermochemical provenance row",
                component.label()
            ),
        })
}

fn element_totals(initial_moles: &[f64], element_composition: &DMatrix<f64>) -> Vec<f64> {
    (0..element_composition.ncols())
        .map(|element| {
            initial_moles
                .iter()
                .enumerate()
                .map(|(component, amount)| amount * element_composition[(component, element)])
                .sum()
        })
        .collect()
}

fn activity_model_for(model: PhaseModel) -> PhaseActivityModel {
    match model {
        PhaseModel::IdealGas => PhaseActivityModel::IdealGas,
        PhaseModel::PureCondensed => PhaseActivityModel::IdealSolution,
    }
}
