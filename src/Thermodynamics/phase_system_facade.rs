//! Public multi-phase facade over the crate-private `PhaseSystem` engine.
//!
//! This module owns the ergonomic API for constructing, inspecting, and
//! evaluating a multi-phase system. The engine itself remains in the parent
//! module so cache ownership and derived-state invalidation stay centralized.

use std::collections::HashMap;
use std::fmt;

use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::SubsDataResult;

use super::{
    MoleNumberSnapshot, NumericThermoCacheContext, PhaseComposition, PhaseEvaluationRequest,
    PhaseFunction, PhaseSpec, PhaseSystem, PhaseThermoPropertyView, ResolvedPhaseSystem,
    ResolvedPhaseSystemReport, SubstancePhaseMapping, SystemLayout, ThermoEvaluationConditions,
    ThermoResultSnapshot, ThermoStateSnapshot, build_multi_cache_snapshot,
    build_thermo_result_snapshot, build_thermo_state_snapshot,
};

#[cfg(test)]
use super::PhaseSymbolicLayout;

/// Multi-phase thermodynamic system managing multiple phases or solutions.
///
/// Each phase is identified by an `Option<String>` key where `Some(name)` represents
/// named phases such as `gas`, `liquid`, or `solid`.
///
/// # Example
/// ```rust
/// use KiThe::Thermodynamics::User_PhaseOrSolution::PhaseOrSolution;
///
/// let mut system = PhaseOrSolution::new();
/// // Typed builders resolve phase payloads into the shared `PhaseSystem`.
/// // Read-only inspection uses `system.phase_data_view()`.
/// ```
///
/// The engine is intentionally private to external callers: raw mutation would
/// bypass layout revisioning and cache invalidation. Use a typed builder or a
/// documented transition method instead.
///
/// ```compile_fail
/// use KiThe::Thermodynamics::User_PhaseOrSolution::PhaseOrSolution;
///
/// let system = PhaseOrSolution::new();
/// let _engine = system.phase_system;
/// ```
#[derive(Clone)]
pub struct PhaseOrSolution {
    pub(crate) substance_phase_mapping: Option<SubstancePhaseMapping>,
    pub(crate) phase_system: PhaseSystem,
}

impl PhaseOrSolution {
    /// Creates a new empty multi-phase system.
    pub fn new() -> Self {
        Self {
            substance_phase_mapping: None,
            phase_system: PhaseSystem::new(),
        }
    }

    /// Returns the current structural revision of the phase system.
    pub fn layout_revision(&self) -> usize {
        self.phase_system.layout_revision()
    }

    /// Canonical resolved phase declarations in solver order.
    pub fn phase_specs(&self) -> &[PhaseSpec] {
        self.phase_system.phase_specs()
    }

    /// Returns the canonical layout when this instance came from a typed
    /// `ResolvedPhaseSystem` rather than manual legacy assembly.
    pub fn resolved_layout(&self) -> Option<&SystemLayout> {
        self.phase_system.resolved_layout()
    }

    /// Read-only payload access for diagnostics and adapters. Mutating raw
    /// phase data is restricted to crate-private transition helpers so the
    /// engine can invalidate caches before a component layout changes.
    pub fn phase_data_view(&self) -> &HashMap<Option<String>, SubsData> {
        self.phase_system.phase_data()
    }

    #[inline]
    pub(crate) fn phase_data(&self) -> &HashMap<Option<String>, SubsData> {
        self.phase_system.phase_data()
    }

    /// Replaces transitional raw phase data only in focused tests. Production
    /// construction consumes a fully validated `ResolvedPhaseSystem` through
    /// `install_resolved_system`.
    #[cfg(test)]
    pub(crate) fn replace_phase_data(&mut self, phase_data: HashMap<Option<String>, SubsData>) {
        self.phase_system.replace_phase_data(phase_data);
    }

    /// Atomically installs the payloads, canonical layout, and provenance
    /// produced by typed phase resolution.
    pub(crate) fn install_resolved_system(&mut self, resolved: ResolvedPhaseSystem) {
        self.phase_system.install_resolved_system(resolved);
    }

    /// Executes one transitional raw-data edit, then invalidates from the
    /// final phase-key set. A closure is required so adding or removing a
    /// phase cannot leave cache bundles keyed to the pre-mutation map.
    #[cfg(test)]
    pub(crate) fn with_phase_data_mut<R>(
        &mut self,
        update: impl FnOnce(&mut HashMap<Option<String>, SubsData>) -> R,
    ) -> R {
        self.phase_system.with_phase_data_mut(update)
    }

    /// Adds or replaces one named phase while invalidating all layout-derived
    /// caches. Typed `PhaseSpec` resolution remains the public construction
    /// API; this helper supports crate-level adapters and focused fixtures.
    #[cfg(test)]
    pub(crate) fn insert_phase_data(&mut self, phase_name: String, subs_data: SubsData) {
        self.with_phase_data_mut(|phase_data| {
            phase_data.insert(Some(phase_name), subs_data);
        });
    }

    pub(crate) fn current_layout(&self) -> SystemLayout {
        self.phase_system.current_layout()
    }

    /// Shared symbolic bookkeeping owned by `PhaseSystem`. This crate-private
    /// view exists for adapters and tests during the facade migration.
    #[cfg(test)]
    pub(crate) fn symbolic_layout_view(&self) -> &PhaseSymbolicLayout {
        self.phase_system.symbolic_layout()
    }

    /// Evaluates Gibbs energies without publishing a `PhaseSystem` cache.
    /// `SubsData` may refresh its own coefficient cache, but this method leaves
    /// the phase-system result snapshot unchanged until a caller explicitly
    /// chooses to publish the returned values.
    pub fn evaluate_gibbs(
        &self,
        temperature: f64,
        pressure: f64,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.evaluate_gibbs_at(
            ThermoEvaluationConditions::new(temperature, pressure)?,
            compositions,
        )
    }

    /// Evaluates Gibbs energy from a validated request without mutating this
    /// system or publishing a cache. Lower-level calculators receive local
    /// snapshots because their coefficient preparation remains stateful.
    pub fn evaluate_gibbs_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system
            .evaluate_gibbs_at(conditions, compositions)
    }

    /// Evaluates Gibbs energy from one self-contained request. This is the
    /// preferred pure-evaluation entry point for new phase-system consumers.
    pub fn evaluate_gibbs_request(
        &self,
        request: &PhaseEvaluationRequest,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system.evaluate_gibbs_request(request)
    }

    /// Evaluates entropies without publishing a `PhaseSystem` cache.
    pub fn evaluate_entropy(
        &self,
        temperature: f64,
        pressure: f64,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.evaluate_entropy_at(
            ThermoEvaluationConditions::new(temperature, pressure)?,
            compositions,
        )
    }

    /// Evaluates entropy from a validated request without publishing a cache.
    pub fn evaluate_entropy_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system
            .evaluate_entropy_at(conditions, compositions)
    }

    /// Entropy counterpart of `evaluate_gibbs_request`.
    pub fn evaluate_entropy_request(
        &self,
        request: &PhaseEvaluationRequest,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system.evaluate_entropy_request(request)
    }

    /// Borrowed Gibbs energy cache for read-only consumers.
    pub fn dG_view(&self) -> PhaseThermoPropertyView<'_, f64> {
        self.phase_system.dG_view()
    }

    /// Borrowed Gibbs function cache for read-only consumers.
    pub fn dG_fun_view(&self) -> PhaseThermoPropertyView<'_, PhaseFunction> {
        self.phase_system.dG_fun_view()
    }

    /// Borrowed symbolic Gibbs cache for read-only consumers.
    pub fn dG_sym_view(&self) -> PhaseThermoPropertyView<'_, Expr> {
        self.phase_system.dG_sym_view()
    }

    /// Borrowed entropy cache for read-only consumers.
    pub fn dS_view(&self) -> PhaseThermoPropertyView<'_, f64> {
        self.phase_system.dS_view()
    }

    /// Borrowed entropy function cache for read-only consumers.
    pub fn dS_fun_view(&self) -> PhaseThermoPropertyView<'_, PhaseFunction> {
        self.phase_system.dS_fun_view()
    }

    /// Borrowed symbolic entropy cache for read-only consumers.
    pub fn dS_sym_view(&self) -> PhaseThermoPropertyView<'_, Expr> {
        self.phase_system.dS_sym_view()
    }

    /// Typed lookup provenance for the resolved payloads, when this system
    /// came from the typed resolution pipeline.
    pub fn resolution_report(&self) -> Option<&ResolvedPhaseSystemReport> {
        self.phase_system.resolution_report()
    }

    /// Returns the request that validates the published numeric Gibbs cache.
    pub fn gibbs_cache_context(&self) -> Option<&NumericThermoCacheContext> {
        self.phase_system.gibbs_cache_context()
    }

    /// Returns whether the published Gibbs result is valid for this exact request.
    pub fn has_current_gibbs_cache(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> bool {
        self.phase_system
            .has_current_gibbs_cache(conditions, compositions)
    }

    pub fn has_current_gibbs_request(&self, request: &PhaseEvaluationRequest) -> bool {
        self.phase_system.has_current_gibbs_request(request)
    }

    /// Returns the request that validates the published numeric entropy cache.
    pub fn entropy_cache_context(&self) -> Option<&NumericThermoCacheContext> {
        self.phase_system.entropy_cache_context()
    }

    /// Returns whether the published entropy result is valid for this exact request.
    pub fn has_current_entropy_cache(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> bool {
        self.phase_system
            .has_current_entropy_cache(conditions, compositions)
    }

    pub fn has_current_entropy_request(&self, request: &PhaseEvaluationRequest) -> bool {
        self.phase_system.has_current_entropy_request(request)
    }

    #[cfg(test)]
    pub(crate) fn debug_replace_dG(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, f64>>,
    ) {
        self.phase_system.replace_dG(values);
    }

    #[cfg(test)]
    pub(crate) fn debug_replace_dG_fun(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, PhaseFunction>>,
    ) {
        self.phase_system.replace_dG_fun(values);
    }

    #[cfg(test)]
    pub(crate) fn debug_replace_dG_sym(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, Expr>>,
    ) {
        self.phase_system.replace_dG_sym(values);
    }

    #[cfg(test)]
    pub(crate) fn debug_replace_dS(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, f64>>,
    ) {
        self.phase_system.replace_dS(values);
    }

    #[cfg(test)]
    pub(crate) fn debug_replace_dS_fun(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, PhaseFunction>>,
    ) {
        self.phase_system.replace_dS_fun(values);
    }

    /// Returns a typed snapshot of the current state and caches.
    pub fn state_snapshot(&self) -> ThermoStateSnapshot<'_> {
        let layout_revision = self.phase_system.layout_revision();
        build_thermo_state_snapshot(
            layout_revision,
            build_multi_cache_snapshot(layout_revision, self.phase_system.dG_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dG_fun_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dG_sym_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dS_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dS_fun_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dS_sym_view()),
        )
    }

    /// Returns a typed thermodynamic result snapshot aligned to the canonical phase layout.
    pub fn result_snapshot(
        &self,
        temperature: Option<f64>,
        pressure: Option<f64>,
    ) -> ThermoResultSnapshot<'_> {
        let layout_revision = self.phase_system.layout_revision();
        build_thermo_result_snapshot(
            layout_revision,
            self.current_layout(),
            temperature,
            pressure,
            build_multi_cache_snapshot(layout_revision, self.phase_system.dG_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dG_fun_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dG_sym_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dS_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dS_fun_view()),
            build_multi_cache_snapshot(layout_revision, self.phase_system.dS_sym_view()),
        )
    }

    /// Normalizes sparse phase mole input into a snapshot that carries the
    /// exact `SystemLayout` used for all ordered component vectors.
    pub fn normalize_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
    ) -> SubsDataResult<MoleNumberSnapshot> {
        self.phase_system
            .normalize_mole_numbers(non_zero_number_of_moles)
    }

    fn bump_layout_revision(&mut self) {
        self.phase_system.bump_layout_revision();
    }

    /// Sets the substance-phase mapping for the system.
    /// Used to track which substances belong to which phases.
    pub fn set_substance_phase_mapping(&mut self, mapping: SubstancePhaseMapping) {
        self.substance_phase_mapping = Some(mapping);
        self.bump_layout_revision();
    }

    pub fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
        >,
        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    > {
        self.phase_system
            .build_numeric_lagrange_equations(A, G_fun, Tm)
    }
}

impl fmt::Debug for PhaseOrSolution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PhaseOrSolution")
            .field("phase_data", self.phase_data())
            .field("substance_phase_mapping", &self.substance_phase_mapping)
            .field("symbolic_layout", self.phase_system.symbolic_layout())
            .field("dG", &self.phase_system.dG_view())
            .field("dG_sym", &self.phase_system.dG_sym_view())
            .field("dS", &self.phase_system.dS_view())
            .field("dS_sym", &self.phase_system.dS_sym_view())
            .finish()
    }
}
