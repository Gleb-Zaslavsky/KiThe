//! Shared phase-engine core for single- and multi-phase thermodynamics.
//!
//! ```text
//! HIERARCHY
//!     
//!                                                      PhaseSystem
//!                           /                                     |                                    \                                  
//!                     phase_bundles                           phase_specs                        symbolic_layout                   resolved_layout
//!       
//!                /             \                                     |                                    |                                    |       
//!           GIBBS             ENTROPY                                 PhaseSpec                      PhaseSymbolicLayout                    SystemLayout
//!             |                  |                                 /     |           \                      |                             /      |         \
//!        PhasePropertyCache      PhasePropertyCache        PhaseId PhasePhysicalState PhaseModel     SystemLayout                   phases   components  phase_ranges
//!        /      |       \          /          |    \                                             /      |         \                    /         |             \
//!   NUMERIC  FUNCTION  SYMBOLIC   NUMERIC  FUNCTION  SYMBOLIC                              phases   components  phase_ranges              PhaseComponentId
//!               |                             |
//!        PhaseFunction                   PhaseFunction    
//! ```

//! This module owns the mutable payload map, canonical layout metadata,
//! transactional updates, and cached Gibbs/entropy result families. Public
//! facades remain thin and delegate orchestration here.

use std::collections::HashMap;
use std::f64;

use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;

use super::{
    IndexedMoleVariables, LegacyMoleNumberSnapshot, MoleNumberSnapshot, NumericPhaseProperty,
    NumericThermoCacheContext, PhaseComposition, PhaseDataView, PhaseEvaluationRequest,
    PhaseFunction, PhaseLagrangeFunction, PhaseSpec, PhaseSymbolicLayout, PhaseThermoCacheBundle,
    PhaseThermoPropertyView, ResolvedPhaseSystem, ResolvedPhaseSystemReport,
    ThermoEvaluationConditions, build_indexed_mole_variables, element_composition_and_molar_mass,
    evaluate_numeric_phase_property, normalize_mole_numbers, stage_entropy_functions,
    stage_gibbs_functions, stage_phase_data_operation, stage_symbolic_entropy,
    stage_symbolic_gibbs, substitute_exprs, validate_phase_evaluation_request,
};

use super::phase_cache::boxify_fun_map;

#[cfg(test)]
use super::phase_cache::{
    ThermoStateSnapshot, build_multi_cache_snapshot, build_thermo_state_snapshot,
};

/// Shared cache/state core for phase thermodynamics wrappers.
///
/// This is crate-private deliberately: public callers work through typed phase
/// specifications and facade views, rather than publishing inconsistent cache
/// state directly.
#[derive(Clone)]
pub(crate) struct PhaseSystem {
    phase_data: HashMap<Option<String>, SubsData>,
    phase_bundles: HashMap<Option<String>, PhaseThermoCacheBundle>,
    symbolic_layout: PhaseSymbolicLayout,
    phase_specs: Vec<PhaseSpec>,
    resolved_layout: Option<SystemLayout>,
    resolution_report: Option<ResolvedPhaseSystemReport>,
    layout_revision: usize,
    configuration_revision: usize,
    gibbs_cache_context: Option<NumericThermoCacheContext>,
    entropy_cache_context: Option<NumericThermoCacheContext>,
}

impl PhaseSystem {
    pub fn new() -> Self {
        Self {
            phase_data: HashMap::new(),
            phase_bundles: HashMap::new(),
            symbolic_layout: PhaseSymbolicLayout::default(),
            phase_specs: Vec::new(),
            resolved_layout: None,
            resolution_report: None,
            layout_revision: 0,
            configuration_revision: 0,
            gibbs_cache_context: None,
            entropy_cache_context: None,
        }
    }

    pub fn new_single_phase() -> Self {
        let mut system = Self::new();
        system.phase_data.insert(None, SubsData::new());
        system
            .phase_bundles
            .insert(None, PhaseThermoCacheBundle::new());
        system
    }

    #[inline]
    pub fn layout_revision(&self) -> usize {
        self.layout_revision
    }

    #[inline]
    pub fn bump_layout_revision(&mut self) {
        self.layout_revision = self.layout_revision.wrapping_add(1);
        self.configuration_revision = self.configuration_revision.wrapping_add(1);
        // A new layout/configuration revision makes every cached representation
        // stale. Rebuild bundles from canonical payload keys so borrowed views
        // cannot expose values calculated for the preceding state.
        self.phase_bundles = self
            .phase_data
            .keys()
            .cloned()
            .map(|phase| (phase, PhaseThermoCacheBundle::new()))
            .collect();
        self.gibbs_cache_context = None;
        self.entropy_cache_context = None;
    }

    /// Canonical declarations installed by the typed resolution path.
    pub fn phase_specs(&self) -> &[PhaseSpec] {
        &self.phase_specs
    }

    /// Ordered component layout installed by the typed resolution path.
    /// Manually assembled legacy wrappers intentionally return `None` until
    /// they are migrated through an explicit resolved specification.
    pub fn resolved_layout(&self) -> Option<&SystemLayout> {
        self.resolved_layout.as_ref()
    }

    /// Returns the resolved layout when available, otherwise derives the
    /// deterministic compatibility layout from the canonical payload map.
    pub(crate) fn current_layout(&self) -> SystemLayout {
        self.resolved_layout
            .clone()
            .unwrap_or_else(|| PhaseDataView::multi(&self.phase_data).layout())
    }

    /// Canonical raw payloads aligned to the phase engine. The map remains a
    /// compatibility/data-acquisition shape; solver-facing order is supplied
    /// exclusively by `SystemLayout`.
    pub(crate) fn phase_data(&self) -> &HashMap<Option<String>, SubsData> {
        &self.phase_data
    }

    pub fn resolution_report(&self) -> Option<&ResolvedPhaseSystemReport> {
        self.resolution_report.as_ref()
    }

    /// Returns the invariant-preserving one-phase payload. `new_single_phase`,
    /// replacement, and resolution all install the `None` entry before a
    /// `OnePhase` facade is constructed.
    pub(crate) fn single_phase_data(&self) -> &SubsData {
        self.phase_data
            .get(&None)
            .expect("OnePhase requires a canonical None phase payload")
    }

    #[cfg(test)]
    pub(crate) fn single_phase_data_mut(&mut self) -> &mut SubsData {
        self.phase_data
            .get_mut(&None)
            .expect("OnePhase requires a canonical None phase payload")
    }

    /// Evaluates Gibbs energy from cloned phase payloads without publishing a
    /// cache. Both one- and multi-phase facades use this one engine path.
    pub(crate) fn evaluate_gibbs_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        let request = PhaseEvaluationRequest::new(conditions, compositions.clone());
        self.evaluate_gibbs_request(&request)
    }

    /// Evaluates Gibbs energy from one atomic request without publishing a
    /// cache. New callers should prefer this over parallel input arguments.
    pub(crate) fn evaluate_gibbs_request(
        &self,
        request: &PhaseEvaluationRequest,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        evaluate_numeric_phase_property(
            PhaseDataView::multi(&self.phase_data),
            request,
            NumericPhaseProperty::Gibbs,
        )
    }

    /// Evaluates entropy from cloned phase payloads without publishing a
    /// cache. Keeping this beside Gibbs prevents one-phase behaviour drifting
    /// from the multi-phase calculation contract.
    pub(crate) fn evaluate_entropy_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        let request = PhaseEvaluationRequest::new(conditions, compositions.clone());
        self.evaluate_entropy_request(&request)
    }

    /// Entropy counterpart of `evaluate_gibbs_request`.
    pub(crate) fn evaluate_entropy_request(
        &self,
        request: &PhaseEvaluationRequest,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        evaluate_numeric_phase_property(
            PhaseDataView::multi(&self.phase_data),
            request,
            NumericPhaseProperty::Entropy,
        )
    }

    /// Builds the solver-facing element matrix from cloned payloads. The
    /// result is pure: asking for element data must not quietly mutate one
    /// facade while leaving the other facade unchanged.
    pub(crate) fn element_composition_and_molar_mass(
        &self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)> {
        element_composition_and_molar_mass(&self.phase_data, &self.current_layout(), groups)
    }

    /// Generates all symbolic mole variables from one ordered phase layout.
    /// A one-phase system deliberately receives the same `Np0`/`n0_i` naming
    /// scheme as a multi-phase system, so expressions remain composable when
    /// a single phase later becomes part of a larger model.
    pub(crate) fn indexed_moles_variables(&mut self) -> SubsDataResult<IndexedMoleVariables> {
        let layout = self.current_layout();
        self.replace_symbolic_layout(build_indexed_mole_variables(&layout)?);
        self.bump_layout_revision();
        Ok(self.symbolic_layout.legacy_projection())
    }

    /// Normalizes sparse mole input without exposing `HashMap` iteration order
    /// to solver vectors. Absent phases retain the legacy meaning of "not
    /// supplied"; supplied phases receive every declared component with zero
    /// inserted where necessary.
    pub(crate) fn normalize_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
    ) -> SubsDataResult<MoleNumberSnapshot> {
        normalize_mole_numbers(self.current_layout(), non_zero_number_of_moles)
    }

    /// Compatibility projection for callers that have not yet migrated from
    /// the tuple-of-maps representation.
    pub(crate) fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
    ) -> SubsDataResult<LegacyMoleNumberSnapshot> {
        self.normalize_mole_numbers(non_zero_number_of_moles)
            .map(MoleNumberSnapshot::into_legacy)
    }

    /// Delegates symbolic stationarity assembly to the pure equilibrium
    /// adapter while keeping layout/cache ownership inside the phase engine.
    pub(crate) fn build_symbolic_lagrange_equations(
        &self,
        element_matrix: DMatrix<f64>,
        reference_temperature: f64,
    ) -> SubsDataResult<Vec<Expr>> {
        crate::Thermodynamics::phase_equilibrium_adapter::build_symbolic_lagrange_equations(
            &self.current_layout(),
            self.dG_sym_view(),
            &element_matrix,
            reference_temperature,
        )
    }

    /// Builds numeric stationarity equations from an explicitly supplied
    /// legacy closure map. Closure ownership stays at this narrow adapter
    /// boundary until equilibrium callers migrate to typed compositions.
    pub(crate) fn build_numeric_lagrange_equations(
        &self,
        element_matrix: DMatrix<f64>,
        gibbs_functions: crate::Thermodynamics::phase_equilibrium_adapter::LegacyPhaseGibbsFunctions,
        reference_temperature: f64,
    ) -> SubsDataResult<PhaseLagrangeFunction> {
        crate::Thermodynamics::phase_equilibrium_adapter::build_numeric_lagrange_equations(
            &self.current_layout(),
            gibbs_functions,
            element_matrix,
            reference_temperature,
        )
    }

    /// Rebuilds Gibbs closures and then constructs numeric stationarity
    /// equations through the same phase-engine path for one or many phases.
    pub(crate) fn calculate_lagrange_equations_fun2(
        &mut self,
        element_matrix: DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_temperature: f64,
    ) -> SubsDataResult<PhaseLagrangeFunction> {
        let (staged_phase_data, phase_functions) =
            stage_gibbs_functions(&self.phase_data, temperature, pressure)?;
        let equations = self.build_numeric_lagrange_equations(
            element_matrix,
            phase_functions
                .iter()
                .map(|(phase, functions)| (phase.clone(), boxify_fun_map(functions)))
                .collect(),
            reference_temperature,
        )?;

        // Equation construction can still reject a matrix or reference
        // temperature. Publish both the payload refresh and the function cache
        // only after that fallible assembly step has succeeded.
        self.validate_phase_property_publication(&phase_functions, "Gibbs closure publication")?;
        self.phase_data = staged_phase_data;
        self.commit_phase_property(phase_functions, |bundle, property| {
            bundle.gibbs.functions = property
        });
        Ok(equations)
    }

    fn try_apply_to_phase_data<R>(
        &mut self,
        update: impl FnMut(&Option<String>, &mut SubsData) -> SubsDataResult<R>,
    ) -> SubsDataResult<HashMap<Option<String>, R>> {
        let (updated_phase_data, results) = stage_phase_data_operation(&self.phase_data, update)?;
        self.phase_data = updated_phase_data;
        Ok(results)
    }

    #[cfg(test)]
    pub(crate) fn debug_try_apply_to_phase_data(
        &mut self,
        update: impl FnMut(&Option<String>, &mut SubsData) -> SubsDataResult<()>,
    ) -> SubsDataResult<()> {
        self.try_apply_to_phase_data(update).map(|_| ())
    }

    /// Prepares all phase-local thermal coefficients for a temperature.
    pub(crate) fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.try_apply_to_phase_data(|_, subs_data| {
            subs_data.extract_all_thermal_coeffs(temperature)
        })?;
        self.bump_layout_revision();
        Ok(())
    }

    /// Rebuilds numeric thermodynamic property maps for every phase.
    pub(crate) fn calculate_therm_map_of_properties(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<()> {
        self.try_apply_to_phase_data(|_, subs_data| {
            subs_data.calculate_therm_map_of_properties(temperature)
        })?;
        self.bump_layout_revision();
        Ok(())
    }

    /// Rebuilds symbolic thermodynamic property maps for every phase.
    pub(crate) fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()> {
        self.try_apply_to_phase_data(|_, subs_data| subs_data.calculate_therm_map_of_sym())?;
        self.bump_layout_revision();
        Ok(())
    }

    /// Ensures every phase has valid coefficients for the requested temperature.
    pub(crate) fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        let results = self.try_apply_to_phase_data(|_, subs_data| {
            subs_data.extract_coeffs_if_current_coeffs_not_valid_for_all_subs(temperature)
        })?;
        self.bump_layout_revision();
        Ok(results.into_values().flatten().collect())
    }

    /// Applies shared pressure and molar-mass configuration to all phases.
    pub(crate) fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        self.try_apply_to_phase_data(|_, subs_data| {
            subs_data.set_P(pressure, pressure_unit.clone())?;
            subs_data.set_M(molar_masses.clone(), mass_unit.clone())?;
            Ok(())
        })?;
        self.bump_layout_revision();
        Ok(())
    }

    /// Runs the phase-local NIST fallback policy across every phase.
    pub(crate) fn fetch_missing_from_nist(&mut self) -> SubsDataResult<()> {
        self.try_apply_to_phase_data(|_, subs_data| subs_data.if_not_found_go_NIST())?;
        self.bump_layout_revision();
        Ok(())
    }

    /// Computes and publishes numeric Gibbs values with their exact request
    /// context. Both facades therefore share one cache-validity contract.
    pub(crate) fn calculate_gibbs_free_energy(
        &mut self,
        temperature: f64,
        pressure: f64,
        compositions: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.calculate_gibbs_free_energy_request(PhaseEvaluationRequest::from_numeric(
            temperature,
            pressure,
            compositions,
        )?)
    }

    pub(crate) fn calculate_gibbs_free_energy_request(
        &mut self,
        request: PhaseEvaluationRequest,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        let values = self.evaluate_gibbs_request(&request)?;
        self.publish_gibbs_evaluation_request(request, values.clone())?;
        Ok(values)
    }

    /// Computes and publishes numeric entropy values with their exact request
    /// context.
    pub(crate) fn calculate_entropy(
        &mut self,
        temperature: f64,
        pressure: f64,
        compositions: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<()> {
        self.calculate_entropy_request(PhaseEvaluationRequest::from_numeric(
            temperature,
            pressure,
            compositions,
        )?)
    }

    pub(crate) fn calculate_entropy_request(
        &mut self,
        request: PhaseEvaluationRequest,
    ) -> SubsDataResult<()> {
        let values = self.evaluate_entropy_request(&request)?;
        self.publish_entropy_evaluation_request(request, values)?;
        Ok(())
    }

    /// Builds symbolic Gibbs expressions against the canonical phase layout.
    /// Missing symbolic mole variables mean a standard-state expression and are
    /// represented uniformly as `(None, None)` for every phase.
    pub(crate) fn calculate_gibbs_sym(&mut self, temperature: f64) -> SubsDataResult<()> {
        let (phase_data, values) =
            stage_symbolic_gibbs(&self.phase_data, &self.symbolic_layout, temperature)?;
        self.validate_phase_property_publication(&values, "symbolic Gibbs publication")?;
        self.phase_data = phase_data;
        self.commit_phase_property(values, |bundle, property| bundle.gibbs.symbolic = property);
        Ok(())
    }

    /// Substitutes pressure in the already published symbolic Gibbs cache.
    pub(crate) fn set_pressure_in_gibbs_sym(&mut self, pressure: f64) {
        for bundle in self.phase_bundles.values_mut() {
            substitute_exprs(bundle.gibbs.symbolic.values_mut(), "P", pressure);
        }
        self.bump_cache_revision();
    }

    /// Substitutes temperature in the already published symbolic Gibbs cache.
    pub(crate) fn set_temperature_in_gibbs_sym(&mut self, temperature: f64) {
        for bundle in self.phase_bundles.values_mut() {
            substitute_exprs(bundle.gibbs.symbolic.values_mut(), "T", temperature);
        }
        self.bump_cache_revision();
    }

    /// Builds and publishes Gibbs closures for all canonical phase payloads.
    pub(crate) fn calculate_gibbs_fun(
        &mut self,
        temperature: f64,
        pressure: f64,
    ) -> SubsDataResult<()> {
        let (phase_data, values) = stage_gibbs_functions(&self.phase_data, temperature, pressure)?;
        self.validate_phase_property_publication(&values, "Gibbs closure publication")?;
        self.phase_data = phase_data;
        self.commit_phase_property(values, |bundle, property| bundle.gibbs.functions = property);
        Ok(())
    }

    /// Builds symbolic entropy expressions against the canonical phase layout.
    pub(crate) fn calculate_entropy_sym(&mut self, temperature: f64) -> SubsDataResult<()> {
        let (phase_data, values) =
            stage_symbolic_entropy(&self.phase_data, &self.symbolic_layout, temperature)?;
        self.validate_phase_property_publication(&values, "symbolic entropy publication")?;
        self.phase_data = phase_data;
        self.commit_phase_property(values, |bundle, property| {
            bundle.entropy.symbolic = property
        });
        Ok(())
    }

    /// Builds and publishes entropy closures for all canonical phase payloads.
    pub(crate) fn calculate_entropy_fun(
        &mut self,
        temperature: f64,
        pressure: f64,
    ) -> SubsDataResult<()> {
        let (phase_data, values) =
            stage_entropy_functions(&self.phase_data, temperature, pressure)?;
        self.validate_phase_property_publication(&values, "entropy closure publication")?;
        self.phase_data = phase_data;
        self.commit_phase_property(values, |bundle, property| {
            bundle.entropy.functions = property
        });
        Ok(())
    }

    /// Replaces legacy/raw payloads as one transaction. The post-replacement
    /// phase keys, rather than stale pre-edit keys, determine the empty cache
    /// bundles installed by the invalidation boundary.
    #[cfg(test)]
    pub(crate) fn replace_phase_data(&mut self, phase_data: HashMap<Option<String>, SubsData>) {
        self.phase_data = phase_data;
        self.resolution_report = None;
        self.prepare_for_legacy_data_change();
    }

    #[cfg(test)]
    pub(crate) fn replace_single_phase_data(&mut self, subs_data: SubsData) {
        self.phase_data = HashMap::from([(None, subs_data)]);
        self.resolution_report = None;
        self.prepare_for_legacy_data_change();
    }

    /// Performs an internal raw-data edit and invalidates using the final map.
    /// This is the only mutation escape hatch during the legacy-consumer
    /// migration; it prevents phase insertion/removal from desynchronising the
    /// cache bundle keys.
    #[cfg(test)]
    pub(crate) fn with_phase_data_mut<R>(
        &mut self,
        update: impl FnOnce(&mut HashMap<Option<String>, SubsData>) -> R,
    ) -> R {
        let result = update(&mut self.phase_data);
        self.prepare_for_legacy_data_change();
        result
    }

    #[cfg(test)]
    pub(crate) fn with_single_phase_data_mut<R>(
        &mut self,
        update: impl FnOnce(&mut SubsData) -> R,
    ) -> R {
        let result = update(self.single_phase_data_mut());
        self.prepare_for_legacy_data_change();
        result
    }

    /// Installs the immutable metadata shared by one- and multi-phase wrappers.
    /// Existing computed caches are cleared because they cannot outlive a new
    /// resolved component layout.
    #[cfg(test)]
    pub(crate) fn configure_resolved_phases(
        &mut self,
        phase_specs: Vec<PhaseSpec>,
        layout: SystemLayout,
    ) {
        self.rebuild_phase_bundles(layout.phases().iter().map(|phase| phase.as_option()));
        self.phase_specs = phase_specs;
        self.resolved_layout = Some(layout);
        self.resolution_report = None;
        self.reset_symbolic_layout();
        self.bump_revision_and_clear_contexts();
    }

    /// Installs one already-validated resolved system as a single state
    /// transition. Payloads, layout metadata, cache bundles, and lookup
    /// provenance become visible together, so no facade can observe an
    /// intermediate payload-only configuration.
    pub(crate) fn install_resolved_system(&mut self, resolved: ResolvedPhaseSystem) {
        let (phase_specs, phase_data, layout, report) = resolved.into_parts();

        self.phase_data = phase_data;
        self.rebuild_phase_bundles(layout.phases().iter().map(|phase| phase.as_option()));
        self.phase_specs = phase_specs;
        self.resolved_layout = Some(layout);
        self.resolution_report = Some(report);
        self.reset_symbolic_layout();
        self.bump_layout_revision();
    }

    /// Invalidates typed resolution metadata before a transitional raw-data
    /// mutation. This remains crate-private: public callers must use the
    /// specification and resolution pipeline instead of editing payloads.
    #[cfg(test)]
    pub(crate) fn prepare_for_legacy_data_change(&mut self) {
        let phase_keys: Vec<Option<String>> = self.phase_data.keys().cloned().collect();
        self.rebuild_phase_bundles(phase_keys.iter());
        self.clear_resolved_layout();
        self.reset_symbolic_layout();
        self.bump_layout_revision();
    }

    /// Rebuilds the cache bundle map from the provided canonical phase keys.
    /// The helper keeps bundle publication tied to one final key set instead of
    /// repeating ad hoc reconstruction logic in each mutation path.
    fn rebuild_phase_bundles<'a, I>(&mut self, phases: I)
    where
        I: IntoIterator<Item = &'a Option<String>>,
    {
        self.phase_bundles = phases
            .into_iter()
            .cloned()
            .map(|phase| (phase, PhaseThermoCacheBundle::new()))
            .collect();
    }

    /// Clears typed layout metadata that cannot survive a legacy raw-data edit.
    #[cfg(test)]
    fn clear_resolved_layout(&mut self) {
        self.phase_specs.clear();
        self.resolved_layout = None;
        self.resolution_report = None;
    }

    /// Resets the derived symbolic layout to its empty canonical state.
    fn reset_symbolic_layout(&mut self) {
        self.symbolic_layout = PhaseSymbolicLayout::default();
    }

    pub(crate) fn symbolic_layout(&self) -> &PhaseSymbolicLayout {
        &self.symbolic_layout
    }

    /// Replaces all symbolic variables atomically.  The caller owns the
    /// surrounding revision policy because indexed-variable construction is a
    /// derived-layout operation in both one- and multi-phase facades.
    pub(crate) fn replace_symbolic_layout(&mut self, layout: PhaseSymbolicLayout) {
        self.symbolic_layout = layout;
    }

    fn bump_cache_revision(&mut self) {
        self.layout_revision = self.layout_revision.wrapping_add(1);
    }

    /// Bumps the structural revisions and invalidates cache contexts without
    /// rebuilding bundle keys from the raw payload map. This is used when a
    /// caller has already supplied the canonical bundle key set directly from
    /// a typed layout.
    #[cfg(test)]
    fn bump_revision_and_clear_contexts(&mut self) {
        self.layout_revision = self.layout_revision.wrapping_add(1);
        self.configuration_revision = self.configuration_revision.wrapping_add(1);
        self.gibbs_cache_context = None;
        self.entropy_cache_context = None;
    }

    #[cfg(test)]
    #[inline]
    fn bundle_or_insert(&mut self, phase: Option<String>) -> &mut PhaseThermoCacheBundle {
        self.phase_bundles.entry(phase).or_default()
    }

    pub(crate) fn phase_bundles(&self) -> &HashMap<Option<String>, PhaseThermoCacheBundle> {
        &self.phase_bundles
    }

    pub(crate) fn single_dG_view(&self) -> &HashMap<String, f64> {
        &self.single_bundle().gibbs.numeric
    }

    pub(crate) fn single_dG_fun_view(&self) -> &HashMap<String, PhaseFunction> {
        &self.single_bundle().gibbs.functions
    }

    pub(crate) fn single_dG_sym_view(&self) -> &HashMap<String, Expr> {
        &self.single_bundle().gibbs.symbolic
    }

    pub(crate) fn single_dS_view(&self) -> &HashMap<String, f64> {
        &self.single_bundle().entropy.numeric
    }

    pub(crate) fn single_dS_fun_view(&self) -> &HashMap<String, PhaseFunction> {
        &self.single_bundle().entropy.functions
    }

    pub(crate) fn single_dS_sym_view(&self) -> &HashMap<String, Expr> {
        &self.single_bundle().entropy.symbolic
    }

    fn single_bundle(&self) -> &PhaseThermoCacheBundle {
        self.phase_bundles
            .get(&None)
            .expect("single-phase bundle must exist")
    }

    #[cfg(test)]
    #[inline]
    fn single_bundle_mut(&mut self) -> &mut PhaseThermoCacheBundle {
        self.phase_bundles
            .entry(None)
            .or_insert_with(PhaseThermoCacheBundle::new)
    }

    #[inline]
    pub fn dG_view(&self) -> PhaseThermoPropertyView<'_, f64> {
        PhaseThermoPropertyView::new(self.phase_bundles(), |bundle| &bundle.gibbs.numeric)
    }

    #[inline]
    pub fn dG_fun_view(&self) -> PhaseThermoPropertyView<'_, PhaseFunction> {
        PhaseThermoPropertyView::new(self.phase_bundles(), |bundle| &bundle.gibbs.functions)
    }

    #[inline]
    pub fn dG_sym_view(&self) -> PhaseThermoPropertyView<'_, Expr> {
        PhaseThermoPropertyView::new(self.phase_bundles(), |bundle| &bundle.gibbs.symbolic)
    }

    #[inline]
    pub fn dS_view(&self) -> PhaseThermoPropertyView<'_, f64> {
        PhaseThermoPropertyView::new(self.phase_bundles(), |bundle| &bundle.entropy.numeric)
    }

    #[inline]
    pub fn dS_fun_view(&self) -> PhaseThermoPropertyView<'_, PhaseFunction> {
        PhaseThermoPropertyView::new(self.phase_bundles(), |bundle| &bundle.entropy.functions)
    }

    #[inline]
    pub fn dS_sym_view(&self) -> PhaseThermoPropertyView<'_, Expr> {
        PhaseThermoPropertyView::new(self.phase_bundles(), |bundle| &bundle.entropy.symbolic)
    }

    #[cfg(test)]
    /// Test/debug escape hatch for single-phase numeric Gibbs fixtures.
    /// Production publication must use the validated cache commit paths.
    pub(crate) fn dG_mut(&mut self) -> &mut HashMap<String, f64> {
        &mut self.single_bundle_mut().gibbs.numeric
    }

    #[cfg(test)]
    /// Test/debug escape hatch for single-phase Gibbs closure fixtures.
    pub(crate) fn dG_fun_mut(&mut self) -> &mut HashMap<String, PhaseFunction> {
        &mut self.single_bundle_mut().gibbs.functions
    }

    #[cfg(test)]
    /// Test/debug escape hatch for single-phase symbolic Gibbs fixtures.
    pub(crate) fn dG_sym_mut(&mut self) -> &mut HashMap<String, Expr> {
        &mut self.single_bundle_mut().gibbs.symbolic
    }

    #[cfg(test)]
    /// Test/debug escape hatch for single-phase entropy closure fixtures.
    pub(crate) fn dS_fun_mut(&mut self) -> &mut HashMap<String, PhaseFunction> {
        &mut self.single_bundle_mut().entropy.functions
    }

    #[cfg(test)]
    #[cfg(test)]
    pub(crate) fn replace_dG(&mut self, values: HashMap<Option<String>, HashMap<String, f64>>) {
        self.replace_phase_property_for_test(values, |bundle, property| {
            bundle.gibbs.numeric = property
        });
        self.gibbs_cache_context = None;
        self.bump_cache_revision();
    }

    #[cfg(test)]
    #[cfg(test)]
    pub(crate) fn replace_dG_fun(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, PhaseFunction>>,
    ) {
        self.replace_phase_property_for_test(values, |bundle, property| {
            bundle.gibbs.functions = property
        });
        self.bump_cache_revision();
    }

    #[cfg(test)]
    #[cfg(test)]
    pub(crate) fn replace_dG_sym(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, Expr>>,
    ) {
        self.replace_phase_property_for_test(values, |bundle, property| {
            bundle.gibbs.symbolic = property
        });
        self.bump_cache_revision();
    }

    #[cfg(test)]
    #[cfg(test)]
    pub(crate) fn replace_dS(&mut self, values: HashMap<Option<String>, HashMap<String, f64>>) {
        self.replace_phase_property_for_test(values, |bundle, property| {
            bundle.entropy.numeric = property
        });
        self.entropy_cache_context = None;
        self.bump_cache_revision();
    }

    #[cfg(test)]
    #[cfg(test)]
    pub(crate) fn replace_dS_fun(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, PhaseFunction>>,
    ) {
        self.replace_phase_property_for_test(values, |bundle, property| {
            bundle.entropy.functions = property
        });
        self.bump_cache_revision();
    }

    #[cfg(test)]
    #[cfg(test)]
    pub(crate) fn replace_dS_sym(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, Expr>>,
    ) {
        self.replace_phase_property_for_test(values, |bundle, property| {
            bundle.entropy.symbolic = property
        });
        self.bump_cache_revision();
    }

    /// Publishes a numeric Gibbs result together with the request that makes it valid.
    #[cfg(test)]
    pub fn publish_gibbs_evaluation(
        &mut self,
        conditions: ThermoEvaluationConditions,
        compositions: HashMap<Option<String>, PhaseComposition>,
        values: HashMap<Option<String>, HashMap<String, f64>>,
    ) -> SubsDataResult<()> {
        self.publish_gibbs_evaluation_request(
            PhaseEvaluationRequest::new(conditions, compositions),
            values,
        )
    }

    /// Publishes a numeric Gibbs result together with the exact owned request.
    pub fn publish_gibbs_evaluation_request(
        &mut self,
        request: PhaseEvaluationRequest,
        values: HashMap<Option<String>, HashMap<String, f64>>,
    ) -> SubsDataResult<()> {
        self.validate_numeric_cache_publication(&request, &values, "Gibbs free energy")?;
        self.commit_phase_property(values, |bundle, property| bundle.gibbs.numeric = property);
        self.gibbs_cache_context = Some(NumericThermoCacheContext::new(
            request,
            self.configuration_revision,
        ));
        Ok(())
    }

    /// Entropy counterpart of `publish_gibbs_evaluation_request`.
    pub fn publish_entropy_evaluation_request(
        &mut self,
        request: PhaseEvaluationRequest,
        values: HashMap<Option<String>, HashMap<String, f64>>,
    ) -> SubsDataResult<()> {
        self.validate_numeric_cache_publication(&request, &values, "entropy")?;
        self.commit_phase_property(values, |bundle, property| bundle.entropy.numeric = property);
        self.entropy_cache_context = Some(NumericThermoCacheContext::new(
            request,
            self.configuration_revision,
        ));
        Ok(())
    }

    pub fn gibbs_cache_context(&self) -> Option<&NumericThermoCacheContext> {
        self.gibbs_cache_context.as_ref()
    }

    pub fn entropy_cache_context(&self) -> Option<&NumericThermoCacheContext> {
        self.entropy_cache_context.as_ref()
    }

    pub fn has_current_gibbs_cache(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> bool {
        self.gibbs_cache_context.as_ref().is_some_and(|context| {
            context.matches_parts(conditions, compositions, self.configuration_revision)
        })
    }

    pub fn has_current_gibbs_request(&self, request: &PhaseEvaluationRequest) -> bool {
        self.gibbs_cache_context
            .as_ref()
            .is_some_and(|context| context.matches_request(request, self.configuration_revision))
    }

    pub fn has_current_entropy_cache(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> bool {
        self.entropy_cache_context.as_ref().is_some_and(|context| {
            context.matches_parts(conditions, compositions, self.configuration_revision)
        })
    }

    pub fn has_current_entropy_request(&self, request: &PhaseEvaluationRequest) -> bool {
        self.entropy_cache_context
            .as_ref()
            .is_some_and(|context| context.matches_request(request, self.configuration_revision))
    }

    #[cfg(test)]
    pub(crate) fn state_snapshot(&self, multi_phase: bool) -> ThermoStateSnapshot<'_> {
        debug_assert!(!multi_phase, "PhaseSystem test snapshots are layout-driven");
        build_thermo_state_snapshot(
            self.layout_revision,
            build_multi_cache_snapshot(self.layout_revision, self.dG_view()),
            build_multi_cache_snapshot(self.layout_revision, self.dG_fun_view()),
            build_multi_cache_snapshot(self.layout_revision, self.dG_sym_view()),
            build_multi_cache_snapshot(self.layout_revision, self.dS_view()),
            build_multi_cache_snapshot(self.layout_revision, self.dS_fun_view()),
            build_multi_cache_snapshot(self.layout_revision, self.dS_sym_view()),
        )
    }

    /// Validates a numeric cache publication before it can replace any visible
    /// cache family. The request, phase keys, components, and values must all
    /// match the current canonical layout.
    fn validate_numeric_cache_publication(
        &self,
        request: &PhaseEvaluationRequest,
        values: &HashMap<Option<String>, HashMap<String, f64>>,
        property_name: &str,
    ) -> SubsDataResult<()> {
        let expected_component_counts = self
            .phase_data
            .iter()
            .map(|(phase, data)| (phase.clone(), data.substances.len()))
            .collect::<Vec<_>>();
        validate_phase_evaluation_request(
            expected_component_counts,
            request.conditions(),
            request.compositions(),
            &format!("{property_name} cache publication"),
        )?;

        let layout = self.current_layout();
        if values.len() != layout.phases().len() {
            return Err(SubsDataError::calculation_failed(
                property_name.to_string(),
                "numeric cache publication",
                format!(
                    "cache contains {} phases; expected {}",
                    values.len(),
                    layout.phases().len()
                ),
            ));
        }
        for phase in layout.phases() {
            let phase_name = phase.as_option();
            let phase_values =
                values
                    .get(phase_name)
                    .ok_or_else(|| SubsDataError::MissingData {
                        field: format!("{property_name} cache phase"),
                        substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                    })?;
            let components =
                layout
                    .components_for_phase(phase)
                    .ok_or_else(|| SubsDataError::MissingData {
                        field: "phase component layout".to_string(),
                        substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                    })?;
            if phase_values.len() != components.len() {
                return Err(SubsDataError::calculation_failed(
                    phase_name.clone().unwrap_or_else(|| "None".to_string()),
                    "numeric cache publication",
                    format!(
                        "{property_name} cache contains {} components; expected {}",
                        phase_values.len(),
                        components.len()
                    ),
                ));
            }
            for component in components {
                let value = phase_values.get(&component.substance).ok_or_else(|| {
                    SubsDataError::MissingData {
                        field: format!("{property_name} cache component"),
                        substance: component.label(),
                    }
                })?;
                if !value.is_finite() {
                    return Err(SubsDataError::calculation_failed(
                        component.label(),
                        "numeric cache publication",
                        format!("{property_name} cache contains a non-finite value"),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Validates that a derived property covers exactly the resolved phase
    /// payloads and their declared components before it can enter a cache.
    fn validate_phase_property_publication<T>(
        &self,
        values: &HashMap<Option<String>, HashMap<String, T>>,
        operation: &str,
    ) -> SubsDataResult<()> {
        if values.len() != self.phase_data.len() {
            return Err(SubsDataError::calculation_failed(
                operation.to_string(),
                "phase cache publication",
                format!(
                    "cache contains {} phases; expected {}",
                    values.len(),
                    self.phase_data.len()
                ),
            ));
        }
        for (phase, data) in &self.phase_data {
            let phase_values = values.get(phase).ok_or_else(|| {
                SubsDataError::calculation_failed(
                    operation.to_string(),
                    "phase cache publication",
                    format!("missing cache values for phase {phase:?}"),
                )
            })?;
            if phase_values.len() != data.substances.len() {
                return Err(SubsDataError::calculation_failed(
                    operation.to_string(),
                    "phase cache publication",
                    format!(
                        "phase {phase:?} contains {} cached components; expected {}",
                        phase_values.len(),
                        data.substances.len()
                    ),
                ));
            }
            for substance in &data.substances {
                if !phase_values.contains_key(substance) {
                    return Err(SubsDataError::calculation_failed(
                        substance.clone(),
                        "phase cache publication",
                        format!("{operation} is missing phase {phase:?}"),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Commits a property map that has already passed the canonical
    /// phase/component validation. No entry creation is permitted here.
    fn commit_phase_property<T>(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, T>>,
        mut assign: impl FnMut(&mut PhaseThermoCacheBundle, HashMap<String, T>),
    ) {
        for (phase, property) in values {
            let bundle = self
                .phase_bundles
                .get_mut(&phase)
                .expect("validated phase cache publication must have a canonical bundle");
            assign(bundle, property);
        }
        self.bump_cache_revision();
    }

    /// Test-only raw cache injection. Fixtures occasionally need to construct
    /// malformed or partially configured legacy states; production paths must
    /// use `validate_phase_property_publication` plus `commit_phase_property`.
    #[cfg(test)]
    fn replace_phase_property_for_test<T>(
        &mut self,
        values: HashMap<Option<String>, HashMap<String, T>>,
        mut assign: impl FnMut(&mut PhaseThermoCacheBundle, HashMap<String, T>),
    ) where
        T: Clone,
    {
        let mut phases = self.phase_bundles.keys().cloned().collect::<Vec<_>>();
        for phase in values.keys().cloned() {
            if !phases.contains(&phase) {
                phases.push(phase);
            }
        }
        for phase in phases {
            let bundle = self.bundle_or_insert(phase.clone());
            let property = values.get(&phase).cloned().unwrap_or_default();
            assign(bundle, property);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::phase_layout::PhaseId;
    use crate::Thermodynamics::physical_state::PhysicalState;

    #[test]
    fn derived_property_publication_rejects_unknown_phase_before_commit() {
        let mut system = PhaseSystem::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string()];
        system.replace_phase_data(HashMap::from([(Some("gas".to_string()), gas)]));
        let revision_before = system.layout_revision();

        let error = system
            .validate_phase_property_publication(
                &HashMap::from([(
                    Some("liquid".to_string()),
                    HashMap::from([("O2".to_string(), Expr::Const(1.0))]),
                )]),
                "symbolic Gibbs publication",
            )
            .expect_err("an unknown phase must not enter a production cache");

        assert!(error.to_string().contains("missing cache values for phase"));
        assert_eq!(system.layout_revision(), revision_before);
        assert!(system.dG_sym_view().get(&Some("gas".to_string())).is_some());
        assert!(
            system
                .dG_sym_view()
                .get(&Some("liquid".to_string()))
                .is_none()
        );
    }

    #[test]
    fn resolved_layout_rebuilds_bundles_and_legacy_edits_clear_typed_metadata() {
        let mut system = PhaseSystem::new();

        let gas_spec = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string()],
        )
        .unwrap();
        let liquid_spec = PhaseSpec::pure_condensed(
            PhaseId::new(Some("liquid".to_string())),
            vec!["H2O".to_string()],
            PhysicalState::Liquid,
        )
        .unwrap();

        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        let phase_data = HashMap::from([
            (Some("gas".to_string()), gas),
            (Some("liquid".to_string()), liquid),
        ]);
        let layout = SystemLayout::from_phase_map(&phase_data);

        system.configure_resolved_phases(vec![gas_spec.clone(), liquid_spec.clone()], layout);

        let mut bundle_keys: Vec<_> = system.phase_bundles().keys().cloned().collect();
        bundle_keys.sort();
        assert_eq!(
            bundle_keys,
            vec![Some("gas".to_string()), Some("liquid".to_string())]
        );
        assert_eq!(system.phase_specs(), &[gas_spec, liquid_spec]);
        assert!(system.resolved_layout().is_some());
        assert!(system.resolution_report().is_none());

        system.replace_phase_data(phase_data);

        let mut rebuilt_keys: Vec<_> = system.phase_bundles().keys().cloned().collect();
        rebuilt_keys.sort();
        assert_eq!(
            rebuilt_keys,
            vec![Some("gas".to_string()), Some("liquid".to_string())]
        );
        assert!(system.phase_specs().is_empty());
        assert!(system.resolved_layout().is_none());
        assert!(system.resolution_report().is_none());
    }
}
