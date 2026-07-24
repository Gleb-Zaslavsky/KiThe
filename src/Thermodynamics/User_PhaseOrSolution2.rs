//! # Single-Phase Thermodynamics Module
//!
//! This module provides the single-phase facade over the shared `PhaseSystem`
//! thermodynamic engine.
//!
//! ## Key Structure
//!
//! - [`OnePhase`]: Manages one `None`-keyed phase through the same canonical
//!   payload, cache, and pure-evaluation paths used by multi-phase systems.
//!
//! ## None Key Convention
//!
//! To maintain API compatibility with multi-phase systems, this module uses `None`
//! as the phase key in all HashMap returns. This allows seamless integration with
//! the unified [`CustomSubstance`] interface.
//!
//! ## Example Usage
//!
//! ```rust
//! use KiThe::Thermodynamics::User_PhaseOrSolution2::OnePhase;
//! use std::collections::HashMap;
//!
//! // Single-phase gas system
//! let system = OnePhase::new();
//! // Inspect resolved data with `system.subs_data_view()`.
//!
//! // Calculate properties - results use None key for compatibility
//! // let result = system.calculate_Gibbs_sym(298.15)?;
//! // assert!(result.contains_key(&None));
//! ```
//!
//! `OnePhase` is intentionally a small facade, not a separate state model.

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::fmt;

use crate::Thermodynamics::User_PhaseOrSolution::{
    MoleNumberSnapshot, PhaseComposition, PhaseDataPreparation, PhaseEquilibriumAssembly,
    PhaseEvaluationRequest, PhaseFunction, PhaseLagrangeFunction, PhaseSpec,
    PhaseSymbolicPropertyBuilder, ResolvedPhaseSystem, ResolvedPhaseSystemReport,
    SubstanceSystemFactoryError, SubstancesContainer, ThermoEvaluationConditions,
    ThermodynamicsCalculatorTrait,
};
use std::collections::HashMap;
use std::f64;

/// Single-phase thermodynamic system with optimized data access.
///
/// Provides direct access to substance data without phase-level indirection.
/// All HashMap returns use `None` as the key for API compatibility with multi-phase systems.
///
/// # Example
/// ```rust
/// use KiThe::Thermodynamics::User_PhaseOrSolution2::OnePhase;
///
/// let mut system = OnePhase::new();
/// // Inspect resolved data through `system.subs_data_view()`.
/// // Typed builders or crate-level transition helpers install components.
/// ```
#[derive(Clone)]
pub struct OnePhase {
    phase_system: crate::Thermodynamics::User_PhaseOrSolution::PhaseSystem,
}
impl OnePhase {
    /// Creates a new single-phase system with empty data.
    pub fn new() -> Self {
        OnePhase {
            phase_system:
                crate::Thermodynamics::User_PhaseOrSolution::PhaseSystem::new_single_phase(),
        }
    }

    /// Returns the current revision of the visible layout/cache state.
    pub fn layout_revision(&self) -> usize {
        self.phase_system.layout_revision()
    }

    /// Returns the resolved single-phase declaration when the system was
    /// created through the typed specification/factory path.
    pub fn phase_spec(&self) -> Option<&PhaseSpec> {
        self.phase_system.phase_specs().first()
    }

    /// Canonical declarations retained by the shared phase engine. The slice
    /// is empty for manually assembled legacy payloads until they are resolved
    /// through a typed `PhaseSpec`.
    pub fn phase_specs(&self) -> &[PhaseSpec] {
        self.phase_system.phase_specs()
    }

    /// Returns the canonical layout when this instance came from the typed
    /// resolution path rather than direct legacy field assembly.
    pub fn resolved_layout(&self) -> Option<&SystemLayout> {
        self.phase_system.resolved_layout()
    }

    /// Typed lookup provenance for the resolved payloads, when this system
    /// came from the typed resolution pipeline.
    pub fn resolution_report(&self) -> Option<&ResolvedPhaseSystemReport> {
        self.phase_system.resolution_report()
    }

    #[cfg(test)]
    pub(crate) fn set_phase_spec(&mut self, phase_spec: PhaseSpec) {
        let layout = self.current_layout();
        self.phase_system
            .configure_resolved_phases(vec![phase_spec], layout);
    }

    /// Installs a typed one-phase resolution in one commit. A multi-phase
    /// resolved payload is rejected here rather than silently becoming a
    /// partially represented `OnePhase` facade.
    pub(crate) fn install_resolved_system(
        &mut self,
        resolved: ResolvedPhaseSystem,
    ) -> Result<(), SubstanceSystemFactoryError> {
        let specs = resolved.phase_specs();
        if specs.len() != 1 || specs[0].id().as_option().is_some() {
            return Err(SubstanceSystemFactoryError::InvalidSpecification {
                field: "resolved phases".to_string(),
                reason: "OnePhase requires exactly one unnamed resolved phase".to_string(),
            });
        }
        self.phase_system.install_resolved_system(resolved);
        Ok(())
    }

    /// Borrowed access to the underlying single-phase database payload.
    pub fn subs_data_view(&self) -> &SubsData {
        self.phase_system.single_phase_data()
    }

    #[inline]
    fn subs_data(&self) -> &SubsData {
        self.phase_system.single_phase_data()
    }

    /// Replaces the one-phase payload through the crate-private migration
    /// boundary. Derived caches and typed metadata are discarded before a new
    /// phase specification can be installed.
    #[cfg(test)]
    pub(crate) fn replace_subs_data(&mut self, subs_data: SubsData) {
        self.phase_system.replace_single_phase_data(subs_data);
    }

    /// Executes one transitional raw-data edit, then invalidates the one-phase
    /// layout and all derived caches. It remains crate-private so external
    /// callers cannot desynchronise raw data from the resolved declaration.
    #[cfg(test)]
    pub(crate) fn with_subs_data_mut<R>(&mut self, update: impl FnOnce(&mut SubsData) -> R) -> R {
        self.phase_system.with_single_phase_data_mut(update)
    }

    /// Replaces the one-phase component list and invalidates all dependent
    /// layout and thermodynamic caches. Full production assembly normally
    /// enters through a typed `PhaseSpec`; this keeps crate fixtures explicit.
    #[cfg(test)]
    pub(crate) fn set_substances(&mut self, substances: Vec<String>) {
        self.with_subs_data_mut(|subs_data| {
            subs_data.substances = substances;
        });
    }

    /// Returns the canonical engine layout used for every one-phase derived
    /// vector and symbolic variable family.
    pub(crate) fn current_layout(&self) -> SystemLayout {
        self.phase_system.current_layout()
    }

    /// Shared symbolic bookkeeping owned by `PhaseSystem`. This crate-private
    /// view lets migration tests assert the common representation directly.
    #[cfg(test)]
    pub(crate) fn symbolic_layout_view(
        &self,
    ) -> &crate::Thermodynamics::User_PhaseOrSolution::PhaseSymbolicLayout {
        self.phase_system.symbolic_layout()
    }

    /// Returns a typed snapshot of the current state and caches.
    pub fn state_snapshot(
        &self,
    ) -> crate::Thermodynamics::User_PhaseOrSolution::ThermoStateSnapshot<'_> {
        use crate::Thermodynamics::User_PhaseOrSolution::{
            build_single_cache_snapshot, build_thermo_state_snapshot,
        };
        let layout_revision = self.layout_revision();
        build_thermo_state_snapshot(
            layout_revision,
            build_single_cache_snapshot(layout_revision, self.dG_view()),
            build_single_cache_snapshot(layout_revision, self.dG_fun_view()),
            build_single_cache_snapshot(layout_revision, self.dG_sym_view()),
            build_single_cache_snapshot(layout_revision, self.dS_view()),
            build_single_cache_snapshot(layout_revision, self.dS_fun_view()),
            build_single_cache_snapshot(layout_revision, self.dS_sym_view()),
        )
    }

    /// Returns a typed thermodynamic result snapshot aligned to a single-phase layout.
    pub fn result_snapshot(
        &self,
        temperature: Option<f64>,
        pressure: Option<f64>,
    ) -> crate::Thermodynamics::User_PhaseOrSolution::ThermoResultSnapshot<'_> {
        use crate::Thermodynamics::User_PhaseOrSolution::{
            build_single_cache_snapshot, build_thermo_result_snapshot,
        };
        let layout_revision = self.layout_revision();
        build_thermo_result_snapshot(
            layout_revision,
            self.current_layout(),
            temperature,
            pressure,
            build_single_cache_snapshot(layout_revision, self.dG_view()),
            build_single_cache_snapshot(layout_revision, self.dG_fun_view()),
            build_single_cache_snapshot(layout_revision, self.dG_sym_view()),
            build_single_cache_snapshot(layout_revision, self.dS_view()),
            build_single_cache_snapshot(layout_revision, self.dS_fun_view()),
            build_single_cache_snapshot(layout_revision, self.dS_sym_view()),
        )
    }

    /// Normalizes sparse mole input into a layout-carrying snapshot. New
    /// callers should prefer this over the legacy tuple returned by the broad
    /// compatibility trait.
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

    /// Borrowed Gibbs energy cache for read-only consumers.
    pub fn dG_view(&self) -> &HashMap<String, f64> {
        self.phase_system.single_dG_view()
    }

    /// Borrowed Gibbs function cache for read-only consumers.
    pub fn dG_fun_view(&self) -> &HashMap<String, PhaseFunction> {
        self.phase_system.single_dG_fun_view()
    }

    /// Borrowed symbolic Gibbs cache for read-only consumers.
    pub fn dG_sym_view(&self) -> &HashMap<String, Expr> {
        self.phase_system.single_dG_sym_view()
    }

    /// Borrowed entropy cache for read-only consumers.
    pub fn dS_view(&self) -> &HashMap<String, f64> {
        self.phase_system.single_dS_view()
    }

    /// Borrowed entropy function cache for read-only consumers.
    pub fn dS_fun_view(&self) -> &HashMap<String, PhaseFunction> {
        self.phase_system.single_dS_fun_view()
    }

    /// Borrowed symbolic entropy cache for read-only consumers.
    pub fn dS_sym_view(&self) -> &HashMap<String, Expr> {
        self.phase_system.single_dS_sym_view()
    }

    /// Returns the request that validates the published numeric Gibbs cache.
    pub fn gibbs_cache_context(
        &self,
    ) -> Option<&crate::Thermodynamics::User_PhaseOrSolution::NumericThermoCacheContext> {
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
    pub fn entropy_cache_context(
        &self,
    ) -> Option<&crate::Thermodynamics::User_PhaseOrSolution::NumericThermoCacheContext> {
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
    pub(crate) fn debug_dG_mut(&mut self) -> &mut HashMap<String, f64> {
        self.phase_system.dG_mut()
    }

    #[cfg(test)]
    pub(crate) fn debug_dG_fun_mut(&mut self) -> &mut HashMap<String, PhaseFunction> {
        self.phase_system.dG_fun_mut()
    }

    #[cfg(test)]
    pub(crate) fn debug_dG_sym_mut(&mut self) -> &mut HashMap<String, Expr> {
        self.phase_system.dG_sym_mut()
    }

    #[cfg(test)]
    pub(crate) fn debug_dS_fun_mut(&mut self) -> &mut HashMap<String, PhaseFunction> {
        self.phase_system.dS_fun_mut()
    }

    /// Evaluates Gibbs energy without publishing a phase-system cache.
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

    /// Delegates the pure query to the shared phase-system engine.
    pub fn evaluate_gibbs_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system
            .evaluate_gibbs_at(conditions, compositions)
    }

    /// Evaluates Gibbs energy from one self-contained request without cache publication.
    pub fn evaluate_gibbs_request(
        &self,
        request: &PhaseEvaluationRequest,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system.evaluate_gibbs_request(request)
    }

    /// Evaluates entropy without publishing a phase-system cache.
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

    /// Delegates the pure query to the shared phase-system engine.
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

    pub fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,

        Tm: f64,
    ) -> Result<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
        SubsDataError,
    > {
        self.phase_system
            .build_numeric_lagrange_equations(A, HashMap::from([(None, G_fun)]), Tm)
    }
}

impl fmt::Debug for OnePhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OnePhase")
            .field("subs_data", self.subs_data())
            .field("symbolic_layout", self.phase_system.symbolic_layout())
            .field("dG", self.phase_system.single_dG_view())
            .field("dG_sym", self.phase_system.single_dG_sym_view())
            .field("dS", self.phase_system.single_dS_view())
            .field("dS_sym", self.phase_system.single_dS_sym_view())
            .finish()
    }
}

impl PhaseDataPreparation for OnePhase {
    fn prepare_thermal_coefficients(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.phase_system.extract_all_thermal_coeffs(temperature)
    }

    fn rebuild_numeric_thermodynamic_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.phase_system
            .calculate_therm_map_of_properties(temperature)
    }

    fn rebuild_symbolic_thermodynamic_properties(&mut self) -> SubsDataResult<()> {
        self.phase_system.calculate_therm_map_of_sym()
    }

    fn ensure_coefficients_for_temperature(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        self.phase_system
            .extract_coeffs_if_current_coeffs_not_valid(temperature)
    }

    fn configure_phase_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        self.phase_system.configure_system_properties(
            pressure,
            pressure_unit,
            molar_masses,
            mass_unit,
        )
    }

    fn fetch_missing_thermochemistry_from_nist(&mut self) -> SubsDataResult<()> {
        self.phase_system.fetch_missing_from_nist()
    }
}

impl PhaseSymbolicPropertyBuilder for OnePhase {
    fn build_symbolic_gibbs(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.phase_system.calculate_gibbs_sym(temperature)
    }

    fn build_gibbs_functions(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()> {
        self.phase_system.calculate_gibbs_fun(temperature, pressure)
    }

    fn substitute_pressure_in_symbolic_gibbs(&mut self, pressure: f64) {
        self.phase_system.set_pressure_in_gibbs_sym(pressure);
    }

    fn substitute_temperature_in_symbolic_gibbs(&mut self, temperature: f64) {
        self.phase_system.set_temperature_in_gibbs_sym(temperature);
    }

    fn build_symbolic_entropy(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.phase_system.calculate_entropy_sym(temperature)
    }

    fn build_entropy_functions(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()> {
        self.phase_system
            .calculate_entropy_fun(temperature, pressure)
    }
}

impl PhaseEquilibriumAssembly for OnePhase {
    fn build_symbolic_lagrange_equations(
        &mut self,
        element_matrix: DMatrix<f64>,
        reference_temperature: f64,
    ) -> SubsDataResult<Vec<Expr>> {
        self.phase_system
            .build_symbolic_lagrange_equations(element_matrix, reference_temperature)
    }

    fn build_numeric_lagrange_equations(
        &mut self,
        element_matrix: DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_temperature: f64,
    ) -> SubsDataResult<PhaseLagrangeFunction> {
        self.phase_system.calculate_lagrange_equations_fun2(
            element_matrix,
            temperature,
            pressure,
            reference_temperature,
        )
    }
}

/// Implementation of ThermodynamicsCalculatorTrait for single-phase systems.
/// All methods maintain compatibility with multi-phase interface using None keys.
impl ThermodynamicsCalculatorTrait for OnePhase {
    /// Extracts thermal coefficients for the given temperature.
    fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        PhaseDataPreparation::prepare_thermal_coefficients(self, temperature)
    }

    /// Calculates thermodynamic property maps at given temperature.
    fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        PhaseDataPreparation::rebuild_numeric_thermodynamic_properties(self, temperature)
    }

    /// Creates symbolic expressions for thermodynamic properties.
    fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()> {
        PhaseDataPreparation::rebuild_symbolic_thermodynamic_properties(self)
    }

    /// Validates and extracts coefficients for temperature range.
    fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        PhaseDataPreparation::ensure_coefficients_for_temperature(self, temperature)
    }

    /// Calculates Gibbs free energy for given conditions and mole numbers.
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system.calculate_gibbs_free_energy(T, P, n)
    }

    /// Creates symbolic expressions for Gibbs free energy.
    fn calculate_Gibbs_sym(&mut self, T: f64) -> SubsDataResult<()> {
        PhaseSymbolicPropertyBuilder::build_symbolic_gibbs(self, T)
    }

    /// Sets pressure value in symbolic Gibbs expressions.
    fn set_P_to_sym_in_G_sym(&mut self, P: f64) {
        PhaseSymbolicPropertyBuilder::substitute_pressure_in_symbolic_gibbs(self, P);
    }

    /// Sets temperature value in symbolic Gibbs expressions.
    fn set_T_to_sym_in_G_sym(&mut self, T: f64) {
        PhaseSymbolicPropertyBuilder::substitute_temperature_in_symbolic_gibbs(self, T);
    }

    /// Creates Gibbs free energy functions for given conditions.
    fn calculate_Gibbs_fun(&mut self, T: f64, P: f64) -> SubsDataResult<()> {
        PhaseSymbolicPropertyBuilder::build_gibbs_functions(self, T, P)
    }

    /// Calculates entropy for given conditions and mole numbers.
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<()> {
        self.phase_system.calculate_entropy(T, P, n)
    }

    /// Creates symbolic expressions for entropy.
    fn calculate_S_sym(&mut self, T: f64) -> SubsDataResult<()> {
        PhaseSymbolicPropertyBuilder::build_symbolic_entropy(self, T)
    }

    /// Creates entropy functions for given conditions.
    fn calculate_S_fun(&mut self, T: f64, P: f64) -> SubsDataResult<()> {
        PhaseSymbolicPropertyBuilder::build_entropy_functions(self, T, P)
    }

    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        PhaseDataPreparation::configure_phase_properties(
            self,
            pressure,
            pressure_unit,
            molar_masses,
            mass_unit,
        )
    }
    fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()> {
        PhaseDataPreparation::fetch_missing_thermochemistry_from_nist(self)
    }
    fn extract_SubstancesContainer(&self) -> SubsDataResult<SubstancesContainer> {
        let substances = self.subs_data().substances.clone();
        Ok(SubstancesContainer::SinglePhase(substances))
    }
    fn get_all_substances(&self) -> Vec<String> {
        self.subs_data().substances.clone()
    }
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)> {
        self.phase_system.element_composition_and_molar_mass(groups)
    }
    fn indexed_moles_variables(
        &mut self,
    ) -> Result<
        (
            HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
            Vec<Expr>,
            Vec<Expr>,
            HashMap<Option<String>, HashMap<String, Expr>>,
        ),
        SubsDataError,
    > {
        self.phase_system.indexed_moles_variables()
    }
    /////////////////////////////equations for G->min///////////////////////////////////////
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>,
        Tm: f64,
    ) -> Result<Vec<Expr>, SubsDataError> {
        PhaseEquilibriumAssembly::build_symbolic_lagrange_equations(self, A, Tm)
    }

    fn calculate_Lagrange_equations_fun2(
        &mut self,
        A: DMatrix<f64>,
        T: f64,
        P: f64,
        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    > {
        PhaseEquilibriumAssembly::build_numeric_lagrange_equations(self, A, T, P, Tm)
    }
    //////////////////////////////equations for Keq///////////////////////////////////////
    /*
        fn calculate_K_eq_sym(
            &mut self,
            stolich_matrix: DMatrix<f64>,

            map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>>,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let T = Expr::Var("T".to_string());

            let subs_data = self.calculate_Gibbs_sym_one_phase(298.0, None, None);
            let mut K_vec: Vec<Expr> = Vec::with_capacity(r);
            for j in 0..r {
                let mut dG_reaction_j = Expr::Const(0.0);
                for (i, subs_i) in substances.iter().enumerate() {
                    let dG_i = subs_data.get(subs_i).unwrap().clone();
                    let nu_ij = stolich_matrix[(i, j)];
                    dG_reaction_j += Expr::Const(nu_ij) * dG_i;
                }
                let ln_K = -dG_reaction_j / (R_sym * T.clone());
                K_vec.push(ln_K);
            }
        }

        fn calculate_K_eq_fun(
            &mut self,
            stolich_matrix: DMatrix<f64>,
            G_fun: HashMap<
                Option<String>,
                HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
            >,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let P = self.P.unwrap_or(101325.0);
            let dG = self.calculate_Gibbs_fun_one_phase(P, 298.0);

            let stol = stolich_matrix; // move into closure
            let K_vec: Box<dyn Fn(f64) -> Vec<f64> + 'static> = Box::new(move |T: f64| -> Vec<f64> {
                let mut K_vec: Vec<f64> = Vec::with_capacity(r);
                for j in 0..r {
                    let mut dG_reaction_j = 0.0;
                    for (i, subs_i) in substances.iter().enumerate() {
                        let dG_i = dG
                            .get(subs_i)
                            .expect("Missing Gibbs function for substance");
                        let nu_ij = stol[(i, j)];
                        dG_reaction_j += nu_ij * dG_i(T, None, None);
                    }
                    let ln_K = -dG_reaction_j / (R * T);
                    K_vec.push(ln_K);
                }
                K_vec
            });
        }
        fn calculate_K_eq_equation_sym(
            &mut self,
            stolich_matrix: DMatrix<f64>,

            map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>>,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let T = Expr::Var("T".to_string());
            let map_of_vars = map_of_var_each_substance.get(&None).unwrap().clone();
            let subs_data = self.calculate_Gibbs_sym_one_phase(298.0, None, None);
            let mut residual: Vec<Expr> = Vec::new();
            for j in 0..r {
                let mut dG_reaction_j = Expr::Const(0.0);
                let mut sum_ln_n = Expr::Const(0.0);
                for (i, subs_i) in substances.iter().enumerate() {
                    let dG_i = subs_data.get(subs_i).unwrap().clone();
                    let nu_ij = Expr::Const(stolich_matrix[(i, j)]);
                    dG_reaction_j += nu_ij.clone() * dG_i;

                    let n_i = map_of_vars.get(subs_i).unwrap();
                    sum_ln_n += nu_ij * n_i.ln();
                }
                let ln_K = -dG_reaction_j / (R_sym * T.clone());
                let residual_j = sum_ln_n - sum_ln_n_phase - ln_K;
                residual.push(ln_K);
            }
        }
        fn calculate_K_eq_equation_fun(
            &mut self,
            stolich_matrix: DMatrix<f64>,
            G_fun: HashMap<
                Option<String>,
                HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
            >,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let P = self.P.unwrap_or(101325.0);
            let dG = self.calculate_Gibbs_fun_one_phase(P, 298.0);

            let stol = stolich_matrix; // move into closure
            let K_vec: Box<dyn Fn(f64, Vec<f64>, Option<f64>) -> Vec<f64> + 'static> =
                Box::new(move |T: f64, n: Vec<f64>, Np: Option<f64>| -> Vec<f64> {
                    let mut K_vec: Vec<f64> = Vec::with_capacity(r);
                    let mut sum_ln_n = 0.0;
                    for j in 0..r {
                        let mut dG_reaction_j = 0.0;
                        for (i, subs_i) in substances.iter().enumerate() {
                            let dG_i = dG
                                .get(subs_i)
                                .expect("Missing Gibbs function for substance");
                            let nu_ij = stol[(i, j)];
                            dG_reaction_j += nu_ij * dG_i(T, None, None);

                            sum_ln_n += nu_ij * n[i].ln();
                        }
                        let ln_K = -dG_reaction_j / (R * T.clone());
                        let residual_j = sum_ln_n - sum_ln_n_phase - ln_K;
                        residual.push(ln_K);
                    }
                    K_vec
                });
        }
    */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
        // physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
    ) -> Result<
        (
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
            HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
            HashMap<String, f64>,
        ),
        SubsDataError,
    > {
        self.phase_system
            .create_full_map_of_mole_numbers(non_zero_number_of_moles)
    }
}
