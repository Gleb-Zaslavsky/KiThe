//! Compatibility facade over single- and multi-phase thermodynamic systems.
//!
//! `CustomSubstance` deliberately owns no calculation state. It normalizes
//! read-only access across `OnePhase` and `PhaseOrSolution` while the shared
//! `PhaseSystem` remains the engine beneath both variants.

use std::collections::HashMap;

use RustedSciThe::symbolic::symbolic_engine::Expr;
use enum_dispatch::enum_dispatch;
use nalgebra::DMatrix;

use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;

use super::{
    MoleNumberSnapshot, NestedPhaseCacheView, PhaseDataPreparation, PhaseEquilibriumAssembly,
    PhaseLagrangeFunction, PhaseLayoutAccess, PhaseOrSolution, PhaseSymbolicPropertyBuilder,
    ResolvedPhaseSystemReport, SubstancesContainer, ThermoCacheSnapshot, ThermoResultSnapshot,
    ThermoStateSnapshot, ThermodynamicsCalculatorTrait, build_cache_snapshot,
};

/// Unified compatibility interface for single- and multi-phase systems.
///
/// New code should prefer the narrow phase traits; this enum preserves one
/// read-only facade for existing consumers that must accept either shape.
#[derive(Debug, Clone)]
#[enum_dispatch(ThermodynamicsCalculatorTrait)]
pub enum CustomSubstance {
    OnePhase(OnePhase),
    PhaseOrSolution(PhaseOrSolution),
}

impl CustomSubstance {
    pub fn system_layout(&self) -> SystemLayout {
        match self {
            Self::OnePhase(one_phase) => one_phase.system_layout(),
            Self::PhaseOrSolution(phase_or_solution) => phase_or_solution.system_layout(),
        }
    }

    pub fn layout_revision(&self) -> usize {
        match self {
            Self::OnePhase(one_phase) => one_phase.layout_revision(),
            Self::PhaseOrSolution(phase_or_solution) => phase_or_solution.layout_revision(),
        }
    }

    pub fn state_snapshot(&self) -> ThermoStateSnapshot<'_> {
        match self {
            Self::OnePhase(one_phase) => one_phase.state_snapshot(),
            Self::PhaseOrSolution(phase_or_solution) => phase_or_solution.state_snapshot(),
        }
    }

    pub fn result_snapshot(
        &self,
        temperature: Option<f64>,
        pressure: Option<f64>,
    ) -> ThermoResultSnapshot<'_> {
        match self {
            Self::OnePhase(one_phase) => one_phase.result_snapshot(temperature, pressure),
            Self::PhaseOrSolution(phase_or_solution) => {
                phase_or_solution.result_snapshot(temperature, pressure)
            }
        }
    }

    pub fn dG_snapshot(&self) -> ThermoCacheSnapshot<'_, f64> {
        build_cache_snapshot(self.layout_revision(), self.get_dG_view())
    }

    pub fn dG_sym_snapshot(&self) -> ThermoCacheSnapshot<'_, Expr> {
        build_cache_snapshot(self.layout_revision(), self.get_dG_sym_view())
    }

    pub fn dS_snapshot(&self) -> ThermoCacheSnapshot<'_, f64> {
        build_cache_snapshot(self.layout_revision(), self.get_dS_view())
    }

    pub fn dS_sym_snapshot(&self) -> ThermoCacheSnapshot<'_, Expr> {
        build_cache_snapshot(self.layout_revision(), self.get_dS_sym_view())
    }

    /// Typed lookup provenance for the resolved payloads, if this facade came
    /// from the typed resolution pipeline.
    pub fn resolution_report(&self) -> Option<&ResolvedPhaseSystemReport> {
        match self {
            Self::OnePhase(one_phase) => one_phase.resolution_report(),
            Self::PhaseOrSolution(phase_or_solution) => phase_or_solution.resolution_report(),
        }
    }

    /// Returns legacy structural input without flattening phase identity.
    pub fn extract_SubstancesContainer(&self) -> SubsDataResult<SubstancesContainer> {
        match self {
            Self::OnePhase(subs_data) => Ok(SubstancesContainer::SinglePhase(
                subs_data.subs_data_view().substances.clone(),
            )),
            Self::PhaseOrSolution(phase_or_solution) => {
                let layout = phase_or_solution.system_layout();
                let mut phase_substances = HashMap::new();
                for phase in layout.phases() {
                    let phase_key =
                        phase
                            .as_option()
                            .clone()
                            .ok_or_else(|| SubsDataError::MissingData {
                                field: "phase name".to_string(),
                                substance: "unknown".to_string(),
                            })?;
                    let components = layout.components_for_phase(phase).ok_or_else(|| {
                        SubsDataError::MissingData {
                            field: "phase component layout".to_string(),
                            substance: phase_key.clone(),
                        }
                    })?;
                    phase_substances.insert(
                        phase_key,
                        components
                            .iter()
                            .map(|component| component.substance.clone())
                            .collect(),
                    );
                }
                Ok(SubstancesContainer::MultiPhase(phase_substances))
            }
        }
    }

    pub fn get_ordered_component_labels(&self) -> SubsDataResult<Vec<(Option<String>, String)>> {
        Ok(match self {
            Self::OnePhase(one_phase) => one_phase
                .subs_data_view()
                .substances
                .iter()
                .cloned()
                .map(|substance| (None, substance))
                .collect(),
            Self::PhaseOrSolution(phase_or_solution) => phase_or_solution
                .system_layout()
                .components()
                .iter()
                .map(|component| {
                    let phase_key = component.phase.as_option().clone().ok_or_else(|| {
                        SubsDataError::MissingData {
                            field: "phase name".to_string(),
                            substance: "unknown".to_string(),
                        }
                    })?;
                    Ok((Some(phase_key), component.substance.clone()))
                })
                .collect::<SubsDataResult<Vec<_>>>()?,
        })
    }

    pub fn get_ordered_substances(&self) -> SubsDataResult<Vec<String>> {
        Ok(self
            .get_ordered_component_labels()?
            .into_iter()
            .map(|(_, substance)| substance)
            .collect())
    }

    pub fn get_ordered_result_labels(&self) -> SubsDataResult<Vec<String>> {
        Ok(match self {
            Self::OnePhase(one_phase) => one_phase.subs_data_view().substances.clone(),
            Self::PhaseOrSolution(phase_or_solution) => {
                phase_or_solution.system_layout().component_labels()
            }
        })
    }

    pub fn normalize_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
    ) -> SubsDataResult<MoleNumberSnapshot> {
        match self {
            Self::OnePhase(one_phase) => one_phase.normalize_mole_numbers(non_zero_number_of_moles),
            Self::PhaseOrSolution(phase_or_solution) => {
                phase_or_solution.normalize_mole_numbers(non_zero_number_of_moles)
            }
        }
    }

    /// Compatibility projection of the typed mole-number snapshot.
    pub fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
        HashMap<String, f64>,
    )> {
        match self {
            Self::OnePhase(one_phase) => {
                one_phase.create_full_map_of_mole_numbers(non_zero_number_of_moles)
            }
            Self::PhaseOrSolution(phase_or_solution) => {
                phase_or_solution.create_full_map_of_mole_numbers(non_zero_number_of_moles)
            }
        }
    }

    pub fn get_dG_sym_view(&self) -> NestedPhaseCacheView<'_, Expr> {
        match self {
            Self::OnePhase(one_phase) => NestedPhaseCacheView::Single(one_phase.dG_sym_view()),
            Self::PhaseOrSolution(phase_or_solution) => {
                NestedPhaseCacheView::Multi(phase_or_solution.dG_sym_view())
            }
        }
    }

    pub fn get_dG_sym(&self) -> NestedPhaseCacheView<'_, Expr> {
        self.get_dG_sym_view()
    }

    pub fn get_dG_view(&self) -> NestedPhaseCacheView<'_, f64> {
        match self {
            Self::OnePhase(one_phase) => NestedPhaseCacheView::Single(one_phase.dG_view()),
            Self::PhaseOrSolution(phase_or_solution) => {
                NestedPhaseCacheView::Multi(phase_or_solution.dG_view())
            }
        }
    }

    pub fn get_dG(&self) -> NestedPhaseCacheView<'_, f64> {
        self.get_dG_view()
    }

    pub fn get_dS_sym_view(&self) -> NestedPhaseCacheView<'_, Expr> {
        match self {
            Self::OnePhase(one_phase) => NestedPhaseCacheView::Single(one_phase.dS_sym_view()),
            Self::PhaseOrSolution(phase_or_solution) => {
                NestedPhaseCacheView::Multi(phase_or_solution.dS_sym_view())
            }
        }
    }

    pub fn get_dS_view(&self) -> NestedPhaseCacheView<'_, f64> {
        match self {
            Self::OnePhase(one_phase) => NestedPhaseCacheView::Single(one_phase.dS_view()),
            Self::PhaseOrSolution(phase_or_solution) => {
                NestedPhaseCacheView::Multi(phase_or_solution.dS_view())
            }
        }
    }
}

// The enum owns no separate state. These implementations are intentionally
// transparent forwarding so callers can use narrow capabilities without
// matching on the compatibility facade themselves.
impl PhaseDataPreparation for CustomSubstance {
    fn prepare_thermal_coefficients(&mut self, temperature: f64) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.prepare_thermal_coefficients(temperature),
            Self::PhaseOrSolution(value) => value.prepare_thermal_coefficients(temperature),
        }
    }

    fn rebuild_numeric_thermodynamic_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.rebuild_numeric_thermodynamic_properties(temperature),
            Self::PhaseOrSolution(value) => {
                value.rebuild_numeric_thermodynamic_properties(temperature)
            }
        }
    }

    fn rebuild_symbolic_thermodynamic_properties(&mut self) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.rebuild_symbolic_thermodynamic_properties(),
            Self::PhaseOrSolution(value) => value.rebuild_symbolic_thermodynamic_properties(),
        }
    }

    fn ensure_coefficients_for_temperature(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        match self {
            Self::OnePhase(value) => value.ensure_coefficients_for_temperature(temperature),
            Self::PhaseOrSolution(value) => value.ensure_coefficients_for_temperature(temperature),
        }
    }

    fn configure_phase_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => {
                value.configure_phase_properties(pressure, pressure_unit, molar_masses, mass_unit)
            }
            Self::PhaseOrSolution(value) => {
                value.configure_phase_properties(pressure, pressure_unit, molar_masses, mass_unit)
            }
        }
    }

    fn fetch_missing_thermochemistry_from_nist(&mut self) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.fetch_missing_thermochemistry_from_nist(),
            Self::PhaseOrSolution(value) => value.fetch_missing_thermochemistry_from_nist(),
        }
    }
}

impl PhaseSymbolicPropertyBuilder for CustomSubstance {
    fn build_symbolic_gibbs(&mut self, temperature: f64) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.build_symbolic_gibbs(temperature),
            Self::PhaseOrSolution(value) => value.build_symbolic_gibbs(temperature),
        }
    }

    fn build_gibbs_functions(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.build_gibbs_functions(temperature, pressure),
            Self::PhaseOrSolution(value) => value.build_gibbs_functions(temperature, pressure),
        }
    }

    fn substitute_pressure_in_symbolic_gibbs(&mut self, pressure: f64) {
        match self {
            Self::OnePhase(value) => value.substitute_pressure_in_symbolic_gibbs(pressure),
            Self::PhaseOrSolution(value) => value.substitute_pressure_in_symbolic_gibbs(pressure),
        }
    }

    fn substitute_temperature_in_symbolic_gibbs(&mut self, temperature: f64) {
        match self {
            Self::OnePhase(value) => value.substitute_temperature_in_symbolic_gibbs(temperature),
            Self::PhaseOrSolution(value) => {
                value.substitute_temperature_in_symbolic_gibbs(temperature)
            }
        }
    }

    fn build_symbolic_entropy(&mut self, temperature: f64) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.build_symbolic_entropy(temperature),
            Self::PhaseOrSolution(value) => value.build_symbolic_entropy(temperature),
        }
    }

    fn build_entropy_functions(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()> {
        match self {
            Self::OnePhase(value) => value.build_entropy_functions(temperature, pressure),
            Self::PhaseOrSolution(value) => value.build_entropy_functions(temperature, pressure),
        }
    }
}

impl PhaseEquilibriumAssembly for CustomSubstance {
    fn build_symbolic_lagrange_equations(
        &mut self,
        element_matrix: DMatrix<f64>,
        reference_temperature: f64,
    ) -> SubsDataResult<Vec<Expr>> {
        match self {
            Self::OnePhase(value) => {
                value.build_symbolic_lagrange_equations(element_matrix, reference_temperature)
            }
            Self::PhaseOrSolution(value) => {
                value.build_symbolic_lagrange_equations(element_matrix, reference_temperature)
            }
        }
    }

    fn build_numeric_lagrange_equations(
        &mut self,
        element_matrix: DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_temperature: f64,
    ) -> SubsDataResult<PhaseLagrangeFunction> {
        match self {
            Self::OnePhase(value) => value.build_numeric_lagrange_equations(
                element_matrix,
                temperature,
                pressure,
                reference_temperature,
            ),
            Self::PhaseOrSolution(value) => value.build_numeric_lagrange_equations(
                element_matrix,
                temperature,
                pressure,
                reference_temperature,
            ),
        }
    }
}
