//! Compatibility API for historical all-in-one thermodynamics consumers.
//!
//! New code should use the narrow traits from `phase_interfaces`. This module
//! preserves the broad legacy surface while routing it through those focused
//! contracts and the shared `PhaseSystem` engine.

use std::collections::{HashMap, HashSet};

use RustedSciThe::symbolic::symbolic_engine::Expr;
use enum_dispatch::enum_dispatch;
use nalgebra::DMatrix;

use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::ordered_phase_entries;

use super::{
    CustomSubstance, PhaseComposition, PhaseDataPreparation, PhaseEquilibriumAssembly,
    PhaseOrSolution, PhaseSymbolicPropertyBuilder, SubstancesContainer,
};

/// Historical wide calculator API retained for compatibility.
///
/// The spelling of established methods remains unchanged here deliberately.
/// New solvers should depend on narrow typed capabilities instead.
#[enum_dispatch]
pub trait ThermodynamicsCalculatorTrait {
    fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()>;
    fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>>;
    fn calcutate_Gibbs_free_energy(
        &mut self,
        temperature: f64,
        pressure: f64,
        composition: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>>;
    fn calculate_Gibbs_sym(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn calculate_Gibbs_fun(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()>;
    fn set_P_to_sym_in_G_sym(&mut self, pressure: f64);
    fn set_T_to_sym_in_G_sym(&mut self, temperature: f64);
    fn calculate_S(
        &mut self,
        temperature: f64,
        pressure: f64,
        composition: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<()>;
    fn calculate_S_sym(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn calculate_S_fun(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()>;
    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()>;
    fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()>;
    fn extract_SubstancesContainer(&self) -> SubsDataResult<SubstancesContainer>;
    fn get_all_substances(&self) -> Vec<String>;
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)>;
    fn indexed_moles_variables(
        &mut self,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
        Vec<Expr>,
        Vec<Expr>,
        HashMap<Option<String>, HashMap<String, Expr>>,
    )>;
    fn calculate_Lagrange_equations_sym(
        &mut self,
        element_matrix: DMatrix<f64>,
        reference_temperature: f64,
    ) -> SubsDataResult<Vec<Expr>>;
    fn calculate_Lagrange_equations_fun2(
        &mut self,
        element_matrix: DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_temperature: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    >;
    fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
        HashMap<String, f64>,
    )>;
}

impl ThermodynamicsCalculatorTrait for PhaseOrSolution {
    fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.prepare_thermal_coefficients(temperature)
    }

    fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.rebuild_numeric_thermodynamic_properties(temperature)
    }

    fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()> {
        self.rebuild_symbolic_thermodynamic_properties()
    }

    fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        self.ensure_coefficients_for_temperature(temperature)
    }

    fn calcutate_Gibbs_free_energy(
        &mut self,
        temperature: f64,
        pressure: f64,
        composition: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.phase_system
            .calculate_gibbs_free_energy(temperature, pressure, composition)
    }

    fn calculate_Gibbs_sym(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.build_symbolic_gibbs(temperature)
    }

    fn calculate_Gibbs_fun(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()> {
        self.build_gibbs_functions(temperature, pressure)
    }

    fn set_P_to_sym_in_G_sym(&mut self, pressure: f64) {
        self.substitute_pressure_in_symbolic_gibbs(pressure);
    }

    fn set_T_to_sym_in_G_sym(&mut self, temperature: f64) {
        self.substitute_temperature_in_symbolic_gibbs(temperature);
    }

    fn calculate_S(
        &mut self,
        temperature: f64,
        pressure: f64,
        composition: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<()> {
        self.phase_system
            .calculate_entropy(temperature, pressure, composition)
    }

    fn calculate_S_sym(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.build_symbolic_entropy(temperature)
    }

    fn calculate_S_fun(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()> {
        self.build_entropy_functions(temperature, pressure)
    }

    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        self.configure_phase_properties(pressure, pressure_unit, molar_masses, mass_unit)
    }

    fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()> {
        self.fetch_missing_thermochemistry_from_nist()
    }

    fn extract_SubstancesContainer(&self) -> SubsDataResult<SubstancesContainer> {
        let mut phase_substances = HashMap::new();
        for (phase_name, subs_data) in ordered_phase_entries(self.phase_data()) {
            let phase_name = phase_name.0.ok_or_else(|| SubsDataError::MissingData {
                field: "phase name".to_string(),
                substance: "unknown".to_string(),
            })?;
            phase_substances.insert(phase_name, subs_data.substances.clone());
        }
        Ok(SubstancesContainer::MultiPhase(phase_substances))
    }

    fn get_all_substances(&self) -> Vec<String> {
        let mut substances = self
            .current_layout()
            .components()
            .iter()
            .map(|component| component.substance.clone())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        substances.sort();
        substances
    }

    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)> {
        self.phase_system.element_composition_and_molar_mass(groups)
    }

    fn indexed_moles_variables(
        &mut self,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
        Vec<Expr>,
        Vec<Expr>,
        HashMap<Option<String>, HashMap<String, Expr>>,
    )> {
        self.phase_system.indexed_moles_variables()
    }

    fn calculate_Lagrange_equations_sym(
        &mut self,
        element_matrix: DMatrix<f64>,
        reference_temperature: f64,
    ) -> SubsDataResult<Vec<Expr>> {
        self.build_symbolic_lagrange_equations(element_matrix, reference_temperature)
    }

    fn calculate_Lagrange_equations_fun2(
        &mut self,
        element_matrix: DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_temperature: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    > {
        self.build_numeric_lagrange_equations(
            element_matrix,
            temperature,
            pressure,
            reference_temperature,
        )
    }

    fn create_full_map_of_mole_numbers(
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
        self.phase_system
            .create_full_map_of_mole_numbers(non_zero_number_of_moles)
    }
}
