//! Narrow capability interfaces for phase-system consumers.
//!
//! Solvers rarely need the historical all-in-one calculator API. These traits
//! describe the smaller contracts for layout inspection, pure evaluation,
//! payload preparation, symbolic construction, and equilibrium assembly.

use std::collections::HashMap;

use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
use crate::Thermodynamics::User_substances_error::SubsDataResult;
use crate::Thermodynamics::phase_layout::SystemLayout;

use super::{
    PhaseComposition, PhaseDataView, PhaseLagrangeFunction, PhaseOrSolution, PhaseSpec,
    ThermoEvaluationConditions,
};

/// Read-only phase identity and solver component order.
pub trait PhaseLayoutAccess {
    /// Active component ordering. Typed resolution takes precedence; legacy
    /// payloads use a deterministic derived layout as a compatibility path.
    fn system_layout(&self) -> SystemLayout;

    /// Typed phase declarations when construction went through resolution.
    fn resolved_phase_specs(&self) -> &[PhaseSpec];
}

/// Pure numeric thermodynamic property evaluation.
pub trait PhasePropertyEvaluator {
    fn evaluate_gibbs_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>>;

    fn evaluate_entropy_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>>;
}

/// Preparation and lookup operations that intentionally mutate phase payloads.
pub trait PhaseDataPreparation {
    fn prepare_thermal_coefficients(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn rebuild_numeric_thermodynamic_properties(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn rebuild_symbolic_thermodynamic_properties(&mut self) -> SubsDataResult<()>;
    fn ensure_coefficients_for_temperature(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>>;
    fn configure_phase_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()>;
    fn fetch_missing_thermochemistry_from_nist(&mut self) -> SubsDataResult<()>;
}

/// Symbolic thermodynamic property construction, distinct from solver assembly.
pub trait PhaseSymbolicPropertyBuilder {
    fn build_symbolic_gibbs(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn build_gibbs_functions(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()>;
    fn substitute_pressure_in_symbolic_gibbs(&mut self, pressure: f64);
    fn substitute_temperature_in_symbolic_gibbs(&mut self, temperature: f64);
    fn build_symbolic_entropy(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn build_entropy_functions(&mut self, temperature: f64, pressure: f64) -> SubsDataResult<()>;
}

/// Adapter-facing Lagrange-equation assembly.
pub trait PhaseEquilibriumAssembly {
    fn build_symbolic_lagrange_equations(
        &mut self,
        element_matrix: DMatrix<f64>,
        reference_temperature: f64,
    ) -> SubsDataResult<Vec<Expr>>;
    fn build_numeric_lagrange_equations(
        &mut self,
        element_matrix: DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_temperature: f64,
    ) -> SubsDataResult<PhaseLagrangeFunction>;
}

impl PhaseLayoutAccess for PhaseOrSolution {
    fn system_layout(&self) -> SystemLayout {
        self.current_layout()
    }

    fn resolved_phase_specs(&self) -> &[PhaseSpec] {
        self.phase_specs()
    }
}

impl PhasePropertyEvaluator for PhaseOrSolution {
    fn evaluate_gibbs_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        PhaseOrSolution::evaluate_gibbs_at(self, conditions, compositions)
    }

    fn evaluate_entropy_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        PhaseOrSolution::evaluate_entropy_at(self, conditions, compositions)
    }
}

impl PhaseDataPreparation for PhaseOrSolution {
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

impl PhaseSymbolicPropertyBuilder for PhaseOrSolution {
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

impl PhaseEquilibriumAssembly for PhaseOrSolution {
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

impl PhaseLayoutAccess for OnePhase {
    fn system_layout(&self) -> SystemLayout {
        self.resolved_layout()
            .cloned()
            .unwrap_or_else(|| PhaseDataView::single(self.subs_data_view()).layout())
    }

    fn resolved_phase_specs(&self) -> &[PhaseSpec] {
        self.phase_specs()
    }
}

impl PhasePropertyEvaluator for OnePhase {
    fn evaluate_gibbs_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        OnePhase::evaluate_gibbs_at(self, conditions, compositions)
    }

    fn evaluate_entropy_at(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        OnePhase::evaluate_entropy_at(self, conditions, compositions)
    }
}
