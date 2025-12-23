use super::ClassicalThermodynamics::{Thermodynamics, ThermodynamicsError};
use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;

use std::collections::HashMap;

impl Thermodynamics {
    //////////////////////////////////////////dG//////////////////////////////////////////////
    pub fn extract_all_thermal_coeffs(
        &mut self,
        temperature: f64,
    ) -> Result<(), ThermodynamicsError> {
        let sd = &mut self.subdata;
        sd.extract_all_thermal_coeffs(temperature)?;
        Ok(())
    }
    pub fn calculate_therm_map_of_properties(
        &mut self,
        temperature: f64,
    ) -> Result<(), ThermodynamicsError> {
        let sd = &mut self.subdata;
        sd.calculate_therm_map_of_properties(temperature)?;
        Ok(())
    }
    pub fn calculate_therm_map_of_sym(&mut self) -> Result<(), ThermodynamicsError> {
        let sd = &mut self.subdata;
        sd.calculate_therm_map_of_sym()?;
        Ok(())
    }
    pub fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> Result<Vec<String>, ThermodynamicsError> {
        let sd = &mut self.subdata;
        let v = sd.extract_coeffs_if_current_coeffs_not_valid(temperature)?;
        Ok(v)
    }
    /// Function for calculating Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_Gibbs_free_energy(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        sd.calcutate_Gibbs_free_energy(T, P, n).unwrap();
    }

    /// Symbolic function for calculating Gibbs free energy of a given mixure of substances (symbolic function at given T, P, concentration)
    pub fn calculate_Gibbs_sym(&mut self, T: f64) {
        self.T = T;

        let sd = &mut self.subdata;
        sd.calculate_Gibbs_sym(T).unwrap();
    }
    /// Function for making Gibbs free energy closure of a given mixure of substances (as function of given T, P, concentration)
    pub fn calculate_Gibbs_fun(&mut self, T: f64) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        sd.calculate_Gibbs_fun(T, P);
    }
    ////////////////////////////////////////////S - ENTHROPY///////////////////////////////////////////////////////////
    /// Function for calculating enthropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_S(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
         sd.calculate_S(T, P, n).unwrap();
    }

    ///Symbolic  function for calculating enthropy of a given mixure of substances ( symbolic result at given T, P, concentration)
    pub fn calculate_S_sym(&mut self, T: f64) {
        self.T = T;

        let sd = &mut self.subdata;
        sd.calculate_S_sym(T).unwrap();
    }

    /// Function for making enthropy closure of a given mixure of substances (as function of given T, P, concentration)
    pub fn calculate_S_fun(&mut self, T: f64) {
        self.T = T;
        let sd = &mut self.subdata;
        sd.calculate_S_fun(T, self.P);
    }
}
