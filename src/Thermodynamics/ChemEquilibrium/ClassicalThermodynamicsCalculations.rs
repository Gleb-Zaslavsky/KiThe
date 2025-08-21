use super::ClassicalThermodynamics::Thermodynamics;
use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;

impl Thermodynamics {
    //////////////////////////////////////////dG//////////////////////////////////////////////

    /// Function for calculating Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_Gibbs_free_energy(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        self.dG = sd.calcutate_Gibbs_free_energy(T, P, n).unwrap();
    }

    /// Symbolic function for calculating Gibbs free energy of a given mixure of substances (symbolic function at given T, P, concentration)
    pub fn calculate_Gibbs_sym(&mut self, T: f64) {
        self.T = T;
        let n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)> =
            self.symbolic_vars.clone();
        let sd = &mut self.subdata;
        self.dG_sym = sd.calculate_Gibbs_sym(T, n).unwrap();
    }
    /// Function for making Gibbs free energy closure of a given mixure of substances (as function of given T, P, concentration)
    pub fn calculate_Gibbs_fun(&mut self, T: f64) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        self.dG_fun = sd.calculate_Gibbs_fun(T, P);
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
        self.dS = sd.calculate_S(T, P, n).unwrap();
    }

    ///Symbolic  function for calculating enthropy of a given mixure of substances ( symbolic result at given T, P, concentration)
    pub fn calculate_S_sym(&mut self, T: f64) {
        self.T = T;
        let n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)> =
            self.symbolic_vars.clone();
        let sd = &mut self.subdata;
        self.dS_sym = sd.calculate_S_sym(T, n).unwrap();
    }

    /// Function for making enthropy closure of a given mixure of substances (as function of given T, P, concentration)
    pub fn calculate_S_fun(&mut self, T: f64) {
        self.T = T;
        let sd = &mut self.subdata;
        self.dS_fun = sd.calculate_S_fun(T, self.P);
    }
}
