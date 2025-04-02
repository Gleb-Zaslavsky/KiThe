// module containing basic formulae of chemical thermodynamics
use crate::Kinetics::molmass::create_elem_composition_matrix;
use crate::Thermodynamics::User_substances::{DataType, SubsData};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;
pub struct Thermodynamics {
    pub vec_of_substances: Vec<String>,
    pub T: f64,
    pub P: f64,

    pub subdata: SubsData,
    pub dmu: HashMap<DataType, f64>,
    pub dmu_fun: HashMap<DataType, Box<dyn Fn(f64) -> f64>>,
    pub dmu_sym: HashMap<DataType, Expr>,
}

impl Thermodynamics {
    pub fn new() -> Self {
        Self {
            vec_of_substances: Vec::new(),
            T: 298.15,
            P: 1e5,
            subdata: SubsData::new(),
            dmu: HashMap::new(),
            dmu_fun: HashMap::new(),
            dmu_sym: HashMap::new(),
        }
    }
    pub fn set_T(&mut self, T: f64) {
        self.T = T;
    }
    pub fn set_P(&mut self, P: f64) {
        self.P = P;
    }
}

pub fn calculate_Gibbs_free_energy(
    T: f64,
    P: f64,
    vec_of_dH: Vec<f64>,
    vec_of_dS: Vec<f64>,
) -> f64 {
    let n = vec_of_dH.len();
    let mut Gibbs_free_energy = 0.0;
    for i in 0..n {
        Gibbs_free_energy += vec_of_dH[i] - T * vec_of_dS[i];
    }
    return Gibbs_free_energy;
}
