// SECTION UNDER CONSTRUCTION

use crate::Kinetics::molmass::calculate_molar_mass;
use std::collections::HashMap;

enum ErrorEnum {
    ParseError,
    MolarMassError,
    UnknownElement,
    UnknownLibrary,
}

/*
fn parse_and_calculate_molar_mass_wrapper(formula: &str) -> Result<(f64, HashMap<String, usize>), ErrorEnum> {
    let (molmass, atomic_composition) = calculate_molar_mass(formula.to_string());
    match (molmass,atomic_composition) {
        Ok( (molmass,atomic_composition) )  => Ok((molmass, atomic_composition)),


        Err(_) => Err(ErrorEnum::MolarMassError),
    }
}

 */
////////////////////////////////////////////////////////////////////////////////
#[derive(Debug, Clone)]
pub struct UserSubstances {
    pub vec_of_substances: Vec<String>,
    pub map_of_substances: HashMap<String, IndividualSubstance>,
    pub atomic_matrix: Vec<Vec<f64>>,
    pub parse_errors: HashMap<String, String>,
}
impl UserSubstances {
    pub fn new() -> Self {
        Self {
            vec_of_substances: Vec::new(),
            map_of_substances: HashMap::new(),
            atomic_matrix: Vec::new(),
            parse_errors: HashMap::new(),
        }
    }

    pub fn set_vec_of_substances(&mut self, vec_of_substances: Vec<String>) {
        self.vec_of_substances = vec_of_substances;
    }
    pub fn set_map_of_substances(&mut self) {
        for substance in &self.vec_of_substances {
            let mut individual_substance = IndividualSubstance::new();
            individual_substance.set_formula(substance.to_string());
            //  individual_substance.parse_and_calculate_molar_mass();
            self.map_of_substances
                .insert(substance.to_string(), individual_substance);
        }
    }
}
#[derive(Debug, Clone)]
// thermal data library
enum ThermalDataLib {
    NASA,
    NIST,
}

impl ThermalDataLib {
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "NASA" => Some(Self::NASA),
            "NIST" => Some(Self::NIST),
            _ => None,
        }
    }

    fn as_str(&self) -> &str {
        match self {
            Self::NASA => "NASA",
            Self::NIST => "NIST",
        }
    }
}

#[derive(Debug, Clone)]
enum TransportLib {
    Aramco,
    CEA,
}

impl TransportLib {
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "Aramco" => Some(Self::Aramco),
            "CEA" => Some(Self::CEA),
            _ => None,
        }
    }
    fn as_str(&self) -> &str {
        match self {
            Self::Aramco => "Aramco",
            Self::CEA => "CEA",
        }
    }
}

// struct for individual substance
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct IndividualSubstance {
    pub formula: String,
    pub molar_mass: f64,
    pub composition: HashMap<String, usize>,
    pub Cp: HashMap<ThermalDataLib, Vec<String>>, // heat capacity polynomial coefficients for substance
    pub dH: HashMap<ThermalDataLib, Vec<String>>, // heat of formation polynomial coefficients for substance
    pub dS: HashMap<ThermalDataLib, Vec<String>>, // entropy polynomial coefficients for substance
    pub dG: HashMap<ThermalDataLib, Vec<String>>, // Gibbs free energy polynomial coefficients for substance
    pub Lambda: HashMap<TransportLib, Vec<String>>, // thermal conductivity polynomial coefficients for substance
}
impl IndividualSubstance {
    pub fn new() -> Self {
        Self {
            formula: String::new(),
            molar_mass: 0.0,
            composition: HashMap::new(),
            Cp: HashMap::new(),
            dH: HashMap::new(),
            dS: HashMap::new(),
            dG: HashMap::new(),
            Lambda: HashMap::new(),
        }
    }
    pub fn set_formula(&mut self, formula: String) {
        self.formula = formula;
    }
    /*
    pub fn parse_and_calculate_molar_mass(&mut self) {
        match parse_and_calculate_molar_mass_wrapper (&self.formula) {
            Ok((molar_mass, composition)) => {
                self.molar_mass = molar_mass;
                self.composition = composition;
            }
            Err(_) => {
                self.molar_mass = 0.0;
                self.composition = HashMap::new();
            }
        }// end of match

    } // end of parse_and_calculate_molar_mass
    */
    pub fn collct_data_from_library(&mut self, _lib: ThermalDataLib) {}
} // end of IndividualSubstance
