#![allow(warnings)]
pub mod kinetics;
mod mechfinder;

use RustedSciThe::symbolic::symbolic_engine::Expr;
/// The module is equipped with a library of kinetic parameters of chemical reactions obtained as a result of parsing publicly available databases
/// The module takes as input the name of the library and the vector of substances and then produces the following data:
/// 1) all reactions of starting substances with each other, and all their possible products with each other.
/// 2) HashMap with kinetic data of all found reactions
/// ru
/// Модуль снабжен библиотекой кинетических параметров химических реакций, полученной в результате парсинга общедоступныхбаз данных
/// Модуль берет на вход название библиотеки и вектор веществ а затем выдает следующие данные:
/// 1) все реакции исходных веществ между собой, и всех их возможных продуктов между собой.
/// 2) HashMap с кинетическими данными всех найденных реакций
/// ----------------------------------------------------------------
// git config --global pack.windowmemory 80000000
use kinetics::{ElementaryStruct, FalloffStruct, PressureStruct, ThreeBodyStruct};
use log::{error, info, warn};
use serde::de::{self, Deserializer, MapAccess, Visitor};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Number, Value};
use std::collections::{HashMap, HashSet};
use std::f64;
use std::fmt;
/// enum for types of chemical kinetics rate contant functions
#[derive(Debug, PartialEq, Serialize, Clone)]
#[serde(rename_all = "lowercase")]
pub enum ReactionType {
    Elem,
    Falloff,
    #[serde(rename = "pres")]
    Pressure,
    #[serde(rename = "three-body")]
    ThreeBody,
    Empirical,
}
impl<'de> Deserialize<'de> for ReactionType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.as_str() {
            "elem" => Ok(ReactionType::Elem),
            "falloff" => Ok(ReactionType::Falloff),
            "pressure" | "pres" => Ok(ReactionType::Pressure),
            "three-body" | "threebody" => Ok(ReactionType::ThreeBody),
            "empirical" => Ok(ReactionType::Empirical),
            _ => Err(serde::de::Error::custom(format!(
                "Unknown reaction type: {}",
                s
            ))),
        }
    }
}
/// struct for reaction data
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ReactionData {
    #[serde(rename = "type")]
    pub reaction_type: ReactionType,
    pub eq: String,
    pub react: Option<HashMap<String, f64>>,
    #[serde(flatten)]
    pub data: ReactionKinetics,
}
/// enum for structs of different types of kinetics
#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(untagged)]
pub enum ReactionKinetics {
    Falloff(FalloffStruct),
    Pressure(PressureStruct),
    ThreeBody(ThreeBodyStruct),
    Elementary(ElementaryStruct),
}
impl ReactionData {
    ///function that checks if the reaction_type field matches the variant of the data field in a ReactionData instance.
    /// If they don't match, the function will panic
    pub fn validate_reaction_type(&self) {
        match (&self.reaction_type, &self.data) {
            (ReactionType::Elem, ReactionKinetics::Elementary(_)) => {}
            (ReactionType::Falloff, ReactionKinetics::Falloff(_)) => {}
            (ReactionType::Pressure, ReactionKinetics::Pressure(_)) => {}
            (ReactionType::ThreeBody, ReactionKinetics::ThreeBody(_)) => {}
            _ => panic!(
                "Mismatch between reaction_type ({:?}) and data variant ({:?})",
                self.reaction_type,
                std::mem::discriminant(&self.data)
            ),
        }
    }
    /// function to calculate reaction rate constant (under thr hood there is a function that calculates the rate constant for each type of reaction)
    /// ATTENTION! don't forget to use absolute temperature in Kelvin!
    pub fn K_const(
        &self,
        temp: f64,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
    ) -> f64 {
        match self.data.clone() {
            ReactionKinetics::Elementary(data) => data.K_const(temp),
            ReactionKinetics::Falloff(data) => data.K_const(temp, concentrations.clone().unwrap()),
            ReactionKinetics::Pressure(data) => data.K_const(temp, pres.clone().unwrap()),
            ReactionKinetics::ThreeBody(data) => {
                data.K_const(temp, concentrations.clone().unwrap())
            }
        }
    }
    ///function to calculate reaction rate constant for the range of temperatures: from T0 to Tend, number of points is n
    /// ATTENTION! don't forget to use absolute temperature in Kelvin!
    pub fn K_const_for_T_range(
        &self,
        T0: f64,
        Tend: f64,
        n: usize,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
    ) -> Vec<f64> {
        let mut K_const_values = Vec::new();
        let T: Vec<f64> = (0..n)
            .map(|i| T0 + i as f64 * (Tend - T0) / n as f64)
            .collect();
        for Ti in T {
            let k = self.K_const(Ti, pres, concentrations.clone());
            K_const_values.push(k);
        }
        K_const_values
    }
    /// generate the symbolic representation of the reaction rate constant
    pub fn K_sym(&self, pres: Option<f64>, concentrations: Option<HashMap<String, Expr>>) -> Expr {
        match self.data.clone() {
            ReactionKinetics::Elementary(data) => data.K_expr(),
            ReactionKinetics::Falloff(data) => data.K_expr(concentrations.clone().unwrap()),
            ReactionKinetics::Pressure(data) => data.K_expr(pres.clone().unwrap()),
            ReactionKinetics::ThreeBody(data) => data.K_expr(concentrations.clone().unwrap()),
        }
    }
}
pub fn parse_kinetic_data(
    big_mech: &str,
    vec_of_reactions: &[String],
    vec_of_reaction_value: Vec<Value>,
) -> (Map<String, Value>, Vec<String>) {
    let mut reaction_data_hash = Map::new();
    let mut equations = Vec::new();
    info!("______________PARCING REACTION DATA INTO STRUCTS________");
    for (reaction_record, reaction_id) in vec_of_reaction_value.iter().zip(vec_of_reactions) {
        info!("reaction_record {:#?} \n \n ", reaction_record);
        let react_code = format!("{}_{}", big_mech, reaction_id);
        if let Ok(mut reactiondata) =
            serde_json::from_value::<ReactionData>(reaction_record.clone())
        {
            equations.push(reactiondata.eq.clone());
            // let reacttype =  &reactiondata.reaction_type;
            reactiondata.validate_reaction_type(); //hecks if the reaction_type field matches the variant of the data field in a ReactionData instance. If they don't match, the function will panic
            let value = serde_json::to_value(&reactiondata).unwrap();
            reaction_data_hash.insert(react_code, value);
        } else {
            error!("Error parsing reaction: {}", reaction_record);
            panic!("Error parsing reaction: {}", reaction_record);
        }
    }
    info!("______________PARCING REACTION DATA INTO STRUCTS ENDED________");
    (reaction_data_hash, equations)
}
/// parse Vec of serde Values with reaction data
pub fn parse_kinetic_data_vec(
    vec_of_reaction_value: Vec<Value>,
) -> (Vec<ReactionData>, Vec<String>) {
    println!("\n \n______________PARCING REACTION DATA INTO STRUCTS________");
    let mut reaction_dat = Vec::new();
    let mut equations = Vec::new();

    for reaction_record in vec_of_reaction_value.iter() {
        println!("reaction_record {:#?} \n \n ", reaction_record);
        if let Ok(mut reactiondata) =
            serde_json::from_value::<ReactionData>(reaction_record.clone())
        {
            equations.push(reactiondata.eq.clone());
            println!("parsed into {:#?} \n", reactiondata);
            reactiondata.validate_reaction_type(); //hecks if the reaction_type field matches the variant of the data field in a ReactionData instance. If they don't match, the function will panic
            reaction_dat.push(reactiondata);
        } else {
            info!("Error parsing reaction: {}", reaction_record);
            panic!("Error parsing reaction: {}", reaction_record);
        }
    }
    info!("______________PARCING REACTION DATA INTO STRUCTS ENDED________");
    (reaction_dat, equations)
}
/// struct for chemical mechanism construction
#[derive(Debug)]
pub struct Mechanism_search {
    pub task_substances: Vec<String>,
    pub task_library: String,
    pub mechanism: Vec<String>,
    pub reactants: Vec<String>,
    pub vec_of_reactions: Vec<String>,
    pub reactdata: Vec<ReactionData>,
}
// set the task to construct mechanism
impl Mechanism_search {
    pub fn new(task_substances: Vec<String>, task_library: String) -> Self {
        Self {
            task_substances: task_substances,
            task_library: task_library,
            mechanism: Vec::new(),
            reactants: Vec::new(),
            vec_of_reactions: Vec::new(),
            reactdata: Vec::new(),
        }
    }

    pub fn default() -> Self {
        Self {
            task_substances: Vec::new(),
            task_library: String::new(),
            mechanism: Vec::new(),
            reactants: Vec::new(),
            vec_of_reactions: Vec::new(),
            reactdata: Vec::new(),
        }
    }
    /// find chemical mechanism using mechfinder API
    pub fn mechfinder_api(&mut self) -> (Vec<String>, Vec<String>, Vec<String>) {
        /*
        let tuple = [ "O", "NH3", "NO", "O2", "N2", "N2O", "CO", "C"];
        O,NH3,NO,O2,N2,N2O,CO,C
        let big_mech = "NUIG".to_string();
            let vec: Vec<&str> =  tuple.into_iter().collect();
        */

        let vec: Vec<&str> = self.task_substances.iter().map(|s| s.as_str()).collect();
        let big_mech = self.task_library.clone();
        info!("задание {:?}, библиотека {:?}", &big_mech, &vec);
        let (mechanism, reactants, vec_of_reactions, vec_of_reaction_value) =
            mechfinder::mechfinder(&big_mech, vec);
        // info!("mechanism {:?}", &mechanism);
        let (mut reactdata, vec_of_equations) =
            parse_kinetic_data_vec(vec_of_reaction_value.clone()); // парсим данные о реакциях
        self.mechanism = mechanism;
        self.reactants = reactants;
        self.vec_of_reactions = vec_of_reactions;
        self.reactdata = reactdata.clone();
        return (
            self.mechanism.to_owned(),
            self.reactants.to_owned(),
            self.vec_of_reactions.to_owned(),
        );
    }
}

//tests
const ELEM_TESTING_JSON: &str = r#"{"type": "elem",
                 "eq": "NAPH+C2H3<=>NAPHV+C2H4",
                  "Arrenius": [0.408, 4.02, 36822.949]}"#;
const FALOFF_TESTING_JSON: &str = r#" {"type": "falloff",
                 "eq": "C4H71-3+CH3(+M)<=>C5H10-2(+M)",
                 "low_rate": [3.91e+60, -12.81, 26143.75],
                "high_rate": [100000000000000.0, -0.32, -1097.2009],
                 "eff": {"H2": 2.0, "H2O": 6.0, "CH4": 2.0, "CO": 1.5, "CO2": 2.0, "C2H6": 3.0, "AR": 0.7},
                 "troe": [0.104, 1606.0, 60000.0, 6118.0]} "#;
const PRES_TESTING_JSON: &str = r#"{"type": "pres", "eq": "SC4H9<=>C3H6+CH3",
               'Arrenius': {"0.001": [2.89e+40, -9.76, 140552.983],
                             "0.01": [1.8e+44, -10.5, 154800.281],
                             "0.1": [2.51e+46, -10.73, 168311.37099999998],
                             "1.0": [4.74e+44, -9.85, 175020.903], 
                             "10.0": [3.79e+37, -7.44, 169846.532],
                             "100.0": [4.79e+26, -4.01, 154344.334] }}"#;
const THREE_BODY_TESTING_JSON: &str = r#"{"type": "threebody",
              "eq": "H2+M<=>H+H+M",
              "Arrenius": [4.577e+19, -1.4, 436705.19999999995],
             "eff": {"H2": 2.5, "H2O": 12.0, "CO": 1.9, "CO2": 3.8, "HE": 0.83, "CH4": 2.0, "C2H6": 3.0} }"#;
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mechfinder_api() {
        let mut mech_search = Mechanism_search::new(
            vec!["O".to_string(), "NH3".to_string(), "NO".to_string()],
            "NUIG".to_string(),
        );

        let (mechanism, reactants, vec_of_reactions) = mech_search.mechfinder_api();

        assert!(!mechanism.is_empty());
        assert!(!reactants.is_empty());
        assert!(!vec_of_reactions.is_empty());
    }

    #[test]
    fn test_default_values() {
        let mech_search = Mechanism_search::default();

        assert!(mech_search.task_substances.is_empty());
        assert!(mech_search.task_library.is_empty());
        assert!(mech_search.mechanism.is_empty());
        assert!(mech_search.reactants.is_empty());
        assert!(mech_search.vec_of_reactions.is_empty());
    }

    #[test]

    fn test_ELEM_parse_kinetic_data() {
        let big_mech: &str = "NUIG";
        /*
        let test_data = [ELEM_TESTING_JSON, FALOFF_TESTING_JSON, PRES_TESTING_JSON, THREE_BODY_TESTING_JSON];
        let   test_reactions_numbers = vec!("1", "2532", "1736", "5");
        let vec_of_reactions: Vec<String> =  test_reactions_numbers.iter().map(|&s| s.trim().to_string()).collect();
        let vec_of_reaction_value: Vec<Value> = test_data.iter().map(|&s| serde_json::from_str(&s).unwrap()).collect();
        */
        let vec_of_reactions = vec!["1".to_string()];
        let reaction = ELEM_TESTING_JSON;
        let vec_of_reaction_value: Vec<Value> = vec![serde_json::from_str(reaction).unwrap()];
        let (ReactionDataHash, _) =
            parse_kinetic_data(big_mech, &vec_of_reactions, vec_of_reaction_value);

        assert!(!ReactionDataHash.is_empty());
        //   let elem_saved_to_hash = ReactionDataHash[test_reactions_numbers[0]];
        let key = format!("{}_{}", big_mech, &vec_of_reactions[0]);
        let elem_react_testing_instance: ElementaryStruct =
            serde_json::from_value::<ElementaryStruct>(ReactionDataHash[&key].clone()).unwrap();
        println!("K_const {:?}", elem_react_testing_instance.K_const(298.15));
        assert!(elem_react_testing_instance.K_const(298.15) > 0.0);
    }
    #[test]

    fn test_THREEBODY_parse_kinetic_data() {
        let ThreeBodyStruct_test_data: &str = r#"{"Arrenius": [4.577e+19, -1.4, 436705.19999999995],
       "eff": {"H2": 2.5, "H2O": 12.0, "CO": 1.9, "CO2": 3.8, "HE": 0.83, "CH4": 2.0, "C2H6": 3.0} }"#;
        let big_mech: &str = "NUIG";
        let vec_of_reactions = vec!["2".to_string()];
        let reaction = THREE_BODY_TESTING_JSON;
        let vec_of_reaction_value: Vec<Value> = vec![serde_json::from_str(reaction).unwrap()];
        let (ReactionDataHash, _) =
            parse_kinetic_data(big_mech, &vec_of_reactions, vec_of_reaction_value);

        assert!(!ReactionDataHash.is_empty());
        println!("ReactionDataHash: {:?} \n \n", ReactionDataHash);
        let key = format!("{}_{}", big_mech, &vec_of_reactions[0]);

        let threebody_react_testing_instance: ThreeBodyStruct =
            serde_json::from_value::<ThreeBodyStruct>(
                serde_json::from_str(ThreeBodyStruct_test_data).unwrap(),
            )
            .unwrap();
        println!(
            "threebody_react_testing_instance {:?}",
            threebody_react_testing_instance
        );
        let mut Concentrations: HashMap<String, f64> = HashMap::new();
        Concentrations.insert("H".to_string(), 0.5);
        Concentrations.insert("O".to_string(), 0.5);
        assert!(threebody_react_testing_instance.K_const(298.15, Concentrations) > 0.0);
        // assert!(elem_react_testing_instance.K_const(298.15) > 0.0);
    }
    #[test]
    fn test_THREEBODY_from_lib() {
        use crate::Kinetics::User_reactions::KinData;
        let mut kinetics = KinData::new();
        let C1_react = Some(vec!["C1".to_string()]);
        kinetics.shortcut_reactions = C1_react.clone();

        kinetics.get_reactions_from_shortcuts();

        kinetics.reactdata_parsing();
        assert!(kinetics.vec_of_reaction_data.iter().len() > 0);
    }
    #[test]

    fn test_FALOFF_parse_kinetic_data() {
        let big_mech: &str = "NUIG";
        let vec_of_reactions = vec!["3".to_string()];
        let reaction = FALOFF_TESTING_JSON;
        let vec_of_reaction_value: Vec<Value> = vec![serde_json::from_str(reaction).unwrap()];
        let (ReactionDataHash, _) =
            parse_kinetic_data(big_mech, &vec_of_reactions, vec_of_reaction_value);
        println!("ReactionDataHash: {:#?}", ReactionDataHash);
        assert!(!ReactionDataHash.is_empty());
        let key = format!("{}_{}", big_mech, &vec_of_reactions[0]);
        let falloff_react_testing_instance: FalloffStruct =
            serde_json::from_value::<FalloffStruct>(ReactionDataHash[&key].clone()).unwrap();
        let mut Concentrations: HashMap<String, f64> = HashMap::new();
        Concentrations.insert("H".to_string(), 0.5);
        Concentrations.insert("O".to_string(), 0.5);
        assert!(falloff_react_testing_instance.K_const(298.15, Concentrations) > 0.0);
        // assert!(elem_react_testing_instance.K_const(298.15) > 0.0);
    }

    #[test]
    fn test_three_body_deserialization() {
        use serde_json::json;
        let three_body_json = json!({
            "type": "three-body",
            "eq": "H2+M<=>H+H+M",
            "Arrenius": [4.577e+19, -1.4, 436705.19999999995],
            "eff": {
                "H2": 2.5,
                "H2O": 12.0,
                "CO": 1.9,
                "CO2": 3.8,
                "HE": 0.83,
                "CH4": 2.0,
                "C2H6": 3.0
            }
        });

        let reaction_data: ReactionData = serde_json::from_value(three_body_json).unwrap();
        println!("reaction_data: {:#?}", reaction_data);
        assert_eq!(
            reaction_data.reaction_type,
            ReactionType::ThreeBody,
            "wrong reaction type!"
        );
        if let ReactionKinetics::ThreeBody(_) = reaction_data.data {
            // Success
        } else {
            panic!("Expected ThreeBody variant");
        }
    }
    #[test]
    fn test_pres_deserialization() {
        let reaction_data: ReactionData =
            serde_json::from_str(PRES_TESTING_JSON).expect("Error parsing JSON: {err:?}");
        println!("reaction_data: {:#?}", reaction_data);
        assert_eq!(
            reaction_data.reaction_type,
            ReactionType::Pressure,
            "wrong reaction type!"
        );
        if let ReactionKinetics::Pressure(_) = reaction_data.data {
            // Success
        } else {
            panic!("Expected Pressure variant");
        }
    }
    #[test]
    fn test_pres_data_deserialization() {
        const PRES_TESTING_JSON: &str = r#"{
                        "Arrenius":{"0.01": [2.89e+40, -9.76, 140552.983],
                                    "0.1": [1.8e+44, -10.5, 154800.281]
                                    }}"#;

        let pres_val =
            serde_json::from_str(PRES_TESTING_JSON).expect("Error parsing JSON: {err:?}");
        println!("val: {:#?}", pres_val);
        let pres_data = serde_json::from_value::<PressureStruct>(pres_val);
        if let Ok(pres) = pres_data {
            println!("pres_data: {:#?}", pres);
            // Success
        } else {
            panic!("Expected Pressure variant");
        }
    }
}
