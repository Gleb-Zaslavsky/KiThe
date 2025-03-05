//mod mechfinder_api;
use crate::Kinetics::kinetics_lib_api::KineticData;
use crate::Kinetics::mechfinder_api::{ parse_kinetic_data_vec, ReactionData, Mechanism_search, ReactionKinetics};
use crate::Kinetics::parsetask::{decipher_vector_of_shortcuts, decipher_vector_of_shortcuts_to_pairs};
use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
use std::fs::File;
use std::io::Write;
use serde_json::json;
use serde_json::Value;
use std::collections::HashMap;
use prettytable::{Table, Row, Cell, row};

use log::{ info, warn};
/// THE STRUCT KinData COLLECTS ALL THE INFORMATION ABOUT SPECIFIC REACTIONS, WHICH ARE NEEDED FOR FURTHER CALCULATIONS.
///  so this is API for allmost all features of Kinetics module
/// Not all features can be used simultaneously. But the list of features is as follows:
///
/// 1) collectiong reactions from library by their library adress or by shortcut names
/// 2) Stoicheometric  data structures: matrix of stoicheometric coefficients, matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the
/// kinetic function, G_matrix. As a rule, the degrees of concentration in the kinetic function coincide with the stoicheometric coefficients of
/// the substances in the reaction; however, for empirical reactions they may differ from the stoicheometric coefficients.
/// 3) calculate molar mass of a substances of reacants and products
/// 4) calculate matrix of elemets
/// In short, the KinData structure is used for processing of user-chosen reactions.
///  So you define what reactions you need by using the constructor of mechanism from the mechfinder_api module or
/// manually with help of kinetics_lib_api module. Now you can

///structure to store user task and reaction data
#[derive(Debug, Clone)]
pub struct KinData {
    pub shortcut_reactions: Option<Vec<String>>, // vector of reaction shortcut names
    pub map_of_reactions: Option<HashMap<String, Vec<String>>>, // full "address" of reaction {'library':"id of reaction"} group reactions by library names {'library':[reaction_ids in that library]}
    pub vec_of_pairs: Option<Vec<(String, String)>>, // vector of pairs of reaction library names and reaction ids
    pub vec_of_reaction_Values: Option<Vec<Value>>, // vector of reaction data in the form of serde Values
    pub vec_of_reaction_data: Option<Vec<ReactionData>>, // data of all reactions
    pub vec_of_equations: Vec<String>,                        // vector of equations of reactions
    pub substances: Vec<String>,                              // vector of substance names
    pub groups: Option<HashMap<String, HashMap<String, usize>>>, // Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
    pub stecheodata: StoichAnalyzer, // matrix of stoichiometric coefficients and other matrices
}

impl KinData {
    pub fn new() -> Self {
        Self {
            shortcut_reactions: None,
            map_of_reactions: None,
            vec_of_pairs: None,
            vec_of_reaction_Values: None,
            vec_of_reaction_data: None,
            vec_of_equations: Vec::new(),
            substances: Vec::new(),
            groups: None,
            stecheodata: StoichAnalyzer::new(),
        }
    }
   /////////////////////////////////SETTING REACTIONS/////////////////////////////////////////// 
    /// set reactions directly
    pub fn set_reactions_directly(
        &mut self,
        reactions: Vec<String>,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> () {
        self.shortcut_reactions = Some(reactions);
        self.groups = groups;
    }
    /// turns range of shortcuts String like C1..20 to Vec of strings С1, C2, C3...etc...C20
    pub fn set_reactions_from_shortcut_range(&mut self, shortcut_range: String) -> Vec<String> {
        let parts: Vec<&str> = shortcut_range.split("..").collect(); //splits the string by the ".." to get the starting and ending parts.

        if parts.len() != 2 {
            panic!("Invalid shortcut range: {}", shortcut_range);
    
        }

        let prefix = parts[0]
            .chars()
            .take_while(|c| c.is_alphabetic())
            .collect::<String>(); //extracts the prefix (letters) and the numeric portion from the start.
   
        let start: usize = parts[0][prefix.len()..].parse().unwrap_or(0); // then parses the numeric values and generates the range of strings.
        let end: usize = if let Some(last_char) = parts[1].chars().last() { //  checks if there's a last character in parts[1]
            if last_char.is_numeric() { // the last character is numerit
                parts[1].chars().rev() // Reverses the string
                    .take_while(|c| c.is_numeric())// Takes characters while they are numeric
                    .collect::<String>() //Collects these into a string
                    .chars().rev().collect::<String>()// Reverses this string again to get the number in the correct order.
                    .parse().unwrap_or(0) //Parses this string into a usize.
            } else {
                panic!("last symbol must be numeric: {}", shortcut_range);
            }
        } else {
            panic!("last symbol must be numeric: {}", shortcut_range);
        };
        let res:Vec<String> = (start..=end).map(|i| format!("{}_{}", prefix, i)).collect(); //formatted strings are collected into a Vec<String>.
        info!("task includes reactions as follows: {:#?}", &res);
        self.shortcut_reactions = Some(res.clone());
        return res;
    }

    /// decifer vector of shortcuts to full reaction names and store them in map {'library':"id of reaction"} and Vec<("library, reaction_id")>
    pub fn get_reactions_from_shortcuts(&mut self) -> () {
        if let Some(shortcut_reactions) = &self.shortcut_reactions {
            let vec: Vec<&str> = shortcut_reactions.iter().map(|s| s.as_str()).collect();
            self.map_of_reactions = Some(decipher_vector_of_shortcuts(vec.clone()));
            let vec_of_pairs = decipher_vector_of_shortcuts_to_pairs(vec);
            self.vec_of_pairs = Some(vec_of_pairs.clone());
            // let's now open the library and get reactions from it
            let mut vec_of_reaction_values = Vec::new();
                // instance of KineticData with opened library json files
                let mut kin_instance = KineticData::new();
                for (lib, reaction_id) in vec_of_pairs.iter() {
                    // collecting reaction data for each library name
                    kin_instance.open_json_files(lib);
    
                    let reaction_data_Value =
                    kin_instance.search_reactdata_by_reaction_id(&reaction_id);
                    
                    vec_of_reaction_values.push(reaction_data_Value);
                  
                    // now we have a vector of json objects with reaction data and a vector of reaction ids
                    // lets parse it to map of structures
          
                }
                self.vec_of_reaction_Values = Some(vec_of_reaction_values);

        } else {
            warn!("KinData::create_map_of_reactions: shortcut_reactions is None");
        }
    }
    ///construct reaction mechanism
    pub fn construct_mechanism(&mut self,    task_substances: Vec<String>, task_library: String,) {

        let mut found_mech = Mechanism_search::new(   task_substances,task_library.clone());
        found_mech.mechfinder_api();
        self.vec_of_equations = found_mech.vec_of_reactions;
  
      //  println!("mechanism found: reaction data: \n {:?}", found_mech.reactdata);
      //  println!("mechanism found: mechanism: \n {:#?}", found_mech.mechanism);
      //  println!("mechanism found: reactants: \n {:#?}", found_mech.reactants);
      //  self.vec_of_equations = found_mech.vec_of_reactions;
        self.vec_of_reaction_data = Some(found_mech.reactdata);
        let reactions = found_mech.mechanism;
        let mut full_addres = Vec::new();
        for reaction in reactions.iter() {
            let addres = format!("{}_{}", task_library, reaction);
            full_addres.push(addres);
        }
        self.shortcut_reactions = Some(full_addres);
        self.substances = found_mech.reactants;
    

    }
    /////////////////////////////////REACTIONS MANIPULATIONS//////////////////////////////////////
    /// add manially reaction data as serde Value
    pub fn append_reaction(&mut self, reactions:Vec<Value>) {
        let mut old_reactions = self.vec_of_reaction_Values.as_mut().unwrap().clone();
        old_reactions.extend(reactions);
        self.vec_of_reaction_Values = Some(old_reactions);
    }
    pub fn remove_by_index(&mut self, index: usize) {
    
            let mut vec_of_equations = self.vec_of_equations.clone();
            vec_of_equations.remove(index);
            self.vec_of_equations = vec_of_equations;
            if let Some(vec_of_reaction_Values) = self.vec_of_reaction_Values.clone().as_mut() {
                vec_of_reaction_Values.remove(index);
            self.vec_of_reaction_Values = Some(vec_of_reaction_Values.clone());
            }
            if let Some(vec_of_reaction_data) = self.vec_of_reaction_data.clone().as_mut() {
                vec_of_reaction_data.remove(index);
            self.vec_of_reaction_data = Some(vec_of_reaction_data.clone());
            }  
            if let Some(vec_of_pairs) = self.vec_of_pairs.clone().as_mut() {
                vec_of_pairs.remove(index);
            self.vec_of_pairs = Some(vec_of_pairs.clone());
            }  
    }
    pub fn remove_reaction_by_name(&mut self, reaction_name: &str) {
       let i=  &self.vec_of_equations.clone().iter().position(|eq_i| *eq_i == reaction_name);
       if let Some(index) = i {
           self.remove_by_index(*index);
       }
    }
  
    /////////////////////////////////COMPUTING AND PARSING REACTION DATA///////////////////////////////////////////
    /// parse reaction libraries and extracts data into structs
    pub fn reactdata_parsing(&mut self) -> () {
        if let Some(vec_of_reaction_values ) = & self.vec_of_reaction_Values  {

            let (vec_ReactionData, vec_of_equations) =
                    parse_kinetic_data_vec( vec_of_reaction_values.clone());
          //  info!("map_of_reaction_data: {:?}", &vec_of_reaction_data);
         //   info!("vec_of_equations: {:?}", &vec_of_equations);
            self.vec_of_reaction_data = Some(vec_ReactionData);
            self.vec_of_equations = vec_of_equations;
        } else {
            warn!("KinData::reactdata_from_shortcuts: map_of_reactions is None");
        }
    }
    /// generates  Stoicheometric  data structures: matrix of stoicheometric coefficients, matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the
    /// kinetic function, G_matrix. As a rule, the degrees of concentration in the kinetic function coincide with the stoicheometric coefficients of
    /// the substances in the reaction; however, for empirical reactions they may differ from the stoicheometric coefficients.
    pub fn analyze_reactions(&mut self) -> () {
        // iniciate instance of ReactionAnalyzer
        let mut StoichAnalyzer_instance = StoichAnalyzer::new();
        // copy vector of reactions to ReactionAnalyzer_instance
        StoichAnalyzer_instance.reactions = self.vec_of_equations.clone();
        StoichAnalyzer_instance.groups = self.groups.clone();
        // parse to find substance names
        if self.substances.is_empty()| (self.substances.len() == 0) { // no need to parse for substances if they are already known
        StoichAnalyzer_instance.search_substances();
        self.substances = StoichAnalyzer_instance.substances.clone();
        } else {StoichAnalyzer_instance.substances = self.substances.clone()};
        //find stoichiometric matrix and other matrices
        StoichAnalyzer_instance.analyse_reactions();
        StoichAnalyzer_instance.create_matrix_of_elements();
        self.stecheodata = StoichAnalyzer_instance;
    }
     /// parsing reaction data into structures and stoichometric calculations under one hood
    pub fn kinetic_main(&mut self) {
        self.reactdata_parsing();
        self.analyze_reactions();
    }
    ///////////////////////////INPUT/OUTPUT/////////////////////////////////////////////////////////
    /// print the chosen reactions 
    pub fn print_raw_reactions(&self) -> Result<(), std::io::Error> {
        if let Some(vec_of_reaction_Values) = & self.vec_of_reaction_Values {
             // Convert the vector of Values to a JSON array
            let json_array = json!(vec_of_reaction_Values);
            // Write the JSON array to a file
            let mut file = File::create("raw_reactions.json")?;
            file.write_all(serde_json::to_string_pretty(&json_array)?.as_bytes())?;
            info!("Raw reactions have been written to raw_reactions.json");
            Ok(())
    
        } else {
            warn!("KinData::reactdata_from_shortcuts: map_of_reactions is None");
            Ok(())
        }
    }
   
    pub fn pretty_print_substances_verbose(&self) -> Result<(), std::io::Error> {
        let subs = &self.substances;
        let reacts = &self.vec_of_equations;
        let sd: &StoichAnalyzer = &self.stecheodata.clone();
        let elem_matrix = sd.matrix_of_elements.clone().unwrap().clone();
        let unique_vec_of_elems = sd.unique_vec_of_elems.clone().unwrap();
        let sm = sd.stecheo_matrx.clone();
        let n_react_from_matrix = sm.len();
        assert_eq!(n_react_from_matrix, reacts.len());
        let n_subs_from_matrix = sm[0].len();
        assert_eq!(n_subs_from_matrix, subs.len());
        let (nrows, ncols) = elem_matrix.shape();
        assert_eq!(nrows, subs.len());
        assert_eq!(ncols, unique_vec_of_elems.len());
         // Code  table of such structure 1) first row: first element - String: "Reactions/Substances" all other elements of the first row taken
    //  from vector subs 2) nest rows look like as follows first element taken from vec reacts, and other elements of row taken from stecheomatrix
    //  vec of vectors
        let mut table_stecheo = Table::new();
        // Create the header row
        let mut header_row = vec![Cell::new("Reactions/Substances")];
        for sub in subs {
            header_row.push(Cell::new(sub));
        }
        table_stecheo.add_row(Row::new(header_row));
    
        // Create the subsequent rows
        for (i, react) in reacts.iter().enumerate() {
            let mut row = vec![Cell::new(react)];
            for j in 0..subs.len() {
                row.push(Cell::new(&format!("{:.4}", sm[i][j])));
            }
            table_stecheo.add_row(Row::new(row));
        }
    
        // Print the table
        table_stecheo.printstd();
    

        let mut elem_table = Table::new();
        let mut header_row = vec![Cell::new("Substances/Elements")];
        for elems in unique_vec_of_elems.clone() {
            header_row.push(Cell::new(&elems));
        }
        elem_table.add_row(Row::new(header_row));
        for (i, sub) in subs.iter().enumerate() {
            let mut row = vec![Cell::new(sub)];
            for j in 0..unique_vec_of_elems.len() {
                row.push(Cell::new(&format!("{:.4}", elem_matrix.get((i, j)).unwrap()    )));
            }
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        Ok(())
    }

    // print table of kinetical information 
    pub fn pretty_print_reaction_data(&self) -> Result<(), std::io::Error> {
        if let Some(vec_of_reaction_data) = &self.vec_of_reaction_data {
            let mut table = Table::new();

            // Add the header row
            table.add_row(row![
                "Reaction number",
                "Reaction Type",
                "Equation",
                "Reactants",
                "Kinetics Data"
            ]);

            // Add rows for each ReactionData
            for (i, reaction_data) in vec_of_reaction_data.iter().enumerate() {
                let reaction_type = format!("{:?}", reaction_data.reaction_type);
                let equation = &reaction_data.eq;
                let reactants = match &reaction_data.react {
                    Some(react) => format!("{:?}", react),
                    None => "None".to_string(),
                };
                let kinetics_data = match &reaction_data.data {
                    ReactionKinetics::Elementary(data) => format!("{:?}", data),
                    ReactionKinetics::Falloff(data) => format!("{:?}", data),
                    ReactionKinetics::Pressure(data) => format!("{:?}", data),
                    ReactionKinetics::ThreeBody(data) => format!("{:?}", data),
                };

                table.add_row(Row::new(vec![
                    Cell::new(i.to_string().as_str()),
                    Cell::new(&reaction_type),
                    Cell::new(equation),
                    Cell::new(&reactants),
                    Cell::new(&kinetics_data),
                ]));
            }

            // Print the table
            table.printstd();
        } else {
            warn!("No reaction data available.");
        }

        Ok(())
    }

    pub fn pretty_print_kindata(&self) {
        let _ = self.pretty_print_reaction_data();
        let _ = self.pretty_print_substances_verbose();
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_user_reactions_new() {
        let shortcut_reactions = Some(vec![
            "C_1".to_string(),
            "C_2".to_string(),
            "C_3".to_string(),
        ]);
        let mut user_reactions = KinData::new();
        // feel shortcut_reactions: Vec<String>, get: map_of_reactions: HashMap<String, String>,grouped_map_of_reactions: HashMap<String, Vec<String>>,
        //  map_of_reaction_data: HashMap<String, Value>, vec_of_equations: Vec<String>, substances: Vec<String>, stecheodata: ReactionAnalyze
        user_reactions.shortcut_reactions = shortcut_reactions.clone();
        assert_eq!(
            user_reactions.shortcut_reactions,
            Some(vec![
                "C_1".to_string(),
                "C_2".to_string(),
                "C_3".to_string()
            ])
        );
        let mut map_of_reactions = HashMap::new();
        map_of_reactions.insert(
            "Cantera".to_string(),
            vec!["1".to_string(), "2".to_string(), "3".to_string()]
                .into_iter()
                .collect(),
        );

        user_reactions.get_reactions_from_shortcuts();
        println!(
            "map_of_reactions: {:?} \n \n ",
            &user_reactions.map_of_reactions
        );
        assert_eq!(user_reactions.map_of_reactions, Some(map_of_reactions));

        user_reactions.reactdata_parsing();
        println!(
            "map_of_reactions data: {:?} \n \n ",
            &user_reactions.clone().vec_of_reaction_data
        );
        assert_eq!(
            user_reactions
                .clone()
                .vec_of_reaction_data
                .unwrap()
                .is_empty(),
            false
        );
        assert_eq!(user_reactions.vec_of_equations.is_empty(), false);

        user_reactions.analyze_reactions();
        let subs = user_reactions.substances;
        assert_eq!(subs.is_empty(), false);
        let S = user_reactions.stecheodata.stecheo_matrx;
        let nunber_of_reactions = S.len();
        let number_of_substances = S[0].len();
        assert_eq!(S.is_empty(), false);
        assert_eq!(nunber_of_reactions, shortcut_reactions.unwrap().len());
        assert_eq!(number_of_substances, subs.len());
    }

    #[test]
    fn test_generate_strings_with_single_character_prefix() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("A1..5".to_string());
        let expected = vec!["A_1", "A_2", "A_3", "A_4", "A_5"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_generate_strings_with_multi_character_prefix() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("Cat5..8".to_string());
        let expected = vec!["Cat_5", "Cat_6", "Cat_7", "Cat_8"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_generate_strings_with_large_numbers() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("Meow20..25".to_string());
        let expected = vec!["Meow_20", "Meow_21", "Meow_22", "Meow_23", "Meow_24", "Meow_25"];
        assert_eq!(result, expected);
    }
}
