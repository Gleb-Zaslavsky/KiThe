//mod mechfinder_api;
use crate::Kinetics::kinetics_lib_api::KineticData;
use crate::Kinetics::mechfinder_api::{ parse_kinetic_data_vec, ReactionData, Mechanism_search};
use crate::Kinetics::parsetask::{decipher_vector_of_shortcuts, decipher_vector_of_shortcuts_to_pairs};
use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
use std::fs::File;
use std::io::Write;
use serde_json::json;
use serde_json::Value;
use std::collections::HashMap;
use prettytable::{Table, Row, Cell};
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
            return vec![]; // Return empty if format is incorrect
        }

        let prefix = parts[0]
            .chars()
            .take_while(|c| c.is_alphabetic())
            .collect::<String>(); //extracts the prefix (letters) and the numeric portion from the start.
        let start: usize = parts[0][prefix.len()..].parse().unwrap_or(0); // then parses the numeric values and generates the range of strings.
        let end: usize = parts[1].parse().unwrap_or(0);

        (start..=end).map(|i| format!("{}{}", prefix, i)).collect() //formatted strings are collected into a Vec<String>.
    }

    /// decifer vector of shortcuts to full reaction names and store them in map {'library':"id of reaction"} and Vec<("library, reaction_id")>
    pub fn set_reactions_from_shortcuts(&mut self) -> () {
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
            println!("KinData::create_map_of_reactions: shortcut_reactions is None");
        }
    }
    ///construct reaction mechanism
    pub fn construct_mechanism(&mut self,    task_substances: Vec<String>, task_library: String,) {

        let found_mech = Mechanism_search::new(   task_substances,task_library.clone());
      //  self.vec_of_equations = found_mech.vec_of_reactions;
        self.vec_of_reaction_data = Some(found_mech.reactdata);
        let reactions = found_mech.mechanism;
        let mut full_addres = Vec::new();
        for reaction in reactions.iter() {
            let addres = format!("{}_{}", task_library, reaction);
            full_addres.push(addres);
        }
        self.shortcut_reactions = Some(full_addres);

    }
    /////////////////////////////////GETTING REACTION DATA///////////////////////////////////////////
    /// parse reaction libraries and extracts data into structs
    pub fn reactdata_parsing(&mut self) -> () {
        if let Some(vec_of_reaction_values ) = & self.vec_of_reaction_Values  {

            let (vec_ReactionData, vec_of_equations) =
                    parse_kinetic_data_vec( vec_of_reaction_values.clone());
          //  println!("map_of_reaction_data: {:?}", &vec_of_reaction_data);
         //   println!("vec_of_equations: {:?}", &vec_of_equations);
            self.vec_of_reaction_data = Some(vec_ReactionData);
            self.vec_of_equations = vec_of_equations;
        } else {
            println!("KinData::reactdata_from_shortcuts: map_of_reactions is None");
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
        StoichAnalyzer_instance.search_substances();
        self.substances = StoichAnalyzer_instance.substances.clone();
        //find stoichiometric matrix and other matrices
        StoichAnalyzer_instance.analyse_reactions();
        StoichAnalyzer_instance.create_matrix_of_elements();
        self.stecheodata = StoichAnalyzer_instance;
    }
    ///////////////////////////INPUT/OUTPUT/////////////////////////////////////////////////////////
    /// printlns the chosen reactions 
    pub fn print_raw_reactions(&self) -> Result<(), std::io::Error> {
        if let Some(vec_of_reaction_Values) = & self.vec_of_reaction_Values {
             // Convert the vector of Values to a JSON array
            let json_array = json!(vec_of_reaction_Values);
            // Write the JSON array to a file
            let mut file = File::create("raw_reactions.json")?;
            file.write_all(serde_json::to_string_pretty(&json_array)?.as_bytes())?;
            println!("Raw reactions have been written to raw_reactions.json");
            Ok(())
    
        } else {
            println!("KinData::reactdata_from_shortcuts: map_of_reactions is None");
            Ok(())
        }
    }// print raw_reactions

    pub fn pretty_print_kindata(&self) -> Result<(),  std::io::Error> {
        if let Some(vec_of_reaction_Values) = &self.vec_of_reaction_Values {
            // Create a new table
            let mut table = Table::new();
    
            // Assuming each Value is an object with the same keys
            if let Some(Value::Object(first_obj)) = vec_of_reaction_Values.get(0) {
                // Add the header row
                let header: Vec<Cell> = first_obj.keys().map(|k| Cell::new(k)).collect();
                table.add_row(Row::new(header));
            }
    
            // Add rows for each Value
            for value in vec_of_reaction_Values {
                if let Value::Object(obj) = value {
                    let row: Vec<Cell> = obj.values().map(|v| Cell::new(&v.to_string())).collect();
                    table.add_row(Row::new(row));
                }
            }
    
            // Print the table to stdout
            table.printstd();
    
            println!("Raw reactions have been written to raw_reactions.json");
            Ok(())
        } else {
            println!("KinData::reactdata_from_shortcuts: map_of_reactions is None");
            Ok(())
        }
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

        user_reactions.set_reactions_from_shortcuts();
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
        let expected = vec!["A1", "A2", "A3", "A4", "A5"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_generate_strings_with_multi_character_prefix() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("Cat5..8".to_string());
        let expected = vec!["Cat5", "Cat6", "Cat7", "Cat8"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_generate_strings_with_large_numbers() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("Meow20..25".to_string());
        let expected = vec!["Meow20", "Meow21", "Meow22", "Meow23", "Meow24", "Meow25"];
        assert_eq!(result, expected);
    }
}
