//! # Kinetics Library API Module
//!
//! ## Purpose
//! This module provides a comprehensive API for accessing and searching chemical reaction databases.
//! It enables users to load kinetic libraries, search for reactions by various criteria, and retrieve
//! detailed reaction data including Arrhenius parameters and reaction mechanisms.
//!
//! ## Main Data Structures
//! - `KineticData`: Core structure that manages reaction libraries and provides search functionality
//!   - `LibKineticData`: HashMap storing reaction ID -> reaction data mappings
//!   - `HashMapOfReactantsAndProducts`: Nested HashMap for substance-based searches
//!   - `AllLibraries`: Vector of available reaction library names
//!   - Search result vectors for found reactions by products/reagents
//!
//! ## Key Logic Implementation
//! 1. **Library Loading**: Reads JSON files containing reaction databases (Reactbase.json, dict_reaction.json)
//! 2. **Multi-criteria Search**: Supports searching by reaction ID, equation, substances, or field values
//! 3. **Substance Matching**: Uses HashSet operations for efficient subset matching of reagents/products
//! 4. **Data Retrieval**: Returns structured reaction data as serde_json::Value objects
//!
//! ## Usage Pattern
//! ```rust
//! use KiThe::Kinetics::kinetics_lib_api::KineticData;
//! let mut kin_instance = KineticData::new();
//! kin_instance.open_json_files("NUIG");  // Load library
//! kin_instance.print_all_reactions();    // Get all reactions
//! kin_instance.search_reaction_by_reagents_and_products(vec!["CO".to_string()]);
//! println!("Found reactions: {:?}", kin_instance.FoundReactionsByProducts);
//! ```
//!
//! ## Interesting Features
//! - **Dual JSON Structure**: Uses separate files for kinetic parameters and substance mappings
//! - **Flexible Search**: Supports both exact matches and subset searches for reaction participants
//! - **Library Agnostic**: Can work with multiple reaction databases (NUIG, Cantera, Aramco, etc.)
//! - **Efficient Indexing**: Pre-builds equation-to-ID mappings for fast lookups

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;

use serde_json::Value;
#[derive(Debug, Clone)]
pub struct KineticData {
    pub HashMapOfReactantsAndProducts: HashMap<String, HashMap<String, HashSet<String>>>, // {'reaction ID':{'reagent'/"product": HashSet[substance]}}
    pub lib_name: String,
    pub LibKineticData: HashMap<String, Value>, // {'reaction ID':{data of reaction}}
    pub AllLibraries: Vec<String>,              // all reaction libraries
    pub AllEquations: Vec<String>,              // all reaction equations in library
    pub EquationReactionMap: HashMap<String, String>, // {'equation':'reaction ID'}
    pub UserEquationReactionMap: HashMap<String, String>,
    pub FoundReactionsByProducts: Vec<String>,
    pub FoundReactionsByReagents: Vec<String>,
    pub FoundReactionDatasByIDs: Vec<Value>,
}
impl KineticData {
    pub fn new() -> Self {
        Self {
            HashMapOfReactantsAndProducts: HashMap::new(),
            lib_name: String::new(),
            LibKineticData: HashMap::new(),
            AllLibraries: Vec::new(),
            AllEquations: Vec::new(),
            EquationReactionMap: HashMap::new(),
            UserEquationReactionMap: HashMap::new(),
            FoundReactionsByProducts: Vec::new(),
            FoundReactionsByReagents: Vec::new(),
            FoundReactionDatasByIDs: Vec::new(),
        }
    }
    /// Load the reaction library from a JSON file by the name of the reaction library
    pub fn open_json_files(&mut self, big_mech: &str) -> () {
        self.lib_name = big_mech.to_owned();
        let mut file = File::open("Reactbase.json").unwrap();
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).unwrap();
        // библиотека:{ номер реакции :{ данные реакции }}
        let reactlibrary: HashMap<String, HashMap<String, Value>> =
            serde_json::from_str::<HashMap<String, HashMap<String, _>>>(&file_contents).unwrap();
        self.AllLibraries = reactlibrary.keys().map(|k| k.to_string()).collect();
        let library_of_kinetic_parameters = reactlibrary.get(big_mech).unwrap();
        self.LibKineticData = library_of_kinetic_parameters.to_owned();
        //
        let mut file = File::open("dict_reaction.json").unwrap();
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).unwrap();
        // библиотека:{ номер реакции :{ "reagents"/"products":{}: }}
        let reaction_db: HashMap<String, HashMap<String, HashMap<String, HashSet<String>>>> =
            serde_json::from_str(&file_contents).unwrap();
        let library_of_reagents_and_products: &HashMap<String, HashMap<String, HashSet<String>>> =
            reaction_db.get(big_mech).unwrap();
        self.HashMapOfReactantsAndProducts = library_of_reagents_and_products.to_owned();
    }

    /// returns vector of all reaction equations and HashMap {reaction equation : reaction ID}
    pub fn print_all_reactions(&mut self) -> () {
        let all_reactions: Vec<String> =
            self.LibKineticData.keys().map(|k| k.to_string()).collect();
        let all_equations: Vec<String> = self
            .LibKineticData
            .keys()
            .map(|k| {
                self.LibKineticData
                    .get(k)
                    .unwrap()
                    .get("eq")
                    .unwrap()
                    .as_str()
                    .unwrap()
                    .to_string()
            })
            .collect();
        self.AllEquations = all_equations.clone();
        let EquationReactionIDMap: HashMap<String, String> = all_equations
            .iter()
            .zip(all_reactions.iter())
            .map(|(eq, r)| (eq.to_string(), r.to_string()))
            .collect();
        self.EquationReactionMap = EquationReactionIDMap;
    } //end print_all_reactions
    /// returns reaction ID and reaction data (parsed from json) for given reaction equation
    pub fn search_reaction_by_equation(&mut self, equation: &str) -> (String, Value) {
        let reaction_id = self.EquationReactionMap.get(equation).unwrap();
        let reaction = self.LibKineticData.get(reaction_id).unwrap();
        return (reaction_id.clone(), reaction.clone());
    }

    // /returns reaction ID and reaction data (parsed from json) for given reagents/products
    fn search_reaction_by_substances(
        &mut self,
        substances: Vec<String>,
        SubstanceType: &str,
    ) -> Vec<String> {
        let substances: HashSet<String> = substances.iter().cloned().collect();
        let mut found_reactions: Vec<String> = Vec::new();
        let db_object = &self.HashMapOfReactantsAndProducts;
        let db_object_mut = &mut db_object.clone();
        let all_reaction_ID: Vec<&str> = self.LibKineticData.keys().map(|k| k.as_str()).collect();
        for r_id in all_reaction_ID {
            let base_substances: &mut HashSet<String> = db_object_mut
                .get_mut(r_id)
                .expect("REASON")
                .get_mut(SubstanceType)
                .expect("REASON");
            if substances.is_subset(&base_substances) {
                found_reactions.push(r_id.to_string());
            }
        }
        return found_reactions;
    }
    /// search reactions where these substances are reagents/products
    pub fn search_reaction_by_reagents_and_products(&mut self, reagents: Vec<String>) -> () {
        let found_reactions_by_reagents =
            self.search_reaction_by_substances(reagents.clone(), "reagents");
        self.FoundReactionsByReagents = found_reactions_by_reagents.clone();
        let found_reactions_by_products = self.search_reaction_by_substances(reagents, "products");
        self.FoundReactionsByProducts = found_reactions_by_products.clone();
    }
    /// get concrete reaction data
    pub fn search_reactdata_by_reaction_id(&mut self, reaction_id: &str) -> Value {
        let reaction = self.LibKineticData.get(reaction_id).unwrap();
        return reaction.clone();
    }
    /// get concrete reaction data for vector of IDs
    pub fn search_reactdata_for_vector_of_IDs(&mut self, reaction_ids: Vec<String>) -> Vec<Value> {
        let mut vec_of_reactions: Vec<Value> = Vec::new();
        for r_id in reaction_ids {
            let reaction = self.LibKineticData.get(&r_id).unwrap();
            vec_of_reactions.push(reaction.clone());
        }
        self.FoundReactionDatasByIDs = vec_of_reactions.clone();
        return vec_of_reactions;
    }
    //TODO! test
    /// get reaction data for a certain value of json field
    pub fn get_reactions_by_field(&self, field: &str, field_value: &str) -> Vec<Value> {
        self.LibKineticData
            .values()
            .filter(|reaction| {
                reaction
                    .get(field)
                    .and_then(|v| v.as_str())
                    .map_or(false, |v| v == field_value)
            })
            .cloned()
            .collect()
    }
    #[allow(dead_code)]
    /// form user chosen data
    pub fn form_user_chosen_data(
        &self,
        what_to_get: SearchVariants,
        mech_name: Option<String>,
    ) -> HashMap<String, HashMap<String, Value>> {
        let big_mech = if let Some(mech_name) = { mech_name } {
            mech_name
        } else {
            self.lib_name.clone()
        };
        let reaction_ids = self.get_search_results(what_to_get);
        if reaction_ids.len() == 0 {
            return HashMap::new();
        }
        let mut hash_map = HashMap::new();
        let mut hash_map_id_data = HashMap::new();
        for r_id in reaction_ids {
            let reaction = self.LibKineticData.get(&r_id).unwrap();
            hash_map_id_data.insert(r_id.clone(), reaction.clone());
        }
        hash_map.insert(big_mech.clone(), hash_map_id_data.clone());
        return hash_map;
    }
    ///
    #[allow(dead_code)]
    pub fn get_search_results(&self, what_to_get: SearchVariants) -> Vec<String> {
        match what_to_get {
            SearchVariants::FoundReactionsByProducts => self.FoundReactionsByProducts.clone(),
            SearchVariants::FoundReactionsByReagents => self.FoundReactionsByReagents.clone(),
        }
    }
}

pub enum SearchVariants {
    FoundReactionsByProducts,
    FoundReactionsByReagents,
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TESTS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kinetic_data_new() {
        let mut kin_instance = KineticData::new();
        let lib = "NUIG";
        let LibVec: HashSet<String> = vec![
            "Beckstead".to_string(),
            "NUIG".to_string(),
            "Cantera".to_string(),
            "Beckstead_c".to_string(),
            "Aramco".to_string(),
        ]
        .into_iter()
        .collect();

        // collecting reaction data for library name lib
        kin_instance.open_json_files(lib);
        // veiew all reactions in library
        kin_instance.print_all_reactions();
        assert_eq!(kin_instance.HashMapOfReactantsAndProducts.is_empty(), false);
        assert_eq!(kin_instance.LibKineticData.is_empty(), false);
        let all_lib: HashSet<String> = kin_instance.AllLibraries.clone().into_iter().collect();
        assert_eq!(all_lib, LibVec);
        assert_eq!(kin_instance.AllEquations.is_empty(), false);

        // search reaction by equation
        let equation = kin_instance
            .LibKineticData
            .get("1")
            .unwrap()
            .get("eq")
            .unwrap()
            .as_str()
            .unwrap()
            .to_string();
        // returns reaction ID and reaction data (parsed from json) for given reaction equation
        let (reaction_id, reaction) = kin_instance.search_reaction_by_equation(&equation);
        assert_eq!(reaction_id, "1".to_string());
        assert_eq!(reaction.get("eq").unwrap().as_str().unwrap(), equation);
        let substances: Vec<String> = vec!["CO".to_string(), "O".to_string()];
        // search reactions by substances
        kin_instance.search_reaction_by_reagents_and_products(substances);
        assert_eq!(kin_instance.FoundReactionsByProducts.is_empty(), false);
        assert_eq!(kin_instance.FoundReactionsByReagents.is_empty(), false);
        kin_instance.search_reactdata_for_vector_of_IDs(vec!["1".to_string(), "2".to_string()]);
        assert_eq!(kin_instance.FoundReactionDatasByIDs.is_empty(), false);
    }

    #[test]
    fn test_search_value() {
        let mut kin_instance = KineticData::new();
        let lib = "NUIG";

        // collecting reaction data for library name lib
        kin_instance.open_json_files(lib);
        // veiew all reactions in library
        kin_instance.print_all_reactions();

        // search reaction by equation
        let data_of_reacion = kin_instance.LibKineticData.get("1").unwrap();
        println!(" data_of_reacion {:?}", data_of_reacion);
    }

    #[test]
    fn test_search_data_by_field() {
        let mut kin_instance = KineticData::new();
        kin_instance.open_json_files("NUIG");

        let reactions_by_field = kin_instance.get_reactions_by_field("type", "falloff");
        println!(" reactions_by_field {:?}", reactions_by_field);
        assert_eq!(reactions_by_field.len(), 109);
    }

    #[test]
    fn test_form_user_chosen_data() {
        let mut kin_instance = KineticData::new();
        kin_instance.open_json_files("NUIG");

        // Prepare some test data
        kin_instance.FoundReactionsByProducts = vec!["1".to_string(), "2".to_string()];
        kin_instance.FoundReactionsByReagents = vec!["3".to_string(), "4".to_string()];

        // Test with FoundReactionsByProducts and default mechanism name
        let result =
            kin_instance.form_user_chosen_data(SearchVariants::FoundReactionsByProducts, None);
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("NUIG"));
        let nuig_data = result.get("NUIG").unwrap();
        assert_eq!(nuig_data.len(), 2);
        assert!(nuig_data.contains_key("1"));
        assert!(nuig_data.contains_key("2"));

        // Test with FoundReactionsByReagents and specified mechanism name
        let result = kin_instance.form_user_chosen_data(
            SearchVariants::FoundReactionsByReagents,
            Some("CustomMech".to_string()),
        );
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("CustomMech"));
        let custom_mech_data = result.get("CustomMech").unwrap();
        assert_eq!(custom_mech_data.len(), 2);
        assert!(custom_mech_data.contains_key("3"));
        assert!(custom_mech_data.contains_key("4"));

        // Test with empty results
        kin_instance.FoundReactionsByProducts = vec![];
        let result =
            kin_instance.form_user_chosen_data(SearchVariants::FoundReactionsByProducts, None);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_form_user_chosen_data_integration() {
        let mut kin_instance = KineticData::new();
        kin_instance.open_json_files("NUIG");

        // Perform a search to populate FoundReactionsByProducts
        kin_instance
            .search_reaction_by_reagents_and_products(vec!["O".to_string(), "H".to_string()]);

        // Test with actual search results
        let result =
            kin_instance.form_user_chosen_data(SearchVariants::FoundReactionsByProducts, None);
        assert!(!result.is_empty());
        assert!(result.contains_key("NUIG"));
        let nuig_data = result.get("NUIG").unwrap();
        assert!(!nuig_data.is_empty());

        // Verify that the returned data matches the LibKineticData
        for (reaction_id, reaction_data) in nuig_data {
            assert_eq!(
                reaction_data,
                kin_instance.LibKineticData.get(reaction_id).unwrap()
            );
        }
    }
    /*
    #[test]
    fn test_kinetic_data_open_json_files() {
        let mut kinetic_data = KineticData::new();

        // Replace "big_mech" with the actual file name
        kinetic_data.open_json_files("big_mech");

        // Add assertions to verify the content of kinetic_data.HashMapOfReactantsAndProducts,
        // kinetic_data.LibKineticData, kinetic_data.AllLibraries, etc.
    }



    */
    // Add more tests for other methods as needed
}
