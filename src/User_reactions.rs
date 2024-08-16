//mod mechfinder_api;
use crate::mechfinder_api::parse_kinetic_data;
use crate::parsetask::decipher_vector_of_shortcuts;
use crate::kinetics_lib_api::KineticData;
use crate::reaction_analyzer::ReactionAnalyzer;

use std::collections::HashMap;
use serde_json::{Value, Map, Number};
/// 0.1.2
/// processing of user-chosen reactions.
///  So you define what reactions you need by using the constructor of mechanism from the mechfinder_api module or
/// manually with help of kinetics_lib_api module. Now you can 


//structure to store user task and reaction data
#[derive(Debug, Clone)]
pub struct UserReactions {
    pub shortcut_reactions: Vec<String>, // vector of reaction shortcut names
    pub  map_of_reactions: HashMap<String, Vec<String>>, // full "address" of reaction {'library':"id of reaction"} group reactions by library names {'library':[reaction_ids in that library]}
   
    pub map_of_reaction_data: HashMap<String, Value>, // data of all reactions
    pub vec_of_equations: Vec<String>, // vector of equations of reactions
    pub substances: Vec<String>, // vector of substance names
    pub stecheodata: ReactionAnalyzer, // matrix of stoichiometric coefficients and other matrices

}

impl UserReactions {
    pub fn new() -> Self {
        Self {
            shortcut_reactions: Vec::new(),
            map_of_reactions: HashMap::new(),
            map_of_reaction_data: HashMap::new(),
            vec_of_equations: Vec::new(),
            substances: Vec::new(), 
            stecheodata:ReactionAnalyzer::new(),
        }
    }
    // decifer vector of shortcuts to full reaction names and store them in map {'library':"id of reaction"}
    pub fn create_map_of_reactions(&mut self) -> () {
        let vec:Vec<&str> = self.shortcut_reactions.iter().map(|s| s.as_str()).collect();
        self.map_of_reactions = decipher_vector_of_shortcuts(vec);

    }

    

    pub fn create_map_of_reaction_data(&mut self) -> () {
        let mut vec_of_equations = Vec::new();
        let mut ReactionDataHash = HashMap::new();
        // instance of KineticData with opened library json files
        let mut kin_instance = KineticData::new();
        for (lib, reaction_id_vector) in self.map_of_reactions.iter() {
            // collecting reaction data for each library name
            kin_instance.open_json_files(lib);
            let mut vec_of_reactions_value: Vec<Value> = Vec::new();
            let mut vec_of_reactions: Vec<String> = Vec::new();
            for reaction_id in reaction_id_vector.iter() {
                let reaction_data_Value = kin_instance.search_reactdata_by_reaction_id( &reaction_id);
                vec_of_reactions_value.push(reaction_data_Value);
                vec_of_reactions.push(reaction_id.to_string());
            }// for
            // now we have a vector of json objects with reaction data and a vector of reaction ids
            // lets parse it to map of structures
            let (mut ReactionDataHash_for_lib, vec_of_equations_for_lib) = parse_kinetic_data( lib, &vec_of_reactions, vec_of_reactions_value);  
            vec_of_equations.extend(vec_of_equations_for_lib);
            ReactionDataHash.extend(ReactionDataHash_for_lib);

        }
        println!("map_of_reaction_data: {:?}", &ReactionDataHash);
        println!("vec_of_equations: {:?}", &vec_of_equations);
        self.map_of_reaction_data = ReactionDataHash;
        self.vec_of_equations = vec_of_equations;

    }
    pub fn analyze_reactions(&mut self) -> () {
        // iniciate instance of ReactionAnalyzer 
        let mut ReactionAnalyzer_instance =ReactionAnalyzer::new();
        // copy vector of reactions to ReactionAnalyzer_instance
        ReactionAnalyzer_instance.reactions = self.vec_of_equations.clone();
        // parse to find substance names
        ReactionAnalyzer_instance.search_substances();
        self .substances = ReactionAnalyzer_instance.substances.clone();
        //find stoichiometric matrix and other matrices
        ReactionAnalyzer_instance.analyse_reactions();

        self.stecheodata = ReactionAnalyzer_instance;

        
    }
    
}
 

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_user_reactions_new() {
        let shortcut_reactions = vec!["C_1".to_string(), "C_2".to_string(), "C_3".to_string()];
        let mut user_reactions = UserReactions::new();
        // feel shortcut_reactions: Vec<String>, get: map_of_reactions: HashMap<String, String>,grouped_map_of_reactions: HashMap<String, Vec<String>>,
       //  map_of_reaction_data: HashMap<String, Value>, vec_of_equations: Vec<String>, substances: Vec<String>, stecheodata: ReactionAnalyze
       user_reactions.shortcut_reactions = shortcut_reactions.clone();
        assert_eq!(user_reactions.shortcut_reactions, vec!["C_1".to_string(), "C_2".to_string(), "C_3".to_string()]   );
        let mut map_of_reactions = HashMap::new();
        map_of_reactions.insert("Cantera".to_string(),  vec! ["1".to_string(), "2".to_string(), "3".to_string()].into_iter().collect());    

     
        user_reactions.create_map_of_reactions();
        assert_eq!(user_reactions.map_of_reactions, map_of_reactions);
        
        user_reactions.create_map_of_reaction_data();
        assert_eq!(user_reactions.map_of_reaction_data.is_empty(), false); 
        assert_eq!(user_reactions.vec_of_equations.is_empty(), false);

        user_reactions.analyze_reactions();
        let subs = user_reactions.substances;
        assert_eq!(subs.is_empty(), false);
        let S = user_reactions.stecheodata.stecheo_matrx;
        let nunber_of_reactions = S.len();
        let number_of_substances = S[0].len();
        assert_eq!(S.is_empty(), false);
        assert_eq!( nunber_of_reactions,  shortcut_reactions.len());
        assert_eq!( number_of_substances,  subs.len());

        /* 
      //  assert_eq!(&stecheodata.is_empty(), &false);
        //assert_eq!(user_reactions.stecheodata. is_empty(), false);
       
       
        assert_eq!(user_reactions.map_of_reactions, HashMap::new());
        assert_eq!(user_reactions.grouped_map_of_reactions, HashMap::new());
        assert_eq!(user_reactions.map_of_reaction_data, HashMap::new());
        assert_eq!(user_reactions.vec_of_equations, Vec::new());
        assert_eq!(user_reactions.substances, Vec::new());
        assert_eq!(user_reactions.stecheodata, ReactionAnalyzer::new());
        */
    }


    
}
    
    