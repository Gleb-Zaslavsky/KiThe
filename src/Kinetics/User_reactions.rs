//mod mechfinder_api;
use crate::Kinetics::kinetics_lib_api::KineticData;
use crate::Kinetics::mechfinder_api::{
    Mechanism_search, ReactionData, ReactionKinetics, parse_kinetic_data_vec,
};
use crate::Kinetics::parsetask::{
    decipher_vector_of_shortcuts, decipher_vector_of_shortcuts_to_pairs,
};
use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
use prettytable::{Cell, Row, Table, row};
use serde_json::json;
use serde_json::{Value, from_reader, to_writer_pretty};
use std::collections::HashMap;

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;

use log::{info, warn};
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
    pub shortcut_reactions: Option<Vec<String>>,
    /// vector of reaction shortcut names
    pub map_of_reactions: Option<HashMap<String, Vec<String>>>,
    /// full "address" of reaction {'library':"id of reaction"} group reactions by library names {'library':[reaction_ids in that library]}
    pub vec_of_pairs: Option<Vec<(String, String)>>,
    /// vector of pairs of reaction library names and reaction ids
    pub vec_of_reaction_Values: Option<Vec<Value>>,
    /// vector of reaction data in the form of serde Values
    pub vec_of_reaction_data: Option<Vec<ReactionData>>, // data of all reactions
    pub vec_of_equations: Vec<String>,
    /// vector of equations of reactions
    pub substances: Vec<String>,
    /// vector of substance names
    pub groups: Option<HashMap<String, HashMap<String, usize>>>,
    /// Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
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
        let end: usize = if let Some(last_char) = parts[1].chars().last() {
            //  checks if there's a last character in parts[1]
            if last_char.is_numeric() {
                // the last character is numerit
                parts[1]
                    .chars()
                    .rev() // Reverses the string
                    .take_while(|c| c.is_numeric()) // Takes characters while they are numeric
                    .collect::<String>() //Collects these into a string
                    .chars()
                    .rev()
                    .collect::<String>() // Reverses this string again to get the number in the correct order.
                    .parse()
                    .unwrap_or(0) //Parses this string into a usize.
            } else {
                panic!("last symbol must be numeric: {}", shortcut_range);
            }
        } else {
            panic!("last symbol must be numeric: {}", shortcut_range);
        };
        let res: Vec<String> = (start..=end).map(|i| format!("{}_{}", prefix, i)).collect(); //formatted strings are collected into a Vec<String>.
        info!("task includes reactions as follows: {:#?}", &res);
        self.shortcut_reactions = Some(res.clone());
        return res;
    }

    /// decifer vector of shortcuts to full reaction names and store them in map {'library':"id of reaction"} and Vec<("library, reaction_id")>
    /// then open kinetic libraries, get reaction info corresponding to the guven shortcuts and store it in structure
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
    pub fn construct_mechanism(&mut self, task_substances: Vec<String>, task_library: String) {
        let mut found_mech = Mechanism_search::new(task_substances, task_library.clone());
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
        self.vec_of_pairs = Some(
            reactions
                .clone()
                .iter()
                .map(|s| (task_library.clone(), s.to_owned()))
                .collect(),
        );
        self.shortcut_reactions = Some(full_addres);
        self.substances = found_mech.reactants;
    }
    /////////////////////////////////REACTIONS MANIPULATIONS//////////////////////////////////////
    /// add manually reaction data as serde Value
    pub fn append_reaction(&mut self, reactions: Vec<Value>) {
        let mut old_reactions = self.vec_of_reaction_Values.as_mut().unwrap().clone();
        old_reactions.extend(reactions);
        self.vec_of_reaction_Values = Some(old_reactions);
    }
    /// add manually reaction data with there shortcut names
    pub fn append_reaction_with_shortcut(&mut self, reactions: Vec<Value>, shortcuts: Vec<String>) {
        assert_eq!(reactions.len(), shortcuts.len());
        self.append_reaction(reactions);
        let mut old_shortcuts = self.shortcut_reactions.clone().unwrap();
        old_shortcuts.extend(shortcuts);
        self.shortcut_reactions = Some(old_shortcuts);
    }
    ///remove manually reaction data by its index
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
    pub fn remove_reaction_by_eq(&mut self, reaction_name: &str) {
        let i = &self
            .vec_of_equations
            .clone()
            .iter()
            .position(|eq_i| *eq_i == reaction_name);
        if let Some(index) = i {
            self.remove_by_index(*index);
        }
    }

    /////////////////////////////////COMPUTING AND PARSING REACTION DATA///////////////////////////////////////////
    /// parse serde Values with kinetic info into structs and store it in KinData structure
    pub fn reactdata_parsing(&mut self) -> () {
        if let Some(vec_of_reaction_values) = &self.vec_of_reaction_Values {
            let (vec_ReactionData, vec_of_equations) =
                parse_kinetic_data_vec(vec_of_reaction_values.clone());
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
        if self.substances.is_empty() | (self.substances.len() == 0) {
            // no need to parse for substances if they are already known
            StoichAnalyzer_instance.search_substances();
            self.substances = StoichAnalyzer_instance.substances.clone();
        } else {
            StoichAnalyzer_instance.substances = self.substances.clone()
        };
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
    /// save the chosen reaction kinetic data to json
    pub fn save_raw_reactions(&self, file_name: &str) -> Result<(), std::io::Error> {
        if let Some(vec_of_reaction_Values) = &self.vec_of_reaction_Values {
            // Convert the vector of Values to a JSON array
            let json_array = json!(vec_of_reaction_Values);
            // Write the JSON array to a file
            let name = format!("{} {}", file_name, "{}.json");
            let mut file = File::create(name)?;
            file.write_all(serde_json::to_string_pretty(&json_array)?.as_bytes())?;
            info!("Raw reactions have been written to raw_reactions.json");
            Ok(())
        } else {
            warn!("KinData::reactdata_from_shortcuts: map_of_reactions is None");
            Ok(())
        }
    }
    /// save the chosen reactions kinetic data with their shortcuts to json as dictionaries {"shortcut":"kinetic data"}
    pub fn save_reactions_with_shortcuts(&self, file_name: &str) -> Result<(), std::io::Error> {
        if let Some(reaction_values) = &self.vec_of_reaction_Values {
            if let Some(shortcuts) = &self.shortcut_reactions {
                if shortcuts.len() != reaction_values.len() {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "Vectors must have the same length",
                    ));
                }

                let map: HashMap<String, Value> = shortcuts
                    .iter()
                    .cloned()
                    .zip(reaction_values.iter().cloned())
                    .collect();

                let file = File::create(file_name)?;
                to_writer_pretty(file, &map)?;
            } else {
                warn!("KinData: no vector of shortcuts");
            }
            Ok(())
        } else {
            warn!("KinData: no vector of Values");
            Ok(())
        }
    }
    ///  algorithm: 1) takes the name of the file being created as an argument.
    ///  2) checks if there is a "KINETICS" or "REACTIONS" header there; if not, it adds them, and under this header, it adds records to a
    /// form available for use by serde (i.e., of the json type): "library_name":{"reaction_name":{ ...reaction data

    pub fn create_kinetics_document(
        &self,
        file_name: &str,
    ) -> std::io::Result<HashMap<String, HashMap<String, Value>>> {
        let path = Path::new(file_name);
        let file_exists = path.exists();

        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .append(true)
            .open(file_name)?;

        if !file_exists {
            writeln!(file, "KINETICS")?;
        } else {
            let reader = BufReader::new(File::open(file_name)?);
            let has_header = reader.lines().any(|line| {
                line.as_ref().unwrap().clone().trim().to_uppercase() == "KINETICS"
                    || line.as_ref().unwrap().clone().trim().to_uppercase() == "REACTIONS"
            });

            if !has_header {
                writeln!(file, "\nKINETICS")?;
            }
        }
        let mut hashmap_to_save: HashMap<String, HashMap<String, Value>> = HashMap::new();
        // if got info of pairs - use pairs
        if let Some(pairs) = self.vec_of_pairs.clone() {
            for (i, (library_name, reaction_name)) in pairs.iter().enumerate() {
                let reactdata = self.vec_of_reaction_Values.as_ref().unwrap()[i].clone();
                if hashmap_to_save.contains_key(library_name) {
                    hashmap_to_save
                        .get_mut(library_name)
                        .unwrap()
                        .insert(reaction_name.clone(), reactdata);
                } else {
                    hashmap_to_save.insert(
                        library_name.clone(),
                        HashMap::from([(reaction_name.clone(), reactdata)]),
                    );
                }
            }
        } else if let Some(map_of_reactions) = self.map_of_reactions.clone() {
            for (i, (library_name, reactions)) in map_of_reactions.iter().enumerate() {
                let reactdata = self.vec_of_reaction_Values.as_ref().unwrap()[i].clone();
                for reaction_id in reactions.iter() {
                    if hashmap_to_save.contains_key(library_name) {
                        hashmap_to_save
                            .get_mut(library_name)
                            .unwrap()
                            .insert(reaction_id.clone(), reactdata.clone());
                    } else {
                        hashmap_to_save.insert(
                            library_name.clone(),
                            HashMap::from([(reaction_id.clone(), reactdata.clone())]),
                        );
                    }
                } // 
            } // for library_name
        } //else if 

        writeln!(file, "{}", serde_json::to_string_pretty(&hashmap_to_save)?)?;

        Ok(hashmap_to_save)
    }
    /// load kinetic data of reactions with their shortcuts from json, where they are stored as dictionaries {"shortcut":"kinetic data"}
    /// and save them into vector of kinetic data serde values and vector of shortcuts
    pub fn load_reactions_from_json(&mut self, file_name: &str) -> Result<(), std::io::Error> {
        let file = File::open(file_name)?;
        let res: Result<HashMap<String, Value>, serde_json::Error> = from_reader(file);
        if let Ok(map) = res {
            let (shortcuts, reaction_values): (Vec<String>, Vec<Value>) = map.into_iter().unzip();

            self.shortcut_reactions = Some(shortcuts);
            self.vec_of_reaction_Values = Some(reaction_values);

            Ok(())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Vectors must have the same length",
            ))
        }
    }

    /// print stecheometric matrix and matrix of elemental composition in the pretty human-readable format
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
                row.push(Cell::new(&format!(
                    "{:.4}",
                    elem_matrix.get((i, j)).unwrap()
                )));
            }
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        Ok(())
    }

    /// print table of kinetical information in pretty format
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
    /// both reaction data and matrces are pretty printed
    pub fn pretty_print_kindata(&self) {
        let _ = self.pretty_print_reaction_data();
        let _ = self.pretty_print_substances_verbose();
    }

    /////////////////////KINETIC CONTANT FUNCTIONS //////////////////////
    ///  Returns the kinetic constant for a reaction with number reaction_id
    pub fn calc_K_const_for_1_react(
        &self,
        reaction_id: usize,
        temp: f64,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
    ) -> Result<f64, std::io::Error> {
        if let Some(reactions) = &self.vec_of_reaction_data {
            let reaction_data = &reactions[reaction_id];
            let K = reaction_data.K_const(temp, pres, concentrations);

            Ok(K)
        } else {
            warn!("no reaction data available");
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "no vector of reaction data available",
            ))
        }
    }
    /// 1) Returns the kinetic constant for all reactions found 2) sort reactions by reaction rate constant and regroup  3)
    /// if  save_rearranged = Some(true) rearrange all data in structure according to reaction rate constant magnitude
    /// 4) pretty print table with headera    "Reaction number",  "Equation",  "log10 of Constant at given T"
    pub fn calc_K_const_for_all_reactions(
        &mut self,
        temp: f64,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
        save_rearranged: Option<bool>,
    ) -> Result<Vec<f64>, std::io::Error> {
        if let Some(reactions) = &self.vec_of_reaction_data {
            let mut vec_of_K_const = Vec::new();
            for reaction in reactions.iter() {
                let K = reaction.K_const(temp, pres, concentrations.clone());
                vec_of_K_const.push(K);
            }
            // Let us sort vec_of_K_const and rearrange vec_of_equations and vec_of_shortcuts in the same order,
            // Create a vector of indices
            let mut indices: Vec<usize> = (0..vec_of_K_const.len()).collect();
            // Sort indices based on K_const values (in descending order)
            indices.sort_by(|&i, &j| {
                vec_of_K_const[j]
                    .partial_cmp(&vec_of_K_const[i])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            // Use the sorted indices to rearrange all vectors
            vec_of_K_const = indices.iter().map(|&i| vec_of_K_const[i]).collect();
            let vec_of_equations = self.rearrange_vectors(indices, save_rearranged);

            // create pretty table
            Self::pretty_print_constants(vec_of_equations, vec_of_K_const.clone());
            Ok(vec_of_K_const)
        } else {
            warn!("no reaction data available");
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "no vector of reaction data available",
            ))
        }
    }
    /// the same but for every reaction for given range of T calculated K values and maximum values compared to those for another reaction
    pub fn calc_K_const_for_all_reactions_forTrange(
        &mut self,
        T0: f64,
        Tend: f64,
        n: usize,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
        save_rearranged: Option<bool>,
    ) -> Result<Vec<f64>, std::io::Error> {
        if let Some(reactions) = &self.vec_of_reaction_data {
            let mut vec_of_K_const = Vec::new();
            for reaction in reactions.iter() {
                let K_vec = reaction.K_const_for_T_range(T0, Tend, n, pres, concentrations.clone());

                let max_K = K_vec.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
                vec_of_K_const.push(*max_K);
            }
            // Let us sort vec_of_K_const and rearrange vec_of_equations and vec_of_shortcuts in the same order,
            // Create a vector of indices
            let mut indices: Vec<usize> = (0..vec_of_K_const.len()).collect();
            // Sort indices based on K_const values (in descending order)
            indices.sort_by(|&i, &j| {
                vec_of_K_const[j]
                    .partial_cmp(&vec_of_K_const[i])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            // Use the sorted indices to rearrange all vectors
            vec_of_K_const = indices.iter().map(|&i| vec_of_K_const[i]).collect();
            let vec_of_equations = self.rearrange_vectors(indices, save_rearranged);
            // println!("reaction {:?}",vec_of_K_const);
            // create pretty table
            Self::pretty_print_constants(vec_of_equations, vec_of_K_const.clone());
            Ok(vec_of_K_const)
        } else {
            warn!("no reaction data available");
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "no vector of reaction data available",
            ))
        }
    }
    /// rearrange all vector data accorfing to new order
    fn rearrange_vectors(
        &mut self,
        indices: Vec<usize>,
        save_rearranged: Option<bool>,
    ) -> Vec<String> {
        // If any of the required Option values are None, return current equations
        if self.vec_of_equations.is_empty() {
            return Vec::new();
        }

        let mut vec_of_equations = self.vec_of_equations.clone();
        if let (Some(shortcuts), Some(pairs), Some(values), Some(data)) = (
            &self.shortcut_reactions,
            &self.vec_of_pairs,
            &self.vec_of_reaction_Values,
            &self.vec_of_reaction_data,
        ) {
            let mut vec_of_shortcuts = shortcuts.clone();
            let mut vec_of_pairs = pairs.clone();
            let mut vec_of_reaction_Values = values.clone();
            let mut vec_of_reaction_data = data.clone();

            // Rearrange all vectors according to indices
            vec_of_equations = indices
                .iter()
                .map(|&i| vec_of_equations[i].clone())
                .collect();
            vec_of_shortcuts = indices
                .iter()
                .map(|&i| vec_of_shortcuts[i].clone())
                .collect();
            vec_of_pairs = indices.iter().map(|&i| vec_of_pairs[i].clone()).collect();
            vec_of_reaction_Values = indices
                .iter()
                .map(|&i| vec_of_reaction_Values[i].clone())
                .collect();
            vec_of_reaction_data = indices
                .iter()
                .map(|&i| vec_of_reaction_data[i].clone())
                .collect();

            // Only update the struct fields if save_rearranged is Some(true)
            if save_rearranged == Some(true) {
                self.vec_of_equations = vec_of_equations.clone();
                self.shortcut_reactions = Some(vec_of_shortcuts);
                self.vec_of_pairs = Some(vec_of_pairs);
                self.vec_of_reaction_Values = Some(vec_of_reaction_Values);
                self.vec_of_reaction_data = Some(vec_of_reaction_data);
            }
        }

        vec_of_equations
    }
    fn pretty_print_constants(vec_of_equations: Vec<String>, vec_of_K_const: Vec<f64>) {
        let mut table = Table::new();
        table.add_row(row![
            "Reaction number",
            "Equation",
            "log10 of Constant at given T"
        ]);
        for (i, K) in vec_of_K_const.iter().enumerate() {
            table.add_row(row![
                i.to_string().as_str(),
                vec_of_equations[i],
                format!("{:.4}", f64::log10(*K))
            ]);
        }
        table.printstd();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::mechfinder_api::ReactionType;
    use crate::Kinetics::mechfinder_api::kinetics::ElementaryStruct;
    use approx::assert_relative_eq;
    use std::fs;
    use tempfile::NamedTempFile;
    const R: f64 = 8.314;
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
        let expected = vec![
            "Meow_20", "Meow_21", "Meow_22", "Meow_23", "Meow_24", "Meow_25",
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_calc_K_const_for_1_react() {
        let mut kin_data = KinData::new();
        // Populate kin_data with some test reaction data
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A + B -> C".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1e13, 0.0, 20000.0],
            }),
        }]);

        let result = kin_data.calc_K_const_for_1_react(0, 1000.0, None, None);
        assert!(result.is_ok());
        let k = result.unwrap();
        let k_expected = 1e13 * f64::exp(-20000.0 / (R * 1000.0));
        assert_relative_eq!(k, k_expected, epsilon = 1e0);
    }

    #[test]
    fn test_calc_K_const_for_all_reactions() {
        let mut kin_data = KinData::new();
        // Populate kin_data with some test reaction data
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A + B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e13, 0.0, 20000.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "C -> D".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e12, 0.0, 20000.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A + B -> C".to_string(), "C -> D".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);

        let result = kin_data.calc_K_const_for_all_reactions(1000.0, None, None, Some(true));
        assert!(result.is_ok());
        let k_values = result.unwrap();
        assert_eq!(k_values.len(), 2);
        let k_expected1 = 1e13 * f64::exp(-20000.0 / (R * 1000.0));
        let k_expected2 = 1e12 * f64::exp(-20000.0 / (R * 1000.0));
        assert_relative_eq!(k_values[0], k_expected1, epsilon = 1e0);
        assert_relative_eq!(k_values[1], k_expected2, epsilon = 1e0);
    }

    #[test]
    fn test_calc_K_const_for_all_reactions_forTrange() {
        let mut kin_data = KinData::new();
        // Populate kin_data with some test reaction data
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A + B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e13, 0.0, 20000.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "C -> D".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e12, 0.0, 20000.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A + B -> C".to_string(), "C -> D".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);

        let result = kin_data.calc_K_const_for_all_reactions_forTrange(
            800.0,
            1200.0,
            5,
            None,
            None,
            Some(true),
        );
        assert!(result.is_ok());
        let max_k_values = result.unwrap();
        assert_eq!(max_k_values.len(), 2);
        let k_expected1 = 1e13 * f64::exp(-20000.0 / (R * 1200.0));
        let k_expected2 = 1e12 * f64::exp(-20000.0 / (R * 1200.0));
        println!(
            "max_k_values {} {}",
            f64::log10(max_k_values[0]),
            f64::log10(max_k_values[1])
        );
        println!(
            "k expected {}, {}",
            f64::log10(k_expected1),
            f64::log10(k_expected2)
        );
        assert_relative_eq!(max_k_values[0], k_expected1, epsilon = 1e0);
        assert_relative_eq!(max_k_values[1], k_expected2, epsilon = 1e0);
    }
    #[test]
    fn test_create_kinetics_document() {
        let mut kin_data = KinData::new();

        // Test with vec_of_pairs
        kin_data.vec_of_pairs = Some(vec![
            ("Library1".to_string(), "Reaction1".to_string()),
            ("Library1".to_string(), "Reaction2".to_string()),
            ("Library2".to_string(), "Reaction3".to_string()),
        ]);
        kin_data.vec_of_reaction_Values = Some(vec![
            json!({"data": "value1"}),
            json!({"data": "value2"}),
            json!({"data": "value3"}),
        ]);

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());

        let content = fs::read_to_string(file_path).unwrap();
        assert!(content.contains("KINETICS"));
        assert!(content.contains("Library1"));
        assert!(content.contains("Library2"));
        assert!(content.contains("Reaction1"));
        assert!(content.contains("Reaction2"));
        assert!(content.contains("Reaction3"));
        assert!(content.contains("value1"));
        assert!(content.contains("value2"));
        assert!(content.contains("value3"));

        // Test with map_of_reactions
        kin_data.vec_of_pairs = None;
        kin_data.map_of_reactions = Some(HashMap::from([
            (
                "Library3".to_string(),
                vec!["Reaction4".to_string(), "Reaction5".to_string()],
            ),
            ("Library4".to_string(), vec!["Reaction6".to_string()]),
        ]));
        kin_data.vec_of_reaction_Values =
            Some(vec![json!({"data": "value4"}), json!({"data": "value5"})]);

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());

        let content = fs::read_to_string(file_path).unwrap();
        assert!(content.contains("KINETICS"));
        assert!(content.contains("Library3"));
        assert!(content.contains("Library4"));
        assert!(content.contains("Reaction4"));
        assert!(content.contains("Reaction5"));
        assert!(content.contains("Reaction6"));
        assert!(content.contains("value4"));
        assert!(content.contains("value5"));

        // Test with existing file and header
        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();
        fs::write(file_path, "KINETICS\nSome existing content\n").unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());

        let content = fs::read_to_string(file_path).unwrap();
        assert_eq!(content.matches("KINETICS").count(), 1);
        assert!(content.contains("Some existing content"));
        assert!(content.contains("Library3"));
        assert!(content.contains("Library4"));
    }
}
