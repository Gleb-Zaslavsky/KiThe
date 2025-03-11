use prettytable::{Table, row};
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
// /Basis functionality to search in library of thermodymical and heat mass transfer data
pub struct ThermoData {
    pub VecOfSubsAdresses: Vec<(String, String)>,
    /// adresses of all substances (lib:sustance)
    pub LibThermoData: HashMap<String, HashMap<String, Value>>,
    /// library of thermodymical and heat mass transfer data  lib:{ substance :{ data }}
    pub AllLibraries: Vec<String>,
    /// all sub-libraries in big library
    pub subs_to_search: Vec<String>,
    /// all substances to search in libraries
    pub hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>>,
    /// search result {substance: {lib: data}}
    pub AllSubstances: Vec<String>, // all substances in the library
}
impl ThermoData {
    pub fn new() -> Self {
        let mut file = File::open("all_keys_substance.json").unwrap();
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).unwrap();
        // (lib:sustance)
        let lib_sub_pairs_vec: Vec<(String, String)> =
            serde_json::from_str(&file_contents).unwrap();

        let mut file = File::open("substance_base_v2.json").unwrap();
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).unwrap();
        // lib:{ substance :{ data }}
        let LibThermoData: HashMap<String, HashMap<String, Value>> =
            serde_json::from_str::<HashMap<String, HashMap<String, _>>>(&file_contents).unwrap();
        let AllLibraries: Vec<String> = LibThermoData.keys().map(|k| k.to_string()).collect();

        Self {
            VecOfSubsAdresses: lib_sub_pairs_vec.clone(),
            LibThermoData: LibThermoData.clone(),
            AllLibraries: AllLibraries.clone(),
            AllSubstances: Vec::new(),
            subs_to_search: Vec::new(),
            hashmap_of_thermo_data: HashMap::new(),
        }
    }
    //////////////////////////////////////////////DATA MANIPULATION//////////////////////////////////////////////////////////
    /// add substances in the form {substanse:{type_of_data*:{data} }}
    /// *type of data may be 1) NASA format of thermodynamic data 2) NIST format of thermodynamic data 3) "transport" format of heat mass transfer data
    pub fn append_substances_data(
        &mut self,
        hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>>,
    ) {
        self.hashmap_of_thermo_data
            .extend(hashmap_of_thermo_data.clone());
        self.subs_to_search
            .extend(hashmap_of_thermo_data.keys().cloned());
    }
    pub fn subs_of_this_lib(&mut self, lib: &str) -> Vec<String> {
        let subs_of_this_lib: Vec<String> = self
            .VecOfSubsAdresses
            .iter()
            .filter(|(k, _)| k == lib)
            .map(|(_, v)| v.clone())
            .collect();
        subs_of_this_lib
    }

    ///  remove substance from hashmap and vector of substances
    pub fn remove_by_name(&mut self, name: String) {
        self.subs_to_search.remove(
            self.subs_to_search
                .iter()
                .position(|x| x == &name)
                .expect("substance not found"),
        );
        self.hashmap_of_thermo_data.remove(&name);
    }
    /////////////////////////////////PARSE DATA FROM LIBRARIES////////////////////////////////////////////////
    /// Search for a specific substance in all reaction libraries and return all libraries where it is found
    pub fn search_libs(
        &mut self,
        substance: &str,
        allowed_libs: Option<Vec<String>>,
    ) -> Vec<String> {
        let libs_for_this_sub: Vec<String> = self
            .VecOfSubsAdresses
            .iter()
            .filter(|(_, subs)| subs == substance)
            .map(|(library, _)| Self::filter_library(library, allowed_libs.clone()))
            .collect();
        libs_for_this_sub
    }

    fn filter_library(library: &str, allowed_libs: Option<Vec<String>>) -> String {
        if let Some(allowed_libs) = allowed_libs {
            if allowed_libs.contains(&library.to_string()) {
                library.to_string()
            } else {
                String::new()
            }
        } else {
            library.to_string()
        }
    }
    /// search a substances from vector in all libraries and making a hashmap of information in the form {substanse:{type_of_data*:{data} }}
    /// *type of data may be 1) NASA format of thermodynamic data 2) NIST format of thermodynamic data 3) "transport" format of heat mass transfer data
    pub fn search_libs_for_subs(
        &mut self,
        substances: Vec<String>,
        allowed_libs: Option<Vec<String>>,
    ) -> Vec<String> {
        let mut hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>> = HashMap::new();
        let mut libs_for_this_subs = Vec::new();

        for sub in substances {
            // for this sub let us find all libraries where this sub is present
            let libs_for_this_sub: Vec<String> = self.search_libs(&sub, allowed_libs.clone());
            let hashmap_lib_data: HashMap<String, Value> = libs_for_this_sub // iterating through all libraries where this sub is present and finding data for this sub
                .iter()
                .filter_map(|lib| {
                    // search data for this sub in this libraries
                    self.LibThermoData
                        .get(Self::rename(&lib)) // get data for this sub from this library and save it in hashmap {lib:substance}
                        .and_then(|lib_data| lib_data.get(&sub).map(|v| ((lib.clone()), v.clone())))
                })
                .collect();
            hashmap_of_thermo_data.insert(sub.clone(), hashmap_lib_data);
            libs_for_this_subs.extend(libs_for_this_sub);
        }

        self.hashmap_of_thermo_data = hashmap_of_thermo_data.clone();
        libs_for_this_subs
    }
    pub fn rename(lib_name: &str) -> &str {
        match lib_name {
            "Cantera_nasa_base_gas" => "NASA_gas",
            _ => lib_name,
        }
    }
    ///////////////////INPUT/OUTPUT/////////////////////////////////////////////////
    ///
    pub fn save_thermo_data_to_json(
        &self,
        filename: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let json = serde_json::to_string(&self.hashmap_of_thermo_data)?;
        let mut file = File::create(filename)?;
        file.write_all(json.as_bytes())?;
        Ok(())
    }

    pub fn load_thermo_data_from_json(
        &mut self,
        filename: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::open(filename)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        self.hashmap_of_thermo_data = serde_json::from_str(&contents)?;
        self.subs_to_search = self
            .hashmap_of_thermo_data
            .keys()
            .map(|k| k.to_string())
            .collect();
        Ok(())
    }
    pub fn pretty_print_thermo_data(&self) {
        let data = self.hashmap_of_thermo_data.clone();
        let mut table = Table::new();
        table.add_row(row!["Substance", "Library", "Data",]);
        for (sub, hashmap_lib_data) in data.iter() {
            for (lib, data) in hashmap_lib_data.iter() {
                table.add_row(row![sub, lib, data]);
            }
        }
        table.printstd();
    }
    pub fn what_handler_to_use(lib_name: &str) -> String {
        match lib_name {
            "Cantera_nasa_base_gas"
            | "NASA_gas"
            | "NASA_cond"
            | "Cantera_nasa_base_cond"
            | "nuig_thermo"
            | "NASA"
            | "NASA7" => "NASA".to_string(),
            "CEA" => "CEA".to_string(),
            "NIST" | "NIST9" => "NIST".to_string(),
            "Aramco_transpot" => "transport".to_string(),
            _ => lib_name.to_string(),
        }
    }
}
