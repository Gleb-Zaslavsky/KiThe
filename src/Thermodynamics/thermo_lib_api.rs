use prettytable::{Table, row};
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};

use std::fs::OpenOptions;

use crate::library_manager::with_library_manager;
use std::io::{BufRead, BufReader};
use std::path::Path;
// /Basis functionality to search in library of thermodymical and heat mass transfer data
#[derive(Debug, Clone)]
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
    pub thermo_libs: Vec<String>,
    pub transport_libs: Vec<String>,
}
impl ThermoData {
    pub fn new() -> Self {
        let mut file =
            with_library_manager(|manager| File::open(manager.all_keys_substance_path())).unwrap();
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).unwrap();
        // (lib:sustance)
        let lib_sub_pairs_vec: Vec<(String, String)> =
            serde_json::from_str(&file_contents).unwrap();

        let mut file =
            with_library_manager(|manager| File::open(manager.substance_base_path())).unwrap();
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).unwrap();
        // lib:{ substance :{ data }}
        let LibThermoData: HashMap<String, HashMap<String, Value>> =
            serde_json::from_str::<HashMap<String, HashMap<String, _>>>(&file_contents).unwrap();
        let AllLibraries: Vec<String> = LibThermoData.keys().map(|k| k.to_string()).collect();

        let mut thermo_libs = Vec::new();
        let mut transport_libs = Vec::new();
        /*
        Iterates through all library names
        For each library, takes the first record from its data
        Checks if the record contains a "Cp" key
        If "Cp" exists, adds the library to thermo_libs
        Otherwise, adds it to transport_libs
        */
        for lib_name in &AllLibraries {
            if let Some(lib_data) = LibThermoData.get(lib_name) {
                if let Some((_, record)) = lib_data.iter().next() {
                    if record.get("Cp").is_some() {
                        thermo_libs.push(lib_name.clone());
                    } else {
                        transport_libs.push(lib_name.clone());
                    }
                }
            }
        }

        Self {
            VecOfSubsAdresses: lib_sub_pairs_vec.clone(),
            LibThermoData: LibThermoData.clone(),
            AllLibraries: AllLibraries.clone(),
            AllSubstances: Vec::new(),
            subs_to_search: Vec::new(),
            hashmap_of_thermo_data: HashMap::new(),
            thermo_libs,
            transport_libs,
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

    pub fn create_substance_document(
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
            writeln!(file, "SUBSTANCES DATA")?;
        } else {
            let reader = BufReader::new(File::open(file_name)?);
            let has_header = reader.lines().any(|line| {
                line.as_ref().unwrap().clone().trim().to_uppercase() == "SUBSTANCES DATA"
                    || line.as_ref().unwrap().clone().trim().to_uppercase() == "SUBS DATA"
            });

            if !has_header {
                writeln!(file, "\nSUBSTANCES DATA")?;
            }
        }
        let hashmap_to_save: HashMap<String, HashMap<String, Value>> =
            self.hashmap_of_thermo_data.clone();
        if hashmap_to_save.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "hashmap_to_save is empty",
            ));
        }

        writeln!(file, "{}", serde_json::to_string_pretty(&hashmap_to_save)?)?;

        Ok(hashmap_to_save)
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_create_substance_document_new_file() {
        let mut thermo_data = ThermoData::new();
        let mut test_data = HashMap::new();
        let mut substance_data = HashMap::new();
        substance_data.insert(
            "test_lib".to_string(),
            serde_json::json!({"test_key": "test_value"}),
        );
        test_data.insert("test_substance".to_string(), substance_data);
        thermo_data.hashmap_of_thermo_data = test_data;

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = thermo_data.create_substance_document(file_path);
        println!("result: {:#?}", result);
        assert!(result.is_ok());

        let mut file_content = String::new();
        temp_file
            .as_file()
            .read_to_string(&mut file_content)
            .unwrap();

        assert!(file_content.contains("SUBSTANCES DATA"));
        assert!(file_content.contains("test_substance"));
        assert!(file_content.contains("test_lib"));
        assert!(file_content.contains("test_key"));
        assert!(file_content.contains("test_value"));
    }
    #[test]
    fn NamedTempFiletest() {
        let text = "Brian was here. Briefly.";

        // Create a file inside of `env::temp_dir()`.
        let mut named_temp_file = NamedTempFile::new().unwrap();

        // Re-open it.
        let mut file2 = named_temp_file.reopen().unwrap();

        // Write some test data to the first handle.
        named_temp_file.write_all(text.as_bytes()).unwrap();

        // Read the test data using the second handle.
        let mut buf = String::new();
        file2.read_to_string(&mut buf).unwrap();
        assert_eq!(buf, text);
    }
    #[test]
    fn test_create_substance_document_existing_file() {
        let mut thermo_data = ThermoData::new();
        let mut test_data = HashMap::new();
        let mut substance_data = HashMap::new();
        substance_data.insert(
            "test_lib".to_string(),
            serde_json::json!({"test_key": "test_value"}),
        );
        test_data.insert("test_substance".to_string(), substance_data);
        thermo_data.hashmap_of_thermo_data = test_data;

        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        temp_file.flush().unwrap();
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        println!("content \n {:#?}", file_content);
        assert!(file_content.contains("Existing content"));

        let file_path = temp_file.path().to_str().unwrap();
        let result = thermo_data.create_substance_document(file_path);
        assert!(result.is_ok());

        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        println!("content \n \n {:#?}", file_content);
        //    assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
        assert!(file_content.contains("test_substance"));
        assert!(file_content.contains("test_lib"));
        assert!(file_content.contains("test_key"));
        assert!(file_content.contains("test_value"));
    }

    #[test]
    fn test_create_substance_document_empty_data() {
        let thermo_data = ThermoData::new();
        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = thermo_data.create_substance_document(file_path);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().kind(), std::io::ErrorKind::InvalidInput);
    }
    #[test]
    fn test_with_real_data() {
        let mut td = ThermoData::new();
        let subdata = td.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = subdata.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        td.search_libs_for_subs(
            vec!["CO".to_owned()],
            Some(vec!["Cantera_nasa_base_gas".to_owned()]),
        );
        println!(
            "td.hashmap_of_thermo_data: {:#?}",
            td.hashmap_of_thermo_data
        );
        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        let file_path = temp_file.path().to_str().unwrap();
        let result = td.create_substance_document(file_path);
        assert!(result.is_ok());
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        println!("content \n \n {:#?}", file_content);
        assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
    }
    #[test]
    fn test_with_real_data2() {
        let mut td = ThermoData::new();
        let subdata = td.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = subdata.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        td.search_libs_for_subs(
            vec!["CO".to_owned(), "O2".to_owned()],
            Some(vec![
                "Cantera_nasa_base_gas".to_owned(),
                "Aramco_transport".to_owned(),
            ]),
        );
        println!(
            "td.hashmap_of_thermo_data: {:#?}",
            td.hashmap_of_thermo_data
        );
        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        let file_path = temp_file.path().to_str().unwrap();
        let result = td.create_substance_document(file_path);
        assert!(result.is_ok());
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        // println!("content \n \n {:#?}", file_content);
        assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
    }

    #[test]
    fn test_iterate_all_libraries() {
        let thermo_data = ThermoData::new();

        for library in &thermo_data.AllLibraries {
            println!("Library: {}", library);

            if let Some(lib_data) = thermo_data.LibThermoData.get(library) {
                if let Some((substance, data)) = lib_data.iter().next() {
                    println!("  Sample substance: {}", substance);
                    println!("  Data: {}", data);
                }
            }
            println!();
        }
    }

    #[test]
    fn test_library_classification() {
        let thermo_data = ThermoData::new();

        println!("Thermo libraries: {:?}", thermo_data.thermo_libs);
        println!("Transport libraries: {:?}", thermo_data.transport_libs);
        println!("All libraries: {:?}", thermo_data.AllLibraries);
        assert!(!thermo_data.thermo_libs.is_empty() || !thermo_data.transport_libs.is_empty());

        for lib in &thermo_data.thermo_libs {
            if let Some(lib_data) = thermo_data.LibThermoData.get(lib) {
                if let Some((_, record)) = lib_data.iter().next() {
                    assert!(
                        record.get("Cp").is_some(),
                        "Thermo lib {} should have Cp key",
                        lib
                    );
                }
            }
        }

        for lib in &thermo_data.transport_libs {
            if let Some(lib_data) = thermo_data.LibThermoData.get(lib) {
                if let Some((_, record)) = lib_data.iter().next() {
                    assert!(
                        record.get("Cp").is_none(),
                        "Transport lib {} should not have Cp key",
                        lib
                    );
                }
            }
        }
    }
}
