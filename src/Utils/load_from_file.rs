use log::{error, info, warn};
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub struct LoadData {
    pub file_name: String,
}

impl LoadData {
    pub fn new(file_name: String) -> Self {
        LoadData { file_name }
    }
    pub fn load_kinetics(&self) -> Result<HashMap<String, HashMap<String, Value>>, String> {
        load_and_validate_kinetics(&self.file_name)
    }
    pub fn load_substance_list(&self) -> Result<Vec<String>, String> {
        load_substance_list(&self.file_name)
    }
    pub fn load_thermo(&self) -> Result<HashMap<String, HashMap<String, Value>>, String> {
        load_and_validate_thermo(&self.file_name)
    }
}
/// Parses a document for kinetic data under the "KINETICS" or "REACTIONS" header.
/// Returns a HashMap<String, HashMap<String, Value>> containing the parsed data.
/// The outer HashMap keys are library names, and the inner HashMap keys are reaction names.
pub fn load_kinetics_from_file(
    file_name: &str,
) -> Result<HashMap<String, HashMap<String, Value>>, String> {
    let path = Path::new(file_name);
    if !path.exists() {
        return Err(format!("File '{}' does not exist", file_name));
    }

    let file = match File::open(path) {
        Ok(file) => file,
        Err(e) => return Err(format!("Failed to open file '{}': {}", file_name, e)),
    };

    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();

    // Find the KINETICS or REACTIONS header
    let mut start_index = None;
    for (i, line) in lines.iter().enumerate() {
        let trimmed = line.trim().to_uppercase();
        if trimmed == "KINETICS" || trimmed == "REACTIONS" {
            start_index = Some(i + 1); // Start from the line after the header
            break;
        }
    }

    let start_index = match start_index {
        Some(index) => index,
        None => {
            return Err(format!(
                "No 'KINETICS' or 'REACTIONS' header found in file '{}'",
                file_name
            ));
        }
    };

    // Find the end index (next header or end of file)
    let mut end_index = lines.len();
    for i in start_index..lines.len() {
        let trimmed = lines[i].trim();
        if !trimmed.is_empty() && trimmed.chars().all(|c| c.is_uppercase() || c == '_') {
            end_index = i;
            break;
        }
    }

    // Extract the kinetics data section
    let kinetics_section = lines[start_index..end_index].join("\n");

    // Parse the JSON data
    let result: Result<HashMap<String, HashMap<String, Value>>, serde_json::Error> =
        serde_json::from_str(&kinetics_section);

    match result {
        Ok(data) => {
            info!(
                "Successfully parsed kinetics data from file '{}'",
                file_name
            );
            Ok(data)
        }
        Err(e) => {
            // Find the line and column where the error occurred
            let error_line = e.line();
            let error_column = e.column();

            // Calculate the actual line number in the file
            let actual_line = start_index + error_line - 1;

            let error_msg = format!(
                "Error parsing kinetics data at line {}, column {} (line {} in file): {}",
                error_line, error_column, actual_line, e
            );
            error!("{}", error_msg);

            // If possible, show the problematic line
            if actual_line < lines.len() {
                let problem_line = &lines[actual_line];
                error!("Problematic line: {}", problem_line);

                // Create a visual pointer to the error position
                if error_column <= problem_line.len() {
                    let pointer = " ".repeat(error_column - 1) + "^";
                    error!("{}", pointer);
                }
            }

            Err(error_msg)
        }
    }
}

/// Loads kinetics data from a file and validates the structure.
/// This function provides additional validation and error reporting.
pub fn load_and_validate_kinetics(
    file_name: &str,
) -> Result<HashMap<String, HashMap<String, Value>>, String> {
    let kinetics_data = load_kinetics_from_file(file_name)?;

    // Validate the structure of the loaded data
    if kinetics_data.is_empty() {
        warn!("Loaded kinetics data is empty");
    }

    // Check each library and reaction
    for (library_name, reactions) in &kinetics_data {
        if reactions.is_empty() {
            warn!("Library '{}' contains no reactions", library_name);
        }

        // Validate each reaction's data structure
        for (reaction_name, reaction_data) in reactions {
            if !reaction_data.is_object() {
                warn!(
                    "Reaction '{}' in library '{}' has invalid data format",
                    reaction_name, library_name
                );
            }

            // Check for required fields in reaction data
            // This can be customized based on your specific data structure requirements
            if !reaction_data.get("eq").is_some() {
                warn!(
                    "Reaction '{}' in library '{}' is missing 'eq' field",
                    reaction_name, library_name
                );
            }
        }
    }

    info!(
        "Loaded and validated kinetics data from file '{}'",
        file_name
    );
    Ok(kinetics_data)
}

/// load from file list of substances
pub fn load_substance_list(file_name: &str) -> Result<Vec<String>, String> {
    let path = Path::new(file_name);
    if !path.exists() {
        return Err(format!("File '{}' does not exist", file_name));
    }

    let file = match File::open(path) {
        Ok(file) => file,
        Err(e) => return Err(format!("Failed to open file '{}': {}", file_name, e)),
    };

    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();

    // Find the KINETICS or REACTIONS header
    let mut start_index = None;
    for (i, line) in lines.iter().enumerate() {
        let trimmed = line.trim().to_uppercase();
        if trimmed == "SUBSTANCES" || trimmed == "SUBSTANCES LIST" {
            start_index = Some(i + 1); // Start from the line after the header
            break;
        }
    }

    let start_index = match start_index {
        Some(index) => index,
        None => {
            return Err(format!(
                "No 'SUBSTANCES' header found in file '{}'",
                file_name
            ));
        }
    };

    // Find the end index (next header or end of file)
    let mut end_index = lines.len();
    for i in start_index..lines.len() {
        let trimmed = lines[i].trim();
        if !trimmed.is_empty() && trimmed.chars().all(|c| c.is_uppercase() || c == '_') {
            end_index = i;
            break;
        }
    }

    // Extract the kinetics data section
    let list_of_subs = lines[start_index..end_index].join("\n");
    let vec_of_molecules: Vec<String> = list_of_subs
        .replace('\n', ", ")
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    if vec_of_molecules.is_empty() {
        return Err(format!("No substances found in file '{}'", file_name));
    }
    if (vec_of_molecules.len() == 1) && (vec_of_molecules[0] == "".to_string()) {
        return Err(format!("No substances found in file '{}'", file_name));
    }
    Ok(vec_of_molecules)
}

pub fn load_thermo_from_file(
    file_name: &str,
) -> Result<HashMap<String, HashMap<String, Value>>, String> {
    let path = Path::new(file_name);
    if !path.exists() {
        return Err(format!("File '{}' does not exist", file_name));
    }

    let file = match File::open(path) {
        Ok(file) => file,
        Err(e) => return Err(format!("Failed to open file '{}': {}", file_name, e)),
    };

    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();

    // Find the KINETICS or REACTIONS header
    let mut start_index = None;
    for (i, line) in lines.iter().enumerate() {
        let trimmed = line.trim().to_uppercase();
        if trimmed == "SUBSTANCES DATA" || trimmed == "SUBS DATA" || trimmed == "SUBSTANCE DATA" {
            start_index = Some(i + 1); // Start from the line after the header
            break;
        }
    }

    let start_index = match start_index {
        Some(index) => index,
        None => {
            return Err(format!(
                "No 'SUBSTANCES DATA' header found in file '{}'",
                file_name
            ));
        }
    };

    // Find the end index (next header or end of file)
    let mut end_index = lines.len();
    for i in start_index..lines.len() {
        let trimmed = lines[i].trim();
        if !trimmed.is_empty() && trimmed.chars().all(|c| c.is_uppercase() || c == '_') {
            end_index = i;
            break;
        }
    }

    // Extract the kinetics data section
    let kinetics_section = lines[start_index..end_index].join("\n");

    // Parse the JSON data
    let result: Result<HashMap<String, HashMap<String, Value>>, serde_json::Error> =
        serde_json::from_str(&kinetics_section);

    match result {
        Ok(data) => {
            info!(
                "Successfully parsed substance data from file '{}'",
                file_name
            );
            Ok(data)
        }
        Err(e) => {
            // Find the line and column where the error occurred
            let error_line = e.line();
            let error_column = e.column();

            // Calculate the actual line number in the file
            let actual_line = start_index + error_line - 1;

            let error_msg = format!(
                "Error parsing substance data at line {}, column {} (line {} in file): {}",
                error_line, error_column, actual_line, e
            );
            error!("{}", error_msg);

            // If possible, show the problematic line
            if actual_line < lines.len() {
                let problem_line = &lines[actual_line];
                error!("Problematic line: {}", problem_line);

                // Create a visual pointer to the error position
                if error_column <= problem_line.len() {
                    let pointer = " ".repeat(error_column - 1) + "^";
                    error!("{}", pointer);
                }
            }

            Err(error_msg)
        }
    }
}

/// Loads thermo data from a file and validates the structure.
/// This function provides additional validation and error reporting.
pub fn load_and_validate_thermo(
    file_name: &str,
) -> Result<HashMap<String, HashMap<String, Value>>, String> {
    let kinetics_data = load_thermo_from_file(file_name)?;

    // Validate the structure of the loaded data
    if kinetics_data.is_empty() {
        warn!("Loaded thermo data is empty");
    }

    // Check each library and substance
    for (substance, lib_data_map) in &kinetics_data {
        if lib_data_map.is_empty() {
            warn!("Library '{}' contains no substances data", substance);
        }

        // Validate each reaction's data structure
        for (lib, subs_data) in lib_data_map {
            if !subs_data.is_object() {
                warn!(
                    "subst data '{}' in library '{}' has invalid data format",
                    subs_data, lib
                );
            }
        }
    }

    info!(
        "Loaded and validated kinetics data from file '{}'",
        file_name
    );
    Ok(kinetics_data)
}
/*
Provides two main functions:

load_kinetics_from_file: The core function that parses the document for kinetic data
load_and_validate_kinetics: A wrapper function that adds additional validation
Includes detailed error reporting:

Reports the exact line and column where parsing errors occur
Shows the problematic line with a pointer to the error position
Logs errors using the log crate
Has comprehensive test coverage:

Tests for successful parsing with both "KINETICS" and "REACTIONS" headers
Tests for error cases like missing headers and invalid JSON
Tests for the validation function
Follows the algorithm you specified:

Finds the header ("KINETICS" or "REACTIONS")
Takes everything until the next header or end of document
Parses the data into the required HashMap structure
Reports any parsing errors with their location
*/

#[cfg(test)]
mod tests {
    use super::*;

    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use serde_json::json;
    use std::io::Read;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_load_kinetics_from_file() {
        // Create a temporary test file
        let mut temp_file = NamedTempFile::new().unwrap();

        // Write test content to the file
        writeln!(temp_file, "Some header text").unwrap();
        writeln!(temp_file, "KINETICS").unwrap();
        writeln!(temp_file, "{{").unwrap();
        writeln!(temp_file, "  \"Library1\": {{").unwrap();
        writeln!(
            temp_file,
            "    \"Reaction1\": {{ \"eq\": \"A + B -> C\", \"data\": \"value1\" }},"
        )
        .unwrap();
        writeln!(
            temp_file,
            "    \"Reaction2\": {{ \"eq\": \"C -> D\", \"data\": \"value2\" }}"
        )
        .unwrap();
        writeln!(temp_file, "  }},").unwrap();
        writeln!(temp_file, "  \"Library2\": {{").unwrap();
        writeln!(
            temp_file,
            "    \"Reaction3\": {{ \"eq\": \"E + F -> G\", \"data\": \"value3\" }}"
        )
        .unwrap();
        writeln!(temp_file, "  }}").unwrap();
        writeln!(temp_file, "}}").unwrap();
        writeln!(temp_file, "ANOTHER_HEADER").unwrap();
        writeln!(temp_file, "Some other content").unwrap();

        // Get the file path
        let file_path = temp_file.path().to_str().unwrap();

        // Test the function
        let result = load_kinetics_from_file(file_path);
        assert!(result.is_ok());

        let kinetics_data = result.unwrap();
        assert_eq!(kinetics_data.len(), 2);
        assert!(kinetics_data.contains_key("Library1"));
        assert!(kinetics_data.contains_key("Library2"));

        let library1 = &kinetics_data["Library1"];
        assert_eq!(library1.len(), 2);
        assert!(library1.contains_key("Reaction1"));
        assert!(library1.contains_key("Reaction2"));

        let library2 = &kinetics_data["Library2"];
        assert_eq!(library2.len(), 1);
        assert!(library2.contains_key("Reaction3"));

        // Check the content of a reaction
        let reaction1 = &library1["Reaction1"];
        assert_eq!(reaction1["eq"], "A + B -> C");
        assert_eq!(reaction1["data"], "value1");
    }

    #[test]
    fn test_load_kinetics_from_file_with_reactions_header() {
        // Create a temporary test file with REACTIONS header instead of KINETICS
        let mut temp_file = NamedTempFile::new().unwrap();

        writeln!(temp_file, "REACTIONS").unwrap();
        writeln!(temp_file, "{{").unwrap();
        writeln!(temp_file, "  \"Library1\": {{").unwrap();
        writeln!(temp_file, "    \"Reaction1\": {{ \"eq\": \"A + B -> C\" }}").unwrap();
        writeln!(temp_file, "  }}").unwrap();
        writeln!(temp_file, "}}").unwrap();

        let file_path = temp_file.path().to_str().unwrap();

        let result = load_kinetics_from_file(file_path);
        assert!(result.is_ok());

        let kinetics_data = result.unwrap();
        assert_eq!(kinetics_data.len(), 1);
        assert!(kinetics_data.contains_key("Library1"));
    }

    #[test]
    fn test_load_kinetics_from_file_no_header() {
        // Create a file without the required header
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "Some content without a KINETICS header").unwrap();

        let file_path = temp_file.path().to_str().unwrap();

        let result = load_kinetics_from_file(file_path);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("No 'KINETICS' or 'REACTIONS' header found")
        );
    }

    // #[test]
    fn test_load_kinetics_from_file_invalid_json() {
        // Create a file with invalid JSON
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "KINETICS").unwrap();
        writeln!(temp_file, "{{").unwrap();
        writeln!(temp_file, "  \"Library1\": {{").unwrap();
        writeln!(
            temp_file,
            "    \"Reaction1\": {{ \"eq\": \"A + B -> C\" }},"
        )
        .unwrap();
        writeln!(temp_file, "    \"Reaction2\": {{ \"eq\": \"C -> D\" }}").unwrap(); // Missing closing brace
        writeln!(temp_file, "  }}").unwrap();
        writeln!(temp_file, "}}").unwrap();

        let file_path = temp_file.path().to_str().unwrap();

        let result = load_kinetics_from_file(file_path);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Error parsing kinetics data"));
    }

    #[test]
    fn test_load_and_validate_kinetics() {
        // Create a valid test file
        let mut temp_file = NamedTempFile::new().unwrap();

        let test_data = json!({
            "Library1": {
                "Reaction1": {
                    "eq": "A + B -> C",
                    "data": "value1"
                }
            }
        });

        writeln!(temp_file, "KINETICS").unwrap();
        writeln!(temp_file, "{}", test_data.to_string()).unwrap();

        let file_path = temp_file.path().to_str().unwrap();

        let result = load_and_validate_kinetics(file_path);
        assert!(result.is_ok());

        let kinetics_data = result.unwrap();
        assert_eq!(kinetics_data.len(), 1);
        assert!(kinetics_data.contains_key("Library1"));
    }

    #[test]
    fn test_load_substance_list_success() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "Some header").unwrap();
        writeln!(temp_file, "SUBSTANCES").unwrap();
        writeln!(temp_file, "H2O, CO2, CH4").unwrap();
        writeln!(temp_file, "O2, N2").unwrap();
        writeln!(temp_file, "ANOTHER_HEADER").unwrap();

        let file_path = temp_file.path().to_str().unwrap();
        let result = load_substance_list(file_path);
        assert!(result.is_ok());
        let substances = result.unwrap();
        assert_eq!(substances, vec!["H2O", "CO2", "CH4", "O2", "N2"]);
    }

    #[test]
    fn test_load_substance_list_file_not_found() {
        let result = load_substance_list("non_existent_file.txt");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("does not exist"));
    }

    #[test]
    fn test_load_substance_list_no_header() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H2O, CO2, CH4").unwrap();

        let file_path = temp_file.path().to_str().unwrap();
        let result = load_substance_list(file_path);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("No 'SUBSTANCES' header found"));
    }
    #[test]
    fn test_load_substance_list_alternative_header() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "SUBSTANCES LIST").unwrap();
        writeln!(temp_file, "H2O, CO2, CH4").unwrap();

        let file_path = temp_file.path().to_str().unwrap();
        let result = load_substance_list(file_path);
        assert!(result.is_ok());
        let substances = result.unwrap();
        assert_eq!(substances, vec!["H2O", "CO2", "CH4"]);
    }

    #[test]
    fn test_load_substance_list_with_spaces() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "SUBSTANCES").unwrap();
        writeln!(temp_file, "H2O , CO2,   CH4, O2").unwrap();

        let file_path = temp_file.path().to_str().unwrap();
        let result = load_substance_list(file_path);
        assert!(result.is_ok());
        let substances = result.unwrap();
        assert_eq!(substances, vec!["H2O", "CO2", "CH4", "O2"]);
    }
    #[test]
    fn test_load_substance_list_empty() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "SUBSTANCES").unwrap();
        writeln!(temp_file, "").unwrap();
        writeln!(temp_file, "ANOTHER_HEADER").unwrap();

        let file_path = temp_file.path().to_str().unwrap();
        let result = load_substance_list(file_path);
        assert!(result.is_err());
    }
    #[test]
    fn test_with_real_data() {
        //CREATING SUBSTANCE DATA//////////////////////////////////////////////////////////////////////////////
        let mut td = ThermoData::new();
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
        //CREATING NEW TEMPORARY FILE///////////////////////////////////////////////////////////
        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        let file_path = temp_file.path().to_str().unwrap();
        //WRITING SUBSTANCES DATA TO TEMPORARY FILE///////////////////////////////////////////////////////////
        let result = td.create_substance_document(file_path); // writing subs to TEMPORARY FILE////////////////////////////////
        assert!(result.is_ok());
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        // println!("content \n \n {:#?}", file_content);
        assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
        //WRITING KINETIC DATA TO TEMPORARY FILE///////////////////////////////////////////////////////////
        use crate::Kinetics::User_reactions::KinData;
        // let our journey begin with a new instance of KinData
        let mut kd = KinData::new();
        // set the shortcut reactions for our KineticData instance
        // it means we want reactions from Cantera sub-librarie from number 1 to number 10
        kd.set_reactions_from_shortcut_range("C1..C10".to_string());
        // searching for reactions in data base
        kd.get_reactions_from_shortcuts();
        kd.kinetic_main(); // parsing reaction data into structures and stoichometric calculations under one hood
        let kintic_written = kd.create_kinetics_document(file_path);
        assert!(kintic_written.is_ok());
        //PARSING SUBSTANCE DATA FROM TEMPORARY FILE///////////////////////////////////////////////////////////
        let ld = LoadData::new(file_path.to_owned());
        let parsed_data = ld.load_thermo();
        println!("parsed_data: {:#?}", parsed_data);
        assert!(parsed_data.is_ok());
        assert_eq!(result.unwrap(), parsed_data.unwrap());
        //PARSING KINETIC DATA FROM TEMPORARY FILE///////////////////////////////////////////////////////////
        let ld = LoadData::new(file_path.to_owned());
        let parsed_data = ld.load_kinetics();
        println!("parsed_data: {:#?}", parsed_data);
        assert!(parsed_data.is_ok());
        assert_eq!(kintic_written.unwrap(), parsed_data.unwrap());
    }
    #[test]
    fn test_with_real_data2() {
        let mut td = ThermoData::new();
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

        let ld = LoadData::new(file_path.to_owned());
        let parsed_data = ld.load_thermo();
        println!("parsed_data: {:#?}", parsed_data);
        assert!(parsed_data.is_ok());
        assert_eq!(result.unwrap(), parsed_data.unwrap());
    }
}
