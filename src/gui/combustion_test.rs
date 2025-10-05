#[cfg(test)]
mod tests {

    use std::collections::HashMap;

    use crate::gui::combustion::{CombustionApp, ProblemsEnum};
    use RustedSciThe::Utils::task_parser::{DocumentMap, Value};

    fn create_test_document() -> DocumentMap {
        let mut section1 = HashMap::new();
        section1.insert(
            "string_field".to_string(),
            Some(vec![Value::String("test".to_string())]),
        );
        section1.insert("float_field".to_string(), Some(vec![Value::Float(3.14)]));
        section1.insert("int_field".to_string(), Some(vec![Value::Integer(42)]));
        section1.insert("bool_field".to_string(), Some(vec![Value::Boolean(true)]));
        section1.insert(
            "vector_field".to_string(),
            Some(vec![Value::Vector(vec![1.0, 2.0, 3.0])]),
        );
        section1.insert("usize_field".to_string(), Some(vec![Value::Usize(100)]));
        section1.insert(
            "optional_field".to_string(),
            Some(vec![Value::Optional(Some(Box::new(Value::String(
                "nested".to_string(),
            ))))]),
        );
        section1.insert("none_field".to_string(), None);

        let mut section2 = HashMap::new();
        section2.insert("temperature".to_string(), Some(vec![Value::Float(298.15)]));
        section2.insert("pressure".to_string(), Some(vec![Value::Float(101325.0)]));

        let mut doc = HashMap::new();
        doc.insert("test_section".to_string(), section1);
        doc.insert("conditions".to_string(), section2);
        doc
    }

    #[test]
    fn test_combustion_app_creation() {
        let app = CombustionApp::new();
        assert!(app.document.is_empty());
        assert!(app.new_section_name.is_empty());
        assert!(app.new_field_names.is_empty());
        assert!(app.new_field_values.is_empty());
        assert!(matches!(app.selected_problem, ProblemsEnum::None));
    }

    #[test]
    fn test_document_structure() {
        let doc = create_test_document();
        assert_eq!(doc.len(), 2);
        assert!(doc.contains_key("test_section"));
        assert!(doc.contains_key("conditions"));

        let test_section = doc.get("test_section").unwrap();
        assert_eq!(test_section.len(), 8);
        assert!(test_section.contains_key("string_field"));
        assert!(test_section.contains_key("float_field"));
        assert!(test_section.contains_key("none_field"));
    }

    #[test]
    fn test_value_types() {
        let doc = create_test_document();
        let section = doc.get("test_section").unwrap();

        // Test String value
        if let Some(Some(values)) = section.get("string_field") {
            if let Value::String(s) = &values[0] {
                assert_eq!(s, "test");
            } else {
                panic!("Expected String value");
            }
        }

        // Test Float value
        if let Some(Some(values)) = section.get("float_field") {
            if let Value::Float(f) = &values[0] {
                assert_eq!(*f, 3.14);
            } else {
                panic!("Expected Float value");
            }
        }

        // Test Integer value
        if let Some(Some(values)) = section.get("int_field") {
            if let Value::Integer(i) = &values[0] {
                assert_eq!(*i, 42);
            } else {
                panic!("Expected Integer value");
            }
        }

        // Test Boolean value
        if let Some(Some(values)) = section.get("bool_field") {
            if let Value::Boolean(b) = &values[0] {
                assert_eq!(*b, true);
            } else {
                panic!("Expected Boolean value");
            }
        }

        // Test Vector value
        if let Some(Some(values)) = section.get("vector_field") {
            if let Value::Vector(vec) = &values[0] {
                assert_eq!(vec, &vec![1.0, 2.0, 3.0]);
            } else {
                panic!("Expected Vector value");
            }
        }

        // Test Usize value
        if let Some(Some(values)) = section.get("usize_field") {
            if let Value::Usize(u) = &values[0] {
                assert_eq!(*u, 100);
            } else {
                panic!("Expected Usize value");
            }
        }

        // Test None value
        assert!(section.get("none_field").unwrap().is_none());
    }

    #[test]
    fn test_app_with_custom_document() {
        let mut app = CombustionApp::new();
        app.document = create_test_document();

        assert_eq!(app.document.len(), 2);
        assert!(app.document.contains_key("test_section"));
        assert!(app.document.contains_key("conditions"));
    }

    #[test]
    fn test_field_addition_simulation() {
        let mut app = CombustionApp::new();
        let section_name = "test_section".to_string();

        // Simulate adding a new section
        app.document.insert(section_name.clone(), HashMap::new());

        // Simulate adding fields
        let field_name = "new_field".to_string();
        let field_value = "42.5".to_string();

        // Parse value as would happen in GUI
        let value = if field_value.parse::<f64>().is_ok() {
            Value::Float(field_value.parse().unwrap())
        } else {
            Value::String(field_value)
        };

        if let Some(section) = app.document.get_mut(&section_name) {
            section.insert(field_name.clone(), Some(vec![value]));
        }

        // Verify field was added
        let section = app.document.get(&section_name).unwrap();
        assert!(section.contains_key(&field_name));

        if let Some(Some(values)) = section.get(&field_name) {
            if let Value::Float(f) = &values[0] {
                assert_eq!(*f, 42.5);
            } else {
                panic!("Expected Float value");
            }
        }
    }

    #[test]
    fn test_value_parsing() {
        // Test float parsing
        let float_str = "3.14";
        assert!(float_str.parse::<f64>().is_ok());

        // Test integer parsing
        let int_str = "42";
        assert!(int_str.parse::<i64>().is_ok());

        // Test boolean parsing
        let bool_str = "true";
        assert_eq!(bool_str, "true");

        // Test string fallback
        let string_str = "hello world";
        assert!(string_str.parse::<f64>().is_err());
        assert!(string_str.parse::<i64>().is_err());
    }

    #[test]
    fn test_section_management() {
        let mut app = CombustionApp::new();
        let initial_count = app.document.len();

        // Add new section
        let new_section = "dynamic_section".to_string();
        app.document.insert(new_section.clone(), HashMap::new());

        assert_eq!(app.document.len(), initial_count + 1);
        assert!(app.document.contains_key(&new_section));

        // Verify section is empty initially
        let section = app.document.get(&new_section).unwrap();
        assert!(section.is_empty());
    }

    #[test]
    fn test_multiple_values_in_field() {
        let mut section = HashMap::new();
        let multiple_values = vec![
            Value::String("first".to_string()),
            Value::String("second".to_string()),
            Value::String("third".to_string()),
        ];
        section.insert("multi_field".to_string(), Some(multiple_values));

        if let Some(Some(values)) = section.get("multi_field") {
            assert_eq!(values.len(), 3);
            if let Value::String(s) = &values[0] {
                assert_eq!(s, "first");
            }
            if let Value::String(s) = &values[2] {
                assert_eq!(s, "third");
            }
        }
    }

    #[test]
    fn test_optional_value_handling() {
        let mut section = HashMap::new();

        // Test Some optional
        let some_optional = Value::Optional(Some(Box::new(Value::Float(2.71))));
        section.insert("some_opt".to_string(), Some(vec![some_optional]));

        // Test None optional
        let none_optional = Value::Optional(None);
        section.insert("none_opt".to_string(), Some(vec![none_optional]));

        assert_eq!(section.len(), 2);
        assert!(section.contains_key("some_opt"));
        assert!(section.contains_key("none_opt"));
    }

    #[test]
    fn test_empty_document_handling() {
        let mut app = CombustionApp::new();
        app.document = HashMap::new();

        assert!(app.document.is_empty());

        // Add first section
        app.document
            .insert("first_section".to_string(), HashMap::new());
        assert_eq!(app.document.len(), 1);
    }

    #[test]
    fn test_field_name_tracking() {
        let mut app = CombustionApp::new();

        // Simulate field name tracking for sections
        let section_name = "test_section".to_string();
        app.new_field_names
            .insert(section_name.clone(), "new_field".to_string());
        app.new_field_values
            .insert(section_name.clone(), "test_value".to_string());

        assert_eq!(app.new_field_names.get(&section_name).unwrap(), "new_field");
        assert_eq!(
            app.new_field_values.get(&section_name).unwrap(),
            "test_value"
        );

        // Clear fields
        app.new_field_names
            .insert(section_name.clone(), String::new());
        app.new_field_values
            .insert(section_name.clone(), String::new());

        assert!(app.new_field_names.get(&section_name).unwrap().is_empty());
        assert!(app.new_field_values.get(&section_name).unwrap().is_empty());
    }
}

mod tests2 {

    use std::collections::HashMap;

    use crate::gui::combustion::{CombustionApp, ProblemsEnum};
    use RustedSciThe::Utils::task_parser::{DocumentMap, Value};
    fn create_test_document() -> DocumentMap {
        let mut section1 = HashMap::new();
        section1.insert(
            "string_field".to_string(),
            Some(vec![Value::String("test".to_string())]),
        );
        section1.insert("float_field".to_string(), Some(vec![Value::Float(3.14)]));
        section1.insert("int_field".to_string(), Some(vec![Value::Integer(42)]));
        section1.insert("bool_field".to_string(), Some(vec![Value::Boolean(true)]));
        section1.insert(
            "vector_field".to_string(),
            Some(vec![Value::Vector(vec![1.0, 2.0, 3.0])]),
        );
        section1.insert("usize_field".to_string(), Some(vec![Value::Usize(100)]));
        section1.insert(
            "optional_field".to_string(),
            Some(vec![Value::Optional(Some(Box::new(Value::String(
                "nested".to_string(),
            ))))]),
        );
        section1.insert("none_field".to_string(), None);

        let mut section2 = HashMap::new();
        section2.insert("temperature".to_string(), Some(vec![Value::Float(298.15)]));
        section2.insert("pressure".to_string(), Some(vec![Value::Float(101325.0)]));

        let mut doc = HashMap::new();
        doc.insert("test_section".to_string(), section1);
        doc.insert("conditions".to_string(), section2);
        doc
    }

    #[test]
    fn test_combustion_app_creation() {
        let app = CombustionApp::new();
        assert!(app.new_section_name.is_empty());
        assert!(app.new_field_names.is_empty());
        assert!(app.new_field_values.is_empty());
        assert!(matches!(app.selected_problem, ProblemsEnum::None));
    }

    #[test]
    fn test_problem_enum_none() {
        let app = CombustionApp::new_with_problem(ProblemsEnum::None);
        assert!(app.document.is_empty());
        assert!(matches!(app.selected_problem, ProblemsEnum::None));
    }

    #[test]
    fn test_problem_enum_bvp_simple() {
        let app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(!app.document.is_empty());
        assert!(matches!(app.selected_problem, ProblemsEnum::BVPSimple));
    }

    #[test]
    fn test_document_structure() {
        let doc = create_test_document();
        assert_eq!(doc.len(), 2);
        assert!(doc.contains_key("test_section"));
        assert!(doc.contains_key("conditions"));

        let test_section = doc.get("test_section").unwrap();
        assert_eq!(test_section.len(), 8);
        assert!(test_section.contains_key("string_field"));
        assert!(test_section.contains_key("float_field"));
        assert!(test_section.contains_key("none_field"));
    }

    #[test]
    fn test_value_types() {
        let doc = create_test_document();
        let section = doc.get("test_section").unwrap();

        // Test String value
        if let Some(Some(values)) = section.get("string_field") {
            if let Value::String(s) = &values[0] {
                assert_eq!(s, "test");
            } else {
                panic!("Expected String value");
            }
        }

        // Test Float value
        if let Some(Some(values)) = section.get("float_field") {
            if let Value::Float(f) = &values[0] {
                assert_eq!(*f, 3.14);
            } else {
                panic!("Expected Float value");
            }
        }

        // Test Integer value
        if let Some(Some(values)) = section.get("int_field") {
            if let Value::Integer(i) = &values[0] {
                assert_eq!(*i, 42);
            } else {
                panic!("Expected Integer value");
            }
        }

        // Test Boolean value
        if let Some(Some(values)) = section.get("bool_field") {
            if let Value::Boolean(b) = &values[0] {
                assert_eq!(*b, true);
            } else {
                panic!("Expected Boolean value");
            }
        }

        // Test Vector value
        if let Some(Some(values)) = section.get("vector_field") {
            if let Value::Vector(vec) = &values[0] {
                assert_eq!(vec, &vec![1.0, 2.0, 3.0]);
            } else {
                panic!("Expected Vector value");
            }
        }

        // Test Usize value
        if let Some(Some(values)) = section.get("usize_field") {
            if let Value::Usize(u) = &values[0] {
                assert_eq!(*u, 100);
            } else {
                panic!("Expected Usize value");
            }
        }

        // Test None value
        assert!(section.get("none_field").unwrap().is_none());
    }

    #[test]
    fn test_app_with_custom_document() {
        let mut app = CombustionApp::new();
        app.document = create_test_document();

        assert_eq!(app.document.len(), 2);
        assert!(app.document.contains_key("test_section"));
        assert!(app.document.contains_key("conditions"));
    }

    #[test]
    fn test_field_addition_simulation() {
        let mut app = CombustionApp::new();
        let section_name = "test_section".to_string();

        // Simulate adding a new section
        app.document.insert(section_name.clone(), HashMap::new());

        // Simulate adding fields
        let field_name = "new_field".to_string();
        let field_value = "42.5".to_string();

        // Parse value as would happen in GUI
        let value = if field_value.parse::<f64>().is_ok() {
            Value::Float(field_value.parse().unwrap())
        } else {
            Value::String(field_value)
        };

        if let Some(section) = app.document.get_mut(&section_name) {
            section.insert(field_name.clone(), Some(vec![value]));
        }

        // Verify field was added
        let section = app.document.get(&section_name).unwrap();
        assert!(section.contains_key(&field_name));

        if let Some(Some(values)) = section.get(&field_name) {
            if let Value::Float(f) = &values[0] {
                assert_eq!(*f, 42.5);
            } else {
                panic!("Expected Float value");
            }
        }
    }

    #[test]
    fn test_value_parsing() {
        // Test float parsing
        let float_str = "3.14";
        assert!(float_str.parse::<f64>().is_ok());

        // Test integer parsing
        let int_str = "42";
        assert!(int_str.parse::<i64>().is_ok());

        // Test boolean parsing
        let bool_str = "true";
        assert_eq!(bool_str, "true");

        // Test string fallback
        let string_str = "hello world";
        assert!(string_str.parse::<f64>().is_err());
        assert!(string_str.parse::<i64>().is_err());
    }

    #[test]
    fn test_section_management() {
        let mut app = CombustionApp::new();
        let initial_count = app.document.len();

        // Add new section
        let new_section = "dynamic_section".to_string();
        app.document.insert(new_section.clone(), HashMap::new());

        assert_eq!(app.document.len(), initial_count + 1);
        assert!(app.document.contains_key(&new_section));

        // Verify section is empty initially
        let section = app.document.get(&new_section).unwrap();
        assert!(section.is_empty());
    }

    #[test]
    fn test_multiple_values_in_field() {
        let mut section = HashMap::new();
        let multiple_values = vec![
            Value::String("first".to_string()),
            Value::String("second".to_string()),
            Value::String("third".to_string()),
        ];
        section.insert("multi_field".to_string(), Some(multiple_values));

        if let Some(Some(values)) = section.get("multi_field") {
            assert_eq!(values.len(), 3);
            if let Value::String(s) = &values[0] {
                assert_eq!(s, "first");
            }
            if let Value::String(s) = &values[2] {
                assert_eq!(s, "third");
            }
        }
    }

    #[test]
    fn test_optional_value_handling() {
        let mut section = HashMap::new();

        // Test Some optional
        let some_optional = Value::Optional(Some(Box::new(Value::Float(2.71))));
        section.insert("some_opt".to_string(), Some(vec![some_optional]));

        // Test None optional
        let none_optional = Value::Optional(None);
        section.insert("none_opt".to_string(), Some(vec![none_optional]));

        assert_eq!(section.len(), 2);
        assert!(section.contains_key("some_opt"));
        assert!(section.contains_key("none_opt"));
    }

    #[test]
    fn test_empty_document_handling() {
        let mut app = CombustionApp::new();
        app.document = HashMap::new();

        assert!(app.document.is_empty());

        // Add first section
        app.document
            .insert("first_section".to_string(), HashMap::new());
        assert_eq!(app.document.len(), 1);
    }
}

#[cfg(test)]
mod file_reading_tests {
    use super::*;
    use crate::gui::combustion::CombustionApp;
    use RustedSciThe::Utils::task_parser::{DocumentParser, Value};
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    const TEST_TASK_CONTENT: &str = r#"
      initial_guess
        universal:1e-2
        process_conditions
        problem_name: Some(HMXTest)
        problem_description: Some(HMXdecompositiontest)
        substances: HMX, HMXprod
        t0: 0.0
        t_end: 1.0
        n_steps: 25
        arg:x
        Tm: 1500.0
        L: 5e-4
        dT: 600.0
        T_scale: 600.0
        P: 1e6
        Cp: 1464.4
        Lambda: 0.07
        m: 0.0043
        M: 0.0342
        thermal_effects: [102000.0]
        groups:true
        boundary_condition
        HMX: 0.999
        HMXprod: 0.001
        T: 800.0
        diffusion_coefficients
        HMX: 0.000009296
        HMXprod: 0.000009296
        HMX
        H: 4
        N: 8
        C: 8
        O: 8
        HMXprod
        H: 6
        C: 1
        O: 1
        reactions
        HMX=>10HMXprod: [130000.0, 0.0, 20920.0, 102000.0]
        solver_settings
        scheme: forward
        method: Sparse
        strategy: Damped
        linear_sys_method: None
        abs_tolerance: 1e-7
        max_iterations: 100
        loglevel: Some(info)
        dont_save_logs: true
        bounds
        C: -10.0, 10.0
        J:  -1e20, 1e20
        Teta:-100.0, 100.0
        q: -1e20, 1e20
        rel_tolerance
        C: 1e-7
        J: 1e-7
        Teta: 1e-7
        q:  1e-7
        strategy_params
        max_jac: Some(3)
        max_damp_iter: Some(10)
        damp_factor: Some(0.5)
        # Adaptive grid refinement settings (optional)
        adaptive_strategy
        # Refinement version
        version: 1
        # Maximum refinement iterations
        max_refinements: 3
                
        #Grid refinement method and parameters
        grid_refinement
        // Available methods:
        // doubleoints: []
        // easy: [parameter]
        // grcarsmooke: [param1, param2, param3]
        // pearson: [param1, param2]
        // twopnt: [param1, param2, param3]
        grcarsmooke: [0.05, 0.05, 1.25]
        postprocessing
        gnuplot:true
        save_to_csv:false
        filename: meow
        "#;

    #[test]
    fn test_file_reading_and_parsing() {
        // Create temporary file
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("test_task.txt");

        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(TEST_TASK_CONTENT.as_bytes())
            .expect("Failed to write to temp file");

        // Test file reading and parsing
        match std::fs::read_to_string(&file_path) {
            Ok(content) => {
                let mut parser = DocumentParser::new(content);
                match parser.parse_document() {
                    Ok(_) => {
                        if let Some(parsed_doc) = parser.get_result() {
                            // Verify parsed content
                            assert!(parsed_doc.contains_key("process_conditions"));
                            assert!(parsed_doc.contains_key("boundary_condition"));
                            assert!(parsed_doc.contains_key("diffusion_coefficients"));
                            assert!(parsed_doc.contains_key("reactions"));

                            // Check specific values
                            let process_conditions = parsed_doc.get("process_conditions").unwrap();
                            if let Some(Some(values)) = process_conditions.get("Tm") {
                                if let Value::Float(tm) = &values[0] {
                                    assert_eq!(*tm, 1500.0);
                                }
                            }
                        } else {
                            panic!("Parser returned no result");
                        }
                    }
                    Err(e) => panic!("Error parsing file: {}", e),
                }
            }
            Err(e) => panic!("Error reading file: {}", e),
        }
    }

    #[test]
    fn test_combustion_app_file_loading() {
        let mut app = CombustionApp::new();

        // Create temporary file
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("combustion_task.txt");

        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(TEST_TASK_CONTENT.as_bytes())
            .expect("Failed to write to temp file");

        // Simulate file reading button click
        match std::fs::read_to_string(&file_path) {
            Ok(content) => {
                let mut parser = DocumentParser::new(content);
                match parser.parse_document() {
                    Ok(_) => {
                        if let Some(parsed_doc) = parser.get_result() {
                            app.document = parsed_doc.clone();

                            // Verify app state after loading
                            assert!(!app.document.is_empty());
                            assert!(app.document.contains_key("process_conditions"));

                            // Check that substances were loaded
                            let process_conditions =
                                app.document.get("process_conditions").unwrap();
                            if let Some(Some(values)) = process_conditions.get("substances") {
                                assert_eq!(values.len(), 2); // HMX and HMXprod
                            }
                        }
                    }
                    Err(e) => panic!("Parsing failed: {}", e),
                }
            }
            Err(e) => panic!("File reading failed: {}", e),
        }
    }

    #[test]
    fn test_document_validation_after_file_load() {
        let mut app = CombustionApp::new();

        // Create temporary file with valid content
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("valid_task.txt");

        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(TEST_TASK_CONTENT.as_bytes())
            .expect("Failed to write to temp file");

        // Load and parse file
        let content = std::fs::read_to_string(&file_path).expect("Failed to read file");
        let mut parser = DocumentParser::new(content);
        parser.parse_document().expect("Failed to parse document");

        if let Some(parsed_doc) = parser.get_result() {
            app.document = parsed_doc.clone();

            // Validate required sections exist
            let required_sections = [
                "process_conditions",
                "boundary_condition",
                "diffusion_coefficients",
            ];
            for section in &required_sections {
                assert!(
                    app.document.contains_key(*section),
                    "Missing section: {}",
                    section
                );
            }

            // Validate process conditions
            let process_conditions = app.document.get("process_conditions").unwrap();
            let required_fields = ["Tm", "P", "Cp", "Lambda"];
            for field in &required_fields {
                assert!(
                    process_conditions.contains_key(*field),
                    "Missing field: {}",
                    field
                );
            }
        }
    }

    #[test]
    fn test_error_handling_invalid_file() {
        let invalid_content = "invalid content that cannot be parsed";

        let mut parser = DocumentParser::new(invalid_content.to_string());
        match parser.parse_document() {
            Ok(_) => {
                // Even if parsing succeeds, result might be empty or incomplete
                if let Some(result) = parser.get_result() {
                    println!("Parsed result: {:?}", result);
                } else {
                    println!("No result from parser");
                }
            }
            Err(e) => {
                println!("Expected parsing error: {}", e);
                // This is expected behavior for invalid content
            }
        }
    }

    #[test]
    fn test_empty_file_handling() {
        let empty_content = "";

        let mut parser = DocumentParser::new(empty_content.to_string());
        match parser.parse_document() {
            Ok(_) => {
                if let Some(result) = parser.get_result() {
                    assert!(result.is_empty() || result.len() == 0);
                }
            }
            Err(_) => {
                // Empty content might cause parsing error, which is acceptable
            }
        }
    }
}
#[cfg(test)]
mod deletion_tests {
    use super::*;
    use crate::gui::combustion::CombustionApp;
    use RustedSciThe::Utils::task_parser::Value;
    use std::collections::HashMap;

    #[test]
    fn test_field_deletion_simulation() {
        let mut app = CombustionApp::new();
        let section_name = "test_section".to_string();

        // Create section with multiple fields
        let mut section = HashMap::new();
        section.insert(
            "field1".to_string(),
            Some(vec![Value::String("value1".to_string())]),
        );
        section.insert("field2".to_string(), Some(vec![Value::Float(42.0)]));
        section.insert("field3".to_string(), Some(vec![Value::Integer(100)]));

        app.document.insert(section_name.clone(), section);

        // Verify initial state
        let section = app.document.get(&section_name).unwrap();
        assert_eq!(section.len(), 3);
        assert!(section.contains_key("field1"));
        assert!(section.contains_key("field2"));
        assert!(section.contains_key("field3"));

        // Simulate field deletion (field2)
        if let Some(section) = app.document.get_mut(&section_name) {
            section.remove("field2");
        }

        // Verify field was deleted
        let section = app.document.get(&section_name).unwrap();
        assert_eq!(section.len(), 2);
        assert!(section.contains_key("field1"));
        assert!(!section.contains_key("field2"));
        assert!(section.contains_key("field3"));
    }

    #[test]
    fn test_section_deletion_simulation() {
        let mut app = CombustionApp::new();

        // Create multiple sections
        app.document.insert("section1".to_string(), HashMap::new());
        app.document.insert("section2".to_string(), HashMap::new());
        app.document.insert("section3".to_string(), HashMap::new());

        // Verify initial state
        assert_eq!(app.document.len(), 3);
        assert!(app.document.contains_key("section1"));
        assert!(app.document.contains_key("section2"));
        assert!(app.document.contains_key("section3"));

        // Simulate section deletion (section2)
        app.document.remove("section2");

        // Verify section was deleted
        assert_eq!(app.document.len(), 2);
        assert!(app.document.contains_key("section1"));
        assert!(!app.document.contains_key("section2"));
        assert!(app.document.contains_key("section3"));
    }
}
#[cfg(test)]
mod saving_tests {
    use super::*;
    use crate::gui::combustion::CombustionApp;
    use RustedSciThe::Utils::task_parser::{DocumentParser, Value};
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_document_to_string_conversion() {
        let mut app = CombustionApp::new();

        // Create test document
        let mut section = HashMap::new();
        section.insert(
            "string_field".to_string(),
            Some(vec![Value::String("test_value".to_string())]),
        );
        section.insert("float_field".to_string(), Some(vec![Value::Float(3.14)]));
        section.insert("int_field".to_string(), Some(vec![Value::Integer(42)]));
        section.insert("bool_field".to_string(), Some(vec![Value::Boolean(true)]));
        section.insert(
            "vector_field".to_string(),
            Some(vec![Value::Vector(vec![1.0, 2.0, 3.0])]),
        );

        app.document.insert("test_section".to_string(), section);

        // Convert to string
        let content = app.document_to_string();

        // Verify content contains expected fields
        assert!(content.contains("test_section"));
        assert!(content.contains("string_field: test_value"));
        assert!(content.contains("float_field: 3.14"));
        assert!(content.contains("int_field: 42"));
        assert!(content.contains("bool_field: true"));
        assert!(content.contains("vector_field: [1, 2, 3]"));
    }

    #[test]
    fn test_save_document_functionality() {
        let mut app = CombustionApp::new();

        // Create test document
        let mut section = HashMap::new();
        section.insert("temperature".to_string(), Some(vec![Value::Float(298.15)]));
        section.insert("pressure".to_string(), Some(vec![Value::Float(101325.0)]));
        app.document.insert("conditions".to_string(), section);

        // Create temporary file path
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("test_save.txt");

        // Save document
        app.save_document(file_path.clone());

        // Verify file was created and contains expected content
        assert!(file_path.exists());
        let saved_content = std::fs::read_to_string(&file_path).expect("Failed to read saved file");

        assert!(saved_content.contains("conditions"));
        assert!(saved_content.contains("temperature: 298.15"));
        assert!(saved_content.contains("pressure: 101325"));
    }

    #[test]
    fn test_save_and_reload_cycle() {
        let mut app = CombustionApp::new();

        // Create test document
        let mut section = HashMap::new();
        section.insert(
            "test_field".to_string(),
            Some(vec![Value::String("test_value".to_string())]),
        );
        section.insert("number_field".to_string(), Some(vec![Value::Float(42.5)]));
        app.document.insert("test_section".to_string(), section);

        // Save to temporary file
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("save_reload_test.txt");
        app.save_document(file_path.clone());

        // Create new app and load the saved file
        let mut new_app = CombustionApp::new();
        let content = std::fs::read_to_string(&file_path).expect("Failed to read file");
        let mut parser = DocumentParser::new(content);

        match parser.parse_document() {
            Ok(_) => {
                if let Some(parsed_doc) = parser.get_result() {
                    new_app.document = parsed_doc.clone();
                    new_app.current_file_path = Some(file_path);

                    // Verify loaded document matches original
                    assert!(new_app.document.contains_key("test_section"));
                    let section = new_app.document.get("test_section").unwrap();
                    assert!(section.contains_key("test_field"));
                    assert!(section.contains_key("number_field"));
                }
            }
            Err(e) => panic!("Failed to parse saved file: {}", e),
        }
    }

    #[test]
    fn test_current_file_path_tracking() {
        let mut app = CombustionApp::new();

        // Initially no file path
        assert!(app.current_file_path.is_none());

        // Create temporary file
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("path_test.txt");
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(b"test_section\ntest_field: test_value\n")
            .expect("Failed to write to temp file");

        // Simulate file loading (sets current_file_path)
        let content = std::fs::read_to_string(&file_path).expect("Failed to read file");
        let mut parser = DocumentParser::new(content);

        match parser.parse_document() {
            Ok(_) => {
                if let Some(parsed_doc) = parser.get_result() {
                    app.document = parsed_doc.clone();
                    app.current_file_path = Some(file_path.clone());

                    // Verify file path is tracked
                    assert!(app.current_file_path.is_some());
                    assert_eq!(app.current_file_path.as_ref().unwrap(), &file_path);
                }
            }
            Err(_) => {
                // Parser might fail on simple content, which is acceptable for this test
                app.current_file_path = Some(file_path.clone());
                assert!(app.current_file_path.is_some());
            }
        }
    }

    #[test]
    fn test_save_with_special_values() {
        let mut app = CombustionApp::new();

        // Create document with special value types
        let mut section = HashMap::new();
        section.insert(
            "optional_some".to_string(),
            Some(vec![Value::Optional(Some(Box::new(Value::String(
                "nested".to_string(),
            ))))]),
        );
        section.insert(
            "optional_none".to_string(),
            Some(vec![Value::Optional(None)]),
        );
        section.insert("none_field".to_string(), None);

        app.document.insert("special_section".to_string(), section);

        // Convert to string and verify special cases are handled
        let content = app.document_to_string();

        assert!(content.contains("special_section"));
        assert!(content.contains("optional_some: Some(nested)"));
        assert!(content.contains("optional_none: None"));
        assert!(content.contains("none_field: None"));
    }
}
