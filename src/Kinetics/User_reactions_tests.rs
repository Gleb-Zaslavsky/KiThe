/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TESTS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use crate::Kinetics::User_reactions::KinData;
    use crate::Kinetics::mechfinder_api::ReactionType;
    use crate::Kinetics::mechfinder_api::kinetics::ElementaryStruct;
    use crate::Kinetics::mechfinder_api::{ReactionData, ReactionKinetics};
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;
    use serde_json::json;
    use std::collections::HashMap;
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

    #[test]
    fn test_kindata_workflow_shortcut_range() {
        let mut kd = KinData::new();

        // Test set_reactions_from_shortcut_range
        let shortcuts = kd.set_reactions_from_shortcut_range("C1..C3".to_string());
        assert_eq!(shortcuts, vec!["C_1", "C_2", "C_3"]);
        assert_eq!(
            kd.shortcut_reactions,
            Some(vec![
                "C_1".to_string(),
                "C_2".to_string(),
                "C_3".to_string()
            ])
        );

        // Test get_reactions_from_shortcuts
        kd.get_reactions_from_shortcuts();
        assert!(kd.map_of_reactions.is_some());
        let map = kd.map_of_reactions.as_ref().unwrap();
        assert!(map.contains_key("Cantera"));
        assert_eq!(map["Cantera"], vec!["1", "2", "3"]);

        // Test kinetic_main (combines parsing and analysis)
        kd.kinetic_main();
        assert!(kd.vec_of_reaction_data.is_some());
        assert!(!kd.vec_of_equations.is_empty());
        assert!(!kd.substances.is_empty());
        assert!(!kd.stecheodata.stecheo_matrx.is_empty());

        // Verify data consistency
        let reactions = kd.vec_of_reaction_data.as_ref().unwrap();
        assert_eq!(reactions.len(), kd.vec_of_equations.len());
        assert_eq!(
            kd.stecheodata.stecheo_matrx.len(),
            kd.vec_of_equations.len()
        );
    }

    #[test]
    fn test_kindata_calc_sym_constants_workflow() {
        let mut kd = KinData::new();

        // Setup with known reactions
        kd.set_reactions_from_shortcut_range("C1..C2".to_string());
        kd.get_reactions_from_shortcuts();
        kd.kinetic_main();
        let conc = HashMap::from([
            ("C_1".to_string(), Expr::Const(1.0)),
            ("C_2".to_string(), Expr::Const(2.0)),
        ]);
        // Test calc_sym_constants for all reactions
        kd.calc_sym_constants(None, Some(conc), None);
    }

    #[test]
    fn test_kindata_pretty_print_integration() {
        let mut kd = KinData::new();

        // Setup complete workflow
        kd.set_reactions_from_shortcut_range("C1..C2".to_string());
        kd.get_reactions_from_shortcuts();
        kd.kinetic_main();

        // Test pretty_print_kindata (should not panic)
        kd.pretty_print_kindata();

        // Verify all required data is present for pretty printing
        assert!(kd.vec_of_reaction_data.is_some());
        assert!(!kd.vec_of_equations.is_empty());
        assert!(!kd.substances.is_empty());
        assert!(kd.stecheodata.matrix_of_elements.is_some());
    }

    #[test]
    fn test_kindata_complete_workflow_with_constants() {
        let mut kd = KinData::new();

        // Complete workflow test
        kd.set_reactions_from_shortcut_range("C1..C3".to_string());
        kd.get_reactions_from_shortcuts();
        kd.kinetic_main();
        let conc = HashMap::from([("H2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
        // Test single reaction constant
        let k_result = kd
            .clone()
            .calc_K_const_for_1_react(0, 1000.0, None, Some(conc.clone()));
        assert!(k_result.is_ok());
        let k_value = k_result.unwrap();
        assert!(k_value > 0.0);
        // Test all reactions constants
        let all_k_result = kd.calc_K_const_for_all_reactions(1000.0, None, Some(conc), Some(false));
        assert!(all_k_result.is_ok());
        let all_k = all_k_result.unwrap();
        assert_eq!(all_k.len(), kd.vec_of_equations.len());
        let mut kd = kd.clone();
        // Test symbolic constants
        let conc = HashMap::from([
            ("C_1".to_string(), Expr::Const(1.0)),
            ("C_2".to_string(), Expr::Const(2.0)),
        ]);
        kd.calc_sym_constants(None, Some(conc), None);
        assert!(kd.K_sym_vec.is_some());
        let sym_vec = kd.K_sym_vec.as_ref().unwrap();
        assert_eq!(sym_vec.len(), kd.vec_of_equations.len());

        // Final pretty print test
        kd.pretty_print_kindata();
    }

    #[test]
    fn test_kindata_error_handling() {
        let mut kd = KinData::new();

        // Test calc_K_const_for_1_react with no data
        let result = kd.calc_K_const_for_1_react(0, 1000.0, None, None);
        assert!(result.is_err());

        // Test calc_sym_constants with no data (should not panic)
        kd.calc_sym_constants(None, None, None);
        assert!(kd.K_sym_vec.is_none());

        // Test calc_K_const_for_all_reactions with no data
        let all_result = kd.calc_K_const_for_all_reactions(1000.0, None, None, None);
        assert!(all_result.is_err());
    }

    #[test]
    fn test_kindata_different_shortcut_ranges() {
        let test_cases = vec![
            ("A1..A5", vec!["A_1", "A_2", "A_3", "A_4", "A_5"]),
            ("NUIG10..NUIG12", vec!["NUIG_10", "NUIG_11", "NUIG_12"]),
            ("Test1..Test1", vec!["Test_1"]),
        ];

        for (input, expected) in test_cases {
            let mut kd = KinData::new();
            let result = kd.set_reactions_from_shortcut_range(input.to_string());
            assert_eq!(result, expected, "Failed for input: {}", input);

            // Verify shortcut_reactions is set correctly
            let expected_shortcuts: Vec<String> = expected.iter().map(|s| s.to_string()).collect();
            assert_eq!(kd.shortcut_reactions, Some(expected_shortcuts));
        }
    }
}
