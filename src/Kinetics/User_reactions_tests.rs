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

    fn parse_kinetics_document(
        content: &str,
    ) -> HashMap<String, HashMap<String, serde_json::Value>> {
        let json_start = content
            .find('{')
            .expect("expected JSON payload in kinetics document");
        serde_json::from_str(&content[json_start..]).expect("expected valid kinetics JSON")
    }

    fn make_every_reaction(
        library: &str,
        reaction_name: &str,
        equation: &str,
        a: f64,
    ) -> crate::Kinetics::User_reactions::EveryReaction {
        crate::Kinetics::User_reactions::EveryReaction {
            reaction_id: crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: library.to_string(),
                reaction_id: reaction_name.to_string(),
            },
            shortcut: Some(format!("{}_{}", library, reaction_name)),
            lib_and_id: Some((library.to_string(), reaction_name.to_string())),
            reaction: ReactionData {
                reaction_type: ReactionType::Elem,
                eq: equation.to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![a, 0.0, 0.0],
                }),
            },
            equation: equation.to_string(),
            K_sym: Some(Expr::Const(a)),
        }
    }

    fn assert_normalized_lengths(kin_data: &KinData) {
        let len = kin_data
            .vec_of_reaction_data
            .as_ref()
            .map(|v| v.len())
            .unwrap();
        assert_eq!(kin_data.vec_of_equations.len(), len);
        assert_eq!(kin_data.reaction_ids.as_ref().map(|v| v.len()), Some(len));
        assert_eq!(
            kin_data.shortcut_reactions.as_ref().map(|v| v.len()),
            Some(len)
        );
        assert_eq!(kin_data.vec_of_pairs.as_ref().map(|v| v.len()), Some(len));
        if let Some(values) = kin_data.vec_of_reaction_Values.as_ref() {
            assert_eq!(values.len(), len);
        }
        assert_eq!(kin_data.K_sym_vec.as_ref().map(|v| v.len()), Some(len));
        assert_eq!(kin_data.every_reaction.as_ref().map(|v| v.len()), Some(len));
    }

    fn assert_raw_lengths_match(kin_data: &KinData) {
        if let Some(values) = kin_data.vec_of_reaction_Values.as_ref() {
            assert_eq!(
                kin_data.reaction_ids.as_ref().map(|ids| ids.len()),
                Some(values.len())
            );
        }
        if let Some(shortcuts) = kin_data.shortcut_reactions.as_ref() {
            assert_eq!(
                kin_data.reaction_ids.as_ref().map(|ids| ids.len()),
                Some(shortcuts.len())
            );
        }
        if let Some(pairs) = kin_data.vec_of_pairs.as_ref() {
            assert_eq!(
                kin_data.reaction_ids.as_ref().map(|ids| ids.len()),
                Some(pairs.len())
            );
        }
        if !kin_data.vec_of_equations.is_empty() {
            assert_eq!(
                kin_data.reaction_ids.as_ref().map(|ids| ids.len()),
                Some(kin_data.vec_of_equations.len())
            );
        }
    }

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

        let _ = user_reactions.get_reactions_from_shortcuts();
        let derived_map = user_reactions.reaction_map();
        println!("map_of_reactions: {:?} \n \n ", &derived_map);
        assert_eq!(derived_map, Some(map_of_reactions));
        assert_eq!(
            user_reactions.state(),
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );
        assert_eq!(user_reactions.workflow_state(), user_reactions.state());
        assert!(user_reactions.reaction_map().is_some());

        user_reactions.reactdata_parsing().unwrap();
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

        user_reactions.analyze_reactions().unwrap();
        let subs = user_reactions.substances.clone();
        assert_eq!(subs.is_empty(), false);
        let s = user_reactions
            .stoichiometric_analyzer()
            .stecheo_matrx
            .clone();
        let nunber_of_reactions = s.len();
        let number_of_substances = s[0].len();
        assert_eq!(s.is_empty(), false);
        assert_eq!(nunber_of_reactions, shortcut_reactions.unwrap().len());
        assert_eq!(number_of_substances, subs.len());
        assert_eq!(
            user_reactions.stoichiometric_analyzer().reactions.len(),
            user_reactions.equations().len()
        );
        assert!(user_reactions.is_analyzed());
    }

    #[test]
    fn test_reactdata_parsing_requires_values() {
        let mut user_reactions = KinData::new();
        let err = user_reactions.reactdata_parsing().unwrap_err();
        assert!(
            err.to_string().contains("no reaction values available"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_equations_from_reactdata_requires_data() {
        let mut kin_data = KinData::new();
        let err = kin_data.equations_from_reactdata().unwrap_err();
        assert!(
            err.to_string().contains("no reaction data available"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_generate_strings_with_single_character_prefix() {
        let mut user_reactions = KinData::new();
        let result = user_reactions
            .set_reactions_from_shortcut_range("A1..5".to_string())
            .unwrap();
        let expected = vec!["A_1", "A_2", "A_3", "A_4", "A_5"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_generate_strings_with_multi_character_prefix() {
        let mut user_reactions = KinData::new();
        let result = user_reactions
            .set_reactions_from_shortcut_range("Cat5..8".to_string())
            .unwrap();
        let expected = vec!["Cat_5", "Cat_6", "Cat_7", "Cat_8"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_generate_strings_with_large_numbers() {
        let mut user_reactions = KinData::new();
        let result = user_reactions
            .set_reactions_from_shortcut_range("Meow20..25".to_string())
            .unwrap();
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
    fn test_calc_k_const_for_1_react_rejects_out_of_bounds_index() {
        let kin_data = KinData::new();
        let err = kin_data
            .calc_K_const_for_1_react(0, 1000.0, None, None)
            .unwrap_err();

        assert!(
            err.to_string()
                .contains("no vector of reaction data available")
        );
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
                    Arrenius: vec![1e12, 0.0, 20000.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "C -> D".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e13, 0.0, 20000.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A + B -> C".to_string(), "C -> D".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.vec_of_pairs = Some(vec![
            ("Library1".to_string(), "Reaction1".to_string()),
            ("Library1".to_string(), "Reaction2".to_string()),
        ]);
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);

        let result = kin_data.calc_K_const_for_all_reactions(1000.0, None, None, Some(true));
        assert!(result.is_ok());
        let k_values = result.unwrap();
        assert_eq!(k_values.len(), 2);
        let k_expected1 = 1e13 * f64::exp(-20000.0 / (R * 1000.0));
        let k_expected2 = 1e12 * f64::exp(-20000.0 / (R * 1000.0));
        assert_relative_eq!(k_values[0], k_expected1, epsilon = 1e0);
        assert_relative_eq!(k_values[1], k_expected2, epsilon = 1e0);
        let every_reaction = kin_data.every_reaction.as_ref().unwrap();
        assert_eq!(every_reaction[0].equation, "C -> D");
        assert_eq!(every_reaction[1].equation, "A + B -> C");
        assert_eq!(every_reaction[0].K_sym, Some(Expr::Const(2.0)));
        assert_eq!(every_reaction[1].K_sym, Some(Expr::Const(1.0)));
    }

    #[test]
    fn test_every_reaction_rebuild_prefers_canonical_ids_over_legacy_metadata() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibA".to_string(),
                reaction_id: "R1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibB".to_string(),
                reaction_id: "R2".to_string(),
            },
        ]);
        kin_data.shortcut_reactions = Some(vec!["STALE_1".to_string(), "STALE_2".to_string()]);
        kin_data.vec_of_pairs = Some(vec![
            ("Legacy".to_string(), "Old1".to_string()),
            ("Legacy".to_string(), "Old2".to_string()),
        ]);
        kin_data.K_sym_vec = Some(vec![Expr::Const(10.0), Expr::Const(20.0)]);

        let values = kin_data
            .calc_K_const_for_all_reactions(1000.0, None, None, Some(true))
            .unwrap();

        assert_eq!(values.len(), 2);
        let every_reaction = kin_data.every_reaction.as_ref().unwrap();
        assert_eq!(
            every_reaction[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibB".to_string(),
                reaction_id: "R2".to_string(),
            }
        );
        assert_eq!(
            every_reaction[1].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibA".to_string(),
                reaction_id: "R1".to_string(),
            }
        );
        assert!(
            every_reaction
                .iter()
                .all(|reaction| reaction.shortcut.is_none())
        );
        assert_eq!(
            every_reaction[0].lib_and_id,
            Some(("LibB".to_string(), "R2".to_string()))
        );
        assert_eq!(
            every_reaction[1].lib_and_id,
            Some(("LibA".to_string(), "R1".to_string()))
        );
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
    fn test_equations_from_reactdata_rebuilds_cache() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "H2 + O2 -> HO2".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "HO2 -> H + O2".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["stale equation".to_string()];

        kin_data.equations_from_reactdata().unwrap();

        assert_eq!(
            kin_data.vec_of_equations,
            vec!["H2 + O2 -> HO2".to_string(), "HO2 -> H + O2".to_string()]
        );
    }

    #[test]
    fn test_reactdata_parsing_refreshes_normalized_views_and_state() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_Values = Some(vec![
            json!({
                "type": "elem",
                "eq": "A + B -> C",
                "Arrenius": [1.0, 0.0, 0.0]
            }),
            json!({
                "type": "elem",
                "eq": "C -> D",
                "Arrenius": [2.0, 0.0, 0.0]
            }),
        ]);

        kin_data.reactdata_parsing().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
        assert_eq!(
            kin_data.vec_of_equations,
            vec!["A + B -> C".to_string(), "C -> D".to_string()]
        );
        assert_eq!(kin_data.reaction_count(), 2);
        assert!(kin_data.every_reaction.is_some());
        assert_eq!(
            kin_data.normalized_reactions().unwrap()[0].equation,
            "A + B -> C"
        );

        // Rebuild the equation cache from parsed data to verify the shared mutation tail stays coherent.
        kin_data.vec_of_equations = vec!["stale equation".to_string()];
        kin_data.equations_from_reactdata().unwrap();
        assert_eq!(
            kin_data.vec_of_equations,
            vec!["A + B -> C".to_string(), "C -> D".to_string()]
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
    }

    #[test]
    fn test_calc_sym_constants_refreshes_state_and_preserves_cache_alignment() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.vec_of_equations = vec!["A -> B".to_string()];
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
        ]);

        kin_data.calc_sym_constants(None, None, None).unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
        assert_eq!(kin_data.symbolic_constants().unwrap().len(), 1);
        assert!(kin_data.every_reaction.is_some());
        assert_eq!(
            kin_data.normalized_reactions().unwrap()[0].equation,
            "A -> B"
        );
    }

    #[test]
    fn test_append_reaction_invalidates_every_reaction() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_Values = Some(vec![json!({"data": "value1"})]);
        kin_data.shortcut_reactions = Some(vec!["R1".to_string()]);
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.vec_of_equations = vec!["A -> B".to_string()];
        kin_data.vec_of_pairs = Some(vec![("Lib1".to_string(), "Reaction1".to_string())]);
        kin_data.shortcut_reactions = Some(vec!["R1".to_string()]);
        kin_data.substances = vec!["A".to_string(), "B".to_string()];
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0)]);
        kin_data.every_reaction = Some(vec![crate::Kinetics::User_reactions::EveryReaction {
            reaction_id: crate::Kinetics::User_reactions::ReactionIdentity::Shortcut(
                "R1".to_string(),
            ),
            shortcut: Some("R1".to_string()),
            lib_and_id: Some(("Library1".to_string(), "Reaction1".to_string())),
            reaction: ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            equation: "A -> B".to_string(),
            K_sym: Some(Expr::Const(1.0)),
        }]);

        kin_data
            .append_reaction(vec![json!({"data": "value2"})])
            .unwrap();

        assert!(kin_data.every_reaction.is_none());
        assert!(kin_data.vec_of_reaction_data.is_none());
        assert!(kin_data.vec_of_pairs.is_none());
        assert!(kin_data.shortcut_reactions.is_none());
        assert!(kin_data.K_sym_vec.is_none());
        assert!(kin_data.substances.is_empty());
        assert_eq!(
            kin_data.vec_of_reaction_Values.as_ref().unwrap(),
            &vec![json!({"data": "value1"}), json!({"data": "value2"})]
        );
    }

    #[test]
    fn test_kindata_state_consistency_after_mutations() {
        let mut raw = KinData::builder()
            .direct_reactions(vec!["A -> B".to_string(), "B -> C".to_string()])
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            raw.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(raw.reaction_ids.as_ref().unwrap().len(), 2);
        assert_eq!(raw.vec_of_equations.len(), 2);
        assert!(raw.vec_of_reaction_data.is_none());
        assert!(raw.every_reaction.is_none());

        raw.append_reaction(vec![json!({"data": "value3"})])
            .unwrap();
        assert_eq!(
            raw.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(raw.vec_of_reaction_Values.as_ref().unwrap().len(), 1);
        assert_eq!(raw.reaction_ids.as_ref().unwrap().len(), 1);
        assert_eq!(
            raw.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "raw".to_string(),
                index: 0,
            }
        );
        assert!(
            raw.vec_of_reaction_Values
                .as_ref()
                .unwrap()
                .last()
                .is_some()
        );
        assert!(raw.vec_of_equations.is_empty());
        assert!(raw.vec_of_reaction_data.is_none());
        assert!(raw.every_reaction.is_none());

        let mut normalized = KinData::new();
        normalized.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        normalized.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        normalized.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        normalized.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        normalized.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string()),
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string()),
        ]);
        let conc = HashMap::from([
            ("A".to_string(), Expr::Const(1.0)),
            ("B".to_string(), Expr::Const(1.0)),
            ("C".to_string(), Expr::Const(1.0)),
        ]);
        normalized
            .calc_sym_constants(None, Some(conc), None)
            .unwrap();
        assert_normalized_lengths(&normalized);

        normalized.remove_by_index(0).unwrap();
        assert_normalized_lengths(&normalized);
        assert_eq!(normalized.vec_of_equations, vec!["B -> C".to_string()]);
    }

    #[test]
    fn test_state_refresh_uses_canonical_ids_and_richer_views() {
        let mut kin_data = KinData::new();
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string()),
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string()),
        ]);

        kin_data.refresh_state_from_available_data().unwrap();
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );

        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.shortcut_reactions = None;
        kin_data.vec_of_pairs = None;

        kin_data.refresh_state_from_available_data().unwrap();
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
    }

    #[test]
    fn test_state_refresh_prefers_shortcuts_over_direct_raw_payload() {
        let mut kin_data = KinData::new();
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string()),
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string()),
        ]);
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.vec_of_reaction_Values =
            Some(vec![json!({"data": "value1"}), json!({"data": "value2"})]);

        // Raw values are present, but the canonical shortcut ids must still win.
        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
    }

    #[test]
    fn test_state_refresh_prefers_library_pairs_over_direct_raw_payload() {
        let mut kin_data = KinData::new();
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib1".to_string(),
                reaction_id: "Reaction1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib2".to_string(),
                reaction_id: "Reaction2".to_string(),
            },
        ]);
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        kin_data.vec_of_reaction_Values =
            Some(vec![json!({"data": "value1"}), json!({"data": "value2"})]);

        // Raw values are present, but the canonical library ids must still win.
        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );
    }

    #[test]
    fn test_state_refresh_treats_document_only_ids_as_direct_branch() {
        let mut kin_data = KinData::new();
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "raw".to_string(),
                index: 0,
            },
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "append".to_string(),
                index: 1,
            },
        ]);

        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
    }

    #[test]
    fn test_state_refresh_treats_analysis_artifacts_as_analyzed() {
        let mut kin_data = KinData::new();
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0)]);

        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
    }

    #[test]
    fn test_state_refresh_rebuilds_missing_ids_from_metadata() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);

        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "Lib1".to_string(),
                    reaction_id: "Reaction1".to_string(),
                },
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "Lib2".to_string(),
                    reaction_id: "Reaction2".to_string(),
                },
            ]
        );
    }

    #[test]
    fn test_state_refresh_keeps_canonical_ids_even_when_legacy_metadata_disagrees() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.vec_of_equations = vec!["A -> B".to_string()];
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "canonical".to_string(),
                index: 0,
            },
        ]);
        kin_data.shortcut_reactions = Some(vec!["STALE_1".to_string()]);
        kin_data.vec_of_pairs = Some(vec![(
            "LegacyLib".to_string(),
            "LegacyReaction".to_string(),
        )]);

        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::Document {
                    source: "canonical".to_string(),
                    index: 0,
                }
            ]
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
    }

    #[test]
    fn test_state_refresh_marks_partially_normalized_data_as_parsed() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.vec_of_equations = vec!["A -> B".to_string()];
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "doc".to_string(),
                index: 0,
            },
        ]);
        // Canonical ids and reaction data are coherent, but no normalized export view exists yet.
        kin_data.shortcut_reactions = Some(vec![]);

        kin_data.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
    }

    #[test]
    fn test_reaction_ids_are_rebuilt_from_metadata_when_missing() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);

        let conc = HashMap::from([
            ("A".to_string(), Expr::Const(1.0)),
            ("B".to_string(), Expr::Const(1.0)),
            ("C".to_string(), Expr::Const(1.0)),
        ]);
        kin_data.calc_sym_constants(None, Some(conc), None).unwrap();

        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "Lib1".to_string(),
                    reaction_id: "Reaction1".to_string(),
                },
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "Lib2".to_string(),
                    reaction_id: "Reaction2".to_string(),
                },
            ]
        );
        assert_eq!(
            kin_data.every_reaction.as_ref().unwrap()[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib1".to_string(),
                reaction_id: "Reaction1".to_string(),
            }
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
    }

    #[test]
    fn test_reaction_ids_keep_canonical_shortcuts_over_legacy_metadata() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string()),
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string()),
        ]);

        // Canonical ids should survive even when older metadata points elsewhere.
        let conc = HashMap::from([
            ("A".to_string(), Expr::Const(1.0)),
            ("B".to_string(), Expr::Const(1.0)),
            ("C".to_string(), Expr::Const(1.0)),
        ]);
        kin_data.calc_sym_constants(None, Some(conc), None).unwrap();

        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string()),
                crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string()),
            ]
        );
        assert_eq!(
            kin_data.every_reaction.as_ref().unwrap()[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string())
        );
        assert_eq!(
            kin_data.every_reaction.as_ref().unwrap()[1].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string())
        );
    }

    #[test]
    fn test_append_reaction_with_shortcut_initializes_empty_state() {
        let mut kin_data = KinData::new();

        kin_data
            .append_reaction_with_shortcut(vec![json!({"data": "value1"})], vec!["R1".to_string()])
            .unwrap();

        assert_eq!(
            kin_data.vec_of_reaction_Values.as_ref().unwrap(),
            &vec![json!({"data": "value1"})]
        );
        assert_eq!(
            kin_data.shortcut_reactions.as_ref().unwrap(),
            &vec!["R1".to_string()]
        );
        assert!(kin_data.every_reaction.is_none());
    }

    #[test]
    fn test_append_reaction_with_shortcut_replaces_stale_branch_state() {
        let mut kin_data = KinData::new();
        kin_data.shortcut_reactions = Some(vec!["OLD".to_string()]);
        kin_data.vec_of_pairs = Some(vec![("Legacy".to_string(), "1".to_string())]);
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "OLD -> OLD".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.vec_of_equations = vec!["OLD -> OLD".to_string()];
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("OLD".to_string()),
        ]);
        kin_data.substances = vec!["OLD".to_string()];
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0)]);
        kin_data.every_reaction = Some(vec![crate::Kinetics::User_reactions::EveryReaction {
            reaction_id: crate::Kinetics::User_reactions::ReactionIdentity::Shortcut(
                "OLD".to_string(),
            ),
            shortcut: Some("OLD".to_string()),
            lib_and_id: Some(("Legacy".to_string(), "1".to_string())),
            reaction: ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "OLD -> OLD".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            equation: "OLD -> OLD".to_string(),
            K_sym: Some(Expr::Const(1.0)),
        }]);

        kin_data
            .append_reaction_with_shortcut(vec![json!({"data": "value1"})], vec!["R1".to_string()])
            .unwrap();

        assert_eq!(kin_data.shortcut_reactions, Some(vec!["R1".to_string()]));
        assert_eq!(kin_data.shortcut_names(), Some(vec!["R1".to_string()]));
        assert!(kin_data.vec_of_pairs.is_none());
        assert!(kin_data.vec_of_reaction_data.is_none());
        assert!(kin_data.vec_of_equations.is_empty());
        assert!(kin_data.substances.is_empty());
        assert!(kin_data.K_sym_vec.is_none());
        assert!(kin_data.every_reaction.is_none());
    }

    #[test]
    fn test_append_and_remove_shortcut_branch_keep_facade_in_sync() {
        let mut kin_data = KinData::new();
        kin_data
            .append_reaction_with_shortcut(
                vec![json!({"data": "value1"}), json!({"data": "value2"})],
                vec!["R1".to_string(), "R2".to_string()],
            )
            .unwrap();

        assert_eq!(
            kin_data.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert!(kin_data.has_shortcuts());
        assert_eq!(
            kin_data.shortcut_names(),
            Some(vec!["R1".to_string(), "R2".to_string()])
        );
        assert_eq!(kin_data.reaction_values().unwrap().len(), 2);
        assert_eq!(kin_data.reaction_count(), 2);

        let err = kin_data.remove_by_index(0).unwrap_err();
        assert!(err.to_string().contains("out of bounds"));
        assert_eq!(
            kin_data.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(
            kin_data.shortcut_names(),
            Some(vec!["R1".to_string(), "R2".to_string()])
        );
        assert_eq!(kin_data.reaction_count(), 2);
    }

    #[test]
    fn test_append_reaction_with_shortcut_rejects_length_mismatch() {
        let mut kin_data = KinData::new();
        let err = kin_data
            .append_reaction_with_shortcut(
                vec![json!({"data": "value1"})],
                vec!["R1".to_string(), "R2".to_string()],
            )
            .unwrap_err();

        assert!(err.to_string().contains("append_reaction_with_shortcut"));
    }

    #[test]
    fn test_append_reaction_from_map_appends_all_entries() {
        let mut kin_data = KinData::new();
        kin_data
            .append_reaction_from_map(vec![
                HashMap::from([(String::from("A"), vec![1.0, 2.0])]),
                HashMap::from([(String::from("B"), vec![3.0, 4.0])]),
            ])
            .unwrap();

        assert_eq!(kin_data.vec_of_reaction_Values.as_ref().unwrap().len(), 2);
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert!(kin_data.every_reaction.is_none());
    }

    #[test]
    fn test_append_reaction_extends_canonical_ids() {
        let mut kin_data = KinData::new();

        kin_data
            .append_reaction(vec![json!({"data": "value1"})])
            .unwrap();
        kin_data
            .append_reaction(vec![json!({"data": "value2"})])
            .unwrap();

        assert_eq!(
            kin_data.vec_of_reaction_Values.as_ref().unwrap(),
            &vec![json!({"data": "value1"}), json!({"data": "value2"})]
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::Document {
                    source: "append".to_string(),
                    index: 0,
                },
                crate::Kinetics::User_reactions::ReactionIdentity::Document {
                    source: "append".to_string(),
                    index: 1,
                },
            ]
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
    }

    #[test]
    fn test_append_reaction_from_direct_branch_resets_document_source() {
        let mut kin_data = KinData::from_direct_reactions(vec!["A -> B".to_string()]).unwrap();

        kin_data
            .append_reaction(vec![json!({"data": "value1"})])
            .unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(
            kin_data.vec_of_reaction_Values.as_ref().unwrap(),
            &vec![json!({"data": "value1"})]
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::Document {
                    source: "raw".to_string(),
                    index: 0,
                }
            ]
        );
        assert!(kin_data.shortcut_reactions.is_none());
    }

    #[test]
    fn test_append_reaction_from_shortcut_branch_starts_fresh_document_source() {
        let mut kin_data = KinData::new();
        kin_data
            .append_reaction_with_shortcut(vec![json!({"data": "value1"})], vec!["R1".to_string()])
            .unwrap();

        kin_data
            .append_reaction(vec![json!({"data": "value2"})])
            .unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(
            kin_data.vec_of_reaction_Values.as_ref().unwrap(),
            &vec![json!({"data": "value2"})]
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::Document {
                    source: "raw".to_string(),
                    index: 0,
                }
            ]
        );
        assert!(kin_data.shortcut_reactions.is_none());
    }

    #[test]
    fn test_save_raw_reactions_writes_json_file() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_Values = Some(vec![json!({"data": "value1"})]);

        let temp_file = NamedTempFile::new().unwrap();
        let file_base = temp_file.path().to_str().unwrap();
        kin_data.save_raw_reactions(file_base).unwrap();

        let saved_path = format!("{}.json", file_base);
        let content = fs::read_to_string(saved_path).unwrap();
        assert!(content.contains("value1"));
    }

    #[test]
    fn test_save_reactions_with_shortcuts_requires_shortcuts() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_Values = Some(vec![json!({"data": "value1"})]);
        let temp_file = NamedTempFile::new().unwrap();
        let err = kin_data
            .save_reactions_with_shortcuts(temp_file.path().to_str().unwrap())
            .unwrap_err();

        assert!(err.to_string().contains("no vector of shortcuts available"));
    }

    #[test]
    fn test_save_reactions_with_shortcuts_requires_reaction_values() {
        let kin_data = KinData::new();
        let temp_file = NamedTempFile::new().unwrap();
        let err = kin_data
            .save_reactions_with_shortcuts(temp_file.path().to_str().unwrap())
            .unwrap_err();

        assert!(err.to_string().contains("no reaction values available"));
    }

    #[test]
    fn test_load_reactions_from_json_initializes_ids_and_state() {
        let mut file = NamedTempFile::new().unwrap();
        let payload = json!({
            "R1": {"data": "value1"},
            "R2": {"data": "value2"}
        });
        serde_json::to_writer_pretty(file.as_file_mut(), &payload).unwrap();

        let mut kin_data = KinData::new();
        kin_data
            .load_reactions_from_json(file.path().to_str().unwrap())
            .unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(kin_data.shortcut_reactions.as_ref().unwrap().len(), 2);
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string())
        );
        assert_eq!(kin_data.vec_of_reaction_Values.as_ref().unwrap().len(), 2);
    }

    #[test]
    fn test_load_reactions_from_json_sorts_shortcuts() {
        let mut file = NamedTempFile::new().unwrap();
        let payload = json!({
            "R2": {"data": "value2"},
            "R1": {"data": "value1"}
        });
        serde_json::to_writer_pretty(file.as_file_mut(), &payload).unwrap();

        let mut kin_data = KinData::new();
        kin_data
            .load_reactions_from_json(file.path().to_str().unwrap())
            .unwrap();

        assert_eq!(
            kin_data.shortcut_reactions,
            Some(vec!["R1".to_string(), "R2".to_string()])
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string())
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap()[1],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string())
        );
    }

    #[test]
    fn test_remove_by_index_updates_all_reaction_views() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.vec_of_reaction_Values = Some(vec![json!({"data": "v1"}), json!({"data": "v2"})]);
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);
        kin_data.remove_by_index(0).unwrap();

        assert_eq!(kin_data.vec_of_equations, vec!["B -> C".to_string()]);
        assert_eq!(
            kin_data.shortcut_reactions.as_ref().unwrap(),
            &vec!["R2".to_string()]
        );
        assert_eq!(
            kin_data.vec_of_pairs.as_ref().unwrap(),
            &vec![("Lib2".to_string(), "Reaction2".to_string())]
        );
        assert_eq!(
            kin_data.K_sym_vec.as_ref().unwrap(),
            &vec![Expr::Const(2.0)]
        );
        let every_reaction = kin_data.every_reaction.as_ref().unwrap();
        assert_eq!(every_reaction.len(), 1);
        assert_eq!(every_reaction[0].equation, "B -> C");
        assert!(every_reaction[0].shortcut.is_none());
        assert_eq!(
            every_reaction[0].lib_and_id,
            Some(("Lib2".to_string(), "Reaction2".to_string()))
        );
        assert_eq!(
            every_reaction[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib2".to_string(),
                reaction_id: "Reaction2".to_string()
            }
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::Sorted
        );
    }

    #[test]
    fn test_remove_by_index_preserves_library_state_without_parsed_data() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib1".to_string(),
                reaction_id: "Reaction1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib2".to_string(),
                reaction_id: "Reaction2".to_string(),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);

        kin_data.remove_by_index(0).unwrap();

        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );
        assert_eq!(
            kin_data.vec_of_pairs.as_ref().unwrap(),
            &vec![("Lib2".to_string(), "Reaction2".to_string())]
        );
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "Lib2".to_string(),
                    reaction_id: "Reaction2".to_string(),
                }
            ]
        );
    }

    #[test]
    fn test_remove_by_index_rejects_out_of_bounds_index() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_equations = vec!["A -> B".to_string()];

        let err = kin_data.remove_by_index(1).unwrap_err();
        assert!(
            err.to_string().contains("out of bounds"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_remove_reaction_by_eq_rejects_missing_equation() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_equations = vec!["A -> B".to_string()];

        let err = kin_data.remove_reaction_by_eq("B -> C").unwrap_err();
        assert!(
            err.to_string().contains("was not found"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_sorting_keeps_canonical_ids_in_sync() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "C -> D".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![10.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "C -> D".to_string()];
        kin_data.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "doc".to_string(),
                index: 0,
            },
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "doc".to_string(),
                index: 1,
            },
        ]);
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);

        let result = kin_data.calc_K_const_for_all_reactions(1000.0, None, None, Some(true));
        assert!(result.is_ok());
        let every_reaction = kin_data.every_reaction.as_ref().unwrap();
        assert_eq!(every_reaction[0].equation, "C -> D");
        assert_eq!(
            every_reaction[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "doc".to_string(),
                index: 1
            }
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::Sorted
        );
    }

    #[test]
    fn test_sorting_keeps_library_ids_in_sync() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "C -> D".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![10.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "C -> D".to_string()];
        kin_data.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib1".to_string(),
                reaction_id: "Reaction1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib2".to_string(),
                reaction_id: "Reaction2".to_string(),
            },
        ]);
        kin_data.vec_of_reaction_Values =
            Some(vec![json!({"data": "value1"}), json!({"data": "value2"})]);
        kin_data.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);

        let result = kin_data.calc_K_const_for_all_reactions(1000.0, None, None, Some(true));
        assert!(result.is_ok());
        let every_reaction = kin_data.every_reaction.as_ref().unwrap();
        assert_eq!(
            every_reaction[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib2".to_string(),
                reaction_id: "Reaction2".to_string(),
            }
        );
        assert_eq!(
            every_reaction[1].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Lib1".to_string(),
                reaction_id: "Reaction1".to_string(),
            }
        );
        assert_eq!(
            kin_data.vec_of_reaction_Values.as_ref().unwrap(),
            &vec![json!({"data": "value2"}), json!({"data": "value1"})]
        );
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::Sorted
        );
    }
    #[test]
    fn test_create_kinetics_document() {
        let mut kin_data = KinData::new();

        // Test with normalized every_reaction view
        let reaction_1 = make_every_reaction("Library1", "Reaction1", "A -> B", 1.0);
        let reaction_2 = make_every_reaction("Library1", "Reaction2", "B -> C", 2.0);
        let reaction_3 = make_every_reaction("Library2", "Reaction3", "C -> D", 3.0);
        kin_data.every_reaction = Some(vec![
            reaction_1.clone(),
            reaction_2.clone(),
            reaction_3.clone(),
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
        let saved = parse_kinetics_document(&content);
        assert_eq!(
            saved["Library1"]["Reaction1"],
            reaction_1.reaction.to_serde_value().unwrap()
        );
        assert_eq!(
            saved["Library1"]["Reaction2"],
            reaction_2.reaction.to_serde_value().unwrap()
        );
        assert_eq!(
            saved["Library2"]["Reaction3"],
            reaction_3.reaction.to_serde_value().unwrap()
        );

        // Test with a different normalized every_reaction ordering
        let reaction_4 = make_every_reaction("Library3", "Reaction4", "D -> E", 4.0);
        let reaction_5 = make_every_reaction("Library3", "Reaction5", "E -> F", 5.0);
        let reaction_6 = make_every_reaction("Library4", "Reaction6", "F -> G", 6.0);
        kin_data.every_reaction = Some(vec![
            reaction_4.clone(),
            reaction_5.clone(),
            reaction_6.clone(),
        ]);
        kin_data.vec_of_reaction_Values = Some(vec![
            json!({"data": "value4"}),
            json!({"data": "value5"}),
            json!({"data": "value6"}),
        ]);

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());

        let content = fs::read_to_string(file_path).unwrap();
        assert!(content.contains("KINETICS"));
        let saved = parse_kinetics_document(&content);
        assert_eq!(
            saved["Library3"]["Reaction4"],
            reaction_4.reaction.to_serde_value().unwrap()
        );
        assert_eq!(
            saved["Library3"]["Reaction5"],
            reaction_5.reaction.to_serde_value().unwrap()
        );
        assert_eq!(
            saved["Library4"]["Reaction6"],
            reaction_6.reaction.to_serde_value().unwrap()
        );

        // Test with existing file and header
        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();
        fs::write(file_path, "KINETICS\nSome existing content\n").unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());

        let content = fs::read_to_string(file_path).unwrap();
        assert_eq!(content.matches("KINETICS").count(), 1);
        assert!(content.contains("Some existing content"));
        let saved = parse_kinetics_document(&content);
        assert_eq!(
            saved["Library3"]["Reaction4"],
            reaction_4.reaction.to_serde_value().unwrap()
        );
        assert_eq!(
            saved["Library4"]["Reaction6"],
            reaction_6.reaction.to_serde_value().unwrap()
        );
    }

    #[test]
    fn test_create_kinetics_document_fills_missing_reaction_ids_from_raw_payload() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "D -> E".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![4.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "E -> F".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![5.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "F -> G".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![6.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_equations = vec![
            "D -> E".to_string(),
            "E -> F".to_string(),
            "F -> G".to_string(),
        ];
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 1,
            },
        ]);

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());
        let saved = parse_kinetics_document(&fs::read_to_string(file_path).unwrap());

        assert_eq!(
            saved["direct"]["0"],
            kin_data.vec_of_reaction_data.as_ref().unwrap()[0]
                .to_serde_value()
                .unwrap()
        );
        assert_eq!(
            saved["direct"]["1"],
            kin_data.vec_of_reaction_data.as_ref().unwrap()[1]
                .to_serde_value()
                .unwrap()
        );
        assert_eq!(
            saved["raw"]["2"],
            kin_data.vec_of_reaction_data.as_ref().unwrap()[2]
                .to_serde_value()
                .unwrap()
        );
    }

    #[test]
    fn test_refresh_state_keeps_canonical_ids_and_fills_missing_raw_slots() {
        let mut kin_data = KinData::new();
        kin_data.vec_of_reaction_Values =
            Some(vec![json!({"data": "value1"}), json!({"data": "value2"})]);
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string()),
        ]);

        kin_data.refresh_state_from_available_data().unwrap();

        let reaction_ids = kin_data.reaction_ids.as_ref().unwrap();
        assert_eq!(reaction_ids.len(), 2);
        assert_eq!(
            reaction_ids[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string())
        );
        assert_eq!(
            reaction_ids[1],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "raw".to_string(),
                index: 1
            }
        );
    }

    #[test]
    fn test_create_kinetics_document_rejects_missing_equations_in_normalized_view() {
        let mut kin_data = KinData::new();
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
        ]);
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.vec_of_reaction_Values = Some(vec![json!({"data": "value1"})]);

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(
            err,
            crate::Kinetics::error::KineticsError::LengthMismatch { .. }
        ));
        assert!(
            err.to_string()
                .contains("vec_of_equations vs vec_of_reaction_data")
        );
    }

    #[test]
    fn test_create_kinetics_document_rebuilds_legacy_views() {
        let mut kin_data = KinData::new();
        kin_data.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LegacyLib".to_string(),
                reaction_id: "1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LegacyLib".to_string(),
                reaction_id: "2".to_string(),
            },
        ]);
        kin_data.vec_of_pairs = Some(vec![
            ("LegacyLib".to_string(), "1".to_string()),
            ("LegacyLib".to_string(), "2".to_string()),
        ]);
        kin_data.shortcut_reactions =
            Some(vec!["LegacyLib_1".to_string(), "LegacyLib_2".to_string()]);
        kin_data.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kin_data.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kin_data.vec_of_reaction_Values = Some(vec![json!({"legacy": 1}), json!({"legacy": 2})]);
        kin_data.every_reaction = None;

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path().to_str().unwrap();

        // The export path should rebuild the normalized view instead of relying on stale caches.
        let result = kin_data.create_kinetics_document(file_path);
        assert!(result.is_ok());
        let saved = parse_kinetics_document(&fs::read_to_string(file_path).unwrap());

        assert_eq!(
            saved["LegacyLib"]["1"],
            kin_data.vec_of_reaction_data.as_ref().unwrap()[0]
                .to_serde_value()
                .unwrap()
        );
        assert_eq!(
            saved["LegacyLib"]["2"],
            kin_data.vec_of_reaction_data.as_ref().unwrap()[1]
                .to_serde_value()
                .unwrap()
        );
    }

    #[test]
    fn test_kindata_workflow_shortcut_range() {
        let mut kd = KinData::new();

        // Test set_reactions_from_shortcut_range
        let shortcuts = kd
            .set_reactions_from_shortcut_range("C1..C3".to_string())
            .unwrap();
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
        let _ = kd.get_reactions_from_shortcuts();
        assert!(kd.reaction_map().is_some());
        let map = kd.reaction_map().unwrap();
        assert!(map.contains_key("Cantera"));
        assert_eq!(map["Cantera"], vec!["1", "2", "3"]);

        // Test kinetic_main (combines parsing and analysis)
        let _ = kd.kinetic_main();
        assert!(kd.vec_of_reaction_data.is_some());
        assert!(!kd.vec_of_equations.is_empty());
        assert!(!kd.substances.is_empty());
        assert!(!kd.stecheodata.stecheo_matrx.is_empty());
        assert!(kd.every_reaction.is_some());

        // Verify data consistency
        let reactions = kd.vec_of_reaction_data.as_ref().unwrap();
        assert_eq!(reactions.len(), kd.vec_of_equations.len());
        assert_eq!(
            kd.stecheodata.stecheo_matrx.len(),
            kd.vec_of_equations.len()
        );
        let every_reaction = kd.every_reaction.as_ref().unwrap();
        assert_eq!(every_reaction.len(), kd.vec_of_equations.len());
        assert_eq!(every_reaction[0].equation, kd.vec_of_equations[0]);
        assert_eq!(every_reaction[0].reaction.eq, kd.vec_of_equations[0]);
    }

    #[test]
    fn test_kindata_calc_sym_constants_workflow() {
        let mut kd = KinData::new();

        // Setup with known reactions
        kd.set_reactions_from_shortcut_range("C1..C2".to_string())
            .unwrap();
        let _ = kd.get_reactions_from_shortcuts();
        let _ = kd.kinetic_main();
        let conc = HashMap::from([
            ("C_1".to_string(), Expr::Const(1.0)),
            ("C_2".to_string(), Expr::Const(2.0)),
        ]);
        // Test calc_sym_constants for all reactions
        kd.calc_sym_constants(None, Some(conc), None).unwrap();
        assert!(kd.K_sym_vec.is_some());
        assert!(kd.every_reaction.is_some());
        assert_eq!(
            kd.every_reaction.as_ref().unwrap().len(),
            kd.K_sym_vec.as_ref().unwrap().len()
        );
        assert!(kd.every_reaction.as_ref().unwrap()[0].K_sym.is_some());
    }

    #[test]
    fn test_kindata_pretty_print_integration() {
        let mut kd = KinData::new();

        // Setup complete workflow
        kd.set_reactions_from_shortcut_range("C1..C2".to_string())
            .unwrap();
        let _ = kd.get_reactions_from_shortcuts();
        let _ = kd.kinetic_main();

        // Test pretty_print_kindata (should not panic)
        assert!(kd.pretty_print_kindata().is_ok());

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
        kd.set_reactions_from_shortcut_range("C1..C3".to_string())
            .unwrap();
        let _ = kd.get_reactions_from_shortcuts();
        let _ = kd.kinetic_main();
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
        kd.calc_sym_constants(None, Some(conc), None).unwrap();
        assert!(kd.K_sym_vec.is_some());
        let sym_vec = kd.K_sym_vec.as_ref().unwrap();
        assert_eq!(sym_vec.len(), kd.vec_of_equations.len());

        // Final pretty print test
        assert!(kd.pretty_print_kindata().is_ok());
    }

    #[test]
    fn test_kindata_error_handling() {
        let mut kd = KinData::new();

        // Test calc_K_const_for_1_react with no data
        let result = kd.calc_K_const_for_1_react(0, 1000.0, None, None);
        assert!(result.is_err());

        // Test calc_sym_constants with no data (should not panic)
        let err = kd.calc_sym_constants(None, None, None).unwrap_err();
        assert!(err.to_string().contains("no reaction data available"));
        assert!(kd.K_sym_vec.is_none());

        let analyze_err = kd.analyze_reactions().unwrap_err();
        assert!(matches!(
            analyze_err,
            crate::Kinetics::error::KineticsError::InvalidState(message)
            if message.contains("reaction analysis")
        ));

        let shortcut_err = kd.get_reactions_from_shortcuts().unwrap_err();
        assert!(
            shortcut_err
                .to_string()
                .contains("no shortcut reactions available")
        );

        // Test calc_K_const_for_all_reactions with no data
        let all_result = kd.calc_K_const_for_all_reactions(1000.0, None, None, None);
        assert!(all_result.is_err());

        // Pretty printing should now return an error instead of panicking.
        assert!(kd.pretty_print_kindata().is_err());
    }

    #[test]
    fn test_load_reactions_from_json_rejects_invalid_json() {
        let mut file = NamedTempFile::new().unwrap();
        std::io::Write::write_all(file.as_file_mut(), b"not-json").unwrap();

        let mut kin_data = KinData::new();
        let err = kin_data
            .load_reactions_from_json(file.path().to_str().unwrap())
            .unwrap_err();
        assert!(err.to_string().contains("JSON error"));
    }

    #[test]
    fn test_load_reactions_from_json_replaces_stale_state() {
        let mut file = NamedTempFile::new().unwrap();
        let payload = json!({
            "R2": {"data": "value2"},
            "R1": {"data": "value1"}
        });
        serde_json::to_writer_pretty(file.as_file_mut(), &payload).unwrap();

        let mut kin_data = KinData::new();
        kin_data.shortcut_reactions = Some(vec!["legacy".to_string()]);
        kin_data.vec_of_pairs = Some(vec![("legacy".to_string(), "1".to_string())]);
        kin_data.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "legacy -> data".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kin_data.every_reaction = Some(vec![make_every_reaction(
            "legacy",
            "1",
            "legacy -> data",
            1.0,
        )]);

        kin_data
            .load_reactions_from_json(file.path().to_str().unwrap())
            .unwrap();

        // Loading raw JSON should replace stale cache-like fields instead of extending them.
        assert_eq!(
            kin_data.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(
            kin_data.shortcut_reactions,
            Some(vec!["R1".to_string(), "R2".to_string()])
        );
        assert!(kin_data.vec_of_pairs.is_none());
        assert!(kin_data.every_reaction.is_none());
        assert_eq!(kin_data.vec_of_reaction_Values.as_ref().unwrap().len(), 2);
        assert_eq!(
            kin_data.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R1".to_string())
        );
    }

    #[test]
    fn test_kindata_builder_creates_direct_state() {
        let kd = KinData::builder()
            .direct_reactions(vec!["A -> B".to_string(), "B -> C".to_string()])
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(kd.reaction_ids.as_ref().unwrap().len(), 2);
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0
            }
        );
        assert_eq!(
            kd.vec_of_equations,
            vec!["A -> B".to_string(), "B -> C".to_string()]
        );
        assert!(kd.shortcut_reactions.is_none());
    }

    #[test]
    fn test_kindata_builder_last_write_wins() {
        let kd = KinData::builder()
            .direct_reactions(vec!["A -> B".to_string()])
            .unwrap()
            .shortcut_reactions(vec!["C_1".to_string(), "C_2".to_string()])
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(
            kd.shortcut_reactions,
            Some(vec!["C_1".to_string(), "C_2".to_string()])
        );
        assert!(kd.vec_of_equations.is_empty());
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("C_1".to_string())
        );
    }

    #[test]
    fn test_kindata_builder_library_reactions_resets_previous_branch() {
        let kd = KinData::builder()
            .direct_reactions(vec!["A -> B".to_string()])
            .unwrap()
            .library_reactions(
                "LibX".to_string(),
                vec!["R1".to_string(), "R2".to_string(), "R3".to_string()],
            )
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );
        assert!(kd.vec_of_equations.is_empty());
        assert!(kd.vec_of_reaction_data.is_none());
        assert!(kd.vec_of_reaction_Values.is_none());
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "LibX".to_string(),
                    reaction_id: "R1".to_string(),
                },
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "LibX".to_string(),
                    reaction_id: "R2".to_string(),
                },
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "LibX".to_string(),
                    reaction_id: "R3".to_string(),
                },
            ]
        );
        assert_eq!(
            kd.vec_of_pairs.as_ref().unwrap(),
            &vec![
                ("LibX".to_string(), "R1".to_string()),
                ("LibX".to_string(), "R2".to_string()),
                ("LibX".to_string(), "R3".to_string()),
            ]
        );
    }

    #[test]
    fn test_direct_reactions_workflow_keeps_equations_for_analysis() {
        let mut kd = KinData::new();
        kd.set_reactions_directly(
            vec!["H2 + O2 -> H2O2".to_string(), "H2O2 -> H2O + O".to_string()],
            None,
        )
        .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(
            kd.vec_of_equations,
            vec!["H2 + O2 -> H2O2".to_string(), "H2O2 -> H2O + O".to_string()]
        );
        assert!(kd.shortcut_reactions.is_none());

        kd.kinetic_main().unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
        assert_eq!(kd.vec_of_equations.len(), 2);
        assert!(!kd.substances.is_empty());
    }

    #[test]
    fn test_set_reactions_directly_preserves_groups() {
        let mut kd = KinData::new();
        let groups = HashMap::from([(
            "Me".to_string(),
            HashMap::from([("C".to_string(), 1usize), ("H".to_string(), 3usize)]),
        )]);

        kd.set_reactions_directly(
            vec!["Me + O2 -> CO2 + H2O".to_string()],
            Some(groups.clone()),
        )
        .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(kd.groups, Some(groups));
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0
            }
        );
    }

    #[test]
    fn test_set_reactions_directly_resets_stale_state_before_installing_new_reactions() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("OLD".to_string()),
        ]);
        kd.shortcut_reactions = Some(vec!["OLD".to_string()]);
        kd.vec_of_pairs = Some(vec![("Legacy".to_string(), "1".to_string())]);
        kd.vec_of_reaction_Values = Some(vec![json!({"old": true})]);
        kd.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "OLD -> OLD".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kd.vec_of_equations = vec!["OLD -> OLD".to_string()];
        kd.substances = vec!["OLD".to_string()];
        kd.K_sym_vec = Some(vec![Expr::Const(1.0)]);
        kd.every_reaction = Some(vec![make_every_reaction("Legacy", "1", "OLD -> OLD", 1.0)]);
        kd.state = crate::Kinetics::User_reactions::KinDataState::Sorted;

        kd.set_reactions_directly(vec!["A -> B".to_string()], None)
            .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0
            }
        );
        assert_eq!(kd.vec_of_equations, vec!["A -> B".to_string()]);
        assert_eq!(kd.shortcut_reactions, None);
        assert!(kd.vec_of_pairs.is_none());
        assert!(kd.vec_of_reaction_Values.is_none());
        assert!(kd.vec_of_reaction_data.is_none());
        assert!(kd.substances.is_empty());
        assert!(kd.K_sym_vec.is_none());
        assert!(kd.every_reaction.is_none());
    }

    #[test]
    fn test_shortcut_range_resets_stale_branch_data_before_installing_shortcuts() {
        let mut kd = KinData::new();
        kd.vec_of_reaction_Values = Some(vec![json!({"stale": true})]);
        kd.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "STALE -> DATA".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kd.vec_of_equations = vec!["STALE -> DATA".to_string()];
        kd.vec_of_pairs = Some(vec![("Legacy".to_string(), "1".to_string())]);
        kd.shortcut_reactions = Some(vec!["OLD".to_string()]);
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "old".to_string(),
                index: 0,
            },
        ]);
        kd.every_reaction = Some(vec![make_every_reaction(
            "Legacy",
            "1",
            "STALE -> DATA",
            1.0,
        )]);
        kd.state = crate::Kinetics::User_reactions::KinDataState::Analyzed;

        let shortcuts = kd
            .set_reactions_from_shortcut_range("C1..C2".to_string())
            .unwrap();

        assert_eq!(shortcuts, vec!["C_1".to_string(), "C_2".to_string()]);
        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(
            kd.shortcut_reactions,
            Some(vec!["C_1".to_string(), "C_2".to_string()])
        );
        assert!(kd.vec_of_reaction_Values.is_none());
        assert!(kd.vec_of_reaction_data.is_none());
        assert!(kd.vec_of_pairs.is_none());
        assert!(kd.vec_of_equations.is_empty());
        assert!(kd.every_reaction.is_none());
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("C_1".to_string())
        );
    }

    #[test]
    fn test_validate_state_contract_rejects_sorted_without_normalized_view() {
        let kd = KinData::new();
        let err = kd
            .validate_state_contract(crate::Kinetics::User_reactions::KinDataState::Sorted)
            .unwrap_err();

        assert!(
            err.to_string()
                .contains("coherent normalized reaction view")
        );
    }

    #[test]
    fn test_validate_state_contract_rejects_empty_direct_branch() {
        let kd = KinData::new();
        let err = kd
            .validate_state_contract(
                crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded,
            )
            .unwrap_err();

        assert!(
            err.to_string()
                .contains("raw reaction values or document ids")
        );
    }

    #[test]
    fn test_validate_state_contract_rejects_analyzed_without_analysis_cache() {
        let mut kd = KinData::new();
        kd.state = crate::Kinetics::User_reactions::KinDataState::Analyzed;

        let err = kd
            .validate_state_contract(crate::Kinetics::User_reactions::KinDataState::Analyzed)
            .unwrap_err();

        assert!(err.to_string().contains("Analyzed requires"));
        assert!(err.to_string().contains("analyzed artifacts"));
    }

    #[test]
    fn test_kindata_builder_supports_shortcut_range() {
        let kd = KinData::builder()
            .with_shortcut_range("C1..C2".to_string())
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(
            kd.shortcut_names(),
            Some(vec!["C_1".to_string(), "C_2".to_string()])
        );
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("C_1".to_string())
        );
    }

    #[test]
    fn test_kindata_facade_names_track_current_branch() {
        let kd = KinData::builder()
            .with_shortcut_range("C1..C2".to_string())
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert!(kd.has_shortcuts());
        assert!(kd.shortcut_names().is_some());
        assert!(kd.reaction_map().is_some());
        assert!(kd.library_pairs().is_none());
        assert!(!kd.has_normalized_view());
        assert_eq!(kd.state(), kd.workflow_state());
    }

    #[test]
    fn test_kindata_builder_facade_exposes_canonical_accessors() {
        let kd = KinData::builder()
            .with_direct_reactions(vec!["A -> B".to_string(), "B -> C".to_string()])
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(
            kd.equations(),
            &["A -> B".to_string(), "B -> C".to_string()]
        );
        assert_eq!(
            kd.canonical_reaction_ids().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0
            }
        );
        assert!(kd.shortcut_names().is_none());
        assert!(kd.library_pairs().is_none());
        assert!(kd.reaction_values().is_none());
        assert!(kd.reaction_data().is_none());
        assert!(kd.symbolic_constants().is_none());
        assert!(kd.normalized_reactions().is_none());
        assert!(!kd.has_normalized_view());
    }

    #[test]
    fn test_kindata_from_direct_reactions_and_queries() {
        let kd = KinData::from_direct_reactions(vec!["X -> Y".to_string()]).unwrap();

        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_eq!(kd.reaction_count(), 1);
        assert!(!kd.has_shortcuts());
        assert!(!kd.has_library_pairs());
        assert!(!kd.is_normalized());
        assert_eq!(kd.state(), kd.workflow_state());
        assert!(!kd.has_reaction_data());
        assert!(!kd.has_reaction_values());
        assert!(!kd.has_symbolic_constants());
    }

    #[test]
    fn test_kindata_set_reaction_data_directly_discovers_substances() {
        let mut kd = KinData::new();
        kd.set_reaction_data_directly(
            vec![
                ReactionData {
                    reaction_type: ReactionType::Elem,
                    eq: "A=>B".to_string(),
                    react: None,
                    data: ReactionKinetics::Elementary(ElementaryStruct {
                        Arrenius: vec![1.0e10, 0.0, 50_000.0],
                    }),
                },
                ReactionData {
                    reaction_type: ReactionType::Elem,
                    eq: "B=>A+C".to_string(),
                    react: None,
                    data: ReactionKinetics::Elementary(ElementaryStruct {
                        Arrenius: vec![1.0e8, 0.5, 30_000.0],
                    }),
                },
            ],
            None,
        )
        .unwrap();

        // Typed reaction input should leave KinData coherent before full molar-mass analysis.
        assert_eq!(kd.equations(), &["A=>B".to_string(), "B=>A+C".to_string()]);
        assert_eq!(
            kd.substances(),
            &["A".to_string(), "B".to_string(), "C".to_string()]
        );
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ReactionDataParsed
        );
        assert!(kd.has_reaction_data());
        assert!(kd.normalized_reactions().is_some());
    }

    #[test]
    fn test_kindata_workflow_state_moves_with_the_data_branch() {
        let mut kd = KinData::from_direct_reactions(vec!["A -> B".to_string()]).unwrap();

        // The facade should expose the branch that is currently loaded, not a manually set state.
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );

        kd.kinetic_main().unwrap();

        // After analysis the state should follow the data flow automatically.
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
        assert!(kd.is_analyzed());
        assert!(!kd.is_parsed());
        assert!(!kd.is_sorted());

        kd.refresh_state_from_available_data().unwrap();
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
    }

    #[test]
    fn test_kindata_mutation_paths_keep_state_invariants() {
        let mut kd = KinData::from_direct_reactions(vec!["A -> B".to_string()]).unwrap();
        assert_raw_lengths_match(&kd);

        kd.append_reaction(vec![json!({"eq": "C -> D"})]).unwrap();
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert_raw_lengths_match(&kd);
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "raw".to_string(),
                index: 0
            }
        );

        let mut file = NamedTempFile::new().unwrap();
        serde_json::to_writer_pretty(
            file.as_file_mut(),
            &json!({"R1": {"data": "value1"}, "R2": {"data": "value2"}}),
        )
        .unwrap();
        kd.load_reactions_from_json(file.path().to_str().unwrap())
            .unwrap();
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_raw_lengths_match(&kd);
        assert!(kd.has_shortcuts());
        assert!(kd.reaction_values().is_some());

        let mut normalized = KinData::new();
        normalized.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        normalized.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        normalized.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 1,
            },
        ]);
        normalized.shortcut_reactions = Some(vec!["R1".to_string(), "R2".to_string()]);
        normalized.vec_of_pairs = Some(vec![
            ("Lib1".to_string(), "Reaction1".to_string()),
            ("Lib2".to_string(), "Reaction2".to_string()),
        ]);
        normalized.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);

        normalized.remove_by_index(0).unwrap();
        assert_eq!(
            normalized.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::Sorted
        );
        assert_normalized_lengths(&normalized);
    }

    #[test]
    fn test_kindata_from_shortcuts_and_library_helpers() {
        let shortcuts = KinData::from_shortcuts(vec!["R1".to_string(), "R2".to_string()]).unwrap();
        assert_eq!(
            shortcuts.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert!(shortcuts.has_shortcuts());
        assert_eq!(shortcuts.reaction_count(), 2);

        let library = KinData::from_library_reactions(
            "Cantera".to_string(),
            vec!["1".to_string(), "2".to_string()],
        )
        .unwrap();
        assert_eq!(
            library.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );
        assert!(library.has_library_pairs());
        assert_eq!(library.reaction_count(), 2);
    }

    #[test]
    fn test_kindata_pure_output_helpers_return_data() {
        let mut normalized = KinData::new();
        normalized.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        normalized.vec_of_equations = vec!["A -> B".to_string()];
        normalized.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
        ]);
        normalized.vec_of_reaction_Values = Some(vec![json!({"eq": "A -> B"})]);
        normalized.shortcut_reactions = None;
        normalized.vec_of_pairs = None;

        let reaction_table = normalized.reaction_data_table().unwrap();
        let reaction_table_text = format!("{}", reaction_table);
        assert!(reaction_table_text.contains("Reaction number"));

        let normalized_table_text =
            format!("{}", normalized.normalized_reaction_data_table().unwrap());
        assert!(normalized_table_text.contains("Canonical ID"));

        let mut kd = KinData::from_direct_reactions(vec!["A -> B".to_string()]).unwrap();
        kd.kinetic_main().unwrap();
        let (substances_table, elements_table) = kd.substances_verbose_tables().unwrap();
        let substances_table_text = format!("{}", substances_table);
        let elements_table_text = format!("{}", elements_table);
        assert!(substances_table_text.contains("Reactions/Substances"));
        assert!(elements_table_text.contains("Substances/Elements"));
    }

    #[test]
    fn test_kindata_normalized_print_helpers_do_not_mutate_state() {
        let mut normalized = KinData::from_direct_reactions(vec!["A -> B".to_string()]).unwrap();
        normalized.kinetic_main().unwrap();

        // Materialize the normalized export view first; the pure print helpers
        // are expected to read this cache without changing it.
        normalized.every_reaction = Some(vec![make_every_reaction("direct", "0", "A -> B", 1.0)]);

        let before_state = normalized.state.clone();
        normalized
            .pretty_print_reaction_data_from_normalized()
            .unwrap();
        normalized.pretty_print_kindata_from_normalized().unwrap();

        assert_eq!(normalized.state, before_state);
    }

    #[test]
    fn test_kindata_pretty_print_reaction_data_does_not_mutate_live_state() {
        let mut normalized = KinData::new();
        normalized.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        normalized.vec_of_equations = vec!["A -> B".to_string()];
        normalized.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
        ]);
        normalized.every_reaction = Some(vec![make_every_reaction("direct", "0", "A -> B", 1.0)]);
        normalized.state = crate::Kinetics::User_reactions::KinDataState::Sorted;

        let before_state = normalized.state;
        let before_equations = normalized.vec_of_equations.clone();

        normalized.pretty_print_reaction_data().unwrap();

        assert_eq!(normalized.state, before_state);
        assert_eq!(normalized.vec_of_equations, before_equations);
    }

    #[test]
    fn test_kindata_kinetics_document_map_builds_without_writing() {
        let mut kd = KinData::builder()
            .with_direct_reactions(vec!["A -> B".to_string()])
            .unwrap()
            .build()
            .unwrap();

        kd.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
        ]);
        kd.vec_of_equations = vec!["A -> B".to_string()];
        kd.every_reaction = Some(vec![make_every_reaction("direct", "0", "A -> B", 1.0)]);

        let payload = kd.kinetics_document_map().unwrap();
        assert!(payload.contains_key("direct"));
        assert!(payload["direct"].contains_key("0"));

        // Repeating a read-only export should not perturb the already normalized state.
        let before_state = kd.state;
        let before_every_reaction = kd.every_reaction.clone();
        let before_every_reaction_len = before_every_reaction.as_ref().map(|values| values.len());
        let second_payload = kd.kinetics_document_map().unwrap();

        assert_eq!(payload, second_payload);
        assert_eq!(kd.state, before_state);
        assert_eq!(
            kd.every_reaction.as_ref().map(|values| values.len()),
            before_every_reaction_len
        );
        assert_eq!(
            kd.every_reaction.as_ref().map(|values| {
                values
                    .iter()
                    .map(|entry| entry.reaction_id.clone())
                    .collect::<Vec<_>>()
            }),
            before_every_reaction.as_ref().map(|values| {
                values
                    .iter()
                    .map(|entry| entry.reaction_id.clone())
                    .collect::<Vec<_>>()
            })
        );
    }

    #[test]
    fn test_materialized_every_reaction_exports_without_assignment_cache() {
        let mut kd = KinData::new();
        kd.every_reaction = Some(vec![make_every_reaction("LibA", "R1", "A -> B", 1.0)]);

        let payload = kd.kinetics_document_map().unwrap();

        assert_eq!(payload["LibA"]["R1"]["eq"], "A -> B");
        assert!(kd.reaction_ids.is_none());
        assert!(kd.vec_of_reaction_data.is_none());
    }

    #[test]
    fn test_refresh_recovers_assignment_ids_from_materialized_result() {
        let mut kd = KinData::new();
        kd.every_reaction = Some(vec![make_every_reaction("LibA", "R1", "A -> B", 1.0)]);

        kd.refresh_state_from_available_data().unwrap();

        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::Sorted
        );
        assert_eq!(
            kd.canonical_reaction_ids().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibA".to_string(),
                reaction_id: "R1".to_string(),
            }
        );
    }

    #[test]
    fn test_shortcut_assignment_survives_failed_database_resolution() {
        let mut kd = KinData::from_shortcuts(vec!["NoSuchLib_1".to_string()]).unwrap();
        let before_ids = kd.canonical_reaction_ids().unwrap().to_vec();

        let err = kd.get_reactions_from_shortcuts().unwrap_err();

        assert!(matches!(
            err,
            crate::Kinetics::error::KineticsError::MissingLibrary(library)
                if library == "nosuchlib"
        ));
        assert_eq!(kd.canonical_reaction_ids().unwrap(), before_ids.as_slice());
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert!(kd.reaction_values().is_none());
        assert!(kd.normalized_reactions().is_none());
    }

    #[test]
    fn test_direct_unknown_formula_analysis_is_transactional() {
        let mut kd = KinData::from_direct_reactions(vec!["A@ -> B".to_string()]).unwrap();

        let err = kd.kinetic_main().unwrap_err();

        assert!(matches!(
            err,
            crate::Kinetics::error::KineticsError::InvalidReactionData(_)
        ));
        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded
        );
        assert!(kd.substances().is_empty());
        assert!(kd.stoichiometric_analyzer().reactions.is_empty());
        assert!(kd.normalized_reactions().is_none());
    }

    #[test]
    fn test_kindata_kinetics_document_map_groups_multiple_library_reactions() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Cantera".to_string(),
                reaction_id: "1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Cantera".to_string(),
                reaction_id: "2".to_string(),
            },
        ]);
        kd.vec_of_reaction_Values = Some(vec![
            json!({"eq": "H2 + O2 -> H2O2"}),
            json!({"eq": "H2O2 -> H2O + O"}),
        ]);
        kd.every_reaction = Some(vec![
            make_every_reaction("Cantera", "1", "H2 + O2 -> H2O2", 1.0),
            make_every_reaction("Cantera", "2", "H2O2 -> H2O + O", 2.0),
        ]);

        let payload = kd.kinetics_document_map_from_normalized().unwrap();
        assert_eq!(payload.len(), 1);
        assert!(payload["Cantera"].contains_key("1"));
        assert!(payload["Cantera"].contains_key("2"));
    }

    #[test]
    fn test_kindata_kinetics_document_map_from_normalized_rejects_empty_view() {
        let mut kd = KinData::new();
        kd.every_reaction = Some(Vec::new());

        let err = kd.kinetics_document_map_from_normalized().unwrap_err();
        assert!(err.to_string().contains("normalized reaction view"));
    }

    #[test]
    fn test_kindata_kinetics_document_map_ignores_stale_raw_values_when_normalized_view_exists() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
        ]);
        kd.vec_of_equations = vec!["A -> B".to_string()];
        kd.vec_of_reaction_Values = Some(vec![
            json!({"eq": "A -> B"}),
            json!({"eq": "stale extra value"}),
        ]);
        kd.every_reaction = Some(vec![make_every_reaction("direct", "0", "A -> B", 1.0)]);

        let payload = kd.kinetics_document_map_from_normalized().unwrap();
        assert_eq!(payload.len(), 1);
        assert!(payload["direct"].contains_key("0"));
        assert_eq!(payload["direct"]["0"]["eq"], "A -> B");
    }

    #[test]
    fn test_kindata_every_reaction_rebuild_uses_reaction_data_equations_and_canonical_ids() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibA".to_string(),
                reaction_id: "R1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("R2".to_string()),
        ]);
        kd.vec_of_pairs = Some(vec![
            ("LegacyLib".to_string(), "Old1".to_string()),
            ("LegacyLib".to_string(), "Old2".to_string()),
        ]);
        kd.shortcut_reactions = Some(vec!["LEGACY_1".to_string(), "LEGACY_2".to_string()]);
        kd.vec_of_equations = vec!["stale eq 1".to_string(), "stale eq 2".to_string()];
        kd.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kd.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);

        let payload = kd.kinetics_document_map().unwrap();
        assert_eq!(payload["LibA"]["R1"]["eq"], "A -> B");
        assert_eq!(payload["SHORTCUTS"]["R2"]["eq"], "B -> C");

        let every_reaction = kd.every_reaction.as_ref().unwrap();
        assert_eq!(every_reaction[0].equation, "A -> B");
        assert_eq!(
            every_reaction[0].lib_and_id,
            Some(("LibA".to_string(), "R1".to_string()))
        );
        assert!(every_reaction[0].shortcut.is_none());
        assert_eq!(every_reaction[1].equation, "B -> C");
        assert_eq!(every_reaction[1].shortcut, Some("R2".to_string()));
        assert!(every_reaction[1].lib_and_id.is_none());
        assert_eq!(every_reaction[0].K_sym, Some(Expr::Const(1.0)));
        assert_eq!(every_reaction[1].K_sym, Some(Expr::Const(2.0)));
    }

    #[test]
    fn test_kindata_kinetics_document_map_rejects_stale_compatibility_cache_lengths() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 0,
            },
            crate::Kinetics::User_reactions::ReactionIdentity::Document {
                source: "direct".to_string(),
                index: 1,
            },
        ]);
        kd.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1.0, 0.0, 0.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![2.0, 0.0, 0.0],
                }),
            },
        ]);
        kd.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kd.shortcut_reactions = Some(vec!["stale_shortcut".to_string()]);
        kd.vec_of_pairs = Some(vec![("Legacy".to_string(), "1".to_string())]);

        let err = kd.kinetics_document_map().unwrap_err();
        assert!(
            err.to_string()
                .contains("shortcut_reactions vs vec_of_reaction_data")
        );
    }

    #[test]
    fn test_kindata_sort_then_remove_keeps_canonical_views_aligned() {
        let mut kd = KinData::new();
        kd.vec_of_reaction_data = Some(vec![
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "A -> B".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e12, 0.0, 20000.0],
                }),
            },
            ReactionData {
                reaction_type: ReactionType::Elem,
                eq: "B -> C".to_string(),
                react: None,
                data: ReactionKinetics::Elementary(ElementaryStruct {
                    Arrenius: vec![1e13, 0.0, 20000.0],
                }),
            },
        ]);
        kd.vec_of_equations = vec!["A -> B".to_string(), "B -> C".to_string()];
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibA".to_string(),
                reaction_id: "R1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibB".to_string(),
                reaction_id: "R2".to_string(),
            },
        ]);
        kd.shortcut_reactions = Some(vec!["LibA_R1".to_string(), "LibB_R2".to_string()]);
        kd.vec_of_pairs = Some(vec![
            ("LibA".to_string(), "R1".to_string()),
            ("LibB".to_string(), "R2".to_string()),
        ]);
        kd.K_sym_vec = Some(vec![Expr::Const(1.0), Expr::Const(2.0)]);

        let k_values = kd
            .calc_K_const_for_all_reactions(1000.0, None, None, Some(true))
            .unwrap();
        assert_eq!(k_values.len(), 2);
        assert_eq!(kd.every_reaction.as_ref().unwrap()[0].equation, "B -> C");

        kd.remove_by_index(1).unwrap();

        assert_eq!(
            kd.reaction_ids.as_ref().unwrap(),
            &vec![
                crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                    library: "LibB".to_string(),
                    reaction_id: "R2".to_string(),
                }
            ]
        );
        assert_eq!(kd.vec_of_equations, vec!["B -> C".to_string()]);
        assert_eq!(kd.shortcut_reactions, None);
        assert_eq!(
            kd.vec_of_pairs.as_ref().unwrap(),
            &vec![("LibB".to_string(), "R2".to_string())]
        );
        assert_eq!(kd.K_sym_vec.as_ref().unwrap(), &vec![Expr::Const(2.0)]);
        assert_eq!(kd.every_reaction.as_ref().unwrap().len(), 1);
        assert_eq!(kd.every_reaction.as_ref().unwrap()[0].equation, "B -> C");
        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::Sorted
        );
    }

    #[test]
    fn test_kindata_builder_aliases_match_legacy_names() {
        let kd = KinData::builder()
            .shortcut_reactions(vec!["C_1".to_string(), "C_2".to_string()])
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.workflow_state(),
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );
        assert_eq!(
            kd.shortcut_names(),
            Some(vec!["C_1".to_string(), "C_2".to_string()])
        );
    }

    #[test]
    fn test_kindata_builder_supports_library_pairs() {
        let kd = KinData::builder()
            .library_reactions(
                "Cantera".to_string(),
                vec!["1".to_string(), "2".to_string()],
            )
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );
        assert_eq!(
            kd.vec_of_pairs,
            Some(vec![
                ("Cantera".to_string(), "1".to_string()),
                ("Cantera".to_string(), "2".to_string())
            ])
        );
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[1],
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Cantera".to_string(),
                reaction_id: "2".to_string()
            }
        );
    }

    #[test]
    fn test_kindata_shortcut_views_ignore_stale_legacy_cache_when_canonical_ids_exist() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("C_1".to_string()),
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("C_2".to_string()),
        ]);
        kd.shortcut_reactions = Some(vec!["legacy_shortcut".to_string()]);
        kd.vec_of_pairs = Some(vec![("LegacyLib".to_string(), "99".to_string())]);

        assert_eq!(
            kd.shortcut_names(),
            Some(vec!["C_1".to_string(), "C_2".to_string()])
        );
        assert_eq!(
            kd.reaction_map(),
            Some(HashMap::from([(
                "Cantera".to_string(),
                vec!["1".to_string(), "2".to_string()]
            )]))
        );
        assert!(kd.library_pairs().is_none());
    }

    #[test]
    fn test_kindata_library_pairs_ignore_stale_legacy_cache_when_canonical_ids_exist() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryA".to_string(),
                reaction_id: "R1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryB".to_string(),
                reaction_id: "R2".to_string(),
            },
        ]);
        kd.shortcut_reactions = Some(vec!["legacy_shortcut".to_string()]);
        kd.vec_of_pairs = Some(vec![("LegacyLib".to_string(), "99".to_string())]);

        assert_eq!(
            kd.library_pairs(),
            Some(vec![
                ("LibraryA".to_string(), "R1".to_string()),
                ("LibraryB".to_string(), "R2".to_string()),
            ])
        );
        let reaction_map = kd.reaction_map().unwrap();
        assert_eq!(reaction_map.get("LibraryA"), Some(&vec!["R1".to_string()]));
        assert_eq!(reaction_map.get("LibraryB"), Some(&vec!["R2".to_string()]));
        assert!(!reaction_map.contains_key("LegacyLib"));
    }

    #[test]
    fn test_get_reactions_from_shortcuts_prefers_canonical_ids_over_stale_cache() {
        let mut kd = KinData::new();
        kd.set_reactions_from_shortcut_range("C1..C2".to_string())
            .unwrap();

        // Leave behind a stale legacy cache to make sure resolution follows the canonical ids.
        kd.shortcut_reactions = Some(vec!["legacy_shortcut".to_string()]);

        kd.get_reactions_from_shortcuts().unwrap();

        assert_eq!(
            kd.vec_of_pairs,
            Some(vec![
                ("Cantera".to_string(), "1".to_string()),
                ("Cantera".to_string(), "2".to_string())
            ])
        );
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[0],
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Cantera".to_string(),
                reaction_id: "1".to_string()
            }
        );
        assert_eq!(
            kd.reaction_ids.as_ref().unwrap()[1],
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Cantera".to_string(),
                reaction_id: "2".to_string()
            }
        );
    }

    #[test]
    fn test_get_reactions_from_shortcuts_rejects_non_shortcut_canonical_ids() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryA".to_string(),
                reaction_id: "R1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryB".to_string(),
                reaction_id: "R2".to_string(),
            },
        ]);
        kd.shortcut_reactions = Some(vec!["legacy_shortcut".to_string()]);

        let err = kd.get_reactions_from_shortcuts().unwrap_err();
        assert!(format!("{err:?}").contains("canonical reaction ids are not shortcut ids"));
    }

    #[test]
    fn test_kindata_public_views_require_canonical_ids() {
        let mut kd = KinData::new();
        kd.shortcut_reactions = Some(vec!["legacy_shortcut".to_string()]);
        kd.vec_of_pairs = Some(vec![("LegacyLib".to_string(), "99".to_string())]);

        assert!(kd.shortcut_names().is_none());
        assert!(kd.library_pairs().is_none());
        assert!(kd.reaction_map().is_none());
    }

    #[test]
    fn test_kindata_public_views_reject_mixed_canonical_ids() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::Shortcut("C_1".to_string()),
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryA".to_string(),
                reaction_id: "R1".to_string(),
            },
        ]);

        assert!(kd.shortcut_names().is_none());
        assert!(kd.library_pairs().is_none());
        assert!(kd.reaction_map().is_none());
    }

    #[test]
    fn test_kindata_state_tracks_workflow_progress() {
        let mut kd = KinData::new();

        kd.set_reactions_from_shortcut_range("C1..C1".to_string())
            .unwrap();
        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::ShortcutsSelected
        );

        let _ = kd.get_reactions_from_shortcuts();
        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::LibraryResolved
        );

        let _ = kd.kinetic_main();
        assert_eq!(
            kd.state,
            crate::Kinetics::User_reactions::KinDataState::Analyzed
        );
        assert!(kd.every_reaction.is_some());
        assert_eq!(
            kd.every_reaction.as_ref().unwrap()[0].reaction_id,
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "Cantera".to_string(),
                reaction_id: "1".to_string()
            }
        );
    }

    #[test]
    fn test_kindata_state_machine_rejects_analysis_from_empty_state() {
        let mut kd = KinData::new();
        kd.vec_of_equations = vec!["A -> B".to_string()];

        let result = kd.analyze_reactions();
        assert!(matches!(
            result,
            Err(crate::Kinetics::error::KineticsError::InvalidState(message))
            if message.contains("reaction analysis")
        ));
    }

    #[test]
    fn test_kindata_state_machine_rejects_symbolic_constants_before_analysis() {
        let mut kd = KinData::new();
        kd.state = crate::Kinetics::User_reactions::KinDataState::DirectReactionsLoaded;
        kd.vec_of_reaction_data = Some(vec![ReactionData {
            reaction_type: ReactionType::Elem,
            eq: "A -> B".to_string(),
            react: None,
            data: ReactionKinetics::Elementary(ElementaryStruct {
                Arrenius: vec![1.0, 0.0, 0.0],
            }),
        }]);

        let result = kd.calc_sym_constants(None, None, None);
        assert!(matches!(
            result,
            Err(crate::Kinetics::error::KineticsError::InvalidState(message))
            if message.contains("symbolic constant calculation")
        ));
    }

    #[test]
    fn test_kindata_calc_k_const_for_trange_empty_temperature_grid() {
        let mut kd = KinData::new();
        kd.set_reactions_from_shortcut_range("C1..C1".to_string())
            .unwrap();
        let _ = kd.get_reactions_from_shortcuts();
        let _ = kd.kinetic_main();

        let result =
            kd.calc_K_const_for_all_reactions_forTrange(300.0, 400.0, 0, None, None, Some(false));

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "reaction data validation failed: n must be at least 2 to include both boundary values"
        );
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
            let result = kd
                .set_reactions_from_shortcut_range(input.to_string())
                .unwrap();
            assert_eq!(result, expected, "Failed for input: {}", input);

            // Verify shortcut_reactions is set correctly
            let expected_shortcuts: Vec<String> = expected.iter().map(|s| s.to_string()).collect();
            assert_eq!(kd.shortcut_reactions, Some(expected_shortcuts));
        }
    }

    #[test]
    fn test_generate_strings_with_invalid_range() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("C1..X".to_string());
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "invalid shortcut range `C1..X`"
        );
    }

    #[test]
    fn test_generate_strings_rejects_missing_range_separator() {
        let mut user_reactions = KinData::new();
        let result = user_reactions.set_reactions_from_shortcut_range("C1-C5".to_string());

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "invalid shortcut range `C1-C5`"
        );
    }

    #[test]
    fn test_get_reactions_from_shortcuts_rejects_unknown_library_prefix() {
        let mut kd = KinData::new();
        kd.shortcut_reactions = Some(vec!["invalid".to_string()]);

        let err = kd.get_reactions_from_shortcuts().unwrap_err();
        assert!(matches!(
            err,
            crate::Kinetics::error::KineticsError::MissingLibrary(library)
                if library == "Unknown"
        ));
    }

    #[test]
    fn test_get_reactions_from_shortcuts_rejects_missing_reaction_id() {
        let mut kd = KinData::new();
        kd.shortcut_reactions = Some(vec!["C999999".to_string()]);

        let err = kd.get_reactions_from_shortcuts().unwrap_err();
        assert!(matches!(
            err,
            crate::Kinetics::error::KineticsError::MissingReaction(reaction_id)
                if reaction_id == "999999"
        ));
    }

    #[test]
    fn test_kindata_library_pairs_are_derived_from_canonical_ids() {
        let mut kd = KinData::new();
        kd.reaction_ids = Some(vec![
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryA".to_string(),
                reaction_id: "R1".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryA".to_string(),
                reaction_id: "R2".to_string(),
            },
            crate::Kinetics::User_reactions::ReactionIdentity::LibraryReaction {
                library: "LibraryB".to_string(),
                reaction_id: "R3".to_string(),
            },
        ]);
        kd.shortcut_reactions = None;
        kd.vec_of_pairs = None;

        assert_eq!(
            kd.library_pairs(),
            Some(vec![
                ("LibraryA".to_string(), "R1".to_string()),
                ("LibraryA".to_string(), "R2".to_string()),
                ("LibraryB".to_string(), "R3".to_string()),
            ])
        );
    }
}
