#[cfg(test)]
mod tests {
    use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;
    use crate::Thermodynamics::User_PhaseOrSolution::{
        CustomSubstance, NestedPhaseCacheView, PhaseComposition, PhaseDataPreparation,
        PhaseDataView, PhaseEquilibriumAssembly, PhaseEvaluationRequest, PhaseLayoutAccess,
        PhaseModel, PhaseOrSolution, PhasePhysicalState, PhasePropertyEvaluator, PhaseSpec,
        PhaseSymbolicPropertyBuilder, PhaseSystem, ResolvedPhaseSystem, SubstancePhaseMapping,
        SubstanceSystemFactory, SubstanceSystemSpec, SubstancesContainer,
        ThermoEvaluationConditions,
    };
    use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
    use crate::Thermodynamics::phase_layout::PhaseId;

    use crate::Thermodynamics::User_substances::SubsData;
    use crate::Thermodynamics::User_substances_error::SubsDataError;
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;
    use nalgebra::DMatrix;
    use std::collections::{HashMap, HashSet};
    use std::sync::Arc;
    use tempfile::tempdir;
    #[test]
    fn test_substance_phase_mapping_basic() {
        let subs = vec!["A".to_string(), "B".to_string()];
        let phases = vec![Some("p1".to_string()), Some("p2".to_string())];
        let mapping = SubstancePhaseMapping::new(subs.clone(), phases.clone()).unwrap();
        assert_eq!(mapping.all_substances, subs);
        assert_eq!(mapping.substance_to_phase, phases);
        assert_eq!(
            mapping.get_phase_for_substance(0).unwrap(),
            &Some("p1".to_string())
        );
    }

    #[test]
    fn test_substance_phase_mapping_rejects_length_mismatch() {
        let subs = vec!["A".to_string()];
        let phases = vec![Some("p1".to_string()), Some("p2".to_string())];
        let err = SubstancePhaseMapping::new(subs, phases).unwrap_err();
        assert!(format!("{}", err).contains("substance-phase mapping"));
    }

    #[test]
    fn test_create_system_single_and_multi_phase() {
        // Single phase
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        let container = SubstancesContainer::SinglePhase(subs.clone());
        let res = SubstanceSystemFactory::create_system(
            container,
            None,
            vec!["NASA_gas".to_string(), "NASA_cond".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );
        assert!(res.is_ok());
        match res.unwrap() {
            crate::Thermodynamics::User_PhaseOrSolution::CustomSubstance::OnePhase(op) => {
                assert_eq!(op.subs_data_view().substances, subs);
                assert_eq!(op.phase_spec().unwrap().model(), PhaseModel::IdealGas);
                assert_eq!(
                    op.resolved_layout().unwrap().component_labels(),
                    vec!["CO".to_string(), "CO2".to_string()]
                );
            }
            _ => panic!("expected one phase"),
        }

        // Multi phase
        let mut phases = HashMap::new();
        phases.insert("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]);
        phases.insert("liquid".to_string(), vec!["H2O".to_string()]);
        let container = SubstancesContainer::MultiPhase(phases);
        let res = SubstanceSystemFactory::create_system(
            container,
            Some(HashMap::from([
                (
                    "gas".to_string(),
                    crate::Thermodynamics::User_substances::Phases::Gas,
                ),
                (
                    "liquid".to_string(),
                    crate::Thermodynamics::User_substances::Phases::Liquid,
                ),
            ])),
            vec!["NASA_gas".to_string(), "NASA_cond".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );
        assert!(res.is_ok());
        match res.unwrap() {
            crate::Thermodynamics::User_PhaseOrSolution::CustomSubstance::PhaseOrSolution(p) => {
                assert!(p.phase_data_view().contains_key(&Some("gas".to_string())));
                assert!(
                    p.phase_data_view()
                        .contains_key(&Some("liquid".to_string()))
                );
                assert_eq!(p.phase_specs().len(), 2);
                assert_eq!(p.resolved_layout().unwrap().component_count(), 3);
            }
            _ => panic!("expected phase or solution"),
        }
    }

    #[test]
    fn test_substance_system_spec_builder_resolves_single_phase() {
        let spec = crate::Thermodynamics::User_PhaseOrSolution::SubstanceSystemSpec::builder(
            SubstancesContainer::SinglePhase(vec!["CO".to_string(), "CO2".to_string()]),
        )
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();

        let resolved = spec.resolve().unwrap();
        match resolved {
            CustomSubstance::OnePhase(one) => {
                assert_eq!(
                    one.subs_data_view().substances,
                    vec!["CO".to_string(), "CO2".to_string()]
                );
                let phase = one
                    .phase_spec()
                    .expect("factory must retain the resolved phase spec");
                assert_eq!(phase.id().as_option(), &None);
                assert_eq!(phase.model(), PhaseModel::IdealGas);
            }
            _ => panic!("expected one-phase system"),
        }
    }

    #[test]
    fn test_substance_system_spec_rejects_empty_single_phase() {
        let err = SubstanceSystemSpec::builder(SubstancesContainer::SinglePhase(Vec::new()))
            .build()
            .unwrap_err();

        assert!(format!("{}", err).contains("each phase must contain at least one substance"));
    }

    #[test]
    fn test_phase_system_new_single_phase_seeds_canonical_none_caches() {
        let system = PhaseSystem::new_single_phase();

        assert_eq!(system.layout_revision(), 0);
        assert!(system.dG_view().contains_key(&None));
        assert!(system.dG_fun_view().contains_key(&None));
        assert!(system.dG_sym_view().contains_key(&None));
        assert!(system.dS_view().contains_key(&None));
        assert!(system.dS_fun_view().contains_key(&None));
        assert!(system.dS_sym_view().contains_key(&None));
        assert!(
            system
                .dG_view()
                .get(&None)
                .expect("canonical single-phase Gibbs cache missing")
                .is_empty()
        );
        assert!(
            system
                .dS_view()
                .get(&None)
                .expect("canonical single-phase entropy cache missing")
                .is_empty()
        );
    }

    #[test]
    fn test_phase_system_replace_methods_bump_revision_and_keep_snapshot_typed() {
        let mut system = PhaseSystem::new_single_phase();
        let start_revision = system.layout_revision();

        system.replace_dG(HashMap::from([(
            None,
            HashMap::from([("A".to_string(), 1.25)]),
        )]));

        assert!(system.layout_revision() > start_revision);
        let snapshot = system.state_snapshot(false);
        assert_eq!(snapshot.layout_revision, system.layout_revision());
        match snapshot.dG.values {
            crate::Thermodynamics::User_PhaseOrSolution::NestedPhaseCacheView::Multi(view) => {
                assert_eq!(
                    view.get(&None).and_then(|inner| inner.get("A")).copied(),
                    Some(1.25)
                );
            }
            _ => panic!("expected multi-phase snapshot view"),
        }
    }

    #[test]
    fn test_phase_system_clone_preserves_cache_representations_without_aliasing_maps() {
        let mut original = PhaseSystem::new_single_phase();
        let numeric = HashMap::from([(None, HashMap::from([("A".to_string(), -12.5)]))]);
        let symbolic = HashMap::from([(
            None,
            HashMap::from([("A".to_string(), Expr::Var("T".to_string()))]),
        )]);
        let functions = HashMap::from([(
            None,
            HashMap::from([(
                "A".to_string(),
                Arc::new(|temperature, _, _| temperature + 2.0)
                    as crate::Thermodynamics::User_PhaseOrSolution::PhaseFunction,
            )]),
        )]);

        original.replace_dG(numeric.clone());
        original.replace_dS(numeric.clone());
        original.replace_dG_sym(symbolic.clone());
        original.replace_dS_sym(symbolic.clone());
        original.replace_dG_fun(functions.clone());
        original.replace_dS_fun(functions);

        let mut cloned = original.clone();
        assert_eq!(cloned.dG_view().to_owned_map(), numeric);
        assert_eq!(cloned.dS_view().to_owned_map(), numeric);
        assert_eq!(cloned.dG_sym_view().to_owned_map(), symbolic);
        assert_eq!(cloned.dS_sym_view().to_owned_map(), symbolic);
        assert_eq!(
            cloned
                .dG_fun_view()
                .get(&None)
                .and_then(|phase| phase.get("A"))
                .expect("cloned Gibbs function missing")(300.0, None, None),
            302.0
        );
        assert_eq!(
            cloned
                .dS_fun_view()
                .get(&None)
                .and_then(|phase| phase.get("A"))
                .expect("cloned entropy function missing")(300.0, None, None),
            302.0
        );

        // Replacing a clone's numeric cache must not mutate the original map.
        cloned.replace_dG(HashMap::from([(
            None,
            HashMap::from([("A".to_string(), 7.0)]),
        )]));
        assert_eq!(
            original
                .dG_view()
                .get(&None)
                .and_then(|phase| phase.get("A")),
            Some(&-12.5)
        );
        assert_eq!(
            cloned.dG_view().get(&None).and_then(|phase| phase.get("A")),
            Some(&7.0)
        );
    }

    #[test]
    fn test_phase_payload_updates_are_transactional_when_a_later_phase_fails() {
        let mut system = PhaseSystem::new();
        system.replace_phase_data(HashMap::from([
            (Some("gas".to_string()), SubsData::new()),
            (Some("liquid".to_string()), SubsData::new()),
        ]));
        let revision_before = system.layout_revision();

        let error = system
            .debug_try_apply_to_phase_data(|phase_name, subs_data| {
                subs_data.set_P(101_325.0, Some("Pa".to_string()))?;
                if phase_name.as_deref() == Some("liquid") {
                    return Err(SubsDataError::MissingData {
                        field: "test rollback trigger".to_string(),
                        substance: "liquid".to_string(),
                    });
                }
                Ok(())
            })
            .expect_err("one phase must reject the transactional update");

        assert!(error.to_string().contains("test rollback trigger"));
        assert_eq!(system.layout_revision(), revision_before);
        assert!(
            system
                .phase_data()
                .values()
                .all(|subs_data| subs_data.pressure().is_none())
        );
    }

    #[test]
    fn test_numeric_cache_context_requires_exact_request_and_current_configuration() {
        let mut system = PhaseSystem::new_single_phase();
        system.with_single_phase_data_mut(|data| data.substances = vec!["A".to_string()]);
        let conditions = ThermoEvaluationConditions::new(300.0, 101_325.0).unwrap();
        let compositions =
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![1.0])))]);
        let values = HashMap::from([(None, HashMap::from([("A".to_string(), -10.0)]))]);

        system
            .publish_gibbs_evaluation(conditions, compositions.clone(), values)
            .unwrap();
        assert!(system.has_current_gibbs_cache(conditions, &compositions));
        assert_eq!(
            system
                .gibbs_cache_context()
                .expect("published Gibbs cache must retain its request")
                .conditions(),
            conditions
        );

        let different_composition =
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![0.8])))]);
        assert!(!system.has_current_gibbs_cache(conditions, &different_composition));

        system.bump_layout_revision();
        assert!(!system.has_current_gibbs_cache(conditions, &compositions));
        assert!(system.gibbs_cache_context().is_none());
    }

    #[test]
    fn test_phase_evaluation_request_is_the_atomic_cache_identity() {
        let mut system = PhaseSystem::new_single_phase();
        system.with_single_phase_data_mut(|data| data.substances = vec!["A".to_string()]);
        let request = PhaseEvaluationRequest::from_numeric(
            300.0,
            101_325.0,
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![1.0])))]),
        )
        .unwrap();
        let values = HashMap::from([(None, HashMap::from([("A".to_string(), -10.0)]))]);

        system
            .publish_gibbs_evaluation_request(request.clone(), values)
            .unwrap();
        let context = system
            .gibbs_cache_context()
            .expect("published Gibbs cache must retain its typed request");
        assert_eq!(context.request(), &request);
        assert!(system.has_current_gibbs_request(&request));

        let different_request =
            PhaseEvaluationRequest::from_numeric(301.0, 101_325.0, request.compositions().clone())
                .unwrap();
        assert!(!system.has_current_gibbs_request(&different_request));
        assert!(
            PhaseEvaluationRequest::from_numeric(f64::NAN, 101_325.0, HashMap::new(),).is_err()
        );
    }

    #[test]
    fn test_invalid_numeric_cache_publication_keeps_previous_cache_and_context() {
        let mut system = PhaseSystem::new_single_phase();
        system.with_single_phase_data_mut(|data| data.substances = vec!["A".to_string()]);
        let request = PhaseEvaluationRequest::from_numeric(
            300.0,
            101_325.0,
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![1.0])))]),
        )
        .unwrap();
        system
            .publish_gibbs_evaluation_request(
                request.clone(),
                HashMap::from([(None, HashMap::from([("A".to_string(), -10.0)]))]),
            )
            .unwrap();
        let revision_before = system.layout_revision();

        let error = system
            .publish_gibbs_evaluation_request(
                request.clone(),
                HashMap::from([(
                    Some("unknown".to_string()),
                    HashMap::from([("A".to_string(), -20.0)]),
                )]),
            )
            .expect_err("an unknown cache phase must not be published");

        assert!(error.to_string().contains("cache phase"));
        assert_eq!(system.layout_revision(), revision_before);
        assert_eq!(
            system
                .dG_view()
                .get(&None)
                .and_then(|values| values.get("A")),
            Some(&-10.0)
        );
        assert!(system.has_current_gibbs_request(&request));
    }

    #[test]
    fn test_symbolic_gibbs_substitution_publishes_a_new_cache_revision() {
        let mut system = PhaseSystem::new_single_phase();
        system.with_single_phase_data_mut(|data| data.substances = vec!["A".to_string()]);
        system.replace_dG_sym(HashMap::from([(
            None,
            HashMap::from([(
                "A".to_string(),
                Expr::Add(
                    Box::new(Expr::Var("T".to_string())),
                    Box::new(Expr::Var("P".to_string())),
                ),
            )]),
        )]));
        let before_pressure = system.layout_revision();

        system.set_pressure_in_gibbs_sym(101_325.0);
        assert!(system.layout_revision() > before_pressure);
        let after_pressure = system.layout_revision();
        assert!(
            format!(
                "{}",
                system.dG_sym_view().get(&None).unwrap().get("A").unwrap()
            )
            .contains("101325")
        );

        system.set_temperature_in_gibbs_sym(300.0);
        assert!(system.layout_revision() > after_pressure);
        let symbolic_values = system.dG_sym_view();
        let expression = symbolic_values.get(&None).unwrap().get("A").unwrap();
        assert!(!format!("{expression}").contains("T"));
        assert!(!format!("{expression}").contains("P"));
    }

    #[test]
    fn test_phase_composition_try_new_numeric_validates_sum_and_component_values() {
        let phase = Some("gas".to_string());
        let ok = PhaseComposition::try_new_numeric(
            Some(1.0),
            Some(vec![0.4, 0.6]),
            &phase,
            "Gibbs free energy",
        );
        assert!(ok.is_ok());

        let bad_sum = PhaseComposition::try_new_numeric(
            Some(1.0),
            Some(vec![0.2, 0.3]),
            &phase,
            "Gibbs free energy",
        );
        assert!(
            bad_sum
                .unwrap_err()
                .to_string()
                .contains("does not match component sum")
        );

        let bad_component = PhaseComposition::try_new_numeric(
            Some(1.0),
            Some(vec![0.5, -0.1, 0.6]),
            &phase,
            "Gibbs free energy",
        );
        assert!(
            bad_component
                .unwrap_err()
                .to_string()
                .contains("invalid component amount")
        );

        for non_finite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let invalid_total = PhaseComposition::try_new_numeric(
                Some(non_finite),
                Some(vec![1.0]),
                &phase,
                "Gibbs free energy",
            );
            assert!(invalid_total.is_err());

            let invalid_component = PhaseComposition::try_new_numeric(
                Some(1.0),
                Some(vec![non_finite]),
                &phase,
                "Gibbs free energy",
            );
            assert!(invalid_component.is_err());
        }
    }

    #[test]
    fn public_typed_phase_boundaries_return_errors_without_unwinding() {
        let outcome = std::panic::catch_unwind(|| {
            assert!(
                SubstancePhaseMapping::new(
                    vec!["A".to_string()],
                    vec![Some("gas".to_string()), Some("liquid".to_string())],
                )
                .is_err()
            );
            assert!(ThermoEvaluationConditions::new(f64::NAN, 101_325.0).is_err());
            assert!(
                PhaseComposition::try_new_numeric(
                    Some(1.0),
                    Some(vec![-1.0]),
                    &None,
                    "public boundary test",
                )
                .is_err()
            );

            let missing_states = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(
                HashMap::from([("gas".to_string(), vec!["A".to_string()])]),
            ))
            .build();
            assert!(missing_states.is_err());

            let one_phase = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
            assert!(
                one_phase
                    .evaluate_gibbs(f64::NAN, 101_325.0, &HashMap::new())
                    .is_err()
            );

            let multi_phase = PhaseOrSolution::new();
            assert!(
                multi_phase
                    .evaluate_entropy(300.0, 0.0, &HashMap::new())
                    .is_err()
            );
        });

        assert!(outcome.is_ok(), "a public typed boundary panicked");
    }

    #[test]
    fn test_failed_evaluation_does_not_publish_or_overwrite_cached_results() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.debug_dG_mut().insert("A".to_string(), 7.0);
        let invalid_single =
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![0.8])))]);
        assert!(
            one.evaluate_gibbs(300.0, 101325.0, &invalid_single)
                .is_err()
        );
        assert_eq!(one.dG_view().get("A"), Some(&7.0));

        let mut multi = PhaseOrSolution::new();
        multi.debug_replace_dS(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("A".to_string(), 9.0)]),
        )]));
        let invalid_multi = HashMap::from([(
            Some("gas".to_string()),
            PhaseComposition::new(Some(1.0), Some(vec![0.2])),
        )]);
        assert!(
            multi
                .evaluate_entropy(300.0, 101325.0, &invalid_multi)
                .is_err()
        );
        assert_eq!(
            multi
                .dS_view()
                .get(&Some("gas".to_string()))
                .and_then(|values| values.get("A")),
            Some(&9.0)
        );
    }

    #[test]
    fn test_numeric_evaluation_rejects_noncanonical_phase_requests_before_cache_publish() {
        let mut system = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["A".to_string(), "B".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["A".to_string()];
        system.insert_phase_data("gas".to_string(), gas);
        system.insert_phase_data("liquid".to_string(), liquid);
        system.debug_replace_dS(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("A".to_string(), 11.0)]),
        )]));

        let missing_phase = HashMap::from([(
            Some("gas".to_string()),
            PhaseComposition::new(Some(1.0), Some(vec![0.4, 0.6])),
        )]);
        let missing_phase_error = system
            .evaluate_entropy(300.0, 101_325.0, &missing_phase)
            .unwrap_err();
        assert!(
            missing_phase_error
                .to_string()
                .contains("must exactly match the resolved phase set")
        );

        let extra_phase = HashMap::from([
            (
                Some("gas".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![0.4, 0.6])),
            ),
            (
                Some("liquid".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![1.0])),
            ),
            (
                Some("unexpected".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![1.0])),
            ),
        ]);
        let extra_phase_error = system
            .evaluate_entropy(300.0, 101_325.0, &extra_phase)
            .unwrap_err();
        assert!(
            extra_phase_error
                .to_string()
                .contains("must exactly match the resolved phase set")
        );

        let wrong_component_count = HashMap::from([
            (
                Some("gas".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![1.0])),
            ),
            (
                Some("liquid".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![1.0])),
            ),
        ]);
        let wrong_component_error = system
            .evaluate_entropy(300.0, 101_325.0, &wrong_component_count)
            .unwrap_err();
        assert!(
            wrong_component_error
                .to_string()
                .contains("component amount vector has length 1; expected 2")
        );

        let complete_composition = HashMap::from([
            (
                Some("gas".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![0.4, 0.6])),
            ),
            (
                Some("liquid".to_string()),
                PhaseComposition::new(Some(1.0), Some(vec![1.0])),
            ),
        ]);
        assert!(
            system
                .evaluate_entropy(f64::NAN, 101_325.0, &complete_composition)
                .is_err()
        );
        assert!(
            system
                .evaluate_entropy(300.0, 0.0, &complete_composition)
                .is_err()
        );
        assert_eq!(
            system
                .dS_view()
                .get(&Some("gas".to_string()))
                .and_then(|values| values.get("A")),
            Some(&11.0)
        );
    }

    #[test]
    fn test_single_phase_evaluation_requires_none_key_and_declared_component_count() {
        let mut system = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        system.set_substances(vec!["A".to_string(), "B".to_string()]);

        let wrong_phase = HashMap::from([(
            Some("gas".to_string()),
            PhaseComposition::new(Some(1.0), Some(vec![0.4, 0.6])),
        )]);
        assert!(
            system
                .evaluate_gibbs(300.0, 101_325.0, &wrong_phase)
                .unwrap_err()
                .to_string()
                .contains("must exactly match the resolved phase set")
        );

        let wrong_component_count =
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![1.0])))]);
        assert!(
            system
                .evaluate_gibbs(300.0, 101_325.0, &wrong_component_count)
                .unwrap_err()
                .to_string()
                .contains("component amount vector has length 1; expected 2")
        );
    }

    #[test]
    fn test_single_phase_symbolic_layout_view_returns_atomic_bundle() {
        let one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        let layout = one.symbolic_layout_view();
        assert!(layout.is_empty());
        assert!(one.current_layout().components().is_empty());
    }

    #[test]
    fn test_one_phase_uses_canonical_layout_for_symbolic_and_mole_snapshots() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.set_substances(vec!["B".to_string(), "A".to_string()]);

        let (indexed, amounts, phase_totals, by_substance) = one.indexed_moles_variables().unwrap();
        assert_eq!(
            amounts,
            vec![Expr::Var("n0_0".to_string()), Expr::Var("n0_1".to_string())]
        );
        assert_eq!(phase_totals, vec![Expr::Var("Np0".to_string())]);
        assert_eq!(
            indexed.get(&None),
            Some(&(
                Some(Expr::Var("Np0".to_string())),
                Some(vec![
                    Expr::Var("n0_0".to_string()),
                    Expr::Var("n0_1".to_string()),
                ]),
            ))
        );
        assert_eq!(
            by_substance.get(&None).and_then(|vars| vars.get("B")),
            Some(&Expr::Var("n0_0".to_string()))
        );
        let snapshot = one.symbolic_layout_view();
        assert_eq!(
            one.current_layout().component_labels(),
            vec!["B".to_string(), "A".to_string()]
        );
        assert_eq!(snapshot.component_variables(), amounts.as_slice());
        assert_eq!(snapshot.phase_totals(), phase_totals.as_slice());
        assert_eq!(snapshot.phase_variables(), &indexed);
        assert_eq!(snapshot.variables_by_phase_and_substance(), &by_substance);

        let input = HashMap::from([(
            None,
            (Some(2.0), Some(HashMap::from([("A".to_string(), 1.5)]))),
        )]);
        let (full, vectors, totals) = one.create_full_map_of_mole_numbers(input).unwrap();
        assert_eq!(full[&None].1.as_ref().unwrap()["B"], 0.0);
        assert_eq!(vectors[&None].1.as_ref().unwrap(), &vec![0.0, 1.5]);
        assert_eq!(totals["A"], 1.5);
        assert_eq!(totals["B"], 0.0);
    }

    #[test]
    fn test_ideal_gas_gibbs_numeric_closure_and_symbolic_paths_agree() {
        // This is a representation-contract test. The same resolved gas
        // composition must retain its thermodynamic meaning across the pure
        // numeric evaluator, cached closure family, and symbolic expression.
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let mut one = match spec.resolve().unwrap() {
            CustomSubstance::OnePhase(one) => one,
            _ => panic!("single-phase specification must resolve to OnePhase"),
        };

        let temperature = 700.0;
        let pressure = 150_000.0;
        let total_amount = 1.0;
        let component_amounts = vec![0.7, 0.3];
        let compositions = HashMap::from([(
            None,
            PhaseComposition::new(Some(total_amount), Some(component_amounts.clone())),
        )]);

        let numeric = one
            .evaluate_gibbs(temperature, pressure, &compositions)
            .unwrap();
        let numeric = numeric.get(&None).unwrap();

        one.indexed_moles_variables().unwrap();
        PhaseSymbolicPropertyBuilder::build_gibbs_functions(&mut one, temperature, pressure)
            .unwrap();
        PhaseSymbolicPropertyBuilder::build_symbolic_gibbs(&mut one, temperature).unwrap();
        PhaseSymbolicPropertyBuilder::substitute_pressure_in_symbolic_gibbs(&mut one, pressure);
        PhaseSymbolicPropertyBuilder::substitute_temperature_in_symbolic_gibbs(
            &mut one,
            temperature,
        );

        let closure_values = one.dG_fun_view();
        let symbolic_values = one.dG_sym_view();
        let symbolic_arguments = ["n0_0", "n0_1", "Np0"];
        let symbolic_inputs = [component_amounts[0], component_amounts[1], total_amount];

        for substance in ["N2", "O2"] {
            let numeric_value = numeric[substance];
            let closure_value = closure_values[substance](
                temperature,
                Some(component_amounts.clone()),
                Some(total_amount),
            );
            let symbolic_value = symbolic_values[substance]
                .clone()
                .lambdify_borrowed_thread_safe(&symbolic_arguments)(
                &symbolic_inputs
            );

            assert_relative_eq!(closure_value, numeric_value, epsilon = 1e-6);
            assert_relative_eq!(symbolic_value, numeric_value, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_phase_lifecycle_resolves_evaluates_invalidates_and_rebuilds() {
        // A cache is valid only for the exact resolved payload that created it.
        // This story makes the lifecycle visible instead of testing its stages
        // in unrelated fixtures.
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let mut one = match spec.clone().resolve().unwrap() {
            CustomSubstance::OnePhase(one) => one,
            _ => panic!("single-phase specification must resolve to OnePhase"),
        };
        let temperature = 650.0;
        let pressure = 101_325.0;
        let composition =
            HashMap::from([(None, PhaseComposition::new(Some(1.0), Some(vec![0.6, 0.4])))]);

        assert_eq!(
            one.resolved_layout().unwrap().component_labels(),
            vec!["N2".to_string(), "O2".to_string()]
        );
        assert!(
            one.evaluate_gibbs(temperature, pressure, &composition)
                .unwrap()
                .get(&None)
                .unwrap()
                .values()
                .all(|value| value.is_finite())
        );
        assert!(one.dG_view().is_empty());

        one.calcutate_Gibbs_free_energy(temperature, pressure, composition)
            .unwrap();
        assert!(!one.dG_view().is_empty());
        let revision_before_mutation = one.layout_revision();

        one.with_subs_data_mut(|data| data.substances.push("CO2".to_string()));
        assert!(one.dG_view().is_empty());
        assert!(one.resolved_layout().is_none());
        assert!(one.resolution_report().is_none());
        assert!(one.layout_revision() > revision_before_mutation);

        let rebuilt = match spec.resolve().unwrap() {
            CustomSubstance::OnePhase(one) => one,
            _ => panic!("rebuilding the same specification must remain one-phase"),
        };
        assert_eq!(
            rebuilt.resolved_layout().unwrap().component_labels(),
            vec!["N2".to_string(), "O2".to_string()]
        );
        assert!(rebuilt.resolution_report().is_some());
        assert!(rebuilt.dG_view().is_empty());
    }

    #[test]
    fn test_multiphase_lifecycle_resolves_publishes_invalidates_and_rebuilds() {
        // The multi-phase facade must obey the same lifecycle contract as the
        // one-phase facade: changing one payload invalidates every derived
        // system-wide artifact rather than retaining a partly valid cache.
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["N2".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            (
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            ),
            (
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let mut system = match spec.clone().resolve().unwrap() {
            CustomSubstance::PhaseOrSolution(system) => system,
            _ => panic!("two named phases must resolve to PhaseOrSolution"),
        };
        let temperature = 650.0;
        let pressure = 101_325.0;
        let composition = HashMap::from([
            (
                Some("gas".to_string()),
                PhaseComposition::new(Some(0.4), Some(vec![0.4])),
            ),
            (
                Some("liquid".to_string()),
                PhaseComposition::new(Some(0.6), Some(vec![0.6])),
            ),
        ]);

        assert!(
            system
                .evaluate_gibbs(temperature, pressure, &composition)
                .unwrap()
                .values()
                .flatten()
                .all(|(_, value)| value.is_finite())
        );
        assert!(system.dG_view().is_empty());

        system
            .calcutate_Gibbs_free_energy(temperature, pressure, composition)
            .unwrap();
        assert_eq!(system.dG_view().len(), 2);
        assert!(system.resolution_report().is_some());
        let revision_before_mutation = system.layout_revision();

        system.with_phase_data_mut(|phase_data| {
            phase_data
                .get_mut(&Some("gas".to_string()))
                .expect("resolved gas phase must exist")
                .substances
                .push("CO2".to_string());
        });
        assert!(system.dG_view().is_empty());
        assert!(system.resolved_layout().is_none());
        assert!(system.resolution_report().is_none());
        assert!(system.layout_revision() > revision_before_mutation);

        let rebuilt = match spec.resolve().unwrap() {
            CustomSubstance::PhaseOrSolution(system) => system,
            _ => panic!("rebuilding the same spec must remain multi-phase"),
        };
        assert_eq!(
            rebuilt.resolved_layout().unwrap().component_labels(),
            vec!["gas::N2".to_string(), "liquid::H2O".to_string()]
        );
        assert!(rebuilt.resolution_report().is_some());
        assert!(rebuilt.dG_view().is_empty());
    }

    #[test]
    fn test_pure_condensed_gibbs_numeric_closure_and_symbolic_paths_agree() {
        // A condensed phase may reuse a NASA source record, but its phase model
        // must suppress the ideal-gas concentration correction consistently in
        // every representation.
        let spec =
            SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([(
                "liquid".to_string(),
                vec!["H2O".to_string()],
            )])))
            .with_phase_natures(Some(HashMap::from([(
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            )])))
            .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
            .with_permitted_libraries(vec!["NIST".to_string()])
            .with_search_in_nist(false)
            .build()
            .unwrap();
        let mut system = match spec.resolve().unwrap() {
            CustomSubstance::PhaseOrSolution(system) => system,
            _ => panic!("named condensed phase must use the multi-phase facade"),
        };
        let temperature = 700.0;
        let pressure = 250_000.0;
        let phase_key = Some("liquid".to_string());
        let compositions = HashMap::from([(
            phase_key.clone(),
            PhaseComposition::new(Some(1.0), Some(vec![1.0])),
        )]);
        let numeric = system
            .evaluate_gibbs(temperature, pressure, &compositions)
            .unwrap();
        let numeric_value = numeric[&phase_key]["H2O"];

        system.indexed_moles_variables().unwrap();
        PhaseSymbolicPropertyBuilder::build_gibbs_functions(&mut system, temperature, pressure)
            .unwrap();
        PhaseSymbolicPropertyBuilder::build_symbolic_gibbs(&mut system, temperature).unwrap();
        PhaseSymbolicPropertyBuilder::substitute_pressure_in_symbolic_gibbs(&mut system, pressure);
        PhaseSymbolicPropertyBuilder::substitute_temperature_in_symbolic_gibbs(
            &mut system,
            temperature,
        );

        let closure_value = system.dG_fun_view().get(&phase_key).unwrap()["H2O"](
            temperature,
            Some(vec![1.0]),
            Some(1.0),
        );
        let symbolic_value = system.dG_sym_view().get(&phase_key).unwrap()["H2O"]
            .clone()
            .lambdify_borrowed_thread_safe(&[])(&[]);

        assert_relative_eq!(closure_value, numeric_value, epsilon = 1e-6);
        assert_relative_eq!(symbolic_value, numeric_value, epsilon = 1e-6);
    }

    #[test]
    fn test_ideal_gas_entropy_numeric_closure_and_symbolic_paths_agree() {
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let mut one = match spec.resolve().unwrap() {
            CustomSubstance::OnePhase(one) => one,
            _ => panic!("single-phase specification must resolve to OnePhase"),
        };
        let temperature = 700.0;
        let pressure = 150_000.0;
        let total_amount = 1.0;
        let component_amounts = vec![0.7, 0.3];
        let compositions = HashMap::from([(
            None,
            PhaseComposition::new(Some(total_amount), Some(component_amounts.clone())),
        )]);
        let numeric = one
            .evaluate_entropy(temperature, pressure, &compositions)
            .unwrap();
        let numeric = numeric.get(&None).unwrap();

        one.indexed_moles_variables().unwrap();
        PhaseSymbolicPropertyBuilder::build_entropy_functions(&mut one, temperature, pressure)
            .unwrap();
        PhaseSymbolicPropertyBuilder::build_symbolic_entropy(&mut one, temperature).unwrap();

        let closure_values = one.dS_fun_view();
        let symbolic_values = one.dS_sym_view();
        let symbolic_arguments = ["T", "P", "n0_0", "n0_1", "Np0"];
        let symbolic_inputs = [
            temperature,
            pressure,
            component_amounts[0],
            component_amounts[1],
            total_amount,
        ];

        for substance in ["N2", "O2"] {
            let numeric_value = numeric[substance];
            let closure_value = closure_values[substance](
                temperature,
                Some(component_amounts.clone()),
                Some(total_amount),
            );
            let symbolic_value = symbolic_values[substance]
                .clone()
                .lambdify_borrowed_thread_safe(&symbolic_arguments)(
                &symbolic_inputs
            );

            assert_relative_eq!(closure_value, numeric_value, epsilon = 1e-6);
            assert_relative_eq!(symbolic_value, numeric_value, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_pure_condensed_entropy_numeric_closure_and_symbolic_paths_agree() {
        let spec =
            SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([(
                "liquid".to_string(),
                vec!["H2O".to_string()],
            )])))
            .with_phase_natures(Some(HashMap::from([(
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            )])))
            .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
            .with_permitted_libraries(vec!["NIST".to_string()])
            .with_search_in_nist(false)
            .build()
            .unwrap();
        let mut system = match spec.resolve().unwrap() {
            CustomSubstance::PhaseOrSolution(system) => system,
            _ => panic!("named condensed phase must use the multi-phase facade"),
        };
        let temperature = 700.0;
        let pressure = 250_000.0;
        let phase_key = Some("liquid".to_string());
        let compositions = HashMap::from([(
            phase_key.clone(),
            PhaseComposition::new(Some(1.0), Some(vec![1.0])),
        )]);
        let numeric_value = system
            .evaluate_entropy(temperature, pressure, &compositions)
            .unwrap()[&phase_key]["H2O"];

        system.indexed_moles_variables().unwrap();
        PhaseSymbolicPropertyBuilder::build_entropy_functions(&mut system, temperature, pressure)
            .unwrap();
        PhaseSymbolicPropertyBuilder::build_symbolic_entropy(&mut system, temperature).unwrap();

        let closure_value = system.dS_fun_view().get(&phase_key).unwrap()["H2O"](
            temperature,
            Some(vec![1.0]),
            Some(1.0),
        );
        let symbolic_value = system.dS_sym_view().get(&phase_key).unwrap()["H2O"]
            .clone()
            .lambdify_borrowed_thread_safe(&["T"])(&[temperature]);

        assert_relative_eq!(closure_value, numeric_value, epsilon = 1e-6);
        assert_relative_eq!(symbolic_value, numeric_value, epsilon = 1e-6);
    }

    #[test]
    fn test_one_phase_element_snapshot_is_pure_and_layout_aligned() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.set_substances(vec!["O2".to_string(), "H2".to_string()]);

        let (matrix, molar_masses, elements) =
            one.calculate_elem_composition_and_molar_mass(None).unwrap();

        assert_eq!(matrix.nrows(), 2);
        assert!(matrix.ncols() >= 2);
        assert!(molar_masses.contains_key("O2"));
        assert!(molar_masses.contains_key("H2"));
        assert!(elements.contains(&"O".to_string()));
        assert!(elements.contains(&"H".to_string()));
        assert!(one.subs_data_view().element_composition_matrix().is_none());
    }

    #[test]
    fn test_multi_phase_symbolic_layout_view_returns_atomic_bundle() {
        let pos = PhaseOrSolution::new();
        let layout = pos.symbolic_layout_view();
        assert!(layout.is_empty());
        assert!(pos.current_layout().components().is_empty());
    }

    #[test]
    fn test_single_phase_result_snapshot_carries_layout_and_conditions() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.set_substances(vec!["CO2".to_string(), "H2O".to_string()]);
        let snapshot = one.result_snapshot(Some(1200.0), Some(101325.0));
        assert_eq!(snapshot.layout_revision, one.layout_revision());
        assert_eq!(snapshot.temperature, Some(1200.0));
        assert_eq!(snapshot.pressure, Some(101325.0));
        assert_eq!(snapshot.phase_count(), 1);
        assert_eq!(snapshot.component_count(), 2);
        assert_eq!(
            snapshot.layout.component_labels(),
            vec!["CO2".to_string(), "H2O".to_string()]
        );
    }

    #[test]
    fn test_multi_phase_result_snapshot_carries_layout_and_conditions() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "H2".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.insert_phase_data("liquid".to_string(), liquid);
        let snapshot = pos.result_snapshot(Some(900.0), Some(2.5e5));
        assert_eq!(snapshot.layout_revision, pos.layout_revision());
        assert_eq!(snapshot.temperature, Some(900.0));
        assert_eq!(snapshot.pressure, Some(2.5e5));
        assert_eq!(snapshot.phase_count(), 2);
        assert_eq!(snapshot.component_count(), 3);
        assert_eq!(
            snapshot.layout.component_labels(),
            vec![
                "gas::O2".to_string(),
                "gas::H2".to_string(),
                "liquid::H2O".to_string()
            ]
        );
    }

    #[test]
    fn test_same_substance_in_two_phases_produces_two_distinct_components() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["H2".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.insert_phase_data("liquid".to_string(), liquid);

        let layout = crate::Thermodynamics::phase_layout::SystemLayout::from_phase_map(
            pos.phase_data_view(),
        );
        assert_eq!(layout.component_count(), 2);
        assert_eq!(
            layout.component_labels(),
            vec!["gas::H2".to_string(), "liquid::H2".to_string()]
        );

        let (indexed, vars, np_vec, map_each) = pos.indexed_moles_variables().unwrap();
        assert_eq!(
            vars,
            vec![Expr::Var("n0_0".to_string()), Expr::Var("n1_0".to_string())]
        );
        assert_eq!(
            np_vec,
            vec![Expr::Var("Np0".to_string()), Expr::Var("Np1".to_string())]
        );
        assert_eq!(
            indexed.get(&Some("gas".to_string())).unwrap(),
            &(
                Some(Expr::Var("Np0".to_string())),
                Some(vec![Expr::Var("n0_0".to_string())])
            )
        );
        assert_eq!(
            indexed.get(&Some("liquid".to_string())).unwrap(),
            &(
                Some(Expr::Var("Np1".to_string())),
                Some(vec![Expr::Var("n1_0".to_string())])
            )
        );
        assert_eq!(
            map_each.get(&Some("gas".to_string())).unwrap().get("H2"),
            Some(&Expr::Var("n0_0".to_string()))
        );
        assert_eq!(
            map_each.get(&Some("liquid".to_string())).unwrap().get("H2"),
            Some(&Expr::Var("n1_0".to_string()))
        );
    }

    #[test]
    fn test_resolved_shared_substance_keeps_phase_qualified_identity_and_totals() {
        // H2O appears in two physical phases. It is one chemical composition,
        // but two independent solver components with different phase models.
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["H2O".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            (
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            ),
            (
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let mut system = match spec.resolve().unwrap() {
            CustomSubstance::PhaseOrSolution(system) => system,
            _ => panic!("two named phases must resolve to PhaseOrSolution"),
        };

        assert_eq!(system.phase_specs().len(), 2);
        assert_eq!(system.phase_specs()[0].model(), PhaseModel::IdealGas);
        assert_eq!(system.phase_specs()[1].model(), PhaseModel::PureCondensed);
        assert_eq!(
            system.resolved_layout().unwrap().component_labels(),
            vec!["gas::H2O".to_string(), "liquid::H2O".to_string()]
        );

        let (_, vectors, totals) = system
            .create_full_map_of_mole_numbers(HashMap::from([
                (
                    Some("gas".to_string()),
                    (Some(0.25), Some(HashMap::from([("H2O".to_string(), 0.25)]))),
                ),
                (
                    Some("liquid".to_string()),
                    (Some(0.75), Some(HashMap::from([("H2O".to_string(), 0.75)]))),
                ),
            ]))
            .unwrap();
        assert_eq!(
            vectors[&Some("gas".to_string())].1.as_ref().unwrap(),
            &vec![0.25]
        );
        assert_eq!(
            vectors[&Some("liquid".to_string())].1.as_ref().unwrap(),
            &vec![0.75]
        );
        assert_relative_eq!(totals["H2O"], 1.0, epsilon = 1e-12);

        let (_, variables, phase_totals, _) = system.indexed_moles_variables().unwrap();
        assert_eq!(
            variables,
            vec![Expr::Var("n0_0".to_string()), Expr::Var("n1_0".to_string())]
        );
        assert_eq!(
            phase_totals,
            vec![Expr::Var("Np0".to_string()), Expr::Var("Np1".to_string())]
        );
    }

    #[test]
    fn test_substance_system_spec_builder_resolves_multi_phase() {
        let mut phases = HashMap::new();
        phases.insert("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]);
        phases.insert("liquid".to_string(), vec!["H2O".to_string()]);

        let spec = crate::Thermodynamics::User_PhaseOrSolution::SubstanceSystemSpec::builder(
            SubstancesContainer::MultiPhase(phases),
        )
        .with_phase_natures(Some(HashMap::from([
            (
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            ),
            (
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();

        let resolved = spec.resolve().unwrap();
        match resolved {
            CustomSubstance::PhaseOrSolution(pos) => {
                assert!(pos.phase_data_view().contains_key(&Some("gas".to_string())));
                assert!(
                    pos.phase_data_view()
                        .contains_key(&Some("liquid".to_string()))
                );
                assert_eq!(pos.phase_specs().len(), 2);
                let gas = pos
                    .phase_specs()
                    .iter()
                    .find(|phase| phase.id().as_option().as_deref() == Some("gas"))
                    .expect("resolved gas phase missing");
                let liquid = pos
                    .phase_specs()
                    .iter()
                    .find(|phase| phase.id().as_option().as_deref() == Some("liquid"))
                    .expect("resolved liquid phase missing");
                assert_eq!(gas.model(), PhaseModel::IdealGas);
                assert_eq!(liquid.model(), PhaseModel::PureCondensed);
            }
            _ => panic!("expected multi-phase system"),
        }
    }

    #[test]
    fn test_phase_resolution_reuses_one_injected_repository_for_all_phase_payloads() {
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            (
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            ),
            (
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let repository =
            crate::Thermodynamics::thermo_lib_api::ThermoData::try_default_repository()
                .expect("bundled test catalog must be available");

        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            spec,
            Arc::clone(&repository),
        )
        .unwrap();
        let gas_repository = &resolved
            .phase_data()
            .get(&Some("gas".to_string()))
            .expect("resolved gas payload missing")
            .thermo_data()
            .repository;
        let liquid_repository = &resolved
            .phase_data()
            .get(&Some("liquid".to_string()))
            .expect("resolved liquid payload missing")
            .thermo_data()
            .repository;

        assert!(Arc::ptr_eq(&repository, gas_repository));
        assert!(Arc::ptr_eq(gas_repository, liquid_repository));

        let report = resolved.report();
        assert!(!report.nist_fallback_enabled());
        assert_eq!(report.phases().len(), 2);
        let gas_summary = report
            .phase(&PhaseId::new(Some("gas".to_string())))
            .expect("gas provenance missing from resolved-system report");
        assert!(gas_summary.search().rows().iter().any(|row| {
            row.substance() == "H2"
                && row.property() == "Thermo"
                && row.state() == "Found"
                && row.library() == "NASA_gas"
                && row.priority() == "Priority"
        }));
    }

    #[test]
    fn test_phase_resolution_reports_mixed_local_library_provenance_offline() {
        // The fixture intentionally has no network path: each phase component
        // exists in exactly one local NASA-compatible library. Resolution must
        // preserve that mixed provenance in the facade report.
        let temp_dir = tempdir().unwrap();
        let all_keys_path = temp_dir.path().join("all_keys.json");
        let substance_base_path = temp_dir.path().join("substance_base.json");
        let elements_path = temp_dir.path().join("elements.json");
        let nasa_coefficients = serde_json::json!([
            200.0, 1000.0, 6000.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0
        ]);

        std::fs::write(
            &all_keys_path,
            serde_json::json!([["NASA_gas", "H2"], ["NASA_cond", "H2O"]]).to_string(),
        )
        .unwrap();
        std::fs::write(
            &substance_base_path,
            serde_json::json!({
                "NASA_gas": {
                    "H2": {
                        "Cp": nasa_coefficients,
                        "composition": "{H: 2}",
                        "model": "NASA7"
                    }
                },
                "NASA_cond": {
                    "H2O": {
                        "Cp": [200.0, 1000.0, 6000.0, 4.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        "composition": "{H: 2, O: 1}",
                        "model": "NASA7"
                    }
                }
            })
            .to_string(),
        )
        .unwrap();
        std::fs::write(
            &elements_path,
            serde_json::json!({
                "H": [["H2", "NASA_gas"], ["H2O", "NASA_cond"]],
                "O": [["H2O", "NASA_cond"]]
            })
            .to_string(),
        )
        .unwrap();

        let repository = crate::Thermodynamics::thermo_lib_api::ThermoData::try_new_from_paths(
            all_keys_path.to_str().unwrap(),
            substance_base_path.to_str().unwrap(),
            elements_path.to_str().unwrap(),
        )
        .unwrap();
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["H2".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            (
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            ),
            (
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();

        let resolved = SubstanceSystemFactory::resolve_spec_with_repository(
            spec,
            Arc::clone(&repository.repository),
        )
        .unwrap();
        let system = match resolved {
            CustomSubstance::PhaseOrSolution(system) => system,
            _ => panic!("named phases must resolve to the multi-phase facade"),
        };
        let report = system
            .resolution_report()
            .expect("resolved facade must retain local lookup provenance");

        assert!(!report.nist_fallback_enabled());
        assert!(
            report
                .phase(&PhaseId::new(Some("gas".to_string())))
                .unwrap()
                .search()
                .rows()
                .iter()
                .any(|row| row.substance() == "H2" && row.library() == "NASA_gas")
        );
        assert!(
            report
                .phase(&PhaseId::new(Some("liquid".to_string())))
                .unwrap()
                .search()
                .rows()
                .iter()
                .any(|row| row.substance() == "H2O" && row.library() == "NASA_cond")
        );
    }

    #[test]
    fn test_resolved_facade_keeps_lookup_provenance_after_resolution() {
        let spec = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            (
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            ),
            (
                "liquid".to_string(),
                crate::Thermodynamics::User_substances::Phases::Liquid,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_permitted_libraries(vec!["NIST".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let repository =
            crate::Thermodynamics::thermo_lib_api::ThermoData::try_default_repository()
                .expect("bundled test catalog must be available");

        let resolved =
            SubstanceSystemFactory::resolve_spec_with_repository(spec, Arc::clone(&repository))
                .unwrap();
        match resolved {
            CustomSubstance::PhaseOrSolution(system) => {
                let report = system
                    .resolution_report()
                    .expect("resolved facade must keep lookup provenance");
                assert!(!report.nist_fallback_enabled());
                assert_eq!(report.phases().len(), 2);
                let gas_summary = report
                    .phase(&PhaseId::new(Some("gas".to_string())))
                    .expect("gas provenance missing from facade report");
                assert!(gas_summary.search().rows().iter().any(|row| {
                    row.substance() == "H2"
                        && row.property() == "Thermo"
                        && row.state() == "Found"
                        && row.library() == "NASA_gas"
                }));
                let snapshot = system.result_snapshot(Some(300.0), Some(101_325.0));
                assert_eq!(snapshot.phase_count(), 2);
                assert_eq!(
                    system
                        .resolution_report()
                        .expect("report must survive result snapshot")
                        .phase(&PhaseId::new(Some("gas".to_string())))
                        .is_some(),
                    true
                );
            }
            _ => panic!("expected multi-phase system"),
        }
    }

    #[test]
    fn test_substance_system_spec_rejects_missing_phase_nature_for_declared_phase() {
        let mut phases = HashMap::new();
        phases.insert("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]);
        phases.insert("liquid".to_string(), vec!["H2O".to_string()]);

        let err = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(phases))
            .with_phase_natures(Some(HashMap::from([(
                "gas".to_string(),
                crate::Thermodynamics::User_substances::Phases::Gas,
            )])))
            .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
            .with_permitted_libraries(vec!["NIST".to_string()])
            .with_search_in_nist(false)
            .build()
            .unwrap_err();

        assert!(
            err.to_string()
                .contains("missing physical state for declared phase")
        );
    }

    #[test]
    fn test_substance_system_spec_rejects_phase_nature_for_unknown_phase() {
        let phases = HashMap::from([("gas".to_string(), vec!["H2".to_string(), "O2".to_string()])]);

        let err = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(phases))
            .with_phase_natures(Some(HashMap::from([
                (
                    "gas".to_string(),
                    crate::Thermodynamics::User_substances::Phases::Gas,
                ),
                (
                    "gas_typo".to_string(),
                    crate::Thermodynamics::User_substances::Phases::Gas,
                ),
            ])))
            .build()
            .expect_err("phase states for unknown phases must be rejected");

        assert!(
            err.to_string()
                .contains("physical state declared for unknown phase 'gas_typo'")
        );
    }

    #[test]
    fn test_substance_system_spec_rejects_implicit_multi_phase_states() {
        let phases = HashMap::from([
            ("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
        ]);

        let error = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(phases))
            .build()
            .expect_err("multi-phase input without physical states must be rejected");

        assert!(
            error
                .to_string()
                .contains("require an explicit physical state")
        );
    }

    #[test]
    fn test_substance_system_spec_rejects_empty_multi_phase() {
        let err = SubstanceSystemSpec::builder(SubstancesContainer::MultiPhase(HashMap::new()))
            .build()
            .unwrap_err();

        assert!(format!("{}", err).contains("a phase system must contain at least one phase"));
    }

    #[test]
    fn test_phase_spec_keeps_state_and_model_as_separate_contracts() {
        let gas = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string(), "N2".to_string()],
        )
        .unwrap();
        let solid = PhaseSpec::pure_condensed(
            PhaseId::new(Some("solid".to_string())),
            vec!["C(s)".to_string()],
            PhasePhysicalState::Solid,
        )
        .unwrap();

        let spec = SubstanceSystemSpec::from_phases(vec![gas.clone(), solid.clone()]).unwrap();
        assert_eq!(spec.phases(), &[gas, solid]);
        assert_eq!(spec.phases()[0].model(), PhaseModel::IdealGas);
        assert_eq!(spec.phases()[1].physical_state(), PhasePhysicalState::Solid);

        let invalid = PhaseSpec::new(
            PhaseId::new(Some("liquid".to_string())),
            vec!["H2O".to_string()],
            PhasePhysicalState::Liquid,
            PhaseModel::IdealGas,
        )
        .unwrap_err();
        assert!(invalid.to_string().contains("not implemented"));
    }

    #[test]
    fn test_phase_spec_derives_legacy_component_states_from_one_canonical_state() {
        let phase = PhaseSpec::pure_condensed(
            PhaseId::new(Some("solid".to_string())),
            vec!["C(s)".to_string(), "Si(s)".to_string()],
            PhasePhysicalState::Solid,
        )
        .unwrap();

        let legacy_states = phase.legacy_component_phase_map();
        assert_eq!(legacy_states.len(), 2);
        assert!(matches!(
            legacy_states.get("C(s)"),
            Some(Some(crate::Thermodynamics::User_substances::Phases::Solid))
        ));
        assert!(matches!(
            legacy_states.get("Si(s)"),
            Some(Some(crate::Thermodynamics::User_substances::Phases::Solid))
        ));
    }

    #[test]
    fn test_phase_spec_rejects_duplicate_components_and_duplicate_phase_ids() {
        let duplicate_component = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string(), "O2".to_string()],
        )
        .unwrap_err();
        assert!(duplicate_component.to_string().contains("duplicated"));

        let phase_a = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string()],
        )
        .unwrap();
        let phase_b = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["N2".to_string()],
        )
        .unwrap();
        let duplicate_phase = SubstanceSystemSpec::from_phases(vec![phase_a, phase_b]).unwrap_err();
        assert!(
            duplicate_phase
                .to_string()
                .contains("declared more than once")
        );
    }

    #[test]
    fn test_resolved_phase_system_aligns_specs_records_and_layout() {
        let gas = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string(), "N2".to_string()],
        )
        .unwrap();
        let solid = PhaseSpec::pure_condensed(
            PhaseId::new(Some("solid".to_string())),
            vec!["C(s)".to_string()],
            PhasePhysicalState::Solid,
        )
        .unwrap();
        let mut gas_data = SubsData::new();
        gas_data.substances = gas.components().to_vec();
        let mut solid_data = SubsData::new();
        solid_data.substances = solid.components().to_vec();

        let resolved = ResolvedPhaseSystem::new(
            vec![solid, gas],
            HashMap::from([
                (Some("gas".to_string()), gas_data),
                (Some("solid".to_string()), solid_data),
            ]),
        )
        .unwrap();
        assert_eq!(resolved.phase_specs().len(), 2);
        assert_eq!(
            resolved.layout().component_labels(),
            vec![
                "gas::O2".to_string(),
                "gas::N2".to_string(),
                "solid::C(s)".to_string(),
            ]
        );

        let wrong_order = ResolvedPhaseSystem::new(
            vec![
                PhaseSpec::ideal_gas(
                    PhaseId::new(Some("gas".to_string())),
                    vec!["O2".to_string(), "N2".to_string()],
                )
                .unwrap(),
            ],
            HashMap::from([(Some("gas".to_string()), {
                let mut data = SubsData::new();
                data.substances = vec!["N2".to_string(), "O2".to_string()];
                data
            })]),
        )
        .unwrap_err();
        assert!(wrong_order.to_string().contains("component order"));
    }

    #[test]
    fn test_resolved_multiphase_install_publishes_payload_layout_and_report_together() {
        let gas = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string()],
        )
        .unwrap();
        let liquid = PhaseSpec::pure_condensed(
            PhaseId::new(Some("liquid".to_string())),
            vec!["H2O".to_string()],
            PhasePhysicalState::Liquid,
        )
        .unwrap();
        let mut gas_data = SubsData::new();
        gas_data.substances = gas.components().to_vec();
        let mut liquid_data = SubsData::new();
        liquid_data.substances = liquid.components().to_vec();
        let resolved = ResolvedPhaseSystem::new(
            vec![gas, liquid],
            HashMap::from([
                (Some("gas".to_string()), gas_data),
                (Some("liquid".to_string()), liquid_data),
            ]),
        )
        .unwrap();

        let mut system = PhaseOrSolution::new();
        let revision_before = system.layout_revision();
        system.install_resolved_system(resolved);

        assert_eq!(system.phase_specs().len(), 2);
        assert_eq!(system.phase_data_view().len(), 2);
        assert_eq!(
            system
                .resolved_layout()
                .expect("install must publish the canonical layout")
                .component_labels(),
            vec!["gas::O2".to_string(), "liquid::H2O".to_string()]
        );
        assert!(system.resolution_report().is_some());
        assert!(system.layout_revision() > revision_before);
        assert_eq!(system.dG_view().len(), 2);
        assert!(system.dG_view().get(&Some("gas".to_string())).is_some());
        assert!(system.dG_view().get(&Some("liquid".to_string())).is_some());
    }

    #[test]
    fn test_resolved_single_phase_install_publishes_the_same_complete_state() {
        let spec = PhaseSpec::ideal_gas(PhaseId::new(None), vec!["CO2".to_string()]).unwrap();
        let mut data = SubsData::new();
        data.substances = spec.components().to_vec();
        let resolved = ResolvedPhaseSystem::new(vec![spec], HashMap::from([(None, data)])).unwrap();

        let mut system = OnePhase::new();
        let revision_before = system.layout_revision();
        system.install_resolved_system(resolved).unwrap();

        assert_eq!(system.phase_specs().len(), 1);
        assert_eq!(system.subs_data_view().substances, vec!["CO2".to_string()]);
        assert_eq!(
            system
                .resolved_layout()
                .expect("install must publish the canonical layout")
                .component_labels(),
            vec!["CO2".to_string()]
        );
        assert!(system.resolution_report().is_some());
        assert!(system.layout_revision() > revision_before);
    }

    #[test]
    fn test_borrowed_phase_cache_iteration_distinguishes_layout_from_population() {
        let mut system = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        system.insert_phase_data("gas".to_string(), gas);
        system.insert_phase_data("liquid".to_string(), liquid);
        system.debug_replace_dG(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("O2".to_string(), -12.0)]),
        )]));

        let cache = system.dG_view();
        assert_eq!(cache.phase_count(), 2);
        assert_eq!(cache.populated_phase_count(), 1);
        let borrowed_phase_names = cache
            .iter_ref()
            .map(|(phase, _)| phase.clone())
            .collect::<std::collections::HashSet<_>>();
        assert_eq!(
            borrowed_phase_names,
            HashSet::from([Some("gas".to_string()), Some("liquid".to_string())])
        );

        let nested = NestedPhaseCacheView::Multi(cache);
        assert_eq!(nested.phase_count(), 2);
        assert_eq!(nested.populated_phase_count(), 1);
        let nested_phase_names = nested
            .iter_ref()
            .map(|(phase, _)| phase.clone())
            .collect::<HashSet<_>>();
        assert_eq!(nested_phase_names, borrowed_phase_names);
    }

    #[test]
    fn test_single_phase_nested_cache_view_retains_its_canonical_none_phase() {
        let mut system = OnePhase::new();

        let empty = NestedPhaseCacheView::Single(system.dG_view());
        assert_eq!(empty.phase_count(), 1);
        assert_eq!(empty.populated_phase_count(), 0);
        assert!(empty.iter_ref().all(|(phase, _)| phase.is_none()));

        system.debug_dG_mut().insert("H2O".to_string(), -18.0);
        let populated = NestedPhaseCacheView::Single(system.dG_view());
        assert_eq!(populated.phase_count(), 1);
        assert_eq!(populated.populated_phase_count(), 1);
        let (phase, values) = populated
            .iter_ref()
            .next()
            .expect("single-phase view must yield its canonical phase");
        assert!(phase.is_none());
        assert_eq!(values.get("H2O"), Some(&-18.0));
    }

    #[test]
    fn test_public_mapping_mutation_clears_all_visible_derived_caches() {
        let mut system = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string()];
        system.replace_phase_data(HashMap::from([(Some("gas".to_string()), gas)]));
        system.debug_replace_dG(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("O2".to_string(), -10.0)]),
        )]));
        let revision_before = system.layout_revision();
        assert_eq!(
            system
                .dG_view()
                .get(&Some("gas".to_string()))
                .and_then(|values| values.get("O2")),
            Some(&-10.0)
        );

        let mapping =
            SubstancePhaseMapping::new(vec!["O2".to_string()], vec![Some("gas".to_string())])
                .unwrap();
        system.set_substance_phase_mapping(mapping);

        assert!(system.layout_revision() > revision_before);
        assert!(system.gibbs_cache_context().is_none());
        assert!(system.entropy_cache_context().is_none());
        assert!(
            system
                .dG_view()
                .get(&Some("gas".to_string()))
                .expect("canonical phase bundle remains present")
                .is_empty()
        );
    }

    #[test]
    fn test_one_phase_rejects_multiphase_resolution_without_state_publication() {
        let gas = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string()],
        )
        .unwrap();
        let liquid = PhaseSpec::pure_condensed(
            PhaseId::new(Some("liquid".to_string())),
            vec!["H2O".to_string()],
            PhasePhysicalState::Liquid,
        )
        .unwrap();
        let mut gas_data = SubsData::new();
        gas_data.substances = gas.components().to_vec();
        let mut liquid_data = SubsData::new();
        liquid_data.substances = liquid.components().to_vec();
        let resolved = ResolvedPhaseSystem::new(
            vec![gas, liquid],
            HashMap::from([
                (Some("gas".to_string()), gas_data),
                (Some("liquid".to_string()), liquid_data),
            ]),
        )
        .unwrap();

        let mut one = OnePhase::new();
        let revision_before = one.layout_revision();
        let error = one
            .install_resolved_system(resolved)
            .expect_err("OnePhase must reject a multi-phase resolved payload");

        assert!(error.to_string().contains("exactly one unnamed"));
        assert_eq!(one.layout_revision(), revision_before);
        assert!(one.phase_specs().is_empty());
        assert!(one.resolved_layout().is_none());
        assert!(one.resolution_report().is_none());
        assert!(one.dG_view().is_empty());
    }

    #[test]
    fn test_create_full_map_of_mole_numbers_phase_or_solution() {
        // Build PhaseOrSolution with two phases
        let mut pos = PhaseOrSolution::new();
        let mut s1 = SubsData::new();
        s1.substances = vec!["A".to_string(), "B".to_string()];
        let mut s2 = SubsData::new();
        s2.substances = vec!["B".to_string(), "C".to_string()];
        pos.insert_phase_data("p1".to_string(), s1);
        pos.insert_phase_data("p2".to_string(), s2);

        // Provide non-zero mole numbers only for p1 with only A present
        let mut input: HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)> =
            HashMap::new();
        input.insert(
            Some("p1".to_string()),
            (Some(1.0), Some(HashMap::from([("A".to_string(), 0.3)]))),
        );

        let (full_map, vec_map, summed) =
            pos.create_full_map_of_mole_numbers(input.clone()).unwrap();

        // Check missing keys inserted for p1 (B should be present with 0.0)
        let p1_inner = full_map.get(&Some("p1".to_string())).unwrap();
        assert!(p1_inner.1.as_ref().unwrap().contains_key("B"));
        assert_eq!(p1_inner.1.as_ref().unwrap().get("B").unwrap(), &0.0);

        // Check vec conversion length
        let p1_vec = vec_map.get(&Some("p1".to_string())).unwrap();
        assert!(p1_vec.1.as_ref().unwrap().len() == 2);

        // Check summed map: A=0.3, B=0.0, C absent -> B should be 0.0 present only if present in at least one phase
        assert!(summed.contains_key("A"));
        assert_eq!(summed.get("A").unwrap(), &0.3);
    }

    #[test]
    fn test_indexed_moles_variables_and_symbolic_vars() {
        let mut pos = PhaseOrSolution::new();
        let mut s1 = SubsData::new();
        s1.substances = vec!["X".to_string(), "Y".to_string()];
        let mut s2 = SubsData::new();
        s2.substances = vec!["Z".to_string()];
        pos.insert_phase_data("p1".to_string(), s1);
        pos.insert_phase_data("p2".to_string(), s2);

        let (_map_indexed, vec_of_n_vars, np_vec, map_of_var_each) =
            pos.indexed_moles_variables().unwrap();

        // Np vector length equals number of phases
        assert_eq!(np_vec.len(), 2);
        // total number of n vars equals 3
        assert_eq!(vec_of_n_vars.len(), 3);
        // Check map_of_var_each_substance contains entries
        let p1_map = map_of_var_each.get(&Some("p1".to_string()));
        assert!(p1_map.is_some());
        assert!(p1_map.unwrap().contains_key("X"));
        assert!(p1_map.unwrap().contains_key("Y"));
    }

    #[test]
    fn test_indexed_moles_variables_use_stable_phase_order_and_unique_names() {
        let mut pos_a = PhaseOrSolution::new();
        let mut phase_b = SubsData::new();
        phase_b.substances = vec!["B1".to_string()];
        let mut phase_a = SubsData::new();
        phase_a.substances = vec!["A1".to_string(), "A2".to_string()];
        pos_a.insert_phase_data("b".to_string(), phase_b.clone());
        pos_a.insert_phase_data("a".to_string(), phase_a.clone());

        let mut pos_b = PhaseOrSolution::new();
        pos_b.insert_phase_data("a".to_string(), phase_a);
        pos_b.insert_phase_data("b".to_string(), phase_b);

        let (_, vars_a, np_a, map_a) = pos_a.indexed_moles_variables().unwrap();
        let (_, vars_b, np_b, map_b) = pos_b.indexed_moles_variables().unwrap();

        assert_eq!(vars_a, vars_b);
        assert_eq!(np_a, np_b);
        assert_eq!(
            vars_a,
            vec!["n0_0", "n0_1", "n1_0"]
                .into_iter()
                .map(|s| Expr::Var(s.to_string()))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            np_a,
            vec![Expr::Var("Np0".to_string()), Expr::Var("Np1".to_string())]
        );
        assert!(
            map_a
                .get(&Some("a".to_string()))
                .unwrap()
                .contains_key("A1")
        );
        assert!(
            map_a
                .get(&Some("b".to_string()))
                .unwrap()
                .contains_key("B1")
        );
        assert_eq!(map_a, map_b);
    }

    #[test]
    fn test_result_snapshot_projects_numeric_caches_into_canonical_dense_order() {
        let mut system = PhaseOrSolution::new();
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "N2".to_string()];
        // Insert in non-canonical order: the snapshot must still follow gas,
        // then liquid, and preserve the declared component sequence per phase.
        system.insert_phase_data("liquid".to_string(), liquid);
        system.insert_phase_data("gas".to_string(), gas);
        system.debug_replace_dG(HashMap::from([
            (
                Some("liquid".to_string()),
                HashMap::from([("H2O".to_string(), -241.8)]),
            ),
            (
                Some("gas".to_string()),
                HashMap::from([("N2".to_string(), -10.0), ("O2".to_string(), -20.0)]),
            ),
        ]));
        system.debug_replace_dS(HashMap::from([
            (
                Some("liquid".to_string()),
                HashMap::from([("H2O".to_string(), 69.9)]),
            ),
            (
                Some("gas".to_string()),
                HashMap::from([("N2".to_string(), 191.5), ("O2".to_string(), 205.1)]),
            ),
        ]));

        let snapshot = system.result_snapshot(Some(298.15), Some(101_325.0));
        assert_eq!(
            snapshot.component_labels(),
            vec![
                "gas::O2".to_string(),
                "gas::N2".to_string(),
                "liquid::H2O".to_string(),
            ]
        );
        assert_eq!(
            snapshot.ordered_gibbs_values().unwrap(),
            vec![-20.0, -10.0, -241.8]
        );
        assert_eq!(
            snapshot.ordered_entropy_values().unwrap(),
            vec![205.1, 191.5, 69.9]
        );
    }

    #[test]
    fn test_result_snapshot_rejects_missing_numeric_component_cache_entry() {
        let mut system = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "N2".to_string()];
        system.insert_phase_data("gas".to_string(), gas);
        system.debug_replace_dG(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("O2".to_string(), -20.0)]),
        )]));

        let error = system
            .result_snapshot(Some(298.15), Some(101_325.0))
            .ordered_gibbs_values()
            .expect_err("every layout component requires a numeric cache value");
        assert!(error.to_string().contains("gas::N2"));
    }

    #[test]
    fn test_result_snapshot_rejects_non_finite_numeric_cache_values() {
        let mut system = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string()];
        system.insert_phase_data("gas".to_string(), gas);
        system.debug_replace_dG(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("O2".to_string(), f64::NAN)]),
        )]));

        let error = system
            .result_snapshot(Some(298.15), Some(101_325.0))
            .ordered_gibbs_values()
            .expect_err("a dense solver vector must reject NaN cache values");
        assert!(error.to_string().contains("non-finite"));
        assert!(error.to_string().contains("gas::O2"));
    }

    #[test]
    fn test_three_phase_layout_and_mole_normalization_ignore_input_map_order() {
        // The phase engine has a canonical ordered layout. This matrix keeps
        // external HashMap insertion order from leaking into solver variables
        // or component mole vectors.
        let phase_orders = [
            ["a", "b", "c"],
            ["a", "c", "b"],
            ["b", "a", "c"],
            ["b", "c", "a"],
            ["c", "a", "b"],
            ["c", "b", "a"],
        ];

        for phase_order in phase_orders {
            let mut system = PhaseOrSolution::new();
            for phase_name in phase_order {
                let mut data = SubsData::new();
                data.substances = match phase_name {
                    "a" => vec!["A1".to_string(), "A2".to_string()],
                    "b" => vec!["B1".to_string(), "B2".to_string()],
                    "c" => vec!["C1".to_string(), "C2".to_string()],
                    _ => unreachable!("the test matrix contains only declared phases"),
                };
                system.insert_phase_data(phase_name.to_string(), data);
            }

            let (_, mole_variables, phase_totals, _) = system.indexed_moles_variables().unwrap();
            assert_eq!(
                mole_variables,
                (0..3)
                    .flat_map(|phase_index| {
                        (0..2).map(move |component_index| {
                            Expr::Var(format!("n{phase_index}_{component_index}"))
                        })
                    })
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                phase_totals,
                (0..3)
                    .map(|phase_index| Expr::Var(format!("Np{phase_index}")))
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                PhaseLayoutAccess::system_layout(&system).component_labels(),
                vec![
                    "a::A1".to_string(),
                    "a::A2".to_string(),
                    "b::B1".to_string(),
                    "b::B2".to_string(),
                    "c::C1".to_string(),
                    "c::C2".to_string(),
                ]
            );

            let mole_input = HashMap::from([
                (
                    Some("c".to_string()),
                    (
                        Some(1.0),
                        Some(HashMap::from([
                            ("C2".to_string(), 0.6),
                            ("C1".to_string(), 0.4),
                        ])),
                    ),
                ),
                (
                    Some("b".to_string()),
                    (
                        Some(1.0),
                        Some(HashMap::from([
                            ("B2".to_string(), 0.8),
                            ("B1".to_string(), 0.2),
                        ])),
                    ),
                ),
                (
                    Some("a".to_string()),
                    (
                        Some(1.0),
                        Some(HashMap::from([
                            ("A2".to_string(), 0.7),
                            ("A1".to_string(), 0.3),
                        ])),
                    ),
                ),
            ]);
            let (_, normalized, _) = system.create_full_map_of_mole_numbers(mole_input).unwrap();
            assert_eq!(
                normalized.get(&Some("a".to_string())).unwrap().1,
                Some(vec![0.3, 0.7])
            );
            assert_eq!(
                normalized.get(&Some("b".to_string())).unwrap().1,
                Some(vec![0.2, 0.8])
            );
            assert_eq!(
                normalized.get(&Some("c".to_string())).unwrap().1,
                Some(vec![0.4, 0.6])
            );
        }
    }

    #[test]
    fn test_create_full_map_of_mole_numbers_respects_declared_component_order() {
        let mut pos = PhaseOrSolution::new();
        let mut phase = SubsData::new();
        phase.substances = vec!["B".to_string(), "A".to_string()];
        pos.insert_phase_data("mix".to_string(), phase);

        let mut input: HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)> =
            HashMap::new();
        input.insert(
            Some("mix".to_string()),
            (Some(2.0), Some(HashMap::from([("A".to_string(), 1.5)]))),
        );

        let (_, vec_map, _) = pos.create_full_map_of_mole_numbers(input).unwrap();
        let ordered = vec_map.get(&Some("mix".to_string())).unwrap();
        assert_eq!(ordered.1.as_ref().unwrap(), &vec![0.0, 1.5]);
    }

    #[test]
    fn test_typed_mole_snapshot_carries_layout_and_preserves_sparse_input() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "N2".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);

        let snapshot = pos
            .normalize_mole_numbers(HashMap::from([(
                Some("gas".to_string()),
                (Some(1.0), Some(HashMap::from([("N2".to_string(), 0.8)]))),
            )]))
            .unwrap();

        assert_eq!(
            snapshot.layout().component_labels(),
            vec!["gas::O2".to_string(), "gas::N2".to_string()]
        );
        let phase = snapshot
            .phase_amounts()
            .get(&Some("gas".to_string()))
            .unwrap();
        assert_eq!(phase.component_amounts().unwrap()["O2"], 0.0);
        assert_eq!(phase.component_amounts().unwrap()["N2"], 0.8);
        assert_eq!(
            snapshot
                .ordered_phase_amounts()
                .get(&Some("gas".to_string()))
                .unwrap()
                .component_amounts(),
            &[0.0, 0.8]
        );
        assert_eq!(snapshot.total_amounts_by_substance()["N2"], 0.8);
    }

    #[test]
    fn test_container_views_split_catalog_and_solver_order() {
        let container = SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["O2".to_string(), "H2".to_string()]),
            (
                "condensed".to_string(),
                vec!["H2".to_string(), "H2O".to_string()],
            ),
        ]));

        let catalog = container.get_all_substances();
        assert_eq!(
            catalog,
            vec!["H2".to_string(), "H2O".to_string(), "O2".to_string()]
        );

        let ordered = container.get_ordered_component_labels();
        assert_eq!(
            ordered,
            vec![
                (Some("condensed".to_string()), "H2".to_string()),
                (Some("condensed".to_string()), "H2O".to_string()),
                (Some("gas".to_string()), "O2".to_string()),
                (Some("gas".to_string()), "H2".to_string()),
            ]
        );
        assert_eq!(
            container.get_ordered_substances(),
            vec![
                "H2".to_string(),
                "H2O".to_string(),
                "O2".to_string(),
                "H2".to_string()
            ]
        );
    }

    #[test]
    fn test_resolved_phase_system_has_stable_layout_for_permuted_input_order() {
        let gas = PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_string())),
            vec!["O2".to_string(), "N2".to_string()],
        )
        .unwrap();
        let liquid = PhaseSpec::pure_condensed(
            PhaseId::new(Some("liquid".to_string())),
            vec!["H2O".to_string()],
            PhasePhysicalState::Liquid,
        )
        .unwrap();
        let solid = PhaseSpec::pure_condensed(
            PhaseId::new(Some("solid".to_string())),
            vec!["C(s)".to_string()],
            PhasePhysicalState::Solid,
        )
        .unwrap();

        let phase_data = |entries: Vec<(&str, Vec<&str>)>| {
            entries
                .into_iter()
                .map(|(phase, components)| {
                    let mut data = SubsData::new();
                    data.substances = components.into_iter().map(str::to_string).collect();
                    (Some(phase.to_string()), data)
                })
                .collect::<HashMap<_, _>>()
        };

        let first = ResolvedPhaseSystem::new(
            vec![solid.clone(), gas.clone(), liquid.clone()],
            phase_data(vec![
                ("solid", vec!["C(s)"]),
                ("gas", vec!["O2", "N2"]),
                ("liquid", vec!["H2O"]),
            ]),
        )
        .unwrap();
        let second = ResolvedPhaseSystem::new(
            vec![liquid, gas, solid],
            phase_data(vec![
                ("liquid", vec!["H2O"]),
                ("gas", vec!["O2", "N2"]),
                ("solid", vec!["C(s)"]),
            ]),
        )
        .unwrap();

        assert_eq!(first.layout(), second.layout());
        assert_eq!(
            first.layout().component_labels(),
            vec![
                "gas::O2".to_string(),
                "gas::N2".to_string(),
                "liquid::H2O".to_string(),
                "solid::C(s)".to_string(),
            ]
        );
        assert_eq!(
            first
                .phase_specs()
                .iter()
                .map(|phase| phase.id().as_option().clone())
                .collect::<Vec<_>>(),
            vec![
                Some("gas".to_string()),
                Some("liquid".to_string()),
                Some("solid".to_string()),
            ]
        );
    }

    #[test]
    fn test_custom_substance_exposes_multiphase_structure_without_flattening() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "H2".to_string()];
        let mut condensed = SubsData::new();
        condensed.substances = vec!["H2".to_string(), "H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.insert_phase_data("condensed".to_string(), condensed);

        let custom = CustomSubstance::PhaseOrSolution(pos);
        let container = custom.extract_SubstancesContainer().unwrap();
        match container {
            SubstancesContainer::MultiPhase(map) => {
                assert_eq!(
                    map.get("gas").unwrap(),
                    &vec!["O2".to_string(), "H2".to_string()]
                );
                assert_eq!(
                    map.get("condensed").unwrap(),
                    &vec!["H2".to_string(), "H2O".to_string()]
                );
            }
            _ => panic!("expected multiphase container"),
        }

        let labels = custom.get_ordered_component_labels().unwrap();
        assert_eq!(
            labels,
            vec![
                (Some("condensed".to_string()), "H2".to_string()),
                (Some("condensed".to_string()), "H2O".to_string()),
                (Some("gas".to_string()), "O2".to_string()),
                (Some("gas".to_string()), "H2".to_string()),
            ]
        );
    }

    #[test]
    fn test_custom_substance_result_labels_follow_layout_order() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "H2".to_string()];
        let mut condensed = SubsData::new();
        condensed.substances = vec!["H2".to_string(), "H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.insert_phase_data("condensed".to_string(), condensed);

        let labels = CustomSubstance::PhaseOrSolution(pos)
            .get_ordered_result_labels()
            .unwrap();
        assert_eq!(
            labels,
            vec![
                "condensed::H2".to_string(),
                "condensed::H2O".to_string(),
                "gas::O2".to_string(),
                "gas::H2".to_string(),
            ]
        );
    }

    #[test]
    fn test_elem_composition_matrix_is_stable_under_phase_insertion_order() {
        let mut pos_a = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "H2".to_string()];
        let mut condensed = SubsData::new();
        condensed.substances = vec!["H2".to_string(), "H2O".to_string()];
        pos_a.insert_phase_data("gas".to_string(), gas.clone());
        pos_a.insert_phase_data("condensed".to_string(), condensed.clone());

        let mut pos_b = PhaseOrSolution::new();
        pos_b.insert_phase_data("condensed".to_string(), condensed);
        pos_b.insert_phase_data("gas".to_string(), gas);

        let (matrix_a, molar_a, elements_a) = pos_a
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        let (matrix_b, molar_b, elements_b) = pos_b
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        assert_eq!(elements_a, elements_b);
        assert_eq!(molar_a, molar_b);
        assert_eq!(matrix_a, matrix_b);
    }

    #[test]
    fn test_calculate_elem_composition_keeps_repeated_substance_rows_for_multiple_phases() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["H2O".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.insert_phase_data("liquid".to_string(), liquid);

        let (matrix, molar_masses, elements) =
            pos.calculate_elem_composition_and_molar_mass(None).unwrap();

        assert_eq!(matrix.nrows(), 2);
        assert!(matrix.ncols() >= 2);
        assert!(molar_masses.contains_key("H2O"));
        assert!(elements.contains(&"H".to_string()));
        assert!(elements.contains(&"O".to_string()));
    }

    #[test]
    fn test_calculate_lagrange_equations_sym_supports_same_substance_in_multiple_phases() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["H2O".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.insert_phase_data("liquid".to_string(), liquid);

        let mut gas_sym = HashMap::new();
        gas_sym.insert("H2O".to_string(), Expr::Const(1.0));
        let mut liquid_sym = HashMap::new();
        liquid_sym.insert("H2O".to_string(), Expr::Const(2.0));
        pos.debug_replace_dG_sym(HashMap::from([
            (Some("gas".to_string()), gas_sym),
            (Some("liquid".to_string()), liquid_sym),
        ]));

        let A = DMatrix::from_vec(2, 1, vec![1.0, 1.0]);
        let eqs = pos.calculate_Lagrange_equations_sym(A, 300.0).unwrap();
        assert_eq!(eqs.len(), 2);
        assert_ne!(format!("{}", eqs[0]), format!("{}", eqs[1]));
    }

    #[test]
    fn test_calculate_lagrange_equations_sym_rejects_wrong_matrix_shape() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.debug_replace_dG_sym(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("H2O".to_string(), Expr::Const(1.0))]),
        )]));
        let revision_before = pos.layout_revision();
        let cached_before = pos
            .dG_sym_view()
            .get(&Some("gas".to_string()))
            .and_then(|values| values.get("H2O"))
            .map(ToString::to_string)
            .expect("fixture must publish symbolic Gibbs data");

        let err = match pos
            .calculate_Lagrange_equations_sym(DMatrix::from_vec(2, 1, vec![1.0, 1.0]), 300.0)
        {
            Ok(_) => panic!("expected shape validation error"),
            Err(err) => err,
        };
        assert!(format!("{}", err).contains("component rows"));
        assert_eq!(pos.layout_revision(), revision_before);
        assert_eq!(
            pos.dG_sym_view()
                .get(&Some("gas".to_string()))
                .and_then(|values| values.get("H2O"))
                .map(ToString::to_string),
            Some(cached_before)
        );
    }

    #[test]
    fn test_calculate_lagrange_equations_fun_rejects_missing_gibbs_data() {
        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["H2O".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        let revision_before = pos.layout_revision();

        let err = match pos.calculate_Lagrange_equations_fun(
            DMatrix::from_vec(1, 1, vec![1.0]),
            HashMap::new(),
            300.0,
        ) {
            Ok(_) => panic!("expected Gibbs data validation error"),
            Err(err) => err,
        };
        assert!(format!("{}", err).contains("phase Gibbs functions"));
        assert_eq!(pos.layout_revision(), revision_before);
        assert!(pos.dG_fun_view().is_empty());
    }

    #[test]
    fn test_one_phase_lagrange_equations_validate_shape_and_missing_gibbs() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.set_substances(vec!["A".to_string(), "B".to_string()]);
        one.debug_dG_sym_mut()
            .insert("A".to_string(), Expr::Const(1.0));

        let err = match one
            .calculate_Lagrange_equations_sym(DMatrix::from_vec(1, 2, vec![1.0, 1.0]), 300.0)
        {
            Ok(_) => panic!("expected shape validation error"),
            Err(err) => err,
        };
        assert!(format!("{}", err).contains("component rows"));

        let mut gfun = HashMap::new();
        gfun.insert(
            "A".to_string(),
            Box::new(|_, _, _| 1.0) as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        let err = match one.calculate_Lagrange_equations_fun(
            DMatrix::from_vec(2, 1, vec![1.0, 1.0]),
            gfun,
            300.0,
        ) {
            Ok(_) => panic!("expected Gibbs data validation error"),
            Err(err) => err,
        };
        assert!(format!("{}", err).contains("Gibbs function"));
    }

    #[test]
    fn test_clone_preserves_function_outputs_for_one_phase_and_multiphase() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.debug_dG_fun_mut().insert(
            "A".to_string(),
            Arc::new(|_, _, _| 7.5) as Arc<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        one.debug_dS_fun_mut().insert(
            "A".to_string(),
            Arc::new(|_, _, _| -2.0) as Arc<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        let one_clone = one.clone();
        assert_eq!(one.dG_fun_view()["A"](300.0, None, None), 7.5);
        assert_eq!(one_clone.dG_fun_view()["A"](300.0, None, None), 7.5);
        assert_eq!(one.dS_fun_view()["A"](300.0, None, None), -2.0);
        assert_eq!(one_clone.dS_fun_view()["A"](300.0, None, None), -2.0);

        let mut pos = PhaseOrSolution::new();
        let mut phase = SubsData::new();
        phase.substances = vec!["A".to_string()];
        pos.insert_phase_data("p".to_string(), phase);
        let mut gfun = HashMap::new();
        gfun.insert(
            "A".to_string(),
            Arc::new(|_, _, _| 3.25) as Arc<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        pos.debug_replace_dG_fun(HashMap::from([(Some("p".to_string()), gfun.clone())]));
        pos.debug_replace_dS_fun(HashMap::from([(Some("p".to_string()), gfun)]));

        let pos_clone = pos.clone();
        assert!(
            pos.dG_fun_view()
                .get(&Some("missing".to_string()))
                .is_none(),
            "a cache view must represent an absent phase without panicking"
        );
        assert_eq!(
            pos.dG_fun_view()
                .get(&Some("p".to_string()))
                .and_then(|phase| phase.get("A"))
                .expect("Gibbs function for phase component must be cached")(
                300.0, None, None
            ),
            3.25
        );
        assert_eq!(
            pos_clone
                .dG_fun_view()
                .get(&Some("p".to_string()))
                .and_then(|phase| phase.get("A"))
                .expect("cloned Gibbs function for phase component must be cached")(
                300.0, None, None
            ),
            3.25
        );
    }

    #[test]
    fn test_transition_payload_accessors_clear_stale_layouts_and_caches() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        let mut one_data = SubsData::new();
        one_data.substances = vec!["CO2".to_string()];
        one.replace_subs_data(one_data);
        one.set_phase_spec(
            PhaseSpec::ideal_gas(PhaseId::new(None), vec!["CO2".to_string()]).unwrap(),
        );
        one.debug_dG_mut().insert("CO2".to_string(), -1.0);
        let one_revision = one.layout_revision();

        one.with_subs_data_mut(|data| data.substances.push("H2O".to_string()));

        assert!(one.phase_spec().is_none());
        assert!(one.resolved_layout().is_none());
        assert!(one.dG_view().is_empty());
        assert!(one.layout_revision() > one_revision);
        assert_eq!(
            PhaseDataView::single(one.subs_data_view())
                .layout()
                .component_labels(),
            vec!["CO2".to_string(), "H2O".to_string()]
        );

        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string()];
        pos.replace_phase_data(HashMap::from([(Some("gas".to_string()), gas)]));
        pos.debug_replace_dG_sym(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("O2".to_string(), Expr::Const(1.0))]),
        )]));
        let pos_revision = pos.layout_revision();

        pos.with_phase_data_mut(|phase_data| {
            phase_data
                .get_mut(&Some("gas".to_string()))
                .unwrap()
                .substances
                .push("N2".to_string());
            let mut liquid = SubsData::new();
            liquid.substances = vec!["H2O".to_string()];
            phase_data.insert(Some("liquid".to_string()), liquid);
        });

        assert!(pos.dG_sym_view().is_empty());
        assert!(pos.resolved_layout().is_none());
        assert!(pos.layout_revision() > pos_revision);
        assert_eq!(
            pos.result_snapshot(None, None).component_labels(),
            vec![
                "gas::O2".to_string(),
                "gas::N2".to_string(),
                "liquid::H2O".to_string(),
            ]
        );
    }

    #[test]
    fn test_narrow_phase_traits_share_layout_and_pure_evaluation_contracts() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.set_substances(vec!["CO2".to_string()]);
        let one_layout = PhaseLayoutAccess::system_layout(&one);
        assert_eq!(one_layout.component_labels(), vec!["CO2".to_string()]);

        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["CO2".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        let pos_layout = PhaseLayoutAccess::system_layout(&pos);
        assert_eq!(pos_layout.component_labels(), vec!["gas::CO2".to_string()]);

        let conditions = ThermoEvaluationConditions::new(300.0, 101_325.0).unwrap();
        assert!(
            PhasePropertyEvaluator::evaluate_gibbs_at(&one, conditions, &HashMap::new())
                .unwrap_err()
                .to_string()
                .contains("must exactly match")
        );
        assert!(
            PhasePropertyEvaluator::evaluate_entropy_at(&pos, conditions, &HashMap::new())
                .unwrap_err()
                .to_string()
                .contains("must exactly match")
        );
    }

    #[test]
    fn test_one_and_many_phase_facades_implement_the_same_narrow_roles() {
        fn assert_phase_roles<T>()
        where
            T: PhaseDataPreparation + PhaseSymbolicPropertyBuilder + PhaseEquilibriumAssembly,
        {
        }

        assert_phase_roles::<crate::Thermodynamics::User_PhaseOrSolution2::OnePhase>();
        assert_phase_roles::<PhaseOrSolution>();
        assert_phase_roles::<CustomSubstance>();
    }

    #[test]
    fn test_phase_or_solution_revision_bumps_on_layout_changes_but_not_on_read_only_lagrange_build()
    {
        let mut pos = PhaseOrSolution::new();
        let mut phase = SubsData::new();
        phase.substances = vec!["A".to_string(), "B".to_string()];
        pos.insert_phase_data("gas".to_string(), phase);

        let start_revision = pos.layout_revision();
        let _ = pos.indexed_moles_variables().unwrap();
        assert!(pos.layout_revision() > start_revision);

        let phase_g = HashMap::from([
            ("A".to_string(), Expr::Const(1.0)),
            ("B".to_string(), Expr::Const(2.0)),
        ]);
        pos.debug_replace_dG_sym(HashMap::from([(Some("gas".to_string()), phase_g)]));
        let revision_after_layout = pos.layout_revision();
        let eqs = pos
            .calculate_Lagrange_equations_sym(DMatrix::from_vec(2, 1, vec![1.0, 1.0]), 300.0)
            .unwrap();
        assert_eq!(eqs.len(), 2);
        assert_eq!(pos.layout_revision(), revision_after_layout);
    }

    #[test]
    fn test_one_phase_revision_bumps_on_layout_changes_but_not_on_read_only_lagrange_build() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();

        let start_revision = one.layout_revision();
        one.configure_system_properties(101325.0, None, HashMap::new(), None)
            .unwrap();
        assert!(one.layout_revision() > start_revision);

        let revision_after_layout = one.layout_revision();
        one.set_substances(vec!["A".to_string(), "B".to_string()]);
        let revision_after_payload_change = one.layout_revision();
        assert!(revision_after_payload_change > revision_after_layout);
        one.debug_dG_sym_mut()
            .insert("A".to_string(), Expr::Const(1.0));
        one.debug_dG_sym_mut()
            .insert("B".to_string(), Expr::Const(2.0));
        let g_fun = HashMap::from([
            (
                "A".to_string(),
                Box::new(|_, _, _| 1.0) as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
            ),
            (
                "B".to_string(),
                Box::new(|_, _, _| 2.0) as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
            ),
        ]);
        let eqs = one
            .calculate_Lagrange_equations_fun(DMatrix::from_vec(2, 1, vec![1.0, 1.0]), g_fun, 300.0)
            .unwrap();
        let sample = eqs(300.0, None, None, vec![1.0]);
        assert_eq!(sample.len(), 2);
        assert_eq!(one.layout_revision(), revision_after_payload_change);
    }

    #[test]
    fn test_custom_substance_typed_snapshots_track_layout_revision() {
        let mut one = crate::Thermodynamics::User_PhaseOrSolution2::OnePhase::new();
        one.debug_dG_mut().insert("A".to_string(), 1.5);
        let one_custom = CustomSubstance::OnePhase(one.clone());
        let one_snapshot = one_custom.state_snapshot();
        assert_eq!(one_snapshot.layout_revision, one_custom.layout_revision());
        match one_snapshot.dG.values {
            crate::Thermodynamics::User_PhaseOrSolution::NestedPhaseCacheView::Single(view) => {
                assert_eq!(view.get("A"), Some(&1.5));
            }
            _ => panic!("expected single-phase view"),
        }

        let mut pos = PhaseOrSolution::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["A".to_string()];
        pos.insert_phase_data("gas".to_string(), gas);
        pos.debug_replace_dS(HashMap::from([(
            Some("gas".to_string()),
            HashMap::from([("A".to_string(), -2.0)]),
        )]));
        let pos_custom = CustomSubstance::PhaseOrSolution(pos.clone());
        let pos_snapshot = pos_custom.state_snapshot();
        assert_eq!(pos_snapshot.layout_revision, pos_custom.layout_revision());
        match pos_snapshot.dS.values {
            crate::Thermodynamics::User_PhaseOrSolution::NestedPhaseCacheView::Multi(view) => {
                assert_eq!(
                    view.get(&Some("gas".to_string()))
                        .and_then(|inner| inner.get("A"))
                        .copied(),
                    Some(-2.0)
                );
            }
            _ => panic!("expected multiphase view"),
        }
    }

    #[test]
    fn test_calculate_lagrange_equations_sym_and_fun() {
        let mut pos = PhaseOrSolution::new();
        let mut s = SubsData::new();
        s.substances = vec!["A".to_string(), "B".to_string()];
        pos.insert_phase_data("mix".to_string(), s);

        // set simple dG_sym: A=1.0, B=2.0
        let mut inner = HashMap::new();
        inner.insert("A".to_string(), Expr::Const(1.0));
        inner.insert("B".to_string(), Expr::Const(2.0));
        let mut outer = HashMap::new();
        outer.insert(Some("mix".to_string()), inner);
        pos.debug_replace_dG_sym(outer);

        // A should be nrows = #substances x #elements (we construct original A as (#substances x #elements))
        // For 2 substances and 1 element, create A as 2x1
        let A = DMatrix::from_vec(2, 1, vec![1.0, 1.0]);

        let eqs = pos
            .calculate_Lagrange_equations_sym(A.clone(), 300.0)
            .unwrap();
        assert_eq!(eqs.len(), 2);

        // Now test functional version
        // build G_fun with simple closures returning constants
        let mut gfun_phase = HashMap::new();
        gfun_phase.insert(
            "A".to_string(),
            Box::new(|_T: f64, _n: Option<Vec<f64>>, _Np: Option<f64>| 1.0f64)
                as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        gfun_phase.insert(
            "B".to_string(),
            Box::new(|_T: f64, _n: Option<Vec<f64>>, _Np: Option<f64>| 2.0f64)
                as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        let mut G_fun = HashMap::new();
        G_fun.insert(Some("mix".to_string()), gfun_phase);

        let f = pos
            .calculate_Lagrange_equations_fun(A.clone(), G_fun, 300.0)
            .unwrap();

        // Call with arbitrary values: T, n, Np, Lambda (Lambda length equals number of elements = 1)
        let vals = f(300.0, None, None, vec![0.5]);
        assert_eq!(vals.len(), 2);
        // For our simple setup, eq_i = sum_by_elements + G_i/(R*Tm)
        // sum_by_elements = 0.5 (since Lambda=0.5 and element count 1)
        // so eq for A should be 0.5 + 1/(R*300)
        let expected_a = 0.5 + 1.0f64 / (crate::Thermodynamics::User_PhaseOrSolution::R * 300.0);
        assert_relative_eq!(vals[0], expected_a, epsilon = 1e-12);
    }

    #[test]
    fn test_set_T_and_P_in_G_sym() {
        let mut pos = PhaseOrSolution::new();
        let mut s = SubsData::new();
        s.substances = vec!["A".to_string()];
        pos.insert_phase_data("p".to_string(), s);

        let mut inner = HashMap::new();
        // construct simple symbolic expression including T and P
        inner.insert(
            "A".to_string(),
            Expr::Var("T".to_string()) + Expr::Var("P".to_string()),
        );
        let mut outer = HashMap::new();
        outer.insert(Some("p".to_string()), inner);
        pos.debug_replace_dG_sym(outer);

        // set P and T
        pos.set_P_to_sym_in_G_sym(101325.0);
        pos.set_T_to_sym_in_G_sym(298.15);

        // now dG_sym should contain constants (simplified)
        let dG_sym_view = pos.dG_sym_view();
        let v = dG_sym_view.get(&Some("p".to_string())).unwrap();
        let expr = v.get("A").unwrap();
        // Expect expression to be simplified to a constant: 298.15 + 101325.0
        let expected = Expr::Const(298.15 + 101325.0);
        assert_eq!(
            format!("{}", expr.simplify()),
            format!("{}", expected.simplify())
        );
    }

    #[test]
    fn test_substance_system_factory_multi_phase() {
        let mut phase_substances = HashMap::new();
        phase_substances.insert("gas".to_string(), vec!["N2".to_string(), "O2".to_string()]);
        phase_substances.insert("gas2".to_string(), vec!["H2O".to_string()]);

        let container = SubstancesContainer::MultiPhase(phase_substances);

        let result = SubstanceSystemFactory::create_system(
            container,
            Some(HashMap::from([
                (
                    "gas".to_string(),
                    crate::Thermodynamics::User_substances::Phases::Gas,
                ),
                (
                    "gas2".to_string(),
                    crate::Thermodynamics::User_substances::Phases::Gas,
                ),
            ])),
            vec!["NASA_gas".to_string(), "NASA_cond".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );

        assert!(result.is_ok());

        // Test that we got a PhaseOrSolution
        let mut pos = match result.unwrap() {
            CustomSubstance::PhaseOrSolution(pos) => pos,
            _ => panic!("Expected PhaseOrSolution"),
        };

        // Test extract_all_thermal_coeffs
        let thermal_result = pos.extract_all_thermal_coeffs(400.0);
        assert!(thermal_result.is_ok());

        // Test calculate_therm_map_of_properties
        let therm_props_result = pos.calculate_therm_map_of_properties(400.0);
        assert!(therm_props_result.is_ok());

        // Test calculate_therm_map_of_sym
        let therm_sym_result = pos.calculate_therm_map_of_sym();
        assert!(therm_sym_result.is_ok());

        // Test get_all_substances
        let all_subs = pos.get_all_substances();
        assert_eq!(all_subs.len(), 3);
        assert!(all_subs.contains(&"N2".to_string()));
        assert!(all_subs.contains(&"O2".to_string()));
        assert!(all_subs.contains(&"H2O".to_string()));

        // Test extract_SubstancesContainer
        let container_result = pos.extract_SubstancesContainer();
        assert!(container_result.is_ok());
        match container_result.unwrap() {
            SubstancesContainer::MultiPhase(phases) => {
                assert!(phases.contains_key("gas"));
                assert!(phases.contains_key("gas2"));
                assert_eq!(phases["gas"].len(), 2);
                assert_eq!(phases["gas2"].len(), 1);
            }
            _ => panic!("Expected MultiPhase container"),
        }

        // Test calculate_elem_composition_and_molar_mass
        let elem_comp_result = pos.calculate_elem_composition_and_molar_mass(None);
        assert!(elem_comp_result.is_ok());
        let (matrix, molar_masses, elements) = elem_comp_result.unwrap();
        assert_eq!(matrix.nrows(), 3); // 3 substances
        assert!(matrix.ncols() >= 2); // At least N, O, H elements
        assert!(molar_masses.contains_key("N2"));
        assert!(molar_masses.contains_key("O2"));
        assert!(molar_masses.contains_key("H2O"));
        assert!(elements.contains(&"N".to_string()));
        assert!(elements.contains(&"O".to_string()));
        assert!(elements.contains(&"H".to_string()));

        // Test indexed_moles_variables
        let indexed_vars_result = pos.indexed_moles_variables();
        assert!(indexed_vars_result.is_ok());
        let (map_indexed, vec_of_n_vars, np_vec, map_of_var_each) = indexed_vars_result.unwrap();
        assert_eq!(vec_of_n_vars.len(), 3); // 3 substances total
        assert_eq!(np_vec.len(), 2); // 2 phases
        assert!(map_indexed.contains_key(&Some("gas".to_string())));
        assert!(map_indexed.contains_key(&Some("gas2".to_string())));
        assert!(map_of_var_each.contains_key(&Some("gas".to_string())));
        assert!(map_of_var_each.contains_key(&Some("gas2".to_string())));

        // Test calculate_Gibbs_sym
        let gibbs_sym_result = pos.calculate_Gibbs_sym(400.0);
        assert!(gibbs_sym_result.is_ok());
        assert!(!pos.dG_sym_view().is_empty());

        // Test calculate_Gibbs_fun
        assert!(pos.calculate_Gibbs_fun(400.0, 101325.0).is_ok());
        assert!(!pos.dG_fun_view().is_empty());

        // Test calculate_S_sym
        let entropy_sym_result = pos.calculate_S_sym(400.0);
        assert!(entropy_sym_result.is_ok());
        assert!(!pos.dS_sym_view().is_empty());

        // Test calculate_S_fun
        assert!(pos.calculate_S_fun(400.0, 101325.0).is_ok());
        assert!(!pos.dS_fun_view().is_empty());

        // Test calcutate_Gibbs_free_energy with mole numbers
        let mut n = HashMap::new();
        let gas_moles = vec![0.4, 0.2];
        let gas2_moles = vec![0.4];
        n.insert(
            Some("gas".to_string()),
            PhaseComposition::new(Some(0.6), Some(gas_moles)),
        );
        n.insert(
            Some("gas2".to_string()),
            PhaseComposition::new(Some(0.4), Some(gas2_moles)),
        );

        let gibbs_energy_result = pos.calcutate_Gibbs_free_energy(400.0, 101325.0, n.clone());
        assert!(gibbs_energy_result.is_ok());
        let gibbs_map = gibbs_energy_result.unwrap();
        assert!(gibbs_map.contains_key(&Some("gas".to_string())));
        assert!(gibbs_map.contains_key(&Some("gas2".to_string())));

        // Test calculate_S with mole numbers
        let entropy_result = pos.calculate_S(400.0, 101325.0, n);
        assert!(entropy_result.is_ok());
        assert!(!pos.dS_view().is_empty());

        // Invalid numeric phase composition must fail before any calculation runs.
        let invalid_gibbs = pos.calcutate_Gibbs_free_energy(
            400.0,
            101325.0,
            HashMap::from([
                (
                    Some("gas".to_string()),
                    PhaseComposition::new(Some(0.5), Some(vec![0.7, -0.2])),
                ),
                (
                    Some("gas2".to_string()),
                    PhaseComposition::new(Some(0.4), Some(vec![0.4])),
                ),
            ]),
        );
        let err = invalid_gibbs.unwrap_err();
        assert!(format!("{}", err).contains("invalid component amount"));

        let invalid_entropy = pos.calculate_S(
            400.0,
            101325.0,
            HashMap::from([
                (
                    Some("gas".to_string()),
                    PhaseComposition::new(Some(0.6), Some(vec![0.4, 0.2])),
                ),
                (
                    Some("gas2".to_string()),
                    PhaseComposition::new(Some(1.0), Some(vec![0.3, 0.4])),
                ),
            ]),
        );
        let err = invalid_entropy.unwrap_err();
        assert!(format!("{}", err).contains("does not match component sum"));

        // Test configure_system_properties
        let molar_masses = HashMap::from([
            ("N2".to_string(), 28.014),
            ("O2".to_string(), 31.998),
            ("H2O".to_string(), 18.015),
        ]);
        let config_result = pos.configure_system_properties(
            101325.0,
            Some("Pa".to_string()),
            molar_masses,
            Some("g/mol".to_string()),
        );
        assert!(config_result.is_ok());

        // Test if_not_found_go_NIST
        let nist_result = pos.if_not_found_go_NIST();
        assert!(nist_result.is_ok());

        // Test extract_coeffs_if_current_coeffs_not_valid
        let coeffs_result = pos.extract_coeffs_if_current_coeffs_not_valid(400.0);
        assert!(coeffs_result.is_ok());

        // Payload/configuration changes invalidate every derived cache. A
        // caller must explicitly rebuild expressions before substituting into
        // them; retaining the old symbolic Gibbs maps would mix coefficient
        // revisions in one solver state.
        assert_eq!(pos.dG_view().populated_phase_count(), 0);
        assert_eq!(pos.dG_fun_view().populated_phase_count(), 0);
        assert_eq!(pos.dG_sym_view().populated_phase_count(), 0);
        assert_eq!(pos.dS_view().populated_phase_count(), 0);
        assert_eq!(pos.dS_fun_view().populated_phase_count(), 0);
        assert_eq!(pos.dS_sym_view().populated_phase_count(), 0);

        // Test create_full_map_of_mole_numbers
        let mut partial_n = HashMap::new();
        let partial_gas_moles = HashMap::from([
            ("N2".to_string(), 0.5),
            // O2 missing - should be filled with 0.0
        ]);
        partial_n.insert(
            Some("gas".to_string()),
            (Some(0.6), Some(partial_gas_moles)),
        );
        partial_n.insert(
            Some("gas2".to_string()),
            (Some(0.4), Some(HashMap::from([("H2O".to_string(), 0.4)]))),
        );

        let full_map_result = pos.create_full_map_of_mole_numbers(partial_n);
        assert!(full_map_result.is_ok());
        let (full_map, _vec_map, summed) = full_map_result.unwrap();

        // Check that missing substances were filled
        let gas_data = full_map.get(&Some("gas".to_string())).unwrap();
        let gas_moles = gas_data.1.as_ref().unwrap();
        assert!(gas_moles.contains_key("O2"));
        assert_relative_eq!(*gas_moles.get("O2").unwrap(), 0.0, epsilon = 1e-10);

        // Check summed map
        assert!(summed.contains_key("N2"));
        assert!(summed.contains_key("O2"));
        assert!(summed.contains_key("H2O"));
        assert_relative_eq!(*summed.get("N2").unwrap(), 0.5, epsilon = 1e-10);
        assert_relative_eq!(*summed.get("O2").unwrap(), 0.0, epsilon = 1e-10);
        assert_relative_eq!(*summed.get("H2O").unwrap(), 0.4, epsilon = 1e-10);

        // Test set_P_to_sym_in_G_sym and set_T_to_sym_in_G_sym
        assert!(pos.calculate_Gibbs_sym(400.0).is_ok());
        pos.set_P_to_sym_in_G_sym(200000.0);
        pos.set_T_to_sym_in_G_sym(350.0);
        // Verify that symbolic expressions have been updated
        assert!(!pos.dG_sym_view().is_empty());
        let dG_sym_view = pos.dG_sym_view();
        let gas_sym = dG_sym_view.get(&Some("gas".to_string())).unwrap();
        let n2_expr = gas_sym.get("N2").unwrap();
        // The expression should now contain constants instead of variables
        assert!(!format!("{}", n2_expr).contains("T"));
        assert!(!format!("{}", n2_expr).contains("P"));

        // Test calculate_Lagrange_equations_sym
        // First need to set up symbolic variables and Gibbs expressions
        let _ = pos.indexed_moles_variables();
        let _ = pos.calculate_Gibbs_sym(400.0);
        let A = DMatrix::from_vec(3, 2, vec![2.0, 0.0, 0.0, 2.0, 1.0, 0.0]); // N2, O2, H2O with N, O elements
        let lagrange_sym_result = pos.calculate_Lagrange_equations_sym(A.clone(), 400.0);
        assert!(lagrange_sym_result.is_ok());
        let lagrange_eqs = lagrange_sym_result.unwrap();
        assert_eq!(lagrange_eqs.len(), 3); // 3 substances

        // Test calculate_Lagrange_equations_fun
        /*
                let lagrange_fun_result = pos.calculate_Lagrange_equations_fun(
                    A.clone(),
                    pos.dG_fun_view().to_owned_map(),
                    400.0,
                );
                assert!(lagrange_fun_result.is_ok());
                let lagrange_fun = lagrange_fun_result.unwrap();
        */
        // Test the function with sample values
        //   let test_lambda = vec![1.0, 2.0]; // Lambda values for N and O elements
        //  let test_vals = lagrange_fun(400.0, None, None, test_lambda);
        // assert_eq!(test_vals.len(), 3); // 3 substances

        // Test calculate_Lagrange_equations_fun2
        let lagrange_fun2_result = pos.calculate_Lagrange_equations_fun2(A, 400.0, 101325.0, 400.0);
        assert!(lagrange_fun2_result.is_ok());
        let lagrange_fun2 = lagrange_fun2_result.unwrap();

        // Test the function with sample values
        let test_vals2 = lagrange_fun2(400.0, None, None, vec![0.5, 1.5]);
        assert_eq!(test_vals2.len(), 3); // 3 substances

        let revision_before_failed_assembly = pos.layout_revision();
        let cache_shape_before_failed_assembly = pos
            .dG_fun_view()
            .iter()
            .map(|(phase, functions)| {
                let mut substances = functions.keys().cloned().collect::<Vec<_>>();
                substances.sort();
                (phase, substances)
            })
            .collect::<HashMap<_, _>>();
        let invalid_matrix_error = match pos.calculate_Lagrange_equations_fun2(
            DMatrix::zeros(2, 2),
            400.0,
            101_325.0,
            400.0,
        ) {
            Ok(_) => panic!("a matrix with too few component rows must be rejected"),
            Err(error) => error,
        };

        assert!(invalid_matrix_error.to_string().contains("component rows"));
        assert_eq!(pos.layout_revision(), revision_before_failed_assembly);
        let cache_shape_after_failed_assembly = pos
            .dG_fun_view()
            .iter()
            .map(|(phase, functions)| {
                let mut substances = functions.keys().cloned().collect::<Vec<_>>();
                substances.sort();
                (phase, substances)
            })
            .collect::<HashMap<_, _>>();
        assert_eq!(
            cache_shape_after_failed_assembly,
            cache_shape_before_failed_assembly
        );
    }
}
