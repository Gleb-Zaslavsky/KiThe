//! Contract tests for the resolved-phase-system equilibrium boundary.
//!
//! These tests exercise the full one-way bridge from resolved phase data to a
//! canonical fixed-set equilibrium solve. They prove deterministic qualified
//! identity, local standard-state extraction, provenance retention, immutable
//! resolution input, and transactional accepted-solution publication.

use std::collections::HashMap;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    equilibrium_logmole_jacobian, equilibrium_logmole_residual,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, PreparedEquilibriumProblem, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RustedSciTheSolver, prepare_rst_symbolic_problem_from_prepared,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::multiphase_equilibrium_residual_generator_sym;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    PhaseEquilibriumBuildRequest, PhaseEquilibriumMetadata, PhaseEquilibriumProblemBundle,
    SupportedPhaseModelPolicy, build_phase_equilibrium_problem,
    build_phase_equilibrium_problem_with_timing,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::User_PhaseOrSolution::{PhaseModel, PhaseSpec, ResolvedPhaseSystem};
use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use RustedSciThe::symbolic::symbolic_engine::Expr;

fn gas_spec() -> PhaseSpec {
    PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2O".to_string(), "O2".to_string()],
    )
    .unwrap()
}

fn liquid_spec() -> PhaseSpec {
    PhaseSpec::pure_condensed(
        PhaseId::new(Some("liquid".to_string())),
        vec!["H2O".to_string()],
        PhysicalState::Liquid,
    )
    .unwrap()
}

fn phase_data(substances: &[&str]) -> SubsData {
    let mut data = SubsData::new();
    data.substances = substances
        .iter()
        .map(|value| (*value).to_string())
        .collect();
    data
}

fn resolved_local_nasa_gas() -> ResolvedPhaseSystem {
    let spec = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()],
    )
    .unwrap();
    let mut data = phase_data(&["H2", "O2", "H2O"]);
    data.set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    data.search_substances().unwrap();
    data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(vec![spec], HashMap::from([(Some("gas".to_string()), data)])).unwrap()
}

/// Real local thermochemistry fixture for the semantic ideal-solution bridge.
///
/// `AL(cr)` and `AL2O3(a)` are deliberately used only to prove that a
/// multicomponent condensed phase preserves every selected local record
/// through the production bridge. This is a boundary-contract fixture, not a
/// claim that the pair forms a physically calibrated solution model.
fn resolved_local_nasa_condensed_ideal_solution() -> ResolvedPhaseSystem {
    let spec = PhaseSpec::ideal_solution(
        PhaseId::new(Some("oxide_solution".to_string())),
        vec!["AL(cr)".to_string(), "AL2O3(a)".to_string()],
        PhysicalState::Solid,
    )
    .unwrap();
    let mut data = phase_data(&["AL(cr)", "AL2O3(a)"]);
    data.set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    data.set_substance_physical_state("AL(cr)".to_string(), PhysicalState::Solid);
    data.set_substance_physical_state("AL2O3(a)".to_string(), PhysicalState::Solid);
    data.search_substances().unwrap();
    data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(
        vec![spec],
        HashMap::from([(Some("oxide_solution".to_string()), data)]),
    )
    .unwrap()
}

fn prepared_local_nasa_gas() -> PhaseEquilibriumProblemBundle {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9]).unwrap();
    build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn bridge_timing_total_encloses_every_recorded_setup_stage() {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9]).unwrap();
    let bundle = build_phase_equilibrium_problem_with_timing(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
        EquilibriumTimingMode::Enabled,
    )
    .unwrap();
    let timing = bundle.timing_report();

    assert!(timing.enabled());
    assert!(timing.total() >= timing.thermochemistry_preparation());
    assert!(timing.total() >= timing.numeric_closure_construction());
    assert!(timing.total() >= timing.symbolic_construction());
    assert!(timing.total() >= timing.equation_construction());
}

fn resolved_local_gas_and_condensed_water() -> ResolvedPhaseSystem {
    resolved_local_gas_and_condensed_water_with_map_order(false)
}

/// Local fixture whose gas phase spans both H and O conservation directions.
/// This makes it suitable for a bounded active-set solve while still retaining
/// H2O as a zero-inventory pure-liquid candidate.
fn resolved_local_gas_oxygen_and_condensed_water() -> ResolvedPhaseSystem {
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2O".to_string(), "O2".to_string()],
    )
    .unwrap();
    let liquid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("liquid".to_string())),
        vec!["H2O".to_string()],
        PhysicalState::Liquid,
    )
    .unwrap();

    let mut gas_data = phase_data(&["H2O", "O2"]);
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.set_substance_physical_state("H2O".to_string(), PhysicalState::Gas);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    let mut condensed_data = phase_data(&["H2O"]);
    condensed_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    condensed_data.set_substance_physical_state("H2O".to_string(), PhysicalState::Liquid);
    condensed_data.search_substances().unwrap();
    condensed_data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(
        vec![gas, liquid],
        HashMap::from([
            (Some("gas".to_string()), gas_data),
            (Some("liquid".to_string()), condensed_data),
        ]),
    )
    .unwrap()
}

fn resolved_local_gas_and_condensed_water_with_map_order(
    reverse_map_order: bool,
) -> ResolvedPhaseSystem {
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2O".to_string()],
    )
    .unwrap();
    let liquid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("liquid".to_string())),
        vec!["H2O".to_string()],
        PhysicalState::Liquid,
    )
    .unwrap();

    let mut gas_data = phase_data(&["H2O"]);
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.set_substance_physical_state("H2O".to_string(), PhysicalState::Gas);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    let mut condensed_data = phase_data(&["H2O"]);
    condensed_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    condensed_data.set_substance_physical_state("H2O".to_string(), PhysicalState::Liquid);
    condensed_data.search_substances().unwrap();
    condensed_data.parse_all_thermal_coeffs().unwrap();

    let entries = [
        (Some("gas".to_string()), gas_data),
        (Some("liquid".to_string()), condensed_data),
    ];
    let mut data = HashMap::new();
    if reverse_map_order {
        for (phase, payload) in entries.into_iter().rev() {
            data.insert(phase, payload);
        }
    } else {
        for (phase, payload) in entries {
            data.insert(phase, payload);
        }
    }

    ResolvedPhaseSystem::new(vec![gas, liquid], data).unwrap()
}

fn resolved_with_missing_last_phase() -> ResolvedPhaseSystem {
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2".to_string()],
    )
    .unwrap();
    let solid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("solid".to_string())),
        vec!["MissingSpecies".to_string()],
        PhysicalState::Solid,
    )
    .unwrap();
    let mut gas_data = phase_data(&["H2"]);
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(
        vec![gas, solid],
        HashMap::from([
            (Some("gas".to_string()), gas_data),
            (Some("solid".to_string()), phase_data(&["MissingSpecies"])),
        ]),
    )
    .unwrap()
}

fn resolved_system(specs: Vec<PhaseSpec>, reverse_map_order: bool) -> ResolvedPhaseSystem {
    let mut data = HashMap::new();
    let entries = [
        (Some("gas".to_string()), phase_data(&["H2O", "O2"])),
        (Some("liquid".to_string()), phase_data(&["H2O"])),
    ];
    if reverse_map_order {
        for (phase, payload) in entries.into_iter().rev() {
            data.insert(phase, payload);
        }
    } else {
        for (phase, payload) in entries {
            data.insert(phase, payload);
        }
    }
    ResolvedPhaseSystem::new(specs, data).unwrap()
}

#[test]
fn metadata_retains_qualified_identity_and_maps_activity_models() {
    let resolved = resolved_system(vec![liquid_spec(), gas_spec()], true);
    let metadata = PhaseEquilibriumMetadata::from_resolved(
        &resolved,
        SupportedPhaseModelPolicy::IdealPhaseModelsV1,
    )
    .unwrap();

    let labels = metadata
        .components()
        .iter()
        .map(|component| component.label())
        .collect::<Vec<_>>();
    assert_eq!(labels, ["gas::H2O", "gas::O2", "liquid::H2O"]);
    assert_eq!(metadata.components()[0].substance(), "H2O");
    assert_eq!(
        metadata.components()[0].activity_model(),
        PhaseActivityModel::IdealGas
    );
    assert_eq!(
        metadata.components()[2].activity_model(),
        PhaseActivityModel::IdealSolution
    );
    assert_eq!(
        metadata.components()[2].phase_model(),
        PhaseModel::PureCondensed
    );
    assert_eq!(
        metadata.components()[2].physical_state(),
        PhysicalState::Liquid
    );

    let gas_h2o = PhaseComponentId::new(PhaseId::new(Some("gas".to_string())), "H2O");
    let liquid_h2o = PhaseComponentId::new(PhaseId::new(Some("liquid".to_string())), "H2O");
    assert_eq!(metadata.component_index(&gas_h2o), Some(0));
    assert_eq!(metadata.component_index(&liquid_h2o), Some(2));
}

#[test]
fn production_bridge_preserves_all_multicomponent_ideal_solution_data() {
    let resolved = resolved_local_nasa_condensed_ideal_solution();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![0.75, 0.25]).unwrap();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(900.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            SupportedPhaseModelPolicy::IdealPhaseModelsV1,
        )
        .unwrap(),
    )
    .unwrap();
    let metadata = bundle.metadata();

    assert_eq!(metadata.components().len(), 2);
    assert_eq!(metadata.phases().len(), 1);
    assert_eq!(metadata.phases()[0].component_range(), 0..2);
    assert_eq!(
        metadata.phases()[0].phase_model(),
        PhaseModel::IdealSolution
    );
    assert_eq!(
        metadata.phases()[0].activity_model(),
        PhaseActivityModel::IdealSolution
    );
    assert!(
        metadata
            .components()
            .iter()
            .all(|component| component.phase_model() == PhaseModel::IdealSolution)
    );
    assert!(
        metadata
            .components()
            .iter()
            .all(|component| component.activity_model() == PhaseActivityModel::IdealSolution)
    );
    assert_eq!(bundle.problem().gibbs().len(), 2);
    assert_eq!(bundle.problem().element_composition().nrows(), 2);
    assert!(bundle.problem().element_composition().ncols() >= 2);
    assert_eq!(bundle.report().components().len(), 2);
}

#[test]
fn hash_map_and_phase_declaration_order_do_not_change_bridge_metadata() {
    let first = resolved_system(vec![gas_spec(), liquid_spec()], false);
    let second = resolved_system(vec![liquid_spec(), gas_spec()], true);

    let first_metadata =
        PhaseEquilibriumMetadata::from_resolved(&first, Default::default()).unwrap();
    let second_metadata =
        PhaseEquilibriumMetadata::from_resolved(&second, Default::default()).unwrap();

    assert_eq!(first_metadata, second_metadata);
    assert_eq!(
        first_metadata
            .phases()
            .iter()
            .map(|phase| phase.index().index())
            .collect::<Vec<_>>(),
        [0, 1]
    );
    assert_eq!(first_metadata.phases()[0].component_range(), 0..2);
    assert_eq!(first_metadata.phases()[1].component_range(), 2..3);
}

#[test]
fn metadata_retains_the_resolution_report_in_canonical_phase_order() {
    let resolved = resolved_system(vec![liquid_spec(), gas_spec()], true);
    let metadata = PhaseEquilibriumMetadata::from_resolved(&resolved, Default::default()).unwrap();

    assert_eq!(metadata.provenance(), resolved.report());
    assert_eq!(metadata.provenance().phases().len(), 2);
    assert_eq!(
        metadata.provenance().phases()[0].phase(),
        &PhaseId::new(Some("gas".to_string()))
    );
    assert_eq!(
        metadata.provenance().phases()[1].phase(),
        &PhaseId::new(Some("liquid".to_string()))
    );
}

#[test]
fn build_request_rejects_composition_from_a_different_layout() {
    let resolved = resolved_system(vec![gas_spec(), liquid_spec()], false);
    let foreign_layout = MultiphaseEquilibriumLayout::new(vec![
        PhaseSpec::ideal_gas(
            PhaseId::new(Some("other".to_string())),
            vec!["H2O".to_string(), "O2".to_string()],
        )
        .unwrap(),
    ])
    .unwrap();
    let foreign_composition =
        MultiphaseInitialComposition::from_dense(&foreign_layout, vec![1.0, 1.0]).unwrap();
    let conditions = EquilibriumConditions::new(1000.0, 101_325.0, 101_325.0).unwrap();

    let error = PhaseEquilibriumBuildRequest::new(
        &resolved,
        conditions,
        foreign_composition,
        TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
        Default::default(),
    )
    .unwrap_err();

    assert!(
        error
            .to_string()
            .contains("composition belongs to a different multiphase layout")
    );
}

#[test]
fn build_request_keeps_physical_zeroes_separate_from_trace_seed_policy() {
    let resolved = resolved_system(vec![gas_spec(), liquid_spec()], false);
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();
    let conditions = EquilibriumConditions::new(1200.0, 202_650.0, 101_325.0).unwrap();
    let trace_policy = TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
        fraction: 1e-20,
        minimum_floor: 1e-30,
    };

    let request = PhaseEquilibriumBuildRequest::new(
        &resolved,
        conditions,
        composition,
        trace_policy,
        Default::default(),
    )
    .unwrap();

    assert_eq!(request.initial_composition().moles(), [2.0, 1.0, 0.0]);
    assert_eq!(request.trace_seed_policy(), trace_policy);
    assert_eq!(request.conditions(), conditions);
    assert_eq!(
        request.metadata().layout_fingerprint(),
        layout.fingerprint()
    );
}

#[test]
fn local_nasa_gas_builds_a_complete_problem_and_retains_provenance() {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();
    let conditions = EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            conditions,
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    assert_eq!(
        bundle.problem().species(),
        ["gas::H2", "gas::O2", "gas::H2O"]
    );
    assert_eq!(bundle.problem().initial_moles(), [2.0, 1.0, 0.0]);
    assert_eq!(bundle.problem().phases().len(), 1);
    assert!(bundle.problem().initial_log_moles().as_slice()[2].is_finite());
    assert!((bundle.problem().initial_log_moles().as_slice()[2] - 1e-30_f64.ln()).abs() < 1e-12);

    let report = bundle.report();
    assert_eq!(report.conditions(), conditions);
    assert_eq!(report.components().len(), 3);
    assert!(
        report
            .components()
            .iter()
            .all(|row| row.standard_gibbs_at_conditions().is_finite())
    );
    assert!(
        report
            .components()
            .iter()
            .all(|row| row.thermo_source().library() == "NASA_gas")
    );

    let totals = report
        .element_labels()
        .iter()
        .cloned()
        .zip(report.element_totals().iter().copied())
        .collect::<HashMap<_, _>>();
    assert_eq!(totals.get("H"), Some(&4.0));
    assert_eq!(totals.get("O"), Some(&2.0));
}

#[test]
fn local_phase_data_solves_through_one_accepted_bridge_bundle() {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();
    let conditions = EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap();
    let prepared = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            conditions,
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    let accepted = prepared.solve().unwrap();

    assert_eq!(accepted.metadata().layout(), resolved.layout());
    assert_eq!(accepted.build_report().conditions(), conditions);
    assert_eq!(accepted.solution().conditions(), conditions);
    assert_eq!(accepted.solution().moles().len(), 3);
    assert!(accepted.solution().moles().iter().all(|moles| *moles > 0.0));
    assert!(accepted.solution().validation().min_moles > 0.0);
    assert!(
        accepted
            .solution()
            .validation()
            .max_abs_element_balance_error
            .is_finite()
    );
    assert!(accepted.solve_report().accepted_attempt().is_some());
    assert!(matches!(
        accepted.solve_report().accepted_backend,
        SolverBackend::RustedSciThe(_)
    ));
    assert_eq!(
        accepted.build_report().components()[2].component().label(),
        "gas::H2O"
    );
    let totals = accepted
        .build_report()
        .element_labels()
        .iter()
        .cloned()
        .zip(accepted.build_report().element_totals().iter().copied())
        .collect::<HashMap<_, _>>();
    // H2O is physically absent here and receives only a numerical trace seed.
    // The conserved inventory must remain exactly 2 H2 + 1 O2.
    assert_eq!(totals.get("H"), Some(&4.0));
    assert_eq!(totals.get("O"), Some(&2.0));
}

#[test]
fn invalid_bridge_solver_settings_fail_without_mutating_resolved_phase_data() {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();
    let before = resolved
        .phase_data()
        .get(&Some("gas".to_string()))
        .expect("gas payload must exist");
    assert!(before.element_composition_matrix().is_none());
    assert!(before.therm_functions().is_empty());

    let prepared = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let error = prepared
        .solve_with(|settings| settings.solver_params.max_iter = 0)
        .unwrap_err();

    assert!(error.to_string().contains("solver_params.max_iter"));
    let after = resolved
        .phase_data()
        .get(&Some("gas".to_string()))
        .expect("gas payload must still exist");
    assert!(after.element_composition_matrix().is_none());
    assert!(after.therm_functions().is_empty());
}

#[test]
fn every_rst_backend_receives_the_same_phase_bridge_problem() {
    for backend in RustedSciTheSolver::recommended_cascade() {
        let result = prepared_local_nasa_gas().solve_with(|settings| {
            settings.solver_policy =
                Some(SolverPolicy::Single(SolverBackend::RustedSciThe(backend)));
        });

        match result {
            Ok(accepted) => {
                assert_eq!(
                    accepted
                        .metadata()
                        .components()
                        .iter()
                        .map(|component| component.label())
                        .collect::<Vec<_>>(),
                    ["gas::H2", "gas::O2", "gas::H2O"]
                );
                assert_eq!(
                    accepted.solve_report().policy,
                    SolverPolicy::Single(SolverBackend::RustedSciThe(backend))
                );
                assert_eq!(
                    accepted.solve_report().accepted_backend,
                    SolverBackend::RustedSciThe(backend)
                );
                assert_eq!(accepted.solution().moles().len(), 3);
                assert!(accepted.solution().moles().iter().all(|moles| *moles > 0.0));
            }
            Err(ReactionExtentError::AllBackendsFailed { attempts }) => {
                assert_eq!(attempts.len(), 1, "{backend:?}");
                assert_eq!(attempts[0].backend, SolverBackend::RustedSciThe(backend));
                assert!(attempts[0].outcome.is_started(), "{backend:?}");
            }
            Err(error) => panic!(
                "{backend:?} must reach bridge-owned RST preparation instead of failing setup: {error}"
            ),
        }
    }
}

#[test]
fn rst_fallback_reuses_one_prepared_phase_bridge_problem() {
    let accepted = prepared_local_nasa_gas()
        .solve_with(|settings| {
            // This test isolates deterministic backend fallback rather than
            // the stricter production acceptance preset.
            settings.solver_params.tol = 1e-5;
            settings.solver_policy = Some(SolverPolicy::Cascade(vec![
                SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
                SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
            ]));
        })
        .unwrap();

    assert_eq!(
        accepted.solve_report().attempt_backends(),
        [
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ]
    );
    assert!(accepted.solve_report().attempts[0].outcome.is_started());
    assert!(accepted.solve_report().attempts[1].outcome.is_accepted());
    assert_eq!(
        accepted.solve_report().accepted_backend,
        SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
    );
    assert_eq!(
        accepted
            .metadata()
            .components()
            .iter()
            .map(|component| component.label())
            .collect::<Vec<_>>(),
        ["gas::H2", "gas::O2", "gas::H2O"]
    );
}

#[test]
fn failed_bridge_backend_cannot_overwrite_an_accepted_snapshot() {
    let accepted = prepared_local_nasa_gas().solve().unwrap();
    let accepted_moles = accepted.solution().moles().to_vec();
    let accepted_sources = accepted
        .build_report()
        .components()
        .iter()
        .map(|component| component.thermo_source().record_key().to_string())
        .collect::<Vec<_>>();

    let error = prepared_local_nasa_gas()
        .solve_with(|settings| {
            settings.solver_policy = Some(SolverPolicy::Single(SolverBackend::RustedSciThe(
                RustedSciTheSolver::NielsenLevenbergMarquardt,
            )));
        })
        .unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::AllBackendsFailed { .. }
    ));
    assert_eq!(accepted.solution().moles(), accepted_moles);
    assert_eq!(
        accepted
            .build_report()
            .components()
            .iter()
            .map(|component| component.thermo_source().record_key().to_string())
            .collect::<Vec<_>>(),
        accepted_sources
    );
}

#[test]
fn bridge_standard_gibbs_matches_direct_subsdata_and_is_pressure_independent() {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let temperature = 900.0;
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();
    let at_reference_pressure = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
            composition.clone(),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let at_double_pressure = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(temperature, 202_650.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    let mut direct_payload = resolved
        .phase_data()
        .get(&Some("gas".to_string()))
        .expect("gas payload must exist")
        .clone();
    direct_payload
        .extract_all_thermal_coeffs(temperature)
        .unwrap();
    let direct_gibbs = direct_payload.calculate_dG0_fun_one_phase().unwrap();
    let (_, direct_compositions, _) =
        SubsData::calculate_elem_composition_and_molar_mass_local(&mut direct_payload, None)
            .unwrap();
    for (index, component) in at_reference_pressure
        .metadata()
        .components()
        .iter()
        .enumerate()
    {
        let expected = direct_gibbs[component.substance()](temperature);
        let reference_value = at_reference_pressure.problem().gibbs()[index](temperature);
        let double_pressure_value = at_double_pressure.problem().gibbs()[index](temperature);
        assert!((reference_value - expected).abs() < 1e-10);
        assert!((double_pressure_value - expected).abs() < 1e-10);
        for (element_index, element) in at_reference_pressure
            .report()
            .element_labels()
            .iter()
            .enumerate()
        {
            assert_eq!(
                at_reference_pressure.problem().element_composition()[(index, element_index)],
                direct_compositions[index]
                    .get(element)
                    .copied()
                    .unwrap_or(0.0)
            );
        }
    }
}

#[test]
fn same_molecule_in_gas_and_condensed_records_keeps_one_formula_and_two_components() {
    let resolved = resolved_local_gas_and_condensed_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0, 0.0]).unwrap();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    assert_eq!(bundle.problem().species(), ["gas::H2O", "liquid::H2O"]);
    assert_eq!(
        bundle.problem().element_composition().row(0),
        bundle.problem().element_composition().row(1)
    );
    assert_eq!(
        bundle.report().components()[0].thermo_source().library(),
        "NASA_gas"
    );
    assert_eq!(
        bundle.report().components()[1].thermo_source().library(),
        "NASA_cond"
    );
}

#[test]
fn direct_pure_water_bridge_reduces_dependent_h_o_constraints_deterministically() {
    let resolved = resolved_local_gas_and_condensed_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0, 0.0]).unwrap();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(300.0, 3_536.717_586_505, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    // The report retains full physical provenance, while the solver receives
    // one independent conservation direction for two copies of H2O.
    assert_eq!(bundle.report().element_labels(), ["H", "O"]);
    assert_eq!(bundle.report().element_totals(), [2.0, 1.0]);
    assert_eq!(bundle.report().solver_element_labels(), ["H"]);
    assert_eq!(bundle.problem().element_composition().ncols(), 1);

    let prepared = PreparedEquilibriumProblem::new(bundle.problem().clone()).unwrap();
    assert_eq!(prepared.reaction_basis().num_reactions, 1);
    assert_eq!(prepared.element_totals(), [2.0]);
}

#[test]
fn direct_pure_water_fixed_pt_solve_needs_no_inert_carrier() {
    let resolved = resolved_local_gas_and_condensed_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0, 1e-12]).unwrap();
    let solution = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
        &resolved,
        EquilibriumConditions::new(300.0, 3_583.115, 101_325.0).unwrap(),
        composition,
    ))
    .expect("direct H2O(g)/H2O(l) fixed-P,T solve must no longer need an O2 carrier");

    assert_eq!(solution.component_moles().len(), 2);
    assert!(
        solution
            .component_moles()
            .iter()
            .all(|moles| moles.is_finite() && *moles >= 0.0)
    );
}

#[test]
fn bridge_phase_control_starts_a_zero_condensed_candidate_as_inactive() {
    let resolved = resolved_local_gas_oxygen_and_condensed_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0]).unwrap();
    let result = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap()
    .solve_with_bounded_phase_control(|_| {}, |_| {})
    .unwrap();

    let phase_control = result
        .phase_control_report()
        .expect("bounded bridge solve must publish phase-control evidence");
    assert_eq!(
        phase_control.initial_phase_set.active_mask(),
        vec![true, false]
    );
    assert_eq!(
        result.phase_status(&PhaseId::new(Some("liquid".to_string()))),
        Some(crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus::Inactive)
    );
}

#[test]
fn bridge_problem_order_is_invariant_to_phase_map_insertion_order() {
    let first = resolved_local_gas_and_condensed_water_with_map_order(false);
    let second = resolved_local_gas_and_condensed_water_with_map_order(true);
    let first_layout = MultiphaseEquilibriumLayout::new(first.phase_specs().to_vec()).unwrap();
    let second_layout = MultiphaseEquilibriumLayout::new(second.phase_specs().to_vec()).unwrap();
    let conditions = EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap();

    let first_bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &first,
            conditions,
            MultiphaseInitialComposition::from_dense(&first_layout, vec![1.0, 0.0]).unwrap(),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let second_bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &second,
            conditions,
            MultiphaseInitialComposition::from_dense(&second_layout, vec![1.0, 0.0]).unwrap(),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    assert_eq!(first_bundle.metadata(), second_bundle.metadata());
    assert_eq!(
        first_bundle.problem().components(),
        second_bundle.problem().components()
    );
    assert_eq!(
        first_bundle.problem().species(),
        second_bundle.problem().species()
    );
    assert_eq!(
        first_bundle.problem().element_composition(),
        second_bundle.problem().element_composition()
    );
    assert_eq!(
        first_bundle.problem().initial_log_moles().as_slice(),
        second_bundle.problem().initial_log_moles().as_slice()
    );
    for (left, right) in first_bundle
        .problem()
        .gibbs()
        .iter()
        .zip(second_bundle.problem().gibbs())
    {
        assert_eq!(
            left(conditions.temperature()),
            right(conditions.temperature())
        );
    }
}

#[test]
fn physical_bridge_numeric_closure_and_symbolic_contracts_agree() {
    // This fixture uses real NASA gas and condensed records so the comparison
    // is not just a synthetic algebra study. It locks the bridge, closure, and
    // symbolic residual/Jacobian paths to the same physical formulation.
    let resolved = resolved_local_gas_oxygen_and_condensed_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0]).unwrap();
    let conditions = EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            conditions,
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    let symbolic_standard_gibbs = bundle.symbolic_standard_gibbs().to_vec();
    let prepared = PreparedEquilibriumProblem::new(bundle.into_problem()).unwrap();
    let log_moles = prepared.problem().initial_log_moles().as_slice().to_vec();
    let numeric_residual = prepared.residual(&log_moles).unwrap();
    let numeric_jacobian = prepared.jacobian(&log_moles).unwrap();

    let closure = equilibrium_logmole_residual(
        prepared.reaction_basis().reactions.clone(),
        prepared.problem().element_composition().clone(),
        prepared.element_totals().to_vec(),
        prepared.problem().gibbs().to_vec(),
        prepared.problem().phases().to_vec(),
        conditions.temperature(),
        conditions.pressure(),
        conditions.reference_pressure(),
        prepared.species_phase().to_vec(),
        1e-30,
        1e-30,
    )
    .unwrap();
    let closure_residual = closure(&log_moles).unwrap();
    assert_eq!(closure_residual, numeric_residual);

    let closure_jacobian = equilibrium_logmole_jacobian(
        &log_moles,
        &prepared.reaction_basis().reactions,
        prepared.problem().element_composition(),
        prepared.species_phase(),
        prepared.phase_stoichiometry(),
        prepared.problem().phases().len(),
        1e-30,
        1e-30,
    )
    .unwrap();
    assert_eq!(closure_jacobian, numeric_jacobian);

    let symbolic = multiphase_equilibrium_residual_generator_sym(
        prepared.reaction_basis().reactions.clone(),
        prepared.problem().element_composition().clone(),
        prepared.element_totals().to_vec(),
        symbolic_standard_gibbs,
        prepared.problem().phases().to_vec(),
        conditions.pressure(),
        conditions.reference_pressure(),
    )
    .unwrap();
    let variables = Expr::IndexedVars(prepared.problem().species().len(), "y").0;
    let names = variables
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let refs = names.iter().map(String::as_str).collect::<Vec<_>>();

    assert_eq!(symbolic.len(), numeric_residual.len());
    for (row, residual) in symbolic.iter().enumerate() {
        let evaluated = residual
            .clone()
            .set_variable("T", conditions.temperature())
            .simplify();
        let evaluate = evaluated.lambdify_borrowed_thread_safe(&refs);
        let symbolic_residual = evaluate(&log_moles);
        assert!(
            (symbolic_residual - numeric_residual[row]).abs() <= 1e-10,
            "physical residual {row} diverged: symbolic={symbolic_residual}, numeric={}",
            numeric_residual[row]
        );

        for (column, variable) in refs.iter().enumerate() {
            let derivative = evaluated
                .diff(variable)
                .set_variable("T", conditions.temperature());
            let evaluate = derivative.lambdify_borrowed_thread_safe(&refs);
            let symbolic_jacobian = evaluate(&log_moles);
            assert!(
                (symbolic_jacobian - numeric_jacobian[(row, column)]).abs() <= 1e-10,
                "physical Jacobian entry ({row}, {column}) diverged: symbolic={symbolic_jacobian}, numeric={}",
                numeric_jacobian[(row, column)]
            );
        }
    }
}

#[test]
fn parameterized_rst_graph_matches_the_canonical_real_nasa_formulation() {
    // The reusable temperature-range graph replaces baked `G0(T)` expressions
    // with parameters. This must be algebraically invisible to a backend.
    let bundle = prepared_local_nasa_gas();
    let symbolic_standard_gibbs = bundle.symbolic_standard_gibbs().to_vec();
    let prepared = PreparedEquilibriumProblem::new(bundle.into_problem()).unwrap();
    let log_moles = prepared.problem().initial_log_moles().as_slice().to_vec();
    let expected_residual = prepared.residual(&log_moles).unwrap();
    let expected_jacobian = prepared.jacobian(&log_moles).unwrap();
    let parameterized =
        prepare_rst_symbolic_problem_from_prepared(&prepared, &symbolic_standard_gibbs).unwrap();

    let actual_residual = parameterized.residual_for_test(&log_moles).unwrap();
    assert_eq!(actual_residual.len(), expected_residual.len());
    for (row, (&actual, &expected)) in actual_residual
        .iter()
        .zip(expected_residual.iter())
        .enumerate()
    {
        assert!(
            (actual - expected).abs() <= 1e-10,
            "parameterized residual row {row} diverged: actual={actual}, expected={expected}"
        );
    }

    let actual_jacobian = parameterized.jacobian_for_test(&log_moles).unwrap();
    assert_eq!(actual_jacobian.shape(), expected_jacobian.shape());
    for row in 0..actual_jacobian.nrows() {
        for column in 0..actual_jacobian.ncols() {
            assert!(
                (actual_jacobian[(row, column)] - expected_jacobian[(row, column)]).abs() <= 1e-10,
                "parameterized Jacobian entry ({row}, {column}) diverged: actual={}, expected={}",
                actual_jacobian[(row, column)],
                expected_jacobian[(row, column)]
            );
        }
    }
}

#[test]
fn missing_last_phase_data_returns_no_bundle_and_does_not_mutate_resolved_input() {
    let resolved = resolved_with_missing_last_phase();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0, 0.0]).unwrap();
    let gas_before = resolved
        .phase_data()
        .get(&Some("gas".to_string()))
        .expect("gas payload must exist");
    assert!(gas_before.element_composition_matrix().is_none());
    assert!(gas_before.therm_functions().is_empty());

    let result = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
            composition,
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    );
    let error = match result {
        Ok(_) => panic!("missing phase thermochemistry must reject bridge construction"),
        Err(error) => error,
    };

    assert!(
        error
            .to_string()
            .contains("equilibrium data preparation failed")
    );
    let gas_after = resolved
        .phase_data()
        .get(&Some("gas".to_string()))
        .expect("gas payload must still exist");
    assert!(gas_after.element_composition_matrix().is_none());
    assert!(gas_after.therm_functions().is_empty());
}
