//! Offline stories for the public fixed-`P,T` multiphase result boundary.
//!
//! These stories check the hypotheses that matter before a future GUI or a
//! phase-control facade consumes a result:
//!
//! - a local NASA gas problem becomes one immutable phase-aware result;
//! - every numeric amount remains attached to its phase-qualified identity;
//! - phase totals and local mole fractions are derived from the accepted
//!   canonical solution rather than independently accumulated state; and
//! - the narrow one-shot facade follows the same bridge and does not mutate
//!   resolved thermochemical source data.

use std::collections::HashMap;

use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumSolverSettings;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    InitialPhaseSet, PhaseManager, PhaseStatus,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    build_phase_equilibrium_problem, PhaseEquilibriumBuildRequest,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    solve_resolved_pt, EquilibriumSolveOptions, PhaseControlPolicy, ResolvedPhaseEquilibriumRequest,
};
use crate::Thermodynamics::User_PhaseOrSolution::{PhaseSpec, ResolvedPhaseSystem};
use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};

fn resolved_local_nasa_gas() -> ResolvedPhaseSystem {
    let phase = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()],
    )
    .unwrap();
    let mut data = SubsData::new();
    data.substances = vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()];
    data.set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    data.search_substances().unwrap();
    data.parse_all_thermal_coeffs().unwrap();
    ResolvedPhaseSystem::new(
        vec![phase],
        HashMap::from([(Some("gas".to_string()), data)]),
    )
    .unwrap()
}

fn resolved_local_nasa_gas_with_two_condensed_candidates() -> ResolvedPhaseSystem {
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()],
    )
    .unwrap();
    let liquid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("liquid".to_string())),
        vec!["H2O".to_string()],
        crate::Thermodynamics::physical_state::PhysicalState::Liquid,
    )
    .unwrap();
    let solid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("solid".to_string())),
        vec!["H2O(s)".to_string()],
        crate::Thermodynamics::physical_state::PhysicalState::Solid,
    )
    .unwrap();

    let mut gas_data = SubsData::new();
    gas_data.substances = vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()];
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    let mut liquid_data = SubsData::new();
    liquid_data.substances = vec!["H2O".to_string()];
    liquid_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    liquid_data.set_substance_physical_state(
        "H2O".to_string(),
        crate::Thermodynamics::physical_state::PhysicalState::Liquid,
    );
    liquid_data.search_substances().unwrap();
    liquid_data.parse_all_thermal_coeffs().unwrap();

    let mut solid_data = SubsData::new();
    solid_data.substances = vec!["H2O(s)".to_string()];
    solid_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    solid_data.set_substance_physical_state(
        "H2O(s)".to_string(),
        crate::Thermodynamics::physical_state::PhysicalState::Solid,
    );
    solid_data.search_substances().unwrap();
    solid_data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(
        vec![gas, liquid, solid],
        HashMap::from([
            (Some("gas".to_string()), gas_data),
            (Some("liquid".to_string()), liquid_data),
            (Some("solid".to_string()), solid_data),
        ]),
    )
    .unwrap()
}

fn resolved_local_nasa_water_phase_pair() -> ResolvedPhaseSystem {
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2O".to_string(), "O2".to_string()],
    )
    .unwrap();
    let liquid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("liquid".to_string())),
        vec!["H2O".to_string()],
        crate::Thermodynamics::physical_state::PhysicalState::Liquid,
    )
    .unwrap();

    let mut gas_data = SubsData::new();
    gas_data.substances = vec!["H2O".to_string(), "O2".to_string()];
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    let mut liquid_data = SubsData::new();
    liquid_data.substances = vec!["H2O".to_string()];
    liquid_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    liquid_data.set_substance_physical_state(
        "H2O".to_string(),
        crate::Thermodynamics::physical_state::PhysicalState::Liquid,
    );
    liquid_data.search_substances().unwrap();
    liquid_data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(
        vec![gas, liquid],
        HashMap::from([
            (Some("gas".to_string()), gas_data),
            (Some("liquid".to_string()), liquid_data),
        ]),
    )
    .unwrap()
}

fn composition_for(resolved: &ResolvedPhaseSystem) -> MultiphaseInitialComposition {
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap()
}

fn conditions() -> EquilibriumConditions {
    EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap()
}

#[test]
fn accepted_fixed_phase_solution_exposes_qualified_amounts_totals_and_summary() {
    let resolved = resolved_local_nasa_gas();
    let result = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            // The bundled ice and liquid records meet at this phase boundary.
            EquilibriumConditions::new(273.15, 101_325.0, 101_325.0).unwrap(),
            composition_for(&resolved),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap()
    .solve()
    .unwrap()
    .into_multiphase_solution()
    .unwrap();

    let gas = PhaseId::new(Some("gas".to_string()));
    let h2 = PhaseComponentId::new(gas.clone(), "H2");
    let h2_amount = result.moles_for(&h2).unwrap();
    let local_fraction_sum = result
        .metadata()
        .layout()
        .components_for_phase(&gas)
        .unwrap()
        .iter()
        .map(|component| result.mole_fraction_for(component).unwrap())
        .sum::<f64>();

    assert!(h2_amount > 0.0);
    assert_eq!(result.phase_status(&gas), Some(PhaseStatus::Active));
    assert!((local_fraction_sum - 1.0).abs() < 1e-12);
    assert!(
        (result.phase_total(&gas).unwrap() - result.component_moles().iter().sum::<f64>()).abs()
            < 1e-12
    );
    assert!(result.aggregate_moles_by_substance().contains_key("H2O"));
    assert!(result.to_string().contains("[component] gas::H2"));
    assert!(result
        .summary_rows()
        .iter()
        .any(|row| row.section == "backend" && row.label == "accepted"));
}

#[test]
fn ideal_gas_facade_matches_direct_canonical_solve() {
    let resolved = resolved_local_nasa_gas();
    let conditions = conditions();
    let composition = composition_for(&resolved);

    let direct = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            conditions,
            composition.clone(),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap()
    .solve()
    .unwrap()
    .into_multiphase_solution()
    .unwrap();

    let facade = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
        &resolved,
        conditions,
        composition,
    ))
    .unwrap();

    assert_eq!(direct.component_moles(), facade.component_moles());
    assert_eq!(
        direct.aggregate_moles_by_substance(),
        facade.aggregate_moles_by_substance()
    );
    assert_eq!(direct.summary_rows(), facade.summary_rows());
    assert_eq!(
        direct.metadata().layout_fingerprint(),
        facade.metadata().layout_fingerprint()
    );
    assert_eq!(
        direct
            .metadata()
            .components()
            .iter()
            .map(|component| component.label())
            .collect::<Vec<_>>(),
        facade
            .metadata()
            .components()
            .iter()
            .map(|component| component.label())
            .collect::<Vec<_>>()
    );
}

#[test]
fn two_independent_pure_condensed_candidates_stay_distinct_in_the_bridge() {
    let resolved = resolved_local_nasa_gas_with_two_condensed_candidates();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(273.15, 101_325.0, 101_325.0).unwrap(),
            MultiphaseInitialComposition::from_dense(
                &MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap(),
                vec![2.0, 1.0, 0.0, 0.0, 0.0],
            )
            .unwrap(),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    assert_eq!(
        bundle.problem().species(),
        [
            "gas::H2",
            "gas::O2",
            "gas::H2O",
            "liquid::H2O",
            "solid::H2O(s)"
        ]
    );
    assert_eq!(
        bundle.report().components()[3].thermo_source().library(),
        "NASA_cond"
    );
    assert_eq!(
        bundle.report().components()[4].thermo_source().library(),
        "NASA_cond"
    );
    assert_eq!(bundle.metadata().provenance().phases().len(), 3);
}

#[test]
fn unsupported_multi_component_condensed_solution_is_rejected_before_bridge_building() {
    let phase = PhaseSpec::new(
        PhaseId::new(Some("bad_condensed".to_string())),
        vec!["A".to_string(), "B".to_string()],
        crate::Thermodynamics::physical_state::PhysicalState::Solid,
        crate::Thermodynamics::User_PhaseOrSolution::PhaseModel::PureCondensed,
    )
    .unwrap();
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["A".to_string(), "B".to_string()],
    )
    .unwrap();

    let error = MultiphaseEquilibriumLayout::new(vec![gas, phase]).unwrap_err();
    assert!(error.to_string().contains("exactly one component"));
}

#[test]
fn one_shot_facade_uses_the_canonical_bridge_without_mutating_resolved_data() {
    let resolved = resolved_local_nasa_gas();
    let before = resolved.phase_data().get(&Some("gas".to_string())).unwrap();
    assert!(before.element_composition_matrix().is_none());
    assert!(before.therm_functions().is_empty());

    let result = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
        &resolved,
        conditions(),
        composition_for(&resolved),
    ))
    .unwrap();

    assert_eq!(
        result
            .metadata()
            .components()
            .iter()
            .map(|component| component.label())
            .collect::<Vec<_>>(),
        ["gas::H2", "gas::O2", "gas::H2O"]
    );
    assert_eq!(result.build_report().components().len(), 3);
    assert!(result.solve_report().accepted_attempt().is_some());
    let after = resolved.phase_data().get(&Some("gas".to_string())).unwrap();
    assert!(after.element_composition_matrix().is_none());
    assert!(after.therm_functions().is_empty());
}

#[test]
fn default_resolved_solve_does_not_run_limited_keq_validator() {
    // K_eq is an independent diagnostic for its supported ideal-gas domain,
    // not a hidden production acceptance gate. The default public request must
    // therefore publish only the canonical solve evidence.
    let resolved = resolved_local_nasa_gas();
    let result = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
        &resolved,
        conditions(),
        composition_for(&resolved),
    ))
    .unwrap();

    assert!(result.keq_validation_status().is_none());
    assert!(!result
        .summary_rows()
        .iter()
        .any(|row| row.section == "keq_validation"));
}

#[test]
fn bounded_phase_control_publishes_its_acceptance_evidence_in_the_same_result() {
    let resolved = resolved_local_nasa_gas();
    let result = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(&resolved, conditions(), composition_for(&resolved))
            .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .unwrap();

    let phase_control = result
        .phase_control_report()
        .expect("bounded solve must retain its transition evidence");
    let acceptance = result
        .acceptance_report()
        .expect("bounded solve must retain its complementarity gate");
    assert_eq!(phase_control.final_phase_set.active_mask(), vec![true]);
    assert!(acceptance.complementarity.satisfied);
    assert!(result
        .summary_rows()
        .iter()
        .any(|row| row.section == "acceptance" && row.label == "complementarity_satisfied"));
}

#[test]
fn bounded_mixed_phase_control_publishes_keq_not_applicable_status() {
    let resolved = resolved_local_nasa_gas_with_two_condensed_candidates();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0, 0.0, 0.0]).unwrap();
    let mut settings = EquilibriumSolverSettings::default();
    settings.keq_validation_mode = EquilibriumConstantValidationMode::WhenApplicable;
    let mut phase_manager = PhaseManager::default();
    // Keep the numerical fixture on the well-conditioned gas-only branch. The
    // declared condensed phases still make K_eq validation inapplicable, so
    // this isolates the contract under test from a separate phase-transition
    // convergence question.
    phase_manager.initial_phase_set = InitialPhaseSet::Explicit {
        active: vec![PhaseIndex::new(0, 3).unwrap()],
        excluded: vec![
            PhaseIndex::new(1, 3).unwrap(),
            PhaseIndex::new(2, 3).unwrap(),
        ],
    };
    let phase_policy = PhaseControlPolicy::new(phase_manager).unwrap();

    let result = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(273.15, 101_325.0, 101_325.0).unwrap(),
            composition,
        )
        .with_solve_options(EquilibriumSolveOptions::from_settings(settings).unwrap())
        .with_phase_control_policy(phase_policy),
    )
    .unwrap();

    assert!(matches!(
        result.keq_validation_status(),
        Some(EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable { .. })
    ));
    assert!(result
        .summary_rows()
        .iter()
        .any(|row| row.section == "keq_validation" && row.value == "not_applicable"));

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let solid = PhaseId::new(Some("solid".to_string()));
    assert_eq!(result.phase_status(&liquid), Some(PhaseStatus::Excluded));
    assert_eq!(result.phase_status(&solid), Some(PhaseStatus::Excluded));
    assert_eq!(result.phase_total(&liquid), Some(0.0));
    assert_eq!(result.phase_total(&solid), Some(0.0));
    assert!(result.numerical_phase_total(&liquid).unwrap() > 0.0);
    assert!(result.numerical_phase_total(&solid).unwrap() > 0.0);
    let liquid_water = PhaseComponentId::new(liquid, "H2O");
    assert_eq!(result.moles_for(&liquid_water), Some(0.0));
    assert!(result.numerical_moles_for(&liquid_water).unwrap() > 0.0);
}

#[test]
fn required_keq_validation_rejects_bounded_mixed_phase_control_explicitly() {
    let resolved = resolved_local_nasa_gas_with_two_condensed_candidates();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0, 0.0, 0.0]).unwrap();
    let mut settings = EquilibriumSolverSettings::default();
    settings.keq_validation_mode = EquilibriumConstantValidationMode::Required;

    let error = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(273.15, 101_325.0, 101_325.0).unwrap(),
            composition,
        )
        .with_solve_options(EquilibriumSolveOptions::from_settings(settings).unwrap())
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::ValidationNotApplicable {
            path: "equilibrium_constant_cross_validation",
            ..
        }
    ));
}

#[test]
fn facade_retains_independent_keq_status_in_the_immutable_result_summary() {
    let resolved = resolved_local_nasa_gas();
    let mut settings = EquilibriumSolverSettings::default();
    settings.keq_validation_mode = EquilibriumConstantValidationMode::WhenApplicable;
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let warm_composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9]).unwrap();

    let result = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(2_500.0, 101_325.0, 101_325.0).unwrap(),
            warm_composition,
        )
        .with_solve_options(EquilibriumSolveOptions::from_settings(settings).unwrap()),
    )
    .unwrap();

    match result.keq_validation_status() {
        Some(EquilibriumConstantCrossValidationStatus::Compared(report)) => {
            assert!(report.accepted);
            assert!(report.max_abs_species_mole_delta <= 1e-5);
            assert!(report.max_abs_species_fraction_delta <= 1e-5);
        }
        other => panic!("expected a compared keq validation status, got {other:?}"),
    }
    assert!(result
        .summary_rows()
        .iter()
        .any(|row| row.section == "keq_validation" && row.label == "status"));
}

#[test]
fn gas_fixture_preserves_component_order_and_solution_across_legacy_and_rst_policies() {
    let resolved = resolved_local_nasa_gas();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0]).unwrap();

    let mut legacy_settings = EquilibriumSolverSettings::default();
    legacy_settings.solver = Solvers::LM;
    legacy_settings.solver_params.max_iter = 200;
    legacy_settings.solver_policy = Some(SolverPolicy::Cascade(vec![
        SolverBackend::Legacy(Solvers::LM),
        SolverBackend::Legacy(Solvers::TR),
        SolverBackend::Legacy(Solvers::NR),
    ]));

    let mut rst_settings = EquilibriumSolverSettings::default();
    rst_settings.solver_params.max_iter = 200;
    rst_settings.solver_policy = Some(SolverPolicy::Single(SolverBackend::RustedSciThe(
        RustedSciTheSolver::TrustRegionLevenbergMarquardt,
    )));

    let legacy = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(&resolved, conditions(), composition.clone())
            .with_solve_options(EquilibriumSolveOptions::from_settings(legacy_settings).unwrap()),
    )
    .unwrap();
    let rst = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(&resolved, conditions(), composition)
            .with_solve_options(EquilibriumSolveOptions::from_settings(rst_settings).unwrap()),
    )
    .unwrap();

    assert_eq!(legacy.metadata().components(), rst.metadata().components());
    assert_eq!(
        legacy.metadata().layout_fingerprint(),
        rst.metadata().layout_fingerprint()
    );
    assert_eq!(legacy.phases(), rst.phases());
    for component in legacy.metadata().components() {
        let legacy_moles = legacy.moles_for(&component.id()).unwrap();
        let rst_moles = rst.moles_for(&component.id()).unwrap();
        assert!(
            (legacy_moles - rst_moles).abs() <= 1e-8,
            "component {} diverged: legacy={}, rst={}",
            component.label(),
            legacy_moles,
            rst_moles
        );
    }
    assert_eq!(
        legacy
            .summary_rows()
            .iter()
            .filter(|row| row.section != "backend" && row.section != "validation")
            .map(|row| (row.section, row.label.as_str()))
            .collect::<Vec<_>>(),
        rst.summary_rows()
            .iter()
            .filter(|row| row.section != "backend" && row.section != "validation")
            .map(|row| (row.section, row.label.as_str()))
            .collect::<Vec<_>>()
    );
}

#[test]
fn water_phase_pair_switches_between_vapor_and_condensed_dominance_with_temperature() {
    let resolved = resolved_local_nasa_water_phase_pair();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0]).unwrap();

    let low = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(350.0, 101_325.0, 101_325.0).unwrap(),
            composition.clone(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .unwrap();

    let high = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(550.0, 101_325.0, 101_325.0).unwrap(),
            composition,
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .unwrap();

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let gas = PhaseId::new(Some("gas".to_string()));
    assert!(low.phase_total(&liquid).unwrap() >= high.phase_total(&liquid).unwrap());
    assert!(high.phase_total(&gas).unwrap() >= low.phase_total(&gas).unwrap());
    assert_eq!(
        low.metadata().layout_fingerprint(),
        high.metadata().layout_fingerprint()
    );
    assert_eq!(
        low.metadata().components().len(),
        layout.system_layout().components.len()
    );
}
