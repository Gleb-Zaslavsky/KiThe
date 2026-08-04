//! Offline bridge tests for real resolved phase systems.
//!
//! These tests pin the production-facing contracts at the resolved-phase
//! boundary:
//!
//! - phase-qualified component identity stays ordered and deterministic;
//! - local NASA gas and condensed records keep distinct provenance;
//! - the same bare molecule can appear in multiple phases without losing its
//!   qualified identity; and
//! - the bounded phase-control path starts from the physical inventory rather
//!   than from a positive numerical trace seed.

use std::collections::HashMap;

use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    build_phase_equilibrium_problem, PhaseEquilibriumBuildRequest,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    solve_resolved_pt, ResolvedPhaseEquilibriumRequest,
};
use crate::Thermodynamics::User_PhaseOrSolution::{PhaseSpec, ResolvedPhaseSystem};
use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};

fn offline_nasa_gas_and_condensed_water() -> ResolvedPhaseSystem {
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

    let mut gas_data = SubsData::new();
    gas_data.substances = vec!["H2O".to_string(), "O2".to_string()];
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.set_substance_physical_state("H2O".to_string(), PhysicalState::Gas);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    let mut liquid_data = SubsData::new();
    liquid_data.substances = vec!["H2O".to_string()];
    liquid_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    liquid_data.set_substance_physical_state("H2O".to_string(), PhysicalState::Liquid);
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

fn offline_nasa_gas_and_condensed_solid_water() -> ResolvedPhaseSystem {
    let gas = PhaseSpec::ideal_gas(
        PhaseId::new(Some("gas".to_string())),
        vec!["H2O".to_string()],
    )
    .unwrap();
    let solid = PhaseSpec::pure_condensed(
        PhaseId::new(Some("solid".to_string())),
        vec!["H2O(s)".to_string()],
        PhysicalState::Solid,
    )
    .unwrap();

    let mut gas_data = SubsData::new();
    gas_data.substances = vec!["H2O".to_string()];
    gas_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);
    gas_data.search_substances().unwrap();
    gas_data.parse_all_thermal_coeffs().unwrap();

    let mut solid_data = SubsData::new();
    solid_data.substances = vec!["H2O(s)".to_string()];
    solid_data
        .set_multiple_library_priorities(vec!["NASA_cond".to_string()], LibraryPriority::Priority);
    solid_data.search_substances().unwrap();
    solid_data.parse_all_thermal_coeffs().unwrap();

    ResolvedPhaseSystem::new(
        vec![gas, solid],
        HashMap::from([
            (Some("gas".to_string()), gas_data),
            (Some("solid".to_string()), solid_data),
        ]),
    )
    .unwrap()
}

fn offline_problem(
) -> crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumProblemBundle
{
    let resolved = offline_nasa_gas_and_condensed_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0]).unwrap();

    build_phase_equilibrium_problem(
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
}

#[test]
fn offline_bridge_keeps_nasa_gas_and_condensed_provenance_separate() {
    let bundle = offline_problem();

    assert_eq!(
        bundle.problem().species(),
        ["gas::H2O", "gas::O2", "liquid::H2O"]
    );
    assert_eq!(
        bundle.report().components()[0].thermo_source().library(),
        "NASA_gas"
    );
    assert_eq!(
        bundle.report().components()[2].thermo_source().library(),
        "NASA_cond"
    );
    assert_eq!(
        bundle.report().components()[0].component().physical_state(),
        PhysicalState::Gas
    );
    assert_eq!(
        bundle.report().components()[2].component().physical_state(),
        PhysicalState::Liquid
    );
    assert_eq!(bundle.metadata().provenance().phases().len(), 2);
}

#[test]
fn offline_bridge_keeps_same_molecule_as_two_phase_qualified_components() {
    let bundle = offline_problem();

    assert_eq!(bundle.problem().components().len(), 3);
    assert_eq!(bundle.problem().species()[0], "gas::H2O");
    assert_eq!(bundle.problem().species()[2], "liquid::H2O");
    assert_eq!(
        bundle.problem().element_composition().row(0),
        bundle.problem().element_composition().row(2)
    );

    let gas_h2o = PhaseComponentId::new(PhaseId::new(Some("gas".to_string())), "H2O");
    let liquid_h2o = PhaseComponentId::new(PhaseId::new(Some("liquid".to_string())), "H2O");
    assert!(bundle.metadata().component_index(&gas_h2o).is_some());
    assert!(bundle.metadata().component_index(&liquid_h2o).is_some());
}

#[test]
fn offline_bridge_keeps_real_nasa_gas_and_condensed_solid_provenance_separate() {
    let resolved = offline_nasa_gas_and_condensed_solid_water();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let bundle = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            &resolved,
            EquilibriumConditions::new(250.0, 101_325.0, 101_325.0).unwrap(),
            MultiphaseInitialComposition::from_dense(&layout, vec![0.75, 0.25]).unwrap(),
            TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
            Default::default(),
        )
        .unwrap(),
    )
    .unwrap();

    assert_eq!(bundle.problem().species(), ["gas::H2O", "solid::H2O(s)"]);
    assert_eq!(
        bundle.report().components()[0].thermo_source().library(),
        "NASA_gas"
    );
    assert_eq!(
        bundle.report().components()[1].thermo_source().library(),
        "NASA_cond"
    );
    assert_eq!(
        bundle.report().components()[1].component().physical_state(),
        PhysicalState::Solid
    );
    assert_eq!(bundle.metadata().provenance().phases().len(), 2);
}

#[test]
fn offline_bounded_bridge_starts_zero_condensed_candidate_inactive() {
    let bundle = offline_problem();
    let result = bundle
        .solve_with_bounded_phase_control(|_| {}, |_| {})
        .unwrap();

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let phase_control = result
        .phase_control_report()
        .expect("bounded solve must retain phase-control evidence");
    assert_eq!(
        phase_control.initial_phase_set.active_mask(),
        vec![true, false]
    );
    assert_eq!(
        result.phase_status(&liquid),
        Some(crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus::Inactive)
    );
    assert!(result
        .summary_rows()
        .iter()
        .any(|row| row.section == "phase_control" && row.label == "iterations"));
    assert!(result.phase_total(&liquid).unwrap() >= 0.0);
}

#[test]
fn offline_bridge_one_shot_facade_preserves_transactional_publication() {
    let resolved = offline_nasa_gas_and_condensed_water();
    let before_gas = resolved
        .phase_data()
        .get(&Some("gas".to_string()))
        .unwrap()
        .clone();
    let before_liquid = resolved
        .phase_data()
        .get(&Some("liquid".to_string()))
        .unwrap()
        .clone();

    let result = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
        &resolved,
        EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0).unwrap(),
        MultiphaseInitialComposition::from_dense(
            &MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap(),
            vec![0.5, 0.25, 0.0],
        )
        .unwrap(),
    ));
    assert!(result.is_err());

    let after_gas = resolved.phase_data().get(&Some("gas".to_string())).unwrap();
    let after_liquid = resolved
        .phase_data()
        .get(&Some("liquid".to_string()))
        .unwrap();
    assert_eq!(before_gas.substances(), after_gas.substances());
    assert_eq!(
        before_gas.library_priorities(),
        after_gas.library_priorities()
    );
    assert_eq!(
        before_gas.explicit_search_map(),
        after_gas.explicit_search_map()
    );
    assert_eq!(before_liquid.substances(), after_liquid.substances());
    assert_eq!(
        before_liquid.library_priorities(),
        after_liquid.library_priorities()
    );
    assert_eq!(
        before_liquid.explicit_search_map(),
        after_liquid.explicit_search_map()
    );
}
