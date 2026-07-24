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

use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    PhaseEquilibriumBuildRequest, build_phase_equilibrium_problem,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::User_PhaseOrSolution::{PhaseSpec, ResolvedPhaseSystem};
use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

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
            conditions(),
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
    assert!(
        result
            .summary_rows()
            .iter()
            .any(|row| row.section == "backend" && row.label == "accepted")
    );
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
