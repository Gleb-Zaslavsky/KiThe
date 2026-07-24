//! Story tests for the typed multiphase equilibrium input boundary.
//!
//! These tests verify that phase-qualified identity is preserved, physical
//! zeroes are not confused with phase exclusion, and unsupported activity
//! models fail before numerical solver construction.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::User_PhaseOrSolution::{PhaseModel, PhaseSpec};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use nalgebra::DMatrix;

fn gas(name: &str, components: &[&str]) -> PhaseSpec {
    PhaseSpec::new(
        PhaseId::new(Some(name.to_string())),
        components
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
        PhysicalState::Gas,
        PhaseModel::IdealGas,
    )
    .unwrap()
}

fn solid(name: &str, component: &str) -> PhaseSpec {
    PhaseSpec::pure_condensed(
        PhaseId::new(Some(name.to_string())),
        vec![component.to_string()],
        PhysicalState::Solid,
    )
    .unwrap()
}

#[test]
fn same_substance_in_gas_and_solid_remains_two_components() {
    let layout =
        MultiphaseEquilibriumLayout::new(vec![gas("gas", &["H2O"]), solid("ice", "H2O")]).unwrap();

    assert_eq!(layout.component_count(), 2);
    assert_eq!(layout.system_layout().components[0].label(), "gas::H2O");
    assert_eq!(layout.system_layout().components[1].label(), "ice::H2O");
}

#[test]
fn rejects_two_ideal_gas_phases_before_solver_setup() {
    let error = MultiphaseEquilibriumLayout::new(vec![gas("gas_a", &["A"]), gas("gas_b", &["B"])])
        .unwrap_err();
    assert!(error.to_string().contains("at most one ideal-gas phase"));
}

#[test]
fn rejects_multicomponent_pure_condensed_phase_before_solver_setup() {
    let phase = PhaseSpec::new(
        PhaseId::new(Some("bad_solid".to_string())),
        vec!["A".to_string(), "B".to_string()],
        PhysicalState::Solid,
        PhaseModel::PureCondensed,
    )
    .unwrap();
    let error = MultiphaseEquilibriumLayout::new(vec![phase]).unwrap_err();
    assert!(error.to_string().contains("exactly one component"));
}

#[test]
fn sparse_and_dense_compositions_are_identical_and_preserve_zero_candidates() {
    let layout =
        MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A", "B"]), solid("solid", "A")])
            .unwrap();
    let dense = MultiphaseInitialComposition::from_dense(&layout, vec![1.0, 0.0, 0.0]).unwrap();
    let sparse = MultiphaseInitialComposition::from_sparse(
        &layout,
        vec![(
            PhaseComponentId::new(PhaseId::new(Some("gas".to_string())), "A"),
            1.0,
        )],
    )
    .unwrap();

    assert_eq!(dense, sparse);
    assert_eq!(sparse.moles(), &[1.0, 0.0, 0.0]);
}

#[test]
fn composition_rejects_duplicates_invalid_amounts_and_all_zero_input() {
    let layout = MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A", "B"])]).unwrap();
    let a = PhaseComponentId::new(PhaseId::new(Some("gas".to_string())), "A");

    assert!(
        MultiphaseInitialComposition::from_sparse(&layout, vec![(a.clone(), 1.0), (a, 2.0)])
            .is_err()
    );
    assert!(MultiphaseInitialComposition::from_dense(&layout, vec![-1.0, 1.0]).is_err());
    assert!(MultiphaseInitialComposition::from_dense(&layout, vec![0.0, 0.0]).is_err());
}

#[test]
fn physical_element_totals_are_computed_before_trace_seed_conversion() {
    let layout =
        MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A", "B"]), solid("solid", "A")])
            .unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 3.0, 0.0]).unwrap();
    let elements = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 2.0, 1.0, 0.0]);

    assert_eq!(
        composition.element_totals(&layout, &elements).unwrap(),
        vec![2.0, 6.0]
    );
}

#[test]
fn layout_fingerprint_is_deterministic_and_rejects_foreign_composition() {
    let first = MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A"])]).unwrap();
    let same = MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A"])]).unwrap();
    let other = MultiphaseEquilibriumLayout::new(vec![gas("gas", &["B"])]).unwrap();
    let composition = MultiphaseInitialComposition::from_dense(&first, vec![1.0]).unwrap();

    assert_eq!(first.fingerprint(), same.fingerprint());
    assert_ne!(first.fingerprint(), other.fingerprint());
    assert!(composition.validate_for(&other).is_err());
}

#[test]
fn composition_rejects_reconstruction_against_a_foreign_layout() {
    let first =
        MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A", "B"]), solid("solid", "A")])
            .unwrap();
    let other =
        MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A", "B"]), solid("solid", "B")])
            .unwrap();
    let composition =
        MultiphaseInitialComposition::from_dense(&first, vec![1.0, 0.0, 0.0]).unwrap();
    let elements = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 2.0, 1.0, 0.0]);

    let error = composition.element_totals(&other, &elements).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("composition belongs to a different multiphase layout")
    );
}

#[test]
fn phase_declaration_order_does_not_change_the_canonical_component_order() {
    let first =
        MultiphaseEquilibriumLayout::new(vec![solid("solid", "A"), gas("gas", &["A", "B"])])
            .unwrap();
    let second =
        MultiphaseEquilibriumLayout::new(vec![gas("gas", &["A", "B"]), solid("solid", "A")])
            .unwrap();

    assert_eq!(
        first.system_layout().components,
        second.system_layout().components
    );
    assert_eq!(first.fingerprint(), second.fingerprint());
}
