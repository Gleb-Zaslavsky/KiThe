//! Offline I1-I4 fixed-P,T validation on shared real pure-phase fixtures.
//!
//! These tests are the P,T peers of `pure_phase_ph_live_data_tests`: every
//! route consumes the same resolved identities, provenance, reaction vector,
//! elemental composition, and physical inventory. The independent route uses
//! `ln(Q)-ln(K)` plus scalar extent mathematics; the production route uses TPD
//! and bounded phase control. Agreement therefore checks two formulations,
//! not two calls into one solver.

use std::fs;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
    TemperatureGrid, TemperatureRangePointPreparation, TemperatureRangeRequest,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::InitialPhaseSet;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
    PurePhaseProductionEvidenceRequest, activation_evidence_from_solution,
    disappearance_evidence_from_solution, stable_inactive_evidence_from_solution,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    CanonicalPurePhaseEvidence, PurePhaseBoundaryStructuralTolerances, PurePhaseBoundaryValidator,
    PurePhaseCrossValidationStatus, PurePhaseCrossValidationTolerances, cross_validate_pure_phase,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    PhaseControlPolicy, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::ChemEquilibrium::real_pure_phase_fixtures::{
    RealPurePhaseFamily, RealPurePhaseGasScenario, RealPurePhaseInventory,
    ResolvedRealPurePhaseFixture,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
use crate::library_manager::with_library_manager;

const PRESSURE_PA: f64 = 101_325.0;

fn local_repository() -> Arc<ThermoRepository> {
    ThermoData::try_default_repository()
        .expect("the bundled offline thermochemistry repository must be available")
}

/// Captures every JSON file used by ordinary local resolution. Real P,T
/// validation is read-only and must leave these bytes unchanged.
fn local_library_snapshot() -> Vec<(String, Vec<u8>)> {
    let paths = with_library_manager(|manager| {
        vec![
            manager.substance_base_path().to_string(),
            manager.all_keys_substance_path().to_string(),
            manager.elements_path().to_string(),
        ]
    });
    paths
        .into_iter()
        .map(|path| {
            let bytes = fs::read(&path)
                .unwrap_or_else(|error| panic!("must read local P,T library '{path}': {error}"));
            (path, bytes)
        })
        .collect()
}

fn water_inventory(
    gas_water_moles: f64,
    oxygen_moles: f64,
    candidate_moles: f64,
) -> RealPurePhaseInventory {
    RealPurePhaseInventory::new(vec![gas_water_moles], candidate_moles)
        .and_then(|inventory| inventory.with_inert("O2", oxygen_moles, vec![0.0, 2.0]))
        .expect("water/O2 pure-phase inventory must validate")
}

fn solve_pt(
    fixture: &ResolvedRealPurePhaseFixture,
    inventory: &RealPurePhaseInventory,
    temperature: f64,
    policy: PhaseControlPolicy,
) -> MultiphaseEquilibriumSolution {
    let initial = fixture
        .initial_composition(inventory)
        .expect("real P,T inventory must match the resolved layout");
    solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            fixture.resolved(),
            EquilibriumConditions::new(temperature, PRESSURE_PA, PRESSURE_PA)
                .expect("real pure-phase P,T conditions must validate"),
            initial,
        )
        .with_phase_control_policy(policy),
    )
    .expect("real pure-phase P,T phase control must solve")
}

/// Compares one immutable production witness with a separately materialized
/// independent boundary problem at the exact recorded gas-only state.
fn assert_complete_cross_validation(
    fixture: &ResolvedRealPurePhaseFixture,
    boundary_gas: &RealPurePhaseGasScenario,
    temperature: f64,
    evidence: &CanonicalPurePhaseEvidence,
) {
    let problem = fixture
        .to_pt_boundary_problem(
            boundary_gas,
            EquilibriumConditions::new(temperature, PRESSURE_PA, PRESSURE_PA)
                .expect("cross-validation P,T conditions must validate"),
        )
        .expect("shared fixture must build the independent P,T problem");
    let reaction_space = problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .expect("real P,T family must retain one strict phase-forming reaction");
    assert_eq!(reaction_space.full_reaction_dimension, 1);
    assert_eq!(reaction_space.gas_only_reaction_dimension, 0);
    let report = cross_validate_pure_phase(
        &problem,
        evidence,
        PurePhaseBoundaryValidator::default(),
        PurePhaseCrossValidationTolerances {
            max_abs_gas_mole_delta: 1e-5,
            max_abs_candidate_mole_delta: 1e-5,
            max_abs_boundary_energy_delta: 1e-3,
            ..PurePhaseCrossValidationTolerances::default()
        },
    )
    .expect("production TPD and independent P,T validation must be comparable");
    assert_eq!(report.status, PurePhaseCrossValidationStatus::Complete);
    assert!(report.accepted, "incomplete real P,T evidence: {report:?}");
    assert_eq!(report.thermodynamic_agreement, Some(true));
    assert_eq!(report.topology_agreement, Some(true));
}

fn recorded_boundary_scenario(
    fixture: &ResolvedRealPurePhaseFixture,
    template: &RealPurePhaseInventory,
    evidence: &CanonicalPurePhaseEvidence,
) -> RealPurePhaseGasScenario {
    fixture
        .gas_scenario_from_canonical_boundary(
            template.gas(),
            &evidence.gas_species,
            &evidence.boundary_gas_moles,
        )
        .expect("production boundary identities must match the shared fixture")
}

#[test]
fn p8_i4_shared_water_pt_fixture_covers_appearance_absence_and_disappearance() {
    let before = local_library_snapshot();
    let fixture = RealPurePhaseFamily::WaterLiquid
        .resolve_offline(local_repository(), &["O2"])
        .expect("local water/liquid P,T fixture must resolve without NIST");
    let gas = PhaseId::new(Some("gas".to_string()));
    let liquid = PhaseId::new(Some("liquid".to_string()));
    let request = PurePhaseProductionEvidenceRequest::new(gas, liquid);

    // I1/I3/I4 appearance from a genuinely absent candidate.
    let absent = water_inventory(0.5, 0.25, 0.0);
    let low = solve_pt(&fixture, &absent, 350.0, PhaseControlPolicy::default());
    let appearance = activation_evidence_from_solution(&low, &request)
        .expect("low-temperature water must publish liquid activation evidence");
    let appearance_boundary = recorded_boundary_scenario(&fixture, &absent, &appearance);
    assert_complete_cross_validation(&fixture, &appearance_boundary, 350.0, &appearance);

    // I1/I3/I4 stable absence at the same material inventory.
    let high = solve_pt(&fixture, &absent, 550.0, PhaseControlPolicy::default());
    let inactive = stable_inactive_evidence_from_solution(&high, &request)
        .expect("high-temperature water must publish stable-inactive liquid evidence");
    let inactive_boundary = recorded_boundary_scenario(&fixture, &absent, &inactive);
    assert_complete_cross_validation(&fixture, &inactive_boundary, 550.0, &inactive);

    // I1/I3/I4 disappearance begins with physical liquid inventory and uses
    // the reduced-boundary restart state recorded by the transition.
    let present = water_inventory(0.5, 0.25, 0.1);
    let policy = PhaseControlPolicy::default()
        .with_initial_phase_set(InitialPhaseSet::AllCandidatePhases)
        .expect("all-candidate water policy must validate");
    let disappeared = solve_pt(&fixture, &present, 550.0, policy);
    let disappearance = disappearance_evidence_from_solution(&disappeared, &request)
        .expect("hot water must publish liquid deactivation evidence");
    let disappearance_boundary = recorded_boundary_scenario(&fixture, &present, &disappearance);
    assert_complete_cross_validation(&fixture, &disappearance_boundary, 550.0, &disappearance);

    assert_eq!(before, local_library_snapshot());
}

#[test]
fn p8_i4_shared_ice_and_boudouard_pt_fixtures_match_independent_boundaries() {
    let before = local_library_snapshot();

    // Water/ice adds a different condensed record identity and a low-
    // temperature phase appearance using the same inventory as the P,H case.
    let ice_fixture = RealPurePhaseFamily::WaterIce
        .resolve_offline(local_repository(), &["O2"])
        .expect("local water/ice P,T fixture must resolve without NIST");
    let ice_inventory = water_inventory(0.5, 0.25, 0.0);
    let ice = solve_pt(
        &ice_fixture,
        &ice_inventory,
        250.0,
        PhaseControlPolicy::default(),
    );
    let ice_evidence = activation_evidence_from_solution(
        &ice,
        &PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some("gas".to_string())),
            PhaseId::new(Some("solid".to_string())),
        ),
    )
    .expect("low-temperature water must publish ice activation evidence");
    let ice_boundary = recorded_boundary_scenario(&ice_fixture, &ice_inventory, &ice_evidence);
    assert_complete_cross_validation(&ice_fixture, &ice_boundary, 250.0, &ice_evidence);

    // Boudouard adds two reactive gas species and graphite, proving that the
    // shared layer is not a water-specific adapter.
    let carbon_fixture = RealPurePhaseFamily::BoudouardCarbon
        .resolve_offline(local_repository(), &[])
        .expect("local Boudouard P,T fixture must resolve without NIST");
    let carbon_inventory = RealPurePhaseInventory::new(vec![1.0, 0.1], 0.0)
        .expect("CO/CO2 gas-only inventory must validate");
    let carbon = solve_pt(
        &carbon_fixture,
        &carbon_inventory,
        700.0,
        PhaseControlPolicy::default(),
    );
    let carbon_evidence = activation_evidence_from_solution(
        &carbon,
        &PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some("gas".to_string())),
            PhaseId::new(Some("solid".to_string())),
        ),
    )
    .expect("low-temperature Boudouard system must publish graphite activation evidence");
    let carbon_boundary =
        recorded_boundary_scenario(&carbon_fixture, &carbon_inventory, &carbon_evidence);
    assert_complete_cross_validation(&carbon_fixture, &carbon_boundary, 700.0, &carbon_evidence);

    assert_eq!(before, local_library_snapshot());
}

#[test]
fn p10_1_shared_water_temperature_sweep_reuses_accepted_phase_control_state() {
    // This is deliberately a short story rather than another range benchmark:
    // it proves that one immutable real fixture can cross an appearance and a
    // disappearance while every continuation seed comes from an accepted point.
    let before = local_library_snapshot();
    let fixture = RealPurePhaseFamily::WaterLiquid
        .resolve_offline(local_repository(), &["O2"])
        .expect("local water fixture must resolve without NIST");
    let inventory = water_inventory(0.5, 0.25, 0.0);
    let initial = fixture
        .initial_composition(&inventory)
        .expect("water range inventory must match the resolved layout");
    let range = TemperatureRangeRequest::new(
        fixture.resolved(),
        initial,
        PRESSURE_PA,
        PRESSURE_PA,
        TemperatureGrid::new(vec![350.0, 450.0, 550.0])
            .expect("water temperature grid must validate"),
    )
    .expect("water range request must validate")
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve()
    .expect("water phase-control range must accept every point transactionally");

    assert_eq!(range.points().len(), 3);
    assert!(range.report().phase_control_enabled());
    assert_eq!(
        range.points()[0].report().preparation(),
        TemperatureRangePointPreparation::InitialFormulation
    );
    for point in &range.points()[1..] {
        assert_eq!(
            point.report().preparation(),
            TemperatureRangePointPreparation::ReusedFormulation
        );
        assert!(point.report().used_continuation_seed());
        assert!(point.report().phase_set_reused());
    }

    let request = PurePhaseProductionEvidenceRequest::new(
        PhaseId::new(Some("gas".to_string())),
        PhaseId::new(Some("liquid".to_string())),
    );
    let appeared = activation_evidence_from_solution(range.points()[0].solution(), &request)
        .expect("first water range point must record liquid appearance");
    let appearance_boundary = recorded_boundary_scenario(&fixture, &inventory, &appeared);
    assert_complete_cross_validation(&fixture, &appearance_boundary, 350.0, &appeared);

    // The phase may deactivate at either continuation point. The final state
    // is nevertheless a valid independent gas-only witness without assuming
    // which neighbouring point owns the transition record.
    let inactive = stable_inactive_evidence_from_solution(range.points()[2].solution(), &request)
        .expect("last water range point must publish stable-inactive liquid evidence");
    let inactive_boundary = recorded_boundary_scenario(&fixture, &inventory, &inactive);
    assert_complete_cross_validation(&fixture, &inactive_boundary, 550.0, &inactive);
    assert_eq!(before, local_library_snapshot());
}
