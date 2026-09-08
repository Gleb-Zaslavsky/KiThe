//! I4 offline real-data tests for the strict pure-phase `P,H` validator.
//!
//! This module deliberately owns only the bridge from an immutable resolved
//! local record to the independent scalar validator. It does not duplicate the
//! broader water/ice phase-control stories in `equilibrium_live_data_tests`.
//! A passing test proves that the same local `G(T)` and `H(T)` capabilities can
//! drive both independent I2 mathematics and the canonical fixed-topology P,H
//! route without a network request or library mutation.

use std::fs;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EquilibriumConstraint, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::MultiphaseInitialComposition;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    PhSolveMode, ResolvedPhaseEnthalpyRequest, solve_resolved_ph,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
    PurePhaseProductionEvidenceRequest, disappearance_evidence_from_solution,
    stable_inactive_evidence_from_solution,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::PurePhaseBoundaryStructuralTolerances;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    PurePhaseBoundaryPrediction, PurePhaseBoundaryTolerances, PurePhaseBoundaryValidator,
    PurePhaseCrossValidationTolerances, cross_validate_pure_phase, evaluate_pure_phase_boundary,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    PhaseControlPolicy, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::ChemEquilibrium::pure_phase_ph_validation::{
    PurePhasePhCanonicalEvidence, PurePhasePhCrossValidationTolerances, PurePhasePhValidator,
    compare_pure_phase_ph_validation,
};
use crate::Thermodynamics::ChemEquilibrium::real_pure_phase_fixtures::{
    RealPurePhaseFamily, RealPurePhaseInventory, ResolvedRealPurePhaseFixture,
};
use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
use crate::library_manager::with_library_manager;

const PRESSURE_PA: f64 = 101_325.0;
const REFERENCE_TEMPERATURE_K: f64 = 350.0;
const INTERIOR_TEMPERATURE_LOWER_K: f64 = 345.0;
const INTERIOR_TEMPERATURE_UPPER_K: f64 = 355.0;

/// Captures every local JSON library ordinary resolution may read. I4 tests
/// compare this byte snapshot after the solve so a passing regression cannot
/// conceal an accidental write workflow.
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
                .unwrap_or_else(|error| panic!("must read local I4 library '{path}': {error}"));
            (path, bytes)
        })
        .collect()
}

fn local_repository() -> Arc<ThermoRepository> {
    ThermoData::try_default_repository()
        .expect("the bundled offline thermochemistry repository must be available")
}

fn local_water_fixture() -> ResolvedRealPurePhaseFixture {
    RealPurePhaseFamily::WaterLiquid
        .resolve_offline(local_repository(), &["O2"])
        .expect("the pinned local water fixture must resolve without NIST")
}

fn local_ice_fixture() -> ResolvedRealPurePhaseFixture {
    RealPurePhaseFamily::WaterIce
        .resolve_offline(local_repository(), &["O2"])
        .expect("the pinned local ice fixture must resolve without NIST")
}

/// Common water inventory used by both the canonical and independent routes.
/// O2 stays explicit: it contributes to gas activities and oxygen inventory,
/// but has exactly zero coefficient in the H2O phase-forming reaction.
fn local_water_scenario() -> RealPurePhaseInventory {
    RealPurePhaseInventory::new(vec![0.5], 0.1)
        .and_then(|scenario| scenario.with_inert("O2", 0.25, vec![0.0, 2.0]))
        .expect("the local water/O2 scenario must be structurally valid")
}

#[test]
fn p9_i4_local_water_fixed_topology_scalar_route_matches_canonical_ph() {
    // I4: all data are resolved from the pinned local NASA gas/condensed
    // catalogs. The P,T point only constructs a reproducible enthalpy target;
    // it is not the P,H reference being asserted below.
    let before = local_library_snapshot();
    let fixture = local_water_fixture();
    let scenario = local_water_scenario();
    assert_eq!(fixture.family(), RealPurePhaseFamily::WaterLiquid);
    let resolved = fixture.resolved();
    let report = resolved.report();
    assert!(!report.nist_fallback_enabled());
    let mut local_sources = report
        .phases()
        .iter()
        .flat_map(|phase| phase.search().rows())
        .filter(|row| row.property() == "Thermo")
        .map(|row| (row.substance().to_string(), row.library().to_string()))
        .collect::<Vec<_>>();
    local_sources.sort();
    assert_eq!(
        local_sources,
        vec![
            ("H2O".to_string(), "NASA_cond".to_string()),
            ("H2O".to_string(), "NASA_gas".to_string()),
            ("O2".to_string(), "NASA_gas".to_string()),
        ],
        "I4 must use the pinned local NASA gas/condensed records"
    );
    assert_eq!(
        resolved.layout().component_labels(),
        ["gas::H2O", "gas::O2", "liquid::H2O"]
    );

    let layout = crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::MultiphaseEquilibriumLayout::new(
        resolved.phase_specs().to_vec(),
    )
    .expect("local water layout must validate");
    let initial = fixture
        .initial_composition(&scenario)
        .expect("fixture must assemble local water moles in the resolved layout");
    let thermochemistry = Arc::new(fixture.thermochemistry().clone());

    let pt = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            resolved,
            EquilibriumConditions::new(REFERENCE_TEMPERATURE_K, PRESSURE_PA, PRESSURE_PA)
                .expect("reference P,T conditions must validate"),
            initial,
        )
        .with_fixed_declared_phases(),
    )
    .expect("fixed local water P,T point must solve");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(pt.component_moles(), REFERENCE_TEMPERATURE_K)
        .expect("accepted local P,T state must have finite enthalpy");

    let independent_problem = fixture
        .to_ph_problem(
            &scenario,
            PRESSURE_PA,
            PRESSURE_PA,
            target_enthalpy,
            TemperatureBounds::new(INTERIOR_TEMPERATURE_LOWER_K, INTERIOR_TEMPERATURE_UPPER_K)
                .expect("local water P,H bounds must validate"),
        )
        .expect("fixture adapter must build the local water independent P,H problem");
    let reaction_space = independent_problem
        .reaction_space(PurePhaseBoundaryStructuralTolerances::default())
        .expect("local water reaction-space evidence must validate");
    assert_eq!(reaction_space.full_reaction_dimension, 1);
    assert_eq!(reaction_space.gas_only_reaction_dimension, 0);
    let independent = PurePhasePhValidator::default()
        .solve(&independent_problem)
        .expect("the local water I2 scalar P,H route must solve");

    // Seed canonical P,H from the accepted P,T composition only to avoid
    // testing an arbitrary numerical starting point. Both paths retain the
    // same original elemental inventory and target enthalpy.
    let canonical_initial =
        MultiphaseInitialComposition::from_dense(&layout, pt.component_moles().to_vec())
            .expect("accepted P,T seed must match the local water layout");
    let canonical = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            canonical_initial,
            EquilibriumConstraint::ph(
                PRESSURE_PA,
                PRESSURE_PA,
                target_enthalpy,
                REFERENCE_TEMPERATURE_K,
            )
            .expect("local water P,H constraint must validate"),
            TemperatureBounds::new(INTERIOR_TEMPERATURE_LOWER_K, INTERIOR_TEMPERATURE_UPPER_K)
                .expect("local water P,H bounds must validate"),
            (*thermochemistry).clone(),
        )
        .expect("canonical local water P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature),
    )
    .expect("canonical fixed-topology local water P,H must solve");
    let canonical_moles = canonical.equilibrium().component_moles();
    let chemical_log_residual = independent_problem
        .chemical_log_residual_for_gas_moles(&canonical_moles[..2], canonical.temperature())
        .expect("canonical local water gas state must remain a valid I1 state")
        .abs();
    let max_abs_element_balance = [
        (2.0 * canonical_moles[0] + 2.0 * canonical_moles[2] - 1.2).abs(),
        (canonical_moles[0] + 2.0 * canonical_moles[1] + canonical_moles[2] - 1.1).abs(),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    let evidence = PurePhasePhCanonicalEvidence::new(
        independent_problem.case_identity(),
        canonical.temperature(),
        canonical_moles[..2].to_vec(),
        canonical_moles[2],
        canonical.calculated_enthalpy(),
        Some(chemical_log_residual),
        Some(max_abs_element_balance),
    )
    .expect("canonical local water evidence must remain finite");
    let comparison = compare_pure_phase_ph_validation(
        &independent_problem,
        &independent,
        &evidence,
        PurePhasePhCrossValidationTolerances {
            max_abs_temperature_delta: 1e-3,
            max_abs_mole_delta: 1e-5,
            max_abs_enthalpy_delta: 1e-2,
            max_abs_chemical_log_residual: 1e-6,
            max_abs_element_balance: 1e-8,
            ..PurePhasePhCrossValidationTolerances::default()
        },
    )
    .expect("local water I4 comparison must be well formed");
    assert!(
        comparison.is_complete_match(),
        "local water P9 I4 fixed-topology comparison: {comparison:?}"
    );
    assert_eq!(before, local_library_snapshot());
}

#[test]
fn p9_i1_i3_local_hot_water_ph_keeps_liquid_inactive_and_matches_boundary_tpd() {
    // I1/I3: derive a physical high-temperature enthalpy target from the
    // immutable local records, then let the production P,H outer loop decide
    // whether the zero-inventory liquid phase stays absent. The independent
    // boundary route shares only G(T), not phase-control code or a solver.
    let before = local_library_snapshot();
    let fixture = local_water_fixture();
    let scenario = RealPurePhaseInventory::new(vec![0.5], 0.0)
        .and_then(|scenario| scenario.with_inert("O2", 0.25, vec![0.0, 2.0]))
        .expect("the inventory-free hot-water scenario must validate");
    let resolved = fixture.resolved();
    let thermochemistry = fixture.thermochemistry();
    let initial = fixture
        .initial_composition(&scenario)
        .expect("hot-water moles must match the resolved layout");
    let reference_temperature = 550.0;
    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            resolved,
            EquilibriumConditions::new(reference_temperature, PRESSURE_PA, PRESSURE_PA)
                .expect("hot-water P,T conditions must validate"),
            initial.clone(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("hot-water P,T phase control must solve");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), reference_temperature)
        .expect("accepted hot-water P,T state must have finite enthalpy");

    let canonical = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            resolved,
            initial,
            EquilibriumConstraint::ph(
                PRESSURE_PA,
                PRESSURE_PA,
                target_enthalpy,
                reference_temperature,
            )
            .expect("hot-water P,H constraint must validate"),
            TemperatureBounds::new(545.0, 555.0).expect("hot-water P,H bounds must validate"),
            thermochemistry.clone(),
        )
        .expect("hot-water P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("hot-water P,H phase control must keep the liquid candidate absent");

    let liquid = crate::Thermodynamics::phase_layout::PhaseId::new(Some("liquid".to_string()));
    assert_eq!(canonical.equilibrium().phase_total(&liquid), Some(0.0));
    let boundary_problem = fixture
        .to_pt_boundary_problem(
            scenario.gas(),
            EquilibriumConditions::new(canonical.temperature(), PRESSURE_PA, PRESSURE_PA)
                .expect("canonical P,H temperature must produce valid P,T conditions"),
        )
        .expect("hot-water independent P,T boundary problem must validate");
    let boundary =
        evaluate_pure_phase_boundary(&boundary_problem, PurePhaseBoundaryTolerances::default())
            .expect("hot-water independent boundary must evaluate");
    assert_eq!(
        boundary.prediction,
        PurePhaseBoundaryPrediction::StableInactive,
        "the independent high-temperature water boundary must reject liquid formation: {boundary:?}"
    );
    let evidence = stable_inactive_evidence_from_solution(
        canonical.equilibrium(),
        &PurePhaseProductionEvidenceRequest::new(
            crate::Thermodynamics::phase_layout::PhaseId::new(Some("gas".to_string())),
            liquid,
        ),
    )
    .expect("accepted P,H result must publish stable-inactive liquid TPD evidence");
    let comparison = cross_validate_pure_phase(
        &boundary_problem,
        &evidence,
        PurePhaseBoundaryValidator::default(),
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect("I1/I3 hot-water cross-validation must be well formed");
    assert!(
        comparison.accepted,
        "canonical P,H stable-inactive liquid evidence must agree with independent boundary physics: {comparison:?}"
    );
    assert_eq!(before, local_library_snapshot());
}

#[test]
fn p9_i1_i3_i4_local_hot_water_ph_deactivates_initial_liquid_transactionally() {
    // This is deliberately distinct from stable absence: liquid begins with a
    // physical positive inventory. The P,H outer loop must reduce to the
    // gas-only boundary, preserve the elemental inventory, and publish the
    // deactivation record that carries its own TPD witness.
    let before = local_library_snapshot();
    let fixture = local_water_fixture();
    let scenario = RealPurePhaseInventory::new(vec![0.5], 0.1)
        .and_then(|scenario| scenario.with_inert("O2", 0.25, vec![0.0, 2.0]))
        .expect("the initially-liquid hot-water scenario must validate");
    let resolved = fixture.resolved();
    let thermochemistry = fixture.thermochemistry();
    let initial = fixture
        .initial_composition(&scenario)
        .expect("initial liquid-water moles must match the resolved layout");
    let reference_temperature = 550.0;
    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            resolved,
            EquilibriumConditions::new(reference_temperature, PRESSURE_PA, PRESSURE_PA)
                .expect("hot-water P,T conditions must validate"),
            initial.clone(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("hot-water P,T boundary recovery must solve");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), reference_temperature)
        .expect("accepted hot-water P,T state must have finite enthalpy");
    let canonical = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            resolved,
            initial,
            EquilibriumConstraint::ph(
                PRESSURE_PA,
                PRESSURE_PA,
                target_enthalpy,
                reference_temperature,
            )
            .expect("hot-water P,H constraint must validate"),
            TemperatureBounds::new(545.0, 555.0).expect("hot-water P,H bounds must validate"),
            thermochemistry.clone(),
        )
        .expect("hot-water P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("hot-water P,H must recover the gas-only boundary");

    let gas = crate::Thermodynamics::phase_layout::PhaseId::new(Some("gas".to_string()));
    let liquid = crate::Thermodynamics::phase_layout::PhaseId::new(Some("liquid".to_string()));
    assert_eq!(canonical.equilibrium().phase_total(&liquid), Some(0.0));
    let lifecycle = canonical
        .equilibrium()
        .phase_control_report()
        .expect("P,H disappearance must retain phase-control evidence");
    assert_eq!(
        lifecycle.transitions.len(),
        1,
        "one initial liquid candidate requires exactly one deactivation: {lifecycle:?}"
    );
    assert_eq!(lifecycle.transitions[0].activated.len(), 0);
    assert_eq!(lifecycle.transitions[0].deactivated.len(), 1);

    let evidence = disappearance_evidence_from_solution(
        canonical.equilibrium(),
        &PurePhaseProductionEvidenceRequest::new(gas, liquid),
    )
    .expect("accepted P,H disappearance must expose immutable boundary evidence");
    // A reduced-boundary deactivation records its own gas-only restart state.
    // Re-evaluate independent thermodynamics at that state, rather than at the
    // original two-phase inventory, because those are different physical
    // points after the candidate material has been projected away.
    let boundary_scenario = RealPurePhaseInventory::new(vec![evidence.boundary_gas_moles[0]], 0.0)
        .and_then(|scenario| {
            scenario.with_inert("O2", evidence.boundary_gas_moles[1], vec![0.0, 2.0])
        })
        .expect("recorded water boundary gas state must remain structurally valid");
    let boundary_problem = fixture
        .to_pt_boundary_problem(
            boundary_scenario.gas(),
            EquilibriumConditions::new(canonical.temperature(), PRESSURE_PA, PRESSURE_PA)
                .expect("canonical P,H temperature must produce valid P,T conditions"),
        )
        .expect("hot-water independent P,T boundary problem must validate");
    let boundary =
        evaluate_pure_phase_boundary(&boundary_problem, PurePhaseBoundaryTolerances::default())
            .expect("hot-water independent boundary must evaluate");
    assert_eq!(
        boundary.prediction,
        PurePhaseBoundaryPrediction::StableInactive
    );
    let comparison = cross_validate_pure_phase(
        &boundary_problem,
        &evidence,
        PurePhaseBoundaryValidator::default(),
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect("I1/I3/I4 hot-water disappearance comparison must be well formed");
    assert!(
        comparison.accepted,
        "canonical P,H disappearance must agree with independent boundary physics: {comparison:?}"
    );
    assert_eq!(before, local_library_snapshot());
}

#[test]
fn p9_i4_local_ice_phase_control_activates_solid_and_matches_independent_ph() {
    // I3/I4: solid water starts with zero physical inventory. The bounded
    // lifecycle must create it on thermodynamic evidence, then publish the
    // same state accepted by the independent scalar I2 calculation.
    let before = local_library_snapshot();
    let fixture = local_ice_fixture();
    // The target below is derived from 0.5 mol total water with no initial
    // solid inventory. The independent problem must carry the identical
    // material inventory; adding the 0.1 mol candidate used by the liquid
    // fixed-topology fixture shifts H by roughly 29.5 kJ and destroys the
    // otherwise valid outer temperature bracket.
    let scenario = RealPurePhaseInventory::new(vec![0.5], 0.0)
        .and_then(|scenario| scenario.with_inert("O2", 0.25, vec![0.0, 2.0]))
        .expect("the gas-only local ice appearance scenario must validate");
    let resolved = fixture.resolved();
    let thermochemistry = Arc::new(fixture.thermochemistry().clone());
    let reference_initial = fixture
        .initial_composition(&scenario)
        .expect("gas-only P,T reference must match the local layout");
    // Use the same immutable resolved payload for the P,T source point and
    // both P,H routes. Re-resolving through a second facade would weaken the
    // I4 guarantee that all compared paths consume identical local records.
    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(250.0, PRESSURE_PA, PRESSURE_PA)
                .expect("reference P,T conditions must validate"),
            reference_initial,
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("phase-controlled ice P,T reference must solve");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), 250.0)
        .expect("accepted phase-controlled ice P,T state must have finite enthalpy");
    let independent_problem = fixture
        .to_ph_problem(
            &scenario,
            PRESSURE_PA,
            PRESSURE_PA,
            target_enthalpy,
            TemperatureBounds::new(245.0, 255.0).expect("local ice P,H bounds must validate"),
        )
        .expect("fixture adapter must build the local ice independent P,H problem");
    let independent = PurePhasePhValidator::default()
        .solve(&independent_problem)
        .expect("the local ice I2 scalar P,H route must solve");

    let gas_only_initial = fixture
        .initial_composition(&scenario)
        .expect("gas-only local ice input must match the layout");
    let canonical = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            gas_only_initial,
            EquilibriumConstraint::ph(PRESSURE_PA, PRESSURE_PA, target_enthalpy, 250.0)
                .expect("local ice P,H constraint must validate"),
            TemperatureBounds::new(245.0, 255.0).expect("local ice P,H bounds must validate"),
            (*thermochemistry).clone(),
        )
        .expect("phase-controlled local ice P,H request must validate")
        // The coupled monolithic attempt may reject this topology change; the
        // production `Auto` route must retain that fact and recover through
        // the guarded nested phase-control path.
        .with_ph_solve_mode(PhSolveMode::Auto)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("local ice P,H phase control must activate solid");
    let phase_control = canonical
        .equilibrium()
        .phase_control_report()
        .expect("bounded P,H solve must publish lifecycle evidence");
    assert_eq!(phase_control.transitions.len(), 1);
    assert!(
        phase_control.transitions[0]
            .activated
            .iter()
            .any(|phase| phase.index() == 1),
        "the only local ice lifecycle transition must activate solid: {phase_control:?}"
    );

    let canonical_moles = canonical.equilibrium().component_moles();
    let chemical_log_residual = independent_problem
        .chemical_log_residual_for_gas_moles(&canonical_moles[..2], canonical.temperature())
        .expect("canonical gas state must remain a valid I1 state")
        .abs();
    let max_abs_element_balance = [
        (2.0 * canonical_moles[0] + 2.0 * canonical_moles[2] - 1.0).abs(),
        (canonical_moles[0] + 2.0 * canonical_moles[1] + canonical_moles[2] - 1.0).abs(),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    let evidence = PurePhasePhCanonicalEvidence::new(
        independent_problem.case_identity(),
        canonical.temperature(),
        canonical_moles[..2].to_vec(),
        canonical_moles[2],
        canonical.calculated_enthalpy(),
        Some(chemical_log_residual),
        Some(max_abs_element_balance),
    )
    .expect("phase-controlled local ice evidence must remain finite");
    let comparison = compare_pure_phase_ph_validation(
        &independent_problem,
        &independent,
        &evidence,
        PurePhasePhCrossValidationTolerances {
            max_abs_temperature_delta: 1e-3,
            max_abs_mole_delta: 1e-5,
            max_abs_enthalpy_delta: 1e-2,
            max_abs_chemical_log_residual: 1e-6,
            max_abs_element_balance: 1e-8,
            ..PurePhasePhCrossValidationTolerances::default()
        },
    )
    .expect("local ice phase-control I4 comparison must be well formed");
    assert!(
        comparison.is_complete_match(),
        "local ice P9 I3/I4 phase-control comparison: {comparison:?}"
    );
    assert_eq!(before, local_library_snapshot());
}
