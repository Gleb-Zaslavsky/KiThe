//! Real Boudouard `P,H` I1/I2/I4 fixed-topology cross-validation.
//!
//! The same independent interior reference anchors fixed-topology agreement,
//! phase-control appearance, stable absence, disappearance, and accepted-state
//! continuation without borrowing a canonical solution to construct `H`.

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::EquilibriumConstraint;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
    PhEnthalpyGrid, PhRangePointPreparation, PhRangeRequest, PhRangeSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    FixedPressureEnthalpySolution, PhSolveMode, ResolvedPhaseEnthalpyRequest, solve_resolved_ph,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
    PurePhaseProductionEvidenceRequest, activation_evidence_from_solution,
    disappearance_evidence_from_solution, stable_inactive_evidence_from_solution,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    PurePhaseBoundaryPrediction, PurePhaseBoundaryTolerances, PurePhaseBoundaryValidator,
    PurePhaseCrossValidationTolerances, cross_validate_pure_phase, evaluate_pure_phase_boundary,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::PhaseControlPolicy;
use crate::Thermodynamics::ChemEquilibrium::pure_phase_ph_validation::{
    PurePhasePhCanonicalEvidence, PurePhasePhCrossValidationTolerances, PurePhasePhProblem,
    PurePhasePhValidator, compare_pure_phase_ph_validation,
};
use crate::Thermodynamics::ChemEquilibrium::real_boudouard_ph_validation::{
    BOUDOUARD_PH_PRESSURE_PA, BOUDOUARD_PH_REFERENCE_TEMPERATURE_K,
    build_boudouard_ph_interior_reference,
};
use crate::Thermodynamics::ChemEquilibrium::real_pure_phase_fixtures::{
    RealPurePhaseFamily, RealPurePhaseInventory, ResolvedRealPurePhaseFixture,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
use crate::library_manager::with_library_manager;

fn local_repository() -> Arc<ThermoRepository> {
    ThermoData::try_default_repository()
        .expect("the bundled offline thermochemistry repository must be available")
}

fn boudouard_fixture() -> ResolvedRealPurePhaseFixture {
    RealPurePhaseFamily::BoudouardCarbon
        .resolve_offline(local_repository(), &[])
        .expect("the pinned local Boudouard fixture must resolve without NIST")
}

fn local_library_snapshot() -> Vec<(String, Vec<u8>)> {
    with_library_manager(|manager| {
        vec![
            manager.substance_base_path().to_string(),
            manager.all_keys_substance_path().to_string(),
            manager.elements_path().to_string(),
        ]
    })
    .into_iter()
    .map(|path| {
        (
            path.clone(),
            fs::read(&path).expect("local Boudouard library must be readable"),
        )
    })
    .collect()
}

fn frozen_janaf_snapshot() -> (Vec<u8>, Vec<u8>) {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src")
        .join("Thermodynamics")
        .join("ChemEquilibrium")
        .join("frozen_reference")
        .join("data")
        .join("janaf");
    (
        fs::read(directory.join("boudouard_reaction_thermodynamics.metadata.json"))
            .expect("frozen Boudouard metadata must be readable"),
        fs::read(directory.join("boudouard_reaction_thermodynamics.rows.json"))
            .expect("frozen Boudouard rows must be readable"),
    )
}

fn canonical_boudouard_evidence(
    problem: &PurePhasePhProblem,
    canonical: &FixedPressureEnthalpySolution,
    initial_moles: &[f64],
) -> PurePhasePhCanonicalEvidence {
    let canonical_moles = canonical.equilibrium().component_moles();
    let chemical_log_residual = problem
        .chemical_log_residual_for_gas_moles(&canonical_moles[..2], canonical.temperature())
        .expect("canonical Boudouard gas state must remain an independent chemical state")
        .abs();
    let initial_carbon = initial_moles[0] + initial_moles[1] + initial_moles[2];
    let initial_oxygen = initial_moles[0] + 2.0 * initial_moles[1];
    let max_abs_element_balance = [
        (canonical_moles[0] + canonical_moles[1] + canonical_moles[2] - initial_carbon).abs(),
        (canonical_moles[0] + 2.0 * canonical_moles[1] - initial_oxygen).abs(),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    PurePhasePhCanonicalEvidence::new(
        problem.case_identity(),
        canonical.temperature(),
        canonical_moles[..2].to_vec(),
        canonical_moles[2],
        canonical.calculated_enthalpy(),
        Some(chemical_log_residual),
        Some(max_abs_element_balance),
    )
    .expect("canonical Boudouard evidence must remain finite")
}

fn boudouard_comparison_tolerances() -> PurePhasePhCrossValidationTolerances {
    PurePhasePhCrossValidationTolerances {
        max_abs_temperature_delta: 1e-3,
        max_abs_mole_delta: 1e-5,
        max_abs_enthalpy_delta: 1e-2,
        max_abs_chemical_log_residual: 1e-6,
        // The canonical Boudouard solve accepts the same element residual
        // scale as its nonlinear acceptance contract. Requiring 1e-8 here
        // would reject an otherwise independently matching state.
        max_abs_element_balance: 2e-6,
        ..PurePhasePhCrossValidationTolerances::default()
    }
}

fn boudouard_ph_request(
    fixture: &ResolvedRealPurePhaseFixture,
    initial: MultiphaseInitialComposition,
    target_enthalpy: f64,
    initial_temperature: f64,
    bounds: TemperatureBounds,
) -> ResolvedPhaseEnthalpyRequest<'_> {
    ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
        fixture.resolved(),
        initial,
        EquilibriumConstraint::ph(
            BOUDOUARD_PH_PRESSURE_PA,
            BOUDOUARD_PH_PRESSURE_PA,
            target_enthalpy,
            initial_temperature,
        )
        .expect("Boudouard P,H constraint must validate"),
        bounds,
        fixture.thermochemistry().clone(),
    )
    .expect("Boudouard P,H request must validate")
}

fn solve_boudouard_ph_with_phase_control(
    fixture: &ResolvedRealPurePhaseFixture,
    initial: MultiphaseInitialComposition,
    target_enthalpy: f64,
    initial_temperature: f64,
    bounds: TemperatureBounds,
) -> FixedPressureEnthalpySolution {
    solve_resolved_ph(
        boudouard_ph_request(
            fixture,
            initial,
            target_enthalpy,
            initial_temperature,
            bounds,
        )
        // Boudouard topology transitions require the guarded nested route.
        // Auto-policy recovery is intentionally tested separately once its
        // monolithic all-active probe has a valid nested fallthrough policy.
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("Boudouard P,H production phase control must solve")
}

fn boudouard_high_temperature_target(
    fixture: &ResolvedRealPurePhaseFixture,
    inventory: &RealPurePhaseInventory,
    temperature: f64,
) -> f64 {
    let initial = fixture
        .initial_composition(inventory)
        .expect("Boudouard high-temperature inventory must match the layout");
    fixture
        .thermochemistry()
        .enthalpy_model()
        .evaluate_total(initial.moles(), temperature)
        .expect("local Boudouard high-temperature enthalpy must be finite")
}

fn solve_boudouard_ph_range(
    fixture: &ResolvedRealPurePhaseFixture,
    initial: MultiphaseInitialComposition,
    targets: PhEnthalpyGrid,
    initial_temperature: f64,
) -> PhRangeSolution {
    PhRangeRequest::from_resolved_thermochemistry(
        fixture.resolved(),
        initial,
        BOUDOUARD_PH_PRESSURE_PA,
        BOUDOUARD_PH_PRESSURE_PA,
        targets,
        TemperatureBounds::new(600.0, 1_450.0).expect("Boudouard P,H sweep bounds must validate"),
        initial_temperature,
        fixture.thermochemistry().clone(),
    )
    .expect("Boudouard P,H range request must validate")
    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve()
    .expect("Boudouard P,H range must publish every accepted target")
}

fn assert_boudouard_range_point_matches_independent_oracle(
    fixture: &ResolvedRealPurePhaseFixture,
    inventory: &RealPurePhaseInventory,
    point: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::PhRangePoint,
) {
    let solution = point.solution();
    let gas = PhaseId::new(Some("gas".to_string()));
    let graphite = PhaseId::new(Some("solid".to_string()));
    let graphite_moles = solution
        .equilibrium()
        .phase_total(&graphite)
        .unwrap_or_default();
    if graphite_moles > 1e-6 {
        // A positive graphite amount is an interior scalar problem, so I2 is
        // available as the independent oracle for the complete continuous
        // state rather than merely for its topology.
        let problem = fixture
            .to_ph_problem(
                inventory,
                BOUDOUARD_PH_PRESSURE_PA,
                BOUDOUARD_PH_PRESSURE_PA,
                point.report().target_enthalpy_joules(),
                TemperatureBounds::new(600.0, 1_450.0)
                    .expect("Boudouard independent sweep bounds must validate"),
            )
            .expect("Boudouard independent P,H sweep point must validate");
        let independent = PurePhasePhValidator::default()
            .solve(&problem)
            .expect("interior Boudouard sweep point must have an independent P,H root");
        let evidence = canonical_boudouard_evidence(
            &problem,
            solution,
            fixture
                .initial_composition(inventory)
                .expect("Boudouard sweep inventory must match the layout")
                .moles(),
        );
        let comparison = compare_pure_phase_ph_validation(
            &problem,
            &independent,
            &evidence,
            boudouard_comparison_tolerances(),
        )
        .expect("Boudouard interior sweep comparison must be well formed");
        assert!(
            comparison.is_complete_match(),
            "Boudouard interior P,H sweep point must match independent P,H: {comparison:?}"
        );
    } else {
        // At the reduced gas-only boundary a finite positive graphite extent
        // need not exist. Here the independent oracle is the driving force,
        // not an artificial interior-root requirement.
        let lifecycle = solution
            .equilibrium()
            .phase_control_report()
            .expect("Boudouard range point must retain lifecycle evidence");
        let has_deactivation = lifecycle
            .transitions
            .iter()
            .any(|transition| !transition.deactivated.is_empty());
        let request = PurePhaseProductionEvidenceRequest::new(gas, graphite);
        let evidence = if has_deactivation {
            disappearance_evidence_from_solution(solution.equilibrium(), &request)
                .expect("Boudouard range deactivation must retain boundary evidence")
        } else {
            stable_inactive_evidence_from_solution(solution.equilibrium(), &request)
                .expect("Boudouard range stable absence must retain boundary evidence")
        };
        let boundary_scenario = fixture
            .gas_scenario_from_canonical_boundary(
                inventory.gas(),
                &evidence.gas_species,
                &evidence.boundary_gas_moles,
            )
            .expect("Boudouard range boundary gas state must match the fixture");
        let boundary_problem = fixture
            .to_pt_boundary_problem(
                &boundary_scenario,
                EquilibriumConditions::new(
                    solution.temperature(),
                    BOUDOUARD_PH_PRESSURE_PA,
                    BOUDOUARD_PH_PRESSURE_PA,
                )
                .expect("Boudouard range P,H temperature must form P,T conditions"),
            )
            .expect("Boudouard range independent boundary must validate");
        let comparison = cross_validate_pure_phase(
            &boundary_problem,
            &evidence,
            PurePhaseBoundaryValidator::default(),
            PurePhaseCrossValidationTolerances::default(),
        )
        .expect("Boudouard boundary sweep comparison must be well formed");
        assert!(
            comparison.accepted,
            "Boudouard gas-only P,H sweep point must match independent boundary physics: {comparison:?}"
        );
    }
}

#[test]
fn p9_i1_i2_i4_boudouard_independent_reference_matches_fixed_topology_ph() {
    let local_before = local_library_snapshot();
    let frozen_before = frozen_janaf_snapshot();
    let fixture = boudouard_fixture();
    let reference = build_boudouard_ph_interior_reference(&fixture)
        .expect("local Boudouard reference state must be a valid interior scalar root");
    assert_eq!(
        reference.problem.conditions().pressure(),
        BOUDOUARD_PH_PRESSURE_PA
    );
    assert_eq!(
        reference.problem.conditions().reference_pressure(),
        BOUDOUARD_PH_PRESSURE_PA
    );
    let gas_only_initial = fixture
        .initial_composition(&reference.inventory)
        .expect("Boudouard gas-only reference inventory must match the layout");
    assert_eq!(gas_only_initial.moles()[2], 0.0);
    assert!(reference.state.candidate_moles > 1e-6);
    assert!(reference.state.gas_moles.iter().all(|moles| *moles > 0.0));
    assert!(reference.state.chemical_log_residual.abs() <= 1e-8);

    let independent = PurePhasePhValidator::default()
        .solve(&reference.problem)
        .expect("independent nested Boudouard P,H must recover its scalar reference state");
    assert!(
        (independent.temperature - BOUDOUARD_PH_REFERENCE_TEMPERATURE_K).abs() <= 1e-3,
        "independent P,H failed to recover the independently selected reference temperature: {independent:?}"
    );

    let resolved = fixture.resolved();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("Boudouard fixed topology layout must validate");
    let canonical_initial = MultiphaseInitialComposition::from_dense(
        &layout,
        vec![
            reference.state.gas_moles[0],
            reference.state.gas_moles[1],
            reference.state.candidate_moles,
        ],
    )
    .expect("independent Boudouard reference moles must match the canonical layout");
    let conditions = reference.problem.conditions();
    let canonical = solve_resolved_ph(
        boudouard_ph_request(
            &fixture,
            canonical_initial,
            reference.target_enthalpy(),
            BOUDOUARD_PH_REFERENCE_TEMPERATURE_K,
            TemperatureBounds::new(
                conditions.lower_temperature(),
                conditions.upper_temperature(),
            )
            .expect("Boudouard canonical P,H bounds must validate"),
        )
        .with_ph_solve_mode(PhSolveMode::NestedTemperature),
    )
    .expect("canonical fixed-topology Boudouard P,H must solve");
    let evidence =
        canonical_boudouard_evidence(&reference.problem, &canonical, gas_only_initial.moles());
    let comparison = compare_pure_phase_ph_validation(
        &reference.problem,
        &independent,
        &evidence,
        boudouard_comparison_tolerances(),
    )
    .expect("Boudouard P,H comparison must be well formed");
    assert!(
        comparison.is_complete_match(),
        "independent Boudouard P,H reference must match canonical fixed topology: {comparison:?}; independent={independent:?}; canonical_moles={:?}",
        canonical.equilibrium().component_moles(),
    );
    assert_eq!(local_before, local_library_snapshot());
    assert_eq!(frozen_before, frozen_janaf_snapshot());
}

#[test]
fn p9_i1_i2_i3_i4_boudouard_gas_only_ph_activates_graphite_and_matches_reference() {
    let local_before = local_library_snapshot();
    let frozen_before = frozen_janaf_snapshot();
    let fixture = boudouard_fixture();
    let reference = build_boudouard_ph_interior_reference(&fixture)
        .expect("independent Boudouard reference state must be available");
    let independent = PurePhasePhValidator::default()
        .solve(&reference.problem)
        .expect("independent Boudouard P,H recovery must solve");
    let resolved = fixture.resolved();
    let graphite = PhaseId::new(Some("solid".to_string()));
    let graphite_phase_index = resolved
        .layout()
        .phase_index(&graphite)
        .expect("Boudouard resolved layout must contain graphite");
    let gas_only_initial = fixture
        .initial_composition(&reference.inventory)
        .expect("gas-only Boudouard input must match the resolved layout");
    let conditions = reference.problem.conditions();
    let production = solve_boudouard_ph_with_phase_control(
        &fixture,
        gas_only_initial.clone(),
        reference.target_enthalpy(),
        BOUDOUARD_PH_REFERENCE_TEMPERATURE_K,
        TemperatureBounds::new(
            conditions.lower_temperature(),
            conditions.upper_temperature(),
        )
        .expect("Boudouard phase-control P,H bounds must validate"),
    );
    assert!(
        production
            .equilibrium()
            .phase_total(&graphite)
            .unwrap_or_default()
            > 1e-6,
        "production Boudouard P,H must publish a positive graphite phase"
    );
    let lifecycle = production
        .equilibrium()
        .phase_control_report()
        .expect("Boudouard production P,H must publish lifecycle evidence");
    assert!(
        lifecycle.transitions.iter().any(|transition| {
            transition
                .activated
                .iter()
                .any(|phase| phase.index() == graphite_phase_index)
        }),
        "gas-only Boudouard P,H must record graphite activation: {lifecycle:?}"
    );
    let activation = activation_evidence_from_solution(
        production.equilibrium(),
        &PurePhaseProductionEvidenceRequest::new(PhaseId::new(Some("gas".to_string())), graphite),
    )
    .expect("graphite activation must retain immutable thermodynamic evidence");
    assert!(
        activation.boundary_minimum_tpd.unwrap_or(f64::NAN) < 0.0,
        "graphite pre-activation TPD must be negative: {activation:?}"
    );
    let evidence =
        canonical_boudouard_evidence(&reference.problem, &production, gas_only_initial.moles());
    let comparison = compare_pure_phase_ph_validation(
        &reference.problem,
        &independent,
        &evidence,
        boudouard_comparison_tolerances(),
    )
    .expect("Boudouard production P,H comparison must be well formed");
    assert!(
        comparison.is_complete_match(),
        "gas-only Boudouard P,H lifecycle must recover the independent state: {comparison:?}"
    );
    assert_eq!(local_before, local_library_snapshot());
    assert_eq!(frozen_before, frozen_janaf_snapshot());
}

#[test]
fn p9_i1_i3_i4_boudouard_hot_gas_ph_keeps_graphite_inactive() {
    // At high temperature the declared CO/CO2 gas branch is graphite-stable
    // without any condensed inventory. I1 evaluates that boundary from local
    // G(T); I3 must publish the same inactive topology and its own TPD report.
    let local_before = local_library_snapshot();
    let frozen_before = frozen_janaf_snapshot();
    let fixture = boudouard_fixture();
    let gas_only = RealPurePhaseInventory::new(vec![1.19, 0.005], 0.0)
        .expect("hot Boudouard gas-only inventory must validate");
    let initial = fixture
        .initial_composition(&gas_only)
        .expect("hot Boudouard gas-only inventory must match the layout");
    let temperature = 1_400.0;
    let target_enthalpy = boudouard_high_temperature_target(&fixture, &gas_only, temperature);
    let canonical = solve_boudouard_ph_with_phase_control(
        &fixture,
        initial,
        target_enthalpy,
        temperature,
        TemperatureBounds::new(1_395.0, 1_405.0).expect("hot Boudouard P,H bounds must validate"),
    );
    let gas = PhaseId::new(Some("gas".to_string()));
    let graphite = PhaseId::new(Some("solid".to_string()));
    assert_eq!(canonical.equilibrium().phase_total(&graphite), Some(0.0));
    let lifecycle = canonical
        .equilibrium()
        .phase_control_report()
        .expect("hot Boudouard P,H must publish lifecycle evidence");
    assert!(
        lifecycle.transitions.is_empty(),
        "graphite-stable gas-only Boudouard state must not activate a phase: {lifecycle:?}"
    );
    let evidence = stable_inactive_evidence_from_solution(
        canonical.equilibrium(),
        &PurePhaseProductionEvidenceRequest::new(gas, graphite),
    )
    .expect("hot Boudouard P,H must expose stable-inactive graphite evidence");
    let boundary_problem = fixture
        .to_pt_boundary_problem(
            gas_only.gas(),
            EquilibriumConditions::new(
                canonical.temperature(),
                BOUDOUARD_PH_PRESSURE_PA,
                BOUDOUARD_PH_PRESSURE_PA,
            )
            .expect("hot Boudouard P,H temperature must form P,T conditions"),
        )
        .expect("hot Boudouard independent boundary must validate");
    let boundary =
        evaluate_pure_phase_boundary(&boundary_problem, PurePhaseBoundaryTolerances::default())
            .expect("hot Boudouard independent boundary must evaluate");
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
    .expect("hot Boudouard I1/I3 cross-validation must be well formed");
    assert!(
        comparison.accepted,
        "hot Boudouard P,H stable-inactive TPD must match independent chemistry: {comparison:?}"
    );
    assert_eq!(local_before, local_library_snapshot());
    assert_eq!(frozen_before, frozen_janaf_snapshot());
}

#[test]
fn p9_i1_i3_i4_boudouard_hot_ph_deactivates_initial_graphite_transactionally() {
    // Preserve the same total C/O inventory as the stable gas-only branch,
    // but start with 0.1 mol of actual graphite. The high-temperature target
    // belongs to a carbon-free boundary; phase control must prove that fact by
    // a recorded active-to-inactive transition rather than by silently
    // overwriting the condensed amount.
    let local_before = local_library_snapshot();
    let frozen_before = frozen_janaf_snapshot();
    let fixture = boudouard_fixture();
    let target_inventory = RealPurePhaseInventory::new(vec![1.19, 0.005], 0.0)
        .expect("hot Boudouard boundary inventory must validate");
    let active_inventory = RealPurePhaseInventory::new(vec![0.99, 0.105], 0.1)
        .expect("initially active Boudouard inventory must validate");
    let initial = fixture
        .initial_composition(&active_inventory)
        .expect("initially active Boudouard inventory must match the layout");
    let temperature = 1_400.0;
    let target_enthalpy =
        boudouard_high_temperature_target(&fixture, &target_inventory, temperature);
    let canonical = solve_boudouard_ph_with_phase_control(
        &fixture,
        initial,
        target_enthalpy,
        temperature,
        TemperatureBounds::new(1_395.0, 1_405.0).expect("hot Boudouard P,H bounds must validate"),
    );
    let gas = PhaseId::new(Some("gas".to_string()));
    let graphite = PhaseId::new(Some("solid".to_string()));
    assert_eq!(canonical.equilibrium().phase_total(&graphite), Some(0.0));
    let lifecycle = canonical
        .equilibrium()
        .phase_control_report()
        .expect("Boudouard P,H disappearance must retain lifecycle evidence");
    assert_eq!(
        lifecycle.transitions.len(),
        1,
        "one initially active graphite phase must yield exactly one deactivation: {lifecycle:?}"
    );
    assert!(lifecycle.transitions[0].activated.is_empty());
    assert_eq!(lifecycle.transitions[0].deactivated.len(), 1);
    let evidence = disappearance_evidence_from_solution(
        canonical.equilibrium(),
        &PurePhaseProductionEvidenceRequest::new(gas, graphite),
    )
    .expect("Boudouard P,H deactivation must expose immutable boundary evidence");
    let boundary_scenario = fixture
        .gas_scenario_from_canonical_boundary(
            target_inventory.gas(),
            &evidence.gas_species,
            &evidence.boundary_gas_moles,
        )
        .expect("recorded Boudouard boundary gas state must match the fixture");
    let boundary_problem = fixture
        .to_pt_boundary_problem(
            &boundary_scenario,
            EquilibriumConditions::new(
                canonical.temperature(),
                BOUDOUARD_PH_PRESSURE_PA,
                BOUDOUARD_PH_PRESSURE_PA,
            )
            .expect("Boudouard P,H temperature must form P,T conditions"),
        )
        .expect("Boudouard independent boundary must validate");
    let boundary =
        evaluate_pure_phase_boundary(&boundary_problem, PurePhaseBoundaryTolerances::default())
            .expect("Boudouard independent boundary must evaluate");
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
    .expect("Boudouard I1/I3/I4 disappearance comparison must be well formed");
    assert!(
        comparison.accepted,
        "Boudouard P,H deactivation TPD must match independent chemistry: {comparison:?}"
    );
    assert_eq!(local_before, local_library_snapshot());
    assert_eq!(frozen_before, frozen_janaf_snapshot());
}

#[test]
fn p9_i1_i2_i3_i4_boudouard_ph_enthalpy_sweep_uses_accepted_continuation_both_ways() {
    // This is intentionally a small physical story rather than a performance
    // grid. It crosses the graphite boundary in both enthalpy directions and
    // proves that every continued point was seeded by a previously accepted
    // production state. Interior points use I2; boundary points use I1 TPD.
    let local_before = local_library_snapshot();
    let frozen_before = frozen_janaf_snapshot();
    let fixture = boudouard_fixture();
    let inventory = RealPurePhaseInventory::new(vec![1.19, 0.005], 0.0)
        .expect("Boudouard sweep gas-only inventory must validate");
    let initial = fixture
        .initial_composition(&inventory)
        .expect("Boudouard sweep inventory must match the layout");
    let mut target_pairs = [650.0, 700.0, 900.0, 1_100.0, 1_400.0]
        .into_iter()
        .map(|temperature| {
            (
                boudouard_high_temperature_target(&fixture, &inventory, temperature),
                temperature,
            )
        })
        .collect::<Vec<_>>();
    target_pairs.sort_by(|left, right| left.0.total_cmp(&right.0));
    let ascending_targets =
        PhEnthalpyGrid::new(target_pairs.iter().map(|(target, _)| *target).collect())
            .expect("Boudouard enthalpy targets must be strictly ascending");
    let ascending = solve_boudouard_ph_range(
        &fixture,
        initial.clone(),
        ascending_targets,
        target_pairs[0].1,
    );
    assert_eq!(ascending.points().len(), target_pairs.len());
    assert_eq!(
        ascending.report().continuation_points(),
        target_pairs.len() - 1
    );
    assert!(
        ascending.report().phase_control_transitions() > 0,
        "Boudouard ascending P,H sweep must cross a graphite lifecycle boundary"
    );
    let graphite = PhaseId::new(Some("solid".to_string()));
    assert!(ascending.points().iter().any(|point| {
        point
            .solution()
            .equilibrium()
            .phase_total(&graphite)
            .unwrap_or_default()
            > 1e-6
    }));
    assert!(
        ascending
            .points()
            .iter()
            .any(|point| { point.solution().equilibrium().phase_total(&graphite) == Some(0.0) })
    );
    for (index, point) in ascending.points().iter().enumerate() {
        assert_eq!(
            point.report().preparation(),
            if index == 0 {
                PhRangePointPreparation::Initial
            } else {
                PhRangePointPreparation::Continued
            }
        );
        if index > 0 {
            assert_eq!(
                point.report().seed_temperature(),
                ascending.points()[index - 1].report().solved_temperature(),
                "continued Boudouard point {index} must seed from the last accepted temperature"
            );
        }
        assert_boudouard_range_point_matches_independent_oracle(&fixture, &inventory, point);
    }

    let descending_pairs = target_pairs.iter().rev().copied().collect::<Vec<_>>();
    let descending = solve_boudouard_ph_range(
        &fixture,
        initial,
        PhEnthalpyGrid::new(descending_pairs.iter().map(|(target, _)| *target).collect())
            .expect("Boudouard enthalpy targets must be strictly descending"),
        descending_pairs[0].1,
    );
    assert_eq!(descending.points().len(), descending_pairs.len());
    assert_eq!(
        descending.report().continuation_points(),
        descending_pairs.len() - 1
    );
    assert!(
        descending.report().phase_control_transitions() > 0,
        "Boudouard descending P,H sweep must cross a graphite lifecycle boundary"
    );
    for (index, point) in descending.points().iter().enumerate() {
        assert_eq!(
            point.report().preparation(),
            if index == 0 {
                PhRangePointPreparation::Initial
            } else {
                PhRangePointPreparation::Continued
            }
        );
        if index > 0 {
            assert_eq!(
                point.report().seed_temperature(),
                descending.points()[index - 1].report().solved_temperature(),
                "reverse continued Boudouard point {index} must use the last accepted temperature"
            );
        }
    }
    // The two endpoints are well outside the hysteresis band. Their topology
    // therefore must be history-independent even though an eventual natural
    // in-band case is deliberately left for a separate evidence story.
    let ascending_cold_graphite = ascending.points()[0]
        .solution()
        .equilibrium()
        .phase_total(&graphite)
        .unwrap_or_default();
    let descending_cold_graphite = descending
        .points()
        .last()
        .unwrap()
        .solution()
        .equilibrium()
        .phase_total(&graphite)
        .unwrap_or_default();
    let ascending_hot_graphite = ascending
        .points()
        .last()
        .unwrap()
        .solution()
        .equilibrium()
        .phase_total(&graphite)
        .unwrap_or_default();
    let descending_hot_graphite = descending.points()[0]
        .solution()
        .equilibrium()
        .phase_total(&graphite)
        .unwrap_or_default();
    assert!(
        (ascending_cold_graphite - descending_cold_graphite).abs() <= 1e-5,
        "cold endpoint must be history-independent outside the hysteresis band: ascending={ascending_cold_graphite:e}, descending={descending_cold_graphite:e}"
    );
    assert!(
        (ascending_hot_graphite - descending_hot_graphite).abs() <= 1e-12,
        "hot endpoint must be carbon-free in both histories: ascending={ascending_hot_graphite:e}, descending={descending_hot_graphite:e}"
    );
    assert_eq!(local_before, local_library_snapshot());
    assert_eq!(frozen_before, frozen_janaf_snapshot());
}

#[test]
#[ignore = "release-oriented real Boudouard P,H lifecycle diagnostic"]
fn p9_i1_i2_i3_i4_boudouard_ph_lifecycle_diagnostic() {
    // Operator-facing characterization only. The strict assertions live in
    // the five focused stories above; this report makes the evidence chain
    // inspectable without dumping nonlinear iteration traces.
    let local_before = local_library_snapshot();
    let frozen_before = frozen_janaf_snapshot();
    let fixture = boudouard_fixture();
    let reference = build_boudouard_ph_interior_reference(&fixture)
        .expect("independent Boudouard interior reference must be available");
    let independent = PurePhasePhValidator::default()
        .solve(&reference.problem)
        .expect("independent Boudouard P,H reference must solve");
    let gas_only_initial = fixture
        .initial_composition(&reference.inventory)
        .expect("Boudouard gas-only reference inventory must match the layout");
    let conditions = reference.problem.conditions();
    let appearance = solve_boudouard_ph_with_phase_control(
        &fixture,
        gas_only_initial.clone(),
        reference.target_enthalpy(),
        BOUDOUARD_PH_REFERENCE_TEMPERATURE_K,
        TemperatureBounds::new(
            conditions.lower_temperature(),
            conditions.upper_temperature(),
        )
        .expect("Boudouard appearance diagnostic bounds must validate"),
    );

    let hot_inventory = RealPurePhaseInventory::new(vec![1.19, 0.005], 0.0)
        .expect("Boudouard hot gas-only inventory must validate");
    let hot_target = boudouard_high_temperature_target(&fixture, &hot_inventory, 1_400.0);
    let stable_absence = solve_boudouard_ph_with_phase_control(
        &fixture,
        fixture
            .initial_composition(&hot_inventory)
            .expect("Boudouard hot inventory must match the layout"),
        hot_target,
        1_400.0,
        TemperatureBounds::new(1_395.0, 1_405.0)
            .expect("Boudouard stable-absence diagnostic bounds must validate"),
    );
    let active_inventory = RealPurePhaseInventory::new(vec![0.99, 0.105], 0.1)
        .expect("initially active Boudouard inventory must validate");
    let disappearance = solve_boudouard_ph_with_phase_control(
        &fixture,
        fixture
            .initial_composition(&active_inventory)
            .expect("active Boudouard inventory must match the layout"),
        hot_target,
        1_400.0,
        TemperatureBounds::new(1_395.0, 1_405.0)
            .expect("Boudouard disappearance diagnostic bounds must validate"),
    );
    let gas = PhaseId::new(Some("gas".to_string()));
    let graphite = PhaseId::new(Some("solid".to_string()));
    let request = PurePhaseProductionEvidenceRequest::new(gas, graphite.clone());
    println!("Boudouard P,H lifecycle diagnostic | P=p0={BOUDOUARD_PH_PRESSURE_PA:.0} Pa");
    println!(
        "route              | H target J      | T K       | CO mol       | CO2 mol      | C(gr) mol    | transitions | boundary TPD"
    );
    println!(
        "independent I2     | {:+.6e} | {:9.4} | {:+.6e} | {:+.6e} | {:+.6e} | -           | chemical={:+.3e}",
        reference.target_enthalpy(),
        independent.temperature,
        independent.gas_moles[0],
        independent.gas_moles[1],
        independent.candidate_moles,
        independent.chemical_log_residual,
    );
    for (route, solution) in [
        ("appearance I3", &appearance),
        ("stable absence", &stable_absence),
        ("disappearance", &disappearance),
    ] {
        let moles = solution.equilibrium().component_moles();
        let transitions = solution
            .equilibrium()
            .phase_control_report()
            .map(|report| report.transitions.len())
            .unwrap_or_default();
        let tpd = if solution
            .equilibrium()
            .phase_total(&graphite)
            .unwrap_or_default()
            > 1e-6
        {
            activation_evidence_from_solution(solution.equilibrium(), &request)
                .ok()
                .and_then(|evidence| evidence.boundary_minimum_tpd)
        } else {
            disappearance_evidence_from_solution(solution.equilibrium(), &request)
                .or_else(|_| {
                    stable_inactive_evidence_from_solution(solution.equilibrium(), &request)
                })
                .ok()
                .and_then(|evidence| evidence.boundary_minimum_tpd)
        };
        println!(
            "{route:18} | {:+.6e} | {:9.4} | {:+.6e} | {:+.6e} | {:+.6e} | {transitions:11} | {:+.6e}",
            solution.target_enthalpy(),
            solution.temperature(),
            moles[0],
            moles[1],
            moles[2],
            tpd.unwrap_or(f64::NAN),
        );
    }

    let target_pairs = [650.0, 700.0, 900.0, 1_100.0, 1_400.0]
        .into_iter()
        .map(|temperature| {
            (
                boudouard_high_temperature_target(&fixture, &reference.inventory, temperature),
                temperature,
            )
        })
        .collect::<Vec<_>>();
    let mut target_pairs = target_pairs;
    target_pairs.sort_by(|left, right| left.0.total_cmp(&right.0));
    let sweep = solve_boudouard_ph_range(
        &fixture,
        gas_only_initial,
        PhEnthalpyGrid::new(target_pairs.iter().map(|(target, _)| *target).collect())
            .expect("Boudouard diagnostic targets must be monotone"),
        target_pairs[0].1,
    );
    println!(
        "sweep index | H target J      | seed K    | solved K  | C(gr) mol    | final transitions | trial phase events | preparation"
    );
    for point in sweep.points() {
        let final_transitions = point
            .solution()
            .equilibrium()
            .phase_control_report()
            .map(|report| report.transitions.len())
            .unwrap_or_default();
        println!(
            "{:11} | {:+.6e} | {:9.4} | {:9.4} | {:+.6e} | {:17} | {:18} | {:?}",
            point.report().index(),
            point.report().target_enthalpy_joules(),
            point.report().seed_temperature(),
            point.report().solved_temperature(),
            point
                .solution()
                .equilibrium()
                .phase_total(&graphite)
                .unwrap_or_default(),
            final_transitions,
            point.report().phase_control_transitions(),
            point.report().preparation(),
        );
    }
    assert_eq!(local_before, local_library_snapshot());
    assert_eq!(frozen_before, frozen_janaf_snapshot());
}
