//! Live-data regression tests for the phase-equilibrium boundary.
//!
//! These tests deliberately use the bundled thermochemistry repository instead
//! of synthetic coefficient fixtures. They check the contracts that matter
//! when the phase subsystem and the equilibrium solver meet in production:
//!
//! - a real local repository resolves named phase systems with explicit lookup
//!   policy;
//! - the resolved phase system keeps provenance and canonical order stable;
//! - the phase-equilibrium pipeline can solve those real records end to end;
//! - solve attempts do not mutate the resolved source data; and
//! - the canonical JSON thermochemistry files remain byte-for-byte unchanged.
//! - bounded phase control remains usable on real local thermochemistry.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::sync::Arc;

use tabled::settings::Style;
use tabled::{Table, Tabled};

use crate::library_manager::with_library_manager;
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EnthalpyScale, EquilibriumConstraint, TemperatureBounds, TotalEnthalpyJoules,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, ReactionExtentErrorKind,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_formulation::PreparedPhFormulation;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
    PhEnthalpyGrid, PhRangeDirection, PhRangeError, PhRangePointPreparation, PhRangeRequest,
    PhRangeSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    solve_resolved_ph, FixedPressureEnthalpySolution, PhSolveMode, PhSolvePath,
    PhTemperatureSolveReport, PhTrialPreparation, ResolvedPhaseEnthalpyRequest,
    ResolvedThermochemistry,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_prepared_runner::PreparedEquilibriumRunner;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, PreparedEquilibriumProblem, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    prepare_baked_rst_symbolic_problem_for_test, prepare_rst_symbolic_problem_from_prepared,
    RustedSciTheSolver,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptFailureKind, SolverAttemptMetrics, SolverAttemptOutcome, SolverAttemptReport,
    SolverBackend, SolverPolicy, SolverTermination,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    postprocess_temperature_range_solution, TemperatureInterpolationPolicy,
    TemperatureInterpolationSpace, TemperaturePostprocessingPolicy, TemperatureResamplingGrid,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
    TemperatureGrid, TemperatureRangeDirection, TemperatureRangePointPreparation,
    TemperatureRangeRequest,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode;
#[allow(deprecated)]
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver_for_T_range;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{InitialPhaseSet, PhaseStatus};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    build_phase_equilibrium_problem_with_timing, PhaseEquilibriumBuildRequest,
    SupportedPhaseModelPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    solve_resolved_pt, EquilibriumSolveOptions, PhaseControlPolicy, PhaseEquilibriumPipelineError,
    PhaseEquilibriumPipelineRequest, ResolvedPhaseEquilibriumRequest,
};
use crate::Thermodynamics::User_PhaseOrSolution::{
    ResolvedPhaseSystem, SubstanceSystemFactory, SubstanceSystemSpec, SubstanceSystemSpecBuilder,
    SubstancesContainer,
};

/// Stable FNV-1a fingerprint for a canonical JSON record.
///
/// This is deliberately a fixture-integrity check, not a cryptographic
/// primitive. A changed value forces a review of the physical scenario rather
/// than silently changing a live regression's inputs.
fn stable_json_fingerprint(value: &serde_json::Value) -> u64 {
    serde_json::to_vec(value)
        .expect("thermochemistry record must serialize canonically")
        .into_iter()
        .fold(0xcbf2_9ce4_8422_2325_u64, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
        })
}

/// Captures the canonical local thermochemistry files used by ordinary lookup.
///
/// Live equilibrium tests must only read these files. Keeping the byte-level
/// snapshot at the test boundary catches an accidental write path even when a
/// resolved in-memory repository remains internally immutable.
fn live_library_file_snapshot() -> Vec<(String, Vec<u8>)> {
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
            let contents = fs::read(&path)
                .unwrap_or_else(|error| panic!("must read live library file '{path}': {error}"));
            (path, contents)
        })
        .collect()
}

fn live_repository() -> Arc<crate::Thermodynamics::thermo_lib_api::ThermoRepository> {
    ThermoData::try_default_repository()
        .expect("bundled thermochemistry repository must be available for live tests")
}

#[test]
fn live_catalog_consistency_report_is_deterministic_and_non_mutating() {
    let before = live_library_file_snapshot();
    let repository = live_repository();
    let report = repository.consistency_report();
    let repeated = repository.consistency_report();
    let after = live_library_file_snapshot();

    assert_eq!(report, repeated);
    assert!(report.indexed_pair_count() > 0);
    assert!(report.payload_pair_count() > 0);
    assert_eq!(
        before, after,
        "catalog diagnostics must not write JSON files"
    );

    println!(
        "catalog consistency: indexed={} unique={} payload={} duplicates={} missing_payload={} orphan_payload={} consistent={}",
        report.indexed_pair_count(),
        report.unique_indexed_pair_count(),
        report.payload_pair_count(),
        report.duplicate_index_pairs().len(),
        report.indexed_without_payload().len(),
        report.payload_without_index().len(),
        report.is_consistent(),
    );
}

#[test]
fn live_exact_element_search_excludes_unrequested_elements() {
    let mut catalog = ThermoData::try_new().expect("bundled catalog must load");
    let requested = HashSet::from(["H".to_string(), "O".to_string()]);
    let found = catalog.search_by_exact_elements(vec!["O".into(), "H".into(), "H".into()]);

    assert!(
        !found.is_empty(),
        "the local catalog must contain H/O records"
    );
    for substance in found {
        let actual: HashSet<String> = catalog
            .get_elements_data()
            .iter()
            .filter(|(_, entries)| entries.iter().any(|pair| pair.first() == Some(&substance)))
            .map(|(element, _)| element.clone())
            .collect();
        assert_eq!(actual, requested, "exact element search leaked {substance}");
    }
}

#[test]
fn live_c_h_o_element_modes_distinguish_allowed_alphabet_from_exact_set() {
    let mut catalog = ThermoData::try_new().expect("bundled catalog must load");
    let mut subset = catalog
        .search_by_elements_only(vec!["C".into(), "H".into(), "O".into()])
        .into_iter()
        .filter(|name| {
            live_repository()
                .LibThermoData
                .get("NASA_gas")
                .is_some_and(|records| records.contains_key(name))
        })
        .collect::<Vec<_>>();
    subset.sort();
    subset.dedup();

    let mut exact = catalog
        .search_by_exact_elements(vec!["C".into(), "H".into(), "O".into()])
        .into_iter()
        .filter(|name| {
            live_repository()
                .LibThermoData
                .get("NASA_gas")
                .is_some_and(|records| records.contains_key(name))
        })
        .collect::<Vec<_>>();
    exact.sort();
    exact.dedup();

    assert!(
        subset.len() >= 50,
        "the local NASA gas catalog must provide at least 50 records using only C/H/O, got {}",
        subset.len()
    );
    assert!(
        exact.len() < subset.len(),
        "strict all-three-element matching must be narrower than the allowed-element search"
    );
    assert!(subset.contains(&"C2H4".to_string()));
    assert!(!exact.contains(&"C2H4".to_string()));
}

fn live_record_fingerprint(library: &str, substance: &str) -> u64 {
    let repository = live_repository();
    let record = repository
        .LibThermoData
        .get(library)
        .and_then(|records| records.get(substance))
        .unwrap_or_else(|| panic!("live fixture requires {library}::{substance}"));
    stable_json_fingerprint(record)
}

fn live_multiphase_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
        ("gas".to_string(), vec!["H2O".to_string(), "O2".to_string()]),
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
    .expect("live multi-phase spec must remain valid against the bundled catalog")
}

fn live_gas_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
        "N2".to_string(),
        "O2".to_string(),
    ]))
    .with_library_priorities(vec!["NASA_gas".to_string()])
    .with_search_in_nist(false)
    .build()
    .expect("live gas spec must remain valid against the bundled catalog")
}

/// The historical O2/O dissociation fixture, expressed through the typed
/// resolved-phase input boundary and the bundled NASA gas records.
fn live_oxygen_dissociation_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
        "O2".to_string(),
        "O".to_string(),
    ]))
    .with_library_priorities(vec!["NASA_gas".to_string()])
    .with_search_in_nist(false)
    .build()
    .expect("live O2/O spec must remain valid against the bundled catalog")
}

fn live_gas_solid_water_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
        ("gas".to_string(), vec!["H2O".to_string(), "O2".to_string()]),
        ("solid".to_string(), vec!["H2O(s)".to_string()]),
    ])))
    .with_phase_natures(Some(HashMap::from([
        (
            "gas".to_string(),
            crate::Thermodynamics::User_substances::Phases::Gas,
        ),
        (
            "solid".to_string(),
            crate::Thermodynamics::User_substances::Phases::Solid,
        ),
    ])))
    .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
    .with_search_in_nist(false)
    .build()
    .expect("live gas/solid water spec must remain valid against the bundled catalog")
}

fn live_gas_solid_carbon_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
        ("gas".to_string(), vec!["CO".to_string(), "CO2".to_string()]),
        ("solid".to_string(), vec!["C(gr)".to_string()]),
    ])))
    .with_phase_natures(Some(HashMap::from([
        (
            "gas".to_string(),
            crate::Thermodynamics::User_substances::Phases::Gas,
        ),
        (
            "solid".to_string(),
            crate::Thermodynamics::User_substances::Phases::Solid,
        ),
    ])))
    .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
    .with_search_in_nist(false)
    .build()
    .expect("live gas/graphite spec must remain valid against the bundled catalog")
}

/// A real local reactive gas system with one H/O reaction degree of freedom.
fn live_reactive_gas_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
        "H2".to_string(),
        "O2".to_string(),
        "H2O".to_string(),
    ]))
    .with_library_priorities(vec!["NASA_gas".to_string()])
    .with_search_in_nist(false)
    .build()
    .expect("live reactive gas spec must remain valid against the bundled catalog")
}

fn live_resolved_system() -> ResolvedPhaseSystem {
    let repository = live_repository();
    let spec = live_multiphase_spec();
    SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)
        .expect("live multi-phase spec must resolve against the bundled repository")
}

fn live_gas_pipeline_request() -> PhaseEquilibriumPipelineRequest {
    PhaseEquilibriumPipelineRequest::new(
        live_gas_spec(),
        vec![0.79, 0.21],
        EquilibriumConditions::new(500.0, 101_325.0, 101_325.0)
            .expect("live temperature and pressure must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
}

fn live_reactive_gas_pipeline_request() -> PhaseEquilibriumPipelineRequest {
    live_reactive_gas_pipeline_request_with_scale(1.0)
}

fn live_reactive_gas_pipeline_request_with_scale(scale: f64) -> PhaseEquilibriumPipelineRequest {
    assert!(scale.is_finite() && scale > 0.0);
    PhaseEquilibriumPipelineRequest::new(
        live_reactive_gas_spec(),
        // This warm start has exactly the same H/O inventory as [2, 1, 0],
        // while avoiding an unnecessarily stiff 1e-30 seed for the dominant
        // product in a strongly reacting live-thermochemistry fixture.
        vec![0.1 * scale, 0.05 * scale, 1.9 * scale],
        EquilibriumConditions::new(2_500.0, 101_325.0, 101_325.0)
            .expect("live reactive gas conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(
        EquilibriumSolveOptions::default()
            .with_keq_validation_mode(EquilibriumConstantValidationMode::WhenApplicable),
    )
}

#[test]
#[ignore = "release-oriented real-data P,H story"]
fn live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas must resolve for the P,H story");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("resolved phase layout must remain valid");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])
        .expect("reactive initial composition must match the layout");
    let pressure = 101_325.0;
    let reference_pressure = 101_325.0;
    let known_temperature = 2_500.0;
    let known_conditions =
        EquilibriumConditions::new(known_temperature, pressure, reference_pressure).unwrap();
    let pt = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(&resolved, known_conditions, initial.clone())
            .with_solve_options(
                EquilibriumSolveOptions::new()
                    .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
                    .expect("legacy NR reference policy must validate"),
            ),
    )
    .expect("the reference P,T solve must be accepted");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live NASA records must provide a complete thermochemistry bundle");
    let monolithic_thermochemistry = thermochemistry.clone();
    let symbolic_monolithic_thermochemistry = thermochemistry.clone();
    let auto_thermochemistry = thermochemistry.clone();
    assert_eq!(thermochemistry.len(), 3);
    assert!(thermochemistry
        .provenance()
        .iter()
        .all(|row| row.library() == "NASA_gas"));
    let enthalpy = thermochemistry.enthalpy_model();
    let target_enthalpy = enthalpy
        .evaluate_total(pt.component_moles(), known_temperature)
        .expect("accepted P,T composition must have finite total enthalpy");

    let constraint = EquilibriumConstraint::ph_joules(
        pressure,
        reference_pressure,
        TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
        1_900.0,
    )
    .unwrap();
    let nested_started = std::time::Instant::now();
    let ph = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            constraint,
            TemperatureBounds::new(1_900.0, 2_900.0).unwrap(),
            thermochemistry,
        )
        .expect("live P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        ),
    )
    .expect("live P,H outer solve must recover the reference state");
    let report_bundle = ph
        .thermochemistry()
        .expect("P,H result must retain thermochemistry provenance");
    assert_eq!(report_bundle.len(), 3);
    assert!(report_bundle.temperature_bounds().lower() <= 1_900.0);
    assert!(report_bundle
        .temperature_bounds()
        .contains(ph.temperature()));

    let max_mole_delta = pt
        .component_moles()
        .iter()
        .zip(ph.equilibrium().component_moles())
        .map(|(expected, recovered)| (expected - recovered).abs())
        .fold(0.0, f64::max);
    println!(
        "live P,H inverse story | reference T={known_temperature:.6} K | recovered T={:.6} K | target H={target_enthalpy:.6e} J | scaled H error={:.3e} | trials={} | max mole delta={max_mole_delta:.3e}",
        ph.temperature(),
        ph.enthalpy_error() / ph.report().enthalpy_scale_joules(),
        ph.report().trials().len(),
    );
    assert!(
        (ph.temperature() - known_temperature).abs() < 1.0e-5,
        "recovered temperature was {} K",
        ph.temperature()
    );
    assert_eq!(ph.report().solve_path(), PhSolvePath::NestedTemperature);
    assert!((ph.enthalpy_error() / ph.report().enthalpy_scale_joules()).abs() < 1.0e-8);
    for (expected, recovered) in pt
        .component_moles()
        .iter()
        .zip(ph.equilibrium().component_moles())
    {
        assert!((expected - recovered).abs() < 1.0e-7);
    }
    assert!(ph.report().trials().len() >= 2);
    assert!(
        ph.report()
            .trials()
            .iter()
            .all(|trial| trial.inner_backend_attempts() > 0),
        "every accepted temperature trial must retain inner backend evidence"
    );
    assert!(
        ph.report()
            .trials()
            .iter()
            .all(|trial| trial.inner_evidence().is_some()),
        "every real nested P,H trial must retain its full immutable inner trace"
    );
    assert!(ph.report().trials().iter().all(|trial| {
        trial
            .inner_evidence()
            .is_some_and(|evidence| evidence.solve_report().accepted_attempt().is_some())
    }));
    assert_eq!(
        ph.report().inner_backend_attempts(),
        ph.report()
            .trials()
            .iter()
            .map(|trial| trial.inner_backend_attempts())
            .sum::<usize>()
    );
    assert!(ph.report().max_temperature_evaluations() >= ph.report().trials().len());
    assert_eq!(
        ph.report().fixed_formulation_builds(),
        1,
        "fixed P,H trials must share one prepared structural formulation"
    );
    assert_eq!(
        ph.report().fixed_formulation_reuses(),
        ph.report().trials().len().saturating_sub(1),
        "every later accepted P,H trial must retarget the shared formulation"
    );
    assert_eq!(
        ph.report().trials()[0].preparation(),
        PhTrialPreparation::FixedFormulationInitial
    );
    assert!(ph.report().trials()[1..]
        .iter()
        .all(|trial| trial.preparation() == PhTrialPreparation::FixedFormulationReused));
    assert!(ph.report().timing().enabled());
    assert!(ph.report().timing().total() > std::time::Duration::ZERO);
    assert!(ph.report().timing().scalar_orchestration() > std::time::Duration::ZERO);
    assert!(ph.report().timing().enthalpy_evaluation() > std::time::Duration::ZERO);
    assert!(ph.report().trials().iter().all(|trial| {
        let timing = trial.timing();
        timing.enabled()
            && timing.total() >= timing.inner_equilibrium()
            && timing.total() >= timing.enthalpy_evaluation()
            && trial.inner_evidence().is_some()
    }));

    // The monolithic system is intentionally compared with the independently
    // safeguarded nested path, not merely with its own residual report.
    let monolithic_started = std::time::Instant::now();
    let monolithic = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            MultiphaseInitialComposition::from_dense(&layout, pt.component_moles().to_vec())
                .expect("monolithic inverse seed must match the reference layout"),
            EquilibriumConstraint::ph_joules(
                pressure,
                reference_pressure,
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                known_temperature,
            )
            .unwrap(),
            TemperatureBounds::new(1_900.0, 2_900.0).unwrap(),
            monolithic_thermochemistry,
        )
        .expect("monolithic live P,H request must validate")
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        ),
    )
    .expect("monolithic live P,H solve must recover the reference state");
    assert!(
        monolithic.report().trials().is_empty(),
        "monolithic P,H must not masquerade as an outer scalar trial"
    );
    assert_eq!(
        monolithic.report().solve_path(),
        PhSolvePath::MonolithicFixedActiveSet
    );
    assert_eq!(monolithic.report().fixed_formulation_builds(), 1);
    assert_eq!(monolithic.report().fixed_formulation_reuses(), 0);
    assert!(
        monolithic
            .report()
            .monolithic_evidence()
            .is_some_and(|evidence| evidence.solve_report().accepted_attempt().is_some()),
        "monolithic P,H must publish the coupled backend trace with its result"
    );
    assert!(
        (monolithic.temperature() - known_temperature).abs() < 1.0e-4,
        "monolithic P,H recovered {} K instead of {known_temperature} K",
        monolithic.temperature()
    );
    assert!(
        (monolithic.enthalpy_error() / monolithic.report().enthalpy_scale_joules()).abs() < 1.0e-8
    );
    for (nested_moles, monolithic_moles) in ph
        .equilibrium()
        .component_moles()
        .iter()
        .zip(monolithic.equilibrium().component_moles())
    {
        assert!((nested_moles - monolithic_moles).abs() < 1.0e-6);
    }
    println!(
        "live monolithic P,H inverse | T={:.6} K | scaled H error={:.3e} | backend={} ",
        monolithic.temperature(),
        monolithic.enthalpy_error() / monolithic.report().enthalpy_scale_joules(),
        monolithic.equilibrium().solve_report().summary(),
    );

    // The selected 1900..2900 K domain is inside one NASA-gas polynomial
    // interval for this fixture. This is the exact real-data counterpart to
    // the synthetic RST symbolic P,H regression: RST owns the full Jacobian
    // for `[ln(n), theta_T]`, while the accepted state is checked against the
    // independently solved P,T reference above.
    let symbolic_monolithic_started = std::time::Instant::now();
    let symbolic_monolithic = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            MultiphaseInitialComposition::from_dense(&layout, pt.component_moles().to_vec())
                .expect("symbolic monolithic seed must match the reference layout"),
            EquilibriumConstraint::ph_joules(
                pressure,
                reference_pressure,
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                known_temperature,
            )
            .unwrap(),
            TemperatureBounds::new(1_900.0, 2_900.0).unwrap(),
            symbolic_monolithic_thermochemistry,
        )
        .expect("symbolic monolithic live P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::Monolithic)
        .with_solve_options(
            EquilibriumSolveOptions::new()
                .with_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
                    RustedSciTheSolver::LevenbergMarquardt,
                )))
                .expect("RST symbolic P,H policy must validate")
                .with_timing_mode(EquilibriumTimingMode::Enabled),
        ),
    )
    .expect("RST symbolic monolithic P,H must recover the reference state");
    assert_eq!(
        symbolic_monolithic.report().solve_path(),
        PhSolvePath::MonolithicFixedActiveSet
    );
    assert!(
        (symbolic_monolithic.temperature() - known_temperature).abs() < 1.0e-4,
        "symbolic monolithic P,H recovered {} K instead of {known_temperature} K",
        symbolic_monolithic.temperature()
    );
    assert!(
        (symbolic_monolithic.enthalpy_error()
            / symbolic_monolithic.report().enthalpy_scale_joules())
        .abs()
            < 1.0e-8
    );
    for (expected, recovered) in monolithic
        .equilibrium()
        .component_moles()
        .iter()
        .zip(symbolic_monolithic.equilibrium().component_moles())
    {
        assert!((expected - recovered).abs() < 1.0e-6);
    }

    let auto_started = std::time::Instant::now();
    let auto = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            MultiphaseInitialComposition::from_dense(&layout, pt.component_moles().to_vec())
                .expect("auto inverse seed must match the reference layout"),
            EquilibriumConstraint::ph_joules(
                pressure,
                reference_pressure,
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                known_temperature,
            )
            .unwrap(),
            TemperatureBounds::new(1_900.0, 2_900.0).unwrap(),
            auto_thermochemistry,
        )
        .expect("auto live P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::Auto),
    )
    .expect("auto live P,H solve must accept the monolithic candidate");
    assert_eq!(
        auto.report().solve_path(),
        PhSolvePath::MonolithicFixedActiveSet
    );
    assert!(auto.report().fallback_reason().is_none());
    assert!((auto.temperature() - known_temperature).abs() < 1.0e-4);
    assert!((auto.enthalpy_error() / auto.report().enthalpy_scale_joules()).abs() < 1.0e-8);
    for (expected, recovered) in monolithic
        .equilibrium()
        .component_moles()
        .iter()
        .zip(auto.equilibrium().component_moles())
    {
        assert!((expected - recovered).abs() < 1.0e-7);
    }

    println!(
        "live real P,H formulation characterization\n{}",
        Table::new(vec![
            live_ph_path_row("nested", &ph, nested_started.elapsed()),
            live_ph_path_row("monolithic", &monolithic, monolithic_started.elapsed()),
            live_ph_path_row(
                "monolithic-rst-lm",
                &symbolic_monolithic,
                symbolic_monolithic_started.elapsed(),
            ),
            live_ph_path_row("auto", &auto, auto_started.elapsed()),
        ])
        .with(Style::rounded())
    );
    assert_eq!(before, live_library_file_snapshot());
}

/// Symbolic P,H deliberately stays exact with respect to the native NASA/NIST
/// polynomial selected for every component. A range that crosses a coefficient
/// switch remains solvable through the analytic route, but must not be handed
/// to RST as one fictitiously global symbolic polynomial.
#[test]
#[ignore = "release-oriented real-data symbolic P,H interval guard"]
fn live_symbolic_ph_rejects_a_native_polynomial_boundary_without_mutating_json() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas must resolve for the symbolic P,H boundary guard");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("resolved phase layout must remain valid");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])
        .expect("reactive initial composition must match the layout");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live NASA records must provide a thermochemistry bundle");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(initial.moles(), 1_000.0)
        .expect("boundary fixture must have finite enthalpy");
    let request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
        &resolved,
        initial,
        EquilibriumConstraint::ph_joules(
            101_325.0,
            101_325.0,
            TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
            1_000.0,
        )
        .unwrap(),
        TemperatureBounds::new(900.0, 1_100.0).unwrap(),
        thermochemistry.clone(),
    )
    .expect("boundary P,H request must validate")
    .with_ph_solve_mode(PhSolveMode::Monolithic)
    .with_solve_options(
        EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
                RustedSciTheSolver::LevenbergMarquardt,
            )))
            .expect("RST policy must validate"),
    );

    assert!(matches!(
        solve_resolved_ph(request),
        Err(ReactionExtentError::UnsupportedBackendCapability {
            backend,
            capability: "monolithic_p_h_symbolic_single_interval",
            ..
        }) if backend == "RustedSciThe"
    ));
    assert_eq!(before, live_library_file_snapshot());
}

/// Release-only bounded monolithic P,H regression using real gas and
/// pure-condensed records.
///
/// The initial liquid inventory is zero. The test therefore exercises the
/// trace-seeded active-set probe and verifies that the shared lifecycle can
/// publish the phase-aware result without mutating the local repository.
#[test]
#[ignore = "release-oriented real-data bounded P,H water story"]
fn live_bounded_water_pt_to_h_to_ph_preserves_phase_evidence_and_json() {
    let before = live_library_file_snapshot();
    let resolved = live_resolved_system();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("resolved water layout must remain valid");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0])
        .expect("live water inventory must match the resolved layout");
    let pressure = 101_325.0;
    let known_temperature = 350.0;
    let pt = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(known_temperature, pressure, pressure)
                .expect("reference P,T conditions must be valid"),
            initial.clone(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("reference bounded water P,T solve must succeed");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live gas/liquid water records must provide enthalpy capabilities");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(pt.component_moles(), known_temperature)
        .expect("reference water enthalpy must be finite");
    let constraint = EquilibriumConstraint::ph_joules(
        pressure,
        pressure,
        TotalEnthalpyJoules::new(target_enthalpy).expect("reference enthalpy must be finite"),
        345.0,
    )
    .expect("P,H constraint must validate");

    let ph = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            constraint,
            TemperatureBounds::new(325.0, 375.0).expect("water P,H bounds must be valid"),
            thermochemistry,
        )
        .expect("bounded live P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::Monolithic)
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        ),
    )
    .expect("bounded live P,H solve must recover the reference water state");

    assert_eq!(
        ph.report().solve_path(),
        PhSolvePath::MonolithicPhaseControl,
        "bounded monolithic P,H must use the shared phase-control lifecycle"
    );
    assert!(
        (ph.temperature() - known_temperature).abs() < 1.0e-4,
        "bounded P,H recovered {} K instead of {known_temperature} K",
        ph.temperature()
    );
    assert_eq!(ph.report().fixed_formulation_builds(), 0);
    assert_eq!(ph.report().fixed_formulation_reuses(), 0);
    assert!(ph.report().trials().is_empty());
    assert!(ph
        .report()
        .monolithic_evidence()
        .is_some_and(|evidence| evidence.phase_control_report().is_some()));
    let phase_control = ph
        .equilibrium()
        .phase_control_report()
        .expect("bounded monolithic P,H must retain lifecycle evidence");
    assert!(
        phase_control
            .transitions
            .iter()
            .any(|transition| { transition.activated.iter().any(|phase| phase.index() == 1) }),
        "zero-inventory liquid must be activated by the trace-seeded P,H probe"
    );
    assert!(
        phase_control
            .transitions
            .iter()
            .all(|transition| transition.transition_duration > std::time::Duration::ZERO),
        "published P,H phase transitions must retain a measured control-pass duration"
    );
    let validation = ph.equilibrium().accepted_solution().validation();
    assert!(validation.residual_l2_norm.is_finite());
    assert!(validation.max_abs_element_balance_error <= 1e-6);
    assert_eq!(before, live_library_file_snapshot());
}

/// Release-only route-equivalence story for a real water/liquid activation.
///
/// Both P,H routes receive the same immutable resolved system, inventory,
/// enthalpy target, bounds, and phase-control policy. This distinguishes a
/// genuine route regression from a difference caused by lookup or setup.
#[test]
#[ignore = "release real-data P,H monolithic/nested phase-control comparison"]
fn live_bounded_water_ph_monolithic_and_nested_routes_agree_after_phase_activation() {
    let before = live_library_file_snapshot();
    let resolved = live_resolved_system();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("resolved water layout must remain valid");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0])
        .expect("live water inventory must match the resolved layout");
    let pressure = 101_325.0;
    let reference_temperature = 350.0;
    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(reference_temperature, pressure, pressure)
                .expect("reference P,T conditions must be valid"),
            initial.clone(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("reference bounded water P,T solve must succeed");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live gas/liquid water records must provide enthalpy capabilities");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), reference_temperature)
        .expect("reference water enthalpy must be finite");
    let constraint = EquilibriumConstraint::ph_joules(
        pressure,
        pressure,
        TotalEnthalpyJoules::new(target_enthalpy).expect("reference enthalpy must be finite"),
        345.0,
    )
    .expect("P,H constraint must validate");
    let bounds = TemperatureBounds::new(325.0, 375.0).expect("water P,H bounds must be valid");
    let options = EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled);

    let monolithic_started = std::time::Instant::now();
    let monolithic = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial.clone(),
            constraint.clone(),
            bounds.clone(),
            thermochemistry.clone(),
        )
        .expect("monolithic P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::Monolithic)
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_solve_options(options.clone()),
    )
    .expect("monolithic P,H phase-control solve must succeed");

    let nested_started = std::time::Instant::now();
    let nested = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            constraint,
            bounds,
            thermochemistry,
        )
        .expect("nested P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_solve_options(options),
    )
    .expect("nested P,H phase-control solve must succeed");

    assert_eq!(
        monolithic.report().solve_path(),
        PhSolvePath::MonolithicPhaseControl
    );
    assert_eq!(nested.report().solve_path(), PhSolvePath::NestedTemperature);
    assert!(monolithic.report().trials().is_empty());
    assert!(!nested.report().trials().is_empty());
    assert!(monolithic.report().monolithic_evidence().is_some());
    assert!(nested.report().monolithic_evidence().is_none());
    assert!(
        (monolithic.temperature() - nested.temperature()).abs() <= 1.0e-4,
        "P,H routes recovered different temperatures: monolithic={} K, nested={} K",
        monolithic.temperature(),
        nested.temperature()
    );

    let liquid_index = 1usize;
    for solution in [&monolithic, &nested] {
        let phase_control = solution
            .equilibrium()
            .phase_control_report()
            .expect("phase-controlled P,H route must retain lifecycle evidence");
        assert!(
            phase_control.transitions.iter().any(|transition| {
                transition
                    .activated
                    .iter()
                    .any(|phase| phase.index() == liquid_index)
            }),
            "zero-inventory liquid must be activated by the P,H lifecycle"
        );
        let validation = solution.equilibrium().accepted_solution().validation();
        assert!(validation.residual_l2_norm.is_finite());
        assert!(validation.max_abs_element_balance_error <= 1.0e-6);
        assert!(solution.enthalpy_error().abs() <= solution.enthalpy_error_limit_joules());
    }

    for (component, (monolithic_moles, nested_moles)) in monolithic
        .equilibrium()
        .component_moles()
        .iter()
        .zip(nested.equilibrium().component_moles())
        .enumerate()
    {
        let tolerance = 1.0e-8 + 1.0e-6 * monolithic_moles.abs().max(nested_moles.abs());
        assert!(
            (monolithic_moles - nested_moles).abs() <= tolerance,
            "component {component} differs between P,H routes: monolithic={monolithic_moles:e}, nested={nested_moles:e}, tolerance={tolerance:e}"
        );
    }

    println!(
        "live bounded water P,H route comparison\n{}",
        Table::new([
            live_ph_path_row(
                "monolithic-phase-control",
                &monolithic,
                monolithic_started.elapsed()
            ),
            live_ph_path_row("nested-phase-control", &nested, nested_started.elapsed()),
        ])
        .with(Style::rounded())
    );
    assert_eq!(before, live_library_file_snapshot());
}

fn live_multiphase_pipeline_request() -> PhaseEquilibriumPipelineRequest {
    PhaseEquilibriumPipelineRequest::new(
        live_multiphase_spec(),
        vec![0.5, 0.25, 0.0],
        // The bundled NASA_cond liquid-water record is valid on 273.15-600 K.
        EquilibriumConditions::new(500.0, 101_325.0, 101_325.0)
            .expect("live temperature and pressure must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
}

#[test]
fn live_repository_resolves_real_nasa_phases_with_explicit_offline_policy() {
    let resolved = live_resolved_system();

    assert_eq!(resolved.phase_specs().len(), 2);
    assert_eq!(resolved.report().nist_fallback_enabled(), false);
    assert_eq!(resolved.report().phases().len(), 2);
    assert_eq!(
        resolved
            .phase_data()
            .get(&Some("gas".to_string()))
            .unwrap()
            .substances(),
        &["H2O".to_string(), "O2".to_string()]
    );
    assert_eq!(
        resolved
            .phase_data()
            .get(&Some("liquid".to_string()))
            .unwrap()
            .substances(),
        &["H2O".to_string()]
    );
}

#[test]
fn live_repository_distinguishes_gaseous_and_solid_water_records() {
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_gas_solid_water_spec(),
        live_repository(),
    )
    .expect("gas/solid water spec must resolve without NIST");

    assert!(!resolved.report().nist_fallback_enabled());
    assert_eq!(resolved.phase_specs().len(), 2);
    assert_eq!(
        resolved
            .phase_data()
            .get(&Some("gas".to_string()))
            .expect("gas phase must be present")
            .substances(),
        &["H2O".to_string(), "O2".to_string()]
    );
    assert_eq!(
        resolved
            .phase_data()
            .get(&Some("solid".to_string()))
            .expect("solid phase must be present")
            .substances(),
        &["H2O(s)".to_string()]
    );

    let report = resolved.report();
    let gas_search = report
        .phases()
        .iter()
        .find(|phase| phase.phase().as_option().as_deref() == Some("gas"))
        .expect("gas lookup report must be present");
    let solid_search = report
        .phases()
        .iter()
        .find(|phase| phase.phase().as_option().as_deref() == Some("solid"))
        .expect("solid lookup report must be present");
    assert!(gas_search.search().rows().iter().any(|row| {
        row.substance() == "H2O" && row.library() == "NASA_gas" && row.record_key() == "H2O"
    }));
    assert!(solid_search.search().rows().iter().any(|row| {
        row.substance() == "H2O(s)" && row.library() == "NASA_cond" && row.record_key() == "H2O(s)"
    }));
}

#[test]
fn live_reactive_fixture_records_are_pinned_for_deliberate_review() {
    let expected = [
        ("H2", 0xf8e3_cc98_a352_016f_u64),
        ("O2", 0xe714_832f_23e0_9998_u64),
        ("H2O", 0x0af2_fe41_48f4_c5ac_u64),
    ];

    for (substance, expected_fingerprint) in expected {
        assert_eq!(
            live_record_fingerprint("NASA_gas", substance),
            expected_fingerprint,
            "the live {substance} fixture changed; review its physical contract before updating this fingerprint"
        );
    }
}

#[test]
fn live_gas_pipeline_solves_without_mutating_the_resolved_repository_view() {
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_gas_spec(),
        Arc::clone(&live_repository()),
    )
    .expect("live gas spec must resolve against the bundled repository");
    let before_gas = resolved
        .phase_data()
        .get(&None)
        .expect("single-phase resolved data must use the canonical None phase key")
        .clone();

    let outcome = live_gas_pipeline_request()
        .solve()
        .expect("live pipeline must solve");

    assert_eq!(outcome.lookup_report(), outcome.resolved().report());
    assert_eq!(
        outcome.solution().build_report().lookup_report(),
        outcome.lookup_report()
    );
    assert_eq!(outcome.resolved().phase_specs().len(), 1);
    assert!(outcome
        .solution()
        .component_moles()
        .iter()
        .all(|value| value.is_finite()));
    assert!(outcome
        .solution()
        .component_moles()
        .iter()
        .sum::<f64>()
        .is_finite());

    assert_eq!(
        before_gas.substances(),
        resolved
            .phase_data()
            .get(&None)
            .expect("single-phase resolved data must remain present")
            .substances()
    );
    assert_eq!(
        before_gas.library_priorities(),
        resolved
            .phase_data()
            .get(&None)
            .expect("single-phase resolved data must remain present")
            .library_priorities()
    );
    assert_eq!(
        before_gas.explicit_search_map(),
        resolved
            .phase_data()
            .get(&None)
            .expect("single-phase resolved data must remain present")
            .explicit_search_map()
    );
}

#[test]
fn live_pipeline_keeps_canonical_thermochemistry_files_byte_for_byte_unchanged() {
    let before = live_library_file_snapshot();

    let outcome = live_multiphase_pipeline_request()
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .expect("read-only live pipeline must solve");

    assert!(outcome
        .solution()
        .component_moles()
        .iter()
        .all(|value| value.is_finite()));
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
fn live_reactive_gas_solves_with_conserved_elements_and_keq_evidence() {
    let outcome = live_reactive_gas_pipeline_request()
        .solve()
        .expect("live H2/O2/H2O equilibrium must solve");
    let validation = outcome.solution().accepted_solution().validation();

    assert!(outcome
        .solution()
        .component_moles()
        .iter()
        .all(|moles| moles.is_finite() && *moles >= 0.0));
    assert!(validation.residual_l2_norm.is_finite());
    assert!(validation.max_abs_element_balance_error <= 1e-8);
    let keq_status = outcome.solution().keq_validation_status();
    assert!(
        matches!(
            keq_status,
            Some(EquilibriumConstantCrossValidationStatus::Compared(report)) if report.accepted
        ),
        "independent live K_eq evidence was not accepted: {keq_status:?}"
    );
}

#[test]
fn live_oxygen_dissociation_migrated_from_classical_fixture() {
    // This deliberately keeps the physical invariant from the old O2/O test,
    // not its mutable solver state or one historical iterate. The typed
    // facade must preserve positive species, finite residual evidence, and the
    // oxygen-element inventory using real local thermochemistry.
    let outcome = PhaseEquilibriumPipelineRequest::new(
        live_oxygen_dissociation_spec(),
        vec![1.0, 1e-5],
        EquilibriumConditions::new(500.0, 101_325.0, 101_325.0)
            .expect("live O2/O conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .solve()
    .expect("migrated live O2/O fixture must solve through the typed facade");

    let moles = outcome.solution().component_moles();
    let validation = outcome.solution().accepted_solution().validation();
    assert_eq!(moles.len(), 2);
    assert!(moles.iter().all(|value| value.is_finite() && *value > 0.0));
    assert!(validation.residual_l2_norm.is_finite());
    assert!(validation.max_abs_element_balance_error <= 1e-8);
}

#[test]
fn live_reactive_gas_is_metamorphic_under_inventory_scaling() {
    // At fixed T and P, scaling every initial amount leaves ideal-gas
    // activities and the equilibrium composition unchanged. The physical
    // amounts must scale, while the acceptance gate must use an
    // absolute-plus-relative conservation tolerance instead of one absolute
    // threshold that is only meaningful near the reference inventory.
    let reference = live_reactive_gas_pipeline_request_with_scale(1.0)
        .solve()
        .expect("reference live reactive gas solve must succeed");
    let reference_moles = reference.solution().component_moles().to_vec();

    for scale in [1e-4, 1e4] {
        let outcome = live_reactive_gas_pipeline_request_with_scale(scale)
            .solve()
            .unwrap_or_else(|error| panic!("scaled live solve failed at {scale:e}: {error:?}"));
        let scaled_moles = outcome.solution().component_moles();
        assert_eq!(scaled_moles.len(), reference_moles.len());
        for (actual, reference) in scaled_moles.iter().zip(&reference_moles) {
            let expected = reference * scale;
            let relative_error = (actual - expected).abs() / expected.abs().max(1e-300);
            assert!(
                relative_error <= 5e-5,
                "scaled live amount drifted: actual={actual:e}, expected={expected:e}, relative_error={relative_error:e}"
            );
        }

        let validation = outcome.solution().accepted_solution().validation();
        let allowed = 1e-7 + 1e-7 * 2.1 * scale;
        assert!(
            validation.max_abs_element_balance_error <= allowed,
            "scale-aware conservation contract failed at {scale:e}: error={:e}, allowed={allowed:e}",
            validation.max_abs_element_balance_error
        );
    }
}

#[test]
fn live_bounded_phase_control_works_on_real_local_thermochemistry() {
    let outcome = live_multiphase_pipeline_request()
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .expect("bounded live pipeline must solve");

    let timing = outcome.timing_report();
    assert!(timing.enabled());
    assert!(timing.repository_lookup() > std::time::Duration::ZERO);
    assert!(timing.numerical_problem_preparation() > std::time::Duration::ZERO);
    assert!(timing.phase_control() > std::time::Duration::ZERO);
    assert!(timing.projection_build() > std::time::Duration::ZERO);
    assert!(timing.total() >= timing.phase_control());

    let components = outcome.solution().build_report().components();
    assert!(components.iter().all(|component| {
        let gibbs = component.standard_gibbs_at_conditions();
        gibbs.is_finite() && gibbs != 0.0
    }));
    let gas_water = components
        .iter()
        .find(|component| component.component().label() == "gas::H2O")
        .expect("gas water bridge row must exist");
    let liquid_water = components
        .iter()
        .find(|component| component.component().label() == "liquid::H2O")
        .expect("liquid water bridge row must exist");
    assert_eq!(gas_water.thermo_source().library(), "NASA_gas");
    assert_eq!(gas_water.thermo_source().record_key(), "H2O");
    assert_eq!(liquid_water.thermo_source().library(), "NASA_cond");
    assert_eq!(liquid_water.thermo_source().record_key(), "H2O(L)");
    assert_ne!(
        gas_water.standard_gibbs_at_conditions(),
        liquid_water.standard_gibbs_at_conditions()
    );

    assert_eq!(outcome.resolved().phase_specs().len(), 2);
    assert!(outcome
        .solution()
        .component_moles()
        .iter()
        .sum::<f64>()
        .is_finite());
    assert!(outcome
        .solution()
        .phase_total(&crate::Thermodynamics::phase_layout::PhaseId::new(Some(
            "liquid".to_string()
        )))
        .is_some());
    assert!(outcome
        .solution()
        .summary_rows()
        .iter()
        .any(|row| row.section == "phase_control"));
}

#[test]
#[ignore = "release live phase-control temperature-range characterization"]
fn live_bounded_phase_control_temperature_range_reuses_accepted_state() {
    let before = live_library_file_snapshot();
    let range = live_multiphase_pipeline_request()
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve_temperature_range(
            TemperatureGrid::new(vec![500.0, 525.0, 550.0])
                .expect("live phase-control temperature grid must validate"),
        )
        .expect("live bounded temperature range must solve");

    assert_eq!(range.points().len(), 3);
    assert!(range.report().phase_control_enabled());
    assert_eq!(range.report().formulation_builds(), 1);
    assert_eq!(range.report().formulation_reuses(), 2);
    assert!(range.report().phase_projection_cache_entries() >= 1);
    assert!(range.report().phase_prepared_cache_entries() >= 1);
    assert!(range.report().phase_rst_cache_entries() >= 1);
    assert!(range.points()[1].report().phase_set_reused());
    assert!(range.points()[2].report().phase_set_reused());
    assert!(range.points()[1].report().symbolic_parameter_reused());
    assert!(range.points()[2].report().symbolic_parameter_reused());
    assert_eq!(
        range.points()[1].report().formulation_build(),
        std::time::Duration::ZERO,
        "an unchanged active set must retarget RST parameters rather than rebuild its graph"
    );
    assert_eq!(
        range.points()[2].report().formulation_build(),
        std::time::Duration::ZERO,
        "an unchanged active set must retarget RST parameters rather than rebuild its graph"
    );
    let postprocessed = postprocess_temperature_range_solution(
        &range,
        &TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 5 },
            interpolation: TemperatureInterpolationPolicy {
                space: TemperatureInterpolationSpace::Linear,
                clamp: true,
            },
        },
    )
    .expect("typed bounded range should feed postprocessing");
    assert_eq!(postprocessed.raw.point_count(), 3);
    assert_eq!(
        postprocessed
            .resampled
            .as_ref()
            .expect("uniform grid should produce resampled data")
            .point_count(),
        5
    );
    assert_eq!(
        postprocessed.raw.labels().len(),
        range.points()[0].solution().component_moles().len()
    );
    assert!(range.points().iter().all(|point| {
        point
            .solution()
            .component_moles()
            .iter()
            .all(|moles| moles.is_finite() && *moles >= 0.0)
    }));
    for point in range.points() {
        let validation = point.solution().accepted_solution().validation();
        let scale = point
            .solution()
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * scale);
    }
    assert_eq!(before, live_library_file_snapshot());

    println!(
        "live bounded typed T-range: points={} transitions={} projections={} prepared={} rst_symbolic={} total={:?}",
        range.points().len(),
        range.report().phase_control_transitions(),
        range.report().phase_projection_cache_entries(),
        range.report().phase_prepared_cache_entries(),
        range.report().phase_rst_cache_entries(),
        range.report().total(),
    );
}

#[test]
#[ignore = "explicit release characterization run"]
fn live_phase_control_timing_report_prints_stage_breakdown() {
    // This is intentionally an ignored characterization test: it exercises
    // the real local repository and prints measurements for human review,
    // while the non-ignored regression above keeps the behavioral contract
    // quiet and deterministic.
    let outcome = live_multiphase_pipeline_request()
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .expect("live timing characterization must solve");
    let timing = outcome.timing_report();

    println!("live phase-control timing report");
    println!(
        "  repository_lookup              {:?}",
        timing.repository_lookup()
    );
    println!(
        "  thermochemistry_preparation   {:?}",
        timing.thermochemistry_preparation()
    );
    println!(
        "  numeric_closure_construction  {:?}",
        timing.numeric_closure_construction()
    );
    println!(
        "  symbolic_construction         {:?}",
        timing.symbolic_construction()
    );
    println!(
        "  equation_construction         {:?}",
        timing.equation_construction()
    );
    println!(
        "  numerical_problem_preparation {:?}",
        timing.numerical_problem_preparation()
    );
    println!(
        "  projection_build               {:?}",
        timing.projection_build()
    );
    println!(
        "  nonlinear_solve                {:?}",
        timing.nonlinear_solve()
    );
    println!(
        "  phase_control                  {:?}",
        timing.phase_control()
    );
    println!("  validation                     {:?}", timing.validation());
    println!(
        "  postprocessing                 {:?}",
        timing.postprocessing()
    );
    println!("  total                          {:?}", timing.total());
}

/// Returns NASA gas records whose formulas use no elements outside C/H/O.
///
/// This is the useful large element-defined inventory contract: a requested
/// element set is an allowed alphabet, so pure C, H, O and binary compounds
/// are valid candidates too. `search_by_exact_elements` is intentionally
/// stricter and is covered by its own regression above.
fn live_element_limited_gas_candidates_up_to(limit: usize) -> Vec<String> {
    assert!(limit > 0, "live candidate limit must be positive");
    let mut catalog = ThermoData::try_new().expect("bundled catalog must load");
    let mut candidates: Vec<String> = catalog
        .search_by_elements_only(vec!["C".into(), "H".into(), "O".into()])
        .into_iter()
        .filter(|name| {
            live_repository()
                .LibThermoData
                .get("NASA_gas")
                .is_some_and(|records| records.contains_key(name))
        })
        .collect();
    candidates.sort();
    candidates.dedup();
    assert!(
        candidates.len() >= limit,
        "the live NASA gas catalog must provide at least {limit} C/H/O-limited candidates, got {}",
        candidates.len(),
    );
    candidates.into_iter().take(limit).collect()
}

fn live_large_element_limited_gas_candidates() -> Vec<String> {
    live_element_limited_gas_candidates_up_to(20)
}

fn live_large_element_limited_gas_request() -> (Vec<String>, PhaseEquilibriumPipelineRequest) {
    let selected = live_large_element_limited_gas_candidates();
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(selected.clone()))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("large C/H/O-limited spec must remain valid");
    let request = PhaseEquilibriumPipelineRequest::new(
        spec,
        vec![1e-3; selected.len()],
        EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0)
            .expect("large live gas conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 });
    (selected, request)
}

fn live_large_element_limited_gas_resolved() -> (Vec<String>, ResolvedPhaseSystem) {
    live_element_limited_gas_resolved_up_to(20)
}

fn live_element_limited_gas_resolved_up_to(count: usize) -> (Vec<String>, ResolvedPhaseSystem) {
    let selected = live_element_limited_gas_candidates_up_to(count);
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(selected.clone()))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("large C/H/O-limited spec must remain valid");
    let resolved =
        SubstanceSystemFactory::resolve_phase_system_with_repository(spec, live_repository())
            .expect("large C/H/O-limited system must resolve against the local repository");
    (selected, resolved)
}

#[test]
#[ignore = "requires the bundled NASA gas catalog and validates the real-data P,H Jacobian"]
fn live_real_gas_ph_jacobian_matches_central_difference_near_interval_boundaries() {
    // Keep the fixture deliberately small: this test validates derivatives,
    // not solver throughput. The selected records are nevertheless resolved
    // through the same local repository and phase-data bridge as production.
    let before = live_library_file_snapshot();
    let (selected, resolved) = live_element_limited_gas_resolved_up_to(5);
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("real phase layout must be valid");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("real thermochemistry bundle must be complete");
    let bounds = TemperatureBounds::new(1_000.0, 1_100.0)
        .expect("real common NASA interval must be non-degenerate");
    let mut checked_entries = 0usize;
    let mut matrix_dimension = 0usize;

    // Exercise the complete block matrix at several extensive scales. The
    // unknowns are logarithmic, so this catches accidental absolute-scale
    // assumptions in the energy row without changing the physical state.
    for mole_scale in [1e-6, 1e-3, 1.0, 1e3] {
        let composition =
            MultiphaseInitialComposition::from_dense(&layout, vec![mole_scale; selected.len()])
                .expect("real composition must match the resolved layout");

        // The logistic coordinate cannot represent an exact endpoint. Values
        // immediately inside each boundary exercise the interval guards
        // without asking the formulation to evaluate out-of-domain data.
        for temperature in [1_000.0001, 1_050.0, 1_099.9999] {
            let molar_enthalpies = thermochemistry
                .evaluate_enthalpy(temperature)
                .expect("real enthalpy must evaluate in the common interval");
            let target = molar_enthalpies
                .iter()
                .zip(composition.moles())
                .map(|(enthalpy, moles)| enthalpy * moles)
                .sum::<f64>();
            let scale =
                EnthalpyScale::from_magnitudes(target, composition.moles(), &molar_enthalpies)
                    .expect("real enthalpy scale must be valid");
            let conditions = EquilibriumConditions::new(temperature, 101_325.0, 101_325.0)
                .expect("real fixed pressure conditions must be valid");
            let build_request = PhaseEquilibriumBuildRequest::new(
                &resolved,
                conditions,
                composition.clone(),
                TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
                SupportedPhaseModelPolicy::default(),
            )
            .expect("real bridge request must validate");
            let bundle = build_phase_equilibrium_problem_with_timing(
                build_request,
                EquilibriumTimingMode::Disabled,
            )
            .expect("real bridge bundle must build");
            let prepared = bundle
                .into_prepared_fixed_active_problem()
                .expect("real fixed active problem must prepare");
            let formulation = PreparedPhFormulation::new(
                prepared.prepared,
                thermochemistry.clone(),
                bounds,
                target,
                scale,
            )
            .expect("real P,H formulation must prepare");
            let unknowns = formulation
                .initial_unknowns(temperature)
                .expect("real P,H seed must be interiorized");
            let analytic = formulation
                .jacobian(&unknowns)
                .expect("real analytic Jacobian must evaluate");
            matrix_dimension = matrix_dimension.max(unknowns.len());
            let step = 1e-6;

            for column in 0..unknowns.len() {
                let mut plus = unknowns.clone();
                let mut minus = unknowns.clone();
                plus[column] += step;
                minus[column] -= step;
                let residual_plus = formulation
                    .residual(&plus)
                    .expect("real plus residual must evaluate");
                let residual_minus = formulation
                    .residual(&minus)
                    .expect("real minus residual must evaluate");
                for row in 0..unknowns.len() {
                    checked_entries += 1;
                    let finite_difference =
                        (residual_plus[row] - residual_minus[row]) / (2.0 * step);
                    let difference = (analytic[(row, column)] - finite_difference).abs();
                    let tolerance = 5e-4
                        * (1.0 + analytic[(row, column)].abs()).max(1.0 + finite_difference.abs());
                    assert!(
                        difference <= tolerance,
                        "real P,H Jacobian mismatch at scale={mole_scale}, T={temperature}, row={row}, column={column}: analytic={}, finite_difference={}, difference={difference}, tolerance={tolerance}",
                        analytic[(row, column)],
                        finite_difference,
                    );
                }
            }
        }
    }
    assert_eq!(before, live_library_file_snapshot());
    println!(
        "live real P,H Jacobian matrix: species={} scales={} temperatures={} dimension={} checked_entries={} status=OK",
        selected.len(),
        4,
        3,
        matrix_dimension,
        checked_entries,
    );
}

fn live_duration_ms(duration: std::time::Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1_000.0)
}

/// One compact row for the expensive real-data backend matrix.
///
/// The table is intentionally assembled at the test boundary. Backend
/// reports remain typed internally, while the operator gets one comparable
/// line per method instead of a stream of nearly identical prefixes.
#[derive(Debug, Clone, PartialEq, Eq, Tabled)]
struct LiveTemperatureRangeBackendRow {
    #[tabled(rename = "Backend")]
    backend: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
    #[tabled(rename = "Wall ms")]
    wall_ms: String,
    #[tabled(rename = "Mean ms")]
    mean_ms: String,
    #[tabled(rename = "Min ms")]
    min_ms: String,
    #[tabled(rename = "Median ms")]
    median_ms: String,
    #[tabled(rename = "Worst ms")]
    worst_ms: String,
    #[tabled(rename = "Worst T K")]
    worst_temperature: String,
    #[tabled(rename = "Worst nonlinear ms")]
    worst_nonlinear_ms: String,
    #[tabled(rename = "Worst symbolic ms")]
    worst_symbolic_ms: String,
    #[tabled(rename = "Worst accepted backend")]
    worst_accepted_backend: String,
    #[tabled(rename = "Worst attempts")]
    worst_attempts: String,
    #[tabled(rename = "Worst iterations")]
    worst_iterations: String,
    #[tabled(rename = "Builds")]
    builds: String,
    #[tabled(rename = "Initial setup ms")]
    initial_setup_ms: String,
    #[tabled(rename = "Initial symbolic ms")]
    initial_symbolic_ms: String,
    #[tabled(rename = "Initial problem prep ms")]
    initial_problem_preparation_ms: String,
    #[tabled(rename = "Formulation build ms")]
    formulation_build_ms: String,
    #[tabled(rename = "Reuses")]
    reuses: String,
    #[tabled(rename = "Symbolic updates")]
    symbolic_updates: String,
    #[tabled(rename = "Max residual")]
    max_residual: String,
    #[tabled(rename = "Max balance")]
    max_balance: String,
    #[tabled(rename = "Error")]
    error: String,
}

/// One accepted point selected for the slow-point release diagnosis.
///
/// The range report already contains all of this evidence. This flat test-only
/// view is intentionally limited to the three slowest points per backend so a
/// 100-point characterization remains readable while exposing the difference
/// between expensive symbolic preparation and time spent inside the opaque RST
/// nonlinear engine.
#[derive(Debug, Clone, PartialEq, Eq, Tabled)]
struct LiveTemperatureRangeSlowPointRow {
    #[tabled(rename = "Backend")]
    backend: String,
    #[tabled(rename = "Point")]
    point_index: String,
    #[tabled(rename = "T K")]
    temperature: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
    #[tabled(rename = "Nonlinear ms")]
    nonlinear_ms: String,
    #[tabled(rename = "Symbolic ms")]
    symbolic_ms: String,
    #[tabled(rename = "Closures ms")]
    numeric_closures_ms: String,
    #[tabled(rename = "Problem prep ms")]
    numerical_problem_ms: String,
    #[tabled(rename = "Validation ms")]
    validation_ms: String,
    #[tabled(rename = "Seed dlog inf")]
    seed_delta_log_inf: String,
    #[tabled(rename = "Seed log range")]
    seed_log_range: String,
    #[tabled(rename = "Result log range")]
    result_log_range: String,
    #[tabled(rename = "Accepted")]
    accepted_backend: String,
    #[tabled(rename = "Termination")]
    termination: String,
    #[tabled(rename = "Iterations")]
    iterations: String,
    #[tabled(rename = "Residual evals")]
    residual_evaluations: String,
    #[tabled(rename = "Jacobian evals")]
    jacobian_evaluations: String,
    #[tabled(rename = "Linear solves")]
    linear_solves: String,
    #[tabled(rename = "Backend ms")]
    backend_ms: String,
    #[tabled(rename = "Residual ms")]
    residual_ms: String,
    #[tabled(rename = "Jacobian ms")]
    jacobian_ms: String,
    #[tabled(rename = "Engine ms")]
    engine_ms: String,
}

/// One started backend attempt preserved through a failed temperature point.
///
/// A range failure used to be stringified at the continuation boundary. The
/// table therefore acts as a regression guard for the new typed wrapper as
/// well as an operator-facing diagnosis of a method that is not suitable for
/// this particular real-data inventory.
#[derive(Debug, Clone, PartialEq, Eq, Tabled)]
struct LiveTemperatureRangeFailureRow {
    #[tabled(rename = "Requested backend")]
    requested_backend: String,
    #[tabled(rename = "Point")]
    point_index: String,
    #[tabled(rename = "T K")]
    temperature: String,
    #[tabled(rename = "Error kind")]
    error_kind: String,
    #[tabled(rename = "Attempt")]
    attempt_index: String,
    #[tabled(rename = "Attempt backend")]
    attempt_backend: String,
    #[tabled(rename = "Outcome")]
    outcome: String,
    #[tabled(rename = "Termination")]
    termination: String,
    #[tabled(rename = "Iterations")]
    iterations: String,
    #[tabled(rename = "Residual evals")]
    residual_evaluations: String,
    #[tabled(rename = "Jacobian evals")]
    jacobian_evaluations: String,
    #[tabled(rename = "Linear solves")]
    linear_solves: String,
    #[tabled(rename = "Backend ms")]
    backend_ms: String,
    #[tabled(rename = "Residual ms")]
    residual_ms: String,
    #[tabled(rename = "Jacobian ms")]
    jacobian_ms: String,
    #[tabled(rename = "Engine ms")]
    engine_ms: String,
    #[tabled(rename = "Detail")]
    detail: String,
}

fn live_log_coordinate_range(values: &[f64]) -> String {
    let Some((&first, rest)) = values.split_first() else {
        return "-".to_string();
    };
    if !first.is_finite() || rest.iter().any(|value| !value.is_finite()) {
        return "non-finite".to_string();
    }
    let (minimum, maximum) = rest
        .iter()
        .fold((first, first), |(minimum, maximum), &value| {
            (minimum.min(value), maximum.max(value))
        });
    format!("[{minimum:.3}, {maximum:.3}]")
}

fn live_log_coordinate_delta_inf(seed: &[f64], result: &[f64]) -> String {
    if seed.len() != result.len() || seed.iter().chain(result).any(|value| !value.is_finite()) {
        return "non-finite".to_string();
    }
    let delta = seed
        .iter()
        .zip(result)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    format!("{delta:.3e}")
}

fn live_attempt_metric<T>(
    attempt: Option<&SolverAttemptReport>,
    render: impl FnOnce(
        &crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptMetrics,
    ) -> T,
) -> String
where
    T: ToString,
{
    attempt
        .and_then(|attempt| attempt.metrics.as_ref())
        .map(render)
        .map(|value| value.to_string())
        .unwrap_or_else(|| "-".to_string())
}

fn live_attempt_evaluation_ms(
    attempt: Option<&SolverAttemptReport>,
    render: impl FnOnce(
        &crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverEvaluationTiming,
    ) -> u128,
) -> String {
    attempt
        .and_then(|attempt| attempt.metrics.as_ref())
        .and_then(|metrics| metrics.evaluation_timing.as_ref())
        .map(render)
        .map(|micros| format!("{:.3}", micros as f64 / 1_000.0))
        .unwrap_or_else(|| "-".to_string())
}

fn live_temperature_range_failure_rows(
    requested_backend: &str,
    error: &ReactionExtentError,
) -> Vec<LiveTemperatureRangeFailureRow> {
    let (point_index, temperature, cause) = match error {
        ReactionExtentError::TemperatureRangePointFailed {
            point_index,
            temperature,
            cause,
        } => (*point_index, *temperature, cause.as_ref()),
        _ => {
            return vec![LiveTemperatureRangeFailureRow {
                requested_backend: requested_backend.to_string(),
                point_index: "-".to_string(),
                temperature: "-".to_string(),
                error_kind: format!("{:?}", error.kind()),
                attempt_index: "-".to_string(),
                attempt_backend: "-".to_string(),
                outcome: "-".to_string(),
                termination: "-".to_string(),
                iterations: "-".to_string(),
                residual_evaluations: "-".to_string(),
                jacobian_evaluations: "-".to_string(),
                linear_solves: "-".to_string(),
                backend_ms: "-".to_string(),
                residual_ms: "-".to_string(),
                jacobian_ms: "-".to_string(),
                engine_ms: "-".to_string(),
                detail: error.to_string(),
            }]
        }
    };
    let attempts = match cause {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => attempts.as_slice(),
        _ => &[],
    };
    if attempts.is_empty() {
        return vec![LiveTemperatureRangeFailureRow {
            requested_backend: requested_backend.to_string(),
            point_index: point_index.to_string(),
            temperature: format!("{temperature:.1}"),
            error_kind: format!("{:?}", cause.kind()),
            attempt_index: "-".to_string(),
            attempt_backend: "-".to_string(),
            outcome: "-".to_string(),
            termination: "-".to_string(),
            iterations: "-".to_string(),
            residual_evaluations: "-".to_string(),
            jacobian_evaluations: "-".to_string(),
            linear_solves: "-".to_string(),
            backend_ms: "-".to_string(),
            residual_ms: "-".to_string(),
            jacobian_ms: "-".to_string(),
            engine_ms: "-".to_string(),
            detail: cause.to_string(),
        }];
    }

    attempts
        .iter()
        .enumerate()
        .map(|(attempt_index, attempt)| LiveTemperatureRangeFailureRow {
            requested_backend: requested_backend.to_string(),
            point_index: point_index.to_string(),
            temperature: format!("{temperature:.1}"),
            error_kind: format!("{:?}", cause.kind()),
            attempt_index: attempt_index.to_string(),
            attempt_backend: format!("{:?}", attempt.backend),
            outcome: live_attempt_outcome_label(attempt),
            termination: live_attempt_metric(Some(attempt), |metrics| {
                format!("{:?}", metrics.termination)
            }),
            iterations: live_attempt_metric(Some(attempt), |metrics| metrics.iterations),
            residual_evaluations: live_attempt_metric(Some(attempt), |metrics| {
                metrics.residual_evaluations
            }),
            jacobian_evaluations: live_attempt_metric(Some(attempt), |metrics| {
                metrics.jacobian_evaluations
            }),
            linear_solves: live_attempt_metric(Some(attempt), |metrics| metrics.linear_solves),
            backend_ms: live_attempt_metric(Some(attempt), |metrics| metrics.elapsed_millis),
            residual_ms: live_attempt_evaluation_ms(Some(attempt), |timing| {
                timing.residual_evaluation_micros
            }),
            jacobian_ms: live_attempt_evaluation_ms(Some(attempt), |timing| {
                timing.jacobian_evaluation_micros
            }),
            engine_ms: live_attempt_evaluation_ms(Some(attempt), |timing| {
                timing.solver_overhead_micros
            }),
            detail: attempt.summary(),
        })
        .collect()
}

fn live_attempt_outcome_label(attempt: &SolverAttemptReport) -> String {
    match &attempt.outcome {
        SolverAttemptOutcome::Accepted => "Accepted".to_string(),
        SolverAttemptOutcome::Failed { kind, .. } => format!("Failed({kind:?})"),
        SolverAttemptOutcome::RejectedCandidate { .. } => "RejectedCandidate".to_string(),
        SolverAttemptOutcome::Skipped { .. } => "Skipped".to_string(),
    }
}

#[test]
fn live_temperature_range_failure_diagnostics_preserve_backend_metrics() {
    let error = ReactionExtentError::TemperatureRangePointFailed {
        point_index: 7,
        temperature: 1_010.0,
        cause: Box::new(ReactionExtentError::AllBackendsFailed {
            attempts: vec![SolverAttemptReport {
                backend: SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
                outcome: SolverAttemptOutcome::Failed {
                    kind: SolverAttemptFailureKind::ResidualEvaluation,
                    reason: "symbolic residual returned NaN".to_string(),
                },
                metrics: Some(SolverAttemptMetrics {
                    termination: SolverTermination::Stagnation,
                    backend_converged: false,
                    iterations: 3,
                    residual_evaluations: 4,
                    jacobian_evaluations: 3,
                    linear_solves: 3,
                    elapsed_millis: 17,
                    evaluation_timing: None,
                }),
            }],
        }),
    };

    let rows = live_temperature_range_failure_rows("rst_nielsen_lm", &error);

    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].point_index, "7");
    assert_eq!(rows[0].temperature, "1010.0");
    assert_eq!(rows[0].error_kind, "AllBackendsFailed");
    assert_eq!(
        rows[0].attempt_backend,
        "RustedSciThe(NielsenLevenbergMarquardt)"
    );
    assert_eq!(rows[0].outcome, "Failed(ResidualEvaluation)");
    assert_eq!(rows[0].termination, "Stagnation");
    assert_eq!(rows[0].iterations, "3");
    assert_eq!(rows[0].backend_ms, "17");
}

/// One compact row for the release-only real-data `P,H` backend matrix.
///
/// `P,H` is an outer scalar solve whose cost is dominated by repeated inner
/// `P,T` equilibrium solves. Keeping both layers in one row makes it clear
/// whether a backend is slow because of its nonlinear work or because the
/// bracket needed more thermodynamic evaluations.
#[derive(Debug, Clone, PartialEq, Eq, Tabled)]
struct LivePhBackendRow {
    #[tabled(rename = "Backend")]
    backend: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Outer trials")]
    outer_trials: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
    #[tabled(rename = "Wall ms")]
    wall_ms: String,
    #[tabled(rename = "Inner nonlinear ms")]
    inner_nonlinear_ms: String,
    #[tabled(rename = "Fixed builds")]
    fixed_builds: String,
    #[tabled(rename = "Fixed reuses")]
    fixed_reuses: String,
    #[tabled(rename = "Inner attempts")]
    inner_attempts: String,
    #[tabled(rename = "Inner iterations")]
    inner_iterations: String,
    #[tabled(rename = "Solved T K")]
    temperature: String,
    #[tabled(rename = "Scaled H error")]
    scaled_enthalpy_error: String,
    #[tabled(rename = "Residual")]
    residual: String,
    #[tabled(rename = "Balance")]
    balance: String,
    #[tabled(rename = "Error")]
    error: String,
}

/// One path-level row for the same real-data P,H release characterization.
///
/// Unlike the backend matrix above, this compares the two P,H formulations
/// themselves. `trials` is intentionally zero for monolithic mode; its
/// coupled backend trace is reported through `monolithic_evidence` instead.
#[derive(Debug, Clone, PartialEq, Eq, Tabled)]
struct LivePhPathRow {
    #[tabled(rename = "Path")]
    path: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Report ms")]
    report_ms: String,
    #[tabled(rename = "Thermo prep ms")]
    thermochemistry_preparation_ms: String,
    #[tabled(rename = "Wall ms")]
    wall_ms: String,
    #[tabled(rename = "Inner solves")]
    inner_solves: String,
    #[tabled(rename = "Trials")]
    trials: String,
    #[tabled(rename = "Attempts")]
    attempts: String,
    #[tabled(rename = "Iterations")]
    iterations: String,
    #[tabled(rename = "Builds")]
    builds: String,
    #[tabled(rename = "Reuses")]
    reuses: String,
    #[tabled(rename = "Transitions")]
    transitions: String,
    #[tabled(rename = "Residual evals")]
    residual_evaluations: String,
    #[tabled(rename = "Jacobian evals")]
    jacobian_evaluations: String,
    #[tabled(rename = "Scaled H error")]
    scaled_enthalpy_error: String,
    #[tabled(rename = "Residual")]
    residual: String,
    #[tabled(rename = "Balance")]
    balance: String,
    #[tabled(rename = "Error")]
    error: String,
}

fn live_ph_metric_totals(report: &PhTemperatureSolveReport) -> (usize, usize, usize) {
    if let Some(evidence) = report.monolithic_evidence() {
        return (
            1,
            evidence.residual_evaluations(),
            evidence.jacobian_evaluations(),
        );
    }

    let mut inner_solves = 0usize;
    let mut residual_evaluations = 0usize;
    let mut jacobian_evaluations = 0usize;
    for trial in report.trials() {
        if let Some(evidence) = trial.inner_evidence() {
            inner_solves += 1;
            residual_evaluations += evidence.residual_evaluations();
            jacobian_evaluations += evidence.jacobian_evaluations();
        }
    }
    (inner_solves, residual_evaluations, jacobian_evaluations)
}

fn live_ph_path_row(
    label: &str,
    solution: &FixedPressureEnthalpySolution,
    wall_time: std::time::Duration,
) -> LivePhPathRow {
    let report = solution.report();
    let (inner_solves, residual_evaluations, jacobian_evaluations) = live_ph_metric_totals(report);
    let validation = solution.equilibrium().accepted_solution().validation();
    LivePhPathRow {
        path: label.to_string(),
        status: "OK".to_string(),
        report_ms: live_duration_ms(report.timing().total()),
        thermochemistry_preparation_ms: live_duration_ms(
            report.inner_timing().thermochemistry_preparation(),
        ),
        wall_ms: live_duration_ms(wall_time),
        inner_solves: inner_solves.to_string(),
        trials: report.trials().len().to_string(),
        attempts: report.inner_backend_attempts().to_string(),
        iterations: report.inner_nonlinear_iterations().to_string(),
        builds: report.fixed_formulation_builds().to_string(),
        reuses: report.fixed_formulation_reuses().to_string(),
        transitions: report.phase_control_transitions().to_string(),
        residual_evaluations: residual_evaluations.to_string(),
        jacobian_evaluations: jacobian_evaluations.to_string(),
        scaled_enthalpy_error: format!(
            "{:.3e}",
            solution.enthalpy_error() / report.enthalpy_scale_joules()
        ),
        residual: format!("{:.3e}", validation.residual_l2_norm),
        balance: format!("{:.3e}", validation.max_abs_element_balance_error),
        error: "-".to_string(),
    }
}

#[test]
#[ignore = "heavy release real-data P,H backend matrix"]
fn live_reactive_gas_ph_backend_matrix() {
    // Construct H* once from an independently accepted P,T state. Every row
    // then receives the exact same resolved records, inventory, energy target,
    // bounds, and strict Single backend policy. A backend failure is useful
    // release evidence and must not be masked by a fallback cascade.
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas must resolve for the P,H backend matrix");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("reactive gas layout must validate");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])
        .expect("reactive gas initial composition must validate");
    let pressure = 101_325.0;
    let known_temperature = 2_500.0;
    let reference_options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
        .expect("reference single backend policy must validate");
    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(known_temperature, pressure, pressure)
                .expect("reference P,T conditions must validate"),
            initial.clone(),
        )
        .with_solve_options(reference_options),
    )
    .expect("reference P,T state must solve with legacy NR");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live NASA records must provide a thermochemistry bundle");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), known_temperature)
        .expect("reference P,T state must have finite total enthalpy");
    let bounds = TemperatureBounds::new(1_900.0, 2_900.0)
        .expect("P,H matrix temperature bounds must validate");
    let backends = [
        (
            "rst_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ),
        (
            "rst_minpack_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
        ),
        (
            "rst_nielsen_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        ),
        (
            "rst_trust_region_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
        ),
        (
            "rst_powell_dogleg",
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg),
        ),
        (
            "rst_damped_newton",
            SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
        ),
        ("legacy_lm", SolverBackend::Legacy(Solvers::LM)),
        ("legacy_nr", SolverBackend::Legacy(Solvers::NR)),
        ("legacy_tr", SolverBackend::Legacy(Solvers::TR)),
    ];

    let mut successful_backends = 0usize;
    let mut failed_backends = Vec::new();
    let mut rows = Vec::with_capacity(backends.len());
    for (name, backend) in backends {
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(backend))
            .expect("single backend policy must validate")
            .with_timing_mode(EquilibriumTimingMode::Enabled);
        let request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial.clone(),
            EquilibriumConstraint::ph_joules(
                pressure,
                pressure,
                TotalEnthalpyJoules::new(target_enthalpy)
                    .expect("reference total enthalpy must remain finite"),
                2_200.0,
            )
            .expect("P,H matrix constraint must validate"),
            bounds,
            thermochemistry.clone(),
        )
        .expect("P,H matrix request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_solve_options(options);
        let started = std::time::Instant::now();

        match solve_resolved_ph(request) {
            Ok(solution) => {
                successful_backends += 1;
                let report = solution.report();
                let validation = solution.equilibrium().accepted_solution().validation();
                assert!(report.trials().len() >= 2);
                assert_eq!(report.fixed_formulation_builds(), 1);
                assert_eq!(
                    report.fixed_formulation_reuses(),
                    report.trials().len().saturating_sub(1)
                );
                assert!(
                    (solution.temperature() - known_temperature).abs() < 1.0e-4,
                    "{name} recovered {} K instead of {known_temperature} K",
                    solution.temperature()
                );
                assert!(validation.residual_l2_norm.is_finite());
                assert!(validation.max_abs_element_balance_error <= 1.0e-6);

                rows.push(LivePhBackendRow {
                    backend: name.to_string(),
                    status: "OK".to_string(),
                    outer_trials: report.trials().len().to_string(),
                    total_ms: live_duration_ms(report.timing().total()),
                    wall_ms: live_duration_ms(started.elapsed()),
                    inner_nonlinear_ms: live_duration_ms(report.inner_timing().nonlinear_solve()),
                    fixed_builds: report.fixed_formulation_builds().to_string(),
                    fixed_reuses: report.fixed_formulation_reuses().to_string(),
                    inner_attempts: report.inner_backend_attempts().to_string(),
                    inner_iterations: report.inner_nonlinear_iterations().to_string(),
                    temperature: format!("{:.6}", solution.temperature()),
                    scaled_enthalpy_error: format!(
                        "{:.3e}",
                        solution.enthalpy_error() / report.enthalpy_scale_joules()
                    ),
                    residual: format!("{:.3e}", validation.residual_l2_norm),
                    balance: format!("{:.3e}", validation.max_abs_element_balance_error),
                    error: "-".to_string(),
                });
            }
            Err(error) => {
                failed_backends.push(name);
                rows.push(LivePhBackendRow {
                    backend: name.to_string(),
                    status: "FAILED".to_string(),
                    outer_trials: "-".to_string(),
                    total_ms: "-".to_string(),
                    wall_ms: live_duration_ms(started.elapsed()),
                    inner_nonlinear_ms: "-".to_string(),
                    fixed_builds: "-".to_string(),
                    fixed_reuses: "-".to_string(),
                    inner_attempts: "-".to_string(),
                    inner_iterations: "-".to_string(),
                    temperature: "-".to_string(),
                    scaled_enthalpy_error: "-".to_string(),
                    residual: "-".to_string(),
                    balance: "-".to_string(),
                    error: format!("{error:?}"),
                });
            }
        }
    }

    assert!(
        successful_backends > 0,
        "the live P,H matrix must have at least one successful backend"
    );
    println!("live real P,H backend matrix (H2/O2/H2O, T*=2500 K)");
    println!("{}", Table::new(rows).with(Style::rounded()));
    println!(
        "summary: attempted={} successful={} failed={failed_backends:?}",
        successful_backends + failed_backends.len(),
        successful_backends,
    );
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
#[ignore = "explicit release characterization run"]
fn live_large_element_limited_gas_timing_report() {
    // Exercise the real candidate-selection path rather than embedding a
    // synthetic species list. The C/H/O-limited alphabet keeps the chemical
    // universe bounded while still providing a materially larger NASA gas
    // system.
    let (selected, request) = live_large_element_limited_gas_request();
    let outcome = request
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        )
        .solve()
        .expect("large C/H/O-limited live gas solve must succeed");
    let moles = outcome.solution().component_moles();
    assert_eq!(moles.len(), selected.len());
    assert!(moles.iter().all(|value| value.is_finite() && *value > 0.0));
    let validation = outcome.solution().accepted_solution().validation();
    let inventory_scale = moles.iter().map(|value| value.abs()).sum::<f64>();
    let balance_limit = 1e-6 + 1e-6 * inventory_scale;
    assert!(validation.residual_l2_norm.is_finite());
    assert!(
        validation.max_abs_element_balance_error <= balance_limit,
        "large real-data timing story violated conservation: balance={} limit={}",
        validation.max_abs_element_balance_error,
        balance_limit
    );
    let timing = outcome.timing_report();

    println!("live large C/H/O-limited gas timing report");
    println!("  selected_species               {}", selected.len());
    println!("  species                        {:?}", selected);
    println!(
        "  repository_lookup              {:?}",
        timing.repository_lookup()
    );
    println!(
        "  thermochemistry_preparation   {:?}",
        timing.thermochemistry_preparation()
    );
    println!(
        "  numeric_closure_construction  {:?}",
        timing.numeric_closure_construction()
    );
    println!(
        "  symbolic_construction         {:?}",
        timing.symbolic_construction()
    );
    println!(
        "  equation_construction         {:?}",
        timing.equation_construction()
    );
    println!(
        "  numerical_problem_preparation {:?}",
        timing.numerical_problem_preparation()
    );
    println!(
        "  nonlinear_solve                {:?}",
        timing.nonlinear_solve()
    );
    println!(
        "  postprocessing                 {:?}",
        timing.postprocessing()
    );
    println!(
        "  validation                     residual={:.3e} balance={:.3e} limit={:.3e}",
        validation.residual_l2_norm, validation.max_abs_element_balance_error, balance_limit,
    );
    println!("  total                          {:?}", timing.total());
}

#[test]
#[ignore = "explicit release solver characterization run"]
fn live_large_element_limited_solver_matrix_compares_rst_and_legacy_backends() {
    // Every backend receives the same real resolved problem and the same
    // initial composition. This is a characterization matrix, not a fallback
    // cascade: failures are reported per backend instead of being hidden by a
    // later successful attempt.
    let (selected, request) = live_large_element_limited_gas_request();
    let backends = [
        (
            "rst_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ),
        (
            "rst_minpack_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
        ),
        (
            "rst_nielsen_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        ),
        (
            "rst_trust_region_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
        ),
        (
            "rst_powell_dogleg",
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg),
        ),
        (
            "rst_damped_newton",
            SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
        ),
        ("legacy_lm", SolverBackend::Legacy(Solvers::LM)),
        ("legacy_nr", SolverBackend::Legacy(Solvers::NR)),
        ("legacy_tr", SolverBackend::Legacy(Solvers::TR)),
    ];

    println!(
        "live large C/H/O-limited solver matrix: species={}",
        selected.len()
    );
    let mut successful = 0;
    for (name, backend) in backends {
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(backend))
            .expect("single backend policy must validate")
            .with_timing_mode(EquilibriumTimingMode::Enabled);
        match request.clone().with_solve_options(options).solve() {
            Ok(outcome) => {
                let solution = outcome.solution();
                let validation = solution.accepted_solution().validation();
                assert!(
                    solution
                        .component_moles()
                        .iter()
                        .all(|value| value.is_finite() && *value > 0.0),
                    "{name} returned invalid component moles"
                );
                assert!(validation.residual_l2_norm.is_finite());
                let inventory_scale = solution
                    .component_moles()
                    .iter()
                    .map(|value| value.abs())
                    .sum::<f64>();
                let balance_limit = 1e-6 + 1e-6 * inventory_scale;
                assert!(
                    validation.max_abs_element_balance_error <= balance_limit,
                    "{name} exceeded the scale-aware live balance contract: error={:e}, limit={balance_limit:e}",
                    validation.max_abs_element_balance_error
                );
                let timing = outcome.timing_report();
                println!(
                    "  {name:24} ok total={:?} nonlinear={:?} residual={:e} balance={:e}",
                    timing.total(),
                    timing.nonlinear_solve(),
                    validation.residual_l2_norm,
                    validation.max_abs_element_balance_error
                );
                successful += 1;
            }
            Err(error) => println!("  {name:24} failed: {error:?}"),
        }
    }
    assert!(
        successful > 0,
        "the real solver matrix must have at least one successful backend"
    );
}

#[test]
#[ignore = "release regression for bounded real-data backend candidates"]
fn live_bounded_backend_regression_avoids_nonfinite_and_iteration_failures() {
    // Keep this fixture materially smaller than the 100 x 50 characterization
    // sweep. It isolates the formerly unstable methods at the same local
    // 1000 K C/H/O starting condition and deliberately forbids fallback: a
    // failure here is direct backend evidence, not an orchestration outcome.
    let before = live_library_file_snapshot();
    let (selected, request) = live_large_element_limited_gas_request();
    let backends = [
        (
            "rst_nielsen_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        ),
        (
            "rst_powell_dogleg",
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg),
        ),
        ("legacy_tr", SolverBackend::Legacy(Solvers::TR)),
    ];

    for (name, backend) in backends {
        let outcome = request
            .clone()
            .with_solve_options(
                EquilibriumSolveOptions::new()
                    .with_solver_policy(SolverPolicy::Single(backend))
                    .expect("single backend policy must validate")
                    .with_timing_mode(EquilibriumTimingMode::Enabled),
            )
            .solve()
            .unwrap_or_else(|error| {
                panic!(
                    "{name} must solve the bounded local {}-species regression fixture: {error:?}",
                    selected.len()
                )
            });
        let validation = outcome.solution().accepted_solution().validation();
        assert!(
            validation.residual_l2_norm.is_finite(),
            "{name} returned a non-finite residual"
        );
        assert!(
            outcome
                .solution()
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles > 0.0),
            "{name} returned an invalid mole vector"
        );
    }
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
#[ignore = "explicit release real-data scaling characterization run"]
fn live_element_limited_release_scaling_matrix() {
    // This intentionally uses the same allowed-element search contract as the
    // 20-species backend matrix, but grows only the real local candidate set.
    // It characterizes solver scaling without changing the catalog or hiding a
    // failed point behind a later fallback backend.
    let before = live_library_file_snapshot();
    for count in [20usize, 50, 100] {
        let selected = live_element_limited_gas_candidates_up_to(count);
        let spec =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(selected.clone()))
                .with_library_priorities(vec!["NASA_gas".to_string()])
                .with_search_in_nist(false)
                .build()
                .expect("real exact-element scaling spec must remain valid");
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
            .expect("single legacy NR policy must validate")
            .with_timing_mode(EquilibriumTimingMode::Enabled);
        let started = std::time::Instant::now();
        let outcome = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![1e-3; selected.len()],
            EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0)
                .expect("real scaling conditions must be valid"),
        )
        .with_repository(live_repository())
        .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
        .with_solve_options(options)
        .solve()
        .unwrap_or_else(|error| panic!("real {count}-species solve failed: {error:?}"));

        let solution = outcome.solution();
        let validation = solution.accepted_solution().validation();
        let inventory_scale = solution
            .component_moles()
            .iter()
            .map(|value| value.abs())
            .sum::<f64>();
        let balance_limit = 1e-6 + 1e-6 * inventory_scale;
        assert_eq!(solution.component_moles().len(), count);
        assert!(solution
            .component_moles()
            .iter()
            .all(|value| value.is_finite() && *value > 0.0));
        assert!(validation.residual_l2_norm.is_finite());
        assert!(
            validation.max_abs_element_balance_error <= balance_limit,
            "{count}-species live balance exceeded limit: error={:e}, limit={balance_limit:e}",
            validation.max_abs_element_balance_error
        );

        println!(
            "live C/H/O-limited scaling: backend=legacy-nr species={} total={:?} nonlinear={:?} residual={:e} balance={:e}",
            count,
            started.elapsed(),
            outcome.timing_report().nonlinear_solve(),
            validation.residual_l2_norm,
            validation.max_abs_element_balance_error,
        );
    }
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
#[ignore = "explicit release typed temperature-range backend characterization run"]
fn live_exact_element_typed_temperature_range_backend_matrix() {
    // This matrix deliberately uses a Single backend policy for every run.
    // A failed method is evidence about that method, not an invitation to
    // silently promote the point to another backend.
    let (selected, resolved) = live_large_element_limited_gas_resolved();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("real resolved phases must produce a canonical range layout");
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; selected.len()])
        .expect("real range composition must match the layout");
    let grid = TemperatureGrid::new(vec![1_000.0, 1_050.0, 1_100.0])
        .expect("real range grid must be valid");
    let backends = [
        (
            "rst_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ),
        (
            "rst_minpack_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
        ),
        (
            "rst_nielsen_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        ),
        (
            "rst_trust_region_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
        ),
        (
            "rst_powell_dogleg",
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg),
        ),
        (
            "rst_damped_newton",
            SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
        ),
        ("legacy_lm", SolverBackend::Legacy(Solvers::LM)),
        ("legacy_nr", SolverBackend::Legacy(Solvers::NR)),
        ("legacy_tr", SolverBackend::Legacy(Solvers::TR)),
    ];

    let mut successful_backends = 0usize;
    for (name, backend) in backends {
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(backend))
            .expect("single backend policy must validate")
            .with_timing_mode(EquilibriumTimingMode::Enabled);
        let result = TemperatureRangeRequest::new(
            &resolved,
            composition.clone(),
            101_325.0,
            101_325.0,
            grid.clone(),
        )
        .expect("real typed range request must validate")
        .with_solve_options(options)
        .solve();

        match result {
            Ok(range) => {
                successful_backends += 1;
                if matches!(backend, SolverBackend::RustedSciThe(_)) {
                    assert!(range.report().symbolic_problem_reused());
                    assert_eq!(
                        range.report().symbolic_parameter_updates(),
                        range.points().len().saturating_sub(1),
                        "parameterized RST graph must be reused after the first accepted point"
                    );
                    assert!(range.points().iter().all(|point| {
                        point.report().formulation_build() == std::time::Duration::ZERO
                    }), "changing a thermochemical coefficient interval must not rebuild the RST graph");
                }
                for point in range.points() {
                    let validation = point.solution().accepted_solution().validation();
                    let inventory_scale = point
                        .solution()
                        .component_moles()
                        .iter()
                        .map(|value| value.abs())
                        .sum::<f64>();
                    let balance_limit = 1e-6 + 1e-6 * inventory_scale;
                    assert_eq!(point.solution().component_moles().len(), selected.len());
                    assert!(point
                        .solution()
                        .component_moles()
                        .iter()
                        .all(|value| value.is_finite() && *value > 0.0));
                    assert!(validation.residual_l2_norm.is_finite());
                    assert!(validation.max_abs_element_balance_error <= balance_limit);
                }
                println!(
                    "live typed T-range backend: backend={} species={} points={} total={:?} mean={:?} worst={:?}",
                    name,
                    selected.len(),
                    range.points().len(),
                    range.report().point_timing().total(),
                    range.report().point_timing().mean(),
                    range.report().point_timing().worst(),
                );
            }
            Err(error) => println!(
                "live typed T-range backend: backend={} species={} failed={error:?}",
                name,
                selected.len(),
            ),
        }
    }
    assert!(
        successful_backends > 0,
        "the real typed temperature-range matrix must have at least one successful backend"
    );
}

#[test]
#[ignore = "release diagnostic: compare reusable and baked RST Damped Newton graphs"]
fn live_damped_newton_parameterized_graph_matches_baked_graph_outcome() {
    // The 100x50 release characterization rejected Damped Newton at 1000 K.
    // Build both graph representations from the identical real 20-species
    // NASA problem and seed.  A difference in acceptance is a regression in
    // the reusable graph; identical rejection is evidence that this is a
    // method-specific limitation rather than a changed thermodynamic model.
    let before = live_library_file_snapshot();
    let (selected, resolved) = live_large_element_limited_gas_resolved();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("real Damped Newton layout must validate");
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; selected.len()])
        .expect("real Damped Newton composition must validate");
    let conditions = EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0)
        .expect("real Damped Newton conditions must validate");

    let build_prepared = || {
        let bundle = build_phase_equilibrium_problem_with_timing(
            PhaseEquilibriumBuildRequest::new(
                &resolved,
                conditions,
                composition.clone(),
                TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 },
                SupportedPhaseModelPolicy::default(),
            )
            .expect("real Damped Newton bridge request must validate"),
            EquilibriumTimingMode::Disabled,
        )
        .expect("real Damped Newton bridge must build");
        let symbolic_standard_gibbs = bundle.symbolic_standard_gibbs().to_vec();
        let prepared = PreparedEquilibriumProblem::new(bundle.into_problem())
            .expect("real Damped Newton problem must prepare");
        let seed = prepared.problem().initial_log_moles().clone();
        (prepared, symbolic_standard_gibbs, seed)
    };
    let policy = SolverPolicy::Single(SolverBackend::RustedSciThe(
        RustedSciTheSolver::DampedNewton,
    ));

    let (parameterized_prepared, parameterized_symbols, parameterized_seed) = build_prepared();
    let parameterized_graph =
        prepare_rst_symbolic_problem_from_prepared(&parameterized_prepared, &parameterized_symbols)
            .expect("parameterized real RST graph must build");
    let mut parameterized_runner =
        PreparedEquilibriumRunner::new(parameterized_prepared, parameterized_symbols)
            .expect("parameterized real RST runner must build");
    parameterized_runner.configure().solver_policy = Some(policy.clone());
    let parameterized =
        parameterized_runner.solve_from_seed_with_rst(parameterized_seed, &parameterized_graph);

    let (baked_prepared, baked_symbols, baked_seed) = build_prepared();
    let baked_graph = prepare_baked_rst_symbolic_problem_for_test(&baked_prepared, &baked_symbols)
        .expect("baked real RST graph must build");
    let mut baked_runner = PreparedEquilibriumRunner::new(baked_prepared, baked_symbols)
        .expect("baked real RST runner must build");
    baked_runner.configure().solver_policy = Some(policy);
    let baked = baked_runner.solve_from_seed_with_rst(baked_seed, &baked_graph);

    match (&parameterized, &baked) {
        (Ok(parameterized), Ok(baked)) => {
            assert!(
                (parameterized.solution.validation().residual_l2_norm
                    - baked.solution.validation().residual_l2_norm)
                    .abs()
                    <= 1e-8,
                "Damped Newton accepted materially different reusable and baked candidates"
            );
        }
        (Err(_), Err(_)) => {
            // Both forms reach the same strict acceptance boundary. The
            // printed reports remain useful release evidence, but there is no
            // parameterization regression to repair.
        }
        (parameterized, baked) => panic!(
            "Damped Newton changed outcome after G0 parameterization: parameterized={parameterized:?}, baked={baked:?}"
        ),
    }
    println!(
        "live Damped Newton graph comparison: species={} parameterized={} baked={}",
        selected.len(),
        if parameterized.is_ok() {
            "accepted"
        } else {
            "rejected"
        },
        if baked.is_ok() {
            "accepted"
        } else {
            "rejected"
        },
    );
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
#[ignore = "heavy release real-data temperature-range backend matrix"]
fn live_100_species_50_point_temperature_range_backend_matrix() {
    // This is an operator-run characterization test. It uses one immutable
    // resolved 100-species system and executes every concrete backend as a
    // strict Single policy over the same 50-point continuation grid. A failed
    // backend is reported explicitly and cannot be hidden by a fallback
    // backend or by publishing a partial range.
    let before = live_library_file_snapshot();
    let (selected, resolved) = live_element_limited_gas_resolved_up_to(100);
    assert_eq!(selected.len(), 100);
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("real 100-species range layout must validate");
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; 100])
        .expect("real 100-species initial composition must validate");
    let temperatures = (0..50)
        .map(|index| 1_000.0 + 10.0 * index as f64)
        .collect::<Vec<_>>();
    let grid = TemperatureGrid::new(temperatures)
        .expect("real 50-point temperature grid must be strictly ascending");
    let backends = [
        (
            "rst_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ),
        (
            "rst_minpack_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
        ),
        (
            "rst_nielsen_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        ),
        (
            "rst_trust_region_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
        ),
        (
            "rst_powell_dogleg",
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg),
        ),
        (
            "rst_damped_newton",
            SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
        ),
        ("legacy_lm", SolverBackend::Legacy(Solvers::LM)),
        ("legacy_nr", SolverBackend::Legacy(Solvers::NR)),
        ("legacy_tr", SolverBackend::Legacy(Solvers::TR)),
    ];

    let mut successful_backends = 0usize;
    let mut failed_backends: Vec<&str> = Vec::new();
    let mut rows = Vec::with_capacity(backends.len());
    let mut slow_point_rows = Vec::new();
    let mut failure_rows = Vec::new();
    for (name, backend) in backends {
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(backend))
            .expect("single backend policy must validate")
            .with_timing_mode(EquilibriumTimingMode::Enabled);
        let started = std::time::Instant::now();
        let result = TemperatureRangeRequest::new(
            &resolved,
            composition.clone(),
            101_325.0,
            101_325.0,
            grid.clone(),
        )
        .expect("real 100-species range request must validate")
        .with_solve_options(options)
        .solve();

        match result {
            Ok(range) => {
                successful_backends += 1;
                assert_eq!(range.points().len(), 50);
                assert_eq!(range.report().point_count(), 50);
                // A NASA record may cross a coefficient interval while the
                // grid advances.  Such a crossing legitimately rebuilds the
                // temperature-dependent formulation, so the contract is
                // accounting completeness rather than exactly one build.
                let formulation_builds = range.report().formulation_builds();
                let formulation_reuses = range.report().formulation_reuses();
                assert!(formulation_builds >= 1);
                assert!(formulation_reuses > 0);
                assert_eq!(formulation_builds + formulation_reuses, 50);

                let mut max_residual = 0.0_f64;
                let mut max_balance = 0.0_f64;
                let mut previous_accepted_log_moles = vec![1e-3_f64.ln(); 100];
                for (index, point) in range.points().iter().enumerate() {
                    let solution = point.solution();
                    let validation = solution.accepted_solution().validation();
                    let inventory_scale = solution
                        .component_moles()
                        .iter()
                        .map(|value| value.abs())
                        .sum::<f64>();
                    let balance_limit = 1e-6 + 1e-6 * inventory_scale;
                    assert_eq!(solution.component_moles().len(), 100);
                    assert!(solution
                        .component_moles()
                        .iter()
                        .all(|value| value.is_finite() && *value > 0.0));
                    assert!(validation.residual_l2_norm.is_finite());
                    assert!(validation.max_abs_element_balance_error <= balance_limit);
                    if index > 0 {
                        assert!(point.report().used_continuation_seed());
                    }
                    max_residual = max_residual.max(validation.residual_l2_norm);
                    max_balance = max_balance.max(validation.max_abs_element_balance_error);
                    let accepted_log_moles = solution.accepted_solution().log_moles();
                    let timing = point.report().timing();
                    let accepted_attempt = solution.solve_report().accepted_attempt();
                    slow_point_rows.push(LiveTemperatureRangeSlowPointRow {
                        backend: name.to_string(),
                        point_index: index.to_string(),
                        temperature: format!("{:.1}", point.report().temperature()),
                        total_ms: live_duration_ms(timing.total()),
                        nonlinear_ms: live_duration_ms(timing.nonlinear_solve()),
                        symbolic_ms: live_duration_ms(timing.symbolic_construction()),
                        numeric_closures_ms: live_duration_ms(
                            timing.numeric_closure_construction(),
                        ),
                        numerical_problem_ms: live_duration_ms(
                            timing.numerical_problem_preparation(),
                        ),
                        validation_ms: live_duration_ms(timing.validation()),
                        seed_delta_log_inf: live_log_coordinate_delta_inf(
                            &previous_accepted_log_moles,
                            accepted_log_moles,
                        ),
                        seed_log_range: live_log_coordinate_range(&previous_accepted_log_moles),
                        result_log_range: live_log_coordinate_range(accepted_log_moles),
                        accepted_backend: format!("{:?}", solution.solve_report().accepted_backend),
                        termination: live_attempt_metric(accepted_attempt, |metrics| {
                            format!("{:?}", metrics.termination)
                        }),
                        iterations: live_attempt_metric(accepted_attempt, |metrics| {
                            metrics.iterations
                        }),
                        residual_evaluations: live_attempt_metric(accepted_attempt, |metrics| {
                            metrics.residual_evaluations
                        }),
                        jacobian_evaluations: live_attempt_metric(accepted_attempt, |metrics| {
                            metrics.jacobian_evaluations
                        }),
                        linear_solves: live_attempt_metric(accepted_attempt, |metrics| {
                            metrics.linear_solves
                        }),
                        backend_ms: live_attempt_metric(accepted_attempt, |metrics| {
                            metrics.elapsed_millis
                        }),
                        residual_ms: live_attempt_evaluation_ms(accepted_attempt, |timing| {
                            timing.residual_evaluation_micros
                        }),
                        jacobian_ms: live_attempt_evaluation_ms(accepted_attempt, |timing| {
                            timing.jacobian_evaluation_micros
                        }),
                        engine_ms: live_attempt_evaluation_ms(accepted_attempt, |timing| {
                            timing.solver_overhead_micros
                        }),
                    });
                    previous_accepted_log_moles = accepted_log_moles.to_vec();
                }
                let min_point_time = range
                    .points()
                    .iter()
                    .map(|point| point.report().timing().total())
                    .min()
                    .unwrap_or_default();
                let worst_point = range
                    .points()
                    .iter()
                    .max_by_key(|point| point.report().timing().total())
                    .expect("successful range must contain a worst point");
                let worst_timing = worst_point.report().timing();
                let worst_solve_report = worst_point.solution().solve_report();
                // Report the complete nonlinear work for the point. This
                // includes every continuation multi-start seed, whereas a
                // maximum over one backend attempt would hide recovery cost.
                let worst_iterations = worst_point.solution().nonlinear_iterations();

                rows.push(LiveTemperatureRangeBackendRow {
                    backend: name.to_string(),
                    status: "OK".to_string(),
                    total_ms: live_duration_ms(range.report().total()),
                    wall_ms: live_duration_ms(started.elapsed()),
                    mean_ms: live_duration_ms(range.report().point_timing().mean()),
                    min_ms: live_duration_ms(min_point_time),
                    median_ms: live_duration_ms(range.report().point_timing().median()),
                    worst_ms: live_duration_ms(range.report().point_timing().worst()),
                    worst_temperature: format!("{:.1}", worst_point.report().temperature()),
                    worst_nonlinear_ms: live_duration_ms(worst_timing.nonlinear_solve()),
                    worst_symbolic_ms: live_duration_ms(worst_timing.symbolic_construction()),
                    worst_accepted_backend: format!("{:?}", worst_solve_report.accepted_backend),
                    worst_attempts: worst_point
                        .solution()
                        .started_backend_attempts()
                        .to_string(),
                    worst_iterations: worst_iterations.to_string(),
                    builds: range.report().formulation_builds().to_string(),
                    initial_setup_ms: live_duration_ms(
                        range.report().initial_formulation_timing().total(),
                    ),
                    initial_symbolic_ms: live_duration_ms(
                        range
                            .report()
                            .initial_formulation_timing()
                            .symbolic_construction(),
                    ),
                    initial_problem_preparation_ms: live_duration_ms(
                        range
                            .report()
                            .initial_formulation_timing()
                            .numerical_problem_preparation(),
                    ),
                    formulation_build_ms: live_duration_ms(
                        range
                            .points()
                            .iter()
                            .map(|point| point.report().formulation_build())
                            .sum(),
                    ),
                    reuses: range.report().formulation_reuses().to_string(),
                    symbolic_updates: range.report().symbolic_parameter_updates().to_string(),
                    max_residual: format!("{max_residual:.3e}"),
                    max_balance: format!("{max_balance:.3e}"),
                    error: "-".to_string(),
                });
            }
            Err(error) => {
                failed_backends.push(name);
                failure_rows.extend(live_temperature_range_failure_rows(name, &error));
                rows.push(LiveTemperatureRangeBackendRow {
                    backend: name.to_string(),
                    status: "FAILED".to_string(),
                    total_ms: "-".to_string(),
                    wall_ms: live_duration_ms(started.elapsed()),
                    mean_ms: "-".to_string(),
                    min_ms: "-".to_string(),
                    median_ms: "-".to_string(),
                    worst_ms: "-".to_string(),
                    worst_temperature: "-".to_string(),
                    worst_nonlinear_ms: "-".to_string(),
                    worst_symbolic_ms: "-".to_string(),
                    worst_accepted_backend: "-".to_string(),
                    worst_attempts: "-".to_string(),
                    worst_iterations: "-".to_string(),
                    builds: "-".to_string(),
                    initial_setup_ms: "-".to_string(),
                    initial_symbolic_ms: "-".to_string(),
                    initial_problem_preparation_ms: "-".to_string(),
                    formulation_build_ms: "-".to_string(),
                    reuses: "-".to_string(),
                    symbolic_updates: "-".to_string(),
                    max_residual: "-".to_string(),
                    max_balance: "-".to_string(),
                    error: format!("{error:?}"),
                });
            }
        }
    }

    assert!(
        successful_backends > 0,
        "the 100-species live temperature-range matrix must have at least one successful backend"
    );
    println!("live real temperature-range backend matrix (species=100, points=50)");
    println!("{}", Table::new(rows).with(Style::rounded()));
    slow_point_rows.sort_by(|left, right| {
        right
            .total_ms
            .parse::<f64>()
            .unwrap_or_default()
            .total_cmp(&left.total_ms.parse::<f64>().unwrap_or_default())
    });
    let mut slowest_per_backend = Vec::new();
    let mut selected_per_backend: HashMap<String, usize> = HashMap::new();
    for row in slow_point_rows {
        let count = selected_per_backend.entry(row.backend.clone()).or_default();
        if *count < 3 {
            *count += 1;
            slowest_per_backend.push(row);
        }
    }
    if !slowest_per_backend.is_empty() {
        println!("three slowest accepted points per successful backend");
        println!("{}", Table::new(slowest_per_backend).with(Style::rounded()));
    }
    if !failure_rows.is_empty() {
        println!("failed first-point/backend attempt diagnostics");
        println!("{}", Table::new(failure_rows).with(Style::rounded()));
    }
    println!(
        "summary: attempted={} successful={} failed={:?}",
        successful_backends + failed_backends.len(),
        successful_backends,
        failed_backends,
    );
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
#[ignore = "release live phase-transition temperature-range characterization"]
fn live_bounded_phase_control_temperature_range_records_real_ice_transition() {
    let before = live_library_file_snapshot();
    let range = PhaseEquilibriumPipelineRequest::new(
        live_gas_solid_water_spec(),
        vec![0.5, 0.25, 0.0],
        EquilibriumConditions::new(250.0, 101_325.0, 101_325.0)
            .expect("live ice range conditions must validate"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(
        EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
    )
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve_temperature_range(
        TemperatureGrid::new(vec![250.0, 260.0, 270.0]).expect("live ice range grid must validate"),
    )
    .expect("live bounded ice range must solve");

    assert!(range.report().phase_control_transitions() > 0);
    assert!(range.report().phase_projection_cache_entries() >= 2);
    assert!(range.report().phase_prepared_cache_entries() >= 2);
    assert!(range.report().phase_rst_cache_entries() >= 2);
    assert!(range.report().initial_formulation_timing().enabled());
    assert!(
        range
            .points()
            .iter()
            .any(|point| point.report().formulation_build() > std::time::Duration::ZERO),
        "the live ice transition must account for at least one reduced formulation build"
    );
    assert!(range.points()[1].report().phase_set_reused());
    assert!(range.points()[1..]
        .iter()
        .any(|point| point.report().symbolic_parameter_reused()));
    let cache_timings = range
        .points()
        .last()
        .expect("real ice range must publish a final point")
        .report()
        .formulation_cache_timings();
    assert_eq!(
        cache_timings.len(),
        range.report().phase_prepared_cache_entries(),
        "the final point must retain one timing snapshot per prepared active-set cache entry"
    );
    assert!(
        cache_timings
            .iter()
            .any(|entry| entry.formulation_build() > std::time::Duration::ZERO),
        "at least one real ice active-set cache entry must have measured build time"
    );
    for point in range.points() {
        let solution = point.solution();
        let validation = solution.accepted_solution().validation();
        let scale = solution
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(solution
            .component_moles()
            .iter()
            .all(|moles| moles.is_finite()));
        assert!(validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * scale);
    }
    assert!(
        range
            .points()
            .iter()
            .flat_map(|point| {
                point
                    .solution()
                    .phase_control_report()
                    .into_iter()
                    .flat_map(|report| report.transitions.iter())
            })
            .all(|transition| transition.transition_duration > std::time::Duration::ZERO),
        "real ice T-range transitions must retain measured control-pass durations"
    );
    assert_eq!(before, live_library_file_snapshot());

    println!(
        "live bounded ice T-range: points={} transitions={} projections={} prepared={} rst_symbolic={} total={:?}",
        range.points().len(),
        range.report().phase_control_transitions(),
        range.report().phase_projection_cache_entries(),
        range.report().phase_prepared_cache_entries(),
        range.report().phase_rst_cache_entries(),
        range.report().total(),
    );
}

/// One row in the release-oriented real phase-transition matrix.
///
/// The matrix intentionally reports evidence instead of asserting a fixed
/// number of transitions. Database updates may move a coexistence boundary,
/// while the physical contracts below remain stable: finite non-negative
/// amounts, conservation, complementarity, and immutable source files.
#[derive(Debug, Tabled)]
struct LivePhaseTransitionMatrixRow {
    #[tabled(rename = "Fixture")]
    fixture: String,
    #[tabled(rename = "T K")]
    temperature: String,
    #[tabled(rename = "Active phases")]
    active_phases: String,
    #[tabled(rename = "Transitions")]
    transitions: String,
    #[tabled(rename = "Max balance")]
    max_balance: String,
    #[tabled(rename = "Complementarity")]
    complementarity: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
}

#[test]
#[ignore = "release real-data phase-transition matrix"]
fn live_real_phase_transition_release_matrix() {
    let before = live_library_file_snapshot();
    let cases = vec![
        (
            "water-gas-ice",
            250.0,
            live_gas_solid_water_spec(),
            vec![0.5, 0.25, 0.0],
        ),
        (
            "water-gas-liquid",
            350.0,
            live_multiphase_spec(),
            vec![0.5, 0.25, 0.0],
        ),
        (
            "water-gas-hot",
            550.0,
            live_multiphase_spec(),
            vec![0.5, 0.25, 0.0],
        ),
        (
            "carbon-graphite",
            700.0,
            live_gas_solid_carbon_spec(),
            vec![1.0, 0.1, 0.0],
        ),
        (
            "carbon-hot",
            1_400.0,
            live_gas_solid_carbon_spec(),
            vec![1.0, 0.1, 0.0],
        ),
    ];
    let mut rows = Vec::with_capacity(cases.len());

    for (fixture, temperature, spec, composition) in cases {
        let outcome = PhaseEquilibriumPipelineRequest::new(
            spec,
            composition,
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0)
                .expect("real phase matrix conditions must validate"),
        )
        .with_repository(live_repository())
        .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .unwrap_or_else(|error| {
            panic!("real phase matrix case {fixture} at {temperature} K failed: {error:?}")
        });

        let solution = outcome.solution();
        let validation = solution.accepted_solution().validation();
        let mole_scale = solution
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(
            solution
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0),
            "real phase matrix case {fixture} produced invalid physical moles"
        );
        assert!(
            validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale,
            "real phase matrix case {fixture} violated conservation: {}",
            validation.max_abs_element_balance_error
        );
        let acceptance = solution
            .acceptance_report()
            .expect("real phase matrix must publish acceptance evidence");
        assert!(
            acceptance.complementarity.satisfied,
            "real phase matrix case {fixture} must satisfy complementarity"
        );
        let active_phases = solution
            .phases()
            .iter()
            .filter(|phase| solution.phase_status(phase.id()) == Some(PhaseStatus::Active))
            .map(|phase| format!("{:?}", phase.id()))
            .collect::<Vec<_>>()
            .join(", ");
        let transitions = solution.phase_control_transitions();
        rows.push(LivePhaseTransitionMatrixRow {
            fixture: fixture.to_string(),
            temperature: format!("{temperature:.1}"),
            active_phases,
            transitions: transitions.to_string(),
            max_balance: format!("{:.3e}", validation.max_abs_element_balance_error),
            complementarity: "OK".to_string(),
            total_ms: live_duration_ms(outcome.timing_report().total()),
        });
    }

    assert_eq!(before, live_library_file_snapshot());
    println!(
        "live real phase-transition matrix (water/ice, water/liquid, graphite)\n{}",
        Table::new(rows).with(Style::rounded())
    );
}

#[test]
fn live_fixed_point_matches_one_point_typed_temperature_range() {
    let (selected, resolved) = live_large_element_limited_gas_resolved();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("real resolved phases must produce a canonical range layout");
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; selected.len()])
        .expect("real fixed/range comparison composition must validate");
    let conditions = EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0)
        .expect("real fixed/range comparison conditions must validate");
    let options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
        .expect("single legacy NR policy must validate");

    let fixed = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, composition.clone())
            .with_solve_options(options.clone()),
    )
    .expect("real fixed point must solve");
    let range = TemperatureRangeRequest::new(
        &resolved,
        composition,
        conditions.pressure(),
        conditions.reference_pressure(),
        TemperatureGrid::new(vec![conditions.temperature()]).expect("one-point grid must validate"),
    )
    .expect("real one-point range request must validate")
    .with_solve_options(options.clone())
    .solve()
    .expect("real one-point range must solve");

    let ranged = range.points()[0].solution();
    assert_eq!(
        fixed.component_moles().len(),
        ranged.component_moles().len()
    );
    for (fixed_moles, ranged_moles) in fixed.component_moles().iter().zip(ranged.component_moles())
    {
        let scale = fixed_moles.abs().max(ranged_moles.abs()).max(1e-30);
        assert!((fixed_moles - ranged_moles).abs() / scale <= 1e-7);
    }
    assert!(range.report().total() >= range.report().point_timing().total());
}

#[test]
#[ignore = "explicit release temperature-range characterization run"]
fn live_large_element_limited_legacy_temperature_range_story() {
    // Keep the compatibility baseline while the typed facade is characterized
    // by the release-oriented story below.
    let (selected, _) = live_large_element_limited_gas_request();
    let initial_moles = vec![1e-3; selected.len()];
    let temperatures: Vec<f64> = (0..5).map(|index| 1_000.0 + 50.0 * index as f64).collect();

    for solver in [Solvers::LM, Solvers::NR, Solvers::TR] {
        let started = std::time::Instant::now();
        #[allow(deprecated)]
        let result = gas_solver_for_T_range(
            selected.clone(),
            initial_moles.clone(),
            101_325.0,
            temperatures[0],
            // The compatibility sweep treats T_end as an exclusive upper
            // bound; use the next grid point to request all five values.
            temperatures.last().unwrap() + 50.0,
            50.0,
            solver,
            None,
            false,
        )
        .unwrap_or_else(|error| panic!("real {solver:?} temperature sweep failed: {error:?}"));

        assert!(
            result.list_of_failed_T.is_empty(),
            "{solver:?} failed points: {:?}",
            result.list_of_failed_T
        );
        assert_eq!(result.moles_for_T_range.len(), temperatures.len());
        for ((actual_temperature, moles), expected_temperature) in
            result.moles_for_T_range.iter().zip(&temperatures)
        {
            assert_eq!(*actual_temperature, *expected_temperature);
            assert_eq!(moles.len(), selected.len());
            assert!(moles.iter().all(|value| value.is_finite() && *value > 0.0));
        }
        assert_eq!(
            result.temperature_solutions.len(),
            temperatures.len(),
            "every accepted compatibility point must retain its solver report"
        );
        assert!(
            result.temperature_solutions.iter().all(|snapshot| matches!(
                snapshot.solve_report.accepted_backend,
                SolverBackend::Legacy(_)
            )),
            "legacy temperature-range compatibility must not silently promote to an RST backend"
        );
        println!(
            "live real T-range: backend={solver:?}, species={}, points={}, elapsed={:?}",
            selected.len(),
            temperatures.len(),
            started.elapsed()
        );
    }
}

#[test]
#[ignore = "explicit release typed temperature-range characterization run"]
fn live_large_element_limited_typed_temperature_range_story() {
    let (selected, resolved) = live_large_element_limited_gas_resolved();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("resolved live phases must produce a canonical layout");
    let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; selected.len()])
        .expect("live initial composition must match the resolved layout");
    let options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
        .expect("legacy characterization policy must validate")
        .with_timing_mode(EquilibriumTimingMode::Enabled);

    for (grid, expected_direction) in [
        (
            TemperatureGrid::new(vec![1_000.0, 1_050.0, 1_100.0]).unwrap(),
            TemperatureRangeDirection::Ascending,
        ),
        (
            TemperatureGrid::new(vec![1_100.0, 1_050.0, 1_000.0]).unwrap(),
            TemperatureRangeDirection::Descending,
        ),
    ] {
        let result = TemperatureRangeRequest::new(
            &resolved,
            composition.clone(),
            101_325.0,
            101_325.0,
            grid,
        )
        .unwrap()
        .with_solve_options(options.clone())
        .solve()
        .expect("typed real temperature range must solve");

        assert_eq!(result.report().direction(), expected_direction);
        assert_eq!(result.report().point_count(), 3);
        assert_eq!(result.report().formulation_builds(), 1);
        assert_eq!(result.report().formulation_reuses(), 2);
        assert_eq!(result.report().symbolic_parameter_updates(), 0);
        assert!(!result.report().symbolic_problem_reused());
        assert!(result.points().iter().all(|point| {
            point
                .solution()
                .component_moles()
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        }));
        assert!(matches!(
            result.points()[0].report().preparation(),
            TemperatureRangePointPreparation::InitialFormulation
        ));
        let postprocessed = postprocess_temperature_range_solution(
            &result,
            &TemperaturePostprocessingPolicy::default(),
        )
        .expect("typed real temperature range should feed raw postprocessing");
        assert_eq!(
            postprocessed.raw.temperatures(),
            &[1_000.0, 1_050.0, 1_100.0]
        );
        assert_eq!(
            postprocessed.raw.labels().len(),
            result.points()[0].solution().component_moles().len()
        );
        assert!(result.points()[1..]
            .iter()
            .all(|point| point.report().used_continuation_seed()));
        let timing = result.report().point_timing();
        println!(
            "live typed T-range: backend=legacy-nr direction={expected_direction:?} \
             species={} points={} setup_builds={} reuses={} \
             point_total={:?} point_mean={:?} point_median={:?} point_worst={:?}",
            selected.len(),
            result.report().point_count(),
            result.report().formulation_builds(),
            result.report().formulation_reuses(),
            timing.total(),
            timing.mean(),
            timing.median(),
            timing.worst(),
        );
    }

    // The default resolved-data policy is symbolic/RST-first. This separate
    // run proves that the range template updates the shared `T` parameter
    // rather than preparing a new symbolic problem for every point.
    let rst_result = TemperatureRangeRequest::new(
        &resolved,
        composition,
        101_325.0,
        101_325.0,
        TemperatureGrid::new(vec![1_000.0, 1_050.0, 1_100.0]).unwrap(),
    )
    .unwrap()
    .with_solve_options(
        EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
    )
    .solve()
    .expect("typed real RST temperature range must solve");
    assert!(rst_result.report().symbolic_problem_reused());
    let point_reuse_count = rst_result
        .points()
        .iter()
        .filter(|point| point.report().symbolic_parameter_reused())
        .count();
    // A real NASA interval boundary may require one symbolic rebuild. The
    // range contract is therefore point-level reuse evidence, not a promise
    // that every temperature uses the first interval's captured constants.
    assert_eq!(
        rst_result.report().symbolic_parameter_updates(),
        point_reuse_count
    );
    assert!(point_reuse_count > 0);
    let timing = rst_result.report().point_timing();
    println!(
        "live typed T-range: backend=rst-default direction=Ascending species={} \
         points={} setup_builds={} reuses={} symbolic_reuses={} \
         point_total={:?} point_mean={:?} point_median={:?} point_worst={:?}",
        selected.len(),
        rst_result.report().point_count(),
        rst_result.report().formulation_builds(),
        rst_result.report().formulation_reuses(),
        point_reuse_count,
        timing.total(),
        timing.mean(),
        timing.median(),
        timing.worst(),
    );
}

#[test]
fn live_low_temperature_water_activates_the_real_ice_phase() {
    let request = PhaseEquilibriumPipelineRequest::new(
        live_gas_solid_water_spec(),
        vec![0.5, 0.25, 0.0],
        EquilibriumConditions::new(250.0, 101_325.0, 101_325.0)
            .expect("ice fixture conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(PhaseControlPolicy::default());

    let outcome = request
        .solve()
        .expect("real low-temperature water fixture must solve");
    let solid = PhaseId::new(Some("solid".to_string()));
    let gas = PhaseId::new(Some("gas".to_string()));
    let solid_moles = outcome
        .solution()
        .phase_total(&solid)
        .expect("solid phase must remain represented");
    let gas_moles = outcome
        .solution()
        .phase_total(&gas)
        .expect("gas phase must remain represented");

    assert!(
        solid_moles > 0.49,
        "ice should contain almost all water at 250 K and 1 atm; solid={solid_moles:e}, gas={gas_moles:e}"
    );
    assert!(gas_moles >= 0.25);
    assert!(outcome
        .solution()
        .phase_control_report()
        .expect("bounded solve must retain its phase-control report")
        .transitions
        .iter()
        .any(|transition| !transition.activated.is_empty()));
    assert!(
        outcome
            .solution()
            .phase_control_report()
            .expect("bounded solve must retain its phase-control report")
            .transitions
            .iter()
            .all(|transition| transition.transition_duration > std::time::Duration::ZERO),
        "real ice activation records must retain a measured control-pass duration"
    );
}

/// Release-only inverse story for a real gas/ice system.
///
/// This keeps a real phase-transition story at the P,H boundary as well as
/// at P,T: the reference state is solved first, its additive enthalpy is
/// reused as the target, and `Auto` must retain the rejected monolithic
/// candidate before publishing the accepted nested recovery.
#[test]
#[ignore = "release-oriented real-data ice P,H Auto story"]
fn live_ice_pt_to_h_to_auto_ph_preserves_phase_transition() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_gas_solid_water_spec(),
        live_repository(),
    )
    .expect("real gas/ice system must resolve");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("real gas/ice layout must remain valid");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0])
        .expect("gas/ice composition must match the resolved layout");
    let pressure = 101_325.0;
    let temperature = 250.0;
    let pt = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(temperature, pressure, pressure)
                .expect("ice P,T conditions must be valid"),
            initial.clone(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("real gas/ice P,T reference must solve");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("gas and ice records must provide P,H thermochemistry");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(pt.component_moles(), temperature)
        .expect("ice reference enthalpy must be finite");
    let constraint = EquilibriumConstraint::ph_joules(
        pressure,
        pressure,
        TotalEnthalpyJoules::new(target_enthalpy).expect("ice target must be finite"),
        temperature,
    )
    .expect("ice P,H constraint must be valid");
    let bounds = TemperatureBounds::new(245.0, 255.0).expect("ice P,H bounds must be valid");

    let monolithic_error = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial.clone(),
            constraint.clone(),
            bounds.clone(),
            thermochemistry.clone(),
        )
        .expect("real gas/ice monolithic request must validate")
        .with_ph_solve_mode(PhSolveMode::Monolithic)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect_err("the direct monolithic ice route must retain its typed failure");
    assert_eq!(
        monolithic_error.kind(),
        ReactionExtentErrorKind::AllBackendsFailed,
        "the reference Auto fallback contract must not hide the direct route failure"
    );

    let nested_started = std::time::Instant::now();
    let nested = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial.clone(),
            constraint.clone(),
            bounds.clone(),
            thermochemistry.clone(),
        )
        .expect("real gas/ice nested request must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("the explicit nested ice route must recover the phase-controlled state");

    let auto_started = std::time::Instant::now();
    let ph = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            constraint,
            bounds,
            thermochemistry,
        )
        .expect("real gas/ice P,H request must validate")
        .with_ph_solve_mode(PhSolveMode::Auto)
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("ice P,H Auto must recover through the nested phase-control path");
    let solid = PhaseId::new(Some("solid".to_string()));
    assert_eq!(nested.report().solve_path(), PhSolvePath::NestedTemperature);
    assert!(nested.report().fallback_reason().is_none());
    assert!((ph.temperature() - temperature).abs() < 1e-4);
    assert!(
        (ph.temperature() - nested.temperature()).abs() <= 1.0e-4,
        "Auto and explicit nested P,H routes recovered different temperatures"
    );
    assert!(ph.equilibrium().phase_total(&solid).unwrap_or_default() > 0.49);
    assert_eq!(ph.report().solve_path(), PhSolvePath::NestedTemperature);
    assert_eq!(
        ph.report()
            .fallback_reason()
            .expect("ice P,H Auto must retain monolithic fallback evidence")
            .error_kind(),
        crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::
            ReactionExtentErrorKind::AllBackendsFailed
    );
    for solution in [&nested, &ph] {
        assert!(solution
            .equilibrium()
            .phase_control_report()
            .expect("ice P,H result must retain phase-control evidence")
            .transitions
            .iter()
            .any(|transition| transition.activated.iter().any(|phase| phase.index() == 1)));
        assert!(solution.report().trials().iter().all(|trial| {
            trial.inner_evidence().is_some_and(|evidence| {
                evidence
                    .phase_control_report()
                    .is_some_and(|report| !report.nonlinear_reports.is_empty())
                    && evidence.acceptance_report().is_some()
            })
        }));
        assert!(solution.enthalpy_error().abs() <= solution.enthalpy_error_limit_joules());
    }
    for (component, (nested_moles, auto_moles)) in nested
        .equilibrium()
        .component_moles()
        .iter()
        .zip(ph.equilibrium().component_moles())
        .enumerate()
    {
        let tolerance = 1.0e-8 + 1.0e-6 * nested_moles.abs().max(auto_moles.abs());
        assert!(
            (nested_moles - auto_moles).abs() <= tolerance,
            "component {component} differs between explicit nested and Auto recovery: nested={nested_moles:e}, auto={auto_moles:e}, tolerance={tolerance:e}"
        );
    }
    assert_eq!(before, live_library_file_snapshot());
    let phase_control = ph
        .equilibrium()
        .phase_control_report()
        .expect("ice P,H result must retain phase-control evidence");
    println!(
        "live real ice P,H route matrix: monolithic=FAILED({:?}) auto_fallback={:?} solid_moles={:.6e} transitions={}\n{}",
        monolithic_error.kind(),
        ph.report().fallback_reason().map(|reason| format!("{:?}", reason.error_kind())),
        ph.equilibrium().phase_total(&solid).unwrap_or_default(),
        phase_control.transitions.len(),
        Table::new([
            live_ph_path_row("nested-phase-control", &nested, nested_started.elapsed()),
            live_ph_path_row("auto-nested-recovery", &ph, auto_started.elapsed()),
        ])
        .with(Style::rounded())
    );
}

#[test]
fn live_high_temperature_water_keeps_inventory_free_liquid_inactive() {
    // AllCandidatePhases means "consider every phase", not "force every phase
    // with zero physical inventory into the positive log-moles solve". At
    // 550 K and one atmosphere the liquid candidate is absent from the input;
    // it must remain inactive while still being available for stability tests.
    let phase_policy = PhaseControlPolicy::default()
        .with_initial_phase_set(InitialPhaseSet::AllCandidatePhases)
        .expect("all-candidate initial policy must be valid");
    let outcome = PhaseEquilibriumPipelineRequest::new(
        live_multiphase_spec(),
        vec![0.5, 0.25, 0.0],
        EquilibriumConditions::new(550.0, 101_325.0, 101_325.0)
            .expect("high-temperature water conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(phase_policy)
    .solve()
    .expect("real high-temperature water lifecycle must solve");

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let report = outcome
        .solution()
        .phase_control_report()
        .expect("bounded solve must publish lifecycle evidence");
    let liquid_index =
        crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex::new(1, 2).unwrap();

    assert_eq!(
        report.initial_phase_set.status(liquid_index),
        PhaseStatus::Inactive
    );
    assert_eq!(
        report.final_phase_set.status(liquid_index),
        PhaseStatus::Inactive
    );
    assert_eq!(outcome.solution().phase_total(&liquid), Some(0.0));
    assert!(report.transitions.is_empty());
}

#[test]
fn live_high_temperature_water_recovers_positive_liquid_disappearance() {
    // This is the complementary boundary case: liquid water is present in the
    // input, but at 550 K and one atmosphere no positive liquid equilibrium is
    // expected. The outer loop must solve the gas-only boundary and publish a
    // typed deactivation record instead of treating the interior failure as a
    // fatal backend error.
    let phase_policy = PhaseControlPolicy::default()
        .with_initial_phase_set(InitialPhaseSet::AllCandidatePhases)
        .expect("all-candidate initial policy must be valid");
    let outcome = PhaseEquilibriumPipelineRequest::new(
        live_multiphase_spec(),
        vec![0.5, 0.25, 0.1],
        EquilibriumConditions::new(550.0, 101_325.0, 101_325.0)
            .expect("high-temperature water conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(phase_policy)
    .solve()
    .expect("real high-temperature water boundary recovery must solve");

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let liquid_index =
        crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex::new(1, 2).unwrap();
    let report = outcome
        .solution()
        .phase_control_report()
        .expect("bounded solve must publish lifecycle evidence");

    assert_eq!(
        report.final_phase_set.status(liquid_index),
        PhaseStatus::Inactive
    );
    assert_eq!(outcome.solution().phase_total(&liquid), Some(0.0));
    assert!(report.transitions.iter().any(|transition| {
        transition.deactivated.iter().any(|phase| phase.index() == 1)
            && matches!(
                transition.reason,
                crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseTransitionReason::BoundaryUnstableActivePhase { .. }
            )
    }));
    assert!(
        report
            .transitions
            .iter()
            .all(|transition| transition.transition_duration > std::time::Duration::ZERO),
        "real phase deactivation records must retain a measured control-pass duration"
    );
}

#[test]
fn live_water_hysteresis_retains_a_marginal_active_liquid_phase() {
    // The broad dead band is deliberate: it makes the test exercise the
    // retention rule itself rather than relying on a particular near-boiling
    // database value. A live NASA liquid record is still used end to end.
    let phase_policy = PhaseControlPolicy::with_explicit_hysteresis(0.05, -1e12, 1e12)
        .expect("finite live hysteresis policy must be valid");
    let outcome = PhaseEquilibriumPipelineRequest::new(
        live_multiphase_spec(),
        vec![0.5, 0.25, 0.1],
        EquilibriumConditions::new(363.0, 101_325.0, 101_325.0)
            .expect("live hysteresis conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(phase_policy)
    .solve()
    .expect("live hysteresis fixture must solve");

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let liquid_index =
        crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex::new(1, 2).unwrap();
    let report = outcome
        .solution()
        .phase_control_report()
        .expect("live hysteresis solve must publish phase-control evidence");
    assert_eq!(
        report.final_phase_set.status(liquid_index),
        PhaseStatus::Active
    );
    assert!(outcome.solution().phase_total(&liquid).unwrap() < 0.05);
    assert_eq!(report.iterations, 1);
    assert!(report.transitions.is_empty());
}

#[test]
fn live_water_budget_failure_rolls_back_without_publishing_partial_solution() {
    // Ice appearance requires a restart. A one-pass budget must therefore
    // fail after the first accepted candidate, while the public transaction
    // returns no partial solution and the local thermochemistry files remain
    // untouched.
    let before = live_library_file_snapshot();
    let phase_policy = PhaseControlPolicy::default()
        .with_max_phase_iterations(1)
        .expect("positive phase-control budget must be valid");
    let result = PhaseEquilibriumPipelineRequest::new(
        live_gas_solid_water_spec(),
        vec![0.5, 0.25, 0.0],
        EquilibriumConditions::new(250.0, 101_325.0, 101_325.0)
            .expect("live budget conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(phase_policy)
    .solve();

    assert!(matches!(
        result,
        Err(PhaseEquilibriumPipelineError::Solve(
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError::PhaseControlDidNotConverge { iterations: 1 }
        ))
    ));
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
fn live_boudouard_system_activates_real_graphite_from_co_co2_gas() {
    let low = PhaseEquilibriumPipelineRequest::new(
        live_gas_solid_carbon_spec(),
        // Both gas chemical potentials are finite at the initial active set.
        // This isolates graphite stability from arbitrary trace-floor values.
        vec![1.0, 0.1, 0.0],
        EquilibriumConditions::new(700.0, 101_325.0, 101_325.0)
            .expect("Boudouard fixture conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve()
    .expect("real Boudouard fixture must solve");
    let high = PhaseEquilibriumPipelineRequest::new(
        live_gas_solid_carbon_spec(),
        vec![1.0, 0.1, 0.0],
        EquilibriumConditions::new(1_400.0, 101_325.0, 101_325.0)
            .expect("high-temperature Boudouard conditions must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve()
    .expect("high-temperature Boudouard fixture must solve");

    let solid = PhaseId::new(Some("solid".to_string()));
    let graphite_moles = low
        .solution()
        .phase_total(&solid)
        .expect("graphite phase must remain represented");
    let high_temperature_graphite = high
        .solution()
        .phase_total(&solid)
        .expect("high-temperature graphite candidate must remain represented");
    assert!(
        graphite_moles > 1e-3,
        "graphite should be stable for the 700 K CO/CO2 fixture; moles={graphite_moles:e}"
    );
    assert!(
        graphite_moles > high_temperature_graphite,
        "graphite inventory must decrease across the 700 K -> 1400 K Boudouard shift; low={graphite_moles:e}, high={high_temperature_graphite:e}"
    );
    assert!(low
        .solution()
        .phase_control_report()
        .expect("bounded solve must retain its phase-control report")
        .transitions
        .iter()
        .any(|transition| !transition.activated.is_empty()));
    assert!(low
        .solution()
        .build_report()
        .components()
        .iter()
        .any(|component| {
            component.component().label() == "solid::C(gr)"
                && component.thermo_source().library() == "NASA_cond"
                && component.thermo_source().record_key() == "C(gr)"
        }));
}

#[test]
fn live_water_pair_shows_temperature_driven_phase_dominance_shift() {
    let low = PhaseEquilibriumPipelineRequest::new(
        live_multiphase_spec(),
        vec![0.5, 0.25, 0.0],
        EquilibriumConditions::new(350.0, 101_325.0, 101_325.0)
            .expect("low temperature must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve()
    .expect("low-temperature live water pair must solve");

    let high = PhaseEquilibriumPipelineRequest::new(
        live_multiphase_spec(),
        vec![0.5, 0.25, 0.0],
        EquilibriumConditions::new(550.0, 101_325.0, 101_325.0)
            .expect("high temperature must be valid"),
    )
    .with_repository(live_repository())
    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1e-30 })
    .with_solve_options(EquilibriumSolveOptions::default())
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve()
    .expect("high-temperature live water pair must solve");

    let liquid = PhaseId::new(Some("liquid".to_string()));
    let gas = PhaseId::new(Some("gas".to_string()));

    assert!(
        low.solution().phase_total(&liquid).unwrap()
            >= high.solution().phase_total(&liquid).unwrap()
    );
    assert!(
        low.solution().phase_total(&liquid).unwrap() > 1e-3,
        "the 350 K fixture must exercise a real liquid-water inventory"
    );
    assert!(
        high.solution().phase_total(&gas).unwrap() >= low.solution().phase_total(&gas).unwrap()
    );
    assert_eq!(
        low.solution().layout_fingerprint(),
        high.solution().layout_fingerprint()
    );
    assert_eq!(
        low.solution().metadata().components().len(),
        high.solution().metadata().components().len()
    );
    assert!(low
        .solution()
        .summary_rows()
        .iter()
        .any(|row| row.section == "phase_control"));
    assert!(high
        .solution()
        .summary_rows()
        .iter()
        .any(|row| row.section == "phase_control"));
}

#[test]
#[ignore = "release-oriented real-data P,H continuation story"]
fn live_reactive_gas_ph_target_range_continues_from_previous_accepted_state() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas must resolve for the P,H target-range story");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("reactive gas layout must validate");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])
        .expect("reactive gas initial composition must validate");
    let pressure = 101_325.0;
    let temperatures = [2_300.0, 2_500.0, 2_700.0];
    let options = EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled);
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live NASA records must provide a thermochemistry bundle");

    let mut target_pairs = temperatures
        .into_iter()
        .map(|temperature| {
            let conditions = EquilibriumConditions::new(temperature, pressure, pressure)
                .expect("reference P,T conditions must validate");
            let solution = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, initial.clone())
                    .with_solve_options(options.clone()),
            )
            .expect("reference P,T state must solve");
            let target = thermochemistry
                .enthalpy_model()
                .evaluate_total(solution.component_moles(), temperature)
                .expect("reference state must have finite total enthalpy");
            (target, temperature)
        })
        .collect::<Vec<_>>();
    target_pairs.sort_by(|left, right| left.0.total_cmp(&right.0));
    let targets = PhEnthalpyGrid::new(target_pairs.iter().map(|pair| pair.0).collect())
        .expect("reference enthalpies must form a strictly monotone grid");
    assert_eq!(targets.direction(), PhRangeDirection::Ascending);

    let range = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial.clone(),
        pressure,
        pressure,
        targets.clone(),
        TemperatureBounds::new(2_100.0, 2_900.0).expect("P,H bounds must validate"),
        target_pairs[0].1,
        thermochemistry.clone(),
    )
    .expect("typed P,H range request must validate")
    .with_solve_options(options.clone())
    .solve()
    .expect("typed P,H range must solve all targets");

    assert_eq!(range.points().len(), target_pairs.len());
    assert_eq!(range.report().continuation_points(), 2);
    assert_eq!(range.report().direction(), PhRangeDirection::Ascending);
    assert_eq!(range.report().formulation_builds(), 1);
    assert_eq!(range.report().formulation_reuses(), 2);
    for (index, (point, (_, expected_temperature))) in
        range.points().iter().zip(target_pairs.iter()).enumerate()
    {
        let point_report = point.report();
        assert_eq!(point_report.index(), index);
        assert_eq!(
            point_report.preparation(),
            if index == 0 {
                PhRangePointPreparation::Initial
            } else {
                PhRangePointPreparation::Continued
            }
        );
        if index > 0 {
            assert_eq!(
                point_report.seed_temperature(),
                range.points()[index - 1].report().solved_temperature()
            );
        }
        assert!(
            (point_report.solved_temperature() - expected_temperature).abs() < 1.0,
            "target {} recovered {} K instead of {} K",
            index,
            point_report.solved_temperature(),
            expected_temperature
        );
        let solution = point.solution().equilibrium();
        let validation = solution.accepted_solution().validation();
        let mole_scale = solution
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(
            validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale,
            "P,H range point {index} violated conservation: {}",
            validation.max_abs_element_balance_error
        );
    }

    // The nested reference route has its own scalar bracket for every target,
    // but fixed-declared inner P,T structure must still be prepared only once
    // for the complete range. This is deliberately separate from bounded
    // phase control, where an active-set lifecycle cannot be shared blindly.
    let nested_options = options
        .clone()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
        .expect("nested reference backend policy must validate");
    let nested = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial.clone(),
        pressure,
        pressure,
        targets.clone(),
        TemperatureBounds::new(2_100.0, 2_900.0).expect("nested P,H bounds must validate"),
        target_pairs[0].1,
        thermochemistry.clone(),
    )
    .expect("nested typed P,H range request must validate")
    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    .with_solve_options(nested_options)
    .solve()
    .expect("nested typed P,H range must solve all targets");
    assert_eq!(nested.points().len(), target_pairs.len());
    assert_eq!(nested.report().formulation_builds(), 1);
    assert!(
        nested.report().formulation_reuses() >= target_pairs.len() - 1,
        "nested P,H range did not reuse its fixed P,T template"
    );

    let descending_pairs = target_pairs.iter().rev().copied().collect::<Vec<_>>();
    let descending = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial,
        pressure,
        pressure,
        PhEnthalpyGrid::new(descending_pairs.iter().map(|pair| pair.0).collect())
            .expect("descending reference enthalpies must validate"),
        TemperatureBounds::new(2_100.0, 2_900.0).expect("P,H bounds must validate"),
        descending_pairs[0].1,
        thermochemistry,
    )
    .expect("descending typed P,H range request must validate")
    .with_solve_options(options)
    .solve()
    .expect("descending typed P,H range must solve all targets");
    assert_eq!(
        descending.report().direction(),
        PhRangeDirection::Descending
    );
    assert_eq!(descending.report().continuation_points(), 2);
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
                descending.points()[index - 1].report().solved_temperature()
            );
        }
    }

    println!(
        "P,H target range | direction=Ascending | points={} | formulation_builds={} | formulation_reuses={}",
        range.points().len(),
        range.report().formulation_builds(),
        range.report().formulation_reuses(),
    );
    println!("index | target J | seed K | solved K | preparation | elapsed");
    for point in range.points() {
        let report = point.report();
        println!(
            "{} | {:.6e} | {:.3} | {:.3} | {:?} | {:?}",
            report.index(),
            report.target_enthalpy_joules(),
            report.seed_temperature(),
            report.solved_temperature(),
            report.preparation(),
            report.elapsed(),
        );
    }
    let nested_timing = nested.report().point_timing();
    println!(
        "nested | points={} | formulation_builds={} | formulation_reuses={} | total={:?} | mean={:?} | worst={:?}",
        nested.points().len(),
        nested.report().formulation_builds(),
        nested.report().formulation_reuses(),
        nested_timing.total(),
        nested_timing.mean(),
        nested_timing.worst(),
    );
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
#[ignore = "heavy release real-data P,H target-range backend matrix"]
fn live_reactive_gas_ph_target_range_backend_matrix() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas must resolve for the P,H range backend matrix");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("reactive gas layout must validate");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])
        .expect("reactive gas initial composition must validate");
    let pressure = 101_325.0;
    let reference_options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
        .expect("reference backend policy must validate");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live NASA records must provide a thermochemistry bundle");
    let target_temperatures = [2_300.0, 2_500.0, 2_700.0];
    let mut target_pairs = target_temperatures
        .into_iter()
        .map(|temperature| {
            let conditions = EquilibriumConditions::new(temperature, pressure, pressure)
                .expect("reference P,T conditions must validate");
            let reference = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, initial.clone())
                    .with_solve_options(reference_options.clone()),
            )
            .expect("reference P,T state must solve");
            let target = thermochemistry
                .enthalpy_model()
                .evaluate_total(reference.component_moles(), temperature)
                .expect("reference state must have finite total enthalpy");
            (target, temperature)
        })
        .collect::<Vec<_>>();
    target_pairs.sort_by(|left, right| left.0.total_cmp(&right.0));

    let backends = [
        (
            "rst_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ),
        (
            "rst_minpack_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
        ),
        (
            "rst_trust_region_lm",
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
        ),
        ("legacy_lm", SolverBackend::Legacy(Solvers::LM)),
        ("legacy_nr", SolverBackend::Legacy(Solvers::NR)),
        ("legacy_tr", SolverBackend::Legacy(Solvers::TR)),
    ];
    let mut successful = 0usize;
    let mut rows = Vec::with_capacity(backends.len());
    for (name, backend) in backends {
        let started = std::time::Instant::now();
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(backend))
            .expect("single backend policy must validate")
            .with_timing_mode(EquilibriumTimingMode::Enabled);
        let request = PhRangeRequest::from_resolved_thermochemistry(
            &resolved,
            initial.clone(),
            pressure,
            pressure,
            PhEnthalpyGrid::new(target_pairs.iter().map(|pair| pair.0).collect())
                .expect("target enthalpies must form a monotone grid"),
            TemperatureBounds::new(2_100.0, 2_900.0).expect("P,H bounds must validate"),
            target_pairs[0].1,
            thermochemistry.clone(),
        )
        .expect("typed P,H range request must validate")
        .with_solve_options(options);

        match request.solve() {
            Ok(solution) => {
                successful += 1;
                let timing = solution.report().point_timing();
                assert_eq!(solution.points().len(), target_pairs.len());
                assert_eq!(solution.report().formulation_builds(), 1);
                assert_eq!(solution.report().formulation_reuses(), 2);
                let mut max_residual = 0.0_f64;
                let mut max_balance = 0.0_f64;
                let mut max_enthalpy_error = 0.0_f64;
                let mut total_attempts = 0usize;
                let mut accepted_backends = std::collections::BTreeSet::new();
                for point in solution.points() {
                    let enthalpy_error = point.solution().enthalpy_error().abs();
                    assert!(enthalpy_error <= point.solution().enthalpy_error_limit_joules());
                    let validation = point
                        .solution()
                        .equilibrium()
                        .accepted_solution()
                        .validation();
                    assert!(validation.residual_l2_norm.is_finite());
                    assert!(validation.max_abs_element_balance_error.is_finite());
                    assert!(validation.max_abs_element_balance_error <= 1e-6);
                    max_enthalpy_error = max_enthalpy_error.max(enthalpy_error);
                    total_attempts += point.solution().report().inner_backend_attempts();
                    accepted_backends.insert(format!(
                        "{:?}",
                        point
                            .solution()
                            .equilibrium()
                            .solve_report()
                            .accepted_backend
                    ));
                    max_residual = max_residual.max(validation.residual_l2_norm);
                    max_balance = max_balance.max(validation.max_abs_element_balance_error);
                }
                rows.push(LivePhTargetRangeBackendRow {
                    backend: name.to_string(),
                    status: "OK".to_string(),
                    points: solution.points().len().to_string(),
                    builds: solution.report().formulation_builds().to_string(),
                    reuses: solution.report().formulation_reuses().to_string(),
                    total_ms: live_duration_ms(timing.total()),
                    wall_ms: live_duration_ms(started.elapsed()),
                    mean_ms: live_duration_ms(timing.mean()),
                    worst_ms: live_duration_ms(timing.worst()),
                    attempts: total_attempts.to_string(),
                    accepted_backend: accepted_backends.into_iter().collect::<Vec<_>>().join(","),
                    validation: "OK".to_string(),
                    max_enthalpy_error: format!("{max_enthalpy_error:.3e}"),
                    max_residual: format!("{max_residual:.3e}"),
                    max_balance: format!("{max_balance:.3e}"),
                    failure_point: "-".to_string(),
                    failure_kind: "-".to_string(),
                    error: "-".to_string(),
                });
            }
            Err(error) => {
                match &error {
                    PhRangeError::Point(point) => {
                        assert_eq!(
                            point.source().kind(),
                            ReactionExtentErrorKind::AllBackendsFailed
                        );
                    }
                    PhRangeError::InvalidProblem(source) => {
                        panic!(
                            "strict backend matrix must fail at a point, not during request construction: {source}"
                        );
                    }
                }
                let (failure_point, failure_kind) = match &error {
                    PhRangeError::Point(point) => (
                        point.index().to_string(),
                        format!("{:?}", point.source().kind()),
                    ),
                    PhRangeError::InvalidProblem(source) => {
                        ("-".to_string(), format!("{:?}", source.kind()))
                    }
                };
                rows.push(LivePhTargetRangeBackendRow {
                    backend: name.to_string(),
                    status: "FAILED".to_string(),
                    points: "-".to_string(),
                    builds: "-".to_string(),
                    reuses: "-".to_string(),
                    total_ms: "-".to_string(),
                    wall_ms: live_duration_ms(started.elapsed()),
                    mean_ms: "-".to_string(),
                    worst_ms: "-".to_string(),
                    attempts: "-".to_string(),
                    accepted_backend: "-".to_string(),
                    validation: "NOT ACCEPTED".to_string(),
                    max_enthalpy_error: "-".to_string(),
                    max_residual: "-".to_string(),
                    max_balance: "-".to_string(),
                    failure_point,
                    failure_kind,
                    error: format!("{error:?}"),
                });
            }
        }
    }
    assert!(
        successful > 0,
        "at least one P,H range backend must succeed"
    );
    println!("live real fixed-phase monolithic P,H target-range backend matrix");
    println!("{}", Table::new(rows).with(Style::rounded()));
    assert_eq!(before, live_library_file_snapshot());
}

/// One row in the release-oriented nested/Auto P,H target-range matrix.
///
/// The route is part of the evidence: a successful nested row proves fixed
/// inner-template reuse, while the Auto water/ice row proves that a real phase
/// transition remains visible after monolithic fallback and scalar recovery.
#[derive(Debug, Tabled)]
struct LivePhRangeMatrixRow {
    #[tabled(rename = "Scenario")]
    scenario: String,
    #[tabled(rename = "Mode")]
    mode: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Points")]
    points: String,
    #[tabled(rename = "Builds")]
    builds: String,
    #[tabled(rename = "Reuses")]
    reuses: String,
    #[tabled(rename = "Transitions")]
    transitions: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
    #[tabled(rename = "Mean ms")]
    mean_ms: String,
    #[tabled(rename = "Worst ms")]
    worst_ms: String,
    #[tabled(rename = "Detail")]
    detail: String,
}

/// One row in the strict single-backend monolithic P,H range matrix.
///
/// A failed backend still gets a row with its wall time and typed error. This
/// keeps backend reliability visible instead of making a successful majority
/// hide a method-specific failure.
#[derive(Debug, Tabled)]
struct LivePhTargetRangeBackendRow {
    #[tabled(rename = "Backend")]
    backend: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Points")]
    points: String,
    #[tabled(rename = "Builds")]
    builds: String,
    #[tabled(rename = "Reuses")]
    reuses: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
    #[tabled(rename = "Wall ms")]
    wall_ms: String,
    #[tabled(rename = "Mean ms")]
    mean_ms: String,
    #[tabled(rename = "Worst ms")]
    worst_ms: String,
    #[tabled(rename = "Attempts")]
    attempts: String,
    #[tabled(rename = "Accepted backend")]
    accepted_backend: String,
    #[tabled(rename = "Validation")]
    validation: String,
    #[tabled(rename = "Max dH error")]
    max_enthalpy_error: String,
    #[tabled(rename = "Max residual")]
    max_residual: String,
    #[tabled(rename = "Max balance")]
    max_balance: String,
    #[tabled(rename = "Failure point")]
    failure_point: String,
    #[tabled(rename = "Failure kind")]
    failure_kind: String,
    #[tabled(rename = "Error")]
    error: String,
}

/// Release characterization row for the nested P,H route's inner P,T
/// cascade. Outer scalar evaluations and inner nonlinear work are reported
/// separately so a fast wall time cannot hide excessive inner work.
#[derive(Debug, Tabled)]
struct LiveNestedPhCascadeRow {
    #[tabled(rename = "Fixture")]
    fixture: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Species")]
    species: String,
    #[tabled(rename = "Points")]
    points: String,
    #[tabled(rename = "Outer evals")]
    outer_evaluations: String,
    #[tabled(rename = "Inner attempts")]
    inner_attempts: String,
    #[tabled(rename = "Inner solve ms")]
    inner_solve_ms: String,
    #[tabled(rename = "Builds")]
    builds: String,
    #[tabled(rename = "Reuses")]
    reuses: String,
    #[tabled(rename = "Total ms")]
    total_ms: String,
    #[tabled(rename = "Mean ms")]
    mean_ms: String,
    #[tabled(rename = "Worst ms")]
    worst_ms: String,
    #[tabled(rename = "Max balance")]
    max_balance: String,
    #[tabled(rename = "Error")]
    error: String,
}

fn ph_range_matrix_success_row(
    scenario: &str,
    mode: &str,
    solution: &PhRangeSolution,
) -> LivePhRangeMatrixRow {
    let timing = solution.report().point_timing();
    LivePhRangeMatrixRow {
        scenario: scenario.to_string(),
        mode: mode.to_string(),
        status: "OK".to_string(),
        points: solution.points().len().to_string(),
        builds: solution.report().formulation_builds().to_string(),
        reuses: solution.report().formulation_reuses().to_string(),
        transitions: solution.report().phase_control_transitions().to_string(),
        total_ms: format!("{:.3}", timing.total().as_secs_f64() * 1_000.0),
        mean_ms: format!("{:.3}", timing.mean().as_secs_f64() * 1_000.0),
        worst_ms: format!("{:.3}", timing.worst().as_secs_f64() * 1_000.0),
        detail: "-".to_string(),
    }
}

/// One accepted point in a release-oriented P,H target-range story.
#[derive(Debug, Tabled)]
struct LivePhTargetRangePointRow {
    #[tabled(rename = "Index")]
    index: String,
    #[tabled(rename = "Target J")]
    target_joules: String,
    #[tabled(rename = "Seed K")]
    seed_temperature: String,
    #[tabled(rename = "Solved K")]
    solved_temperature: String,
    #[tabled(rename = "Path")]
    path: String,
    #[tabled(rename = "Fallback")]
    fallback: String,
    #[tabled(rename = "Transitions")]
    transitions: String,
    #[tabled(rename = "Balance")]
    max_balance: String,
    #[tabled(rename = "Elapsed ms")]
    elapsed_ms: String,
}

fn ph_target_range_story_rows(solution: &PhRangeSolution) -> Vec<LivePhTargetRangePointRow> {
    solution
        .points()
        .iter()
        .map(|point| {
            let report = point.report();
            let validation = point
                .solution()
                .equilibrium()
                .accepted_solution()
                .validation();
            LivePhTargetRangePointRow {
                index: report.index().to_string(),
                target_joules: format!("{:.6e}", report.target_enthalpy_joules()),
                seed_temperature: format!("{:.3}", report.seed_temperature()),
                solved_temperature: format!("{:.3}", report.solved_temperature()),
                path: format!("{:?}", report.solve_path()),
                fallback: report
                    .fallback_reason()
                    .map(|reason| format!("{:?}", reason.error_kind()))
                    .unwrap_or_else(|| "-".to_string()),
                transitions: report.phase_control_transitions().to_string(),
                max_balance: format!("{:.3e}", validation.max_abs_element_balance_error),
                elapsed_ms: format!("{:.3}", report.elapsed().as_secs_f64() * 1_000.0),
            }
        })
        .collect()
}

#[test]
#[ignore = "release nested/Auto P,H target-range lifecycle matrix"]
fn live_reactive_and_water_ph_target_range_route_matrix() {
    let before = live_library_file_snapshot();
    let timing_options =
        || EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled);
    let legacy_nr_options = || {
        timing_options()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
            .expect("legacy-NR matrix policy must validate")
    };
    let pressure = 101_325.0;
    let mut rows = Vec::new();

    // Stable real-gas nested range: every P,H target owns its scalar bracket,
    // while the fixed P,T formulation is shared across the whole batch.
    let reactive_resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas must resolve for the route matrix");
    let reactive_layout =
        MultiphaseEquilibriumLayout::new(reactive_resolved.phase_specs().to_vec())
            .expect("reactive gas layout must validate");
    let reactive_initial =
        MultiphaseInitialComposition::from_dense(&reactive_layout, vec![0.1, 0.05, 1.9])
            .expect("reactive gas composition must validate");
    let reactive_thermochemistry =
        ResolvedThermochemistry::from_resolved_system(&reactive_resolved)
            .expect("reactive gas thermochemistry must resolve");
    let reactive_reference_temperatures = [2_300.0, 2_500.0, 2_700.0];
    let mut reactive_targets = reactive_reference_temperatures
        .into_iter()
        .map(|temperature| {
            let conditions = EquilibriumConditions::new(temperature, pressure, pressure)
                .expect("reactive reference conditions must validate");
            let solution = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    &reactive_resolved,
                    conditions,
                    reactive_initial.clone(),
                )
                .with_solve_options(legacy_nr_options()),
            )
            .expect("reactive reference P,T state must solve");
            let target = reactive_thermochemistry
                .enthalpy_model()
                .evaluate_total(solution.component_moles(), temperature)
                .expect("reactive reference enthalpy must be finite");
            (target, temperature)
        })
        .collect::<Vec<_>>();
    reactive_targets.sort_by(|left, right| left.0.total_cmp(&right.0));
    let reactive_range = PhRangeRequest::from_resolved_thermochemistry(
        &reactive_resolved,
        reactive_initial.clone(),
        pressure,
        pressure,
        PhEnthalpyGrid::new(reactive_targets.iter().map(|pair| pair.0).collect())
            .expect("reactive target grid must be monotone"),
        TemperatureBounds::new(2_100.0, 2_900.0).expect("reactive P,H bounds must validate"),
        reactive_targets[0].1,
        reactive_thermochemistry.clone(),
    )
    .expect("nested reactive range must validate")
    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    .with_solve_options(legacy_nr_options())
    .solve()
    .expect("nested reactive range must solve");
    assert_eq!(reactive_range.report().formulation_builds(), 1);
    assert!(reactive_range.report().formulation_reuses() > 0);
    assert!(reactive_range
        .points()
        .iter()
        .all(|point| point.solution().enthalpy_error().abs()
            <= point.solution().enthalpy_error_limit_joules()));
    assert!(reactive_range.points().iter().all(|point| {
        let equilibrium = point.solution().equilibrium();
        let validation = equilibrium.accepted_solution().validation();
        let mole_scale = equilibrium
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale
            && equilibrium
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0)
    }));
    let mut nested_row = ph_range_matrix_success_row("reactive-gas", "Nested", &reactive_range);
    nested_row.detail = "fixed-template=shared".to_string();
    rows.push(nested_row);

    // Real water/ice Auto range: the monolithic phase-control attempt may be
    // rejected, but Auto must preserve that evidence and publish nested
    // accepted points with the phase transition intact.
    let water_resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_gas_solid_water_spec(),
        live_repository(),
    )
    .expect("water/ice system must resolve for the route matrix");
    let water_layout = MultiphaseEquilibriumLayout::new(water_resolved.phase_specs().to_vec())
        .expect("water/ice layout must validate");
    let water_initial =
        MultiphaseInitialComposition::from_dense(&water_layout, vec![0.5, 0.25, 0.0])
            .expect("water/ice composition must validate");
    let water_thermochemistry = ResolvedThermochemistry::from_resolved_system(&water_resolved)
        .expect("water/ice thermochemistry must resolve");
    let water_reference_temperatures = [250.0, 260.0, 270.0];
    let mut water_targets = water_reference_temperatures
        .into_iter()
        .map(|temperature| {
            let solution = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    &water_resolved,
                    EquilibriumConditions::new(temperature, pressure, pressure)
                        .expect("water reference conditions must validate"),
                    water_initial.clone(),
                )
                .with_solve_options(timing_options())
                .with_phase_control_policy(PhaseControlPolicy::default()),
            )
            .expect("water reference P,T state must solve");
            let target = water_thermochemistry
                .enthalpy_model()
                .evaluate_total(solution.component_moles(), temperature)
                .expect("water reference enthalpy must be finite");
            (target, temperature)
        })
        .collect::<Vec<_>>();
    water_targets.sort_by(|left, right| left.0.total_cmp(&right.0));
    let water_range = PhRangeRequest::from_resolved_thermochemistry(
        &water_resolved,
        water_initial.clone(),
        pressure,
        pressure,
        PhEnthalpyGrid::new(water_targets.iter().map(|pair| pair.0).collect())
            .expect("water target grid must be monotone"),
        TemperatureBounds::new(250.0, 270.0).expect("water P,H bounds must validate"),
        water_targets[0].1,
        water_thermochemistry.clone(),
    )
    .expect("Auto water/ice range must validate")
    .with_ph_solve_mode(PhSolveMode::Auto)
    .with_phase_control_policy(PhaseControlPolicy::default())
    .with_solve_options(timing_options())
    .solve()
    .expect("Auto water/ice range must solve");
    assert!(water_range.report().phase_control_transitions() > 0);
    assert!(water_range.points().iter().any(|point| {
        point.report().solve_path() == PhSolvePath::NestedTemperature
            && point.report().fallback_reason().is_some()
    }));
    for (index, point) in water_range.points().iter().enumerate() {
        let equilibrium = point.solution().equilibrium();
        let validation = equilibrium.accepted_solution().validation();
        let mole_scale = equilibrium
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(
            equilibrium
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0),
            "Auto water/ice point {index} produced invalid moles"
        );
        assert!(
            validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale,
            "Auto water/ice point {index} violated conservation: {}",
            validation.max_abs_element_balance_error
        );
        assert!(
            point.solution().enthalpy_error().abs()
                <= point.solution().enthalpy_error_limit_joules(),
            "Auto water/ice point {index} violated the enthalpy acceptance contract"
        );
    }
    let fallback_points = water_range
        .points()
        .iter()
        .filter(|point| point.report().fallback_reason().is_some())
        .count();
    let mut auto_row = ph_range_matrix_success_row("water-ice", "Auto", &water_range);
    auto_row.detail = format!(
        "nested_points={} fallback_points={fallback_points}",
        water_range
            .points()
            .iter()
            .filter(|point| point.report().solve_path() == PhSolvePath::NestedTemperature)
            .count(),
    );
    rows.push(auto_row);

    // An unreachable later target must abort the complete range. The first
    // accepted point is intentionally not returned as a partial solution, and
    // the typed error must identify the failing point index.
    let rollback_targets = PhEnthalpyGrid::new(vec![reactive_targets[0].0, 1.0e12])
        .expect("rollback target grid must be monotone");
    let rollback_request = PhRangeRequest::from_resolved_thermochemistry(
        &reactive_resolved,
        reactive_initial,
        pressure,
        pressure,
        rollback_targets,
        TemperatureBounds::new(2_100.0, 2_900.0).expect("rollback bounds must validate"),
        reactive_targets[0].1,
        reactive_thermochemistry,
    )
    .expect("rollback range must validate")
    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    .with_solve_options(legacy_nr_options());
    let rollback_started = std::time::Instant::now();
    match rollback_request.solve() {
        Err(PhRangeError::Point(error)) => {
            assert_eq!(error.index(), 1);
            rows.push(LivePhRangeMatrixRow {
                scenario: "unreachable-target".to_string(),
                mode: "Nested".to_string(),
                status: "ROLLBACK".to_string(),
                points: "0".to_string(),
                builds: "-".to_string(),
                reuses: "-".to_string(),
                transitions: "-".to_string(),
                total_ms: format!("{:.3}", rollback_started.elapsed().as_secs_f64() * 1_000.0),
                mean_ms: "-".to_string(),
                worst_ms: "-".to_string(),
                detail: format!("failed_point={}", error.index()),
            });
        }
        Err(other) => panic!("rollback returned the wrong typed error: {other}"),
        Ok(_) => panic!("unreachable target incorrectly published a complete range"),
    }

    assert_eq!(before, live_library_file_snapshot());
    println!(
        "live P,H nested/Auto route matrix (real NASA gas and water/ice)\n{}",
        Table::new(rows).with(Style::rounded())
    );
}

/// Release-only continuation story for a materially larger real gas system.
///
/// Nine enthalpy targets exercise accepted-state hand-off repeatedly while
/// twenty local NASA records make every inner `P,T` solve a nontrivial
/// reaction-basis problem. The fixture deliberately uses the nested route:
/// it characterizes reusable fixed-P,T preparation separately from coupled
/// monolithic P,H backend behavior.
#[test]
#[ignore = "release large real-data nested P,H target-range story"]
fn live_large_element_limited_nested_ph_target_range_story() {
    let before = live_library_file_snapshot();
    let (selected, resolved) = live_large_element_limited_gas_resolved();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("large real gas layout must validate");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; selected.len()])
        .expect("large real gas composition must validate");
    let pressure = 101_325.0;
    let options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
        .expect("large nested P,H policy must validate")
        .with_timing_mode(EquilibriumTimingMode::Enabled);
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("large local NASA gas system must provide P,H capabilities");
    // Keep inverse targets strictly inside the shared NASA interval. Exact
    // endpoints are covered by the dedicated interval/Jacobian story; using
    // them here would turn a continuation benchmark into a comparison of two
    // independent endpoint acceptance tolerances.
    let temperatures = [
        1_005.0, 1_016.25, 1_027.5, 1_038.75, 1_050.0, 1_061.25, 1_072.5, 1_083.75, 1_095.0,
    ];
    let mut targets = temperatures
        .into_iter()
        .map(|temperature| {
            let solution = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    &resolved,
                    EquilibriumConditions::new(temperature, pressure, pressure)
                        .expect("large reference conditions must validate"),
                    initial.clone(),
                )
                .with_solve_options(options.clone()),
            )
            .expect("large reference P,T point must solve");
            let target = thermochemistry
                .enthalpy_model()
                .evaluate_total(solution.component_moles(), temperature)
                .expect("large reference target enthalpy must be finite");
            (target, temperature)
        })
        .collect::<Vec<_>>();
    targets.sort_by(|left, right| left.0.total_cmp(&right.0));

    let range = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial,
        pressure,
        pressure,
        PhEnthalpyGrid::new(targets.iter().map(|pair| pair.0).collect())
            .expect("large P,H targets must be strictly monotone"),
        TemperatureBounds::new(1_000.0, 1_100.0).expect("large P,H bounds must validate"),
        targets[0].1,
        thermochemistry,
    )
    .expect("large nested P,H request must validate")
    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    .with_solve_options(options)
    .solve()
    .expect("large nested P,H target range must solve");

    assert_eq!(range.points().len(), targets.len());
    assert_eq!(range.report().continuation_points(), targets.len() - 1);
    assert_eq!(range.report().formulation_builds(), 1);
    assert!(
        range.report().formulation_reuses() >= targets.len() - 1,
        "the large nested range must reuse its fixed P,T template"
    );
    for (index, (point, (_, expected_temperature))) in
        range.points().iter().zip(targets.iter()).enumerate()
    {
        let equilibrium = point.solution().equilibrium();
        let validation = equilibrium.accepted_solution().validation();
        let mole_scale = equilibrium
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(
            (point.solution().temperature() - expected_temperature).abs() < 1.0e-3,
            "large P,H point {index} recovered {} K instead of {expected_temperature} K",
            point.solution().temperature()
        );
        assert!(
            equilibrium
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0),
            "large P,H point {index} produced invalid moles"
        );
        assert!(
            validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale,
            "large P,H point {index} violated conservation: {}",
            validation.max_abs_element_balance_error
        );
        assert!(
            point.solution().enthalpy_error().abs()
                <= point.solution().enthalpy_error_limit_joules(),
            "large P,H point {index} violated the enthalpy contract"
        );
    }
    println!(
        "live large nested P,H target range: species={} points={} builds={} reuses={} total={:?}\n{}",
        selected.len(),
        range.points().len(),
        range.report().formulation_builds(),
        range.report().formulation_reuses(),
        range.report().total(),
        Table::new(ph_target_range_story_rows(&range)).with(Style::rounded()),
    );
    assert_eq!(before, live_library_file_snapshot());
}

/// Release-only characterization of the nested P,H route on three real-data
/// scales. Each row reuses one prepared inner P,T formulation and uses the
/// ordered legacy fallback cascade; this measures lifecycle behavior and
/// backend evidence, not a production-backend selection decision. RST-first
/// behavior is characterized separately because its acceptance and
/// monotonicity evidence has a different contract.
#[test]
#[ignore = "release nested P,H inner P,T cascade characterization matrix"]
fn live_nested_ph_inner_pt_cascade_release_matrix() {
    let before = live_library_file_snapshot();
    let pressure = 101_325.0;
    let temperatures = [1_005.0, 1_050.0, 1_095.0];
    let options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::legacy_default(Solvers::NR))
        .expect("nested inner P,T cascade policy must validate")
        .with_timing_mode(EquilibriumTimingMode::Enabled);
    let mut rows = Vec::new();

    for (fixture, species_count) in [("small", 5usize), ("medium", 20), ("large", 100)] {
        let (selected, resolved) = live_element_limited_gas_resolved_up_to(species_count);
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
            .expect("nested cascade layout must validate");
        let initial = MultiphaseInitialComposition::from_dense(&layout, vec![1e-3; selected.len()])
            .expect("nested cascade composition must validate");
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("nested cascade thermochemistry must resolve");
        let mut targets = temperatures
            .into_iter()
            .map(|temperature| {
                let solution = solve_resolved_pt(
                    ResolvedPhaseEquilibriumRequest::new(
                        &resolved,
                        EquilibriumConditions::new(temperature, pressure, pressure)
                            .expect("nested cascade reference conditions must validate"),
                        initial.clone(),
                    )
                    .with_solve_options(options.clone()),
                )
                .expect("nested cascade reference P,T point must solve");
                let target = thermochemistry
                    .enthalpy_model()
                    .evaluate_total(solution.component_moles(), temperature)
                    .expect("nested cascade target enthalpy must be finite");
                (target, temperature)
            })
            .collect::<Vec<_>>();
        targets.sort_by(|left, right| left.0.total_cmp(&right.0));

        let range = PhRangeRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            pressure,
            pressure,
            PhEnthalpyGrid::new(targets.iter().map(|pair| pair.0).collect())
                .expect("nested cascade target grid must be monotone"),
            TemperatureBounds::new(1_000.0, 1_100.0)
                .expect("nested cascade temperature bounds must validate"),
            targets[0].1,
            thermochemistry,
        )
        .expect("nested cascade range must validate")
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_solve_options(options.clone())
        .solve()
        .expect("nested cascade range must solve");

        assert_eq!(range.points().len(), temperatures.len());
        assert_eq!(range.report().formulation_builds(), 1);
        assert!(range.report().formulation_reuses() >= temperatures.len() - 1);

        let mut outer_evaluations = 0usize;
        let mut inner_attempts = 0usize;
        let mut inner_solve_ms = 0.0;
        let mut max_balance: f64 = 0.0;
        for point in range.points() {
            outer_evaluations += point.solution().report().iterations();
            inner_attempts += point.solution().report().inner_backend_attempts();
            inner_solve_ms += point
                .solution()
                .report()
                .inner_timing()
                .nonlinear_solve()
                .as_secs_f64()
                * 1_000.0;
            let equilibrium = point.solution().equilibrium();
            let validation = equilibrium.accepted_solution().validation();
            let mole_scale = equilibrium
                .component_moles()
                .iter()
                .map(|moles| moles.abs())
                .sum::<f64>();
            max_balance = max_balance.max(validation.max_abs_element_balance_error);
            assert!(
                validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale,
                "{fixture} nested point violated conservation: {}",
                validation.max_abs_element_balance_error
            );
            assert!(
                point.solution().enthalpy_error().abs()
                    <= point.solution().enthalpy_error_limit_joules(),
                "{fixture} nested point violated enthalpy acceptance"
            );
        }
        assert!(inner_attempts >= range.points().len());

        let timing = range.report().point_timing();
        rows.push(LiveNestedPhCascadeRow {
            fixture: fixture.to_string(),
            status: "OK".to_string(),
            species: selected.len().to_string(),
            points: range.points().len().to_string(),
            outer_evaluations: outer_evaluations.to_string(),
            inner_attempts: inner_attempts.to_string(),
            inner_solve_ms: format!("{inner_solve_ms:.3}"),
            builds: range.report().formulation_builds().to_string(),
            reuses: range.report().formulation_reuses().to_string(),
            total_ms: format!("{:.3}", timing.total().as_secs_f64() * 1_000.0),
            mean_ms: format!("{:.3}", timing.mean().as_secs_f64() * 1_000.0),
            worst_ms: format!("{:.3}", timing.worst().as_secs_f64() * 1_000.0),
            max_balance: format!("{max_balance:.3e}"),
            error: "-".to_string(),
        });
    }

    println!(
        "live nested P,H inner P,T cascade release matrix\n{}",
        Table::new(rows).with(Style::rounded())
    );
    assert_eq!(before, live_library_file_snapshot());
}

/// Release-only denser real water/ice `P,H` range through the `Auto` route.
///
/// The target sequence crosses the same low-temperature phase boundary more
/// densely than the compact route matrix. It is intentionally a route story,
/// not a claim that every point must select the same numerical path: accepted
/// points may remain monolithic while a difficult point uses the classified
/// nested fallback.
#[test]
#[ignore = "release dense real-data water/ice Auto P,H target-range story"]
fn live_water_ice_auto_ph_target_range_story() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_gas_solid_water_spec(),
        live_repository(),
    )
    .expect("water/ice system must resolve for the dense P,H story");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("water/ice layout must validate");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0])
        .expect("water/ice composition must validate");
    let pressure = 101_325.0;
    let options = EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled);
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("water/ice thermochemistry must provide P,H capabilities");
    let temperatures = [250.0, 255.0, 260.0, 265.0, 270.0];
    let mut targets = temperatures
        .into_iter()
        .map(|temperature| {
            let solution = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    &resolved,
                    EquilibriumConditions::new(temperature, pressure, pressure)
                        .expect("water/ice reference conditions must validate"),
                    initial.clone(),
                )
                .with_solve_options(options.clone())
                .with_phase_control_policy(PhaseControlPolicy::default()),
            )
            .expect("water/ice reference P,T point must solve");
            let target = thermochemistry
                .enthalpy_model()
                .evaluate_total(solution.component_moles(), temperature)
                .expect("water/ice reference target enthalpy must be finite");
            (target, temperature)
        })
        .collect::<Vec<_>>();
    targets.sort_by(|left, right| left.0.total_cmp(&right.0));

    let range = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial,
        pressure,
        pressure,
        PhEnthalpyGrid::new(targets.iter().map(|pair| pair.0).collect())
            .expect("water/ice P,H targets must be strictly monotone"),
        TemperatureBounds::new(250.0, 270.0).expect("water/ice P,H bounds must validate"),
        targets[0].1,
        thermochemistry,
    )
    .expect("dense Auto water/ice P,H request must validate")
    .with_ph_solve_mode(PhSolveMode::Auto)
    .with_phase_control_policy(PhaseControlPolicy::default())
    .with_solve_options(options)
    .solve()
    .expect("dense Auto water/ice P,H range must solve");

    assert_eq!(range.points().len(), targets.len());
    assert_eq!(range.report().continuation_points(), targets.len() - 1);
    assert!(
        range.report().phase_control_transitions() > 0,
        "the dense water/ice range must retain a real phase transition"
    );
    assert!(range.points().iter().any(|point| {
        point.report().solve_path() == PhSolvePath::NestedTemperature
            && point.report().fallback_reason().is_some()
    }));
    for (index, (point, (_, expected_temperature))) in
        range.points().iter().zip(targets.iter()).enumerate()
    {
        let equilibrium = point.solution().equilibrium();
        let validation = equilibrium.accepted_solution().validation();
        let mole_scale = equilibrium
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(
            (point.solution().temperature() - expected_temperature).abs() < 1.0e-3,
            "water/ice P,H point {index} recovered {} K instead of {expected_temperature} K",
            point.solution().temperature()
        );
        assert!(
            equilibrium
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0),
            "water/ice P,H point {index} produced invalid moles"
        );
        assert!(
            validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * mole_scale,
            "water/ice P,H point {index} violated conservation: {}",
            validation.max_abs_element_balance_error
        );
        assert!(
            point.solution().enthalpy_error().abs()
                <= point.solution().enthalpy_error_limit_joules(),
            "water/ice P,H point {index} violated the enthalpy contract"
        );
    }
    let fallback_points = range
        .points()
        .iter()
        .filter(|point| point.report().fallback_reason().is_some())
        .count();
    println!(
        "live dense water/ice Auto P,H target range: points={} transitions={} fallback_points={} total={:?}\n{}",
        range.points().len(),
        range.report().phase_control_transitions(),
        fallback_points,
        range.report().total(),
        Table::new(ph_target_range_story_rows(&range)).with(Style::rounded()),
    );
    assert_eq!(before, live_library_file_snapshot());
}

/// Release-only validation matrix for the resolved real-data P,H boundary.
///
/// The matrix deliberately mixes request-construction failures with one
/// solver-level unreachable target. A production facade must reject all of
/// them transactionally and must never turn a failed point into a partial
/// P,H range solution.
#[test]
#[ignore = "release real-data P,H invalid-input and rollback matrix"]
fn live_real_ph_invalid_input_and_rollback_matrix() {
    let before = live_library_file_snapshot();
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        live_reactive_gas_spec(),
        live_repository(),
    )
    .expect("reactive gas system must resolve for the P,H validation matrix");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("reactive gas layout must validate");
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])
        .expect("reactive gas composition must validate");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("reactive gas thermochemistry must provide P,H capabilities");
    let pressure = 101_325.0;

    let mut rows = Vec::new();
    let non_finite = PhEnthalpyGrid::new(vec![f64::NAN]);
    assert!(matches!(
        non_finite,
        Err(ReactionExtentError::InvalidProblem { .. })
    ));
    rows.push(LivePhValidationRow {
        scenario: "non-finite target".to_string(),
        status: "REJECTED".to_string(),
        detail: "typed grid validation".to_string(),
    });

    let reversed_bounds = TemperatureBounds::new(2_900.0, 2_100.0);
    assert!(reversed_bounds.is_err());
    rows.push(LivePhValidationRow {
        scenario: "reversed temperature bounds".to_string(),
        status: "REJECTED".to_string(),
        detail: "typed bounds validation".to_string(),
    });

    let invalid_pressure = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial.clone(),
        0.0,
        pressure,
        PhEnthalpyGrid::new(vec![-1.0e5]).expect("finite pressure test target must validate"),
        TemperatureBounds::new(2_100.0, 2_900.0).expect("real gas bounds must validate"),
        2_500.0,
        thermochemistry.clone(),
    );
    assert!(matches!(
        invalid_pressure,
        Err(PhRangeError::InvalidProblem(_))
    ));
    rows.push(LivePhValidationRow {
        scenario: "non-positive pressure".to_string(),
        status: "REJECTED".to_string(),
        detail: "typed condition validation".to_string(),
    });

    let unreachable = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial,
        pressure,
        pressure,
        PhEnthalpyGrid::new(vec![1.0e30]).expect("finite unreachable target is a valid grid"),
        TemperatureBounds::new(2_100.0, 2_900.0).expect("real gas bounds must validate"),
        2_500.0,
        thermochemistry,
    )
    .expect("unreachable finite target remains structurally valid")
    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    .with_solve_options(
        EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
            .expect("single legacy backend policy must validate"),
    )
    .solve();

    match unreachable {
        Err(PhRangeError::Point(error)) => {
            assert_eq!(error.index(), 0);
            rows.push(LivePhValidationRow {
                scenario: "finite unreachable target".to_string(),
                status: "ROLLED BACK".to_string(),
                detail: "point error; no range published".to_string(),
            });
        }
        Err(error) => panic!("unreachable target returned the wrong error: {error}"),
        Ok(_) => panic!("unreachable target incorrectly published a complete range"),
    }

    println!("live real P,H invalid-input and rollback matrix");
    println!("{}", Table::new(rows).with(Style::rounded()));
    assert_eq!(before, live_library_file_snapshot());
}

/// Compact row for the real-data P,H validation matrix.
#[derive(Debug, Tabled)]
struct LivePhValidationRow {
    #[tabled(rename = "Scenario")]
    scenario: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Detail")]
    detail: String,
}
