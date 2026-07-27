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

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    TemperatureInterpolationPolicy, TemperatureInterpolationSpace, TemperaturePostprocessingPolicy,
    TemperatureResamplingGrid, postprocess_temperature_range_solution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
    TemperatureGrid, TemperatureRangeDirection, TemperatureRangePointPreparation,
    TemperatureRangeRequest,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode;
#[allow(deprecated)]
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver_for_T_range;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{InitialPhaseSet, PhaseStatus};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    EquilibriumSolveOptions, PhaseControlPolicy, PhaseEquilibriumPipelineError,
    PhaseEquilibriumPipelineRequest, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::User_PhaseOrSolution::{
    CustomSubstance, PhaseOrSolution, ResolvedPhaseSystem, SubstanceSystemFactory,
    SubstanceSystemSpec, SubstanceSystemSpecBuilder, SubstancesContainer,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::library_manager::with_library_manager;

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

fn live_resolved_system() -> PhaseOrSolution {
    let repository = live_repository();
    let spec = live_multiphase_spec();
    match spec
        .resolve_with_repository(Arc::clone(&repository))
        .expect("live multi-phase spec must resolve against the bundled repository")
    {
        CustomSubstance::PhaseOrSolution(system) => system,
        other => panic!("expected a multi-phase resolved system, got {other:?}"),
    }
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
    assert_eq!(
        resolved
            .resolution_report()
            .expect("live repository must retain lookup provenance")
            .nist_fallback_enabled(),
        false
    );
    assert_eq!(
        resolved
            .resolution_report()
            .expect("live repository must retain lookup provenance")
            .phases()
            .len(),
        2
    );
    assert_eq!(
        resolved
            .phase_data_view()
            .get(&Some("gas".to_string()))
            .unwrap()
            .substances(),
        &["H2O".to_string(), "O2".to_string()]
    );
    assert_eq!(
        resolved
            .phase_data_view()
            .get(&Some("liquid".to_string()))
            .unwrap()
            .substances(),
        &["H2O".to_string()]
    );
}

#[test]
fn live_repository_distinguishes_gaseous_and_solid_water_records() {
    let resolved = match live_gas_solid_water_spec()
        .resolve_with_repository(live_repository())
        .expect("gas/solid water spec must resolve without NIST")
    {
        CustomSubstance::PhaseOrSolution(system) => system,
        other => panic!("expected a multi-phase gas/solid system, got {other:?}"),
    };

    assert!(
        !resolved
            .resolution_report()
            .expect("resolved facade must retain lookup provenance")
            .nist_fallback_enabled()
    );
    assert_eq!(resolved.phase_specs().len(), 2);
    assert_eq!(
        resolved
            .phase_data_view()
            .get(&Some("gas".to_string()))
            .expect("gas phase must be present")
            .substances(),
        &["H2O".to_string(), "O2".to_string()]
    );
    assert_eq!(
        resolved
            .phase_data_view()
            .get(&Some("solid".to_string()))
            .expect("solid phase must be present")
            .substances(),
        &["H2O(s)".to_string()]
    );

    let report = resolved
        .resolution_report()
        .expect("resolved facade must retain lookup provenance");
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
    let resolved = match live_gas_spec()
        .resolve_with_repository(Arc::clone(&live_repository()))
        .expect("live gas spec must resolve against the bundled repository")
    {
        CustomSubstance::OnePhase(system) => system,
        other => panic!("expected a live single-phase resolved system, got {other:?}"),
    };
    let before_gas = resolved.subs_data_view().clone();

    let outcome = live_gas_pipeline_request()
        .solve()
        .expect("live pipeline must solve");

    assert_eq!(outcome.lookup_report(), outcome.resolved().report());
    assert_eq!(
        outcome.solution().build_report().lookup_report(),
        outcome.lookup_report()
    );
    assert_eq!(outcome.resolved().phase_specs().len(), 1);
    assert!(
        outcome
            .solution()
            .component_moles()
            .iter()
            .all(|value| value.is_finite())
    );
    assert!(
        outcome
            .solution()
            .component_moles()
            .iter()
            .sum::<f64>()
            .is_finite()
    );

    assert_eq!(
        before_gas.substances(),
        resolved.subs_data_view().substances()
    );
    assert_eq!(
        before_gas.library_priorities(),
        resolved.subs_data_view().library_priorities()
    );
    assert_eq!(
        before_gas.explicit_search_map(),
        resolved.subs_data_view().explicit_search_map()
    );
}

#[test]
fn live_pipeline_keeps_canonical_thermochemistry_files_byte_for_byte_unchanged() {
    let before = live_library_file_snapshot();

    let outcome = live_multiphase_pipeline_request()
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .expect("read-only live pipeline must solve");

    assert!(
        outcome
            .solution()
            .component_moles()
            .iter()
            .all(|value| value.is_finite())
    );
    assert_eq!(before, live_library_file_snapshot());
}

#[test]
fn live_reactive_gas_solves_with_conserved_elements_and_keq_evidence() {
    let outcome = live_reactive_gas_pipeline_request()
        .solve()
        .expect("live H2/O2/H2O equilibrium must solve");
    let validation = outcome.solution().accepted_solution().validation();

    assert!(
        outcome
            .solution()
            .component_moles()
            .iter()
            .all(|moles| moles.is_finite() && *moles >= 0.0)
    );
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
    assert!(
        outcome
            .solution()
            .component_moles()
            .iter()
            .sum::<f64>()
            .is_finite()
    );
    assert!(
        outcome
            .solution()
            .phase_total(&crate::Thermodynamics::phase_layout::PhaseId::new(Some(
                "liquid".to_string()
            )))
            .is_some()
    );
    assert!(
        outcome
            .solution()
            .summary_rows()
            .iter()
            .any(|row| row.section == "phase_control")
    );
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

fn live_exact_element_gas_candidates_up_to(limit: usize) -> Vec<String> {
    assert!(limit > 0, "live candidate limit must be positive");
    let mut catalog = ThermoData::try_new().expect("bundled catalog must load");
    let mut candidates: Vec<String> = catalog
        .search_by_exact_elements(vec!["C".into(), "H".into(), "O".into()])
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
        "the live NASA gas catalog must provide at least {limit} exact C/H/O candidates, got {}",
        candidates.len(),
    );
    candidates.into_iter().take(limit).collect()
}

fn live_large_exact_element_gas_candidates() -> Vec<String> {
    live_exact_element_gas_candidates_up_to(20)
}

fn live_large_exact_element_gas_request() -> (Vec<String>, PhaseEquilibriumPipelineRequest) {
    let selected = live_large_exact_element_gas_candidates();
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(selected.clone()))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("large exact-element spec must remain valid");
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

fn live_large_exact_element_gas_resolved() -> (Vec<String>, ResolvedPhaseSystem) {
    let selected = live_large_exact_element_gas_candidates();
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(selected.clone()))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("large exact-element spec must remain valid");
    let resolved =
        SubstanceSystemFactory::resolve_phase_system_with_repository(spec, live_repository())
            .expect("large exact-element system must resolve against the local repository");
    (selected, resolved)
}

#[test]
#[ignore = "explicit release characterization run"]
fn live_large_exact_element_gas_timing_report() {
    // Exercise the real candidate-selection path rather than embedding a
    // synthetic species list. Exact C/H/O matching keeps the chemical universe
    // bounded while still providing a materially larger NASA gas system.
    let (selected, request) = live_large_exact_element_gas_request();
    let outcome = request
        .with_solve_options(
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
        )
        .solve()
        .expect("large exact-element live gas solve must succeed");
    let moles = outcome.solution().component_moles();
    assert_eq!(moles.len(), selected.len());
    assert!(moles.iter().all(|value| value.is_finite() && *value > 0.0));
    let validation = outcome.solution().accepted_solution().validation();
    assert!(validation.residual_l2_norm.is_finite());
    assert!(validation.max_abs_element_balance_error <= 1e-8);
    let timing = outcome.timing_report();

    println!("live large exact C/H/O gas timing report");
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
    println!("  total                          {:?}", timing.total());
}

#[test]
#[ignore = "explicit release solver characterization run"]
fn live_large_exact_element_solver_matrix_compares_rst_and_legacy_backends() {
    // Every backend receives the same real resolved problem and the same
    // initial composition. This is a characterization matrix, not a fallback
    // cascade: failures are reported per backend instead of being hidden by a
    // later successful attempt.
    let (selected, request) = live_large_exact_element_gas_request();
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
        "live large exact C/H/O solver matrix: species={}",
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
#[ignore = "explicit release real-data scaling characterization run"]
fn live_exact_element_release_scaling_matrix() {
    // This intentionally uses the same exact-element search contract as the
    // 20-species backend matrix, but grows only the real local candidate set.
    // It characterizes solver scaling without changing the catalog or hiding a
    // failed point behind a later fallback backend.
    let before = live_library_file_snapshot();
    for count in [20usize, 50, 100] {
        let selected = live_exact_element_gas_candidates_up_to(count);
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
        assert!(
            solution
                .component_moles()
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        );
        assert!(validation.residual_l2_norm.is_finite());
        assert!(
            validation.max_abs_element_balance_error <= balance_limit,
            "{count}-species live balance exceeded limit: error={:e}, limit={balance_limit:e}",
            validation.max_abs_element_balance_error
        );

        println!(
            "live exact C/H/O scaling: backend=legacy-nr species={} total={:?} nonlinear={:?} residual={:e} balance={:e}",
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
    let (selected, resolved) = live_large_exact_element_gas_resolved();
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
                    assert!(
                        point
                            .solution()
                            .component_moles()
                            .iter()
                            .all(|value| value.is_finite() && *value > 0.0)
                    );
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
    .with_phase_control_policy(PhaseControlPolicy::default())
    .solve_temperature_range(
        TemperatureGrid::new(vec![250.0, 260.0, 270.0]).expect("live ice range grid must validate"),
    )
    .expect("live bounded ice range must solve");

    assert!(range.report().phase_control_transitions() > 0);
    assert!(range.report().phase_projection_cache_entries() >= 2);
    assert!(range.report().phase_prepared_cache_entries() >= 2);
    assert!(range.report().phase_rst_cache_entries() >= 2);
    assert!(range.points()[1].report().phase_set_reused());
    assert!(
        range.points()[1..]
            .iter()
            .any(|point| point.report().symbolic_parameter_reused())
    );
    for point in range.points() {
        let solution = point.solution();
        let validation = solution.accepted_solution().validation();
        let scale = solution
            .component_moles()
            .iter()
            .map(|moles| moles.abs())
            .sum::<f64>();
        assert!(
            solution
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite())
        );
        assert!(validation.max_abs_element_balance_error <= 1e-6 + 1e-6 * scale);
    }
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

#[test]
fn live_fixed_point_matches_one_point_typed_temperature_range() {
    let (selected, resolved) = live_large_exact_element_gas_resolved();
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
    .with_solve_options(options)
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
fn live_large_exact_element_legacy_temperature_range_story() {
    // Keep the compatibility baseline while the typed facade is characterized
    // by the release-oriented story below.
    let (selected, _) = live_large_exact_element_gas_request();
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
fn live_large_exact_element_typed_temperature_range_story() {
    let (selected, resolved) = live_large_exact_element_gas_resolved();
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
        assert!(
            result.points()[1..]
                .iter()
                .all(|point| point.report().used_continuation_seed())
        );
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
    assert!(
        outcome
            .solution()
            .phase_control_report()
            .expect("bounded solve must retain its phase-control report")
            .transitions
            .iter()
            .any(|transition| !transition.activated.is_empty())
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
    assert!(
        low.solution()
            .phase_control_report()
            .expect("bounded solve must retain its phase-control report")
            .transitions
            .iter()
            .any(|transition| !transition.activated.is_empty())
    );
    assert!(
        low.solution()
            .build_report()
            .components()
            .iter()
            .any(|component| {
                component.component().label() == "solid::C(gr)"
                    && component.thermo_source().library() == "NASA_cond"
                    && component.thermo_source().record_key() == "C(gr)"
            })
    );
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
    assert!(
        low.solution()
            .summary_rows()
            .iter()
            .any(|row| row.section == "phase_control")
    );
    assert!(
        high.solution()
            .summary_rows()
            .iter()
            .any(|row| row.section == "phase_control")
    );
}
