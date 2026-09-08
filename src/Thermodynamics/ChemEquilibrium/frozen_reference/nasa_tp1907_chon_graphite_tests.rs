//! Offline I5 evidence for NASA TP-1907 Table 11.3E CHON + graphite.

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::{Arc, Mutex};

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
        EquilibriumDiagnosticEvent, EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
        EquilibriumExecutionControl, EquilibriumProgressStage,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_extensive_normalization::ExtensiveNormalization;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_prepared_runner::{
        PreparedEquilibriumRunner, PreparedSolveOutcome,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        LogMolesInitialGuess, PreparedEquilibriumProblem, TraceSpeciesSeedPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
        TemperatureGrid, TemperatureRangePointPreparation, TemperatureRangeRequest,
        TemperatureRangeSolution,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        InitialPhaseSet, PhaseSet, PhaseStatus, compute_phase_stability_reports,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::assertions::{
        AcceptedSolutionContract, assert_bounded_solution_accepted,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenMultiphaseNormalization, FrozenReferenceDataset, FrozenReferenceEvidenceKind,
        NasaTp1907MultiphaseReference,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_interior_seed::{
        ElementFeasibleInteriorSeed, InteriorSeedAvailability, InteriorSeedSettings,
        InteriorSeedTarget, build_element_feasible_interior_seed,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite::{
        ResolvedNasaTp1907ChonGraphiteFixture, TP1907_GAS_PHASE, TP1907_GRAPHITE_PHASE,
        TP1907_LOCAL_GAS_SPECIES, graphite_component, load_nasa_tp1907_chon_graphite_dataset,
        reconstructed_source_feed, validate_published_feed_diagnostics,
        verify_element_equivalent_feed,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
        PurePhaseProductionEvidenceRequest, activation_evidence_from_solution,
        stable_inactive_evidence_from_solution,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
        PhaseEquilibriumBuildRequest, build_phase_equilibrium_problem,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        EquilibriumSolveOptions, ExtensiveNormalizationPolicy, ResolvedPhaseEquilibriumRequest,
        solve_resolved_pt,
    };
    use crate::Thermodynamics::ChemEquilibrium::prelude::PhaseControlPolicy;
    use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
        PreparedPhaseControlOutcome, PreparedPhaseControlRunner,
    };
    use crate::Thermodynamics::phase_layout::PhaseId;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;

    // Regression guards selected after repeated debug/release characterization.
    // They are not NASA TP-1907 source uncertainty or physical acceptance tolerances.
    const MAX_MAJOR_GAS_RELATIVE_ERROR_GUARD: f64 = 0.10;
    const MAX_MAJOR_GAS_RMS_RELATIVE_ERROR_GUARD: f64 = 0.05;
    const MAX_GRAPHITE_680_RELATIVE_ERROR_GUARD: f64 = 0.30;
    const MAX_GRAPHITE_700_ABSOLUTE_SYSTEM_FRACTION_ERROR_GUARD: f64 = 1.0e-3;

    fn dataset() -> FrozenReferenceDataset<NasaTp1907MultiphaseReference> {
        load_nasa_tp1907_chon_graphite_dataset().expect("frozen NASA TP-1907 data must load")
    }

    fn repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository().expect("bundled local repository must be available")
    }

    fn fixture() -> ResolvedNasaTp1907ChonGraphiteFixture {
        ResolvedNasaTp1907ChonGraphiteFixture::resolve_offline(repository())
            .expect("reviewed NASA CHON + graphite records must resolve offline")
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
                fs::read(&path).expect("local library must be readable"),
            )
        })
        .collect()
    }

    fn frozen_snapshot() -> (Vec<u8>, Vec<u8>) {
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nasa_tp1907");
        (
            fs::read(directory.join("chon_graphite_er125_1atm.metadata.json"))
                .expect("frozen metadata must be readable"),
            fs::read(directory.join("chon_graphite_er125_1atm.rows.json"))
                .expect("frozen rows must be readable"),
        )
    }

    fn graphite_tpd(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperature_k: f64,
        evidence_request: &PurePhaseProductionEvidenceRequest,
    ) -> f64 {
        let solution = solve_graphite_pt(fixture, temperature_k);
        graphite_tpd_from_solution(&solution, evidence_request)
    }

    fn solve_graphite_pt(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperature_k: f64,
    ) -> MultiphaseEquilibriumSolution {
        let solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                fixture
                    .conditions_at(temperature_k, 101_325.0)
                    .expect("TPD probe temperature must be inside local NASA coverage"),
                fixture
                    .initial_composition()
                    .expect("TPD probe inventory must validate"),
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .expect("canonical TPD probe must solve");
        solution
    }

    fn solve_graphite_pt_with_trace_floor(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperature_k: f64,
        floor: f64,
    ) -> MultiphaseEquilibriumSolution {
        solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                fixture
                    .conditions_at(temperature_k, 101_325.0)
                    .expect("TP-1907 trace-floor probe temperature must be supported"),
                fixture
                    .initial_composition()
                    .expect("TP-1907 trace-floor inventory must validate"),
            )
            .with_phase_control_policy(PhaseControlPolicy::default())
            .with_solve_options(
                EquilibriumSolveOptions::default()
                    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor }),
            ),
        )
        .expect("canonical bounded TP-1907 trace-floor probe must solve")
    }

    fn prepared_fixed_tp1907_problem(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperature_k: f64,
    ) -> PreparedEquilibriumProblem {
        prepared_fixed_tp1907_problem_at_scale(fixture, temperature_k, 1.0)
    }

    fn prepared_fixed_tp1907_problem_at_scale(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperature_k: f64,
        inventory_scale: f64,
    ) -> PreparedEquilibriumProblem {
        let bundle = build_phase_equilibrium_problem(
            PhaseEquilibriumBuildRequest::new(
                fixture.resolved(),
                fixture
                    .conditions_at(temperature_k, 101_325.0)
                    .expect("TP-1907 fixed-active temperature must be supported"),
                scaled_initial_composition(fixture, inventory_scale),
                TraceSpeciesSeedPolicy::Absolute { floor: 1.0e-30 },
                Default::default(),
            )
            .expect("TP-1907 fixed-active request must validate"),
        )
        .expect("TP-1907 fixed-active problem must build");
        PreparedEquilibriumProblem::new(bundle.into_problem())
            .expect("TP-1907 fixed-active problem must prepare")
    }

    /// Test-only exact extensive normalization of one fixed-`P,T` problem.
    ///
    /// The scale is derived solely from the current physical input inventory;
    /// it is not a fixture constant and does not inspect an accepted solution.
    /// Thermochemistry, conditions, reference pressure, and solver policy are
    /// rebuilt unchanged. Only the physical input `n0` becomes `n0 / s`, with
    /// a normalized active inventory of one mole.
    fn normalized_fixed_problem_from_input(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperature_k: f64,
        physical: &PreparedEquilibriumProblem,
    ) -> (PreparedEquilibriumProblem, f64) {
        let physical_moles = physical.problem().initial_moles();
        let normalization = ExtensiveNormalization::from_physical_moles(physical_moles)
            .expect("TP-1907 physical inventory must define extensive normalization");

        let layout = MultiphaseEquilibriumLayout::new(fixture.resolved().phase_specs().to_vec())
            .expect("TP-1907 normalized layout must remain valid");
        let physical_composition =
            MultiphaseInitialComposition::from_dense(&layout, physical_moles.to_vec())
                .expect("TP-1907 physical composition must remain valid");
        let normalized = normalization
            .normalize_initial_composition(&layout, &physical_composition)
            .expect("normalized TP-1907 composition must remain valid");
        let bundle = build_phase_equilibrium_problem(
            PhaseEquilibriumBuildRequest::new(
                fixture.resolved(),
                fixture
                    .conditions_at(temperature_k, 101_325.0)
                    .expect("TP-1907 normalized temperature must be supported"),
                normalized,
                normalization
                    .normalize_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute {
                        floor: 1.0e-30,
                    })
                    .expect("TP-1907 trace floor must normalize"),
                Default::default(),
            )
            .expect("normalized TP-1907 request must validate"),
        )
        .expect("normalized TP-1907 problem must build");
        (
            PreparedEquilibriumProblem::new(bundle.into_problem())
                .expect("normalized TP-1907 problem must prepare"),
            normalization.physical_inventory_scale(),
        )
    }

    fn scaled_initial_composition(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        factor: f64,
    ) -> MultiphaseInitialComposition {
        let layout = MultiphaseEquilibriumLayout::new(fixture.resolved().phase_specs().to_vec())
            .expect("TP-1907 phase layout must remain valid");
        let moles = fixture
            .initial_composition()
            .expect("TP-1907 inventory must validate")
            .moles()
            .iter()
            .map(|moles| factor * moles)
            .collect();
        MultiphaseInitialComposition::from_dense(&layout, moles)
            .expect("scaled TP-1907 inventory must validate")
    }

    fn active_phase_topology(solution: &MultiphaseEquilibriumSolution) -> String {
        solution
            .phases()
            .iter()
            .filter_map(|phase| {
                matches!(
                    solution.phase_status(phase.id()),
                    Some(PhaseStatus::Active | PhaseStatus::Appeared)
                )
                .then(|| format!("{:?}", phase.id()))
            })
            .collect::<Vec<_>>()
            .join("+")
    }

    fn smallest_active_phase_total(solution: &MultiphaseEquilibriumSolution) -> f64 {
        solution
            .phases()
            .iter()
            .filter(|phase| {
                matches!(
                    solution.phase_status(phase.id()),
                    Some(PhaseStatus::Active | PhaseStatus::Appeared)
                )
            })
            .filter_map(|phase| solution.phase_total(phase.id()))
            .fold(f64::INFINITY, f64::min)
    }

    fn component_moles_by_identity(
        solution: &MultiphaseEquilibriumSolution,
    ) -> BTreeMap<String, f64> {
        solution
            .metadata()
            .components()
            .iter()
            .zip(solution.component_moles())
            .map(|(component, &moles)| (component.id().label(), moles))
            .collect()
    }

    fn gas_mole_fractions_by_identity(
        solution: &MultiphaseEquilibriumSolution,
    ) -> BTreeMap<String, f64> {
        let gas = solution
            .phases()
            .iter()
            .find(|phase| phase.id() == &PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())))
            .expect("TP-1907 gas phase must be represented");
        let range = gas.component_range();
        let moles = &solution.component_moles()[range];
        let total = moles.iter().sum::<f64>();
        assert!(
            total.is_finite() && total > 0.0,
            "active gas total must be positive"
        );
        solution.metadata().components()[gas.component_range()]
            .iter()
            .zip(moles)
            .map(|(component, &moles)| (component.id().label(), moles / total))
            .collect()
    }

    fn max_recovered_relative_error(
        expected: &BTreeMap<String, f64>,
        actual: &BTreeMap<String, f64>,
        factor: f64,
    ) -> f64 {
        assert_eq!(
            expected.keys().collect::<Vec<_>>(),
            actual.keys().collect::<Vec<_>>(),
            "metamorphic comparisons must align components by phase-qualified identity"
        );
        expected
            .iter()
            .map(|(identity, &expected)| {
                let actual = actual[identity];
                let recovered = actual / factor;
                (expected - recovered).abs() / expected.abs().max(recovered.abs()).max(1.0e-30)
            })
            .fold(0.0_f64, f64::max)
    }

    /// Diagnostic-only outcome for a uniform extensive-scaling experiment.
    ///
    /// This is deliberately a test-suite type rather than a production solver
    /// status. A classification requires evidence from both ordinary and
    /// oracle-style seeds; it must not steer a user calculation.
    #[allow(dead_code)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum ScaleConditioningClassification {
        StrictScaleInvariant,
        ThresholdLimited,
        FreshSeedBasinLimited,
        PhysicalScaleRegression,
    }

    #[derive(Debug, Clone, Copy)]
    struct LogSeedStatistics {
        min: f64,
        max: f64,
        mean: f64,
        dynamic_range: f64,
        trace_coordinates: usize,
    }

    fn log_seed_statistics(seed: &LogMolesInitialGuess, trace_floor: f64) -> LogSeedStatistics {
        let values = seed.as_slice();
        let min = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let trace_log = trace_floor.ln();
        let trace_coordinates = values
            .iter()
            .filter(|&&value| (value - trace_log).abs() <= 1.0e-12)
            .count();
        LogSeedStatistics {
            min,
            max,
            mean: values.iter().sum::<f64>() / values.len() as f64,
            dynamic_range: max - min,
            trace_coordinates,
        }
    }

    fn print_seed_statistics(label: &str, statistics: LogSeedStatistics) {
        println!(
            "  {label:27} min={:10.4} max={:10.4} mean={:10.4} range={:10.4} trace={}",
            statistics.min,
            statistics.max,
            statistics.mean,
            statistics.dynamic_range,
            statistics.trace_coordinates,
        );
    }

    /// Builds the proposed Route E from current input alone.
    ///
    /// Positive physical input coordinates are rescaled to the total physical
    /// inventory represented by their ordinary fresh seed. Zero-input entries
    /// retain their explicit numerical trace coordinates. The returned scale
    /// is diagnostic evidence: a value of one proves this proposal is already
    /// equivalent to the canonical `from_moles_with_policy` construction.
    fn input_derived_active_normalized_seed(
        prepared: &PreparedEquilibriumProblem,
    ) -> (LogMolesInitialGuess, f64, f64, f64) {
        let input = prepared.problem().initial_moles();
        let fresh = prepared.problem().initial_log_moles();
        let input_total = input
            .iter()
            .copied()
            .filter(|amount| *amount > 0.0)
            .sum::<f64>();
        let fresh_active_total = input
            .iter()
            .zip(fresh.as_slice())
            .filter_map(|(&input_amount, &log_seed)| (input_amount > 0.0).then_some(log_seed.exp()))
            .sum::<f64>();
        let scale = input_total / fresh_active_total;
        let values = input
            .iter()
            .zip(fresh.as_slice())
            .map(|(&input_amount, &log_seed)| {
                if input_amount > 0.0 {
                    log_seed + scale.ln()
                } else {
                    log_seed
                }
            })
            .collect();
        (
            LogMolesInitialGuess::new(values)
                .expect("input-derived scale-normalized diagnostic seed must be finite"),
            scale,
            input_total,
            fresh_active_total,
        )
    }

    fn print_component_seed_comparison(
        temperature_k: f64,
        factor: f64,
        baseline: &PreparedEquilibriumProblem,
        scaled: &PreparedEquilibriumProblem,
        accepted_unit_solution: &crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
    ) {
        let base_input = baseline.problem().initial_moles();
        let base_seed = baseline.problem().initial_log_moles().as_slice();
        let scaled_input = scaled.problem().initial_moles();
        let scaled_seed = scaled.problem().initial_log_moles().as_slice();
        let expected_shift = factor.ln();
        println!(
            "  seed components at T={temperature_k:.1} K (input -> fresh -> transformed oracle)"
        );
        println!(
            "    component                         role             input n       fresh n    transformed n   fresh shift err   fresh-oracle log"
        );
        for (index, component) in scaled.problem().components().iter().enumerate() {
            let role = if base_input[index] > 0.0 {
                "PhysicalInput"
            } else {
                "NumericalTrace"
            };
            let transformed_log = accepted_unit_solution.log_moles()[index] + expected_shift;
            let shift_error = if base_input[index] > 0.0 {
                scaled_seed[index] - base_seed[index] - expected_shift
            } else {
                f64::NAN
            };
            println!(
                "    {label:32} {role:16} {input:11.3e} {fresh:11.3e} {transformed:14.3e} {shift_error:17.3e} {delta:17.3e}",
                label = component.label(),
                input = scaled_input[index],
                fresh = scaled_seed[index].exp(),
                transformed = transformed_log.exp(),
                delta = scaled_seed[index] - transformed_log,
            );
        }
    }

    fn print_initial_residual_diagnostics(
        label: &str,
        prepared: &PreparedEquilibriumProblem,
        seed: &LogMolesInitialGuess,
    ) {
        let raw = prepared
            .residual(seed.as_slice())
            .expect("diagnostic seed residual must evaluate");
        let scale = prepared
            .residual_scale()
            .expect("diagnostic row scales must evaluate");
        let scaled = prepared
            .scaled_residual(seed.as_slice(), &scale)
            .expect("diagnostic scaled residual must evaluate");
        let reaction_rows = prepared.reaction_basis().reactions.ncols();
        let max_abs = |rows: &[f64]| {
            rows.iter()
                .fold(0.0_f64, |max, &value| max.max(value.abs()))
        };
        let bounds = prepared
            .finite_log_mole_bounds()
            .expect("diagnostic finite log-mole bounds must evaluate");
        let lower = bounds
            .iter()
            .map(|(lower, _)| *lower)
            .fold(f64::INFINITY, f64::min);
        let upper = bounds
            .iter()
            .map(|(_, upper)| *upper)
            .fold(f64::NEG_INFINITY, f64::max);
        println!(
            "  {label:27} reaction raw={:10.3e} element raw={:10.3e} reaction scaled={:10.3e} element scaled={:10.3e}",
            max_abs(&raw[..reaction_rows]),
            max_abs(&raw[reaction_rows..]),
            max_abs(&scaled[..reaction_rows]),
            max_abs(&scaled[reaction_rows..]),
        );
        println!("  finite log-mole bounds      lower={lower:10.4} upper={upper:10.4}");
    }

    fn strict_fixed_active_runner(
        prepared: PreparedEquilibriumProblem,
    ) -> PreparedEquilibriumRunner {
        let mut runner = PreparedEquilibriumRunner::new(prepared, Vec::new())
            .expect("prepared TP-1907 runner must construct");
        let settings = runner.configure();
        settings.solver = Solvers::LM;
        settings.solver_params.tol = 1.0e-10;
        settings.solver_params.max_iter = 500;
        runner
    }

    /// Builds the diagnostic interior seed over exactly the specified physical
    /// phase set. Inactive phases are zero in `physical_moles` and retain only
    /// their ordinary numerical trace coordinate in `log_moles`.
    fn element_feasible_interior_seed(
        prepared: &PreparedEquilibriumProblem,
        active_phase_mask: &[bool],
    ) -> ElementFeasibleInteriorSeed {
        element_feasible_interior_seed_with_settings(
            prepared,
            active_phase_mask,
            InteriorSeedSettings::default(),
        )
    }

    fn element_feasible_interior_seed_with_settings(
        prepared: &PreparedEquilibriumProblem,
        active_phase_mask: &[bool],
        settings: InteriorSeedSettings,
    ) -> ElementFeasibleInteriorSeed {
        let seed = build_element_feasible_interior_seed(
            prepared.problem().element_composition(),
            prepared.element_totals(),
            prepared.species_phase(),
            active_phase_mask,
            prepared.problem().initial_moles(),
            settings,
        )
        .expect("TP-1907 active set must have a structural interior seed");
        assert!(
            seed.max_element_balance_error <= 1.0e-7 * seed.total_inventory_scale.max(1.0),
            "element-feasible diagnostic seed must preserve the physical inventory, error={:e}, scale={:e}",
            seed.max_element_balance_error,
            seed.total_inventory_scale,
        );
        seed
    }

    fn print_interior_seed(label: &str, seed: &ElementFeasibleInteriorSeed) {
        println!(
            "  {label:27} availability={:?} active={} scale={:.6e} requested={:.3e} achieved={:.3e} balance={:.3e}",
            seed.availability,
            seed.active_species_count,
            seed.total_inventory_scale,
            seed.requested_minimum_fraction,
            seed.achieved_minimum_fraction,
            seed.max_element_balance_error,
        );
    }

    /// Solves the physically correct 720 K reduced formulation. The graphite
    /// phase remains a stability candidate, but is not forced into the
    /// positive-log nonlinear system when its accepted TPD is positive.
    ///
    /// This is diagnostic infrastructure only: production enters this state
    /// through its ordinary bounded phase-control lifecycle rather than an
    /// externally supplied accepted-state seed.
    fn solve_reduced_gas_only_from_seed(
        prepared: &PreparedEquilibriumProblem,
        seed: LogMolesInitialGuess,
    ) -> Result<PreparedPhaseControlOutcome, ReactionExtentError> {
        let problem = prepared.problem().clone();
        let conditions = problem.conditions();
        let gibbs = problem.gibbs().to_vec();
        let phase_count = problem.phases().len();
        let gas_only = PhaseSet::from_policy(
            &InitialPhaseSet::Explicit {
                active: vec![PhaseIndex::new(0, phase_count)?],
                excluded: Vec::new(),
            },
            &[true, false],
        )?;
        let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false)?;
        let settings = runner.configure_solver();
        settings.solver = Solvers::LM;
        settings.solver_params.tol = 1.0e-10;
        settings.solver_params.max_iter = 500;
        runner.retarget_numeric(conditions, seed, gibbs)?;
        runner.set_continuation_phase_set(gas_only)?;
        runner.solve()
    }

    fn assert_fixed_active_scaled_solution(
        baseline: &crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        scaled: &crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
        factor: f64,
        temperature_k: f64,
    ) {
        for (index, (&expected, &actual)) in baseline.moles().iter().zip(scaled.moles()).enumerate()
        {
            let recovered = actual / factor;
            let relative =
                (expected - recovered).abs() / expected.abs().max(recovered.abs()).max(1.0e-30);
            assert!(
                relative <= 1.0e-7,
                "TP-1907 {temperature_k} K fixed-active large-scale mismatch at component {index}: baseline={expected:e}, recovered={recovered:e}, relative={relative:e}"
            );
        }
        assert!(
            scaled.validation().max_abs_reaction_affinity <= 1.0e-6,
            "the scaled fixed-active state must retain the physical reaction-affinity contract at {temperature_k} K"
        );
        assert!(
            scaled.validation().max_abs_element_balance_error <= 1.0e-6,
            "the scaled fixed-active state must retain the element-balance contract at {temperature_k} K"
        );
    }

    /// Compares an internally normalized solution after physical extensive
    /// reconstruction with a solution of the original input problem.
    fn reconstructed_mole_and_gas_fraction_errors(
        normalized_moles: &[f64],
        normalization_factor: f64,
        physical_moles: &[f64],
        species_phase: &[usize],
    ) -> (f64, f64) {
        let normalized_gas_total = normalized_moles
            .iter()
            .enumerate()
            .filter_map(|(index, &amount)| (species_phase[index] == 0).then_some(amount))
            .sum::<f64>();
        let physical_gas_total = physical_moles
            .iter()
            .enumerate()
            .filter_map(|(index, &amount)| (species_phase[index] == 0).then_some(amount))
            .sum::<f64>();
        let mut max_mole_error = 0.0_f64;
        let mut max_gas_fraction_error = 0.0_f64;
        for (index, (&normalized, &physical)) in
            normalized_moles.iter().zip(physical_moles).enumerate()
        {
            let reconstructed = normalized * normalization_factor;
            let relative = (reconstructed - physical).abs()
                / reconstructed.abs().max(physical.abs()).max(1.0e-30);
            max_mole_error = max_mole_error.max(relative);
            if species_phase[index] == 0 {
                max_gas_fraction_error = max_gas_fraction_error
                    .max((normalized / normalized_gas_total - physical / physical_gas_total).abs());
            }
        }
        (max_mole_error, max_gas_fraction_error)
    }

    /// Computes the graphite TPD for a prepared gas-only formulation without
    /// building a new physical request. The candidate stays inactive.
    fn graphite_tpd_from_prepared_gas_only(
        prepared: &PreparedEquilibriumProblem,
        log_moles: &[f64],
    ) -> f64 {
        let gas_only = PhaseSet::from_policy(
            &InitialPhaseSet::Explicit {
                active: vec![
                    PhaseIndex::new(0, prepared.problem().phases().len())
                        .expect("gas phase index must be valid"),
                ],
                excluded: Vec::new(),
            },
            &[true, false],
        )
        .expect("gas-only phase set must be valid");
        let stability = compute_phase_stability_reports(
            log_moles,
            prepared.problem().gibbs(),
            prepared.problem().phases(),
            prepared.species_phase(),
            prepared.problem().element_composition(),
            prepared.problem().conditions().temperature(),
            prepared.problem().conditions().pressure(),
            prepared.problem().conditions().reference_pressure(),
            &gas_only,
        )
        .expect("graphite stability evidence must evaluate");
        stability[1]
            .minimum_tpd
            .expect("inactive graphite must have canonical TPD evidence")
    }

    /// Uses a stricter inner nonlinear contract only for the extensivity
    /// metamorphic test. Production defaults deliberately balance accuracy and
    /// throughput; comparing inventories three orders apart needs numerical
    /// error materially below the physical invariant.
    fn strict_extensivity_solve_options() -> EquilibriumSolveOptions {
        EquilibriumSolveOptions::default()
            // The ordinary RST-first production cascade limits every backend
            // to 50 iterations. This invariant needs a deterministic
            // high-accuracy reference solve instead of testing that budget.
            .with_solver_backend(Solvers::LM)
            .with_max_iterations(500)
            .expect("strict extensivity iteration budget must validate")
            .with_tolerance(1.0e-10)
            .expect("strict extensivity tolerance must validate")
    }

    fn solve_prepared_with_legacy_lm(
        prepared: PreparedEquilibriumProblem,
    ) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution {
        let mut runner = PreparedEquilibriumRunner::new(prepared, Vec::new())
            .expect("prepared TP-1907 runner must construct");
        let settings = runner.configure();
        settings.solver = Solvers::LM;
        // Basis invariance is a formulation test, so solve more tightly than
        // the ordinary production tolerance before comparing coordinates.
        settings.solver_params.tol = 1.0e-10;
        settings.solver_params.max_iter = 500;
        runner
            .solve()
            .expect("legacy LM must solve the fixed TP-1907 formulation")
            .solution
    }

    fn graphite_tpd_from_solution(
        solution: &MultiphaseEquilibriumSolution,
        evidence_request: &PurePhaseProductionEvidenceRequest,
    ) -> f64 {
        let evidence = match solution.phase_status(&evidence_request.candidate_phase) {
            Some(PhaseStatus::Inactive) => {
                stable_inactive_evidence_from_solution(&solution, evidence_request)
            }
            Some(PhaseStatus::Active | PhaseStatus::Appeared) => {
                activation_evidence_from_solution(&solution, evidence_request)
            }
            status => panic!("graphite has unsupported status {status:?} during TPD search"),
        }
        .expect("canonical phase-control result must retain graphite evidence");
        evidence
            .boundary_minimum_tpd
            .filter(|value| value.is_finite())
            .expect("canonical graphite evidence must expose a finite TPD")
    }

    /// Reads final-state TPD directly from the acceptance bundle. Unlike the
    /// transition adapter, this remains valid when continuation starts with
    /// graphite already active and therefore records no new activation.
    fn graphite_final_tpd(solution: &MultiphaseEquilibriumSolution) -> f64 {
        let graphite = PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()));
        let phase_index = solution
            .phases()
            .iter()
            .find(|phase| phase.id() == &graphite)
            .expect("graphite descriptor must be represented")
            .index()
            .index();
        solution
            .acceptance_report()
            .expect("bounded solution must retain final acceptance evidence")
            .phase_stability
            .get(phase_index)
            .and_then(|report| report.minimum_tpd)
            .filter(|value| value.is_finite())
            .expect("final graphite stability report must expose finite TPD")
    }

    /// Bisection over the continuous canonical TPD. This is deliberately a
    /// test-local characterization helper: it does not pretend that the four
    /// printed NASA rows provide an externally tabulated transition temperature.
    fn bisect_graphite_tpd_boundary(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        evidence_request: &PurePhaseProductionEvidenceRequest,
    ) -> (f64, f64, usize) {
        let mut lower_temperature_k = 700.0;
        let mut lower_tpd = graphite_tpd(fixture, lower_temperature_k, evidence_request);
        let mut upper_temperature_k = 720.0;
        let upper_tpd = graphite_tpd(fixture, upper_temperature_k, evidence_request);
        assert!(
            lower_tpd < 0.0 && upper_tpd > 0.0,
            "TPD search requires the source topology bracket: lower={lower_tpd:e}, upper={upper_tpd:e}"
        );

        for iteration in 1..=64 {
            let temperature_k = 0.5 * (lower_temperature_k + upper_temperature_k);
            let tpd = graphite_tpd(fixture, temperature_k, evidence_request);
            if tpd.abs() <= 1.0e-7 || (upper_temperature_k - lower_temperature_k) <= 1.0e-6 {
                return (temperature_k, tpd, iteration);
            }
            if lower_tpd.signum() != tpd.signum() {
                upper_temperature_k = temperature_k;
            } else {
                lower_temperature_k = temperature_k;
                lower_tpd = tpd;
            }
        }
        panic!("canonical graphite TPD boundary did not converge within 64 bisection steps");
    }

    fn solve_graphite_temperature_range(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        temperatures: Vec<f64>,
    ) -> TemperatureRangeSolution {
        TemperatureRangeRequest::new(
            fixture.resolved(),
            fixture
                .initial_composition()
                .expect("range inventory must validate"),
            101_325.0,
            101_325.0,
            TemperatureGrid::new(temperatures).expect("range grid must be monotone"),
        )
        .expect("resolved graphite range request must validate")
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .expect("bounded graphite temperature continuation must solve")
    }

    fn system_mole_fraction(
        solution: &MultiphaseEquilibriumSolution,
        component: &crate::Thermodynamics::phase_layout::PhaseComponentId,
    ) -> f64 {
        let total = solution.component_moles().iter().sum::<f64>();
        assert!(
            total.is_finite() && total > 0.0,
            "accepted system total must be positive"
        );
        solution
            .moles_for(component)
            .expect("fixture component must be represented")
            / total
    }

    fn assert_tp1907_external_regression_envelope(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        reference: &NasaTp1907MultiphaseReference,
        solution: &MultiphaseEquilibriumSolution,
    ) {
        let rows = fixture
            .compare_system_composition(reference, solution)
            .expect("accepted state must align to frozen source identities");
        let mut major_gas_relative_errors = Vec::new();
        for row in &rows {
            match row.identity.as_str() {
                "C(gr)" if reference.temperature_k == 680.0 => {
                    let relative = row
                        .relative_error
                        .expect("active graphite must compare")
                        .abs();
                    assert!(
                        relative <= MAX_GRAPHITE_680_RELATIVE_ERROR_GUARD,
                        "NASA TP-1907 external quality regression at 680 K for C(gr): relative error={relative:e} exceeds reviewed guard={MAX_GRAPHITE_680_RELATIVE_ERROR_GUARD:e}; source={:e}, local={:e}",
                        row.source_system_fraction,
                        row.kithe_system_fraction
                            .expect("active graphite must have a local amount"),
                    );
                }
                "C(gr)" if reference.temperature_k == 700.0 => {
                    let absolute = row
                        .absolute_error
                        .expect("active graphite must compare")
                        .abs();
                    assert!(
                        absolute <= MAX_GRAPHITE_700_ABSOLUTE_SYSTEM_FRACTION_ERROR_GUARD,
                        "NASA TP-1907 external quality regression at 700 K for C(gr): absolute system-fraction error={absolute:e} exceeds reviewed guard={MAX_GRAPHITE_700_ABSOLUTE_SYSTEM_FRACTION_ERROR_GUARD:e}; source={:e}, local={:e}",
                        row.source_system_fraction,
                        row.kithe_system_fraction
                            .expect("active graphite must have a local amount"),
                    );
                }
                "C(gr)" => {}
                _ if row.source_system_fraction >= 1.0e-3 => {
                    if let Some(relative) = row.relative_error {
                        major_gas_relative_errors.push(relative);
                    }
                }
                _ => {}
            }
        }
        assert!(
            !major_gas_relative_errors.is_empty(),
            "NASA TP-1907 fixture must retain comparable major gas species at {} K",
            reference.temperature_k,
        );
        let max = major_gas_relative_errors
            .iter()
            .map(|error| error.abs())
            .fold(0.0_f64, f64::max);
        let rms = (major_gas_relative_errors
            .iter()
            .map(|error| error * error)
            .sum::<f64>()
            / major_gas_relative_errors.len() as f64)
            .sqrt();
        assert!(
            max <= MAX_MAJOR_GAS_RELATIVE_ERROR_GUARD,
            "NASA TP-1907 external quality regression at {} K: max major-gas relative error={max:e} exceeds reviewed guard={MAX_MAJOR_GAS_RELATIVE_ERROR_GUARD:e}",
            reference.temperature_k,
        );
        assert!(
            rms <= MAX_MAJOR_GAS_RMS_RELATIVE_ERROR_GUARD,
            "NASA TP-1907 external quality regression at {} K: RMS major-gas relative error={rms:e} exceeds reviewed guard={MAX_MAJOR_GAS_RMS_RELATIVE_ERROR_GUARD:e}",
            reference.temperature_k,
        );
    }

    #[test]
    fn i5_nasa_tp1907_external_regression_envelope_preserves_topology_and_boundary() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let dataset = dataset();
        let fixture = fixture();
        let evidence_request = PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())),
            PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
        );

        for reference in dataset.rows() {
            let solution = solve_graphite_pt(&fixture, reference.temperature_k);
            let status = solution
                .phase_status(&evidence_request.candidate_phase)
                .expect("graphite phase must be represented");
            let tpd = graphite_tpd_from_solution(&solution, &evidence_request);
            let expects_graphite = reference
                .condensed_species
                .iter()
                .find(|entry| entry.phase == "C(gr)")
                .expect("frozen graphite row must exist")
                .system_mole_fraction
                > 0.0;
            assert_eq!(
                matches!(status, PhaseStatus::Active | PhaseStatus::Appeared),
                expects_graphite,
                "NASA TP-1907 topology regression at {} K: status={status:?}",
                reference.temperature_k,
            );
            assert!(
                if expects_graphite {
                    tpd < 0.0
                } else {
                    tpd > 0.0
                },
                "NASA TP-1907 TPD-sign regression at {} K: topology expects graphite={expects_graphite}, TPD={tpd:e}",
                reference.temperature_k,
            );
            assert_tp1907_external_regression_envelope(&fixture, reference, &solution);
        }
        let (boundary_temperature_k, _, _) =
            bisect_graphite_tpd_boundary(&fixture, &evidence_request);
        assert!(
            (700.0..=720.0).contains(&boundary_temperature_k),
            "NASA TP-1907 graphite boundary regression: internal TPD root={boundary_temperature_k} K is outside frozen topology bracket [700, 720] K"
        );
        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    #[test]
    fn tp1907_real_reaction_basis_permutation_and_scaling_preserve_fixed_active_solution() {
        let fixture = fixture();
        let prepared = prepared_fixed_tp1907_problem(&fixture, 700.0);
        let baseline = solve_prepared_with_legacy_lm(prepared.clone());

        let source = prepared.reaction_basis().reactions.clone();
        let mut transformed = nalgebra::DMatrix::zeros(source.nrows(), source.ncols());
        let scales = [0.25, 1.5, 3.0, 0.75, 2.0];
        for target_column in 0..source.ncols() {
            let source_column = source.ncols() - 1 - target_column;
            let scale = scales[target_column % scales.len()];
            for row in 0..source.nrows() {
                transformed[(row, target_column)] = scale * source[(row, source_column)];
            }
        }
        let transformed = prepared
            .with_reaction_matrix_for_test(transformed)
            .expect("permuted and rescaled real reaction basis must remain valid");
        let candidate = solve_prepared_with_legacy_lm(transformed);

        for (index, (&expected, &actual)) in
            baseline.moles().iter().zip(candidate.moles()).enumerate()
        {
            let relative =
                (expected - actual).abs() / expected.abs().max(actual.abs()).max(1.0e-12);
            assert!(
                relative <= 1.0e-7,
                "TP-1907 reaction-basis metamorphic mismatch at component {index}: baseline={expected:e}, transformed={actual:e}, relative={relative:e}"
            );
        }
        assert!(
            candidate.validation().max_abs_reaction_affinity <= 1.0e-6,
            "transformed reaction coordinates must satisfy the physical affinity contract"
        );
        assert!(
            candidate.validation().max_abs_element_balance_error <= 1.0e-6,
            "transformed reaction coordinates must preserve elemental balance"
        );
    }

    #[test]
    fn tp1907_reaction_row_normalization_keeps_dimensionless_residuals_invariant() {
        let fixture = fixture();
        let prepared = prepared_fixed_tp1907_problem(&fixture, 700.0);
        let source = prepared.reaction_basis().reactions.clone();
        let mut probe = solve_prepared_with_legacy_lm(prepared.clone())
            .log_moles()
            .to_vec();
        for (index, value) in probe.iter_mut().enumerate() {
            // An off-solution finite probe makes the raw reaction residual
            // observable. At an exact root the same invariant would only
            // compare numerical zero against numerical zero.
            *value += 0.013 * (index + 1) as f64;
        }
        let baseline_raw = prepared
            .residual(&probe)
            .expect("baseline residual must evaluate");
        let baseline_scale = prepared
            .residual_scale()
            .expect("baseline residual scale must evaluate");
        let reaction_count = source.ncols();

        for factor in [1.0e-8, 1.0e-4, 1.0, 1.0e4, 1.0e8, -1.0] {
            let transformed = prepared
                .with_reaction_matrix_for_test(source.clone() * factor)
                .expect("scaled real reaction basis must remain element-conserving");
            let raw = transformed
                .residual(&probe)
                .expect("transformed residual must evaluate");
            let scale = transformed
                .residual_scale()
                .expect("transformed residual scale must evaluate");

            for index in 0..reaction_count {
                let expected_raw = factor * baseline_raw[index];
                let raw_error = (raw[index] - expected_raw).abs();
                assert!(
                    raw_error <= 1.0e-10 * expected_raw.abs().max(1.0),
                    "reaction row {index} raw residual must follow its normalization: factor={factor:e}, expected={expected_raw:e}, actual={:e}, error={raw_error:e}",
                    raw[index],
                );

                let expected_scale = factor.abs() * baseline_scale[index];
                let scale_error = (scale[index] - expected_scale).abs();
                assert!(
                    scale_error <= 1.0e-10 * expected_scale.max(1.0),
                    "reaction row {index} scale must follow its normalization: factor={factor:e}, expected={expected_scale:e}, actual={:e}, error={scale_error:e}",
                    scale[index],
                );

                let baseline_dimensionless = baseline_raw[index] / baseline_scale[index];
                let dimensionless = raw[index] / scale[index];
                let expected_dimensionless = factor.signum() * baseline_dimensionless;
                assert!(
                    (dimensionless - expected_dimensionless).abs()
                        <= 1.0e-10 * expected_dimensionless.abs().max(1.0),
                    "reaction row {index} dimensionless residual must be invariant up to a sign reversal: factor={factor:e}, expected={expected_dimensionless:e}, actual={dimensionless:e}",
                );
            }
        }
    }

    #[test]
    fn tp1907_fixed_active_large_inventory_accepts_scaled_accepted_seed() {
        let fixture = fixture();
        let baseline_problem = prepared_fixed_tp1907_problem(&fixture, 700.0);
        let baseline = solve_prepared_with_legacy_lm(baseline_problem);
        let factor = 1.0e4_f64;
        let scaled_seed = LogMolesInitialGuess::new(
            baseline
                .log_moles()
                .iter()
                .map(|value| value + factor.ln())
                .collect(),
        )
        .expect("scaled accepted log-mole seed must remain finite");
        let mut runner = PreparedEquilibriumRunner::new(
            prepared_fixed_tp1907_problem_at_scale(&fixture, 700.0, factor),
            Vec::new(),
        )
        .expect("scaled fixed-active runner must construct");
        let settings = runner.configure();
        settings.solver = Solvers::LM;
        settings.solver_params.tol = 1.0e-10;
        settings.solver_params.max_iter = 500;
        let large = runner
            .solve_from_seed(scaled_seed)
            .expect("scaled accepted seed must remain in the fixed-active basin")
            .solution;

        for (index, (&expected, &actual)) in baseline.moles().iter().zip(large.moles()).enumerate()
        {
            let recovered = actual / factor;
            let relative =
                (expected - recovered).abs() / expected.abs().max(recovered.abs()).max(1.0e-30);
            assert!(
                relative <= 1.0e-7,
                "fixed-active large-scale continuation mismatch at component {index}: baseline={expected:e}, recovered={recovered:e}, relative={relative:e}"
            );
        }
    }

    #[test]
    #[ignore = "release diagnostic separating fresh log-seed conditioning from phase control"]
    fn i5_tp1907_large_inventory_seed_basin_matrix_diagnostic() {
        let fixture = fixture();
        let factor = 1.0e4_f64;
        const ABSOLUTE_TRACE_FLOOR: f64 = 1.0e-30;
        const LOW_ABSOLUTE_TRACE_FLOOR: f64 = 1.0e-36;
        let relative_trace_policy = TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
            fraction: 1.0e-30,
            minimum_floor: ABSOLUTE_TRACE_FLOOR,
        };

        println!("TP-1907 fixed-active 1e4 fresh-seed basin matrix");
        println!(
            "routes: A=fresh production, B=transformed accepted oracle, C=low-absolute-trace diagnostic, D=relative-trace diagnostic"
        );
        for temperature_k in [680.0, 700.0, 720.0] {
            if temperature_k == 720.0 {
                // Graphite is physically absent here. An all-active
                // fixed-set solve is therefore not the same formulation as
                // the bounded production state and may have no compatible
                // root. The scale oracle must instead keep graphite inactive
                // while retaining it as a TPD candidate.
                let bounded = solve_graphite_pt(&fixture, temperature_k);
                assert_eq!(
                    bounded.phase_status(&PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()))),
                    Some(PhaseStatus::Inactive),
                    "720 K must remain on the gas-only production topology"
                );
                let baseline_tpd = graphite_tpd_from_solution(
                    &bounded,
                    &PurePhaseProductionEvidenceRequest::new(
                        PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())),
                        PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
                    ),
                );
                assert!(
                    baseline_tpd > 0.0,
                    "the 720 K canonical gas-only state must reject graphite by positive TPD, got {baseline_tpd:e}"
                );

                let scaled_problem =
                    prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, factor);
                let full_problem = scaled_problem.problem();
                let interior_seed = element_feasible_interior_seed(&scaled_problem, &[true, false]);
                print_interior_seed("I element-feasible interior", &interior_seed);
                print_seed_statistics(
                    "I element-feasible interior",
                    log_seed_statistics(&interior_seed.log_moles, ABSOLUTE_TRACE_FLOOR),
                );
                print_initial_residual_diagnostics(
                    "I interior residual blocks",
                    &scaled_problem,
                    &interior_seed.log_moles,
                );
                let transformed_seed = LogMolesInitialGuess::new(
                    bounded
                        .numerical_component_moles()
                        .iter()
                        .enumerate()
                        .map(|(index, &moles)| {
                            // The 720 K baseline has only gas physically
                            // active. Preserve the graphite trace coordinate
                            // from the scaled problem instead of incorrectly
                            // translating it into an active graphite amount.
                            if scaled_problem.species_phase()[index] == 0 {
                                moles.ln() + factor.ln()
                            } else {
                                full_problem.initial_log_moles().as_slice()[index]
                            }
                        })
                        .collect(),
                )
                .expect("scaled reduced gas-only transformed seed must be finite");
                let fresh = solve_resolved_pt(
                    ResolvedPhaseEquilibriumRequest::new(
                        fixture.resolved(),
                        fixture
                            .conditions_at(temperature_k, 101_325.0)
                            .expect("720 K large-inventory conditions must be supported"),
                        scaled_initial_composition(&fixture, factor),
                    )
                    .with_phase_control_policy(PhaseControlPolicy::default())
                    .with_solve_options(strict_extensivity_solve_options()),
                );
                let transformed = solve_reduced_gas_only_from_seed(&scaled_problem, transformed_seed)
                    .expect("the transformed gas-only 720 K seed must prove the scaled physical state exists");
                let interior = solve_reduced_gas_only_from_seed(
                    &scaled_problem,
                    interior_seed.log_moles.clone(),
                );

                assert_eq!(
                    transformed.phase_statuses,
                    vec![PhaseStatus::Active, PhaseStatus::Inactive],
                    "720 K transformed solve must preserve the gas-only topology"
                );
                for (index, (&baseline_moles, &scaled_moles)) in bounded
                    .component_moles()
                    .iter()
                    .zip(transformed.solution.moles())
                    .enumerate()
                {
                    if scaled_problem.species_phase()[index] != 0 {
                        continue;
                    }
                    let recovered = scaled_moles / factor;
                    let relative = (baseline_moles - recovered).abs()
                        / baseline_moles.abs().max(recovered.abs()).max(1.0e-30);
                    assert!(
                        relative <= 1.0e-6,
                        "720 K gas-only large-scale mismatch at component {index}: baseline={baseline_moles:e}, recovered={recovered:e}, relative={relative:e}"
                    );
                }
                let baseline_gas = gas_mole_fractions_by_identity(&bounded);
                let gas_indices = scaled_problem
                    .species_phase()
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &phase)| (phase == 0).then_some(index))
                    .collect::<Vec<_>>();
                let scaled_gas_total = gas_indices
                    .iter()
                    .map(|&index| transformed.solution.moles()[index])
                    .sum::<f64>();
                for (index, component) in full_problem.components().iter().enumerate() {
                    if scaled_problem.species_phase()[index] != 0 {
                        continue;
                    }
                    let actual = transformed.solution.moles()[index] / scaled_gas_total;
                    let expected = baseline_gas[&component.id().label()];
                    assert!(
                        (expected - actual).abs() <= 1.0e-7,
                        "720 K gas composition must be extensive for {}: baseline={expected:e}, scaled={actual:e}",
                        component.id().label(),
                    );
                }
                assert!(
                    transformed.solution.validation().max_abs_reaction_affinity <= 1.0e-6
                        && transformed
                            .solution
                            .validation()
                            .max_abs_element_balance_error
                            <= 1.0e-6,
                    "720 K transformed gas-only solve must retain the accepted nonlinear contract"
                );
                let gas_only = PhaseSet::from_policy(
                    &InitialPhaseSet::Explicit {
                        active: vec![
                            PhaseIndex::new(0, full_problem.phases().len())
                                .expect("gas phase index must be valid"),
                        ],
                        excluded: Vec::new(),
                    },
                    &[true, false],
                )
                .expect("720 K gas-only phase set must be valid");
                let stability = compute_phase_stability_reports(
                    transformed.solution.log_moles(),
                    full_problem.gibbs(),
                    full_problem.phases(),
                    scaled_problem.species_phase(),
                    full_problem.element_composition(),
                    full_problem.conditions().temperature(),
                    full_problem.conditions().pressure(),
                    full_problem.conditions().reference_pressure(),
                    &gas_only,
                )
                .expect("720 K graphite stability evidence must evaluate");
                let scaled_tpd = stability[1]
                    .minimum_tpd
                    .expect("inactive graphite must have canonical TPD evidence");
                let tpd_relative = (baseline_tpd - scaled_tpd).abs()
                    / baseline_tpd.abs().max(scaled_tpd.abs()).max(1.0);
                assert!(
                    // The baseline arrives through ordinary bounded phase
                    // control while this diagnostic enters the same reduced
                    // active set directly. Their accepted TPDs are
                    // intensive and agree to a few parts per million; the
                    // physics-critical invariant is the shared positive
                    // sign, not bitwise equality of two nonlinear routes.
                    scaled_tpd > 0.0 && tpd_relative <= 1.0e-5,
                    "720 K graphite TPD must stay positive and intensive: baseline={baseline_tpd:e}, scaled={scaled_tpd:e}, relative={tpd_relative:e}"
                );
                println!("T={temperature_k:.1} K factor={factor:.1e}");
                println!(
                    "  A fresh production            {}",
                    if fresh.is_ok() { "OK" } else { "FAILED" }
                );
                println!(
                    "  B transformed reduced gas    OK     residual={:.3e} balance={:.3e} TPD_graphite={scaled_tpd:.3e}",
                    transformed.solution.validation().residual_l2_norm,
                    transformed
                        .solution
                        .validation()
                        .max_abs_element_balance_error,
                );
                match &interior {
                    Ok(outcome) => println!(
                        "  I element-feasible interior  OK     residual={:.3e} balance={:.3e} topology={:?}",
                        outcome.solution.validation().residual_l2_norm,
                        outcome.solution.validation().max_abs_element_balance_error,
                        outcome.phase_statuses,
                    ),
                    Err(error) => println!(
                        "  I element-feasible interior  FAILED kind={:?} attempts={} iterations={} error={error}",
                        error.kind(),
                        error.started_backend_attempts(),
                        error.nonlinear_iterations(),
                    ),
                }
                let input_anchored = element_feasible_interior_seed_with_settings(
                    &scaled_problem,
                    &[true, false],
                    InteriorSeedSettings {
                        target: InteriorSeedTarget::InputAnchored,
                        ..InteriorSeedSettings::default()
                    },
                );
                print_interior_seed("I input-anchored interior", &input_anchored);
                print_initial_residual_diagnostics(
                    "I input-anchored residual blocks",
                    &scaled_problem,
                    &input_anchored.log_moles,
                );
                match solve_reduced_gas_only_from_seed(&scaled_problem, input_anchored.log_moles) {
                    Ok(outcome) => println!(
                        "  I input-anchored interior  OK     residual={:.3e} balance={:.3e} topology={:?}",
                        outcome.solution.validation().residual_l2_norm,
                        outcome.solution.validation().max_abs_element_balance_error,
                        outcome.phase_statuses,
                    ),
                    Err(error) => println!(
                        "  I input-anchored interior  FAILED kind={:?} attempts={} iterations={} error={error:?}",
                        error.kind(),
                        error.started_backend_attempts(),
                        error.nonlinear_iterations(),
                    ),
                }
                for fraction in [1.0e-4, 1.0e-8] {
                    let variant = element_feasible_interior_seed_with_settings(
                        &scaled_problem,
                        &[true, false],
                        InteriorSeedSettings {
                            requested_minimum_fraction: fraction,
                            ..InteriorSeedSettings::default()
                        },
                    );
                    let label = format!("I interior fraction {fraction:.0e}");
                    print_interior_seed(&label, &variant);
                    print_initial_residual_diagnostics(
                        &format!("{label} residual blocks"),
                        &scaled_problem,
                        &variant.log_moles,
                    );
                    match solve_reduced_gas_only_from_seed(&scaled_problem, variant.log_moles) {
                        Ok(outcome) => println!(
                            "  {label:27} OK     residual={:.3e} balance={:.3e} topology={:?}",
                            outcome.solution.validation().residual_l2_norm,
                            outcome.solution.validation().max_abs_element_balance_error,
                            outcome.phase_statuses,
                        ),
                        Err(error) => println!(
                            "  {label:27} FAILED kind={:?} attempts={} iterations={} error={error:?}",
                            error.kind(),
                            error.started_backend_attempts(),
                            error.nonlinear_iterations(),
                        ),
                    }
                }
                println!(
                    "  C/D trace variants            NOT APPLICABLE: 680/700 already reject the trace-floor hypothesis"
                );
                let classification = if fresh.is_ok() {
                    ScaleConditioningClassification::StrictScaleInvariant
                } else {
                    ScaleConditioningClassification::FreshSeedBasinLimited
                };
                println!("  classification              {classification:?}");
                continue;
            }
            let baseline_problem = prepared_fixed_tp1907_problem(&fixture, temperature_k);
            let baseline = solve_prepared_with_legacy_lm(baseline_problem);
            let scaled_problem =
                prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, factor);
            let fresh_seed = scaled_problem.problem().initial_log_moles().clone();
            let transformed_seed = LogMolesInitialGuess::new(
                baseline
                    .log_moles()
                    .iter()
                    .map(|value| value + factor.ln())
                    .collect(),
            )
            .expect("scaled accepted log-mole seed must remain finite");
            let relative_trace_seed = LogMolesInitialGuess::from_moles_with_policy(
                scaled_problem.problem().initial_moles(),
                relative_trace_policy,
            )
            .expect("relative-trace diagnostic seed must be valid");
            let low_trace_seed = LogMolesInitialGuess::from_moles(
                scaled_problem.problem().initial_moles(),
                LOW_ABSOLUTE_TRACE_FLOOR,
            )
            .expect("low-absolute-trace diagnostic seed must be valid");
            let (route_e_seed, route_e_scale, input_total, fresh_active_total) =
                input_derived_active_normalized_seed(&scaled_problem);
            let interior_seed = element_feasible_interior_seed(&scaled_problem, &[true, true]);

            println!("T={temperature_k:.1} K factor={factor:.1e}");
            print_seed_statistics(
                "A fresh production",
                log_seed_statistics(&fresh_seed, ABSOLUTE_TRACE_FLOOR),
            );
            print_seed_statistics(
                "B transformed accepted",
                log_seed_statistics(&transformed_seed, ABSOLUTE_TRACE_FLOOR),
            );
            let relative_floor = relative_trace_policy
                .trace_floor_for(scaled_problem.problem().initial_moles())
                .expect("relative trace floor must be valid");
            print_seed_statistics(
                "C low absolute trace",
                log_seed_statistics(&low_trace_seed, LOW_ABSOLUTE_TRACE_FLOOR),
            );
            print_seed_statistics(
                "D relative trace",
                log_seed_statistics(&relative_trace_seed, relative_floor),
            );
            print_seed_statistics(
                "E input-normalized",
                log_seed_statistics(&route_e_seed, ABSOLUTE_TRACE_FLOOR),
            );
            print_interior_seed("I element-feasible interior", &interior_seed);
            print_seed_statistics(
                "I element-feasible interior",
                log_seed_statistics(&interior_seed.log_moles, ABSOLUTE_TRACE_FLOOR),
            );
            println!(
                "  Route E active totals        input={input_total:.6e} fresh={fresh_active_total:.6e} scale={route_e_scale:.6e}"
            );
            for (index, (&ordinary, &route_e)) in fresh_seed
                .as_slice()
                .iter()
                .zip(route_e_seed.as_slice())
                .enumerate()
            {
                assert!(
                    (ordinary - route_e).abs() <= 1.0e-12,
                    "Route E must remain observationally identical to the ordinary seed until instrumentation proves an input-scale mismatch at component {index}"
                );
            }
            assert!(
                (route_e_scale - 1.0).abs() <= 1.0e-12,
                "the current fresh seed already carries the physical active inventory scale; Route E must not claim a missing global correction"
            );
            print_component_seed_comparison(
                temperature_k,
                factor,
                &prepared_fixed_tp1907_problem(&fixture, temperature_k),
                &scaled_problem,
                &baseline,
            );
            print_initial_residual_diagnostics(
                "A fresh residual blocks",
                &scaled_problem,
                &fresh_seed,
            );
            print_initial_residual_diagnostics(
                "B transformed residual blocks",
                &scaled_problem,
                &transformed_seed,
            );
            print_initial_residual_diagnostics(
                "I interior residual blocks",
                &scaled_problem,
                &interior_seed.log_moles,
            );

            let fresh = strict_fixed_active_runner(scaled_problem.clone()).solve();
            let transformed = strict_fixed_active_runner(scaled_problem.clone())
                .solve_from_seed(transformed_seed);
            let interior = strict_fixed_active_runner(scaled_problem.clone())
                .solve_from_seed(interior_seed.log_moles.clone());
            let low_trace =
                strict_fixed_active_runner(scaled_problem.clone()).solve_from_seed(low_trace_seed);
            let relative_trace = strict_fixed_active_runner(scaled_problem.clone())
                .solve_from_seed(relative_trace_seed);

            let describe =
                |label: &str, result: &Result<PreparedSolveOutcome, ReactionExtentError>| {
                    match result {
                        Ok(outcome) => println!(
                            "  {label:27} OK     residual={:.3e} balance={:.3e} attempts={} iterations={}",
                            outcome.solution.validation().residual_l2_norm,
                            outcome.solution.validation().max_abs_element_balance_error,
                            outcome.solve_report.attempts.len(),
                            outcome.solve_report.nonlinear_iterations(),
                        ),
                        Err(error) => println!(
                            "  {label:27} FAILED kind={:?} attempts={} iterations={} error={error:?}",
                            error.kind(),
                            error.started_backend_attempts(),
                            error.nonlinear_iterations(),
                        ),
                    }
                };
            describe("A fresh production", &fresh);
            describe("B transformed accepted", &transformed);
            describe("I element-feasible interior", &interior);
            let input_anchored = element_feasible_interior_seed_with_settings(
                &scaled_problem,
                &[true, true],
                InteriorSeedSettings {
                    target: InteriorSeedTarget::InputAnchored,
                    ..InteriorSeedSettings::default()
                },
            );
            print_interior_seed("I input-anchored interior", &input_anchored);
            print_initial_residual_diagnostics(
                "I input-anchored residual blocks",
                &scaled_problem,
                &input_anchored.log_moles,
            );
            let input_anchored_result = strict_fixed_active_runner(scaled_problem.clone())
                .solve_from_seed(input_anchored.log_moles);
            describe("I input-anchored interior", &input_anchored_result);
            for fraction in [1.0e-4, 1.0e-8] {
                let variant = element_feasible_interior_seed_with_settings(
                    &scaled_problem,
                    &[true, true],
                    InteriorSeedSettings {
                        requested_minimum_fraction: fraction,
                        ..InteriorSeedSettings::default()
                    },
                );
                let label = format!("I interior fraction {fraction:.0e}");
                print_interior_seed(&label, &variant);
                print_initial_residual_diagnostics(
                    &format!("{label} residual blocks"),
                    &scaled_problem,
                    &variant.log_moles,
                );
                let result = strict_fixed_active_runner(scaled_problem.clone())
                    .solve_from_seed(variant.log_moles);
                describe(&label, &result);
            }
            describe("C low absolute trace", &low_trace);
            describe("D relative trace", &relative_trace);
            describe("E input-normalized", &fresh);

            let transformed = transformed.expect(
                "the transformed accepted seed is the physical-existence oracle for this scale diagnostic",
            );
            assert_fixed_active_scaled_solution(
                &baseline,
                &transformed.solution,
                factor,
                temperature_k,
            );
            let classification = if fresh.is_ok() {
                ScaleConditioningClassification::StrictScaleInvariant
            } else {
                ScaleConditioningClassification::FreshSeedBasinLimited
            };
            println!("  classification              {classification:?}");
        }
    }

    #[test]
    #[ignore = "release diagnostic for exact test-only extensive normalization"]
    fn i5_tp1907_large_inventory_extensive_normalization_matrix() {
        let fixture = fixture();
        let inventory_factor = 1.0e4_f64;

        println!("TP-1907 exact extensive-normalization matrix");
        println!(
            "T K    input factor  route  internal total  topology       status   max n/B error  max gas-x error  TPD error"
        );

        for temperature_k in [680.0, 700.0] {
            let unit = solve_prepared_with_legacy_lm(prepared_fixed_tp1907_problem(
                &fixture,
                temperature_k,
            ));
            let physical =
                prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, inventory_factor);
            let transformed_seed = LogMolesInitialGuess::new(
                unit.log_moles()
                    .iter()
                    .map(|value| value + inventory_factor.ln())
                    .collect(),
            )
            .expect("large transformed fixed-active oracle seed must be finite");
            let fresh = strict_fixed_active_runner(physical.clone()).solve();
            assert!(
                fresh.is_err(),
                "the established large fresh fixed-active witness must still fail at {temperature_k} K"
            );
            let oracle = strict_fixed_active_runner(physical.clone())
                .solve_from_seed(transformed_seed)
                .expect("transformed accepted seed must remain the large physical oracle");

            let (normalized, normalization_factor) =
                normalized_fixed_problem_from_input(&fixture, temperature_k, &physical);
            let internal_total = normalized.problem().initial_moles().iter().sum::<f64>();
            assert!(
                (internal_total - 1.0).abs() <= 1.0e-12,
                "the test normalization must map active physical input to one mole"
            );
            let normalized_outcome = strict_fixed_active_runner(normalized.clone())
                .solve()
                .expect("ordinary fresh solve of the normalized equivalent must accept");

            let mut max_mole_error = 0.0_f64;
            let mut max_gas_fraction_error = 0.0_f64;
            let normalized_gas_total = normalized_outcome
                .solution
                .moles()
                .iter()
                .enumerate()
                .filter_map(|(index, &amount)| {
                    (normalized.species_phase()[index] == 0).then_some(amount)
                })
                .sum::<f64>();
            let oracle_gas_total = oracle
                .solution
                .moles()
                .iter()
                .enumerate()
                .filter_map(|(index, &amount)| {
                    (physical.species_phase()[index] == 0).then_some(amount)
                })
                .sum::<f64>();
            for (index, (&normalized_moles, &oracle_moles)) in normalized_outcome
                .solution
                .moles()
                .iter()
                .zip(oracle.solution.moles())
                .enumerate()
            {
                let reconstructed = normalized_moles * normalization_factor;
                let relative = (reconstructed - oracle_moles).abs()
                    / reconstructed.abs().max(oracle_moles.abs()).max(1.0e-30);
                max_mole_error = max_mole_error.max(relative);
                if physical.species_phase()[index] == 0 {
                    let normalized_fraction = normalized_moles / normalized_gas_total;
                    let oracle_fraction = oracle_moles / oracle_gas_total;
                    max_gas_fraction_error =
                        max_gas_fraction_error.max((normalized_fraction - oracle_fraction).abs());
                }
            }
            assert!(
                max_mole_error <= 1.0e-6 && max_gas_fraction_error <= 1.0e-7,
                "normalized {temperature_k} K state must reconstruct the large oracle: moles={max_mole_error:e}, gas_x={max_gas_fraction_error:e}"
            );
            assert!(
                normalized_outcome
                    .solution
                    .validation()
                    .max_abs_reaction_affinity
                    <= 1.0e-6
                    && normalized_outcome
                        .solution
                        .validation()
                        .max_abs_element_balance_error
                        <= 1.0e-8,
                "normalized {temperature_k} K solve must retain the ordinary nonlinear contract"
            );
            println!(
                "{temperature_k:5.1}  {inventory_factor:12.1e}  N      {internal_total:14.6e}  gas+graphite  OK       {max_mole_error:12.3e}  {max_gas_fraction_error:15.3e}  -"
            );
        }

        // At 720 K graphite is a stability candidate, not an active unknown.
        // The normalization experiment must retain that reduced topology.
        let temperature_k = 720.0;
        let physical =
            prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, inventory_factor);
        let unit_bounded = solve_graphite_pt(&fixture, temperature_k);
        let baseline_tpd = graphite_tpd_from_solution(
            &unit_bounded,
            &PurePhaseProductionEvidenceRequest::new(
                PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())),
                PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
            ),
        );
        let oracle_seed = LogMolesInitialGuess::new(
            unit_bounded
                .numerical_component_moles()
                .iter()
                .enumerate()
                .map(|(index, &amount)| {
                    if physical.species_phase()[index] == 0 {
                        amount.ln() + inventory_factor.ln()
                    } else {
                        physical.problem().initial_log_moles().as_slice()[index]
                    }
                })
                .collect(),
        )
        .expect("large gas-only transformed oracle seed must be finite");
        let fresh = solve_reduced_gas_only_from_seed(
            &physical,
            physical.problem().initial_log_moles().clone(),
        );
        assert!(
            fresh.is_err(),
            "the established large fresh gas-only witness must still fail at 720 K"
        );
        let oracle = solve_reduced_gas_only_from_seed(&physical, oracle_seed)
            .expect("transformed gas-only seed must remain the physical oracle");
        let (normalized, normalization_factor) =
            normalized_fixed_problem_from_input(&fixture, temperature_k, &physical);
        let internal_total = normalized.problem().initial_moles().iter().sum::<f64>();
        let normalized_outcome = solve_reduced_gas_only_from_seed(
            &normalized,
            normalized.problem().initial_log_moles().clone(),
        )
        .expect("ordinary fresh normalized gas-only solve must accept");
        assert_eq!(
            normalized_outcome.phase_statuses,
            vec![PhaseStatus::Active, PhaseStatus::Inactive],
            "720 K normalization must retain gas-only topology"
        );

        let mut max_mole_error = 0.0_f64;
        let mut max_gas_fraction_error = 0.0_f64;
        let normalized_gas_total = normalized_outcome
            .solution
            .moles()
            .iter()
            .enumerate()
            .filter_map(|(index, &amount)| {
                (normalized.species_phase()[index] == 0).then_some(amount)
            })
            .sum::<f64>();
        let oracle_gas_total = oracle
            .solution
            .moles()
            .iter()
            .enumerate()
            .filter_map(|(index, &amount)| (physical.species_phase()[index] == 0).then_some(amount))
            .sum::<f64>();
        for (index, (&normalized_moles, &oracle_moles)) in normalized_outcome
            .solution
            .moles()
            .iter()
            .zip(oracle.solution.moles())
            .enumerate()
        {
            if physical.species_phase()[index] != 0 {
                continue;
            }
            let reconstructed = normalized_moles * normalization_factor;
            let relative = (reconstructed - oracle_moles).abs()
                / reconstructed.abs().max(oracle_moles.abs()).max(1.0e-30);
            max_mole_error = max_mole_error.max(relative);
            max_gas_fraction_error = max_gas_fraction_error.max(
                (normalized_moles / normalized_gas_total - oracle_moles / oracle_gas_total).abs(),
            );
        }
        let gas_only = PhaseSet::from_policy(
            &InitialPhaseSet::Explicit {
                active: vec![
                    PhaseIndex::new(0, normalized.problem().phases().len())
                        .expect("gas phase index must be valid"),
                ],
                excluded: Vec::new(),
            },
            &[true, false],
        )
        .expect("normalized 720 K gas-only phase set must be valid");
        let stability = compute_phase_stability_reports(
            normalized_outcome.solution.log_moles(),
            normalized.problem().gibbs(),
            normalized.problem().phases(),
            normalized.species_phase(),
            normalized.problem().element_composition(),
            normalized.problem().conditions().temperature(),
            normalized.problem().conditions().pressure(),
            normalized.problem().conditions().reference_pressure(),
            &gas_only,
        )
        .expect("normalized 720 K graphite stability evidence must evaluate");
        let normalized_tpd = stability[1]
            .minimum_tpd
            .expect("normalized inactive graphite must have TPD evidence");
        let tpd_relative = (normalized_tpd - baseline_tpd).abs()
            / normalized_tpd.abs().max(baseline_tpd.abs()).max(1.0);
        assert!(
            max_mole_error <= 1.0e-6
                && max_gas_fraction_error <= 1.0e-7
                && normalized_tpd > 0.0
                && tpd_relative <= 1.0e-5,
            "normalized 720 K state must reconstruct gas-only oracle: moles={max_mole_error:e}, gas_x={max_gas_fraction_error:e}, TPD={normalized_tpd:e}, relative={tpd_relative:e}"
        );
        println!(
            "{temperature_k:5.1}  {inventory_factor:12.1e}  N      {internal_total:14.6e}  gas            OK       {max_mole_error:12.3e}  {max_gas_fraction_error:15.3e}  {tpd_relative:9.3e}"
        );
        assert!(
            (normalization_factor - physical.problem().initial_moles().iter().sum::<f64>()).abs()
                <= 1.0e-8,
            "normalization factor must derive solely from physical input"
        );
    }

    #[test]
    #[ignore = "release metamorphic matrix for exact extensive normalization"]
    fn i5_tp1907_extensive_normalization_scale_metamorphism() {
        let fixture = fixture();
        let factors = [1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4];

        println!("TP-1907 extensive-normalization scale metamorphism");
        println!(
            "T K    physical factor  internal total  topology       max n/B error  max gas-x error  TPD error"
        );

        for temperature_k in [680.0, 700.0] {
            let unit = solve_prepared_with_legacy_lm(prepared_fixed_tp1907_problem(
                &fixture,
                temperature_k,
            ));
            for factor in factors {
                let physical =
                    prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, factor);
                let oracle_seed = LogMolesInitialGuess::new(
                    unit.log_moles()
                        .iter()
                        .map(|value| value + factor.ln())
                        .collect(),
                )
                .expect("scaled fixed-active oracle seed must be finite");
                let oracle = strict_fixed_active_runner(physical.clone())
                    .solve_from_seed(oracle_seed)
                    .expect("scaled fixed-active oracle must accept");
                let (normalized, normalization_factor) =
                    normalized_fixed_problem_from_input(&fixture, temperature_k, &physical);
                let internal_total = normalized.problem().initial_moles().iter().sum::<f64>();
                let normalized_outcome = strict_fixed_active_runner(normalized.clone())
                    .solve()
                    .expect("ordinary fresh normalized fixed-active solve must accept");
                let (mole_error, gas_fraction_error) = reconstructed_mole_and_gas_fraction_errors(
                    normalized_outcome.solution.moles(),
                    normalization_factor,
                    oracle.solution.moles(),
                    physical.species_phase(),
                );
                assert!(
                    (internal_total - 1.0).abs() <= 1.0e-12
                        && mole_error <= 1.0e-6
                        && gas_fraction_error <= 1.0e-7,
                    "normalized {temperature_k} K factor={factor:e} must reconstruct fixed-active oracle: internal={internal_total:e}, n={mole_error:e}, gas_x={gas_fraction_error:e}"
                );
                println!(
                    "{temperature_k:5.1}  {factor:14.1e}  {internal_total:14.6e}  gas+graphite  {mole_error:12.3e}  {gas_fraction_error:15.3e}  -"
                );
            }
        }

        let temperature_k = 720.0;
        let unit_bounded = solve_graphite_pt(&fixture, temperature_k);
        let baseline_tpd = graphite_tpd_from_solution(
            &unit_bounded,
            &PurePhaseProductionEvidenceRequest::new(
                PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())),
                PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
            ),
        );
        for factor in factors {
            let physical = prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, factor);
            let oracle_seed = LogMolesInitialGuess::new(
                unit_bounded
                    .numerical_component_moles()
                    .iter()
                    .enumerate()
                    .map(|(index, &amount)| {
                        if physical.species_phase()[index] == 0 {
                            amount.ln() + factor.ln()
                        } else {
                            physical.problem().initial_log_moles().as_slice()[index]
                        }
                    })
                    .collect(),
            )
            .expect("scaled gas-only oracle seed must be finite");
            let oracle = solve_reduced_gas_only_from_seed(&physical, oracle_seed)
                .expect("scaled gas-only oracle must accept");
            let (normalized, normalization_factor) =
                normalized_fixed_problem_from_input(&fixture, temperature_k, &physical);
            let internal_total = normalized.problem().initial_moles().iter().sum::<f64>();
            let normalized_outcome = solve_reduced_gas_only_from_seed(
                &normalized,
                normalized.problem().initial_log_moles().clone(),
            )
            .expect("ordinary fresh normalized gas-only solve must accept");
            assert_eq!(
                normalized_outcome.phase_statuses,
                vec![PhaseStatus::Active, PhaseStatus::Inactive],
                "normalized 720 K factor={factor:e} must retain gas-only topology"
            );
            let (mole_error, gas_fraction_error) = reconstructed_mole_and_gas_fraction_errors(
                normalized_outcome.solution.moles(),
                normalization_factor,
                oracle.solution.moles(),
                physical.species_phase(),
            );
            let normalized_tpd = graphite_tpd_from_prepared_gas_only(
                &normalized,
                normalized_outcome.solution.log_moles(),
            );
            let tpd_relative = (normalized_tpd - baseline_tpd).abs()
                / normalized_tpd.abs().max(baseline_tpd.abs()).max(1.0);
            assert!(
                (internal_total - 1.0).abs() <= 1.0e-12
                    && mole_error <= 1.0e-6
                    && gas_fraction_error <= 1.0e-7
                    && normalized_tpd > 0.0
                    && tpd_relative <= 1.0e-5,
                "normalized 720 K factor={factor:e} must reconstruct gas-only oracle: internal={internal_total:e}, n={mole_error:e}, gas_x={gas_fraction_error:e}, tpd={normalized_tpd:e}, tpd_rel={tpd_relative:e}"
            );
            println!(
                "{temperature_k:5.1}  {factor:14.1e}  {internal_total:14.6e}  gas            {mole_error:12.3e}  {gas_fraction_error:15.3e}  {tpd_relative:9.3e}"
            );
        }
    }

    #[test]
    fn tp1907_element_feasible_interior_seed_is_extensive_and_phase_scoped() {
        let fixture = fixture();
        for (temperature_k, active_phase_mask) in [
            (680.0, vec![true, true]),
            (700.0, vec![true, true]),
            (720.0, vec![true, false]),
        ] {
            let baseline_problem =
                prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, 1.0);
            for target in [
                InteriorSeedTarget::Uniform,
                InteriorSeedTarget::InputAnchored,
            ] {
                let settings = InteriorSeedSettings {
                    target,
                    ..InteriorSeedSettings::default()
                };
                let baseline = element_feasible_interior_seed_with_settings(
                    &baseline_problem,
                    &active_phase_mask,
                    settings,
                );
                assert_eq!(
                    baseline.availability,
                    InteriorSeedAvailability::FullInteriorAvailable,
                    "TP-1907 {temperature_k} K must support the initial diagnostic interior fraction for {target:?}"
                );

                for factor in [1.0e-4, 1.0, 1.0e4] {
                    let prepared =
                        prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, factor);
                    let candidate = element_feasible_interior_seed_with_settings(
                        &prepared,
                        &active_phase_mask,
                        settings,
                    );
                    assert_eq!(candidate.availability, baseline.availability);
                    for (index, (&expected, &actual)) in baseline_problem
                        .problem()
                        .initial_moles()
                        .iter()
                        .zip(prepared.problem().initial_moles().iter())
                        .enumerate()
                    {
                        let recovered = actual / factor;
                        let relative = (expected - recovered).abs()
                            / expected.abs().max(recovered.abs()).max(1.0e-30);
                        assert!(
                            relative <= 1.0e-12,
                            "fixture input must itself scale at {temperature_k} K, factor={factor:e}, component {index}: baseline={expected:e}, recovered={recovered:e}, relative={relative:e}"
                        );
                    }
                    for (index, (&expected, &actual)) in baseline
                        .physical_moles
                        .iter()
                        .zip(candidate.physical_moles.iter())
                        .enumerate()
                    {
                        let active = active_phase_mask[prepared.species_phase()[index]];
                        if !active {
                            assert_eq!(
                                actual, 0.0,
                                "inactive phase species {index} must remain physically absent at {temperature_k} K"
                            );
                            continue;
                        }
                        let recovered = actual / factor;
                        let relative = (expected - recovered).abs()
                            / expected.abs().max(recovered.abs()).max(1.0e-30);
                        assert!(
                            relative <= 1.0e-9,
                            "{target:?} interior seed must be extensive at {temperature_k} K, factor={factor:e}, component {index}: baseline={expected:e}, recovered={recovered:e}, relative={relative:e}"
                        );
                        assert!(
                            actual > 0.0,
                            "active component {index} must receive a positive physical interior seed at {temperature_k} K"
                        );
                    }
                    if temperature_k == 720.0 {
                        let graphite = prepared
                            .species_phase()
                            .iter()
                            .position(|&phase| phase == 1)
                            .expect("TP-1907 graphite component must be present");
                        assert_eq!(candidate.physical_moles[graphite], 0.0);
                        assert!(
                            (candidate.log_moles.as_slice()[graphite] - 1.0e-30_f64.ln()).abs()
                                <= 1.0e-12,
                            "inactive graphite must retain only the declared numerical trace"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn tp1907_real_bounded_phase_control_is_invariant_to_small_trace_floors() {
        let fixture = fixture();
        for temperature_k in [700.0, 720.0] {
            let baseline = solve_graphite_pt_with_trace_floor(&fixture, temperature_k, 1.0e-30);
            let baseline_status = baseline
                .phase_status(&PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())))
                .expect("graphite phase must be represented");
            for floor in [1.0e-24, 1.0e-36] {
                let candidate = solve_graphite_pt_with_trace_floor(&fixture, temperature_k, floor);
                assert_eq!(
                    candidate.phase_status(&PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()))),
                    Some(baseline_status),
                    "TP-1907 graphite topology must not depend on the numerical trace floor at {temperature_k} K"
                );
                for (index, (&expected, &actual)) in baseline
                    .component_moles()
                    .iter()
                    .zip(candidate.component_moles())
                    .enumerate()
                {
                    let relative =
                        (expected - actual).abs() / expected.abs().max(actual.abs()).max(1.0e-12);
                    assert!(
                        relative <= 1.0e-7 || (expected - actual).abs() <= 1.0e-12,
                        "TP-1907 trace-floor mismatch at {temperature_k} K, component {index}: baseline={expected:e}, candidate={actual:e}, relative={relative:e}"
                    );
                }
            }
        }
    }

    #[test]
    fn tp1907_real_pt_solution_is_extensive_under_inventory_scaling() {
        let fixture = fixture();
        let solve = |factor| {
            solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    fixture.resolved(),
                    fixture
                        .conditions_at(700.0, 101_325.0)
                        .expect("TP-1907 scaling temperature must be supported"),
                    scaled_initial_composition(&fixture, factor),
                )
                .with_phase_control_policy(PhaseControlPolicy::default())
                .with_solve_options(strict_extensivity_solve_options()),
            )
            .expect("strict TP-1907 scaling solve must accept")
        };
        let baseline = solve(1.0);
        let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);
        assert_bounded_solution_accepted(&baseline, contract);
        let graphite = PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()));
        let baseline_status = baseline
            .phase_status(&graphite)
            .expect("graphite phase must be represented");

        for factor in [1.0e-3, 1.0e3] {
            let candidate = solve(factor);
            assert_bounded_solution_accepted(&candidate, contract);
            assert_eq!(candidate.phase_status(&graphite), Some(baseline_status));
            for (index, (&expected, &actual)) in baseline
                .component_moles()
                .iter()
                .zip(candidate.component_moles())
                .enumerate()
            {
                let recovered = actual / factor;
                let relative =
                    (expected - recovered).abs() / expected.abs().max(recovered.abs()).max(1.0e-12);
                assert!(
                    relative <= 1.0e-7,
                    "TP-1907 extensive P,T mismatch at component {index}, factor={factor:e}: baseline={expected:e}, recovered={recovered:e}, relative={relative:e}"
                );
            }
        }
    }

    #[test]
    #[ignore = "release characterization of real TP-1907 scale and threshold bands"]
    fn i5_tp1907_extensive_pt_scale_matrix_characterization() {
        const TRACE_FLOOR: f64 = 1.0e-30;
        // This is the documented default PhaseControlPolicy threshold. It is
        // printed here as diagnostic context, not exposed or modified by the
        // test. A small phase near it would be a threshold interaction, not a
        // failure of the extensive thermodynamic invariant by itself.
        const DEFAULT_PHASE_EPS: f64 = 1.0e-30;

        let fixture = fixture();
        let evidence_request = PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())),
            PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
        );
        let solve = |temperature_k, factor| {
            solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    fixture.resolved(),
                    fixture
                        .conditions_at(temperature_k, 101_325.0)
                        .expect("TP-1907 scale-matrix temperature must be supported"),
                    scaled_initial_composition(&fixture, factor),
                )
                .with_phase_control_policy(PhaseControlPolicy::default())
                .with_solve_options(
                    strict_extensivity_solve_options()
                        .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute {
                            floor: TRACE_FLOOR,
                        })
                        .with_diagnostics(EquilibriumDiagnosticsOptions::enabled(
                            EquilibriumDiagnosticsMode::Summary,
                        )),
                ),
            )
        };

        println!("NASA TP-1907 CHON+graphite P,T extensive-scale characterization");
        println!(
            "trace_floor={TRACE_FLOOR:e} mol phase_eps={DEFAULT_PHASE_EPS:e} mol; fresh solves, no continuation"
        );
        println!(
            "  T K    factor     topology              min-active mol   graphite TPD  max n/f err  max gas x err  residual    balance     route"
        );

        for temperature_k in [680.0, 700.0, 720.0] {
            let baseline = solve(temperature_k, 1.0)
                .expect("unit-scale TP-1907 characterization solve must accept");
            let baseline_topology = active_phase_topology(&baseline);
            let baseline_tpd = graphite_tpd_from_solution(&baseline, &evidence_request);
            let baseline_component_moles = component_moles_by_identity(&baseline);
            let baseline_gas_x = gas_mole_fractions_by_identity(&baseline);
            let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);
            assert_bounded_solution_accepted(&baseline, contract);

            for factor in [1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4] {
                let candidate = solve(temperature_k, factor).unwrap_or_else(|error| {
                    panic!(
                        "production extensive recovery failed at {temperature_k} K, factor={factor:e}: {error:?}"
                    )
                });
                assert_bounded_solution_accepted(&candidate, contract);
                let amount_error = max_recovered_relative_error(
                    &baseline_component_moles,
                    &component_moles_by_identity(&candidate),
                    factor,
                );
                let gas_error = max_recovered_relative_error(
                    &baseline_gas_x,
                    &gas_mole_fractions_by_identity(&candidate),
                    1.0,
                );
                let tpd = graphite_tpd_from_solution(&candidate, &evidence_request);
                let validation = candidate.accepted_solution().validation();
                let recovery = candidate.extensive_normalization_recovery();
                let recovery_label = recovery.map_or("direct", |evidence| {
                    if evidence.reconstructed_physical_boundary {
                        "reconstructed"
                    } else {
                        "physical-retry"
                    }
                });
                println!(
                    "  {temperature_k:6.1}  {factor:8.1e}  {topology:20}  {min_active:14.6e}  {tpd:12.4e}  {amount_error:11.3e}  {gas_error:13.3e}  {residual:10.3e}  {balance:10.3e}  {recovery_label}",
                    topology = active_phase_topology(&candidate),
                    min_active = smallest_active_phase_total(&candidate),
                    residual = validation.residual_l2_norm,
                    balance = validation.max_abs_element_balance_error,
                );
                assert_eq!(
                    active_phase_topology(&candidate),
                    baseline_topology,
                    "scale-matrix topology drift must be explicit in the characterization"
                );
                assert_eq!(
                    tpd.signum(),
                    baseline_tpd.signum(),
                    "scale-matrix graphite TPD-sign drift must be explicit in the characterization"
                );
                assert!(
                    amount_error <= 1.0e-4,
                    "physical component amounts lost extensive covariance at {temperature_k} K, factor={factor:e}: {amount_error:e}"
                );
                assert!(
                    gas_error <= 1.0e-8,
                    "gas composition lost intensive invariance at {temperature_k} K, factor={factor:e}: {gas_error:e}"
                );
                for (element_index, (&baseline_total, &candidate_total)) in baseline
                    .build_report()
                    .element_totals()
                    .iter()
                    .zip(candidate.build_report().element_totals())
                    .enumerate()
                {
                    let recovered = candidate_total / factor;
                    let relative = (baseline_total - recovered).abs()
                        / baseline_total.abs().max(recovered.abs()).max(1.0e-30);
                    assert!(
                        relative <= 1.0e-12,
                        "physical build-report element total {element_index} was not reconstructed at {temperature_k} K, factor={factor:e}: baseline={baseline_total:e}, candidate={candidate_total:e}"
                    );
                }
                if factor == 1.0e4 {
                    let evidence = recovery.expect(
                        "the known large-inventory basin must retain explicit recovery provenance",
                    );
                    let expected_scale: f64 = scaled_initial_composition(&fixture, factor)
                        .moles()
                        .iter()
                        .filter(|moles| **moles > 0.0)
                        .sum();
                    assert!(
                        (evidence.physical_inventory_scale - expected_scale).abs()
                            <= 1.0e-12 * expected_scale.max(1.0),
                        "recovery evidence must report the physical input inventory scale"
                    );
                    if temperature_k < 720.0 {
                        assert!(evidence.reconstructed_physical_boundary);
                        assert!(evidence.physical_retry_backend.is_none());
                        assert!(evidence.physical_retry_failure.is_some());
                    } else {
                        assert!(!evidence.reconstructed_physical_boundary);
                        assert!(evidence.physical_retry_backend.is_some());
                        assert!(evidence.physical_retry_failure.is_none());
                    }
                    assert!(
                        candidate
                            .diagnostics_report()
                            .expect("enabled diagnostics must retain recovery evidence")
                            .events()
                            .iter()
                            .any(|event| matches!(
                                event,
                                EquilibriumDiagnosticEvent::ExtensiveNormalizationRecoveryAccepted {
                                    reconstructed_physical_boundary,
                                    ..
                                } if *reconstructed_physical_boundary
                                    == evidence.reconstructed_physical_boundary
                            ))
                    );
                } else {
                    assert!(
                        recovery.is_none(),
                        "ordinary scale factor {factor:e} unexpectedly entered recovery"
                    );
                }
            }
        }

        let disabled_error = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                fixture
                    .conditions_at(700.0, 101_325.0)
                    .expect("TP-1907 disabled-recovery temperature must be supported"),
                scaled_initial_composition(&fixture, 1.0e4),
            )
            .with_phase_control_policy(PhaseControlPolicy::default())
            .with_solve_options(
                strict_extensivity_solve_options()
                    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: TRACE_FLOOR })
                    .with_extensive_normalization_policy(ExtensiveNormalizationPolicy::Disabled),
            ),
        )
        .expect_err("disabling recovery must preserve the known physical-coordinate failure");
        assert_eq!(
            disabled_error.kind(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::
                ReactionExtentErrorKind::AllBackendsFailed
        );
    }

    #[test]
    #[ignore = "release characterization of prepared T-range normalization recovery"]
    fn i5_tp1907_extensive_temperature_range_recovers_large_inventory_transactionally() {
        const FACTOR: f64 = 1.0e4;
        const TEMPERATURES: [f64; 3] = [680.0, 700.0, 720.0];
        let fixture = fixture();
        let progress = Arc::new(Mutex::new(Vec::new()));
        let progress_sink = Arc::clone(&progress);
        let execution_control =
            EquilibriumExecutionControl::new().with_progress_sink(move |event| {
                progress_sink
                    .lock()
                    .expect("progress lock must remain available")
                    .push(event);
            });
        let solve_range = |factor, normalization_policy| {
            TemperatureRangeRequest::new(
                fixture.resolved(),
                scaled_initial_composition(&fixture, factor),
                101_325.0,
                101_325.0,
                TemperatureGrid::new(TEMPERATURES.to_vec()).expect("grid must be monotone"),
            )
            .expect("scaled range request must validate")
            .with_phase_control_policy(PhaseControlPolicy::default())
            .with_solve_options(
                strict_extensivity_solve_options()
                    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1.0e-30 })
                    .with_extensive_normalization_policy(normalization_policy)
                    .with_execution_control(execution_control.clone()),
            )
            .solve()
        };

        let baseline = solve_range(1.0, ExtensiveNormalizationPolicy::OnNumericalFailure)
            .expect("unit-scale prepared range must solve");
        progress
            .lock()
            .expect("progress lock must remain available")
            .clear();
        let scaled = solve_range(FACTOR, ExtensiveNormalizationPolicy::OnNumericalFailure)
            .expect("large-inventory prepared range must recover transactionally");
        assert_eq!(baseline.points().len(), TEMPERATURES.len());
        assert_eq!(scaled.points().len(), TEMPERATURES.len());

        let mut recovery_points = 0usize;
        println!("NASA TP-1907 prepared T-range extensive recovery");
        println!(
            "  T K    preparation           topology              max n/f err  graphite TPD  route"
        );
        for ((&temperature_k, baseline_point), scaled_point) in TEMPERATURES
            .iter()
            .zip(baseline.points())
            .zip(scaled.points())
        {
            let baseline_solution = baseline_point.solution();
            let scaled_solution = scaled_point.solution();
            assert_bounded_solution_accepted(
                scaled_solution,
                AcceptedSolutionContract::new(1.0e-6, 1.0e-6),
            );
            let amount_error = max_recovered_relative_error(
                &component_moles_by_identity(baseline_solution),
                &component_moles_by_identity(scaled_solution),
                FACTOR,
            );
            let baseline_tpd = graphite_final_tpd(baseline_solution);
            let scaled_tpd = graphite_final_tpd(scaled_solution);
            assert_eq!(
                active_phase_topology(scaled_solution),
                active_phase_topology(baseline_solution),
                "prepared range topology changed under extensive scaling at {temperature_k} K"
            );
            assert_eq!(
                scaled_tpd.signum(),
                baseline_tpd.signum(),
                "prepared range graphite stability changed under extensive scaling at {temperature_k} K"
            );
            assert!(
                amount_error <= 1.0e-4,
                "prepared range lost extensive covariance at {temperature_k} K: {amount_error:e}"
            );
            for (element_index, (&baseline_total, &scaled_total)) in baseline_solution
                .build_report()
                .element_totals()
                .iter()
                .zip(scaled_solution.build_report().element_totals())
                .enumerate()
            {
                let recovered = scaled_total / FACTOR;
                let relative = (baseline_total - recovered).abs()
                    / baseline_total.abs().max(recovered.abs()).max(1.0e-30);
                assert!(
                    relative <= 1.0e-12,
                    "range build-report element {element_index} remained normalized at {temperature_k} K"
                );
            }

            let recovery = scaled_solution.extensive_normalization_recovery();
            if scaled_point.report().preparation()
                == TemperatureRangePointPreparation::RecoveryFormulation
            {
                recovery_points += 1;
                let evidence = recovery.expect("recovery preparation must retain typed evidence");
                assert!(
                    (evidence.physical_inventory_scale
                        - scaled_initial_composition(&fixture, FACTOR)
                            .moles()
                            .iter()
                            .filter(|moles| **moles > 0.0)
                            .sum::<f64>())
                    .abs()
                        <= 1.0e-8 * evidence.physical_inventory_scale
                );
            } else {
                assert!(
                    recovery.is_none(),
                    "a direct prepared point must not carry stale recovery evidence"
                );
            }
            println!(
                "  {temperature_k:6.1}  {preparation:20?}  {topology:20}  {amount_error:11.3e}  {scaled_tpd:12.4e}  {route}",
                preparation = scaled_point.report().preparation(),
                topology = active_phase_topology(scaled_solution),
                route = if recovery.is_some() {
                    "recovery"
                } else {
                    "direct"
                },
            );
        }
        assert!(
            recovery_points > 0,
            "known 10^4 range must exercise recovery"
        );
        assert_eq!(scaled.report().formulation_builds(), 1 + recovery_points);
        assert_eq!(
            scaled.report().formulation_reuses(),
            scaled
                .points()
                .iter()
                .filter(|point| {
                    point.report().preparation()
                        == TemperatureRangePointPreparation::ReusedFormulation
                })
                .count()
        );
        let scaled_progress = progress
            .lock()
            .expect("progress lock must remain available");
        assert_eq!(
            scaled_progress
                .iter()
                .filter(|event| event.stage() == EquilibriumProgressStage::PointStarted)
                .count(),
            TEMPERATURES.len(),
            "internal point recovery must not duplicate range progress"
        );
        assert_eq!(
            scaled_progress
                .iter()
                .filter(|event| event.stage() == EquilibriumProgressStage::PointAccepted)
                .count(),
            TEMPERATURES.len(),
            "each transactionally published point must emit exactly one acceptance"
        );
        drop(scaled_progress);
        progress
            .lock()
            .expect("progress lock must remain available")
            .clear();

        let disabled = solve_range(FACTOR, ExtensiveNormalizationPolicy::Disabled)
            .expect_err("disabled recovery must preserve the known prepared-range failure");
        let ReactionExtentError::TemperatureRangePointFailed {
            point_index, cause, ..
        } = disabled
        else {
            panic!("range failure must preserve its point coordinate");
        };
        assert_eq!(point_index, 0);
        assert_eq!(
            cause.kind(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::
                ReactionExtentErrorKind::AllBackendsFailed
        );
        let failed_progress = progress
            .lock()
            .expect("progress lock must remain available");
        assert_eq!(
            failed_progress
                .iter()
                .filter(|event| event.stage() == EquilibriumProgressStage::PointStarted)
                .count(),
            1
        );
        assert_eq!(
            failed_progress
                .iter()
                .filter(|event| event.stage() == EquilibriumProgressStage::PointAccepted)
                .count(),
            0,
            "failed range must not publish an accepted point"
        );
    }

    #[test]
    fn tp1907_default_and_legacy_nr_keep_the_same_bounded_pt_state() {
        let fixture = fixture();
        let solve = |options| {
            solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    fixture.resolved(),
                    fixture
                        .conditions_at(700.0, 101_325.0)
                        .expect("TP-1907 backend-matrix temperature must be supported"),
                    fixture
                        .initial_composition()
                        .expect("TP-1907 backend-matrix inventory must validate"),
                )
                .with_phase_control_policy(PhaseControlPolicy::default())
                .with_solve_options(options),
            )
            .expect("TP-1907 bounded backend-matrix solve must accept")
        };
        let default = solve(EquilibriumSolveOptions::default());
        let legacy_nr = solve(EquilibriumSolveOptions::default().with_solver_backend(Solvers::NR));
        let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);
        assert_bounded_solution_accepted(&default, contract);
        assert_bounded_solution_accepted(&legacy_nr, contract);

        for phase in default.phases() {
            assert_eq!(
                default.phase_status(phase.id()),
                legacy_nr.phase_status(phase.id()),
                "backend selection must not change phase {:?} topology",
                phase.id(),
            );
        }
        for (index, (&default_moles, &legacy_moles)) in default
            .component_moles()
            .iter()
            .zip(legacy_nr.component_moles())
            .enumerate()
        {
            let relative = (default_moles - legacy_moles).abs()
                / default_moles.abs().max(legacy_moles.abs()).max(1.0e-12);
            assert!(
                relative <= 1.0e-5,
                "P,T backend matrix differs at component {index}: default={default_moles:e}, legacy_nr={legacy_moles:e}, relative={relative:e}"
            );
        }
    }

    #[test]
    fn source_faithful_single_point_ranges_match_independent_pt_phase_control() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let fixture = fixture();
        for temperature_k in [680.0, 700.0, 720.0, 740.0] {
            let direct = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    fixture.resolved(),
                    fixture
                        .conditions_at(temperature_k, 101_325.0)
                        .expect("fixture temperature must be supported"),
                    fixture
                        .initial_composition()
                        .expect("fixture inventory must validate"),
                )
                .with_phase_control_policy(PhaseControlPolicy::default()),
            )
            .expect("independent source-faithful P,T solve must succeed");
            let range = solve_graphite_temperature_range(&fixture, vec![temperature_k]);
            let ranged = range
                .points()
                .first()
                .expect("one-point range must publish its only point")
                .solution();

            assert_eq!(
                direct.phase_status(&PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()))),
                ranged.phase_status(&PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()))),
                "the first range point must use the same graphite lifecycle result as an independent P,T solve at {temperature_k} K"
            );
            for (index, (&direct_moles, &range_moles)) in direct
                .component_moles()
                .iter()
                .zip(ranged.component_moles())
                .enumerate()
            {
                assert!(
                    (direct_moles - range_moles).abs() <= 1e-10,
                    "component {index} changed between direct P,T and a one-point range at {temperature_k} K: direct={direct_moles:e}, range={range_moles:e}"
                );
            }
            assert!(
                (system_mole_fraction(&direct, &graphite_component())
                    - system_mole_fraction(ranged, &graphite_component()))
                .abs()
                    <= 1e-12,
                "one-point range graphite system fraction must equal direct P,T at {temperature_k} K"
            );
            assert!(
                !range.points()[0].report().used_continuation_seed(),
                "the first range point must not claim physical continuation at {temperature_k} K"
            );
        }
        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    #[test]
    fn i5_nasa_tp1907_rows_preserve_multiphase_system_normalization() {
        let before = frozen_snapshot();
        let dataset = dataset();
        assert_eq!(
            dataset.metadata().evidence_kind,
            FrozenReferenceEvidenceKind::FrozenExternal
        );
        assert_eq!(dataset.rows().len(), 4);
        assert_eq!(
            dataset
                .rows()
                .iter()
                .map(|row| row.temperature_k)
                .collect::<Vec<_>>(),
            [680.0, 700.0, 720.0, 740.0]
        );
        for row in dataset.rows() {
            assert_eq!(
                row.normalization,
                FrozenMultiphaseNormalization::SystemTotalMoleFraction
            );
            let total = row
                .gas_species
                .iter()
                .map(|entry| entry.system_mole_fraction)
                .sum::<f64>()
                + row
                    .condensed_species
                    .iter()
                    .map(|entry| entry.system_mole_fraction)
                    .sum::<f64>();
            assert!(
                (total - 1.0).abs() <= 2.0e-5,
                "{} K source total={total}",
                row.temperature_k
            );
        }
        let graphite = |temperature| {
            dataset
                .rows()
                .iter()
                .find(|row| row.temperature_k == temperature)
                .and_then(|row| {
                    row.condensed_species
                        .iter()
                        .find(|entry| entry.phase == "C(gr)")
                })
                .map(|entry| entry.system_mole_fraction)
                .expect("graphite row must exist")
        };
        assert!(graphite(700.0) > 0.0);
        assert_eq!(graphite(720.0), 0.0);
        assert_eq!(before, frozen_snapshot());
    }

    #[test]
    fn i5_nasa_tp1907_source_faithful_feed_preserves_chon_ar_and_published_scalars() {
        let feed = reconstructed_source_feed()
            .expect("TP-1906/1907 source formula must yield a valid feed");
        let totals = verify_element_equivalent_feed()
            .expect("ordinary molecule basis must preserve all source elements");
        assert_eq!(feed.h_to_c, 2.0);
        assert_eq!(feed.equivalence_ratio, 1.25);
        assert!((feed.alpha - 1.2).abs() <= 1.0e-14);
        assert!((feed.o2_moles - 1.2).abs() <= 1.0e-14);
        assert!((feed.n2_moles - 4.473_104_4).abs() <= 1.0e-14);
        assert!((feed.ar_moles - 0.053_648_16).abs() <= 1.0e-14);
        assert!((feed.co2_moles - 0.001_827_36).abs() <= 1.0e-14);
        assert!((totals["C"] - (1.0 + feed.co2_moles)).abs() <= 1.0e-14);
        assert_eq!(totals.get("H"), Some(&2.0));
        assert!((totals["O"] - 2.0 * (feed.o2_moles + feed.co2_moles)).abs() <= 1.0e-14);
        assert!((totals["N"] - 2.0 * feed.n2_moles).abs() <= 1.0e-14);
        assert!((totals["Ar"] - feed.ar_moles).abs() <= 1.0e-14);
        validate_published_feed_diagnostics(feed)
            .expect("published F/A and chemical-ER diagnostics must remain consistent");
        assert!((feed.fuel_air_mass_ratio() - 0.084_535_289_081_489_9).abs() <= 1.0e-15);
        assert!((feed.chemical_equivalence_ratio() - 1.25).abs() <= 1.0e-14);
    }

    #[test]
    fn i5_nasa_tp1907_exact_chon_ar_graphite_universe_resolves_offline() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let dataset = dataset();
        let fixture = fixture();
        assert_eq!(
            fixture.resolved().layout().component_count(),
            TP1907_LOCAL_GAS_SPECIES.len() + 1
        );
        assert_eq!(fixture.structure().element_count, 5);
        assert_eq!(fixture.structure().element_rank, 5);
        assert_eq!(
            fixture.structure().reaction_dimension,
            TP1907_LOCAL_GAS_SPECIES.len() - 4
        );
        assert_eq!(
            fixture.preflight().len(),
            TP1907_LOCAL_GAS_SPECIES.len() + 1
        );
        assert!(
            fixture
                .preflight()
                .iter()
                .all(|row| !row.record_key.is_empty())
        );
        assert!(
            fixture
                .preflight()
                .iter()
                .any(|row| row.component == graphite_component() && row.library == "NASA_cond")
        );
        for row in dataset.rows() {
            fixture
                .conditions(row)
                .expect("all selected source temperatures must lie in the common local domain");
        }
        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    #[test]
    #[ignore = "release I5 characterization: source-faithful NASA TP-1907 graphite sweep"]
    fn i5_nasa_tp1907_chon_graphite_source_audit_and_sweep() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let dataset = dataset();
        let fixture = fixture();
        let evidence_request = PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())),
            PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
        );
        let feed = reconstructed_source_feed().expect("source-faithful feed must be available");
        println!("NASA TP-1906/1907 source reconstruction");
        println!(
            "  H/C={:.6} ER={:.6} alpha={:.8}",
            feed.h_to_c, feed.equivalence_ratio, feed.alpha
        );
        println!(
            "  dry air / O2: N2={:.7} Ar={:.7} CO2={:.7}",
            feed.n2_moles / feed.o2_moles,
            feed.ar_moles / feed.o2_moles,
            feed.co2_moles / feed.o2_moles,
        );
        println!(
            "  diagnostics: F/A={:.9} (published 0.084535), chemical ER={:.7} (published 1.2496)",
            feed.fuel_air_mass_ratio(),
            feed.chemical_equivalence_ratio(),
        );
        println!("NASA TP-1907 Table 11.3E | source-faithful CH2 + dry-air P=101325 Pa");
        println!(
            "T K    NASA C(gr)    KiThe C(gr)    NASA topology  KiThe topology  TPD              transitions"
        );
        for reference in dataset.rows() {
            let solution = solve_graphite_pt(&fixture, reference.temperature_k);
            let external_graphite = reference
                .condensed_species
                .iter()
                .find(|entry| entry.phase == "C(gr)")
                .expect("frozen graphite row must exist")
                .system_mole_fraction;
            let status = solution
                .phase_status(&evidence_request.candidate_phase)
                .expect("graphite phase must be represented");
            let evidence = match status {
                PhaseStatus::Inactive => {
                    stable_inactive_evidence_from_solution(&solution, &evidence_request).unwrap()
                }
                PhaseStatus::Active | PhaseStatus::Appeared => {
                    activation_evidence_from_solution(&solution, &evidence_request).unwrap()
                }
                other => panic!(
                    "{} K graphite has unsupported final status {other:?}",
                    reference.temperature_k
                ),
            };
            let graphite_moles = solution
                .moles_for(&graphite_component())
                .expect("graphite component must be represented");
            let graphite_system_fraction =
                graphite_moles / solution.component_moles().iter().sum::<f64>();
            println!(
                "{:<6.1} {:<15.6e} {:<15.6e} {:<15} {:<15?} {:<16.6e} {}",
                reference.temperature_k,
                external_graphite,
                graphite_system_fraction,
                if external_graphite > 0.0 {
                    "active"
                } else {
                    "inactive"
                },
                status,
                evidence.boundary_minimum_tpd.unwrap_or(f64::NAN),
                solution.phase_control_transitions()
            );
            println!(
                "  accepted backend={:?} nonlinear_iterations={}",
                solution.solve_report().accepted_backend,
                solution.nonlinear_iterations(),
            );
            if external_graphite > 0.0 {
                assert!(matches!(
                    status,
                    PhaseStatus::Active | PhaseStatus::Appeared
                ));
                assert!(graphite_moles > 0.0);
                assert!(
                    evidence
                        .boundary_minimum_tpd
                        .expect("activation evidence needs TPD")
                        < 0.0
                );
            } else {
                assert_eq!(status, PhaseStatus::Inactive);
                assert_eq!(graphite_moles, 0.0);
                assert!(
                    evidence
                        .boundary_minimum_tpd
                        .expect("inactive evidence needs TPD")
                        > 0.0
                );
            }
            let validation = solution.accepted_solution().validation();
            assert!(validation.residual_l2_norm.is_finite());
            assert!(validation.max_abs_element_balance_error.is_finite());
            assert_tp1907_external_regression_envelope(&fixture, reference, &solution);
            println!(
                "  component     NASA(system)   KiThe(system)  relative       dlog10       note"
            );
            for row in fixture
                .compare_system_composition(reference, &solution)
                .expect("accepted state must align to frozen source identities")
            {
                match row.kithe_system_fraction {
                    Some(local) => println!(
                        "  {:<12} {:<14.6e} {:<14.6e} {:+.3e}   {:+.3e}   -",
                        row.identity,
                        row.source_system_fraction,
                        local,
                        row.relative_error.unwrap_or(f64::NAN),
                        row.delta_log10.unwrap_or(f64::NAN),
                    ),
                    None => println!(
                        "  {:<12} {:<14.6e} excluded       -            -            {}",
                        row.identity,
                        row.source_system_fraction,
                        row.exclusion_reason.as_deref().unwrap_or("-")
                    ),
                }
            }
        }
        println!("accepted-state continuation");
        for (direction, range) in [
            (
                "forward",
                solve_graphite_temperature_range(&fixture, vec![680.0, 700.0, 720.0, 740.0]),
            ),
            (
                "reverse",
                solve_graphite_temperature_range(&fixture, vec![740.0, 720.0, 700.0, 680.0]),
            ),
        ] {
            println!("  direction={direction}");
            for (index, point) in range.points().iter().enumerate() {
                let status = point
                    .solution()
                    .phase_status(&evidence_request.candidate_phase)
                    .expect("continued graphite phase must be represented");
                let graphite = system_mole_fraction(point.solution(), &graphite_component());
                println!(
                    "    T={:6.1} K status={status:?} C(gr)/system={graphite:.6e} continuation={} transitions={} backend={:?} iterations={}",
                    point.report().temperature(),
                    point.report().used_continuation_seed(),
                    point.report().phase_control_transitions(),
                    point.solution().solve_report().accepted_backend,
                    point.solution().nonlinear_iterations(),
                );
                if index == 0 {
                    assert!(!point.report().used_continuation_seed());
                } else {
                    assert!(point.report().used_continuation_seed());
                }
                let expects_graphite = if direction == "forward" {
                    point.report().temperature() <= 700.0
                } else {
                    point.report().temperature() <= 700.0
                };
                assert_eq!(
                    matches!(status, PhaseStatus::Active | PhaseStatus::Appeared),
                    expects_graphite,
                    "{direction} continuation returned wrong graphite topology at {} K",
                    point.report().temperature()
                );
            }
        }
        let (boundary_temperature_k, boundary_tpd, boundary_iterations) =
            bisect_graphite_tpd_boundary(&fixture, &evidence_request);
        println!(
            "canonical graphite TPD boundary: T={boundary_temperature_k:.6} K tpd={boundary_tpd:.3e} iterations={boundary_iterations} NASA bracket=[700, 720] K"
        );
        assert!((700.0..=720.0).contains(&boundary_temperature_k));
        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }
}
