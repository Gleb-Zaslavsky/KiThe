//! Offline I5 evidence for NASA TP-1907 Table 11.3E CHON + graphite.

#[cfg(test)]
mod tests {
    #[path = "../nasa_tp1907_frontend_diagnostics.rs"]
    mod frontend_diagnostics;

    use std::collections::BTreeMap;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::{Arc, Mutex};

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::{
        EquilibriumCalculator, EquilibriumCalculatorOutcome,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
        EquilibriumDiagnosticEvent, EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_element_inventory::{
        ElementInventory, FormalElementCarrier,
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
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
        RustedSciTheSolveContract, RustedSciTheSolver, prepare_rst_symbolic_problem_from_prepared,
        prepare_rst_symbolic_problem_from_prepared_with_frontend, with_test_rst_frontend,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        SolverAttemptOutcome, SolverBackend, SolverPolicy,
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
        TP1907_LOCAL_GAS_SPECIES, TP1907_REFERENCE_PRESSURE_PA, graphite_component,
        load_nasa_tp1907_chon_graphite_dataset, reconstructed_source_feed,
        validate_published_feed_diagnostics, verify_element_equivalent_feed,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
        PurePhaseProductionEvidenceRequest, activation_evidence_from_solution,
        stable_inactive_evidence_from_solution,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
        PhaseEquilibriumBuildRequest, PhaseEquilibriumInputKind, build_phase_equilibrium_problem,
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
    use crate::Thermodynamics::User_PhaseOrSolution::element_composition_and_molar_mass;
    use crate::Thermodynamics::phase_layout::PhaseId;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;
    use RustedSciThe::numerical::Nonlinear_systems::prelude::{
        JacobianProvider, LambdifyExecutionPolicy, NonlinearProblem, SymbolicLambdifyFrontend,
        SymbolicNonlinearProblem, SymbolicProblemOptions,
    };
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use nalgebra::DVector;
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use std::time::Instant;

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
                // This is a seed-only basin hint for the fixed-active
                // diagnostic, not a physical phase threshold. The physical
                // inventory and all frozen thermochemistry remain unchanged.
                TraceSpeciesSeedPolicy::Absolute { floor: 1.0e-20 },
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
                // Keep the normalized diagnostic seed scale-free. The
                // minimum floor is only a lower bound for trace species.
                normalization
                    .normalize_trace_seed_policy(
                        TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                            fraction: 1.0e-12,
                            minimum_floor: 1.0e-20,
                        },
                    )
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

    fn l2_norm(values: &[f64]) -> f64 {
        values.iter().map(|value| value * value).sum::<f64>().sqrt()
    }

    fn stable_equation_hash(equations: &[String]) -> u64 {
        // FNV-1a with length separators keeps the fixture hash independent of
        // Rust's randomized hashers and unambiguous for concatenated strings.
        let mut hash = 0xcbf29ce484222325_u64;
        for equation in equations {
            for byte in (equation.len() as u64).to_le_bytes() {
                hash ^= u64::from(byte);
                hash = hash.wrapping_mul(0x100000001b3_u64);
            }
            for byte in equation.as_bytes() {
                hash ^= u64::from(*byte);
                hash = hash.wrapping_mul(0x100000001b3_u64);
            }
        }
        hash
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
        let mut runner = PreparedEquilibriumRunner::new(prepared)
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
        let mut runner = PreparedPhaseControlRunner::new(problem, false)?;
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
        let mut runner = PreparedEquilibriumRunner::new(prepared)
            .expect("prepared TP-1907 runner must construct");
        let settings = runner.configure();
        settings.solver = Solvers::LM;
        // The metamorphic case deliberately changes reaction-coordinate
        // norms. Apply the canonical row scaling so the reference solver
        // measures the same dimensionless formulation in every basis.
        settings.scaling_flag = true;
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

    fn frozen_element_inventory(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        initial: &MultiphaseInitialComposition,
    ) -> (ElementInventory, Vec<String>, Vec<f64>) {
        let layout = MultiphaseEquilibriumLayout::new(fixture.resolved().phase_specs().to_vec())
            .expect("TP-1907 calculator layout must remain valid");
        let (element_matrix, _, element_labels) = element_composition_and_molar_mass(
            fixture.resolved().phase_data(),
            fixture.resolved().layout(),
            None,
        )
        .expect("TP-1907 element matrix must remain available");
        let totals = initial
            .element_totals(&layout, &element_matrix)
            .expect("TP-1907 reference b must be derivable as A^T n_initial");
        let inventory = ElementInventory::from_amounts(
            element_labels.iter().cloned().zip(totals.iter().copied()),
        )
        .expect("TP-1907 reference elemental inventory must validate");
        (inventory, element_labels, totals)
    }

    fn calculator_tp1907_builder(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        repository: Arc<ThermoRepository>,
    ) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorBuilder
    {
        EquilibriumCalculator::from_phases(fixture.resolved().phase_specs().iter().cloned())
            .with_repository(repository)
            .prefer_libraries(["NASA_gas", "NASA_cond"])
            .offline_only()
            .pressure_pa(TP1907_REFERENCE_PRESSURE_PA)
            .reference_pressure_pa(TP1907_REFERENCE_PRESSURE_PA)
            .phase_control(PhaseControlPolicy::default())
            .production_cascade()
    }

    fn calculator_tp1907_pt(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        repository: Arc<ThermoRepository>,
        temperature_k: f64,
        initial: Option<Vec<f64>>,
        inventory: Option<ElementInventory>,
    ) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorPoint
    {
        let builder = calculator_tp1907_builder(fixture, repository).at_temperature(temperature_k);
        let outcome = match (initial, inventory) {
            (Some(moles), None) => builder.initial_moles(moles).solve(),
            (None, Some(inventory)) => builder.element_inventory(inventory).solve(),
            _ => panic!("TP-1907 calculator test requires exactly one physical input"),
        }
        .expect("TP-1907 calculator P,T route must solve");
        match outcome {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected TP-1907 calculator P,T point, got {other:?}"),
        }
    }

    fn calculator_tp1907_pt_with_element_seed(
        fixture: &ResolvedNasaTp1907ChonGraphiteFixture,
        repository: Arc<ThermoRepository>,
        temperature_k: f64,
        inventory: ElementInventory,
        seed: MultiphaseInitialComposition,
    ) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorPoint
    {
        let outcome = calculator_tp1907_builder(fixture, repository)
            .at_temperature(temperature_k)
            .element_inventory(inventory)
            .element_numerical_seed(seed)
            .solve()
            .expect("TP-1907 seeded elemental calculator P,T route must solve");
        match outcome {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected seeded TP-1907 calculator P,T point, got {other:?}"),
        }
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
    fn tp1907_frozen_pt_molecular_elemental_and_formal_routes_match() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let fixture = fixture();
        let frozen_dataset = dataset();
        let reference = frozen_dataset
            .rows()
            .iter()
            .find(|row| row.temperature_k == 700.0)
            .expect("TP-1907 frozen P,T row at 700 K must remain available");
        let initial = fixture
            .initial_composition()
            .expect("TP-1907 molecular source feed must validate");
        let (inventory, element_labels, reference_b) = frozen_element_inventory(&fixture, &initial);
        let formal_carriers = element_labels
            .iter()
            .zip(reference_b.iter().copied())
            .map(|(element, amount)| {
                FormalElementCarrier::new(element.clone(), amount)
                    .expect("TP-1907 formal element carrier must validate")
            })
            .collect::<Vec<_>>();
        let formal_inventory = ElementInventory::from_formal_carriers(formal_carriers)
            .expect("TP-1907 formal carriers must reconstruct reference b");

        let molecular = calculator_tp1907_pt(
            &fixture,
            repository(),
            reference.temperature_k,
            Some(initial.moles().to_vec()),
            None,
        );
        let elemental = calculator_tp1907_pt(
            &fixture,
            repository(),
            reference.temperature_k,
            None,
            Some(inventory),
        );
        let formal = calculator_tp1907_pt(
            &fixture,
            repository(),
            reference.temperature_k,
            None,
            Some(formal_inventory),
        );

        let expected_components = molecular
            .solution()
            .metadata()
            .components()
            .iter()
            .map(|component| component.id().clone())
            .collect::<Vec<_>>();
        let expected_phases = molecular
            .solution()
            .metadata()
            .phases()
            .iter()
            .map(|phase| phase.id().clone())
            .collect::<Vec<_>>();
        let expected_topology = expected_phases
            .iter()
            .map(|phase| (phase.clone(), molecular.solution().phase_status(phase)))
            .collect::<Vec<_>>();

        for (route, point) in [("elemental", &elemental), ("formal", &formal)] {
            assert_eq!(
                point.resolved().phase_specs(),
                molecular.resolved().phase_specs(),
                "{route} route changed the explicit frozen universe"
            );
            let solution = point.solution();
            assert_eq!(
                solution
                    .metadata()
                    .components()
                    .iter()
                    .map(|component| component.id().clone())
                    .collect::<Vec<_>>(),
                expected_components,
                "{route} route changed component identity/order"
            );
            assert_eq!(
                solution
                    .metadata()
                    .phases()
                    .iter()
                    .map(|phase| phase.id().clone())
                    .collect::<Vec<_>>(),
                expected_phases,
                "{route} route changed phase topology layout"
            );
            assert_eq!(
                solution.build_report().input_kind(),
                PhaseEquilibriumInputKind::ElementInventory
            );
            assert_eq!(solution.build_report().element_labels(), &element_labels);
            for (index, (&expected, &actual)) in solution
                .build_report()
                .element_totals()
                .iter()
                .zip(reference_b.iter())
                .enumerate()
            {
                assert!(
                    (expected - actual).abs() <= 1.0e-12,
                    "{route} reference b mismatch at element {index}: actual={actual:e}, expected={expected:e}"
                );
            }
            for (phase, expected_status) in &expected_topology {
                assert_eq!(
                    solution.phase_status(phase),
                    *expected_status,
                    "{route} phase status changed for {phase:?}"
                );
            }
            for component in solution.metadata().components() {
                assert!(
                    !element_labels.iter().any(|element| {
                        component.id().substance == *element
                            && component.id().phase.as_option().is_none()
                    }),
                    "formal elemental carriers must not become bare thermodynamic components"
                );
            }
            for (index, (&expected, &actual)) in molecular
                .solution()
                .component_moles()
                .iter()
                .zip(solution.component_moles())
                .enumerate()
            {
                let relative =
                    (expected - actual).abs() / expected.abs().max(actual.abs()).max(1.0e-12);
                assert!(
                    relative <= 1.0e-6 || (expected - actual).abs() <= 1.0e-10,
                    "{route} component {index} mismatch: molecular={expected:e}, actual={actual:e}, relative={relative:e}"
                );
            }
            for component in molecular.solution().metadata().components() {
                let expected_fraction = molecular
                    .solution()
                    .mole_fraction_for(component.id())
                    .expect("molecular component fraction must exist");
                let actual_fraction = solution
                    .mole_fraction_for(component.id())
                    .expect("elemental component fraction must exist");
                assert!(
                    (expected_fraction - actual_fraction).abs() <= 1.0e-8,
                    "{route} fraction mismatch for {}: molecular={expected_fraction:e}, actual={actual_fraction:e}",
                    component.id().label()
                );
            }
            let expected_validation = molecular.solution().accepted_solution().validation();
            let actual_validation = solution.accepted_solution().validation();
            assert!(actual_validation.residual_l2_norm < 1.0e-5);
            assert!(
                (expected_validation.residual_l2_norm - actual_validation.residual_l2_norm).abs()
                    <= 1.0e-5,
                "{route} residual quality changed: molecular={:e}, actual={:e}",
                expected_validation.residual_l2_norm,
                actual_validation.residual_l2_norm
            );
            assert!(actual_validation.max_abs_element_balance_error <= 1.0e-6);

            let rows = fixture
                .compare_system_composition(reference, solution)
                .expect("{route} result must align with frozen TP-1907 identities");
            let max_external_discrepancy = rows
                .iter()
                .filter_map(|row| row.absolute_error.map(f64::abs))
                .fold(0.0_f64, f64::max);
            println!(
                "TP-1907 frozen P,T route={route} T={} K b={reference_b:?} max_external_absolute_discrepancy={max_external_discrepancy:e} residual={:e} balance={:e}",
                reference.temperature_k,
                actual_validation.residual_l2_norm,
                actual_validation.max_abs_element_balance_error,
            );
            assert_tp1907_external_regression_envelope(&fixture, reference, solution);
        }

        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    fn assert_tp1907_facade_inventory_scale(scale: f64, class: &str) {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let fixture = fixture();
        let frozen_dataset = dataset();
        let reference = frozen_dataset
            .rows()
            .iter()
            .find(|row| row.temperature_k == 700.0)
            .expect("TP-1907 frozen P,T row at 700 K must remain available");
        let initial = scaled_initial_composition(&fixture, scale);
        let (inventory, element_labels, reference_b) = frozen_element_inventory(&fixture, &initial);
        let formal_inventory = ElementInventory::from_formal_carriers(
            element_labels
                .iter()
                .zip(reference_b.iter().copied())
                .map(|(element, amount)| {
                    FormalElementCarrier::new(element.clone(), amount)
                        .expect("scaled TP-1907 formal carrier must validate")
                }),
        )
        .expect("scaled TP-1907 formal carriers must reconstruct b");

        let molecular = calculator_tp1907_pt(
            &fixture,
            repository(),
            reference.temperature_k,
            Some(initial.moles().to_vec()),
            None,
        );
        let elemental = calculator_tp1907_pt_with_element_seed(
            &fixture,
            repository(),
            reference.temperature_k,
            inventory,
            initial.clone(),
        );
        let formal = calculator_tp1907_pt_with_element_seed(
            &fixture,
            repository(),
            reference.temperature_k,
            formal_inventory,
            initial,
        );

        let expected_components = molecular
            .solution()
            .metadata()
            .components()
            .iter()
            .map(|component| component.id().clone())
            .collect::<Vec<_>>();
        let expected_phases = molecular
            .solution()
            .metadata()
            .phases()
            .iter()
            .map(|phase| phase.id().clone())
            .collect::<Vec<_>>();
        let expected_topology = expected_phases
            .iter()
            .map(|phase| (phase.clone(), molecular.solution().phase_status(phase)))
            .collect::<Vec<_>>();

        for (route, point) in [("elemental", &elemental), ("formal", &formal)] {
            let solution = point.solution();
            assert_eq!(
                point.resolved().phase_specs(),
                molecular.resolved().phase_specs(),
                "{class} {route} route changed the frozen universe"
            );
            assert_eq!(
                solution
                    .metadata()
                    .components()
                    .iter()
                    .map(|component| component.id().clone())
                    .collect::<Vec<_>>(),
                expected_components,
                "{class} {route} route changed component identity/order"
            );
            assert_eq!(
                solution
                    .metadata()
                    .phases()
                    .iter()
                    .map(|phase| phase.id().clone())
                    .collect::<Vec<_>>(),
                expected_phases,
                "{class} {route} route changed phase order"
            );
            assert_eq!(
                solution.build_report().input_kind(),
                PhaseEquilibriumInputKind::ElementInventory
            );
            assert_eq!(solution.build_report().element_labels(), &element_labels);
            for (index, (&expected, &actual)) in solution
                .build_report()
                .element_totals()
                .iter()
                .zip(reference_b.iter())
                .enumerate()
            {
                assert!(
                    (expected - actual).abs() <= 1.0e-10 * expected.abs().max(1.0),
                    "{class} {route} b mismatch at element {index}: actual={actual:e}, expected={expected:e}"
                );
            }
            for (phase, expected_status) in &expected_topology {
                assert_eq!(
                    solution.phase_status(phase),
                    *expected_status,
                    "{class} {route} phase status changed for {phase:?}"
                );
            }
            for (index, (&expected, &actual)) in molecular
                .solution()
                .component_moles()
                .iter()
                .zip(solution.component_moles())
                .enumerate()
            {
                let relative =
                    (expected - actual).abs() / expected.abs().max(actual.abs()).max(1.0e-12);
                assert!(
                    relative <= 1.0e-5 || (expected - actual).abs() <= 1.0e-10,
                    "{class} {route} component {index} mismatch: molecular={expected:e}, actual={actual:e}, relative={relative:e}"
                );
            }
            for component in molecular.solution().metadata().components() {
                let expected_fraction = molecular
                    .solution()
                    .mole_fraction_for(component.id())
                    .expect("molecular fraction must exist");
                let actual_fraction = solution
                    .mole_fraction_for(component.id())
                    .expect("elemental fraction must exist");
                assert!(
                    (expected_fraction - actual_fraction).abs() <= 1.0e-8,
                    "{class} {route} fraction mismatch for {}",
                    component.id().label()
                );
            }
            let validation = solution.accepted_solution().validation();
            println!(
                "TP-1907 frozen facade candidate scale={scale:e} class={class} route={route} backend={:?} attempts={} residual={:.3e} balance={:.3e}",
                solution.solve_report().accepted_backend,
                solution.solve_report().started_attempt_count(),
                validation.residual_l2_norm,
                validation.max_abs_element_balance_error,
            );
            assert!(
                validation.residual_l2_norm < 1.0e-5,
                "{class} {route} residual={:e}",
                validation.residual_l2_norm
            );
            assert!(
                validation.max_abs_element_balance_error <= 1.0e-6,
                "{class} {route} balance={:e}",
                validation.max_abs_element_balance_error
            );
            assert_tp1907_external_regression_envelope(&fixture, reference, solution);
            println!(
                "TP-1907 frozen facade scale={scale:e} class={class} route={route} backend={:?} attempts={} total_b={:.6e} residual={:.3e} balance={:.3e}",
                solution.solve_report().accepted_backend,
                solution.solve_report().started_attempt_count(),
                reference_b.iter().sum::<f64>(),
                validation.residual_l2_norm,
                validation.max_abs_element_balance_error,
            );
        }

        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    fn run_tp1907_upper_backend_acceptance_matrix() -> (usize, usize) {
        let fixture = fixture();
        let frozen_dataset = dataset();
        let reference = frozen_dataset
            .rows()
            .iter()
            .find(|row| row.temperature_k == 700.0)
            .expect("TP-1907 frozen P,T row at 700 K must remain available");
        let initial = scaled_initial_composition(&fixture, 1.0e8);
        let (inventory, _element_labels, reference_b) =
            frozen_element_inventory(&fixture, &initial);
        let (element_matrix, _, _) = element_composition_and_molar_mass(
            fixture.resolved().phase_data(),
            fixture.resolved().layout(),
            None,
        )
        .expect("TP-1907 backend matrix element composition must remain available");
        let oracle_options = EquilibriumSolveOptions::default()
            .with_max_iterations(500)
            .expect("TP-1907 backend matrix oracle iteration budget must validate")
            .with_tolerance(1.0e-12)
            .expect("TP-1907 backend matrix oracle tolerance must validate");
        let molecular = calculator_tp1907_builder(&fixture, repository())
            .solve_options(oracle_options)
            .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::TR)))
            .expect("TP-1907 backend matrix oracle policy must validate")
            .at_temperature(reference.temperature_k)
            .initial_moles(initial.moles().to_vec())
            .solve()
            .expect("TP-1907 Legacy(TR) molecular oracle must solve");
        let molecular = match molecular {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected TP-1907 molecular oracle point, got {other:?}"),
        };
        let molecular_solution = molecular.solution();
        let expected_topology = molecular_solution
            .metadata()
            .phases()
            .iter()
            .map(|phase| {
                (
                    phase.id().clone(),
                    molecular_solution.phase_status(phase.id()),
                )
            })
            .collect::<Vec<_>>();

        let backends = [
            (
                "rst-lm",
                "RST",
                SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
            ),
            (
                "rst-backtracking-lm",
                "RST",
                SolverBackend::RustedSciThe(RustedSciTheSolver::BacktrackingLevenbergMarquardt),
            ),
            (
                "rst-minpack-lm",
                "RST",
                SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
            ),
            (
                "rst-nielsen-lm",
                "RST",
                SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
            ),
            (
                "rst-trust-region-lm",
                "RST",
                SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
            ),
            (
                "rst-damped-newton",
                "RST",
                SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
            ),
            ("legacy-lm", "Legacy", SolverBackend::Legacy(Solvers::LM)),
            ("legacy-nr", "Legacy", SolverBackend::Legacy(Solvers::NR)),
            ("legacy-tr", "Legacy", SolverBackend::Legacy(Solvers::TR)),
        ];

        println!("TP-1907 upper inventory facade backend acceptance matrix");
        println!(
            "| backend | family | status | scaled residual | raw residual | max physical balance | max component rel | max fraction error | topology |"
        );
        println!("|---|---|---|---:|---:|---:|---:|---:|---|");
        let format_metric = |value: f64| {
            if value.is_finite() {
                format!("{value:.3e}")
            } else {
                "-".to_owned()
            }
        };

        let mut rst_passes = 0;
        let mut legacy_passes = 0;
        let options = EquilibriumSolveOptions::default()
            .with_max_iterations(500)
            .expect("TP-1907 backend matrix iteration budget must validate")
            .with_tolerance(1.0e-12)
            .expect("TP-1907 backend matrix tolerance must validate");

        for (name, family, backend) in backends {
            let outcome = calculator_tp1907_builder(&fixture, repository())
                .solve_options(options.clone())
                .solver_policy(SolverPolicy::Single(backend))
                .expect("TP-1907 single-backend policy must validate")
                .at_temperature(reference.temperature_k)
                .element_inventory(inventory.clone())
                .element_numerical_seed(initial.clone())
                .solve();

            let (
                status,
                scaled_residual,
                raw_residual,
                max_balance,
                max_component_relative,
                max_fraction_error,
                topology,
                detail,
            ) = match outcome {
                Ok(EquilibriumCalculatorOutcome::PtPoint(point)) => {
                    let solution = point.solution();
                    let validation = solution.accepted_solution().validation();
                    let moles = solution.component_moles();
                    let max_balance = (0..element_matrix.ncols())
                        .map(|element| {
                            (0..element_matrix.nrows())
                                .map(|species| element_matrix[(species, element)] * moles[species])
                                .sum::<f64>()
                                - reference_b[element]
                        })
                        .map(f64::abs)
                        .fold(0.0_f64, f64::max);
                    let mut max_component_relative = 0.0_f64;
                    let mut component_accuracy = true;
                    for (&expected, &actual) in molecular_solution
                        .component_moles()
                        .iter()
                        .zip(moles.iter())
                    {
                        let absolute = (expected - actual).abs();
                        let relative = absolute / expected.abs().max(actual.abs()).max(1.0e-30);
                        max_component_relative = max_component_relative.max(relative);
                        component_accuracy &= relative <= 1.0e-5 || absolute <= 1.0e-10;
                    }
                    let mut max_fraction_error = 0.0_f64;
                    let mut fraction_accuracy = true;
                    for component in molecular_solution.metadata().components() {
                        let expected_fraction = molecular_solution
                            .mole_fraction_for(component.id())
                            .expect("molecular fraction must exist");
                        let actual_fraction = solution
                            .mole_fraction_for(component.id())
                            .expect("backend fraction must exist");
                        let error = (expected_fraction - actual_fraction).abs();
                        max_fraction_error = max_fraction_error.max(error);
                        fraction_accuracy &= error <= 1.0e-8;
                    }
                    let topology_accuracy = expected_topology
                        .iter()
                        .all(|(phase, expected)| solution.phase_status(phase) == *expected);
                    let finite_non_negative =
                        moles.iter().all(|value| value.is_finite() && *value >= 0.0);
                    let accuracy = validation.residual_l2_norm <= 1.0e-5
                        && validation.raw_residual_l2_norm <= 1.0e-5
                        && max_balance <= 1.0e-6
                        && component_accuracy
                        && fraction_accuracy
                        && topology_accuracy
                        && finite_non_negative;
                    let status = if accuracy {
                        if family == "RST" {
                            rst_passes += 1;
                        } else {
                            legacy_passes += 1;
                        }
                        "PASS"
                    } else {
                        "ACCURACY_FAIL"
                    };
                    (
                        status,
                        validation.residual_l2_norm,
                        validation.raw_residual_l2_norm,
                        max_balance,
                        max_component_relative,
                        max_fraction_error,
                        topology_accuracy,
                        String::new(),
                    )
                }
                Ok(other) => (
                    "HARD_FAIL",
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    false,
                    format!("unexpected facade outcome: {other:?}"),
                ),
                Err(error) => {
                    let attempts = match &error {
                        crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorError::Solve(
                            solve_error,
                        )
                        | crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorError::Preparation(
                            solve_error,
                        ) => solve_error.backend_attempt_reports(),
                        _ => Vec::new(),
                    };
                    let accuracy_failure = attempts.iter().any(|attempt| {
                        matches!(
                            attempt.outcome,
                            SolverAttemptOutcome::RejectedCandidate { .. }
                        )
                    });
                    (
                        if accuracy_failure {
                            "ACCURACY_FAIL"
                        } else {
                            "HARD_FAIL"
                        },
                        f64::NAN,
                        f64::NAN,
                        f64::NAN,
                        f64::NAN,
                        f64::NAN,
                        false,
                        attempts
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                            .join("; "),
                    )
                }
            };
            println!(
                "| {name} | {family} | {status} | {} | {} | {} | {} | {} | {topology} |",
                format_metric(scaled_residual),
                format_metric(raw_residual),
                format_metric(max_balance),
                format_metric(max_component_relative),
                format_metric(max_fraction_error),
            );
            if !detail.is_empty() {
                println!("  {name} detail: {detail}");
            }
        }

        assert!(
            legacy_passes > 0,
            "TP-1907 upper inventory failed: no Legacy backend passed the physical acceptance matrix"
        );
        (rst_passes, legacy_passes)
    }

    #[test]
    fn tp1907_frozen_facade_small_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e-3, "small");
    }

    #[test]
    fn tp1907_frozen_facade_medium_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e-1, "medium");
    }

    #[test]
    fn tp1907_frozen_facade_large_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e1, "large");
    }

    #[test]
    fn tp1907_frozen_facade_very_large_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e2, "very-large");
    }

    #[test]
    fn tp1907_frozen_facade_extreme_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e3, "extreme");
    }

    #[test]
    fn tp1907_frozen_facade_beyond_extreme_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e4, "beyond-extreme");
    }

    #[test]
    fn tp1907_frozen_facade_high_inventory_matches_the_same_oracle() {
        assert_tp1907_facade_inventory_scale(1.0e6, "high");
    }

    #[test]
    fn tp1907_frozen_facade_upper_inventory_matches_the_same_oracle() {
        let (rst_passes, _) = run_tp1907_upper_backend_acceptance_matrix();
        assert!(
            rst_passes > 0,
            "TP-1907 upper inventory failed: no RST backend passed the physical acceptance matrix"
        );
    }

    #[test]
    fn tp1907_atom_view_and_expr_legacy_have_componentwise_callback_parity_at_the_seed() {
        let fixture = fixture();
        let prepared = prepared_fixed_tp1907_problem_at_scale(&fixture, 700.0, 1.0e8);
        let (legacy, _) = prepare_rst_symbolic_problem_from_prepared_with_frontend(
            &prepared,
            crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
            SymbolicLambdifyFrontend::ExprLegacy,
        )
        .expect("TP-1907 ExprLegacy parity problem must prepare");
        let (atom, _) = prepare_rst_symbolic_problem_from_prepared_with_frontend(
            &prepared,
            crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
            SymbolicLambdifyFrontend::AtomViewNative,
        )
        .expect("TP-1907 AtomView parity problem must prepare");
        let seed = prepared.problem().initial_log_moles().as_slice();
        let legacy_residual = legacy.residual_for_test(seed).unwrap();
        let atom_residual = atom.residual_for_test(seed).unwrap();
        let residual_drift = legacy_residual
            .iter()
            .zip(atom_residual.iter())
            .map(|(left, right)| (left - right).abs())
            .fold(0.0_f64, f64::max);
        let residual_scale = legacy_residual
            .iter()
            .chain(atom_residual.iter())
            .map(|value| value.abs())
            .fold(1.0_f64, f64::max);
        let residual_relative_drift = residual_drift / residual_scale;
        let legacy_jacobian = legacy.jacobian_for_test(seed).unwrap();
        let atom_jacobian = atom.jacobian_for_test(seed).unwrap();
        let jacobian_drift = legacy_jacobian
            .iter()
            .zip(atom_jacobian.iter())
            .map(|(left, right)| (left - right).abs())
            .fold(0.0_f64, f64::max);
        let jacobian_scale = legacy_jacobian
            .iter()
            .chain(atom_jacobian.iter())
            .map(|value| value.abs())
            .fold(1.0_f64, f64::max);
        let jacobian_relative_drift = jacobian_drift / jacobian_scale;

        println!(
            "TP-1907 upper seed frontend parity max_residual_abs={residual_drift:.6e} max_residual_relative={residual_relative_drift:.6e} max_jacobian_abs={jacobian_drift:.6e} max_jacobian_relative={jacobian_relative_drift:.6e}"
        );
        assert!(residual_relative_drift <= 2e-9);
        assert!(jacobian_relative_drift <= 1e-12);
        assert!(legacy_residual.len() == atom_residual.len());
        assert!(legacy_jacobian.shape() == atom_jacobian.shape());
    }

    #[test]
    fn tp1907_upper_inventory_facade_atom_view_and_expr_legacy_match_published_equilibrium() {
        let fixture = fixture();
        let temperature_k = 700.0;
        let scale = 1.0e8;
        let initial = scaled_initial_composition(&fixture, scale);
        let (inventory, _element_labels, reference_b) =
            frozen_element_inventory(&fixture, &initial);
        let solve = |frontend| {
            with_test_rst_frontend(frontend, || {
                calculator_tp1907_pt_with_element_seed(
                    &fixture,
                    repository(),
                    temperature_k,
                    inventory.clone(),
                    initial.clone(),
                )
            })
        };
        let expr_point = solve(SymbolicLambdifyFrontend::ExprLegacy);
        let atom_point = solve(SymbolicLambdifyFrontend::AtomViewNative);
        let expr = expr_point.solution();
        let atom = atom_point.solution();
        let expr_components = expr.metadata().components();
        let atom_components = atom.metadata().components();
        assert_eq!(expr_components.len(), atom_components.len());

        let mut max_fraction_difference = 0.0_f64;
        let mut max_component_relative_difference = 0.0_f64;
        for (expr_component, atom_component) in expr_components.iter().zip(atom_components) {
            assert_eq!(expr_component.id(), atom_component.id());
            let expr_fraction = expr
                .mole_fraction_for(expr_component.id())
                .expect("ExprLegacy fraction must exist");
            let atom_fraction = atom
                .mole_fraction_for(atom_component.id())
                .expect("AtomView fraction must exist");
            max_fraction_difference =
                max_fraction_difference.max((expr_fraction - atom_fraction).abs());
        }
        for (&expr_moles, &atom_moles) in expr.component_moles().iter().zip(atom.component_moles())
        {
            let relative = (expr_moles - atom_moles).abs()
                / expr_moles.abs().max(atom_moles.abs()).max(1.0e-30);
            max_component_relative_difference = max_component_relative_difference.max(relative);
        }
        let topology_matches = expr
            .metadata()
            .phases()
            .iter()
            .all(|phase| expr.phase_status(phase.id()) == atom.phase_status(phase.id()));
        let temperature_difference =
            (expr.conditions().temperature() - atom.conditions().temperature()).abs();
        let expr_validation = expr.accepted_solution().validation();
        let atom_validation = atom.accepted_solution().validation();
        let expr_backend = expr.solve_report().accepted_backend;
        let atom_backend = atom.solve_report().accepted_backend;

        println!(
            "TP-1907 facade frontend A/B scale={scale:.1e} ExprLegacy_backend={expr_backend:?} AtomView_backend={atom_backend:?} ExprLegacy_T={:.9} K AtomView_T={:.9} K delta_T={temperature_difference:.3e} max_abs_fraction_delta={max_fraction_difference:.6e} max_component_relative_delta={max_component_relative_difference:.6e} ExprLegacy_balance={:.6e} AtomView_balance={:.6e} topology_matches={topology_matches}",
            expr.conditions().temperature(),
            atom.conditions().temperature(),
            expr_validation.max_abs_element_balance_error,
            atom_validation.max_abs_element_balance_error,
        );
        for component in expr_components {
            let expr_fraction = expr
                .mole_fraction_for(component.id())
                .expect("ExprLegacy fraction must exist");
            let atom_fraction = atom
                .mole_fraction_for(component.id())
                .expect("AtomView fraction must exist");
            println!(
                "  {:32} ExprLegacy={expr_fraction:.12e} AtomView={atom_fraction:.12e} abs_delta={:.6e}",
                component.label(),
                (expr_fraction - atom_fraction).abs(),
            );
        }
        assert_eq!(expr.conditions().temperature(), temperature_k);
        assert_eq!(atom.conditions().temperature(), temperature_k);
        assert_eq!(expr_backend, atom_backend);
        let balance_tolerance = 1.0e-12 * reference_b.iter().sum::<f64>();
        assert!(expr_validation.max_abs_element_balance_error <= balance_tolerance);
        assert!(atom_validation.max_abs_element_balance_error <= balance_tolerance);
        assert!(topology_matches);
        assert!(max_fraction_difference <= 1.0e-8);
        assert!(max_component_relative_difference <= 1.0e-5);
    }

    #[test]
    fn tp1907_upper_inventory_records_frontend_acceptance_parity() {
        let expr = with_test_rst_frontend(SymbolicLambdifyFrontend::ExprLegacy, || {
            run_tp1907_upper_backend_acceptance_matrix()
        });
        let atom = with_test_rst_frontend(SymbolicLambdifyFrontend::AtomViewNative, || {
            run_tp1907_upper_backend_acceptance_matrix()
        });
        println!(
            "TP-1907 upper inventory frontend acceptance ExprLegacy RST_passes={} Legacy_passes={} AtomView RST_passes={} Legacy_passes={}",
            expr.0, expr.1, atom.0, atom.1
        );
        assert!(expr.0 > 0, "ExprLegacy must retain its accepted RST route");
        assert!(
            expr.1 > 0 && atom.1 > 0,
            "legacy controls must pass for both runs"
        );
        assert_eq!(
            expr.0, atom.0,
            "RST frontend acceptance count changed at the upper inventory scale"
        );
        assert!(atom.0 > 0, "AtomView must retain an accepted RST backend");
        assert_eq!(
            expr.1, atom.1,
            "legacy control acceptance changed by frontend"
        );
    }

    /// Exports the accepted KiThe facade state as a deterministic RST parity
    /// fixture. This is intentionally ignored: it is a diagnostic hand-off,
    /// not another production solve or correctness gate.
    #[test]
    #[ignore = "emit the full KiThe/RST parity fixture for the RustedSciThe test suite"]
    fn tp1907_upper_inventory_emits_paste_ready_rst_parity_fixture() {
        let fixture = fixture();
        let inventory_scale = 1.0e8;
        let temperature_k = 700.0;
        let pressure_pa = TP1907_REFERENCE_PRESSURE_PA;
        let initial_composition = scaled_initial_composition(&fixture, inventory_scale);
        let (inventory, element_labels, element_totals_b) =
            frozen_element_inventory(&fixture, &initial_composition);
        let normalization = ExtensiveNormalization::from_element_inventory(&inventory)
            .expect("TP-1907 parity inventory normalization must validate");
        let prepared =
            prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, inventory_scale);
        let rst_problem = prepare_rst_symbolic_problem_from_prepared(&prepared)
            .expect("TP-1907 parity RST problem must prepare");
        let (equations, variables, parameter_values, parameters) =
            rst_problem.symbolic_payload_for_test();
        let initial_log_moles = prepared.problem().initial_log_moles().as_slice().to_vec();
        let bounds = rst_problem.log_mole_bounds_for_test();
        let (element_matrix, _, matrix_element_labels) = element_composition_and_molar_mass(
            fixture.resolved().phase_data(),
            fixture.resolved().layout(),
            None,
        )
        .expect("TP-1907 parity element matrix must remain available");
        assert_eq!(element_labels, matrix_element_labels);
        assert_eq!(variables.len(), 18);
        assert_eq!(equations.len(), 18);
        assert_eq!(initial_log_moles.len(), 18);
        assert_eq!(element_totals_b.len(), element_labels.len());

        let options = EquilibriumSolveOptions::default()
            .with_max_iterations(500)
            .expect("TP-1907 parity iteration budget must validate")
            .with_tolerance(1.0e-12)
            .expect("TP-1907 parity tolerance must validate");
        RustedSciTheSolver::begin_history_capture_for_test();
        let outcome = calculator_tp1907_builder(&fixture, repository())
            .solve_options(options)
            .solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
                RustedSciTheSolver::TrustRegionLevenbergMarquardt,
            )))
            .expect("TP-1907 parity solver policy must validate")
            .at_temperature(temperature_k)
            .element_inventory(inventory.clone())
            .element_numerical_seed(initial_composition.clone())
            .solve()
            .expect("TP-1907 accepting RST facade route must solve");
        let captured_results = RustedSciTheSolver::take_history_capture_for_test();
        let point = match outcome {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected TP-1907 parity point, got {other:?}"),
        };
        let solution = point.solution();
        let accepted = solution.accepted_solution();
        let y_final = accepted.log_moles().to_vec();
        let physical_moles = y_final.iter().map(|value| value.exp()).collect::<Vec<_>>();
        let numerical_moles = accepted.moles().to_vec();
        let published_moles = solution.component_moles().to_vec();
        assert_eq!(y_final.len(), 18);
        assert!(y_final.iter().all(|value| value.is_finite()));
        assert!(
            physical_moles
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        );

        let raw_residual = rst_problem
            .residual_for_test(&y_final)
            .expect("TP-1907 parity raw residual must be evaluable");
        let raw_jacobian = rst_problem
            .jacobian_for_test(&y_final)
            .expect("TP-1907 parity raw Jacobian must be evaluable");
        let row_scale = prepared
            .residual_scale()
            .expect("TP-1907 parity row scaling must remain available");
        let scaled_residual = prepared
            .scaled_residual(&y_final, &row_scale)
            .expect("TP-1907 parity scaled residual must be evaluable");
        let analytic_scaled_jacobian = prepared
            .scaled_jacobian(&y_final, &row_scale)
            .expect("TP-1907 parity scaled Jacobian must be evaluable");
        assert_eq!(raw_residual.len(), 18);
        assert_eq!(raw_jacobian.nrows(), 18);
        assert_eq!(raw_jacobian.ncols(), 18);
        assert_eq!(scaled_residual.len(), 18);
        assert_eq!(analytic_scaled_jacobian.nrows(), 18);
        assert_eq!(analytic_scaled_jacobian.ncols(), 18);

        let raw_jt_f = raw_jacobian.transpose() * &raw_residual;
        let scaled_jt_f = analytic_scaled_jacobian.transpose()
            * nalgebra::DVector::from_vec(scaled_residual.clone());
        let numerical_balance = (0..element_matrix.ncols())
            .map(|element| {
                (0..element_matrix.nrows())
                    .map(|species| element_matrix[(species, element)] * physical_moles[species])
                    .sum::<f64>()
                    - element_totals_b[element]
            })
            .collect::<Vec<_>>();
        let max_abs_element_balance_error = numerical_balance
            .iter()
            .map(|value| value.abs())
            .fold(0.0_f64, f64::max);
        let published_balance = (0..element_matrix.ncols())
            .map(|element| {
                (0..element_matrix.nrows())
                    .map(|species| element_matrix[(species, element)] * published_moles[species])
                    .sum::<f64>()
                    - element_totals_b[element]
            })
            .collect::<Vec<_>>();
        let max_abs_published_element_balance_error = published_balance
            .iter()
            .map(|value| value.abs())
            .fold(0.0_f64, f64::max);
        let validation = accepted.validation();
        let graphite = PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()));
        let graphite_status = solution.phase_status(&graphite);
        let graphite_total = solution.phase_total(&graphite);
        let graphite_tpd = graphite_final_tpd(solution);

        let contract = RustedSciTheSolveContract::new(1.0e-12, 500)
            .expect("TP-1907 parity RST contract must validate");
        let direct_result = RustedSciTheSolver::TrustRegionLevenbergMarquardt
            .solve_with_history_for_test(&rst_problem, &initial_log_moles, &contract)
            .expect("TP-1907 direct RST parity replay must return a result");
        let normalization_scale = normalization.physical_inventory_scale();
        let normalized_result = captured_results
            .iter()
            .filter(|result| format!("{:?}", result.termination) == "Converged")
            .max_by_key(|result| result.iterations)
            .expect("TP-1907 facade recovery must capture a converged normalized RST attempt");
        let normalized_reconstructed_y = normalized_result
            .x
            .iter()
            .map(|value| value + normalization_scale.ln())
            .collect::<Vec<_>>();
        let max_normalized_reconstructed_y_error = y_final
            .iter()
            .zip(normalized_reconstructed_y.iter())
            .map(|(&facade, &reconstructed)| (facade - reconstructed).abs())
            .fold(0.0_f64, f64::max);
        let max_facade_direct_y_error = y_final
            .iter()
            .zip(direct_result.x.iter())
            .map(|(&facade, &direct)| (facade - direct).abs())
            .fold(0.0_f64, f64::max);

        let format_vec = |values: &[f64]| {
            values
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(", ")
        };
        let format_matrix = |matrix: &nalgebra::DMatrix<f64>| {
            matrix
                .row_iter()
                .map(|row| {
                    format!(
                        "    vec![{}]",
                        row.iter()
                            .map(|value| format!("{value:.17e}"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                })
                .collect::<Vec<_>>()
                .join(",\n")
        };
        let phase_status_rows = solution
            .phases()
            .iter()
            .map(|phase| {
                let phase_id = format!("{:?}", phase.id());
                let physical_state = format!("{:?}", phase.physical_state());
                let phase_model = format!("{:?}", phase.phase_model());
                let activity_model = format!("{:?}", phase.activity_model());
                let status = format!("{:?}", solution.phase_status(phase.id()));
                format!(
                    "    ({phase_id:?}.to_string(), {physical_state:?}.to_string(), {phase_model:?}.to_string(), {activity_model:?}.to_string(), {status:?}.to_string(), {:?}),",
                    solution.phase_total(phase.id())
                )
            })
            .collect::<Vec<_>>()
            .join("\n");
        let accepted_backend_name = format!("{:?}", solution.solve_report().accepted_backend);
        let recovery_used = solution.extensive_normalization_recovery().is_some();
        let recovery_discovery_backend = solution
            .extensive_normalization_recovery()
            .map(|evidence| format!("{:?}", evidence.discovery_backend));
        let recovery_physical_retry_backend = solution
            .extensive_normalization_recovery()
            .and_then(|evidence| evidence.physical_retry_backend.as_ref())
            .map(|backend| format!("{:?}", backend));
        let reconstructed_physical_boundary = solution
            .extensive_normalization_recovery()
            .map(|evidence| evidence.reconstructed_physical_boundary)
            .unwrap_or(false);

        println!("TP1907 KiThe/RST strict parity fixture");
        println!("// source: NASA TP-1907 CHON + graphite, accepted fixed P,T facade route");
        println!(
            "// requested policy: elemental ElementInventory -> explicit Single RST Trust-Region LM"
        );
        println!(
            "// actual accepted route: physical attempt -> extensive normalization -> normalized RST discovery -> failed physical retry -> reconstructed physical publication"
        );
        println!(
            "// temperature_k: {temperature_k:.17e}; pressure_pa: {pressure_pa:.17e}; inventory_scale: {inventory_scale:.17e}"
        );
        println!(
            "// The fixture is emitted from the accepted facade state; no second solver determines y_final."
        );
        println!(
            "// Coordinates: y[i] = ln(n[i]); n[i] = exp(y[i]); no dependent-variable elimination."
        );
        println!(
            "// Rows 0..13 are reaction-affinity equations; rows 13..18 are elemental balances A*n-b."
        );
        println!("// F_raw(y) = [reaction-affinity rows; A*n-b].");
        println!(
            "// residual_row_scale_w is denominator D: F_scaled[i] = F_raw[i] / D[i], W=diag(1/D)."
        );
        println!(
            "// Reaction D rows are reaction-norm scales; balance D rows are abs(b[e]) with unit fallback at zero."
        );
        println!(
            "// Bounds are adapter log-mole roundoff bounds; exp(y) is strictly positive before publication."
        );
        println!(
            "// Published inactive-phase amounts may be zero; phase control is external active-set/TPD logic."
        );
        println!(
            "// At this point graphite must be active and its final TPD is checked by the facade acceptance report."
        );
        println!(
            "// This is the same 18-variable/18-row prepared graph exported by the RST adapter."
        );
        println!(
            "// The facade acceptance report may retain normalized-recovery solver norms; the raw/scaled vectors and physical balances below are freshly recomputed on this physical y_final."
        );
        println!("let accepted_backend: &str = {accepted_backend_name:?};");
        println!(
            "let active_topology: &str = {:?};",
            active_phase_topology(solution)
        );
        println!(
            "let started_backend_attempts: usize = {};",
            solution.started_backend_attempts()
        );
        println!("let used_extensive_normalization_recovery: bool = {recovery_used};");
        println!("let normalization_scale: f64 = {normalization_scale:.17e};");
        println!(
            "let recovery_discovery_backend: Option<&str> = {:?};",
            recovery_discovery_backend.as_deref()
        );
        println!(
            "let recovery_physical_retry_backend: Option<&str> = {:?};",
            recovery_physical_retry_backend.as_deref()
        );
        println!("let reconstructed_physical_boundary: bool = {reconstructed_physical_boundary};");
        println!(
            "// accepted backend: {:?}; started attempts: {}; active topology: {}.",
            solution.solve_report().accepted_backend,
            solution.started_backend_attempts(),
            active_phase_topology(solution)
        );
        println!(
            "// extensive normalization/recovery evidence: {:?}",
            solution.extensive_normalization_recovery()
        );
        println!(
            "let equations_fnv1a64: &str = \"0x{:016x}\";",
            stable_equation_hash(&equations)
        );
        println!(
            "let equations: Vec<String> = vec![{}];",
            equations
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let variables: Vec<String> = vec![{}];",
            variables
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let equation_parameters: Vec<String> = vec![{}];",
            parameters
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let equation_parameter_values: Vec<f64> = vec![{}];",
            format_vec(&parameter_values)
        );
        println!(
            "let element_labels: Vec<String> = vec![{}];",
            element_labels
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let element_totals_b: Vec<f64> = vec![{}];",
            format_vec(&element_totals_b)
        );
        println!(
            "let species_order: Vec<String> = vec![{}];",
            prepared
                .problem()
                .species()
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!("let element_composition_a_element_by_species: Vec<Vec<f64>> = vec![");
        for element in 0..element_matrix.ncols() {
            println!(
                "    vec![{}], // A[{}={}]",
                (0..element_matrix.nrows())
                    .map(|species| format!("{:.17e}", element_matrix[(species, element)]))
                    .collect::<Vec<_>>()
                    .join(", "),
                element,
                element_labels[element]
            );
        }
        println!("];");
        println!(
            "let species_phase: Vec<usize> = vec![{}];",
            prepared
                .species_phase()
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let residual_row_scale_w: Vec<f64> = vec![{}];",
            format_vec(&row_scale)
        );
        println!(
            "let initial_log_moles: Vec<f64> = vec![{}];",
            format_vec(&initial_log_moles)
        );
        println!(
            "let log_mole_bounds: Option<Vec<(f64, f64)>> = {};",
            bounds.as_ref().map_or_else(
                || "None".to_owned(),
                |values| format!(
                    "Some(vec![{}])",
                    values
                        .iter()
                        .map(|(lower, upper)| format!("({lower:.17e}, {upper:.17e})"))
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            )
        );
        println!("let y_final: Vec<f64> = vec![{}];", format_vec(&y_final));
        println!(
            "let physical_moles_from_y_final: Vec<f64> = vec![{}];",
            format_vec(&physical_moles)
        );
        println!(
            "let numerical_moles_from_accepted_solution: Vec<f64> = vec![{}];",
            format_vec(&numerical_moles)
        );
        println!(
            "let published_component_moles: Vec<f64> = vec![{}];",
            format_vec(&published_moles)
        );
        println!(
            "let phase_status_rows: Vec<(String, String, String, String, String, Option<f64>)> = vec!["
        );
        print!("{phase_status_rows}");
        println!("];");
        println!(
            "let graphite_status: Option<String> = {:?};",
            graphite_status.map(|status| format!("{status:?}"))
        );
        println!("let graphite_total: Option<f64> = {:?};", graphite_total);
        println!("let graphite_final_minimum_tpd: f64 = {graphite_tpd:.17e};");
        println!(
            "let raw_residual: Vec<f64> = vec![{}];",
            format_vec(raw_residual.as_slice())
        );
        println!(
            "let scaled_residual: Vec<f64> = vec![{}];",
            format_vec(scaled_residual.as_slice())
        );
        println!(
            "let raw_jacobian: Vec<Vec<f64>> = vec![\n{}\n];",
            format_matrix(&raw_jacobian)
        );
        println!(
            "let scaled_jacobian: Vec<Vec<f64>> = vec![\n{}\n];",
            format_matrix(&analytic_scaled_jacobian)
        );
        println!(
            "let jt_f_raw: Vec<f64> = vec![{}];",
            format_vec(raw_jt_f.as_slice())
        );
        println!(
            "let jt_f_scaled: Vec<f64> = vec![{}];",
            format_vec(scaled_jt_f.as_slice())
        );
        println!(
            "let raw_residual_l2: f64 = {:.17e};",
            l2_norm(raw_residual.as_slice())
        );
        println!(
            "let scaled_residual_l2: f64 = {:.17e};",
            l2_norm(scaled_residual.as_slice())
        );
        println!(
            "let raw_element_balance: Vec<f64> = vec![{}];",
            format_vec(&numerical_balance)
        );
        println!("let max_abs_element_balance_error: f64 = {max_abs_element_balance_error:.17e};");
        println!(
            "let published_element_balance: Vec<f64> = vec![{}];",
            format_vec(&published_balance)
        );
        println!(
            "let max_abs_published_element_balance_error: f64 = {max_abs_published_element_balance_error:.17e};"
        );
        println!(
            "// KiThe validation: scaled residual L2 = {:.17e}; raw residual L2 = {:.17e}; max reported physical balance = {:.17e}.",
            validation.residual_l2_norm,
            validation.raw_residual_l2_norm,
            validation.max_abs_element_balance_error
        );
        println!(
            "// Current facade publication acceptance uses published moles: max |A*n_published-b| <= 1e-6. A strict RST exp(y_final) transfer check must use max_abs_element_balance_error separately."
        );
        println!(
            "// Proposed transfer acceptance: raw/scaled residual L2 <= 1e-5, physical max |A*n-b| <= 1e-6, finite/non-negative publication, exact ordering/topology, TPD phase-control checks."
        );
        println!(
            "// In this run exp(y_final) balance is {:.17e}, while published balance is {:.17e}; retain both values to expose log/exp round-trip error rather than silently substituting publication moles.",
            max_abs_element_balance_error, max_abs_published_element_balance_error
        );
        println!(
            "// Direct RST replay starts from the original physical seed, while y_final is the accepted facade state; this is not asserted equal because facade basin preparation/recovery may change coordinates."
        );
        println!("let max_facade_direct_y_error: f64 = {max_facade_direct_y_error:.17e};");
        println!(
            "// Direct raw RST replay termination = {:?}; residual_l2 = {:.17e}.",
            direct_result.termination, direct_result.residual_norm
        );
        println!(
            "let normalized_discovery_reconstructed_y: Vec<f64> = vec![{}];",
            format_vec(&normalized_reconstructed_y)
        );
        println!(
            "let max_normalized_reconstructed_y_error: f64 = {max_normalized_reconstructed_y_error:.17e};"
        );
        println!(
            "// The normalized discovery replay uses the same elemental inventory scaling and numerical seed as the accepted recovery route; its denormalized y is compared with y_final below."
        );
        println!(
            "// RST 0.4.17 history exposes iteration, post-step residual norm, step norm, and accepted flag only."
        );
        println!(
            "// lambda, regularization form, trial alpha, exact residual-before/after for internal retries, and per-iteration termination are not public in this version; they are emitted as None below rather than inferred."
        );
        println!(
            "#[derive(Debug, Clone)] struct IterationTraceRow {{ iteration: usize, residual_norm_before: Option<f64>, lambda: Option<f64>, regularization: Option<&'static str>, trial_alpha: Option<f64>, accepted: bool, residual_norm_after: f64, step_norm: f64, termination: Option<&'static str> }}"
        );
        println!("let iteration_trace_last_20: Vec<IterationTraceRow> = vec![");
        let history_start = normalized_result.history.len().saturating_sub(20);
        for (position, record) in normalized_result
            .history
            .iter()
            .enumerate()
            .skip(history_start)
        {
            let before = if position > 0 {
                Some(normalized_result.history[position - 1].residual_norm)
            } else {
                None
            };
            println!(
                "    IterationTraceRow {{ iteration: {}, residual_norm_before: {:?}, lambda: None, regularization: None, trial_alpha: None, accepted: {}, residual_norm_after: {:.17e}, step_norm: {:.17e}, termination: None }},",
                record.iteration, before, record.accepted, record.residual_norm, record.step_norm
            );
        }
        println!("];");
        println!(
            "let final_termination: &str = {:?};",
            format!("{:?}", normalized_result.termination)
        );
        println!(
            "// Aggregate RST calls: residual={} jacobian={} linear_solves={} iterations={} accepted_steps={} rejected_steps={} total_duration_us={}.",
            normalized_result.statistics.residual_evaluations,
            normalized_result.statistics.jacobian_evaluations,
            normalized_result.statistics.linear_solves,
            normalized_result.statistics.iterations,
            normalized_result.statistics.accepted_steps,
            normalized_result.statistics.rejected_steps,
            normalized_result.statistics.total_duration.as_micros()
        );
    }

    /// Proves that the exact equation strings exported from KiThe preserve
    /// both residual and Jacobian semantics when rebuilt through RST's public
    /// parser. This catches a lossy `Expr::to_string()` hand-off before an RST
    /// reproducer is used to diagnose solver behavior.
    #[test]
    #[ignore = "currently exposes RST parser incompatibility with Expr::to_string() output"]
    fn tp1907_upper_inventory_symbolic_string_round_trip_matches_original_graph() {
        let fixture = fixture();
        let inventory_scale = 1.0e8;
        let temperature_k = 700.0;
        let initial_composition = scaled_initial_composition(&fixture, inventory_scale);
        let (inventory, _, _) = frozen_element_inventory(&fixture, &initial_composition);
        let prepared =
            prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, inventory_scale);
        let rst_problem = prepare_rst_symbolic_problem_from_prepared(&prepared)
            .expect("TP-1907 round-trip source graph must prepare");
        let (equations, variables, parameter_values, parameters) =
            rst_problem.symbolic_payload_for_test();
        let equation_hash = stable_equation_hash(&equations);
        println!("TP1907 equation_hash_fnv1a64=0x{equation_hash:016x}");

        for (index, equation) in equations.iter().enumerate() {
            let parsed = catch_unwind(AssertUnwindSafe(|| Expr::parse_expression(equation)));
            if parsed.is_err() {
                panic!("TP-1907 exported equation {index} is rejected by RST parser: {equation:?}");
            }
        }

        let options = SymbolicProblemOptions::new()
            .with_variables(variables.clone())
            .with_equation_parameters(parameters.clone())
            .with_equation_parameter_values(DVector::from_vec(parameter_values.clone()))
            .with_lambdify_execution_policy(LambdifyExecutionPolicy::Parallel { min_work: 1 })
            .with_lambdify_backend();
        let reparsed_problem =
            SymbolicNonlinearProblem::from_strings_with_options(equations.to_vec(), options)
                .expect("TP-1907 exported equations must parse through public RST string API");

        let solve_options = EquilibriumSolveOptions::default()
            .with_max_iterations(500)
            .expect("TP-1907 round-trip iteration budget must validate")
            .with_tolerance(1.0e-12)
            .expect("TP-1907 round-trip tolerance must validate");
        let outcome = calculator_tp1907_builder(&fixture, repository())
            .solve_options(solve_options)
            .solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
                RustedSciTheSolver::TrustRegionLevenbergMarquardt,
            )))
            .expect("TP-1907 round-trip solver policy must validate")
            .at_temperature(temperature_k)
            .element_inventory(inventory)
            .element_numerical_seed(initial_composition)
            .solve()
            .expect("TP-1907 round-trip facade route must accept the solution");
        let point = match outcome {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected TP-1907 round-trip point, got {other:?}"),
        };
        let y_final = DVector::from_vec(point.solution().accepted_solution().log_moles().to_vec());
        assert_eq!(y_final.len(), equations.len());
        assert_eq!(y_final.len(), variables.len());

        let original_residual = rst_problem
            .residual_for_test(y_final.as_slice())
            .expect("original TP-1907 graph residual must evaluate");
        let original_jacobian = rst_problem
            .jacobian_for_test(y_final.as_slice())
            .expect("original TP-1907 graph Jacobian must evaluate");
        let reparsed_residual = reparsed_problem
            .residual(&y_final)
            .expect("reparsed TP-1907 graph residual must evaluate");
        let reparsed_jacobian = reparsed_problem
            .jacobian(&y_final)
            .expect("reparsed TP-1907 graph Jacobian must evaluate");

        let residual_differences = original_residual
            .iter()
            .zip(reparsed_residual.iter())
            .map(|(original, reparsed)| (original - reparsed).abs())
            .collect::<Vec<_>>();
        let max_abs_residual_difference =
            residual_differences.iter().copied().fold(0.0_f64, f64::max);
        let mut max_abs_jacobian_difference = 0.0_f64;
        let mut first_jacobian_difference = None;
        for row in 0..original_jacobian.nrows() {
            for column in 0..original_jacobian.ncols() {
                let original = original_jacobian[(row, column)];
                let reparsed = reparsed_jacobian[(row, column)];
                let difference = (original - reparsed).abs();
                if difference > max_abs_jacobian_difference {
                    max_abs_jacobian_difference = difference;
                }
                if first_jacobian_difference.is_none() && difference != 0.0 {
                    first_jacobian_difference = Some((row, column, original, reparsed, difference));
                }
            }
        }
        let first_residual_difference = residual_differences
            .iter()
            .enumerate()
            .find(|(_, difference)| **difference != 0.0)
            .map(|(row, difference)| {
                (
                    row,
                    original_residual[row],
                    reparsed_residual[row],
                    *difference,
                )
            });
        println!("TP1907 KiThe/RST symbolic string round-trip");
        println!(
            "case=upper-inventory scale={inventory_scale:.17e} temperature_k={temperature_k:.17e}"
        );
        println!(
            "equation_count={} variable_count={} parameter_count={}",
            equations.len(),
            variables.len(),
            parameters.len()
        );
        println!("equation_hash_fnv1a64=0x{equation_hash:016x}");
        println!("max_abs_residual_difference={max_abs_residual_difference:.17e}");
        println!("max_abs_jacobian_difference={max_abs_jacobian_difference:.17e}");
        println!(
            "original_residual_l2={:.17e}",
            l2_norm(original_residual.as_slice())
        );
        println!(
            "reparsed_residual_l2={:.17e}",
            l2_norm(reparsed_residual.as_slice())
        );
        println!("first_residual_difference={first_residual_difference:?}");
        println!("first_jacobian_difference={first_jacobian_difference:?}");
        println!("variables={variables:?}");
        println!("equation_parameters={parameters:?}");
        println!("equation_parameter_values={parameter_values:?}");
        println!("canonical_equations:");
        for (index, equation) in equations.iter().enumerate() {
            println!("  [{index}] {equation:?}");
        }

        assert!(
            max_abs_residual_difference <= 1.0e-12,
            "TP-1907 Expr::to_string() residual round-trip drift is {max_abs_residual_difference:.17e}; first={first_residual_difference:?}"
        );
        assert!(
            max_abs_jacobian_difference <= 1.0e-12,
            "TP-1907 Expr::to_string() Jacobian round-trip drift is {max_abs_jacobian_difference:.17e}; first={first_jacobian_difference:?}"
        );
    }

    /// Emits a self-contained RST reproducer for the currently fragile upper
    /// inventory scale. This is intentionally ignored: it is a hand-off
    /// artifact for RustedSciThe solver development, not a production gate.
    #[test]
    #[ignore = "paste the emitted symbolic payload into the RustedSciThe test suite"]
    fn tp1907_upper_inventory_emits_paste_ready_rst_symbolic_problem() {
        let fixture = fixture();
        let inventory_scale = 1.0e8;
        let initial_composition = scaled_initial_composition(&fixture, inventory_scale);
        let (inventory, element_labels, element_totals_b) =
            frozen_element_inventory(&fixture, &initial_composition);
        let normalization = ExtensiveNormalization::from_element_inventory(&inventory)
            .expect("TP-1907 elemental inventory must define normalization");
        let prepared = prepared_fixed_tp1907_problem_at_scale(&fixture, 700.0, inventory_scale);
        let rst_problem = prepare_rst_symbolic_problem_from_prepared(&prepared)
            .expect("TP-1907 upper-scale RST problem must prepare");
        let (equations, variables, parameter_values, parameters) =
            rst_problem.symbolic_payload_for_test();
        let initial = prepared.problem().initial_log_moles().as_slice().to_vec();
        let bounds = rst_problem.log_mole_bounds_for_test();
        let (element_matrix, _, matrix_element_labels) = element_composition_and_molar_mass(
            fixture.resolved().phase_data(),
            fixture.resolved().layout(),
            None,
        )
        .expect("TP-1907 element matrix must remain available");
        assert_eq!(element_labels, matrix_element_labels);

        let singular_values = rst_problem
            .jacobian_for_test(&initial)
            .expect("TP-1907 Jacobian must remain evaluable")
            .svd(false, false)
            .singular_values;
        let largest_singular_value = singular_values.iter().copied().fold(0.0_f64, f64::max);
        let rank_tolerance = 1.0e-12 * largest_singular_value.max(1.0);
        let numerical_rank = singular_values
            .iter()
            .filter(|&&value| value > rank_tolerance)
            .count();
        assert_eq!(equations.len(), 18);
        assert_eq!(variables.len(), 18);
        assert_eq!(numerical_rank, 17);

        let oracle_options = EquilibriumSolveOptions::default()
            .with_max_iterations(500)
            .expect("TP-1907 oracle iteration budget must validate")
            .with_tolerance(1.0e-12)
            .expect("TP-1907 oracle tolerance must validate");
        let molecular_oracle = calculator_tp1907_builder(&fixture, repository())
            .solve_options(oracle_options)
            .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::TR)))
            .expect("TP-1907 independent oracle policy must validate")
            .at_temperature(700.0)
            .initial_moles(initial_composition.moles().to_vec())
            .solve()
            .expect("TP-1907 legacy-TR molecular oracle must solve");
        let molecular_oracle = match molecular_oracle {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected TP-1907 molecular oracle point, got {other:?}"),
        };
        let oracle_solution = molecular_oracle.solution();
        let oracle_log_moles = oracle_solution
            .accepted_solution()
            .log_moles()
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let oracle_numerical_moles = oracle_solution
            .accepted_solution()
            .moles()
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let oracle_published_moles = oracle_solution
            .component_moles()
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let oracle_mole_fractions = oracle_solution
            .metadata()
            .components()
            .iter()
            .map(|component| {
                oracle_solution
                    .mole_fraction_for(component.id())
                    .expect("TP-1907 oracle mole fraction must exist")
            })
            .collect::<Vec<_>>();
        let oracle_raw_residual = prepared
            .residual(&oracle_log_moles)
            .expect("TP-1907 oracle raw residual must remain evaluable");
        let row_scale = prepared
            .residual_scale()
            .expect("TP-1907 row scaling contract must remain available");
        let oracle_scaled_residual = prepared
            .scaled_residual(&oracle_log_moles, &row_scale)
            .expect("TP-1907 oracle scaled residual must remain evaluable");
        let oracle_element_balance = (0..element_matrix.ncols())
            .map(|element| {
                (0..element_matrix.nrows())
                    .map(|species| {
                        element_matrix[(species, element)] * oracle_numerical_moles[species]
                    })
                    .sum::<f64>()
                    - element_totals_b[element]
            })
            .collect::<Vec<_>>();
        let phase_manager = PhaseControlPolicy::default().into_phase_manager();
        let (dg_create, dg_keep) = phase_manager
            .thresholds_at(700.0)
            .expect("TP-1907 phase hysteresis thresholds must resolve");

        assert_eq!(equations.len(), variables.len());
        assert!(!equations.is_empty());
        assert_eq!(parameters.len(), parameter_values.len());

        println!("TP1907 RST paste-ready reproducer");
        println!("// source: NASA TP-1907 CHON + graphite, fixed P,T");
        println!(
            "// physical inventory multiplier relative to the frozen seed: {inventory_scale:.17e}"
        );
        println!("// temperature_k: 700.0; pressure_pa: 101325.0");
        println!(
            "// This is the prepared symbolic residual/Jacobian problem in physical coordinates."
        );
        println!(
            "// It intentionally excludes phase control, normalization, and facade cascade code."
        );
        println!("// The reduced graph has no missing scalar equality to insert blindly:");
        println!("// phase control is an outer active-set/complementarity condition.");
        println!("// Equalities inside this graph are reaction-affinity rows plus A*n-b rows.");
        println!(
            "// TPD, active topology, and phase add/remove decisions are external inequalities."
        );
        println!(
            "// At 700 K graphite is active when published graphite total > 0 and final minimum TPD"
        );
        println!(
            "// satisfies the phase-control acceptance tolerance (the oracle reports +1.87993e-9);"
        );
        println!(
            "// an inactive phase must have zero published amount and satisfy its inactive-phase"
        );
        println!(
            "// stability tolerance (usually expressed as a non-negative minimum TPD condition)."
        );
        println!("// Production elemental normalization: s = sum(element_totals_b), b_hat = b/s,");
        println!(
            "// n_hat = n/s, y_hat = log(n_hat) = log(n) - log(s), and n = s*n_hat on publication."
        );
        println!(
            "// measured production normalization scale s: {:.17e}",
            normalization.physical_inventory_scale()
        );
        println!(
            "// This witness deliberately keeps raw b and raw y, so RST can inspect conditioning directly."
        );
        println!(
            "// Solver coordinates here are y=log(n), not n; normalized recovery uses y_hat=log(n/s)."
        );
        println!(
            "// Current RST production adapter binds the raw symbolic residual/Jacobian; it does not consume the runner's scaled f/j closures."
        );
        println!(
            "// Legacy backends consume the scaled closures. The common validation gate recomputes both raw and scaled residuals afterward."
        );
        println!(
            "// Ideal TPD(x) = sum_i x_i * (g_reduced[i] + R*T*ln(x_i)); zero-support terms contribute zero."
        );
        println!(
            "// TPD units: J/mol; active accepted iff TPD <= dg_keep; inactive accepted iff TPD >= dg_create."
        );
        println!(
            "// phase_eps = {:.17e} mol; dg_create = {:.17e} J/mol; dg_keep = {:.17e} J/mol",
            phase_manager.phase_eps, dg_create, dg_keep
        );
        println!("// Activate an inactive candidate when minimum TPD < dg_create.");
        println!(
            "// Deactivate an active phase only when total < phase_eps and minimum TPD > dg_keep."
        );
        println!("// jacobian singular values: {singular_values:?}");
        println!(
            "// numerical rank threshold: {rank_tolerance:.17e}; rank: {numerical_rank}/{}",
            variables.len()
        );
        println!(
            "// Rank 17 is not claimed to be physically expected: this direct reduced graph is intentionally a singular safety regression."
        );
        println!(
            "// Phase control does not add an equality row here; it selects active sets and checks TPD inequalities outside Newton."
        );
        println!("// The source facade acceptance contract is NOT solver termination alone:");
        println!("//   future strict physical gate: raw_residual_l2 < 1e-5");
        println!("//   max_abs_element_balance_error <= 1e-6");
        println!("//   all published moles finite and non-negative");
        println!(
            "//   elemental b, component/phase ordering, and active topology match the molecular oracle"
        );
        println!("//   component relative error <= 1e-5 (or absolute error <= 1e-10)");
        println!("//   mole-fraction error <= 1e-8 and frozen external regression envelope passes");
        println!("let equations: Vec<String> = vec![");
        for equation in &equations {
            println!("    {equation:?}.to_string(),");
        }
        println!("];");
        println!(
            "let variables: Vec<String> = vec![{}];",
            variables
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let equation_parameters: Vec<String> = vec![{}];",
            parameters
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let equation_parameter_values: Vec<f64> = vec![{}];",
            parameter_values
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let element_labels: Vec<String> = vec![{}];",
            element_labels
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let species_order: Vec<String> = vec![{}];",
            prepared
                .problem()
                .species()
                .iter()
                .map(|value| format!("{value:?}.to_string()"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!("let element_composition_a_element_by_species: Vec<Vec<f64>> = vec![");
        for element in 0..element_matrix.ncols() {
            println!(
                "    vec![{}], // A[{}={}]",
                (0..element_matrix.nrows())
                    .map(|species| format!("{:.17e}", element_matrix[(species, element)]))
                    .collect::<Vec<_>>()
                    .join(", "),
                element,
                element_labels[element]
            );
        }
        println!("];");
        println!("// A is printed as A[element, species]; the library stores species-by-element.");
        println!("// Exact balance definition: balance[e] = sum_i A[e,i] * n[i] - b[e].");
        println!(
            "let element_totals_b: Vec<f64> = vec![{}];",
            element_totals_b
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let species_phase: Vec<usize> = vec![{}];",
            prepared
                .species_phase()
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!("// phase order and species membership:");
        for (index, phase) in prepared.problem().phases().iter().enumerate() {
            println!(
                "// phase[{index}] kind={:?} species_indices={:?}",
                phase.kind, phase.species
            );
        }
        println!("// phase descriptors:");
        for phase in prepared.problem().phase_descriptors() {
            println!(
                "// phase index={} id={:?} component_range={:?}",
                phase.index().index(),
                phase.id(),
                phase.component_range()
            );
        }
        println!(
            "let residual_row_scale_w: Vec<f64> = vec![{}];",
            row_scale
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "// Despite the historical variable name, residual_row_scale_w is the denominator vector D, not diag(W)."
        );
        println!(
            "// W = diag(1 / D); F_raw(y) = [reaction-affinity rows; A*n-b]; F_scaled(y) = W * F_raw(y)."
        );
        println!(
            "// Therefore F_scaled[i] = F_raw[i] / D[i] and J_scaled[i,:] = J_raw[i,:] / D[i]."
        );
        println!(
            "// Reaction D rows are ||nu_k||; elemental D rows are abs(b[e]), or 1 when b[e] == 0."
        );
        println!(
            "// Scaling is numerical conditioning metadata and never replaces the physical A*n-b check."
        );
        println!(
            "let initial_log_moles: Vec<f64> = vec![{}];",
            initial
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "let log_mole_bounds: Option<Vec<(f64, f64)>> = {};",
            bounds.as_ref().map_or_else(
                || "None".to_string(),
                |values| {
                    format!(
                        "Some(vec![{}])",
                        values
                            .iter()
                            .map(|(lower, upper)| { format!("({lower:.17e}, {upper:.17e})") })
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                },
            )
        );
        println!(
            "let initial_moles: Vec<f64> = vec![{}];",
            prepared
                .problem()
                .initial_moles()
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!("variables_count = {};", variables.len());
        println!("equations_count = {};", equations.len());
        println!(
            "// Molecular oracle: production molecular-input facade solve at physical scale 1e8."
        );
        println!(
            "// It is independent of the direct standalone RST callback and uses an explicit legacy-TR facade policy."
        );
        println!(
            "// oracle accepted backend: {:?}",
            oracle_solution.solve_report().accepted_backend
        );
        println!(
            "// oracle component totals / published moles (scale 1e8, exact species order): {oracle_published_moles:?}"
        );
        println!(
            "// oracle numerical moles (positive solver coordinates, scale 1e8): {oracle_numerical_moles:?}"
        );
        println!("// oracle mole fractions: {oracle_mole_fractions:?}");
        println!(
            "// oracle active topology: {}; phase statuses: {:?}",
            active_phase_topology(oracle_solution),
            oracle_solution
                .phases()
                .iter()
                .map(|phase| (phase.id(), oracle_solution.phase_status(phase.id())))
                .collect::<Vec<_>>()
        );
        println!(
            "// oracle phase totals (published, physical scale 1e8): {:?}",
            oracle_solution
                .phases()
                .iter()
                .map(|phase| (phase.id(), oracle_solution.phase_total(phase.id())))
                .collect::<Vec<_>>()
        );
        println!(
            "// oracle final graphite minimum TPD: {:.17e}",
            graphite_final_tpd(oracle_solution)
        );
        println!(
            "// oracle raw residual_l2 on this graph: {:.17e}",
            l2_norm(&oracle_raw_residual)
        );
        println!(
            "// oracle scaled residual_l2 on this graph: {:.17e}",
            l2_norm(&oracle_scaled_residual)
        );
        println!("// oracle raw A*n-b: {oracle_element_balance:?}");
        println!(
            "// current common-gate residual metric: scaled_residual_l2 < 1e-5 when scaling is enabled; raw_residual_l2 is diagnostic only."
        );
        println!("// future physical acceptance gate: raw_residual_l2 < 1e-5;");
        println!(
            "// max_abs(raw A*n-b) <= 1e-6; component relative error <= 1e-5 or abs <= 1e-10;"
        );
        println!(
            "// mole-fraction error <= 1e-8; finite/non-negative published moles; exact ordering/topology."
        );
        println!(
            "// Diagnostic component relative error uses abs(candidate-oracle) / max(abs(candidate), abs(oracle), 1e-30)."
        );
        println!(
            "// Published inactive-phase amounts may be exact zero; internal exp(log-moles) candidates must be finite and strictly positive."
        );
        println!("// Required future RST gate methods: Minpack LM, Nielsen LM, Trust-Region LM.");
        println!("// Advisory characterization methods: plain LM, Damped Newton, Powell Dogleg.");
        println!(
            "// On this rank-deficient reduced graph, typed SingularJacobian/StepTooSmall is valid characterization, not physical acceptance."
        );
        println!(
            "// Production pipeline: physical build -> raw RST symbolic or scaled legacy backend -> candidate validation (raw+scaled residual, physical A*n-b) -> TPD/complementarity -> restart/publication; normalized recovery is used only after classified numerical failure."
        );

        let contract = RustedSciTheSolveContract::new(1.0e-12, 500)
            .expect("RST diagnostic solve contract must validate");
        let mut solvers = RustedSciTheSolver::recommended_cascade();
        solvers.push(RustedSciTheSolver::PowellDogleg);

        println!("RST method matrix:");
        for solver in solvers {
            let started = Instant::now();
            match solver.solve(&rst_problem, &initial, &contract) {
                Ok(outcome) => {
                    let residual = rst_problem
                        .residual_for_test(&outcome.solution)
                        .expect("RST residual must remain evaluable")
                        .iter()
                        .map(|value| value * value)
                        .sum::<f64>()
                        .sqrt();
                    let jacobian = rst_problem
                        .jacobian_for_test(&outcome.solution)
                        .expect("RST Jacobian must remain evaluable");
                    let max_abs_jacobian = jacobian
                        .iter()
                        .map(|value| value.abs())
                        .fold(0.0_f64, f64::max);
                    let candidate_physical = prepared.reconstruct_moles(&outcome.solution);
                    let (
                        finite_non_negative,
                        max_balance,
                        max_component_relative,
                        max_fraction_error,
                    ) = match candidate_physical {
                        Ok(candidate_moles) => {
                            let finite_non_negative = candidate_moles
                                .iter()
                                .all(|value| value.is_finite() && *value >= 0.0);
                            let max_balance = (0..element_matrix.ncols())
                                .map(|element| {
                                    (0..element_matrix.nrows())
                                        .map(|species| {
                                            element_matrix[(species, element)]
                                                * candidate_moles[species]
                                        })
                                        .sum::<f64>()
                                        - element_totals_b[element]
                                })
                                .map(f64::abs)
                                .fold(0.0_f64, f64::max);
                            let max_component_relative = candidate_moles
                                .iter()
                                .zip(oracle_numerical_moles.iter())
                                .map(|(&candidate, &oracle)| {
                                    (candidate - oracle).abs()
                                        / candidate.abs().max(oracle.abs()).max(1.0e-30)
                                })
                                .fold(0.0_f64, f64::max);
                            let mut max_fraction_error = 0.0_f64;
                            for phase in prepared.problem().phases() {
                                let total = phase
                                    .species
                                    .iter()
                                    .map(|&species| candidate_moles[species])
                                    .sum::<f64>();
                                if total.is_finite() && total > 0.0 {
                                    for &species in &phase.species {
                                        max_fraction_error = max_fraction_error.max(
                                            (candidate_moles[species] / total
                                                - oracle_mole_fractions[species])
                                                .abs(),
                                        );
                                    }
                                }
                            }
                            (
                                finite_non_negative,
                                max_balance,
                                max_component_relative,
                                max_fraction_error,
                            )
                        }
                        Err(_) => (false, f64::INFINITY, f64::INFINITY, f64::INFINITY),
                    };
                    println!(
                        "  {} status=returned termination={:?} converged={} residual_l2={residual:.17e} max_abs_jacobian={max_abs_jacobian:.17e} finite_non_negative={} max_abs_element_balance={max_balance:.17e} max_component_relative={max_component_relative:.17e} max_mole_fraction_error={max_fraction_error:.17e} iterations={} residual_calls={} jacobian_calls={} linear_solves={} elapsed_ms={} wall_ms={}",
                        solver.name(),
                        outcome.metrics.termination,
                        outcome.metrics.backend_converged,
                        finite_non_negative,
                        outcome.metrics.iterations,
                        outcome.metrics.residual_evaluations,
                        outcome.metrics.jacobian_evaluations,
                        outcome.metrics.linear_solves,
                        outcome.metrics.elapsed_millis,
                        started.elapsed().as_millis(),
                    );
                }
                Err(error) => {
                    println!(
                        "  {} status=error error={error:?} wall_ms={}",
                        solver.name(),
                        started.elapsed().as_millis(),
                    );
                }
            }
        }
    }

    #[test]
    fn tp1907_upper_inventory_analytic_and_lambdified_symbolic_contracts_agree() {
        let fixture = fixture();
        let prepared = prepared_fixed_tp1907_problem_at_scale(&fixture, 700.0, 1.0e8);
        let rst_problem = prepare_rst_symbolic_problem_from_prepared(&prepared)
            .expect("TP-1907 upper-scale RST problem must prepare");
        let initial = prepared.problem().initial_log_moles().as_slice().to_vec();
        let probes = [
            initial.clone(),
            initial
                .iter()
                .enumerate()
                .map(|(index, value)| value + if index % 2 == 0 { 1.0e-3 } else { -1.0e-3 })
                .collect::<Vec<_>>(),
        ];

        let mut max_residual_abs_error = 0.0_f64;
        let mut max_jacobian_abs_error = 0.0_f64;
        let mut max_residual_relative_error = 0.0_f64;
        let mut max_jacobian_relative_error = 0.0_f64;
        for probe in probes {
            let analytic_residual = prepared.residual(&probe).expect("analytic residual");
            let symbolic_residual = rst_problem
                .residual_for_test(&probe)
                .expect("lambdified symbolic residual");
            assert_eq!(analytic_residual.len(), symbolic_residual.len());
            for (row, (&analytic, symbolic)) in analytic_residual
                .iter()
                .zip(symbolic_residual.iter())
                .enumerate()
            {
                let absolute = (analytic - symbolic).abs();
                let relative = absolute / analytic.abs().max(symbolic.abs()).max(1.0);
                max_residual_abs_error = max_residual_abs_error.max(absolute);
                max_residual_relative_error = max_residual_relative_error.max(relative);
                assert!(
                    relative <= 1.0e-11,
                    "TP-1907 residual row {row} differs: analytic={analytic:e}, symbolic={symbolic:e}, relative={relative:e}"
                );
            }

            let analytic_jacobian = prepared.jacobian(&probe).expect("analytic Jacobian");
            let symbolic_jacobian = rst_problem
                .jacobian_for_test(&probe)
                .expect("lambdified symbolic Jacobian");
            assert_eq!(analytic_jacobian.shape(), symbolic_jacobian.shape());
            for (index, (&analytic, symbolic)) in analytic_jacobian
                .iter()
                .zip(symbolic_jacobian.iter())
                .enumerate()
            {
                let absolute = (analytic - symbolic).abs();
                let relative = absolute / analytic.abs().max(symbolic.abs()).max(1.0);
                max_jacobian_abs_error = max_jacobian_abs_error.max(absolute);
                max_jacobian_relative_error = max_jacobian_relative_error.max(relative);
                assert!(
                    relative <= 1.0e-11,
                    "TP-1907 Jacobian entry {index} differs: analytic={analytic:e}, symbolic={symbolic:e}, relative={relative:e}"
                );
            }
        }

        println!(
            "TP-1907 analytic-vs-symbolic upper scale max_residual_abs={max_residual_abs_error:e} max_residual_relative={max_residual_relative_error:e} max_jacobian_abs={max_jacobian_abs_error:e} max_jacobian_relative={max_jacobian_relative_error:e}"
        );
    }

    #[test]
    fn tp1907_upper_inventory_is_backend_invariant() {
        let fixture = fixture();
        let repository = repository();
        let frozen_dataset = dataset();
        let reference = frozen_dataset
            .rows()
            .iter()
            .find(|row| row.temperature_k == 700.0)
            .expect("TP-1907 frozen P,T row at 700 K must remain available");
        let initial = scaled_initial_composition(&fixture, 1.0e8);
        let (inventory, element_labels, reference_b) = frozen_element_inventory(&fixture, &initial);
        let molecular = calculator_tp1907_pt(
            &fixture,
            repository.clone(),
            reference.temperature_k,
            Some(initial.moles().to_vec()),
            None,
        );
        let backends = [
            (
                "rst-lm",
                SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
            ),
            (
                "rst-backtracking-lm",
                SolverBackend::RustedSciThe(RustedSciTheSolver::BacktrackingLevenbergMarquardt),
            ),
            (
                "rst-minpack-lm",
                SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
            ),
            (
                "rst-nielsen-lm",
                SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
            ),
            (
                "rst-trust-region-lm",
                SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
            ),
            (
                "rst-damped-newton",
                SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
            ),
            ("legacy-lm", SolverBackend::Legacy(Solvers::LM)),
            ("legacy-nr", SolverBackend::Legacy(Solvers::NR)),
            ("legacy-tr", SolverBackend::Legacy(Solvers::TR)),
        ];
        let mut failures = Vec::new();
        let mut successes = Vec::new();

        for (name, backend) in backends {
            let options = EquilibriumSolveOptions::default()
                .with_max_iterations(500)
                .expect("backend invariant iteration budget must validate")
                .with_tolerance(1.0e-12)
                .expect("backend invariant tolerance must validate");
            let outcome = match calculator_tp1907_builder(&fixture, repository.clone())
                .solve_options(options)
                .solver_policy(SolverPolicy::Single(backend))
                .expect("single backend policy must validate")
                .at_temperature(reference.temperature_k)
                .element_inventory(inventory.clone())
                .element_numerical_seed(initial.clone())
                .solve()
            {
                Ok(outcome) => outcome,
                Err(error) => {
                    println!("TP-1907 backend={name} scale=1e8 FAILED error={error:?}",);
                    failures.push(name);
                    continue;
                }
            };
            let point = match outcome {
                EquilibriumCalculatorOutcome::PtPoint(point) => point,
                other => panic!("{name} returned unexpected outcome: {other:?}"),
            };
            let solution = point.solution();
            assert_eq!(
                solution.build_report().element_labels(),
                &element_labels,
                "{name} changed the canonical element order"
            );
            for (index, (&expected, &actual)) in solution
                .build_report()
                .element_totals()
                .iter()
                .zip(reference_b.iter())
                .enumerate()
            {
                assert!(
                    (expected - actual).abs() <= 1.0e-10 * expected.abs().max(1.0),
                    "{name} changed b at element {index}: actual={actual:e}, expected={expected:e}"
                );
            }
            assert_eq!(
                active_phase_topology(solution),
                active_phase_topology(molecular.solution()),
                "{name} changed the accepted phase topology"
            );
            assert!(
                solution
                    .component_moles()
                    .iter()
                    .all(|amount| amount.is_finite() && *amount >= 0.0),
                "{name} published invalid physical moles"
            );
            let validation = solution.accepted_solution().validation();
            assert!(
                validation.residual_l2_norm < 1.0e-5,
                "{name} residual={:e}",
                validation.residual_l2_norm
            );
            assert!(
                validation.max_abs_element_balance_error
                    <= 1.0e-12 * reference_b.iter().sum::<f64>(),
                "{name} relative balance={:e}",
                validation.max_abs_element_balance_error / reference_b.iter().sum::<f64>()
            );
            for (index, (&expected, &actual)) in molecular
                .solution()
                .component_moles()
                .iter()
                .zip(solution.component_moles())
                .enumerate()
            {
                let relative =
                    (expected - actual).abs() / expected.abs().max(actual.abs()).max(1.0e-12);
                assert!(
                    relative <= 1.0e-5 || (expected - actual).abs() <= 1.0e-10,
                    "{name} component {index} differs from molecular oracle: relative={relative:e}"
                );
            }
            println!(
                "TP-1907 backend={name} scale=1e8 residual={:.3e} balance={:.3e}",
                validation.residual_l2_norm, validation.max_abs_element_balance_error,
            );
            successes.push(name);
        }
        assert_eq!(
            successes,
            vec!["rst-trust-region-lm", "legacy-lm", "legacy-nr", "legacy-tr",],
            "accepted backend matrix changed at inventory scale 1e8"
        );
        assert_eq!(
            failures,
            vec![
                "rst-lm",
                "rst-backtracking-lm",
                "rst-minpack-lm",
                "rst-nielsen-lm",
                "rst-damped-newton",
            ],
            "isolated backend failure matrix changed at inventory scale 1e8"
        );
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
        let mut runner = PreparedEquilibriumRunner::new(prepared_fixed_tp1907_problem_at_scale(
            &fixture, 700.0, factor,
        ))
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
