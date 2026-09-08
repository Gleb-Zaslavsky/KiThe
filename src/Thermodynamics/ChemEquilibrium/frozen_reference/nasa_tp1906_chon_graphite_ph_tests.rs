//! I5-H/I5-composition lifecycle evidence for NASA TP-1906/1907 CHON graphite.

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TotalEnthalpyJoules;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_extensive_normalization::ExtensiveNormalization;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
        PhEnthalpyGrid, PhRangePointPreparation, PhRangeRequest,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
        PhSolveMode, PhTemperatureSolveOptions, solve_resolved_ph,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::TraceSpeciesSeedPolicy;
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::assertions::{
        AcceptedSolutionContract, assert_bounded_solution_accepted, assert_byte_snapshots_unchanged,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph::{
        ResolvedNasaTp1906ChonGraphitePhFixture, load_nasa_tp1906_chon_graphite_enthalpy_dataset,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        EquilibriumSolveOptions, ExtensiveNormalizationPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::prelude::PhaseControlPolicy;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;

    fn repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository().expect("bundled local repository must be available")
    }

    fn fixture() -> ResolvedNasaTp1906ChonGraphitePhFixture {
        ResolvedNasaTp1906ChonGraphitePhFixture::resolve_offline(repository())
            .expect("reviewed TP-1906/1907 local fixture must resolve offline")
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

    fn frozen_snapshot() -> Vec<(String, Vec<u8>)> {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data");
        [
            root.join("nasa_tp1906/chon_graphite_er125_1atm_heterogeneous_enthalpy.metadata.json"),
            root.join("nasa_tp1906/chon_graphite_er125_1atm_heterogeneous_enthalpy.rows.json"),
            root.join("nasa_tp1907/chon_graphite_er125_1atm.metadata.json"),
            root.join("nasa_tp1907/chon_graphite_er125_1atm.rows.json"),
        ]
        .into_iter()
        .map(|path| {
            let label = path.display().to_string();
            (
                label,
                fs::read(path).expect("frozen source must be readable"),
            )
        })
        .collect()
    }

    fn ph_pt_witness_max_relative_difference(
        fixture: &ResolvedNasaTp1906ChonGraphitePhFixture,
        ph: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::FixedPressureEnthalpySolution,
    ) -> f64 {
        let witness = fixture.solve_pt_witness(ph.temperature()).unwrap();
        assert_eq!(
            fixture.graphite_is_active(ph.equilibrium()),
            fixture.graphite_is_active(&witness),
            "P,H must land on the canonical P,T phase manifold at recovered temperature"
        );
        ph.equilibrium()
            .component_moles()
            .iter()
            .zip(witness.component_moles())
            .map(|(&ph_moles, &pt_moles)| {
                (ph_moles - pt_moles).abs() / ph_moles.abs().max(pt_moles.abs()).max(1.0e-12)
            })
            .fold(0.0_f64, f64::max)
    }

    fn print_ph_pt_witness_delta(
        fixture: &ResolvedNasaTp1906ChonGraphitePhFixture,
        ph: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::FixedPressureEnthalpySolution,
    ) {
        let witness = fixture.solve_pt_witness(ph.temperature()).unwrap();
        let mut rows = ph
            .equilibrium()
            .component_moles()
            .iter()
            .zip(witness.component_moles())
            .enumerate()
            .map(|(index, (&ph_moles, &pt_moles))| {
                let absolute = (ph_moles - pt_moles).abs();
                let relative = absolute / ph_moles.abs().max(pt_moles.abs()).max(1.0e-12);
                (index, ph_moles, pt_moles, absolute, relative)
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| right.3.partial_cmp(&left.3).unwrap());
        println!(
            "    witness: T={:.6} residual={:.3e} balance={:.3e} graphite P,H/P,T={:.6e}/{:.6e}",
            ph.temperature(),
            ph.equilibrium().accepted_solution().validation().residual_l2_norm,
            ph.equilibrium().accepted_solution().validation().max_abs_element_balance_error,
            ph.equilibrium()
                .moles_for(&crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite::graphite_component())
                .unwrap(),
            witness
                .moles_for(&crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite::graphite_component())
                .unwrap(),
        );
        for (index, ph_moles, pt_moles, absolute, relative) in rows.into_iter().take(4) {
            println!(
                "      component[{index}] P,H={ph_moles:.6e} P,T={pt_moles:.6e} abs={absolute:.3e} rel={relative:.3e}"
            );
        }
    }

    fn assert_ph_matches_pt_witness(
        fixture: &ResolvedNasaTp1906ChonGraphitePhFixture,
        ph: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::FixedPressureEnthalpySolution,
    ) {
        let maximum = ph_pt_witness_max_relative_difference(fixture, ph);
        assert!(
            maximum <= 1.0e-6,
            "P,H/P,T witness maximum relative mole mismatch is {maximum:e}"
        );
    }

    fn graphite_history(
        fixture: &ResolvedNasaTp1906ChonGraphitePhFixture,
        range: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::PhRangeSolution,
    ) -> Vec<bool> {
        range
            .points()
            .iter()
            .map(|point| fixture.graphite_is_active(point.solution().equilibrium()))
            .collect()
    }

    fn assert_forward_reverse_agree(
        forward: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::PhRangeSolution,
        reverse: &crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::PhRangeSolution,
    ) {
        let reverse_by_target = reverse
            .points()
            .iter()
            .map(|point| (point.report().target_enthalpy_joules().to_bits(), point))
            .collect::<BTreeMap<_, _>>();
        for forward_point in forward.points() {
            let reverse_point = reverse_by_target
                .get(&forward_point.report().target_enthalpy_joules().to_bits())
                .expect("both continuation directions must publish the same targets");
            assert!(
                (forward_point.solution().temperature() - reverse_point.solution().temperature())
                    .abs()
                    <= 1.0e-3,
                "forward/reverse temperature mismatch at target {:.6e} J: forward={:.12} K, reverse={:.12} K, delta={:.6e} K",
                forward_point.report().target_enthalpy_joules(),
                forward_point.solution().temperature(),
                reverse_point.solution().temperature(),
                forward_point.solution().temperature() - reverse_point.solution().temperature(),
            );
            for (index, (&forward_moles, &reverse_moles)) in forward_point
                .solution()
                .equilibrium()
                .component_moles()
                .iter()
                .zip(reverse_point.solution().equilibrium().component_moles())
                .enumerate()
            {
                let scale = forward_moles.abs().max(reverse_moles.abs()).max(1.0e-12);
                assert!(
                    (forward_moles - reverse_moles).abs() / scale <= 1.0e-4
                        || (forward_moles - reverse_moles).abs() <= 1.0e-7,
                    "forward/reverse mole mismatch at component {index}: forward={forward_moles:e}, reverse={reverse_moles:e}, relative={:e}",
                    (forward_moles - reverse_moles).abs() / scale,
                );
            }
        }
    }

    fn scaled_initial_composition(
        fixture: &ResolvedNasaTp1906ChonGraphitePhFixture,
        factor: f64,
    ) -> MultiphaseInitialComposition {
        let layout =
            MultiphaseEquilibriumLayout::new(fixture.tp1907().resolved().phase_specs().to_vec())
                .expect("TP-1906/1907 phase layout must remain valid");
        let moles = fixture
            .tp1907()
            .initial_composition()
            .expect("TP-1906/1907 inventory must validate")
            .moles()
            .iter()
            .map(|moles| factor * moles)
            .collect();
        MultiphaseInitialComposition::from_dense(&layout, moles)
            .expect("scaled TP-1906/1907 inventory must validate")
    }

    /// Test-only exact extensive normalization for a P,H request.
    ///
    /// It derives the factor from the request's current physical inventory,
    /// not from a frozen reference state. The caller must apply the same
    /// factor to its total enthalpy target before invoking the ordinary nested
    /// P,H workflow.
    fn normalized_initial_composition(
        fixture: &ResolvedNasaTp1906ChonGraphitePhFixture,
        physical: &MultiphaseInitialComposition,
    ) -> (MultiphaseInitialComposition, ExtensiveNormalization) {
        let normalization = ExtensiveNormalization::from_physical_moles(physical.moles())
            .expect("TP-1906/1907 physical inventory must define extensive normalization");
        let layout =
            MultiphaseEquilibriumLayout::new(fixture.tp1907().resolved().phase_specs().to_vec())
                .expect("TP-1906/1907 normalized layout must remain valid");
        let normalized = normalization
            .normalize_initial_composition(&layout, physical)
            .expect("TP-1906/1907 normalized inventory must validate");
        (normalized, normalization)
    }

    /// The canonical production P,H tolerances scale with extensive enthalpy.
    /// This test instead needs a tighter, reproducible scalar contract, while
    /// retaining the standard backend cascade and bounded phase lifecycle.
    fn extensivity_temperature_options() -> PhTemperatureSolveOptions {
        let mut options = PhTemperatureSolveOptions::default();
        options.scaled_enthalpy_tolerance = 1.0e-9;
        options.absolute_enthalpy_tolerance_joules = 1.0e-9;
        options
    }

    fn extensivity_solve_options() -> EquilibriumSolveOptions {
        EquilibriumSolveOptions::default()
            .with_max_iterations(100)
            .expect("extensivity iteration budget must validate")
            .with_tolerance(1.0e-9)
            .expect("extensivity tolerance must validate")
    }

    // The P,H outer root is assembled from independently accepted P,T trials.
    // The controlled two-order inventory probe therefore permits a few
    // microkelvin of numerical movement while still rejecting physical scale
    // dependence. It is deliberately separate from the mole-scaling contract.
    const PH_EXTENSIVITY_TEMPERATURE_TOLERANCE_K: f64 = 1.0e-5;

    // The most temperature-sensitive trace component remains materially below
    // this envelope under the controlled two-order inventory probe. The bound
    // remains two orders tighter than the continuation comparison and detects
    // any physical scale dependence.
    const PH_EXTENSIVITY_MOLE_RELATIVE_TOLERANCE: f64 = 1.0e-4;

    #[test]
    fn tp1906_rows_keep_specific_enthalpy_provenance_and_semantic_tp1907_joins() {
        let before_library = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let dataset = load_nasa_tp1906_chon_graphite_enthalpy_dataset().unwrap();
        assert!(dataset.is_external_evidence());
        assert_eq!(
            dataset.metadata().dataset_id,
            "nasa_tp1906.chon_graphite.er125.1atm.heterogeneous_enthalpy.v1"
        );
        assert_eq!(
            dataset
                .metadata()
                .quantities
                .iter()
                .find(|quantity| quantity.column == "specific_enthalpy_j_g")
                .expect("specific enthalpy metadata is required")
                .unit,
            "J/g"
        );
        assert_eq!(
            dataset
                .rows()
                .iter()
                .map(|row| row.specific_enthalpy_j_g)
                .collect::<Vec<_>>(),
            vec![-2375.5, -2334.1, -2291.4, -2247.1]
        );

        let fixture = fixture();
        let cases = fixture.cases().unwrap();
        assert_eq!(cases.len(), 4);
        assert!(fixture.inventory_mass_g().is_finite() && fixture.inventory_mass_g() > 0.0);
        for case in &cases {
            assert_eq!(case.enthalpy.temperature_k, case.composition.temperature_k);
            assert_eq!(case.enthalpy.pressure_pa, case.composition.pressure_pa);
            assert!(
                (case.target_enthalpy_j
                    - case.source_specific_enthalpy_j_g() * case.inventory_mass_g)
                    .abs()
                    <= 1.0e-12
            );
        }
        assert_byte_snapshots_unchanged(&before_library, &local_library_snapshot());
        assert_byte_snapshots_unchanged(&before_frozen, &frozen_snapshot());
    }

    #[test]
    fn monolithic_h700_lands_on_the_canonical_pt_graphite_branch() {
        let fixture = fixture();
        let case = fixture
            .cases()
            .unwrap()
            .into_iter()
            .find(|case| (case.enthalpy.temperature_k - 700.0).abs() <= f64::EPSILON)
            .expect("the frozen TP-1906 dataset must contain the 700 K row");
        let solution = solve_resolved_ph(
            fixture
                .ph_request(&case, 900.0)
                .unwrap()
                .with_ph_solve_mode(PhSolveMode::Monolithic),
        )
        .expect("monolithic P,H must solve the near-boundary graphite case");

        assert_eq!(
            solution.report().solve_path(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolvePath::MonolithicPhaseControl,
        );
        assert!(fixture.graphite_is_active(solution.equilibrium()));
        assert_ph_matches_pt_witness(&fixture, &solution);
        assert!(
            solution
                .equilibrium()
                .accepted_solution()
                .validation()
                .max_abs_reaction_affinity
                <= 1.0e-6,
            "accepted monolithic state must satisfy the dimensionless affinity contract"
        );
    }

    #[test]
    fn tp1906_nested_monolithic_and_auto_agree_on_h700_physical_state() {
        let fixture = fixture();
        let case = fixture
            .cases()
            .unwrap()
            .into_iter()
            .find(|case| (case.enthalpy.temperature_k - 700.0).abs() <= f64::EPSILON)
            .expect("the frozen TP-1906 dataset must contain the 700 K row");
        let solve = |mode| {
            solve_resolved_ph(
                fixture
                    .ph_request(&case, 900.0)
                    .unwrap()
                    .with_ph_solve_mode(mode),
            )
            .expect("TP-1906 route-matrix solve must accept")
        };
        let nested = solve(PhSolveMode::NestedTemperature);
        let monolithic = solve(PhSolveMode::Monolithic);
        let automatic = solve(PhSolveMode::Auto);
        let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);
        for solution in [&nested, &monolithic, &automatic] {
            assert_bounded_solution_accepted(solution.equilibrium(), contract);
            assert_ph_matches_pt_witness(&fixture, solution);
        }
        assert_eq!(
            nested.report().solve_path(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolvePath::NestedTemperature
        );
        assert_eq!(
            monolithic.report().solve_path(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolvePath::MonolithicPhaseControl
        );
        assert_eq!(
            automatic.report().solve_path(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolvePath::MonolithicPhaseControl,
            "Auto must retain the direct accepted monolithic route instead of an unnecessary fallback"
        );

        for (label, candidate) in [("monolithic", &monolithic), ("auto", &automatic)] {
            assert!(
                (candidate.temperature() - nested.temperature()).abs() <= 1.0e-5,
                "P,H route {label} changed the intensive temperature: nested={:.12} K, candidate={:.12} K",
                nested.temperature(),
                candidate.temperature(),
            );
            for (index, (&nested_moles, &candidate_moles)) in nested
                .equilibrium()
                .component_moles()
                .iter()
                .zip(candidate.equilibrium().component_moles())
                .enumerate()
            {
                let relative = (nested_moles - candidate_moles).abs()
                    / nested_moles.abs().max(candidate_moles.abs()).max(1.0e-12);
                assert!(
                    relative <= 1.0e-5,
                    "P,H route {label} differs at component {index}: nested={nested_moles:e}, candidate={candidate_moles:e}, relative={relative:e}"
                );
            }
        }
    }

    #[test]
    fn tp1906_four_point_ph_lifecycle_preserves_pt_manifold_and_graphite_topology() {
        let before_library = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let fixture = fixture();
        let cases = fixture.cases().unwrap();

        for case in &cases {
            let solution = solve_resolved_ph(
                fixture
                    .ph_request(case, 900.0)
                    .unwrap()
                    .with_ph_solve_mode(PhSolveMode::NestedTemperature),
            )
            .expect("each frozen TP-1906 P,H point must solve through the canonical route");
            assert!(solution.scaled_enthalpy_error().is_finite());
            assert!(
                solution.enthalpy_error().abs() <= solution.enthalpy_error_limit_joules(),
                "TP-1906 P,H point at {} K must satisfy its enthalpy acceptance contract",
                case.enthalpy.temperature_k,
            );
            assert_eq!(
                fixture.graphite_is_active(solution.equilibrium()),
                case.expected_graphite_active(),
                "TP-1906 P,H point at {} K must preserve TP-1907 graphite topology",
                case.enthalpy.temperature_k,
            );
            assert_ph_matches_pt_witness(&fixture, &solution);
        }

        let solve_range = |targets: Vec<f64>| {
            PhRangeRequest::from_resolved_thermochemistry(
                fixture.tp1907().resolved(),
                fixture.tp1907().initial_composition().unwrap(),
                101_325.0,
                101_325.0,
                PhEnthalpyGrid::new(targets).unwrap(),
                fixture.tp1907().thermochemistry().temperature_bounds(),
                900.0,
                fixture.tp1907().thermochemistry().clone(),
            )
            .unwrap()
            .with_phase_control_policy(PhaseControlPolicy::default())
            .with_ph_solve_mode(PhSolveMode::NestedTemperature)
            .solve()
            .expect("frozen TP-1906 P,H continuation must solve")
        };
        let forward = solve_range(cases.iter().map(|case| case.target_enthalpy_j).collect());
        let reverse = solve_range(
            cases
                .iter()
                .rev()
                .map(|case| case.target_enthalpy_j)
                .collect(),
        );
        let expected_forward = cases
            .iter()
            .map(|case| case.expected_graphite_active())
            .collect::<Vec<_>>();
        assert_eq!(graphite_history(&fixture, &forward), expected_forward);
        assert_eq!(
            graphite_history(&fixture, &reverse),
            expected_forward.iter().rev().copied().collect::<Vec<_>>(),
        );
        assert_eq!(
            forward.points()[0].report().preparation(),
            PhRangePointPreparation::Initial
        );
        assert!(
            forward.points()[1..]
                .iter()
                .all(|point| point.report().preparation() == PhRangePointPreparation::Continued)
        );
        assert_eq!(
            reverse.points()[0].report().preparation(),
            PhRangePointPreparation::Initial
        );
        assert!(
            reverse.points()[1..]
                .iter()
                .all(|point| point.report().preparation() == PhRangePointPreparation::Continued)
        );
        assert_forward_reverse_agree(&forward, &reverse);

        assert_eq!(before_library, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    #[test]
    fn tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together() {
        let fixture = fixture();
        let case = fixture
            .cases()
            .unwrap()
            .into_iter()
            .find(|case| (case.enthalpy.temperature_k - 700.0).abs() <= f64::EPSILON)
            .expect("TP-1906 fixture must retain the 700 K graphite row");
        let baseline = solve_resolved_ph(
            fixture
                .ph_request(&case, 900.0)
                .unwrap()
                .with_ph_solve_mode(PhSolveMode::NestedTemperature)
                .with_solve_options(extensivity_solve_options())
                .with_temperature_options(extensivity_temperature_options())
                .expect("TP-1906 extensivity scalar controls must validate"),
        )
        .expect("baseline TP-1906 P,H solve must accept");
        let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);
        assert_bounded_solution_accepted(baseline.equilibrium(), contract);

        for factor in [1.0e-1, 1.0e1] {
            let request = fixture
                .ph_request(&case, 900.0)
                .unwrap()
                .with_initial_composition(scaled_initial_composition(&fixture, factor))
                .expect("scaled TP-1906 P,H composition must validate")
                .with_target_enthalpy_and_seed(
                    TotalEnthalpyJoules::new(factor * case.target_enthalpy_j)
                        .expect("scaled TP-1906 enthalpy target must be finite"),
                    900.0,
                )
                .expect("scaled TP-1906 P,H target must validate")
                .with_ph_solve_mode(PhSolveMode::NestedTemperature)
                .with_solve_options(extensivity_solve_options())
                .with_temperature_options(extensivity_temperature_options())
                .expect("TP-1906 extensivity scalar controls must validate");
            let candidate =
                solve_resolved_ph(request).expect("scaled TP-1906 P,H solve must accept");
            assert_bounded_solution_accepted(candidate.equilibrium(), contract);
            assert!(
                (candidate.temperature() - baseline.temperature()).abs()
                    <= PH_EXTENSIVITY_TEMPERATURE_TOLERANCE_K,
                "P,H temperature must be intensive under factor {factor:e}: baseline={:.12} K, candidate={:.12} K, delta={:.6e} K",
                baseline.temperature(),
                candidate.temperature(),
                candidate.temperature() - baseline.temperature(),
            );
            assert_eq!(
                fixture.graphite_is_active(candidate.equilibrium()),
                fixture.graphite_is_active(baseline.equilibrium()),
                "P,H graphite topology must be intensive under factor {factor:e}"
            );
            for (index, (&expected, &actual)) in baseline
                .equilibrium()
                .component_moles()
                .iter()
                .zip(candidate.equilibrium().component_moles())
                .enumerate()
            {
                let recovered = actual / factor;
                let relative =
                    (expected - recovered).abs() / expected.abs().max(recovered.abs()).max(1.0e-12);
                assert!(
                    relative <= PH_EXTENSIVITY_MOLE_RELATIVE_TOLERANCE,
                    "TP-1906 extensive P,H mismatch at component {index}, factor={factor:e}: baseline={expected:e}, recovered={recovered:e}, relative={relative:e}"
                );
            }
        }
    }

    #[test]
    #[ignore = "release characterization of real TP-1906 P,H scale and threshold bands"]
    fn i5_tp1906_extensive_ph_scale_matrix_characterization() {
        const TRACE_FLOOR: f64 = 1.0e-30;
        const DEFAULT_PHASE_EPS: f64 = 1.0e-30;

        let fixture = fixture();
        let cases = fixture.cases().unwrap();
        println!("NASA TP-1906/1907 CHON+graphite P,H extensive-scale characterization");
        println!(
            "trace_floor={TRACE_FLOOR:e} mol phase_eps={DEFAULT_PHASE_EPS:e} mol; fresh nested P,H solves, no continuation"
        );
        println!(
            "  T source K  factor     T KiThe    delta T K   graphite  max n/f err  H scaled    route                  conditioning"
        );

        for temperature_k in [700.0, 720.0] {
            let case = cases
                .iter()
                .find(|case| (case.enthalpy.temperature_k - temperature_k).abs() <= f64::EPSILON)
                .expect("TP-1906 scale matrix must retain the requested source row");
            let solve = |factor| {
                let request = fixture
                    .ph_request(case, 900.0)
                    .expect("TP-1906 P,H scale-matrix request must validate")
                    .with_initial_composition(scaled_initial_composition(&fixture, factor))
                    .expect("scaled TP-1906 P,H inventory must validate")
                    .with_target_enthalpy_and_seed(
                        TotalEnthalpyJoules::new(factor * case.target_enthalpy_j)
                            .expect("scaled TP-1906 P,H target must be finite"),
                        900.0,
                    )
                    .expect("scaled TP-1906 P,H target and seed must validate")
                    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
                    .with_solve_options(extensivity_solve_options().with_trace_seed_policy(
                        TraceSpeciesSeedPolicy::Absolute { floor: TRACE_FLOOR },
                    ))
                    .with_temperature_options(extensivity_temperature_options())
                    .expect("TP-1906 P,H scale-matrix controls must validate");
                solve_resolved_ph(request)
            };
            let baseline =
                solve(1.0).expect("unit-scale TP-1906 P,H characterization solve must accept");
            let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);
            assert_bounded_solution_accepted(baseline.equilibrium(), contract);

            for factor in [1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4] {
                let candidate = solve(factor).unwrap_or_else(|error| {
                    panic!(
                        "production nested P,H extensive recovery failed at source {temperature_k} K, factor={factor:e}: {error}"
                    )
                });
                assert_bounded_solution_accepted(candidate.equilibrium(), contract);
                let max_mole_error = baseline
                    .equilibrium()
                    .component_moles()
                    .iter()
                    .zip(candidate.equilibrium().component_moles())
                    .map(|(&expected, &actual)| {
                        let recovered = actual / factor;
                        (expected - recovered).abs()
                            / expected.abs().max(recovered.abs()).max(1.0e-30)
                    })
                    .fold(0.0_f64, f64::max);
                let recovery = candidate.equilibrium().extensive_normalization_recovery();
                let recovery_label = recovery.map_or("direct", |evidence| {
                    if evidence.reconstructed_physical_boundary {
                        "reconstructed"
                    } else {
                        "physical-retry"
                    }
                });
                println!(
                    "  {temperature_k:10.3}  {factor:8.1e}  {temperature:9.4}  {delta_temperature:10.3e}  {graphite:>8}  {max_mole_error:11.3e}  {enthalpy:10.3e}  {route:22?}  {recovery_label}",
                    temperature = candidate.temperature(),
                    delta_temperature = candidate.temperature() - baseline.temperature(),
                    graphite = fixture.graphite_is_active(candidate.equilibrium()),
                    enthalpy = candidate.scaled_enthalpy_error(),
                    route = candidate.report().solve_path(),
                );
                assert!(
                    max_mole_error <= 1.0e-4,
                    "physical P,H component amounts lost extensive covariance at source {temperature_k} K, factor={factor:e}: {max_mole_error:e}"
                );
                assert!(
                    (candidate.temperature() - baseline.temperature()).abs() <= 1.0e-4,
                    "P,H temperature lost intensive invariance at source {temperature_k} K, factor={factor:e}"
                );
                for (element_index, (&baseline_total, &candidate_total)) in baseline
                    .equilibrium()
                    .build_report()
                    .element_totals()
                    .iter()
                    .zip(candidate.equilibrium().build_report().element_totals())
                    .enumerate()
                {
                    let recovered = candidate_total / factor;
                    let relative = (baseline_total - recovered).abs()
                        / baseline_total.abs().max(recovered.abs()).max(1.0e-30);
                    assert!(
                        relative <= 1.0e-12,
                        "physical P,H build-report element total {element_index} was not reconstructed at source {temperature_k} K, factor={factor:e}"
                    );
                }
                if factor == 1.0e4 {
                    let evidence = recovery.expect(
                        "the known large-inventory nested P,H basin must retain recovery provenance",
                    );
                    if temperature_k < 720.0 {
                        assert!(evidence.reconstructed_physical_boundary);
                    } else {
                        assert!(!evidence.reconstructed_physical_boundary);
                        assert!(evidence.physical_retry_backend.is_some());
                    }
                } else {
                    assert!(
                        recovery.is_none(),
                        "ordinary P,H scale factor {factor:e} unexpectedly entered recovery"
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "release diagnostic for exact test-only P,H extensive normalization"]
    fn i5_tp1906_extensive_ph_normalization_recovery_matrix() {
        let fixture = fixture();
        let cases = fixture.cases().unwrap();
        let factor = 1.0e4_f64;
        let contract = AcceptedSolutionContract::new(1.0e-6, 1.0e-6);

        println!("NASA TP-1906/1907 P,H exact extensive-normalization recovery");
        println!(
            "T source K  physical factor  internal total  T normalized  delta T K  graphite  max n/base err  H internal J  H physical J  status"
        );

        for temperature_k in [700.0, 720.0] {
            let case = cases
                .iter()
                .find(|case| (case.enthalpy.temperature_k - temperature_k).abs() <= f64::EPSILON)
                .expect("TP-1906 normalization matrix must retain the source row");
            let build_request = |composition, target_enthalpy| {
                fixture
                    .ph_request(case, 900.0)
                    .expect("TP-1906 P,H normalization request must validate")
                    .with_initial_composition(composition)
                    .expect("P,H normalization composition must validate")
                    .with_target_enthalpy_and_seed(
                        TotalEnthalpyJoules::new(target_enthalpy)
                            .expect("P,H normalization target must be finite"),
                        900.0,
                    )
                    .expect("P,H normalization target and seed must validate")
                    .with_ph_solve_mode(PhSolveMode::NestedTemperature)
                    .with_solve_options(extensivity_solve_options())
                    .with_temperature_options(extensivity_temperature_options())
                    .expect("P,H normalization controls must validate")
            };

            // The base solution is an external comparison target only. It is
            // never used to form the normalized request or its initial seed.
            let base = solve_resolved_ph(build_request(
                scaled_initial_composition(&fixture, 1.0),
                case.target_enthalpy_j,
            ))
            .expect("unit-scale P,H reference solve must accept");
            let physical_composition = scaled_initial_composition(&fixture, factor);
            let physical_target = factor * case.target_enthalpy_j;
            // This is the historical fresh-coordinate witness. Disable the
            // production recovery policy here so the test continues to prove
            // that the unnormalized basin is ill-conditioned, while the
            // normalized request below proves the supported recovery route.
            let fresh_physical = solve_resolved_ph(
                build_request(physical_composition.clone(), physical_target).with_solve_options(
                    extensivity_solve_options().with_extensive_normalization_policy(
                        ExtensiveNormalizationPolicy::Disabled,
                    ),
                ),
            );
            assert!(
                fresh_physical.is_err(),
                "the established factor=1e4 fresh P,H witness must still fail at {temperature_k} K"
            );

            let (normalized_composition, normalization) =
                normalized_initial_composition(&fixture, &physical_composition);
            let normalization_factor = normalization.physical_inventory_scale();
            let internal_total = normalized_composition.moles().iter().sum::<f64>();
            let normalized_target = normalization
                .normalize_total_enthalpy(TotalEnthalpyJoules::new(physical_target).unwrap())
                .expect("P,H target enthalpy must normalize");
            let mut normalized_temperature_options = extensivity_temperature_options();
            normalized_temperature_options.absolute_enthalpy_tolerance_joules = normalization
                .normalize_physical_energy_tolerance(
                    normalized_temperature_options.absolute_enthalpy_tolerance_joules,
                )
                .expect("P,H absolute enthalpy tolerance must normalize");
            let normalized = solve_resolved_ph(
                build_request(normalized_composition, normalized_target.joules())
                    .with_solve_options(
                        extensivity_solve_options().with_trace_seed_policy(
                            normalization
                                .normalize_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute {
                                    floor: 1.0e-30,
                                })
                                .expect("P,H trace floor must normalize"),
                        ),
                    )
                    .with_phase_control_policy(
                        PhaseControlPolicy::default()
                            .normalized_for_extensive_representation(normalization)
                            .expect("P,H phase epsilon must normalize"),
                    )
                    .with_temperature_options(normalized_temperature_options)
                    .expect("P,H normalized controls must validate"),
            )
            .expect("ordinary nested P,H solve of normalized equivalent must accept");
            assert_bounded_solution_accepted(normalized.equilibrium(), contract);
            assert!(
                (internal_total - 1.0).abs() <= 1.0e-12,
                "P,H internal input must be normalized to one mole"
            );
            assert!(
                (normalized.temperature() - base.temperature()).abs()
                    <= PH_EXTENSIVITY_TEMPERATURE_TOLERANCE_K,
                "P,H normalization must preserve recovered temperature at {temperature_k} K: base={} normalized={}",
                base.temperature(),
                normalized.temperature(),
            );
            assert_eq!(
                fixture.graphite_is_active(normalized.equilibrium()),
                fixture.graphite_is_active(base.equilibrium()),
                "P,H normalization must preserve graphite topology at {temperature_k} K"
            );
            let max_mole_error = base
                .equilibrium()
                .component_moles()
                .iter()
                .zip(normalized.equilibrium().component_moles())
                .map(|(&base_moles, &normalized_moles)| {
                    let reconstructed = normalized_moles * normalization_factor;
                    let physical_baseline = base_moles * factor;
                    (physical_baseline - reconstructed).abs()
                        / physical_baseline
                            .abs()
                            .max(reconstructed.abs())
                            .max(1.0e-30)
                })
                .fold(0.0_f64, f64::max);
            let physical_enthalpy_error = normalization
                .denormalize_extensive_error(normalized.enthalpy_error())
                .expect("P,H enthalpy error must denormalize");
            println!(
                "  probe {temperature_k:.1} K: base_T={:.6} normalized_T={:.6} base_graphite={} normalized_graphite={} base_total={:.6e} normalized_total={:.6e} scale={normalization_factor:.6e} max_n={max_mole_error:.3e} H_internal={:.3e} H_physical={physical_enthalpy_error:.3e}",
                base.temperature(),
                normalized.temperature(),
                fixture.graphite_is_active(base.equilibrium()),
                fixture.graphite_is_active(normalized.equilibrium()),
                base.equilibrium().component_moles().iter().sum::<f64>(),
                normalized
                    .equilibrium()
                    .component_moles()
                    .iter()
                    .sum::<f64>(),
                normalized.enthalpy_error(),
            );
            assert!(
                max_mole_error <= PH_EXTENSIVITY_MOLE_RELATIVE_TOLERANCE
                    && normalized.enthalpy_error().abs()
                        <= normalized.enthalpy_error_limit_joules()
                    && physical_enthalpy_error.abs()
                        <= normalization
                            .denormalize_extensive_error(normalized.enthalpy_error_limit_joules(),)
                            .expect("P,H enthalpy limit must denormalize"),
                "P,H normalized solution must reconstruct the physical extensive state at {temperature_k} K: moles={max_mole_error:e}, internal_H={:e}, physical_H={physical_enthalpy_error:e}",
                normalized.enthalpy_error(),
            );
            println!(
                "{temperature_k:10.3}  {factor:15.1e}  {internal_total:14.6e}  {temperature:12.6}  {delta_temperature:9.3e}  {graphite:>8}  {max_mole_error:14.3e}  {internal_enthalpy:12.3e}  {physical_enthalpy_error:12.3e}  OK",
                temperature = normalized.temperature(),
                delta_temperature = normalized.temperature() - base.temperature(),
                graphite = fixture.graphite_is_active(normalized.equilibrium()),
                internal_enthalpy = normalized.enthalpy_error(),
            );
        }
    }

    #[test]
    #[ignore = "release story for transactional P,H range normalization recovery"]
    fn i5_tp1906_extensive_ph_target_range_recovery_is_transactional() {
        let fixture = fixture();
        let cases = fixture.cases().unwrap();
        let factor = 1.0e4_f64;
        let selected = [700.0, 720.0]
            .into_iter()
            .map(|temperature_k| {
                cases
                    .iter()
                    .find(|case| {
                        (case.enthalpy.temperature_k - temperature_k).abs() <= f64::EPSILON
                    })
                    .expect("range story must retain both frozen source rows")
            })
            .collect::<Vec<_>>();
        let physical = scaled_initial_composition(&fixture, factor);
        let targets = PhEnthalpyGrid::new(
            selected
                .iter()
                .map(|case| factor * case.target_enthalpy_j)
                .collect(),
        )
        .expect("scaled frozen targets must be strictly monotone");
        let range = PhRangeRequest::from_resolved_thermochemistry(
            fixture.tp1907().resolved(),
            physical,
            101_325.0,
            101_325.0,
            targets,
            fixture.tp1907().thermochemistry().temperature_bounds(),
            900.0,
            fixture.tp1907().thermochemistry().clone(),
        )
        .expect("large-inventory P,H range request must validate")
        .with_solve_options(extensivity_solve_options())
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .solve()
        .expect("normalized recovery must publish the complete P,H range");

        assert_eq!(range.points().len(), 2);
        assert_eq!(range.report().point_count(), 2);
        assert_eq!(range.report().continuation_points(), 1);
        assert_eq!(
            range.points()[0].report().preparation(),
            PhRangePointPreparation::Initial
        );
        assert_eq!(
            range.points()[1].report().preparation(),
            PhRangePointPreparation::Continued
        );
        assert_eq!(
            range.report().extensive_normalization_recoveries(),
            range
                .points()
                .iter()
                .filter(|point| point.report().used_extensive_normalization_recovery())
                .count()
        );
        for (point, case) in range.points().iter().zip(selected) {
            assert!(point.solution().temperature().is_finite());
            assert!(
                point.solution().enthalpy_error().abs()
                    <= point.solution().enthalpy_error_limit_joules()
            );
            assert_eq!(
                fixture.graphite_is_active(point.solution().equilibrium()),
                case.expected_graphite_active()
            );
        }
        println!("NASA TP-1906/1907 transactional P,H range normalization recovery");
        println!(
            "factor={factor:.1e} points={} recoveries={}",
            range.points().len(),
            range.report().extensive_normalization_recoveries()
        );
        for point in range.points() {
            println!(
                "  {:13.6e} {:>10?} T={:10.6} recovery={} transitions={}",
                point.report().target_enthalpy_joules(),
                point.report().preparation(),
                point.solution().temperature(),
                point.report().used_extensive_normalization_recovery(),
                point.report().phase_control_transitions(),
            );
        }
    }

    #[test]
    #[ignore = "release I5 TP-1906/1907 P,H preflight and lifecycle characterization"]
    fn i5_nasa_tp1906_tp1907_chon_graphite_ph_production_lifecycle_characterization() {
        let before_library = local_library_snapshot();
        let before_frozen = frozen_snapshot();
        let fixture = fixture();
        let cases = fixture.cases().unwrap();

        println!("NASA TP-1906/1907 CHON+graphite P,H characterization");
        println!(
            "P=101325 Pa inventory_mass_g={:.9}",
            fixture.inventory_mass_g()
        );
        println!("P,T enthalpy-reference preflight:");
        println!(
            "  T source K    h source J/g    h local J/g    delta J/g   graphite expected/local"
        );
        for case in &cases {
            let witness = fixture
                .solve_pt_witness(case.enthalpy.temperature_k)
                .unwrap();
            let preflight = fixture.enthalpy_preflight(case, &witness).unwrap();
            println!(
                "  {:10.3} {:15.6} {:14.6} {:12.6} {:>8}/{:<8}",
                preflight.temperature_k,
                preflight.source_specific_enthalpy_j_g,
                preflight.local_specific_enthalpy_j_g,
                preflight.delta_specific_enthalpy_j_g,
                case.expected_graphite_active(),
                fixture.graphite_is_active(&witness),
            );
        }

        println!("isolated P,H cases (common numerical seed=900 K):");
        println!(
            "  T source  H target J     T KiThe  delta T  graphite expected/local  H scaled  route/trials"
        );
        for case in &cases {
            let solution = solve_resolved_ph(
                fixture
                    .ph_request(case, 900.0)
                    .unwrap()
                    .with_ph_solve_mode(PhSolveMode::NestedTemperature),
            )
            .unwrap();
            println!(
                "  {:8.3} {:13.6e} {:10.3} {:8.3} {:>8}/{:<8} {:9.3e} {:>18}",
                case.enthalpy.temperature_k,
                case.target_enthalpy_j,
                solution.temperature(),
                solution.temperature() - case.enthalpy.temperature_k,
                case.expected_graphite_active(),
                fixture.graphite_is_active(solution.equilibrium()),
                solution.scaled_enthalpy_error(),
                format!(
                    "{:?}/{}",
                    solution.report().solve_path(),
                    solution.report().trials().len()
                ),
            );
            assert!(solution.scaled_enthalpy_error().is_finite());
            assert!(
                solution.enthalpy_error().abs() <= solution.enthalpy_error_limit_joules(),
                "accepted P,H result must satisfy its own energy contract"
            );
            assert_eq!(
                fixture.graphite_is_active(solution.equilibrium()),
                case.expected_graphite_active(),
                "external topology characterization must retain its full diagnostics before any policy change"
            );
            assert_ph_matches_pt_witness(&fixture, &solution);
            let comparison = fixture
                .tp1907()
                .compare_system_composition(&case.composition, solution.equilibrium())
                .unwrap();
            assert_eq!(
                comparison.len(),
                case.composition.gas_species.len() + case.composition.condensed_species.len(),
                "P,H result must use the established TP-1907 identity comparator"
            );
        }

        let forward_targets =
            PhEnthalpyGrid::new(cases.iter().map(|case| case.target_enthalpy_j).collect()).unwrap();
        let forward = PhRangeRequest::from_resolved_thermochemistry(
            fixture.tp1907().resolved(),
            fixture.tp1907().initial_composition().unwrap(),
            101_325.0,
            101_325.0,
            forward_targets,
            fixture.tp1907().thermochemistry().temperature_bounds(),
            900.0,
            fixture.tp1907().thermochemistry().clone(),
        )
        .unwrap()
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .solve()
        .unwrap();
        assert_eq!(forward.points().len(), cases.len());
        assert_eq!(
            forward.points()[0].report().preparation(),
            PhRangePointPreparation::Initial
        );
        assert!(
            forward.points()[1..]
                .iter()
                .all(|point| point.report().preparation() == PhRangePointPreparation::Continued)
        );
        println!("forward accepted P,H continuation:");
        println!(
            "  H target J       seed/continued  T KiThe  graphite  accepted transitions  trial events"
        );
        for point in forward.points() {
            println!(
                "  {:13.6e} {:>14?} {:8.3} {:>8} {:>10} {:>13}",
                point.report().target_enthalpy_joules(),
                point.report().preparation(),
                point.solution().temperature(),
                fixture.graphite_is_active(point.solution().equilibrium()),
                point.report().phase_control_transitions(),
                point.solution().report().phase_control_transitions(),
            );
        }

        let reverse_targets = PhEnthalpyGrid::new(
            cases
                .iter()
                .rev()
                .map(|case| case.target_enthalpy_j)
                .collect(),
        )
        .unwrap();
        let reverse = PhRangeRequest::from_resolved_thermochemistry(
            fixture.tp1907().resolved(),
            fixture.tp1907().initial_composition().unwrap(),
            101_325.0,
            101_325.0,
            reverse_targets,
            fixture.tp1907().thermochemistry().temperature_bounds(),
            900.0,
            fixture.tp1907().thermochemistry().clone(),
        )
        .unwrap()
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .solve()
        .unwrap();
        assert_eq!(
            reverse.points()[0].report().preparation(),
            PhRangePointPreparation::Initial
        );
        assert!(
            reverse.points()[1..]
                .iter()
                .all(|point| point.report().preparation() == PhRangePointPreparation::Continued)
        );
        let expected_forward = cases
            .iter()
            .map(|case| case.expected_graphite_active())
            .collect::<Vec<_>>();
        let expected_reverse = expected_forward.iter().rev().copied().collect::<Vec<_>>();
        let forward_history = graphite_history(&fixture, &forward);
        let reverse_history = graphite_history(&fixture, &reverse);
        assert_eq!(forward_history, expected_forward);
        assert_eq!(reverse_history, expected_reverse);
        assert_eq!(
            forward_history
                .windows(2)
                .filter(|pair| pair[0] && !pair[1])
                .count(),
            1,
            "forward accepted history must contain one graphite disappearance"
        );
        assert_eq!(
            reverse_history
                .windows(2)
                .filter(|pair| !pair[0] && pair[1])
                .count(),
            1,
            "reverse accepted history must contain one graphite appearance"
        );
        assert_eq!(
            forward_history
                .windows(2)
                .filter(|pair| pair[0] != pair[1])
                .count(),
            1,
            "forward accepted topology must not chatter"
        );
        assert_eq!(
            reverse_history
                .windows(2)
                .filter(|pair| pair[0] != pair[1])
                .count(),
            1,
            "reverse accepted topology must not chatter"
        );
        println!("reverse accepted P,H continuation:");
        println!(
            "  H target J       seed/continued  T KiThe  graphite  accepted transitions  trial events"
        );
        for point in reverse.points() {
            println!(
                "  {:13.6e} {:>14?} {:8.3} {:>8} {:>10} {:>13}",
                point.report().target_enthalpy_joules(),
                point.report().preparation(),
                point.solution().temperature(),
                fixture.graphite_is_active(point.solution().equilibrium()),
                point.report().phase_control_transitions(),
                point.solution().report().phase_control_transitions(),
            );
        }
        assert_forward_reverse_agree(&forward, &reverse);

        assert_eq!(before_library, local_library_snapshot());
        assert_eq!(before_frozen, frozen_snapshot());
    }

    #[test]
    #[ignore = "diagnostic route matrix for the TP-1906/1907 graphite P,H witness contract"]
    fn i5_nasa_tp1906_tp1907_chon_graphite_ph_route_witness_diagnostic() {
        let fixture = fixture();
        let cases = fixture.cases().unwrap();
        println!("NASA TP-1906/1907 P,H route vs canonical P,T witness");
        println!(
            "  T source  profile/route          accepted route              max relative mole delta"
        );
        let mut required_route_successes = 0;
        for case in &cases[..2] {
            for (profile, mode, options) in [
                (
                    "default",
                    PhSolveMode::NestedTemperature,
                    EquilibriumSolveOptions::default(),
                ),
                (
                    "default",
                    PhSolveMode::Monolithic,
                    EquilibriumSolveOptions::default(),
                ),
                (
                    "default",
                    PhSolveMode::Auto,
                    EquilibriumSolveOptions::default(),
                ),
                (
                    "tight-1e-8",
                    PhSolveMode::Auto,
                    EquilibriumSolveOptions::default()
                        .with_tolerance(1.0e-8)
                        .unwrap(),
                ),
            ] {
                let required_route =
                    matches!(mode, PhSolveMode::NestedTemperature | PhSolveMode::Auto);
                let result = solve_resolved_ph(
                    fixture
                        .ph_request(case, 900.0)
                        .unwrap()
                        .with_ph_solve_mode(mode)
                        .with_solve_options(options),
                );
                match result {
                    Ok(solution) => {
                        assert_ph_matches_pt_witness(&fixture, &solution);
                        if required_route {
                            required_route_successes += 1;
                        }
                        println!(
                            "  {:8.3} {:>11}/{:<8?} {:>27?} {:24.6e}",
                            case.enthalpy.temperature_k,
                            profile,
                            mode,
                            solution.report().solve_path(),
                            ph_pt_witness_max_relative_difference(&fixture, &solution),
                        );
                        print_ph_pt_witness_delta(&fixture, &solution);
                    }
                    Err(error) => {
                        assert!(
                            !required_route,
                            "canonical {mode:?} route must solve the witness case at {} K: {error}",
                            case.enthalpy.temperature_k
                        );
                        println!(
                            "  {:8.3} {:>11}/{:<8?} {:>27} {}",
                            case.enthalpy.temperature_k, profile, mode, "FAILED", error
                        );
                    }
                }
            }
        }
        assert_eq!(
            required_route_successes, 6,
            "NestedTemperature and both Auto profiles must retain the two canonical P,H witnesses"
        );
    }
}
