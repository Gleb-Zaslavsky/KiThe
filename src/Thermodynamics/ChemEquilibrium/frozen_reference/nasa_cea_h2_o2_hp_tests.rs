//! Loader and offline-inventory tests for the first NASA CEA full-equilibrium I5 case.

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
        EquilibriumConstraint, TemperatureBounds,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
        PhMonolithicSeedPolicy, PhSolveMode, PhTemperatureSolveOptions,
        ResolvedPhaseEnthalpyRequest, solve_resolved_ph,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::TraceSpeciesSeedPolicy;
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_h2_o2_hp::{
        CEA_COMPONENT_MAPPINGS, CEA_DECLARED_SPECIES, NasaCeaEquilibriumComparisonReport,
        NasaCeaPhRouteDiagnostic, NasaCeaSpeciesMagnitudeClass, ResolvedNasaCeaH2O2HpFixture,
        ResolvedNasaCeaH2O2HpGasFixture, load_nasa_cea_h2_o2_hp_dataset,
        preflight_nasa_cea_h2_o2_hp, reconstruct_nasa_cea_h2_o2_reactants,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::EquilibriumSolveOptions;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;

    // Regression guards selected after repeated debug/release characterization.
    // They are not CEA source uncertainty or physical acceptance tolerances.
    const MAX_TEMPERATURE_DELTA_K_GUARD: f64 = 10.0;
    const MAX_TOTAL_AMOUNT_RELATIVE_ERROR_GUARD: f64 = 0.02;
    const MAX_MAJOR_RELATIVE_ERROR_GUARD: f64 = 0.10;
    const MAX_MAJOR_RMS_RELATIVE_ERROR_GUARD: f64 = 0.075;
    const MAX_MINOR_OR_TRACE_LOG10_ERROR_GUARD: f64 = 0.10;

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("the bundled offline thermochemistry repository must be available")
    }

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
                (
                    path.clone(),
                    fs::read(&path).expect("local database must be readable"),
                )
            })
            .collect()
    }

    fn frozen_dataset_snapshot() -> Vec<(PathBuf, Vec<u8>)> {
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nasa_cea");
        [
            directory.join("h2_o2_hp_scitech_2025.metadata.json"),
            directory.join("h2_o2_hp_scitech_2025.rows.json"),
        ]
        .into_iter()
        .map(|path| {
            let bytes = fs::read(&path).expect("frozen CEA fixture must be readable");
            (path, bytes)
        })
        .collect()
    }

    fn render_comparison(report: &NasaCeaEquilibriumComparisonReport) {
        println!("NASA CEA H2/O2 HP gas-only external characterization");
        println!(
            "dataset={} | route={:?} | P={:.0} Pa (activity reference {:.0} Pa)",
            report.dataset_id,
            report.solve_path,
            report.physical_pressure_pa,
            report.activity_reference_pressure_pa,
        );
        if let Some(reason) = &report.fallback_reason {
            println!("Auto fallback: {reason}");
        }
        println!(
            "reactants T={:.3} K H_target={:.9e} J | T: CEA={:.3} K KiThe={:.3} K delta={:+.3} K",
            report.reactant_temperature_k,
            report.target_enthalpy_j,
            report.cea_temperature_k,
            report.kithe_temperature_k,
            report.temperature_delta_k,
        );
        println!(
            "total: CEA={:.9e} KiThe={:.9e} kgmol/kg",
            report.cea_total_kgmol_per_kg, report.kithe_total_kgmol_per_kg,
        );
        println!(
            "H_error={:.3e} J scaled_H={:.3e} residual={:.3e} balance={:.3e} final_transitions={} trial_events={}",
            report.enthalpy_error_j,
            report.scaled_enthalpy_error,
            report.residual_l2_norm,
            report.max_abs_element_balance_error,
            report.final_phase_transitions,
            report.trial_phase_events,
        );
        let summary = &report.external_summary;
        println!(
            "external deltas: |dT|={:.6} K relT={:+.3e} |dtotal|={:.7e} reltotal={:+.3e}",
            summary.absolute_temperature_delta_k,
            summary.relative_temperature_delta,
            summary.absolute_total_amount_delta_kgmol_per_kg,
            summary.relative_total_amount_delta,
        );
        if let Some(major) = summary.major_relative {
            println!(
                "major rel: n={} max={:.3e} rms={:.3e} mean={:+.3e}",
                major.sample_count,
                major.max_absolute_relative_error,
                major.rms_relative_error,
                major.mean_signed_relative_error,
            );
        }
        for (label, values) in [
            ("minor log10", summary.minor_log10),
            ("trace log10", summary.trace_log10),
        ] {
            if let Some(values) = values {
                println!(
                    "{label}: n={} max={:.3e} rms={:.3e} mean={:+.3e}",
                    values.sample_count,
                    values.max_absolute_delta_log10,
                    values.rms_delta_log10,
                    values.mean_signed_delta_log10,
                );
            }
        }
        println!(
            "{:<9} {:<19} {:>14} {:>14} {:>14} {:>12} {:>12}",
            "CEA", "class", "CEA kgmol/kg", "KiThe", "abs delta", "relative", "dlog10"
        );
        for row in &report.species {
            let local = row
                .kithe_amount_kgmol_per_kg
                .map(|value| format!("{value:.7e}"))
                .unwrap_or_else(|| "excluded".to_owned());
            let absolute = row
                .absolute_error
                .map(|value| format!("{value:+.7e}"))
                .unwrap_or_else(|| "-".to_owned());
            let relative = row
                .relative_error
                .map(|value| format!("{value:+.3e}"))
                .unwrap_or_else(|| "-".to_owned());
            let delta_log10 = row
                .delta_log10
                .map(|value| format!("{value:+.3}"))
                .unwrap_or_else(|| "-".to_owned());
            println!(
                "{:<9} {:<19} {:>14.7e} {:>14} {:>14} {:>12} {:>12}",
                row.cea_identity,
                format!("{:?}", row.magnitude_class),
                row.cea_amount_kgmol_per_kg,
                local,
                absolute,
                relative,
                delta_log10,
            );
            if let Some(reason) = &row.exclusion_reason {
                println!("          note: {reason}");
            }
        }
    }

    fn assert_cea_external_regression_envelope(report: &NasaCeaEquilibriumComparisonReport) {
        let summary = &report.external_summary;
        assert!(
            summary.absolute_temperature_delta_k <= MAX_TEMPERATURE_DELTA_K_GUARD,
            "NASA CEA H2/O2 external quality regression: |delta T|={:e} K exceeds reviewed guard={MAX_TEMPERATURE_DELTA_K_GUARD:e} K",
            summary.absolute_temperature_delta_k,
        );
        assert!(
            summary.relative_total_amount_delta.abs() <= MAX_TOTAL_AMOUNT_RELATIVE_ERROR_GUARD,
            "NASA CEA H2/O2 external quality regression: relative total-amount error={:e} exceeds reviewed guard={MAX_TOTAL_AMOUNT_RELATIVE_ERROR_GUARD:e}",
            summary.relative_total_amount_delta,
        );
        let major = summary
            .major_relative
            .expect("CEA fixture must retain positive major species");
        assert!(
            major.max_absolute_relative_error <= MAX_MAJOR_RELATIVE_ERROR_GUARD,
            "NASA CEA H2/O2 external quality regression: max major relative error={:e} exceeds reviewed guard={MAX_MAJOR_RELATIVE_ERROR_GUARD:e}",
            major.max_absolute_relative_error,
        );
        assert!(
            major.rms_relative_error <= MAX_MAJOR_RMS_RELATIVE_ERROR_GUARD,
            "NASA CEA H2/O2 external quality regression: RMS major relative error={:e} exceeds reviewed guard={MAX_MAJOR_RMS_RELATIVE_ERROR_GUARD:e}",
            major.rms_relative_error,
        );
        for (label, summary) in [
            ("minor", summary.minor_log10),
            ("trace", summary.trace_log10),
        ] {
            if let Some(summary) = summary {
                assert!(
                    summary.max_absolute_delta_log10 <= MAX_MINOR_OR_TRACE_LOG10_ERROR_GUARD,
                    "NASA CEA H2/O2 external quality regression: max {label} |delta log10|={:e} exceeds reviewed guard={MAX_MINOR_OR_TRACE_LOG10_ERROR_GUARD:e}",
                    summary.max_absolute_delta_log10,
                );
            }
        }
    }

    fn render_route_diagnostic(diagnostic: &NasaCeaPhRouteDiagnostic) {
        match &diagnostic.comparison {
            Some(report) => println!(
                "route={:?} seed={:.3} trace={:.0e} status=OK path={:?} T={:.3} scaled_H={:.3e} residual={:.3e} balance={:.3e} trials={} attempts={} iterations={} fallback={}",
                diagnostic.requested_mode,
                diagnostic.temperature_seed_k,
                diagnostic.trace_floor,
                diagnostic.solve_path,
                report.kithe_temperature_k,
                report.scaled_enthalpy_error,
                report.residual_l2_norm,
                report.max_abs_element_balance_error,
                diagnostic.temperature_trials.unwrap_or(0),
                diagnostic.inner_backend_attempts.unwrap_or(0),
                diagnostic.inner_nonlinear_iterations.unwrap_or(0),
                diagnostic.fallback_reason.as_deref().unwrap_or("-"),
            ),
            None => println!(
                "route={:?} seed={:.3} trace={:.0e} status=FAILED kind={:?} attempts={} message={}",
                diagnostic.requested_mode,
                diagnostic.temperature_seed_k,
                diagnostic.trace_floor,
                diagnostic.failure_kind,
                diagnostic.backend_attempts.len(),
                diagnostic.failure_message.as_deref().unwrap_or("-"),
            ),
        }
        for (index, attempt) in diagnostic.temperature_seed_attempts.iter().enumerate() {
            println!(
                "  seed_attempt={} T={:.3} status={} failure={:?} backend_attempts={} iterations={} error={}",
                index,
                attempt.temperature_seed_k(),
                if attempt.accepted() { "OK" } else { "FAILED" },
                attempt.failure_kind(),
                attempt.started_backend_attempts(),
                attempt.nonlinear_iterations(),
                attempt.error().unwrap_or("-"),
            );
        }
        for attempt in &diagnostic.backend_attempts {
            println!(
                "  backend={} outcome={} failure={:?} termination={:?} iter={:?} residual_evals={:?} jacobian_evals={:?} reason={}",
                attempt.backend,
                attempt.outcome,
                attempt.failure_kind,
                attempt.termination,
                attempt.iterations,
                attempt.residual_evaluations,
                attempt.jacobian_evaluations,
                attempt.reason.as_deref().unwrap_or("-"),
            );
        }
    }

    fn assert_internal_route_agreement(
        nested: &NasaCeaEquilibriumComparisonReport,
        monolithic: &NasaCeaEquilibriumComparisonReport,
    ) {
        let temperature_scale = nested.kithe_temperature_k.abs().max(1.0);
        assert!(
            (nested.kithe_temperature_k - monolithic.kithe_temperature_k).abs() / temperature_scale
                <= 1.0e-6,
            "nested and monolithic temperatures must agree for the same local P,H problem"
        );
        let total_scale = nested.kithe_total_kgmol_per_kg.abs().max(1.0e-12);
        assert!(
            (nested.kithe_total_kgmol_per_kg - monolithic.kithe_total_kgmol_per_kg).abs()
                / total_scale
                <= 1.0e-6,
            "nested and monolithic total amounts must agree"
        );
        for (nested_row, monolithic_row) in nested.species.iter().zip(&monolithic.species) {
            assert_eq!(nested_row.cea_identity, monolithic_row.cea_identity);
            match (
                nested_row.kithe_amount_kgmol_per_kg,
                monolithic_row.kithe_amount_kgmol_per_kg,
            ) {
                (Some(left), Some(right)) => {
                    assert!(
                        (left - right).abs() / left.abs().max(1.0e-12) <= 1.0e-6,
                        "nested and monolithic amounts differ for {}",
                        nested_row.cea_identity
                    );
                }
                (None, None) => {}
                _ => panic!(
                    "route comparison changed local participation for {}",
                    nested_row.cea_identity
                ),
            }
        }
        for (name, nested_value, monolithic_value) in [
            (
                "scaled enthalpy error",
                nested.scaled_enthalpy_error,
                monolithic.scaled_enthalpy_error,
            ),
            (
                "residual norm",
                nested.residual_l2_norm,
                monolithic.residual_l2_norm,
            ),
            (
                "element balance",
                nested.max_abs_element_balance_error,
                monolithic.max_abs_element_balance_error,
            ),
        ] {
            assert!(
                nested_value.is_finite() && monolithic_value.is_finite(),
                "{name} must remain finite for both accepted P,H routes"
            );
            assert!(
                nested_value <= 1.0e-6 && monolithic_value <= 1.0e-6,
                "{name} must satisfy the strict internal acceptance envelope: nested={nested_value:.3e}, monolithic={monolithic_value:.3e}"
            );
        }
        assert_eq!(
            nested.final_phase_transitions, monolithic.final_phase_transitions,
            "equivalent fixed-gas P,H routes must publish the same phase lifecycle"
        );
    }

    #[test]
    fn i5_nasa_cea_hp_frozen_case_has_exact_semantic_universe() {
        let before = frozen_dataset_snapshot();
        let dataset =
            load_nasa_cea_h2_o2_hp_dataset().expect("reviewed NASA CEA dataset must load");
        assert!(dataset.is_external_evidence());
        assert_eq!(
            dataset.metadata().dataset_id,
            "nasa_cea.h2_o2.hp.scitech_2025.v1"
        );
        assert_eq!(
            dataset.metadata().source.stable_identifier.as_deref(),
            Some("NASA-20240016039")
        );
        let row = dataset.rows().first().expect("one CEA case must exist");
        assert_eq!(row.problem_type, "HP");
        assert_eq!(row.pressure_pa, 101_325.0);
        assert_eq!(row.reactant_temperature_k, 2_000.0);
        assert_eq!(row.species_amounts.len(), CEA_DECLARED_SPECIES.len());
        assert_eq!(
            row.species_amounts
                .iter()
                .map(|amount| amount.species.as_str())
                .collect::<Vec<_>>(),
            CEA_DECLARED_SPECIES
        );
        for condensed in ["H2O(L)", "H2O(cr)"] {
            assert_eq!(
                row.species_amounts
                    .iter()
                    .find(|amount| amount.species == condensed)
                    .expect("the frozen CEA universe must retain condensed water")
                    .amount_kgmol_per_kg,
                0.0,
                "gas-only execution is valid only while CEA reports {condensed} as absent"
            );
        }
        assert_eq!(CEA_COMPONENT_MAPPINGS.len(), CEA_DECLARED_SPECIES.len());
        assert_eq!(before, frozen_dataset_snapshot());
    }

    #[test]
    fn i5_nasa_cea_hp_preflight_never_substitutes_an_out_of_domain_condensed_record() {
        let before = local_library_snapshot();
        let dataset = load_nasa_cea_h2_o2_hp_dataset().unwrap();
        let reference = dataset.rows().first().unwrap();
        let reactants = reconstruct_nasa_cea_h2_o2_reactants(local_repository(), reference)
            .expect("exact local H2/O2 input preparation must work at 2000 K");
        assert!((reactants.total_mass_kg - 1.0).abs() < f64::EPSILON);
        assert!((reactants.hydrogen_mass_kg + reactants.oxygen_mass_kg - 1.0).abs() < f64::EPSILON);
        assert!(reactants.hydrogen_molar_mass_kg_mol > 0.002);
        assert!(reactants.oxygen_molar_mass_kg_mol > 0.031);
        assert!(reactants.hydrogen_moles > 0.0 && reactants.oxygen_moles > 0.0);
        assert!(reactants.target_enthalpy_j.is_finite());
        let preflight = preflight_nasa_cea_h2_o2_hp(local_repository(), reference).unwrap();
        assert_eq!(preflight.rows().len(), CEA_DECLARED_SPECIES.len());
        let (common_lower, common_upper) = preflight.common_temperature_interval();
        assert_eq!(common_lower, 273.15);
        assert_eq!(common_upper, 273.15);
        assert!(
            preflight.rows().iter().any(|row| {
                row.cea_identity == "H2O(L)" && !row.supports_equilibrium_temperature
            })
        );
        assert!(
            preflight.rows().iter().any(|row| {
                row.cea_identity == "H2O(cr)" && !row.supports_equilibrium_temperature
            })
        );
        match ResolvedNasaCeaH2O2HpFixture::resolve_offline(local_repository(), reference) {
            Ok(fixture) => {
                let reconstruction = fixture.reconstruct_reactants(reference).unwrap();
                let initial = fixture.initial_composition(&reconstruction).unwrap();
                let structure = fixture.structure();
                assert_eq!(fixture.component_inventory().len(), 11);
                assert_eq!(fixture.resolved().layout().component_count(), 11);
                assert_eq!(
                    initial.moles().iter().filter(|&&moles| moles > 0.0).count(),
                    2
                );
                assert_eq!(structure.component_count, 11);
                assert_eq!(structure.element_count, 2);
                assert_eq!(structure.element_rank, 2);
                assert_eq!(structure.reaction_dimension, 9);
                assert!(
                    fixture
                        .thermochemistry()
                        .temperature_bounds()
                        .contains(reference.equilibrium_temperature_k)
                );
                assert!((reconstruction.total_mass_kg - 1.0).abs() < f64::EPSILON);
                assert!(reconstruction.hydrogen_moles > 0.0 && reconstruction.oxygen_moles > 0.0);
                assert!(reconstruction.target_enthalpy_j.is_finite());
                for (row, mapping) in fixture
                    .component_inventory()
                    .iter()
                    .zip(CEA_COMPONENT_MAPPINGS)
                {
                    assert_eq!(row.cea_identity, mapping.cea_identity);
                    assert_eq!(row.component.substance, mapping.local_species);
                    assert_eq!(row.library, mapping.library);
                    assert_eq!(row.physical_state, mapping.physical_state);
                    assert!(row.gibbs_j_mol.is_finite());
                    assert!(row.enthalpy_j_mol.is_finite());
                    assert!(row.heat_capacity_j_mol_k.is_finite());
                }
            }
            Err(ReactionExtentError::ValidationNotApplicable { path, message }) => {
                // Current local H2O(L)/H2O(s) records collapse the shared
                // interval to 273.15 K. This is a data/capability gap, not a
                // license to omit them from the declared CEA universe.
                assert_eq!(path, "nasa_cea_h2_o2_hp_temperature_domain");
                assert!(message.contains("3181.23") || message.contains("2000"));
                assert!(message.contains("H2O(L)"));
                assert!(message.contains("H2O(cr)"));
            }
            Err(error) => panic!("CEA preflight must be typed, got {error:?}"),
        }
        assert_eq!(before, local_library_snapshot());
    }

    #[test]
    fn i5_nasa_cea_hp_gas_subset_is_an_exact_offline_solver_layout() {
        let before = local_library_snapshot();
        let dataset = load_nasa_cea_h2_o2_hp_dataset().unwrap();
        let reference = dataset.rows().first().unwrap();
        let repository = local_repository();
        let reactants = reconstruct_nasa_cea_h2_o2_reactants(repository.clone(), reference)
            .expect("exact local H2/O2 input preparation must work");
        let fixture = ResolvedNasaCeaH2O2HpGasFixture::resolve_offline(repository, reference)
            .expect("the exact nine-component gas universe must resolve offline");
        let initial = fixture.initial_composition(&reactants).unwrap();
        let broad = fixture
            .broad_initial_composition(&reactants, 1.0e-2)
            .unwrap();
        let structure = fixture.structure();

        assert_eq!(fixture.resolved().layout().component_count(), 9);
        assert_eq!(fixture.component_inventory().len(), 9);
        assert_eq!(structure.component_count, 9);
        assert_eq!(structure.element_count, 2);
        assert_eq!(structure.element_rank, 2);
        assert_eq!(structure.reaction_dimension, 7);
        assert_eq!(
            initial.moles().iter().filter(|&&moles| moles > 0.0).count(),
            2
        );
        assert_eq!(
            broad.moles().iter().filter(|&&moles| moles > 0.0).count(),
            9,
            "the diagnostic broad seed must make every exact gas component positive"
        );
        assert!(
            fixture
                .thermochemistry()
                .temperature_bounds()
                .contains(reference.equilibrium_temperature_k)
        );
        assert_eq!(before, local_library_snapshot());
    }

    #[test]
    fn i5_nasa_cea_hp_gas_only_rejects_nonzero_excluded_condensed_reference_rows() {
        let dataset = load_nasa_cea_h2_o2_hp_dataset().unwrap();
        let mut reference = dataset.rows().first().unwrap().clone();
        reference
            .species_amounts
            .iter_mut()
            .find(|row| row.species == "H2O(L)")
            .expect("the frozen declaration includes liquid water")
            .amount_kgmol_per_kg = 1.0e-12;

        let result =
            ResolvedNasaCeaH2O2HpGasFixture::resolve_offline(local_repository(), &reference);
        assert!(matches!(
            result,
            Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_gas_only_eligibility",
                ..
            })
        ));
    }

    #[test]
    #[ignore = "release I5 NASA CEA H2/O2 characterization with conservative software-regression envelope"]
    fn i5_nasa_cea_hp_gas_only_equilibrium_characterization() {
        let before_library = local_library_snapshot();
        let before_frozen = frozen_dataset_snapshot();
        let dataset = load_nasa_cea_h2_o2_hp_dataset().unwrap();
        let reference = dataset.rows().first().unwrap();
        let repository = local_repository();
        let reactants = reconstruct_nasa_cea_h2_o2_reactants(repository.clone(), reference)
            .expect("exact local H2/O2 input preparation must work");
        let fixture = ResolvedNasaCeaH2O2HpGasFixture::resolve_offline(repository, reference)
            .expect("the exact nine-component gas universe must resolve offline");
        let initial = fixture.initial_composition(&reactants).unwrap();
        let bounds = fixture.thermochemistry().temperature_bounds();
        let request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            fixture.resolved(),
            initial,
            EquilibriumConstraint::ph(
                reference.pressure_pa,
                reference.pressure_pa,
                reactants.target_enthalpy_j,
                reference.reactant_temperature_k,
            )
            .unwrap(),
            TemperatureBounds::new(bounds.lower(), bounds.upper()).unwrap(),
            fixture.thermochemistry().clone(),
        )
        .unwrap()
        .with_solve_options(
            EquilibriumSolveOptions::default()
                .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1.0e-30 }),
        )
        .with_ph_solve_mode(PhSolveMode::Auto);

        let solution = solve_resolved_ph(request)
            .expect("production P,H route must solve the exact nine-component gas benchmark");
        let report = NasaCeaEquilibriumComparisonReport::from_solution(
            &dataset,
            &reactants,
            reference.pressure_pa,
            &solution,
        )
        .unwrap();
        render_comparison(&report);

        assert_eq!(report.species.len(), CEA_DECLARED_SPECIES.len());
        assert_eq!(report.excluded_condensed_count(), 2);
        assert_eq!(
            report
                .species
                .iter()
                .filter(|row| row.local_component.is_some())
                .count(),
            9
        );
        assert!(report.kithe_temperature_k.is_finite());
        assert!(report.kithe_total_kgmol_per_kg.is_finite());
        assert!(report.residual_l2_norm.is_finite());
        assert!(report.max_abs_element_balance_error.is_finite());
        assert_cea_external_regression_envelope(&report);
        if report.solve_path
            == crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolvePath::NestedTemperature
        {
            assert!(
                report.fallback_reason.is_some(),
                "an Auto recovery must preserve the rejected monolithic reason"
            );
        }
        assert!(
            report
                .species
                .iter()
                .filter(|row| {
                    row.magnitude_class == NasaCeaSpeciesMagnitudeClass::ExternallyAbsentExcluded
                })
                .all(|row| row.cea_amount_kgmol_per_kg == 0.0
                    && row.kithe_amount_kgmol_per_kg.is_none()
                    && row.exclusion_reason.is_some())
        );
        assert_eq!(before_library, local_library_snapshot());
        assert_eq!(before_frozen, frozen_dataset_snapshot());
    }

    #[test]
    #[ignore = "release diagnostic for direct P,H route and seed sensitivity on NASA CEA H2/O2"]
    fn i5_nasa_cea_hp_gas_only_ph_route_matrix_diagnostic() {
        let before_library = local_library_snapshot();
        let before_frozen = frozen_dataset_snapshot();
        let dataset = load_nasa_cea_h2_o2_hp_dataset().unwrap();
        let reference = dataset.rows().first().unwrap();
        let repository = local_repository();
        let reactants = reconstruct_nasa_cea_h2_o2_reactants(repository.clone(), reference)
            .expect("exact local H2/O2 input preparation must work");
        let fixture = ResolvedNasaCeaH2O2HpGasFixture::resolve_offline(repository, reference)
            .expect("the exact nine-component gas universe must resolve offline");
        let bounds = fixture.thermochemistry().temperature_bounds();

        let run = |mode: PhSolveMode,
                   temperature_seed_k: f64,
                   trace_floor: f64,
                   broad: bool,
                   seed_policy: PhMonolithicSeedPolicy| {
            let initial = if broad {
                fixture
                    .broad_initial_composition(&reactants, 1.0e-2)
                    .expect("broad local gas seed must be valid")
            } else {
                fixture
                    .initial_composition(&reactants)
                    .expect("physical H2/O2 seed must be valid")
            };
            let temperature_options =
                PhTemperatureSolveOptions::default().with_monolithic_seed_policy(seed_policy);
            let request =
                ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
                    fixture.resolved(),
                    initial,
                    EquilibriumConstraint::ph(
                        reference.pressure_pa,
                        reference.pressure_pa,
                        reactants.target_enthalpy_j,
                        temperature_seed_k,
                    )
                    .expect("P,H constraint must be valid"),
                    TemperatureBounds::new(bounds.lower(), bounds.upper())
                        .expect("gas thermochemistry bounds must be valid"),
                    fixture.thermochemistry().clone(),
                )
                .expect("resolved gas P,H request must be valid")
                .with_solve_options(EquilibriumSolveOptions::default().with_trace_seed_policy(
                    TraceSpeciesSeedPolicy::Absolute { floor: trace_floor },
                ))
                .with_ph_solve_mode(mode)
                .with_temperature_options(temperature_options)
                .expect("temperature-seed policy must be valid");
            NasaCeaPhRouteDiagnostic::from_result(
                &dataset,
                &reactants,
                reference.pressure_pa,
                mode,
                temperature_seed_k,
                trace_floor,
                solve_resolved_ph(request),
            )
            .expect("route diagnostics must be constructible")
        };

        println!("NASA CEA H2/O2 HP gas-only P,H route matrix");
        let nested = run(
            PhSolveMode::NestedTemperature,
            reference.reactant_temperature_k,
            1.0e-30,
            false,
            PhMonolithicSeedPolicy::DeterministicBoundedRecovery,
        );
        let monolithic = run(
            PhSolveMode::Monolithic,
            reference.reactant_temperature_k,
            1.0e-30,
            false,
            PhMonolithicSeedPolicy::DeterministicBoundedRecovery,
        );
        let auto = run(
            PhSolveMode::Auto,
            reference.reactant_temperature_k,
            1.0e-30,
            false,
            PhMonolithicSeedPolicy::DeterministicBoundedRecovery,
        );
        for diagnostic in [&nested, &monolithic, &auto] {
            render_route_diagnostic(diagnostic);
        }
        assert!(
            nested.comparison.is_some(),
            "nested route is the reference solver"
        );
        assert!(
            auto.comparison.is_some(),
            "Auto must preserve numerical recovery"
        );
        assert!(
            monolithic.comparison.is_some(),
            "bounded temperature recovery must keep this case on the monolithic route"
        );
        assert_eq!(
            monolithic.solve_path,
            Some(crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolvePath::MonolithicFixedActiveSet)
        );
        assert!(
            matches!(
                monolithic.temperature_seed_attempts.as_slice(),
                [attempt] if attempt.accepted()
                    && (attempt.temperature_seed_k() - reference.reactant_temperature_k).abs()
                        <= f64::EPSILON
            ),
            "the dimensionally correct affinity scaling must let the physical 2000 K seed converge without bounded recovery: {:?}",
            monolithic.temperature_seed_attempts,
        );
        assert_eq!(auto.solve_path, monolithic.solve_path);
        assert!(auto.fallback_reason.is_none());

        println!("monolithic temperature-seed sensitivity");
        let warm_monolithic = run(
            PhSolveMode::Monolithic,
            3_000.0,
            1.0e-30,
            false,
            PhMonolithicSeedPolicy::InitialOnly,
        );
        render_route_diagnostic(&warm_monolithic);
        for temperature_seed_k in [
            reference.reactant_temperature_k,
            3_500.0,
            reference.equilibrium_temperature_k,
        ] {
            render_route_diagnostic(&run(
                PhSolveMode::Monolithic,
                temperature_seed_k,
                1.0e-30,
                false,
                PhMonolithicSeedPolicy::InitialOnly,
            ));
        }
        if let (Some(nested), Some(monolithic)) = (
            nested.comparison.as_ref(),
            warm_monolithic.comparison.as_ref(),
        ) {
            assert_internal_route_agreement(nested, monolithic);
            println!("nested/monolithic warm-seed internal state agreement: OK");
        }
        println!("monolithic trace-floor sensitivity");
        for trace_floor in [1.0e-20, 1.0e-25, 1.0e-30] {
            render_route_diagnostic(&run(
                PhSolveMode::Monolithic,
                reference.reactant_temperature_k,
                trace_floor,
                false,
                PhMonolithicSeedPolicy::InitialOnly,
            ));
        }
        println!("monolithic broad local composition probe");
        render_route_diagnostic(&run(
            PhSolveMode::Monolithic,
            reference.reactant_temperature_k,
            1.0e-30,
            true,
            PhMonolithicSeedPolicy::InitialOnly,
        ));

        assert_eq!(before_library, local_library_snapshot());
        assert_eq!(before_frozen, frozen_dataset_snapshot());
    }

    #[test]
    #[ignore = "diagnostic evidence for the currently blocked full NASA CEA I5 solve"]
    fn i5_nasa_cea_hp_exact_inventory_preflight_diagnostic() {
        let dataset = load_nasa_cea_h2_o2_hp_dataset().unwrap();
        let reference = dataset.rows().first().unwrap();
        let reactants = reconstruct_nasa_cea_h2_o2_reactants(local_repository(), reference)
            .expect("exact local H2/O2 input preparation must work at 2000 K");
        let preflight = preflight_nasa_cea_h2_o2_hp(local_repository(), reference).unwrap();
        assert_eq!(preflight.rows().len(), CEA_DECLARED_SPECIES.len());
        assert!(
            !preflight.supports_required_temperatures(),
            "the full eleven-component route must remain honestly blocked until both condensed-water records cover the HP temperature"
        );
        let (common_lower, common_upper) = preflight.common_temperature_interval();
        assert!(
            !(common_lower.is_finite() && common_upper.is_finite() && common_lower < common_upper),
            "the full exact universe is expected to have no common local temperature interval"
        );
        assert!(
            preflight
                .rows()
                .iter()
                .any(|row| !row.supports_equilibrium_temperature),
            "the blocked status must be traceable to an explicit record-domain gap"
        );
        println!("NASA CEA H2/O2 HP preflight");
        println!(
            "P={:.0} Pa | reactants T={:.2} K | O/F mass={:.5}",
            reference.pressure_pa,
            reference.reactant_temperature_k,
            reference.oxidizer_fuel_mass_ratio
        );
        println!(
            "H2: mass={:.9e} kg mol={:.9e} | O2: mass={:.9e} kg mol={:.9e} | H_target={:.9e} J",
            reactants.hydrogen_mass_kg,
            reactants.hydrogen_moles,
            reactants.oxygen_mass_kg,
            reactants.oxygen_moles,
            reactants.target_enthalpy_j,
        );
        println!(
            "{:<10} {:<18} {:<10} {:<14} {:>10} {:>10} {:>7} {:>7}",
            "CEA", "KiThe component", "library", "record", "Tmin", "Tmax", "2000", "3181"
        );
        for row in preflight.rows() {
            println!(
                "{:<10} {:<18} {:<10} {:<14} {:>10.2} {:>10.2} {:>7} {:>7}",
                row.cea_identity,
                row.component.label(),
                row.library,
                row.record_key,
                row.temperature_lower_k,
                row.temperature_upper_k,
                if row.supports_reactant_temperature {
                    "yes"
                } else {
                    "no"
                },
                if row.supports_equilibrium_temperature {
                    "yes"
                } else {
                    "no"
                },
            );
        }
        println!(
            "common interval=[{common_lower:.2}, {common_upper:.2}] K | full 11-component solve={}",
            if preflight.supports_required_temperatures() {
                "available"
            } else {
                "blocked"
            },
        );
    }
}
