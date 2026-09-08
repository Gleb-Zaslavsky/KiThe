//! Offline I5 tests for the Argonne/STANJAN general CHON fixed-`P,T` case.

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        ArgonneStanjanChonPtReference, FrozenReferenceDataset, FrozenReferenceError,
        FrozenReferenceEvidenceKind,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon::{
        ResolvedArgonneStanjanChonFixture, STANJAN_DECLARED_SPECIES, STANJAN_LOCAL_GAS_SPECIES,
        StanjanChonComparisonReport, StanjanChonMagnitudeClass, load_argonne_stanjan_chon_dataset,
        stanjan_conditions, verify_element_equivalent_feeds,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
    };
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;

    // Regression guards selected after repeated debug/release characterization.
    // They are not STANJAN source uncertainty or physical acceptance tolerances.
    const MAX_MAJOR_RELATIVE_ERROR_GUARD: f64 = 0.075;
    const MAX_MAJOR_RMS_RELATIVE_ERROR_GUARD: f64 = 0.05;
    const MAX_MINOR_LOG10_ERROR_GUARD: f64 = 0.10;
    const MAX_TRACE_LOG10_ERROR_GUARD: f64 = 0.15;

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
                let bytes = fs::read(&path).expect("local database must be readable");
                (path, bytes)
            })
            .collect()
    }

    fn frozen_dataset_snapshot() -> Vec<(PathBuf, Vec<u8>)> {
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/argonne_stanjan");
        [
            directory.join("pentane_methane_air_tp_2500k_35atm.metadata.json"),
            directory.join("pentane_methane_air_tp_2500k_35atm.rows.json"),
        ]
        .into_iter()
        .map(|path| {
            let bytes = fs::read(&path).expect("frozen STANJAN data must be readable");
            (path, bytes)
        })
        .collect()
    }

    fn dataset() -> FrozenReferenceDataset<ArgonneStanjanChonPtReference> {
        load_argonne_stanjan_chon_dataset().expect("reviewed STANJAN dataset must load")
    }

    fn render_comparison(report: &StanjanChonComparisonReport) {
        println!("Argonne/STANJAN CHON fixed-P,T external characterization");
        println!(
            "dataset={} | T={} K | P={} Pa (activity reference {} Pa)",
            report.dataset_id,
            report.temperature_k,
            report.physical_pressure_pa,
            report.activity_reference_pressure_pa,
        );
        println!(
            "structure: reactions={} | accepted_backend={} | transitions={}",
            report.reaction_dimension, report.accepted_backend, report.phase_control_transitions,
        );
        println!(
            "sum: STANJAN={:.9} KiThe={:.9} | residual={:.3e} balance={:.3e}",
            report.external_sum,
            report.kithe_sum,
            report.residual_l2_norm,
            report.max_abs_element_balance_error,
        );
        println!(
            "species   class                           STANJAN         KiThe        abs delta     relative       dlog10"
        );
        for row in &report.species {
            match row.kithe_mole_fraction {
                Some(local) => println!(
                    "{:<9} {:?} {:>15.7e} {:>15.7e} {:>+13.6e} {:>+12.3e} {:>+10.3}",
                    row.external_identity,
                    row.magnitude_class,
                    row.stanjan_mole_fraction,
                    local,
                    row.absolute_error.unwrap_or_default(),
                    row.relative_error.unwrap_or_default(),
                    row.delta_log10.unwrap_or_default(),
                ),
                None => println!(
                    "{:<9} {:?} {:>15.7e} {:>15}  note: {}",
                    row.external_identity,
                    row.magnitude_class,
                    row.stanjan_mole_fraction,
                    "excluded",
                    row.exclusion_reason.as_deref().unwrap_or("-")
                ),
            }
        }
        for (name, summary) in [
            ("major relative", report.major_relative),
            ("minor dlog10", report.minor_log10),
            ("trace dlog10", report.trace_log10),
        ] {
            if let Some(summary) = summary {
                println!(
                    "{name}: n={} max={:.3e} rms={:.3e} mean={:+.3e}",
                    summary.sample_count,
                    summary.maximum_absolute,
                    summary.rms,
                    summary.mean_signed,
                );
            }
        }
    }

    fn assert_stanjan_external_regression_envelope(report: &StanjanChonComparisonReport) {
        let major = report
            .major_relative
            .expect("STANJAN fixture must retain major species");
        assert!(
            major.maximum_absolute <= MAX_MAJOR_RELATIVE_ERROR_GUARD,
            "Argonne/STANJAN external quality regression: max major relative error={:e} exceeds reviewed guard={MAX_MAJOR_RELATIVE_ERROR_GUARD:e}",
            major.maximum_absolute,
        );
        assert!(
            major.rms <= MAX_MAJOR_RMS_RELATIVE_ERROR_GUARD,
            "Argonne/STANJAN external quality regression: RMS major relative error={:e} exceeds reviewed guard={MAX_MAJOR_RMS_RELATIVE_ERROR_GUARD:e}",
            major.rms,
        );
        for (label, summary, guard) in [
            ("minor", report.minor_log10, MAX_MINOR_LOG10_ERROR_GUARD),
            ("trace", report.trace_log10, MAX_TRACE_LOG10_ERROR_GUARD),
        ] {
            if let Some(summary) = summary {
                assert!(
                    summary.maximum_absolute <= guard,
                    "Argonne/STANJAN external quality regression: max {label} |delta log10|={:e} exceeds reviewed guard={guard:e}",
                    summary.maximum_absolute,
                );
            }
        }
    }

    #[test]
    fn i5_argonne_stanjan_chon_dataset_preserves_full_source_evidence() {
        let before = frozen_dataset_snapshot();
        let dataset = dataset();
        let reference = &dataset.rows()[0];

        assert_eq!(
            dataset.metadata().evidence_kind,
            FrozenReferenceEvidenceKind::FrozenExternal
        );
        assert_eq!(reference.problem_type, "TP");
        assert_eq!(reference.temperature_k, 2_500.0);
        assert_eq!(reference.pressure_pa, 35.0 * 101_325.0);
        assert_eq!(
            reference
                .species_mole_fractions
                .iter()
                .map(|row| row.species.as_str())
                .collect::<Vec<_>>(),
            STANJAN_DECLARED_SPECIES
        );
        assert_eq!(
            reference
                .species_mole_fractions
                .iter()
                .find(|row| row.species == "C5H12")
                .expect("source table must retain pentane")
                .mole_fraction,
            0.0,
        );
        assert_eq!(before, frozen_dataset_snapshot());
    }

    #[test]
    fn i5_argonne_stanjan_chon_element_equivalent_feed_is_exact() {
        let dataset = dataset();
        let reference = &dataset.rows()[0];
        let totals = verify_element_equivalent_feeds(reference)
            .expect("published and executable feed must have identical C/H/O/N totals");
        assert_eq!(totals.get("C"), Some(&6.0));
        assert_eq!(totals.get("H"), Some(&16.0));
        assert_eq!(totals.get("O"), Some(&20.0));
        assert_eq!(totals.get("N"), Some(&75.2));
    }

    #[test]
    fn i5_argonne_stanjan_chon_rejects_positive_pentane_in_reduced_universe_source() {
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/argonne_stanjan");
        let metadata =
            fs::read_to_string(directory.join("pentane_methane_air_tp_2500k_35atm.metadata.json"))
                .expect("frozen metadata must be readable");
        let mut rows: serde_json::Value = serde_json::from_slice(
            &fs::read(directory.join("pentane_methane_air_tp_2500k_35atm.rows.json"))
                .expect("frozen rows must be readable"),
        )
        .expect("frozen rows must be valid JSON before mutation");
        rows["rows"][0]["species_mole_fractions"][0]["mole_fraction"] = serde_json::json!(1.0e-8);

        let error = FrozenReferenceDataset::<ArgonneStanjanChonPtReference>::from_json_strs(
            &metadata,
            &serde_json::to_string(&rows).expect("mutated rows must serialize"),
            "pentane_methane_air_tp_2500k_35atm.rows.json",
        )
        .expect_err("a positive external pentane row invalidates the reduced local universe");
        assert!(matches!(
            error,
            FrozenReferenceError::InvalidRows { field, .. }
                if field == "species_mole_fractions.C5H12"
        ));
    }

    #[test]
    fn i5_argonne_stanjan_chon_exact_universe_resolves_offline_with_full_preflight() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_dataset_snapshot();
        let dataset = dataset();
        let reference = &dataset.rows()[0];
        let fixture =
            ResolvedArgonneStanjanChonFixture::resolve_offline(local_repository(), reference)
                .expect("the exact 15-species NASA-gas universe must resolve offline at 2500 K");

        assert_eq!(fixture.resolved().layout().component_count(), 15);
        assert_eq!(fixture.preflight().len(), STANJAN_LOCAL_GAS_SPECIES.len());
        assert!(
            fixture
                .preflight()
                .iter()
                .all(|row| row.library == "NASA_gas" && row.supports_target_temperature)
        );
        assert!(
            fixture
                .preflight()
                .iter()
                .all(|row| row.gibbs_j_mol.is_finite())
        );
        assert_eq!(fixture.structure().component_count, 15);
        assert_eq!(fixture.structure().element_count, 4);
        assert_eq!(fixture.structure().element_rank, 4);
        assert_eq!(fixture.structure().reaction_dimension, 11);
        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_dataset_snapshot());
    }

    #[test]
    #[ignore = "release I5 Argonne/STANJAN characterization with conservative software-regression envelope"]
    fn i5_argonne_stanjan_chon_fixed_pt_equilibrium_characterization() {
        let before_local = local_library_snapshot();
        let before_frozen = frozen_dataset_snapshot();
        let dataset = dataset();
        let reference = &dataset.rows()[0];
        let fixture =
            ResolvedArgonneStanjanChonFixture::resolve_offline(local_repository(), reference)
                .expect("the exact 15-species NASA-gas universe must resolve offline");
        let solution = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
            fixture.resolved(),
            stanjan_conditions(reference).expect("published conditions must be valid"),
            fixture
                .initial_composition()
                .expect("equivalent local feed must map to layout"),
        ))
        .expect("production fixed-P,T general solver must accept the CHON benchmark");
        let report = StanjanChonComparisonReport::from_solution(&dataset, &fixture, &solution)
            .expect("accepted local solution must align to the frozen external identities");
        render_comparison(&report);

        assert_eq!(report.species.len(), STANJAN_DECLARED_SPECIES.len());
        assert_eq!(report.phase_control_transitions, 0);
        assert!(report.kithe_sum.is_finite() && (report.kithe_sum - 1.0).abs() <= 1.0e-10);
        assert!(report.residual_l2_norm.is_finite());
        assert!(report.max_abs_element_balance_error.is_finite());
        assert_stanjan_external_regression_envelope(&report);
        let pentane = report
            .species
            .iter()
            .find(|row| row.external_identity == "C5H12")
            .expect("report must retain the source pentane row");
        assert_eq!(
            pentane.magnitude_class,
            StanjanChonMagnitudeClass::ExternallyZeroReactantNotSolved
        );
        assert!(pentane.kithe_mole_fraction.is_none());
        assert!(pentane.relative_error.is_none());
        assert_eq!(before_local, local_library_snapshot());
        assert_eq!(before_frozen, frozen_dataset_snapshot());
    }
}
