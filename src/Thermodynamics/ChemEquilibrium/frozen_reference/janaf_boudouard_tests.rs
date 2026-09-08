//! First Boudouard I4/I5 evidence: local standard thermochemistry versus
//! frozen primary NIST-JANAF species quantities.
//!
//! This test module must remain free of phase-control calls. A disagreement in
//! standard reaction Gibbs energy belongs to thermochemistry/provenance review
//! first; adding TPD or an active-set solve would obscure that diagnosis.

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenReferenceDataset, JanafBoudouardReference,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_thermochemistry::{
        JanafBoudouardComparisonContract, JanafBoudouardDerived,
        JanafBoudouardThermochemistryReport, JanafBoudouardThermochemistryRow,
        JANAF_STANDARD_PRESSURE_PA,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
        PurePhaseBoundaryStructuralTolerances,
    };
    use crate::Thermodynamics::ChemEquilibrium::real_pure_phase_fixtures::{
        RealPurePhaseFamily, RealPurePhaseGasScenario, ResolvedRealPurePhaseFixture,
    };
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};

    // The JANAF log-Kf columns are published at three decimal places. This
    // checks their arithmetic consistency with Gibbs-derived log10(Kp) while
    // leaving room for that documented source rounding.
    const MAX_JANAF_INTERNAL_LOG10_ROUTE_DELTA: f64 = 0.005;
    // Regression guards selected after repeated debug/release characterization.
    // They are not JANAF source uncertainty or physical acceptance tolerances.
    const MAX_EXTERNAL_DELTA_G_ERROR_J_MOL_GUARD: f64 = 100.0;
    const MAX_EXTERNAL_LOG10_K_ERROR_GUARD: f64 = 0.01;

    fn janaf_directory() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("janaf")
    }

    fn janaf_dataset() -> FrozenReferenceDataset<JanafBoudouardReference> {
        let directory = janaf_directory();
        FrozenReferenceDataset::load(
            directory.join("boudouard_reaction_thermodynamics.metadata.json"),
            directory.join("boudouard_reaction_thermodynamics.rows.json"),
        )
        .expect("reviewed frozen JANAF Boudouard dataset must load")
    }

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("the bundled offline thermochemistry repository must be available")
    }

    fn boudouard_fixture() -> ResolvedRealPurePhaseFixture {
        RealPurePhaseFamily::BoudouardCarbon
            .resolve_offline(local_repository(), &[])
            .expect("the pinned local Boudouard fixture must resolve without NIST")
    }

    fn component_index(fixture: &ResolvedRealPurePhaseFixture, expected_label: &str) -> usize {
        fixture
            .resolved()
            .layout()
            .component_labels()
            .iter()
            .position(|label| label == expected_label)
            .unwrap_or_else(|| panic!("Boudouard layout lacks '{expected_label}'"))
    }

    fn local_reaction_thermochemistry(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature_k: f64,
    ) -> (f64, f64) {
        let co = component_index(fixture, "gas::CO");
        let co2 = component_index(fixture, "gas::CO2");
        let graphite = component_index(fixture, "solid::C(gr)");
        let g = fixture
            .thermochemistry()
            .evaluate_gibbs(temperature_k)
            .unwrap_or_else(|error| {
                panic!("local Boudouard G0 failed at {temperature_k} K: {error}")
            });
        let delta_r_g_j_mol = g[co2] + g[graphite] - 2.0 * g[co];
        let log10_kp = -delta_r_g_j_mol
            / (crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::MOLAR_GAS_CONSTANT
                * temperature_k * std::f64::consts::LN_10);
        (delta_r_g_j_mol, log10_kp)
    }

    fn build_report(
        fixture: &ResolvedRealPurePhaseFixture,
        dataset: &FrozenReferenceDataset<JanafBoudouardReference>,
    ) -> JanafBoudouardThermochemistryReport {
        let bounds = fixture.thermochemistry().temperature_bounds();
        let unsupported = dataset
            .rows()
            .iter()
            .filter(|row| !bounds.contains(row.temperature_k))
            .map(|row| row.temperature_k)
            .collect::<Vec<_>>();
        assert!(
            unsupported.is_empty(),
            "local Boudouard bounds [{}, {}] K exclude frozen JANAF rows {unsupported:?}; do not extrapolate local thermochemistry",
            bounds.lower(),
            bounds.upper(),
        );

        let rows = dataset
            .rows()
            .iter()
            .map(|reference| {
                let (kithe_delta_r_g_j_mol, kithe_log10_kp) =
                    local_reaction_thermochemistry(fixture, reference.temperature_k);
                JanafBoudouardThermochemistryRow {
                    temperature_k: reference.temperature_k,
                    janaf: JanafBoudouardDerived::from_reference(*reference),
                    kithe_delta_r_g_j_mol,
                    kithe_log10_kp,
                }
            })
            .collect();
        JanafBoudouardThermochemistryReport::new(dataset, rows, JANAF_STANDARD_PRESSURE_PA)
            .expect("finite local and frozen Boudouard evidence must build a report")
    }

    fn assert_janaf_external_regression_envelope(report: &JanafBoudouardThermochemistryReport) {
        assert!(
            report.max_delta_g_error_j_mol() <= MAX_EXTERNAL_DELTA_G_ERROR_J_MOL_GUARD,
            "JANAF Boudouard external quality regression: max |delta G|={:e} J/mol exceeds reviewed guard={MAX_EXTERNAL_DELTA_G_ERROR_J_MOL_GUARD:e} J/mol",
            report.max_delta_g_error_j_mol(),
        );
        assert!(
            report.max_log10_k_error() <= MAX_EXTERNAL_LOG10_K_ERROR_GUARD,
            "JANAF Boudouard external quality regression: max |delta log10 K|={:e} exceeds reviewed guard={MAX_EXTERNAL_LOG10_K_ERROR_GUARD:e}",
            report.max_log10_k_error(),
        );
    }

    #[test]
    fn i5_janaf_boudouard_dataset_and_local_fixture_contract_are_valid() {
        let dataset = janaf_dataset();
        let fixture = boudouard_fixture();
        let contract = JanafBoudouardComparisonContract::characterization_only();
        contract.validate_dataset(&dataset).unwrap();
        assert_eq!(contract.standard_pressure_pa(), 100_000.0);
        assert_eq!(dataset.rows().len(), 11);
        assert_eq!(dataset.rows().first().unwrap().temperature_k, 500.0);
        assert_eq!(dataset.rows().last().unwrap().temperature_k, 1_500.0);
        assert_eq!(
            fixture.resolved().layout().component_labels(),
            ["gas::CO", "gas::CO2", "solid::C(gr)"]
        );

        let provenance = fixture.thermochemistry().provenance();
        assert_eq!(provenance.len(), 3);
        assert_eq!(provenance[0].library(), "NASA_gas");
        assert_eq!(provenance[0].record_key(), "CO");
        assert_eq!(provenance[1].library(), "NASA_gas");
        assert_eq!(provenance[1].record_key(), "CO2");
        assert_eq!(provenance[2].library(), "NASA_cond");
        assert_eq!(provenance[2].record_key(), "C(gr)");
        // The bundled NASA records currently have no machine-readable p0
        // field. Keep that fact visible: this G0-only characterization applies
        // no correction and a future pressure-boundary comparison must supply
        // an explicit source convention instead of guessing one from NASA7.
        assert!(
            provenance
                .iter()
                .all(|row| !row.standard_state_pressure().is_declared())
        );
    }

    #[test]
    fn i5_janaf_boudouard_structure_has_no_gas_only_reaction() {
        let fixture = boudouard_fixture();
        let gas = RealPurePhaseGasScenario::new(vec![0.5, 0.25]).unwrap();
        let problem = fixture
            .to_pt_boundary_problem(
                &gas,
                EquilibriumConditions::new(
                    1_000.0,
                    JANAF_STANDARD_PRESSURE_PA,
                    JANAF_STANDARD_PRESSURE_PA,
                )
                .unwrap(),
            )
            .expect("Boudouard independent structure must materialize");
        let reaction_space = problem
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .expect("Boudouard fixture must have one phase-forming direction only");
        assert_eq!(reaction_space.gas_only_reaction_dimension, 0);
        assert_eq!(reaction_space.full_reaction_dimension, 1);
    }

    #[test]
    fn i5_janaf_boudouard_local_standard_gibbs_is_finite_on_full_frozen_grid() {
        let dataset = janaf_dataset();
        let fixture = boudouard_fixture();
        let report = build_report(&fixture, &dataset);
        assert_eq!(report.rows().len(), dataset.rows().len());
        report
            .validate_janaf_internal_routes(MAX_JANAF_INTERNAL_LOG10_ROUTE_DELTA)
            .expect("rounded JANAF Gibbs and log-Kf routes must remain internally consistent");
        assert_janaf_external_regression_envelope(&report);
    }

    #[test]
    #[ignore = "release I5 JANAF Boudouard characterization with conservative software-regression envelope"]
    fn i5_janaf_boudouard_reaction_thermochemistry_diagnostic() {
        let directory = janaf_directory();
        let metadata_path = directory.join("boudouard_reaction_thermodynamics.metadata.json");
        let rows_path = directory.join("boudouard_reaction_thermodynamics.rows.json");
        let metadata_before = fs::read(&metadata_path).unwrap();
        let rows_before = fs::read(&rows_path).unwrap();
        let dataset = janaf_dataset();
        let contract = JanafBoudouardComparisonContract::characterization_only();
        contract.validate_dataset(&dataset).unwrap();
        let fixture = boudouard_fixture();
        let report = build_report(&fixture, &dataset);
        contract.validate_report(&report).unwrap();
        report
            .validate_janaf_internal_routes(MAX_JANAF_INTERNAL_LOG10_ROUTE_DELTA)
            .unwrap_or_else(|error| panic!("frozen JANAF route mismatch: {error}"));
        assert_janaf_external_regression_envelope(&report);
        println!("KiThe standard-state pressure provenance:");
        for row in fixture.thermochemistry().provenance() {
            println!(
                "  {:<16} {}",
                row.component().label(),
                row.standard_state_pressure()
            );
        }
        println!("{report}");
        assert_eq!(fs::read(metadata_path).unwrap(), metadata_before);
        assert_eq!(fs::read(rows_path).unwrap(), rows_before);
    }
}
