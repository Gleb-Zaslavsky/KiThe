//! Frozen NIST evidence for the low-temperature extension of the canonical CO2(g) state.
//!
//! This is intentionally not a production CO2 adapter yet. The NIST table has
//! thermal ideal-gas functions, while the existing KiThe record also encodes a
//! formation-state datum. The future splice must establish that datum once and
//! preserve it across the full gas/solid closure.

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::path::PathBuf;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::ResolvedThermochemistry;
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenReferenceDataset, NistCo2IdealGasReference,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemFactory, SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    use crate::Thermodynamics::User_substances::Phases;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;

    const SOURCE_GAS_CONSTANT_J_MOL_K: f64 = 8.314_459_8;
    const SOURCE_STANDARD_PRESSURE_PA: f64 = 100_000.0;
    const JANAF_100_K_CP_J_MOL_K: f64 = 29.208;
    const JANAF_100_K_ENTROPY_J_MOL_K: f64 = 179.009;
    const JANAF_100_K_H_MINUS_298_15_J_MOL: f64 = -6_456.0;
    const JANAF_200_K_CP_J_MOL_K: f64 = 32.359;
    const JANAF_200_K_ENTROPY_J_MOL_K: f64 = 199.975;
    const JANAF_200_K_H_MINUS_298_15_J_MOL: f64 = -3_414.0;

    fn dataset() -> FrozenReferenceDataset<NistCo2IdealGasReference> {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nist_tashkun_harvey");
        FrozenReferenceDataset::load(
            root.join("co2_ideal_gas_low_temperature.metadata.json"),
            root.join("co2_ideal_gas_low_temperature.rows.json"),
        )
        .expect("frozen Tashkun-Harvey CO2 evidence must load")
    }

    fn local_co2_thermochemistry() -> ResolvedThermochemistry {
        let specification =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
                "CO2".to_owned(),
            ]))
            .with_phase_natures(Some(HashMap::from([("gas".to_owned(), Phases::Gas)])))
            .with_library_priorities(vec!["NASA_gas".to_owned()])
            .with_search_in_nist(false)
            .build()
            .expect("CO2 gas-only source specification must be valid");
        let repository = ThermoData::try_default_repository()
            .expect("bundled offline thermochemistry repository must be available");
        let resolved =
            SubstanceSystemFactory::resolve_phase_system_with_repository(specification, repository)
                .expect("canonical local CO2(g) record must resolve offline");
        ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("local CO2(g) must expose Gibbs, enthalpy, and Cp capabilities")
    }

    #[test]
    fn tashkun_harvey_low_temperature_rows_are_source_faithful_and_ordered() {
        let dataset = dataset();
        assert!(dataset.is_external_evidence());
        assert_eq!(
            dataset.rows().first().map(|row| row.temperature_k),
            Some(100.0)
        );
        assert_eq!(
            dataset.rows().last().map(|row| row.temperature_k),
            Some(298.0)
        );
        assert!(dataset.rows().iter().any(|row| row.temperature_k == 200.0));
        assert!(dataset.rows().iter().any(|row| row.temperature_k == 195.0));
        assert!(
            dataset
                .metadata()
                .source
                .stable_identifier
                .as_deref()
                .is_some_and(|value| value.contains("4e0e057bfd3b92dc"))
        );
    }

    #[test]
    fn tashkun_harvey_conversions_are_explicit_and_physical() {
        let dataset = dataset();
        let row = dataset
            .rows()
            .iter()
            .find(|row| row.temperature_k == 200.0)
            .expect("frozen source must retain the 200 K splice anchor");
        let cp = row.heat_capacity_over_r * SOURCE_GAS_CONSTANT_J_MOL_K;
        let entropy = row.entropy_over_r * SOURCE_GAS_CONSTANT_J_MOL_K;
        let enthalpy = row.enthalpy_over_rt * SOURCE_GAS_CONSTANT_J_MOL_K * row.temperature_k;

        assert!((cp - 32.354_896_930_746).abs() < 1e-9);
        assert!((entropy - 200.067_573_596_023).abs() < 1e-9);
        assert!((enthalpy - 5_952.126_381_014_7).abs() < 1e-9);
        assert_eq!(SOURCE_STANDARD_PRESSURE_PA, 100_000.0);
    }

    /// JANAF is independent, older evidence. These deliberately compare
    /// thermal enthalpy differences, avoiding an assumption about either
    /// source's absolute formation-enthalpy zero.
    #[test]
    fn tashkun_harvey_characterizes_independent_janaf_low_temperature_anchors() {
        let dataset = dataset();
        let row_at = |temperature_k| {
            dataset
                .rows()
                .iter()
                .find(|row| row.temperature_k == temperature_k)
                .copied()
                .expect("frozen source must retain the requested JANAF anchor")
        };
        let row_100 = row_at(100.0);
        let row_200 = row_at(200.0);
        let row_298 = row_at(298.0);
        let values = |row: NistCo2IdealGasReference| {
            (
                row.heat_capacity_over_r * SOURCE_GAS_CONSTANT_J_MOL_K,
                row.entropy_over_r * SOURCE_GAS_CONSTANT_J_MOL_K,
                row.enthalpy_over_rt * SOURCE_GAS_CONSTANT_J_MOL_K * row.temperature_k,
            )
        };
        let (cp_100, s_100, h_100) = values(row_100);
        let (cp_200, s_200, h_200) = values(row_200);
        let (_, _, h_298) = values(row_298);

        println!("NIST-2025 versus independent JANAF low-temperature anchors");
        println!("  T K    dCp J/mol/K    dS J/mol/K    d(H-H298) J/mol");
        println!(
            "  100    {:+.6}       {:+.6}       {:+.3}",
            cp_100 - JANAF_100_K_CP_J_MOL_K,
            s_100 - JANAF_100_K_ENTROPY_J_MOL_K,
            (h_100 - h_298) - JANAF_100_K_H_MINUS_298_15_J_MOL,
        );
        println!(
            "  200    {:+.6}       {:+.6}       {:+.3}",
            cp_200 - JANAF_200_K_CP_J_MOL_K,
            s_200 - JANAF_200_K_ENTROPY_J_MOL_K,
            (h_200 - h_298) - JANAF_200_K_H_MINUS_298_15_J_MOL,
        );

        assert!((cp_100 - JANAF_100_K_CP_J_MOL_K).abs() < 0.02);
        assert!((s_100 - JANAF_100_K_ENTROPY_J_MOL_K).abs() < 0.2);
        assert!(((h_100 - h_298) - JANAF_100_K_H_MINUS_298_15_J_MOL).abs() < 50.0);
        assert!((cp_200 - JANAF_200_K_CP_J_MOL_K).abs() < 0.02);
        assert!((s_200 - JANAF_200_K_ENTROPY_J_MOL_K).abs() < 0.2);
        assert!(((h_200 - h_298) - JANAF_200_K_H_MINUS_298_15_J_MOL).abs() < 50.0);
    }

    /// Establishes the numerical evidence needed before a future gas splice.
    /// This is intentionally characterization-only: the Tashkun--Harvey table
    /// supplies thermal ideal-gas functions while local NASA thermochemistry
    /// carries an absolute formation-state convention that has not yet been
    /// joined and reviewed.
    #[test]
    fn local_nasa_co2_join_boundary_is_characterized_without_second_gas_state() {
        let frozen = dataset();
        let local = local_co2_thermochemistry();
        let source = |temperature_k| {
            frozen
                .rows()
                .iter()
                .find(|row| row.temperature_k == temperature_k)
                .copied()
                .expect("frozen source must retain a join characterization row")
        };

        println!("CO2(g) low-temperature source versus canonical NASA join characterization");
        println!("  T K    dCp J/mol/K   dS J/mol/K    dH(T)-dH(200) J/mol");
        for temperature_k in [200.0, 201.0, 205.0, 298.0] {
            let source_row = source(temperature_k);
            let source_cp = source_row.heat_capacity_over_r * SOURCE_GAS_CONSTANT_J_MOL_K;
            let source_entropy = source_row.entropy_over_r * SOURCE_GAS_CONSTANT_J_MOL_K;
            let source_h = source_row.enthalpy_over_rt
                * SOURCE_GAS_CONSTANT_J_MOL_K
                * source_row.temperature_k;
            let local_cp = local
                .evaluate_heat_capacity(temperature_k)
                .expect("local source temperature must be in range")[0]
                .expect("CO2(g) local record must expose Cp");
            let local_h = local
                .evaluate_enthalpy(temperature_k)
                .expect("local source temperature must be in range")[0];
            let local_g = local
                .evaluate_gibbs(temperature_k)
                .expect("local source temperature must be in range")[0];
            let local_entropy = (local_h - local_g) / temperature_k;
            let source_h_200 = source(200.0).enthalpy_over_rt * SOURCE_GAS_CONSTANT_J_MOL_K * 200.0;
            let local_h_200 = local
                .evaluate_enthalpy(200.0)
                .expect("local 200 K source temperature must be in range")[0];
            println!(
                "  {temperature_k:>3.0}    {:+.6}      {:+.6}      {:+.3}",
                local_cp - source_cp,
                local_entropy - source_entropy,
                (local_h - local_h_200) - (source_h - source_h_200),
            );
        }
    }
}
