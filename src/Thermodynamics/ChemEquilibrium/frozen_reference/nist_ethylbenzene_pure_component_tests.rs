#[cfg(test)]
mod tests {
    use std::fs;

    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_ethylbenzene_pure_component::{
        load_nist_ethylbenzene_pure_component_seed, nist_ethylbenzene_antoine_pressure_pa,
        nist_ethylbenzene_majer_svoboda_vaporization_enthalpy_j_mol,
    };

    #[test]
    fn nist_ethylbenzene_seed_preserves_approved_standard_state_and_phase_change_anchors() {
        let dataset = load_nist_ethylbenzene_pure_component_seed().unwrap();
        assert!(dataset.is_external_evidence());
        assert_eq!(dataset.rows().len(), 1);
        let row = &dataset.rows()[0];

        assert_eq!(row.compound, "ethylbenzene");
        assert_eq!(row.formula, "C8H10");
        assert_eq!(row.cas_registry_number, "100-41-4");
        assert_eq!(row.reference_temperature_k, 298.15);
        assert_eq!(row.liquid_formation_enthalpy_j_mol, -12_500.0);
        assert_eq!(row.liquid_standard_entropy_j_mol_k, 255.01);
        assert_eq!(row.gas_formation_enthalpy_j_mol, 29_800.0);
        assert_eq!(row.liquid_heat_capacity_anchors.len(), 6);
        assert_eq!(row.vaporization_enthalpy_points.len(), 11);
        assert_eq!(row.vaporization_enthalpy_points[0].temperature_k, 294.01);
        assert_eq!(
            row.vaporization_enthalpy_points[0].vaporization_enthalpy_j_mol,
            42_490.0
        );

        // This is a source-anchor consistency check, not a fitted closure or
        // an external VLE acceptance envelope.
        let formation_difference =
            row.gas_formation_enthalpy_j_mol - row.liquid_formation_enthalpy_j_mol;
        assert!((formation_difference - 42_300.0).abs() < 1.0e-12);
        assert!(
            (formation_difference
                - row.vaporization_enthalpy_points[0].vaporization_enthalpy_j_mol)
                .abs()
                < 500.0
        );
    }

    #[test]
    fn pure_component_oracles_are_bounded_finite_and_monotonic_over_ternary_review_window() {
        let dataset = load_nist_ethylbenzene_pure_component_seed().unwrap();
        let row = &dataset.rows()[0];
        let temperatures = [340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0];
        let mut previous_pressure = 0.0_f64;
        for temperature in temperatures {
            let pressure = nist_ethylbenzene_antoine_pressure_pa(row.antoine, temperature).unwrap();
            let vaporization = nist_ethylbenzene_majer_svoboda_vaporization_enthalpy_j_mol(
                row.majer_svoboda,
                temperature,
            )
            .unwrap();
            assert!(
                pressure > previous_pressure,
                "Antoine pressure must rise with T at {temperature} K"
            );
            assert!(vaporization.is_finite() && vaporization > 0.0);
            previous_pressure = pressure;
        }
        assert!(nist_ethylbenzene_antoine_pressure_pa(row.antoine, 329.73).is_err());
        assert!(nist_ethylbenzene_antoine_pressure_pa(row.antoine, 410.28).is_err());
        assert!(
            nist_ethylbenzene_majer_svoboda_vaporization_enthalpy_j_mol(row.majer_svoboda, 294.99)
                .is_err()
        );
        assert!(
            nist_ethylbenzene_majer_svoboda_vaporization_enthalpy_j_mol(row.majer_svoboda, 437.01)
                .is_err()
        );
    }

    #[test]
    fn loading_nist_ethylbenzene_seed_never_rewrites_frozen_source_files() {
        let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("nist_webbook");
        let metadata = directory.join("ethylbenzene_pure_component_seed.metadata.json");
        let rows = directory.join("ethylbenzene_pure_component_seed.rows.json");
        let metadata_before = fs::read(&metadata).unwrap();
        let rows_before = fs::read(&rows).unwrap();

        load_nist_ethylbenzene_pure_component_seed().unwrap();

        assert_eq!(fs::read(metadata).unwrap(), metadata_before);
        assert_eq!(fs::read(rows).unwrap(), rows_before);
    }
}
