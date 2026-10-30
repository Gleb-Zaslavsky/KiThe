#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::Thermodynamics::DBhandlers::NISTdata::NISTdata;

    fn fixture_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nist_webbook")
            .join(name)
    }

    #[test]
    fn local_nist_water_shomate_fixture_matches_the_published_500k_anchor_offline() {
        let value: serde_json::Value = serde_json::from_slice(
            &std::fs::read(fixture_path("water_gas_shomate.json"))
                .expect("frozen NIST Shomate payload must be readable"),
        )
        .expect("frozen NIST Shomate payload must be valid JSON");
        let metadata: serde_json::Value = serde_json::from_slice(
            &std::fs::read(fixture_path("water_gas_shomate.metadata.json"))
                .expect("frozen NIST source metadata must be readable"),
        )
        .expect("frozen NIST source metadata must be valid JSON");

        assert_eq!(
            metadata["fixture_id"],
            "nist-webbook-water-gas-shomate-chase-1998-v1"
        );
        assert_eq!(metadata["record"]["cas_registry_number"], "7732-18-5");

        let mut nist = NISTdata::new();
        nist.set_unit("J").unwrap();
        nist.from_serde(value).unwrap();
        nist.extract_coefficients(500.0).unwrap();
        nist.calculate_cp_dh_ds(500.0).unwrap();

        let anchor = &metadata["reviewed_500_kelvin_anchor"];
        let cp = anchor["cp_j_mol_k"].as_f64().unwrap();
        let enthalpy_increment = anchor["enthalpy_increment_kj_mol"].as_f64().unwrap();
        let entropy = anchor["entropy_j_mol_k"].as_f64().unwrap();
        assert!((nist.Cp - cp).abs() <= 0.02, "Cp={} J/mol/K", nist.Cp);
        assert!(
            ((nist.dh - (-241.8264 * 1000.0)) / 1000.0 - enthalpy_increment).abs() <= 0.02,
            "H-H298={} kJ/mol",
            (nist.dh - (-241.8264 * 1000.0)) / 1000.0
        );
        assert!((nist.ds - entropy).abs() <= 0.1, "S={} J/mol/K", nist.ds);
    }
}
