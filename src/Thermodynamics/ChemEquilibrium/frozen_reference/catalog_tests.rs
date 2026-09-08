//! Completeness checks for the read-only frozen-reference registry.
//!
//! Typed adapters validate the physical meaning of their own rows.  This
//! catalog is deliberately orthogonal: it makes the repository-level contract
//! explicit so a new JSON pair cannot be silently orphaned, duplicated, or
//! mislabeled while every individual adapter test still happens to pass.

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::fs;
    use std::path::{Path, PathBuf};

    use serde::Deserialize;
    use serde_json::Value;

    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FROZEN_REFERENCE_FORMAT_VERSION, FrozenReferenceEvidenceKind, FrozenReferenceMetadata,
    };

    #[derive(Debug, Clone, Copy)]
    struct FrozenReferenceCatalogEntry {
        dataset_id: &'static str,
        evidence_kind: FrozenReferenceEvidenceKind,
        metadata_relative: &'static str,
        expected_row_count: usize,
        /// Owning typed adapter, retained as reviewable provenance for the
        /// registry rather than inferred from an arbitrary filename.
        adapter: &'static str,
    }

    #[derive(Debug, Deserialize)]
    #[serde(deny_unknown_fields)]
    struct RawRowsFile {
        dataset_id: String,
        rows: Vec<Value>,
    }

    const CATALOG: [FrozenReferenceCatalogEntry; 14] = [
        FrozenReferenceCatalogEntry {
            dataset_id: "argonne.stanjan.pentane_methane_air.tp.2500k_35atm.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "argonne_stanjan/pentane_methane_air_tp_2500k_35atm.metadata.json",
            expected_row_count: 1,
            adapter: "argonne_stanjan_chon",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "iapws.water_ice_ih_sublimation.low_pressure.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "iapws/water_ice_sublimation_low_pressure.metadata.json",
            expected_row_count: 8,
            adapter: "iapws_ice_sublimation",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "iapws.water_liquid_saturation.low_pressure.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "iapws/water_liquid_saturation_low_pressure.metadata.json",
            expected_row_count: 10,
            adapter: "iapws_water_saturation",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "janaf.boudouard.reaction_thermodynamics.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "janaf/boudouard_reaction_thermodynamics.metadata.json",
            expected_row_count: 11,
            adapter: "janaf_boudouard_thermochemistry",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nasa_cea.h2_o2.hp.scitech_2025.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nasa_cea/h2_o2_hp_scitech_2025.metadata.json",
            expected_row_count: 1,
            adapter: "nasa_cea_h2_o2_hp",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nasa_cea.rp1311.example3.hydrocarbon_air.hp.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nasa_cea/rp1311_example3_hp.metadata.json",
            expected_row_count: 3,
            adapter: "nasa_cea_rp1311_example3",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nasa_tp1906.chon_graphite.er125.1atm.heterogeneous_enthalpy.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nasa_tp1906/chon_graphite_er125_1atm_heterogeneous_enthalpy.metadata.json",
            expected_row_count: 4,
            adapter: "nasa_tp1906_chon_graphite_ph",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nasa_tp1907.chon_graphite.er125.1atm.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nasa_tp1907/chon_graphite_er125_1atm.metadata.json",
            expected_row_count: 4,
            adapter: "nasa_tp1907_chon_graphite",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nist.tashkun_harvey.co2_ideal_gas_low_temperature.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nist_tashkun_harvey/co2_ideal_gas_low_temperature.metadata.json",
            expected_row_count: 14,
            adapter: "nist_tashkun_harvey_co2",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nist.thermoml.benzene_toluene.vle.353_15k.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nist_thermoml/benzene_toluene_vle_353_15k.metadata.json",
            expected_row_count: 5,
            adapter: "nist_benzene_toluene_vle",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nist.thermoml.toluene_ethylbenzene_chlorobenzene.selected_interior.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nist_thermoml/toluene_ethylbenzene_chlorobenzene_selected.metadata.json",
            expected_row_count: 4,
            adapter: "nist_toluene_ethylbenzene_chlorobenzene_preflight",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nist.webbook.ethylbenzene.pure_component_seed.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nist_webbook/ethylbenzene_pure_component_seed.metadata.json",
            expected_row_count: 1,
            adapter: "nist_ethylbenzene_pure_component",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "nist.webbook.ternary_vle_antoine_gauge.v1",
            evidence_kind: FrozenReferenceEvidenceKind::FrozenExternal,
            metadata_relative: "nist_webbook/toluene_ethylbenzene_chlorobenzene_antoine.metadata.json",
            expected_row_count: 3,
            adapter: "nist_ternary_vle_antoine_gauge",
        },
        FrozenReferenceCatalogEntry {
            dataset_id: "synthetic.temperature_pressure.v1",
            evidence_kind: FrozenReferenceEvidenceKind::SyntheticInfrastructure,
            metadata_relative: "synthetic/temperature_pressure.metadata.json",
            expected_row_count: 3,
            adapter: "loader_tests",
        },
    ];

    fn data_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data")
    }

    fn collect_relative_files(root: &Path, suffix: &str) -> Vec<String> {
        fn visit(directory: &Path, root: &Path, suffix: &str, found: &mut Vec<String>) {
            for entry in fs::read_dir(directory).expect("frozen data directory must be readable") {
                let entry = entry.expect("frozen data directory entry must be readable");
                let path = entry.path();
                if path.is_dir() {
                    visit(&path, root, suffix, found);
                } else if path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.ends_with(suffix))
                {
                    found.push(
                        path.strip_prefix(root)
                            .expect("catalog path must stay below frozen data root")
                            .to_string_lossy()
                            .replace('\\', "/"),
                    );
                }
            }
        }

        let mut found = Vec::new();
        visit(root, root, suffix, &mut found);
        found.sort();
        found
    }

    #[test]
    fn catalog_covers_every_frozen_pair_with_exact_identity_kind_version_and_row_count() {
        let root = data_root();
        let mut identifiers = HashSet::new();
        let mut metadata_paths = HashSet::new();
        let mut row_paths = HashSet::new();

        for entry in CATALOG {
            assert!(
                identifiers.insert(entry.dataset_id),
                "catalog contains duplicate dataset id '{}'",
                entry.dataset_id
            );
            assert!(
                metadata_paths.insert(entry.metadata_relative),
                "catalog contains duplicate metadata path '{}'",
                entry.metadata_relative
            );
            assert!(
                !entry.adapter.trim().is_empty(),
                "catalog entry '{}' must name its typed adapter",
                entry.dataset_id
            );

            let metadata_path = root.join(entry.metadata_relative);
            let metadata: FrozenReferenceMetadata = serde_json::from_slice(
                &fs::read(&metadata_path).expect("catalog metadata file must be readable"),
            )
            .expect("catalog metadata must deserialize");
            assert_eq!(
                metadata.dataset_format_version,
                FROZEN_REFERENCE_FORMAT_VERSION
            );
            assert_eq!(metadata.dataset_id, entry.dataset_id);
            assert_eq!(metadata.evidence_kind, entry.evidence_kind);

            let data_path = metadata_path
                .parent()
                .expect("metadata file must have a parent directory")
                .join(&metadata.data_file);
            let data_relative = data_path
                .strip_prefix(&root)
                .expect("catalog rows path must stay below frozen data root")
                .to_string_lossy()
                .replace('\\', "/");
            assert!(
                row_paths.insert(data_relative.clone()),
                "catalog contains duplicate rows path '{data_relative}'"
            );
            let rows: RawRowsFile = serde_json::from_slice(
                &fs::read(data_path).expect("catalog rows file must be readable"),
            )
            .expect("catalog rows file must deserialize");
            assert_eq!(rows.dataset_id, entry.dataset_id);
            assert_eq!(rows.rows.len(), entry.expected_row_count);
        }

        let mut catalog_metadata_paths = metadata_paths.into_iter().collect::<Vec<_>>();
        catalog_metadata_paths.sort_unstable();
        assert_eq!(
            catalog_metadata_paths,
            collect_relative_files(&root, ".metadata.json"),
            "every metadata file must have exactly one catalog entry"
        );
        let mut catalog_row_paths = row_paths.into_iter().collect::<Vec<_>>();
        catalog_row_paths.sort_unstable();
        assert_eq!(
            catalog_row_paths,
            collect_relative_files(&root, ".rows.json"),
            "every rows file must be the declared payload of exactly one catalog entry"
        );
    }
}
