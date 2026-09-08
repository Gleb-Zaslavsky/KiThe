//! Read-only infrastructure for frozen external reference evidence.
//!
//! The validation ladder used by the equilibrium suite distinguishes:
//! - I1: independent thermodynamic equations;
//! - I2: an independent numerical route;
//! - I3: the production phase lifecycle;
//! - I4: offline thermochemistry resolved from KiThe's local repository;
//! - I5: authoritative external numbers frozen in the test repository.
//!
//! I5 does not imply I1: a published equilibrium table may only be compared
//! with production output. I1 does not imply I5: exact synthetic mathematics
//! is not evidence of agreement with external reference data.
//!
//! This module intentionally has no writer, downloader, cache, tolerance, or
//! thermochemical-format adapter. A frozen dataset is small, read-only test
//! evidence. NASA/NIST handlers remain responsible for KiThe source records;
//! this loader reads already-published numerical reference values.

use serde::Deserialize;
use serde::de::DeserializeOwned;
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use thiserror::Error;

/// Shared assertions used only by frozen-reference regression modules.
///
/// They intentionally stay out of production code: the contracts are test
/// evidence, while physical topology and external-source comparisons remain
/// explicit in each fixture.
#[cfg(test)]
pub(crate) mod assertions;

/// The only frozen-reference metadata schema understood by this build.
///
/// A future schema must be introduced deliberately with an explicit migration
/// reader. Silently accepting a newer positive version would make tests claim
/// provenance they do not actually understand.
pub(crate) const FROZEN_REFERENCE_FORMAT_VERSION: u32 = 1;

/// Distinguishes real I5 evidence from fixtures that only test this loader.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum FrozenReferenceEvidenceKind {
    SyntheticInfrastructure,
    FrozenExternal,
}

/// Optional uncertainty explicitly stated by the source.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ReferenceUncertainty {
    /// Numerical magnitude of the uncertainty.
    pub(crate) value: f64,
    /// Unit of the uncertainty magnitude.
    pub(crate) unit: String,
}

/// Meaning and unit of one typed numerical field.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ReferenceQuantity {
    /// Column name as it appears in the rows JSON.
    pub(crate) column: String,
    /// Human-readable physical meaning of the column.
    pub(crate) meaning: String,
    /// Unit expressed for this column.
    pub(crate) unit: String,
    /// Optional source-stated printing precision.
    #[serde(default)]
    pub(crate) source_precision: Option<String>,
    /// Optional uncertainty explicitly stated by the source.
    #[serde(default)]
    pub(crate) uncertainty: Option<ReferenceUncertainty>,
}

/// Identity of the publication or correlation from which values were copied.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenReferenceSource {
    /// Publishing organization (e.g. NASA, NIST, IAPWS, JANAF).
    pub(crate) organization: String,
    /// Publication or correlation name.
    pub(crate) name: String,
    /// Optional release/version tag of the publication.
    #[serde(default)]
    pub(crate) version: Option<String>,
    /// Bibliographic citation for the copied values.
    pub(crate) citation: String,
    /// Optional table or figure identifier within the source.
    #[serde(default)]
    pub(crate) source_table: Option<String>,
    /// Optional stable identifier (DOI, dataset id) for the source.
    #[serde(default)]
    pub(crate) stable_identifier: Option<String>,
}

/// Provenance stored separately from numerical rows.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenReferenceMetadata {
    /// Version of KiThe's frozen-data metadata schema, not the source release.
    pub(crate) dataset_format_version: u32,
    /// Unique identifier shared between the metadata and rows files.
    pub(crate) dataset_id: String,
    /// Human-readable title of the frozen dataset.
    pub(crate) title: String,
    /// Classification: real external evidence or loader-test fixture.
    pub(crate) evidence_kind: FrozenReferenceEvidenceKind,
    /// Filename of the rows file declared by this metadata.
    pub(crate) data_file: String,
    /// Identity and citation of the original publication.
    pub(crate) source: FrozenReferenceSource,
    /// Free-text description of how the source numbers were transcribed.
    pub(crate) transcription: String,
    /// Explicitly dimensioned columns that define the numeric schema.
    pub(crate) quantities: Vec<ReferenceQuantity>,
}

/// Static expected schema supplied by a typed row implementation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReferenceColumnSchema {
    /// Expected column name in the metadata schema.
    pub(crate) column: &'static str,
    /// Expected unit for the column.
    pub(crate) unit: &'static str,
}

/// Contract implemented separately by every physical reference-table shape.
pub(crate) trait FrozenReferenceRow: DeserializeOwned {
    /// Declares the expected typed column/unit schema for this table shape.
    fn schema() -> &'static [ReferenceColumnSchema];

    /// Validates the physical-table invariants specific to this row type.
    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError>;
}

/// One published NIST ThermoML benzene/toluene P-x datum.
///
/// `liquid_benzene_mole_fraction` is deliberately the only composition field:
/// the selected ThermoML dataset is P-x evidence and does not publish an
/// experimental vapor composition. The independent Raoult layer derives its
/// own vapor composition later; it must not be smuggled into frozen I5 data.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistBenzeneTolueneVleReference {
    /// Isotherm temperature of the frozen P-x datum, K.
    pub(crate) temperature_k: f64,
    /// Liquid-phase benzene mole fraction `x_B`.
    pub(crate) liquid_benzene_mole_fraction: f64,
    /// Experimental total pressure of the binary mixture, Pa.
    pub(crate) experimental_pressure_pa: f64,
}

const NIST_BENZENE_TOLUENE_VLE_SCHEMA: [ReferenceColumnSchema; 3] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "liquid_benzene_mole_fraction",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "experimental_pressure_pa",
        unit: "Pa",
    },
];

impl FrozenReferenceRow for NistBenzeneTolueneVleReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &NIST_BENZENE_TOLUENE_VLE_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() < 5 {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "benzene/toluene VLE evidence requires both pure endpoints and at least three interior rows".to_owned(),
            });
        }
        let mut previous_x = -1.0_f64;
        for (index, row) in rows.iter().enumerate() {
            if row.temperature_k != 353.15 {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "temperature_k".to_owned(),
                    message: format!(
                        "expected the reviewed 353.15 K isotherm, got {} K",
                        row.temperature_k
                    ),
                });
            }
            if !row.liquid_benzene_mole_fraction.is_finite()
                || !(0.0..=1.0).contains(&row.liquid_benzene_mole_fraction)
                || row.liquid_benzene_mole_fraction <= previous_x
            {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "liquid_benzene_mole_fraction".to_owned(),
                    message:
                        "liquid benzene mole fractions must be finite, ordered, and within [0, 1]"
                            .to_owned(),
                });
            }
            if !row.experimental_pressure_pa.is_finite() || row.experimental_pressure_pa <= 0.0 {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "experimental_pressure_pa".to_owned(),
                    message: "experimental pressure must be finite and positive".to_owned(),
                });
            }
            previous_x = row.liquid_benzene_mole_fraction;
        }
        if rows
            .first()
            .is_none_or(|row| row.liquid_benzene_mole_fraction != 0.0)
            || rows
                .last()
                .is_none_or(|row| row.liquid_benzene_mole_fraction != 1.0)
        {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "liquid_benzene_mole_fraction".to_owned(),
                message: "frozen P-x evidence must preserve pure toluene x_B=0 and pure benzene x_B=1 endpoints".to_owned(),
            });
        }
        Ok(())
    }
}

/// One species amount printed by an external equilibrium program.
///
/// Complete-equilibrium publications often present a heterogeneous list of
/// major, minor, trace, and absent condensed components.  A semantic name is
/// therefore part of the frozen datum: positional arrays would make a later
/// comparison vulnerable to a harmless reordering of either solver layout.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenEquilibriumSpeciesAmount {
    /// Semantic identity of the species as printed by the source program.
    pub(crate) species: String,
    /// Equilibrium amount in kg-mol per kg of reactant.
    pub(crate) amount_kgmol_per_kg: f64,
}

/// One published NASA CEA `H2/O2`, `P,H` equilibrium table.
///
/// This is deliberately a single typed row rather than a temperature series.
/// The case owns a physical input convention, one final equilibrium state,
/// and a declared component universe.  The input reactant enthalpy and all
/// KiThe-side reconstructed quantities are intentionally *not* frozen here.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NasaCeaH2O2HpReference {
    /// Problem type; must be the fixed-`P,H` code `HP`.
    pub(crate) problem_type: String,
    /// Temperature at which the reactants enter, K.
    pub(crate) reactant_temperature_k: f64,
    /// Oxidizer-to-fuel mass ratio defining the feed.
    pub(crate) oxidizer_fuel_mass_ratio: f64,
    /// Fixed pressure of the case, Pa.
    pub(crate) pressure_pa: f64,
    /// Published equilibrium temperature, K.
    pub(crate) equilibrium_temperature_k: f64,
    /// Total equilibrium amount per kg of reactant, kg-mol/kg.
    pub(crate) total_amount_kgmol_per_kg: f64,
    /// Per-species equilibrium amounts in kg-mol per kg.
    pub(crate) species_amounts: Vec<FrozenEquilibriumSpeciesAmount>,
}

const NASA_CEA_H2_O2_HP_SCHEMA: [ReferenceColumnSchema; 7] = [
    ReferenceColumnSchema {
        column: "problem_type",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "reactant_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "oxidizer_fuel_mass_ratio",
        unit: "kg/kg",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
    ReferenceColumnSchema {
        column: "equilibrium_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "total_amount_kgmol_per_kg",
        unit: "kg-mol/kg",
    },
    ReferenceColumnSchema {
        column: "species_amounts",
        unit: "kg-mol/kg",
    },
];

impl FrozenReferenceRow for NasaCeaH2O2HpReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &NASA_CEA_H2_O2_HP_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() != 1 {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "the NASA CEA H2/O2 HP dataset must contain exactly one case".to_owned(),
            });
        }
        let row = &rows[0];
        if row.problem_type != "HP" {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: "problem_type".to_owned(),
                message: format!("expected HP, got '{}'", row.problem_type),
            });
        }
        for (field, value) in [
            ("reactant_temperature_k", row.reactant_temperature_k),
            ("oxidizer_fuel_mass_ratio", row.oxidizer_fuel_mass_ratio),
            ("pressure_pa", row.pressure_pa),
            ("equilibrium_temperature_k", row.equilibrium_temperature_k),
            ("total_amount_kgmol_per_kg", row.total_amount_kgmol_per_kg),
        ] {
            validate_positive_finite(dataset_id, 0, field, value)?;
        }

        const DECLARED_SPECIES: [&str; 11] = [
            "H", "H2", "H2O", "H2O2", "HO2", "O", "O2", "O3", "OH", "H2O(L)", "H2O(cr)",
        ];
        if row.species_amounts.len() != DECLARED_SPECIES.len() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: "species_amounts".to_owned(),
                message: format!(
                    "expected {} declared CEA species, got {}",
                    DECLARED_SPECIES.len(),
                    row.species_amounts.len()
                ),
            });
        }
        let mut seen = HashSet::with_capacity(row.species_amounts.len());
        for amount in &row.species_amounts {
            if !DECLARED_SPECIES.contains(&amount.species.as_str()) {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(0),
                    field: "species_amounts.species".to_owned(),
                    message: format!("undeclared CEA identity '{}'", amount.species),
                });
            }
            if !seen.insert(amount.species.as_str()) {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(0),
                    field: "species_amounts.species".to_owned(),
                    message: format!("duplicate CEA identity '{}'", amount.species),
                });
            }
            if !amount.amount_kgmol_per_kg.is_finite() || amount.amount_kgmol_per_kg < 0.0 {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(0),
                    field: "species_amounts.amount_kgmol_per_kg".to_owned(),
                    message: format!(
                        "amount for '{}' must be finite and non-negative",
                        amount.species
                    ),
                });
            }
        }
        for species in DECLARED_SPECIES {
            if !seen.contains(species) {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(0),
                    field: "species_amounts.species".to_owned(),
                    message: format!("missing declared CEA identity '{species}'"),
                });
            }
        }
        Ok(())
    }
}

/// One molecular amount used to state a published equilibrium feed.
///
/// These rows intentionally preserve the source molecular basis even when an
/// executable local fixture uses a different, element-equivalent basis.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenMolarAmount {
    /// Molecular identity of the feed species.
    pub(crate) species: String,
    /// Amount of the species in the source molecular basis, mol.
    pub(crate) amount_mol: f64,
}

/// One explicitly published elemental inventory total.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenElementAmount {
    /// Chemical element symbol.
    pub(crate) element: String,
    /// Explicit elemental inventory total, in mol-atoms.
    pub(crate) amount_mol_atoms: f64,
}

/// One published STANJAN equilibrium mole fraction.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenEquilibriumMoleFraction {
    /// Species identity as published in the STANJAN table.
    pub(crate) species: String,
    /// Published equilibrium mole fraction.
    pub(crate) mole_fraction: f64,
}

/// One gas-species composition value from a heterogeneous external table.
///
/// This is deliberately not reused for pure condensed phases.  A source can
/// normalize gases and condensed amounts differently, so the normalization
/// must be stated by the enclosing physical reference row.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenGasSystemFraction {
    /// Gas species identity from the heterogeneous table.
    pub(crate) species: String,
    /// Species amount normalized by the complete heterogeneous system total.
    pub(crate) system_mole_fraction: f64,
}

/// One pure condensed-phase amount from a heterogeneous external table.
///
/// `phase` is a source identity such as `C(gr)`, rather than a gas species
/// alias. Keeping this separate prevents a later comparator from silently
/// applying a gas-only normalization to a condensed phase amount.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenCondensedSystemFraction {
    /// Condensed phase identity such as `C(gr)`.
    pub(crate) phase: String,
    /// Phase amount normalized by the complete heterogeneous system total.
    pub(crate) system_mole_fraction: f64,
}

/// Declared normalization of an external heterogeneous composition table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum FrozenMultiphaseNormalization {
    /// Every published gas and condensed value is divided by one total amount
    /// for the complete heterogeneous system.
    SystemTotalMoleFraction,
}

/// Frozen NASA TP-1907 multiphase fixed-`P,T` composition row.
///
/// NASA Table 11.3E contains gas species and condensed candidates in the same
/// composition table. The source rows sum to unity (within printed rounding)
/// only when graphite is included, which is retained here as an explicit
/// `SystemTotalMoleFraction` contract. This row intentionally records source
/// conditions separately from the executable element-equivalent local feed.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NasaTp1907MultiphaseReference {
    /// Problem type; must be the fixed-`P,T` code `TP`.
    pub(crate) problem_type: String,
    /// Published equilibrium temperature, K.
    pub(crate) temperature_k: f64,
    /// Fixed total pressure, Pa.
    pub(crate) pressure_pa: f64,
    /// Fuel hydrogen-to-carbon atom ratio.
    pub(crate) fuel_h_to_c_atom_ratio: f64,
    /// Fuel-to-air mass ratio.
    pub(crate) fuel_air_mass_ratio: f64,
    /// Equivalence ratio.
    pub(crate) equivalence_ratio: f64,
    /// Chemical (atom-balance) equivalence ratio.
    pub(crate) chemical_equivalence_ratio: f64,
    /// Whether the published case uses dry air (required by the contract).
    pub(crate) dry_air: bool,
    /// Declared normalization of the heterogeneous composition table.
    pub(crate) normalization: FrozenMultiphaseNormalization,
    /// Gas-species system mole fractions.
    pub(crate) gas_species: Vec<FrozenGasSystemFraction>,
    /// Condensed-phase system mole fractions.
    pub(crate) condensed_species: Vec<FrozenCondensedSystemFraction>,
}

const NASA_TP1907_MULTIPHASE_SCHEMA: [ReferenceColumnSchema; 11] = [
    ReferenceColumnSchema {
        column: "problem_type",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
    ReferenceColumnSchema {
        column: "fuel_h_to_c_atom_ratio",
        unit: "mol-atom/mol-atom",
    },
    ReferenceColumnSchema {
        column: "fuel_air_mass_ratio",
        unit: "kg/kg",
    },
    ReferenceColumnSchema {
        column: "equivalence_ratio",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "chemical_equivalence_ratio",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "dry_air",
        unit: "bool",
    },
    ReferenceColumnSchema {
        column: "normalization",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "gas_species",
        unit: "system-mole-fraction",
    },
    ReferenceColumnSchema {
        column: "condensed_species",
        unit: "system-mole-fraction",
    },
];

impl FrozenReferenceRow for NasaTp1907MultiphaseReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &NASA_TP1907_MULTIPHASE_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        const GAS_SPECIES: [&str; 9] = ["Ar", "CH4", "CO", "CO2", "H2", "H2O", "NH3", "N2", "O2"];
        const CONDENSED_SPECIES: [&str; 3] = ["C(gr)", "H2O(s)", "H2O(l)"];
        const TEMPERATURES: [f64; 4] = [680.0, 700.0, 720.0, 740.0];

        if rows.len() != TEMPERATURES.len() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: format!(
                    "expected {} Table 11.3E rows, got {}",
                    TEMPERATURES.len(),
                    rows.len()
                ),
            });
        }
        for (row_index, row) in rows.iter().enumerate() {
            if row.problem_type != "TP" {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(row_index),
                    field: "problem_type".to_owned(),
                    message: format!("expected TP, got '{}'", row.problem_type),
                });
            }
            if row.temperature_k != TEMPERATURES[row_index] {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(row_index),
                    field: "temperature_k".to_owned(),
                    message: format!(
                        "expected {} K, got {} K",
                        TEMPERATURES[row_index], row.temperature_k
                    ),
                });
            }
            for (field, value, expected) in [
                ("pressure_pa", row.pressure_pa, 101_325.0),
                ("fuel_h_to_c_atom_ratio", row.fuel_h_to_c_atom_ratio, 2.0),
                ("fuel_air_mass_ratio", row.fuel_air_mass_ratio, 0.084535),
                ("equivalence_ratio", row.equivalence_ratio, 1.25),
                (
                    "chemical_equivalence_ratio",
                    row.chemical_equivalence_ratio,
                    1.2496,
                ),
            ] {
                if !value.is_finite() || (value - expected).abs() > 1.0e-12 {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(row_index),
                        field: field.to_owned(),
                        message: format!("expected {expected}, got {value}"),
                    });
                }
            }
            if !row.dry_air
                || row.normalization != FrozenMultiphaseNormalization::SystemTotalMoleFraction
            {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(row_index),
                    field: "normalization".to_owned(),
                    message:
                        "Table 11.3E contract requires dry air and system-total mole fractions"
                            .to_owned(),
                });
            }
            validate_named_amounts(
                dataset_id,
                "gas_species",
                &row.gas_species,
                &GAS_SPECIES,
                |entry| entry.species.as_str(),
                |entry| entry.system_mole_fraction,
                false,
            )?;
            validate_named_amounts(
                dataset_id,
                "condensed_species",
                &row.condensed_species,
                &CONDENSED_SPECIES,
                |entry| entry.phase.as_str(),
                |entry| entry.system_mole_fraction,
                false,
            )?;
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
            if (total - 1.0).abs() > 2.0e-5 {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(row_index),
                    field: "system_total".to_owned(),
                    message: format!(
                        "gas plus condensed system fractions must sum to one within printed rounding, got {total}"
                    ),
                });
            }
            for phase in ["H2O(s)", "H2O(l)"] {
                let amount = row
                    .condensed_species
                    .iter()
                    .find(|entry| entry.phase == phase)
                    .expect("validated condensed list contains water phase")
                    .system_mole_fraction;
                if amount != 0.0 {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(), row: Some(row_index), field: format!("condensed_species.{phase}"),
                        message: "the selected Table 11.3E rows publish this water condensed phase as zero".to_owned(),
                    });
                }
            }
        }
        Ok(())
    }
}

/// Frozen NASA TP-1906 Table 11.3E heterogeneous equilibrium enthalpy row.
///
/// TP-1906 publishes an equilibrium-mixture *specific* enthalpy.  It remains
/// intentionally separate from the TP-1907 composition rows: each technical
/// paper owns its own source values and provenance, while the P,H adapter
/// validates their shared physical-family fields before joining them.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NasaTp1906HeterogeneousEnthalpyReference {
    /// Published equilibrium temperature, K.
    pub(crate) temperature_k: f64,
    /// Fixed total pressure, Pa.
    pub(crate) pressure_pa: f64,
    /// Fuel hydrogen-to-carbon atom ratio.
    pub(crate) fuel_h_to_c_atom_ratio: f64,
    /// Fuel-to-air mass ratio.
    pub(crate) fuel_air_mass_ratio: f64,
    /// Equivalence ratio.
    pub(crate) equivalence_ratio: f64,
    /// Chemical (atom-balance) equivalence ratio.
    pub(crate) chemical_equivalence_ratio: f64,
    /// Whether the published case uses dry air (required by the contract).
    pub(crate) dry_air: bool,
    /// Published heterogeneous-equilibrium specific enthalpy in the source
    /// unit.  This is not a molar quantity.
    pub(crate) specific_enthalpy_j_g: f64,
}

const NASA_TP1906_HETEROGENEOUS_ENTHALPY_SCHEMA: [ReferenceColumnSchema; 8] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
    ReferenceColumnSchema {
        column: "fuel_h_to_c_atom_ratio",
        unit: "mol-atom/mol-atom",
    },
    ReferenceColumnSchema {
        column: "fuel_air_mass_ratio",
        unit: "kg/kg",
    },
    ReferenceColumnSchema {
        column: "equivalence_ratio",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "chemical_equivalence_ratio",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "dry_air",
        unit: "bool",
    },
    ReferenceColumnSchema {
        column: "specific_enthalpy_j_g",
        unit: "J/g",
    },
];

impl FrozenReferenceRow for NasaTp1906HeterogeneousEnthalpyReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &NASA_TP1906_HETEROGENEOUS_ENTHALPY_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        const TEMPERATURES: [f64; 4] = [680.0, 700.0, 720.0, 740.0];
        const ENTHALPIES_J_G: [f64; 4] = [-2375.5, -2334.1, -2291.4, -2247.1];
        if rows.len() != TEMPERATURES.len() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: format!(
                    "expected {} TP-1906 Table 11.3E rows, got {}",
                    TEMPERATURES.len(),
                    rows.len()
                ),
            });
        }
        for (index, row) in rows.iter().enumerate() {
            for (field, value, expected) in [
                ("temperature_k", row.temperature_k, TEMPERATURES[index]),
                ("pressure_pa", row.pressure_pa, 101_325.0),
                ("fuel_h_to_c_atom_ratio", row.fuel_h_to_c_atom_ratio, 2.0),
                ("fuel_air_mass_ratio", row.fuel_air_mass_ratio, 0.084535),
                ("equivalence_ratio", row.equivalence_ratio, 1.25),
                (
                    "chemical_equivalence_ratio",
                    row.chemical_equivalence_ratio,
                    1.2496,
                ),
                (
                    "specific_enthalpy_j_g",
                    row.specific_enthalpy_j_g,
                    ENTHALPIES_J_G[index],
                ),
            ] {
                if !value.is_finite() || (value - expected).abs() > 1.0e-12 {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: field.to_owned(),
                        message: format!("expected {expected}, got {value}"),
                    });
                }
            }
            if !row.dry_air {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "dry_air".to_owned(),
                    message: "TP-1906 Table 11.3E fixture requires dry air".to_owned(),
                });
            }
        }
        Ok(())
    }
}

/// Frozen Argonne/STANJAN fixed-`P,T` CHON equilibrium case.
///
/// The source product table retains `C5H12 = 0`. KiThe's executable fixture
/// deliberately uses the same C/H/O/N inventory with a 15-species local gas
/// universe and therefore reports that row as an external zero reactant that
/// was not locally solved, rather than as a missing record.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ArgonneStanjanChonPtReference {
    /// Problem type; must be the fixed-`P,T` code `TP`.
    pub(crate) problem_type: String,
    /// Published equilibrium temperature, K.
    pub(crate) temperature_k: f64,
    /// Fixed total pressure, Pa.
    pub(crate) pressure_pa: f64,
    /// Feed in the source molecular basis (including `C5H12`).
    pub(crate) original_reactants: Vec<FrozenMolarAmount>,
    /// Explicitly published C/H/O/N elemental inventory totals.
    pub(crate) element_totals: Vec<FrozenElementAmount>,
    /// Published equilibrium mole fractions for the 16 declared species.
    pub(crate) species_mole_fractions: Vec<FrozenEquilibriumMoleFraction>,
}

const ARGONNE_STANJAN_CHON_PT_SCHEMA: [ReferenceColumnSchema; 6] = [
    ReferenceColumnSchema {
        column: "problem_type",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
    ReferenceColumnSchema {
        column: "original_reactants",
        unit: "mol",
    },
    ReferenceColumnSchema {
        column: "element_totals",
        unit: "mol-atoms",
    },
    ReferenceColumnSchema {
        column: "species_mole_fractions",
        unit: "1",
    },
];

impl FrozenReferenceRow for ArgonneStanjanChonPtReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &ARGONNE_STANJAN_CHON_PT_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        const REACTANTS: [&str; 4] = ["C5H12", "CH4", "O2", "N2"];
        const ELEMENTS: [&str; 4] = ["C", "H", "O", "N"];
        const SPECIES: [&str; 16] = [
            "C5H12", "CH4", "O2", "CO2", "H2O", "N2", "N", "O", "NO", "OH", "H", "N2O", "CO", "H2",
            "NO2", "HO2",
        ];
        if rows.len() != 1 {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "the Argonne/STANJAN CHON dataset must contain exactly one case"
                    .to_owned(),
            });
        }
        let row = &rows[0];
        if row.problem_type != "TP" {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: "problem_type".to_owned(),
                message: format!("expected TP, got '{}'", row.problem_type),
            });
        }
        validate_positive_finite(dataset_id, 0, "temperature_k", row.temperature_k)?;
        validate_positive_finite(dataset_id, 0, "pressure_pa", row.pressure_pa)?;
        validate_named_amounts(
            dataset_id,
            "original_reactants",
            &row.original_reactants,
            &REACTANTS,
            |amount| amount.species.as_str(),
            |amount| amount.amount_mol,
            true,
        )?;
        validate_named_amounts(
            dataset_id,
            "element_totals",
            &row.element_totals,
            &ELEMENTS,
            |amount| amount.element.as_str(),
            |amount| amount.amount_mol_atoms,
            true,
        )?;
        validate_named_amounts(
            dataset_id,
            "species_mole_fractions",
            &row.species_mole_fractions,
            &SPECIES,
            |amount| amount.species.as_str(),
            |amount| amount.mole_fraction,
            false,
        )?;
        let pentane = row
            .species_mole_fractions
            .iter()
            .find(|amount| amount.species == "C5H12")
            .expect("validated STANJAN species list must include C5H12");
        if pentane.mole_fraction != 0.0 {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: "species_mole_fractions.C5H12".to_owned(),
                message: "the reduced-universe benchmark requires an externally zero C5H12 row"
                    .to_owned(),
            });
        }
        Ok(())
    }
}

/// Validates that a named-amount list exactly matches an expected identity set,
/// rejects duplicates and undeclared identities, and checks finiteness and sign.
fn validate_named_amounts<T>(
    dataset_id: &str,
    field: &str,
    values: &[T],
    expected: &[&str],
    name: impl Fn(&T) -> &str,
    value: impl Fn(&T) -> f64,
    require_positive: bool,
) -> Result<(), FrozenReferenceError> {
    if values.len() != expected.len() {
        return Err(FrozenReferenceError::InvalidRows {
            dataset_id: dataset_id.to_owned(),
            row: Some(0),
            field: field.to_owned(),
            message: format!("expected {} rows, got {}", expected.len(), values.len()),
        });
    }
    let mut seen = HashSet::with_capacity(values.len());
    for entry in values {
        let identity = name(entry);
        let amount = value(entry);
        if !expected.contains(&identity) {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: field.to_owned(),
                message: format!("undeclared identity '{identity}'"),
            });
        }
        if !seen.insert(identity) {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: field.to_owned(),
                message: format!("duplicate identity '{identity}'"),
            });
        }
        let valid = amount.is_finite()
            && if require_positive {
                amount > 0.0
            } else {
                amount >= 0.0
            };
        if !valid {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: field.to_owned(),
                message: format!(
                    "amount for '{identity}' must be finite and {}",
                    if require_positive {
                        "positive"
                    } else {
                        "non-negative"
                    }
                ),
            });
        }
    }
    for identity in expected {
        if !seen.contains(identity) {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: Some(0),
                field: field.to_owned(),
                message: format!("missing declared identity '{identity}'"),
            });
        }
    }
    Ok(())
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct FrozenReferenceRows<R> {
    dataset_id: String,
    rows: Vec<R>,
}

/// A validated reference table together with provenance and source paths.
#[derive(Debug, Clone)]
pub(crate) struct FrozenReferenceDataset<R> {
    metadata: FrozenReferenceMetadata,
    rows: Vec<R>,
    metadata_path: PathBuf,
    data_path: PathBuf,
}

/// Structured failures from parsing or validating frozen evidence.
#[derive(Debug, Error)]
pub(crate) enum FrozenReferenceError {
    /// The source file could not be read from disk.
    #[error("failed to read frozen reference file '{path}': {message}")]
    Read { path: PathBuf, message: String },
    /// A source JSON file could not be deserialized.
    #[error("failed to parse frozen reference file '{path}': {message}")]
    Parse { path: PathBuf, message: String },
    /// Metadata violates the expected schema or content rules.
    #[error("invalid metadata for dataset '{dataset_id}', field '{field}': {message}")]
    InvalidMetadata {
        dataset_id: String,
        field: String,
        message: String,
    },
    /// The metadata and rows files declare different dataset ids.
    #[error("reference dataset id mismatch: metadata='{metadata_id}', data='{data_id}'")]
    DatasetIdMismatch {
        metadata_id: String,
        data_id: String,
    },
    /// The loaded rows filename does not match the metadata declaration.
    #[error(
        "reference data filename mismatch for dataset '{dataset_id}': metadata='{metadata_file}', loaded='{loaded_file}'"
    )]
    DataFileMismatch {
        dataset_id: String,
        metadata_file: String,
        loaded_file: String,
    },
    /// The typed schema does not match the metadata quantities.
    #[error("schema mismatch for dataset '{dataset_id}', column '{column}': {message}")]
    SchemaMismatch {
        dataset_id: String,
        column: String,
        message: String,
    },
    /// One or more numerical rows violate the physical-table invariants.
    #[error("invalid rows for dataset '{dataset_id}', row {row:?}, field '{field}': {message}")]
    InvalidRows {
        dataset_id: String,
        row: Option<usize>,
        field: String,
        message: String,
    },
}

impl<R: FrozenReferenceRow> FrozenReferenceDataset<R> {
    /// Loads and validates one immutable metadata/rows pair.
    pub(crate) fn load(
        metadata_path: impl AsRef<Path>,
        data_path: impl AsRef<Path>,
    ) -> Result<Self, FrozenReferenceError> {
        let metadata_path = metadata_path.as_ref();
        let data_path = data_path.as_ref();
        let metadata_json = read_reference_file(metadata_path)?;
        let rows_json = read_reference_file(data_path)?;
        Self::from_json(
            &metadata_json,
            &rows_json,
            metadata_path.to_path_buf(),
            data_path.to_path_buf(),
        )
    }

    /// Parses in-memory fixtures for strict failure-path tests.
    pub(crate) fn from_json_strs(
        metadata_json: &str,
        rows_json: &str,
        loaded_data_file: &str,
    ) -> Result<Self, FrozenReferenceError> {
        Self::from_json(
            metadata_json,
            rows_json,
            PathBuf::from("<inline-metadata>"),
            PathBuf::from(loaded_data_file),
        )
    }

    /// Parses and validates an in-memory metadata/rows pair with real paths.
    fn from_json(
        metadata_json: &str,
        rows_json: &str,
        metadata_path: PathBuf,
        data_path: PathBuf,
    ) -> Result<Self, FrozenReferenceError> {
        let metadata =
            serde_json::from_str(metadata_json).map_err(|error| FrozenReferenceError::Parse {
                path: metadata_path.clone(),
                message: error.to_string(),
            })?;
        let row_file: FrozenReferenceRows<R> =
            serde_json::from_str(rows_json).map_err(|error| FrozenReferenceError::Parse {
                path: data_path.clone(),
                message: error.to_string(),
            })?;

        let dataset = Self {
            metadata,
            rows: row_file.rows,
            metadata_path,
            data_path,
        };
        dataset.validate_pair_identity(&row_file.dataset_id)?;
        dataset.validate()?;
        Ok(dataset)
    }

    /// Revalidates metadata, typed schema, and numerical rows.
    pub(crate) fn validate(&self) -> Result<(), FrozenReferenceError> {
        validate_metadata::<R>(&self.metadata)?;
        R::validate_rows(&self.metadata.dataset_id, &self.rows)
    }

    /// Verifies the rows dataset id and filename agree with the metadata.
    fn validate_pair_identity(&self, rows_dataset_id: &str) -> Result<(), FrozenReferenceError> {
        if self.metadata.dataset_id != rows_dataset_id {
            return Err(FrozenReferenceError::DatasetIdMismatch {
                metadata_id: self.metadata.dataset_id.clone(),
                data_id: rows_dataset_id.to_owned(),
            });
        }

        let loaded_file = self
            .data_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default();
        if self.metadata.data_file != loaded_file {
            return Err(FrozenReferenceError::DataFileMismatch {
                dataset_id: self.metadata.dataset_id.clone(),
                metadata_file: self.metadata.data_file.clone(),
                loaded_file: loaded_file.to_owned(),
            });
        }
        Ok(())
    }

    pub(crate) fn metadata(&self) -> &FrozenReferenceMetadata {
        &self.metadata
    }

    pub(crate) fn rows(&self) -> &[R] {
        &self.rows
    }

    pub(crate) fn is_external_evidence(&self) -> bool {
        self.metadata.evidence_kind == FrozenReferenceEvidenceKind::FrozenExternal
    }

    /// Compact context intended for future comparison failures.
    pub(crate) fn provenance_context(&self, row: usize) -> String {
        let version = self
            .metadata
            .source
            .version
            .as_deref()
            .unwrap_or("unspecified");
        format!(
            "dataset={} source={}/{} version={} row={} metadata={} data={}",
            self.metadata.dataset_id,
            self.metadata.source.organization,
            self.metadata.source.name,
            version,
            row,
            self.metadata_path.display(),
            self.data_path.display()
        )
    }
}

/// Reads a frozen evidence file into a UTF-8 string, mapping I/O errors into
/// [`FrozenReferenceError::Read`].
fn read_reference_file(path: &Path) -> Result<String, FrozenReferenceError> {
    fs::read_to_string(path).map_err(|error| FrozenReferenceError::Read {
        path: path.to_path_buf(),
        message: error.to_string(),
    })
}

/// Validates metadata format version, required text fields, quantity columns,
/// and agreement with the typed row schema.
fn validate_metadata<R: FrozenReferenceRow>(
    metadata: &FrozenReferenceMetadata,
) -> Result<(), FrozenReferenceError> {
    if metadata.dataset_format_version != FROZEN_REFERENCE_FORMAT_VERSION {
        return Err(invalid_metadata(
            metadata,
            "dataset_format_version",
            format!(
                "unsupported format version {}; this build supports only {}",
                metadata.dataset_format_version, FROZEN_REFERENCE_FORMAT_VERSION
            ),
        ));
    }
    require_metadata_text(metadata, "dataset_id", &metadata.dataset_id)?;
    require_metadata_text(metadata, "title", &metadata.title)?;
    require_metadata_text(metadata, "data_file", &metadata.data_file)?;
    require_metadata_text(
        metadata,
        "source.organization",
        &metadata.source.organization,
    )?;
    require_metadata_text(metadata, "source.name", &metadata.source.name)?;
    require_metadata_text(metadata, "source.citation", &metadata.source.citation)?;
    require_metadata_text(metadata, "transcription", &metadata.transcription)?;
    validate_optional_metadata_text(
        metadata,
        "source.version",
        metadata.source.version.as_deref(),
    )?;
    validate_optional_metadata_text(
        metadata,
        "source.source_table",
        metadata.source.source_table.as_deref(),
    )?;
    validate_optional_metadata_text(
        metadata,
        "source.stable_identifier",
        metadata.source.stable_identifier.as_deref(),
    )?;

    if metadata.quantities.is_empty() {
        return Err(invalid_metadata(
            metadata,
            "quantities",
            "at least one explicitly dimensioned quantity is required",
        ));
    }

    let mut seen_columns = HashSet::with_capacity(metadata.quantities.len());
    for quantity in &metadata.quantities {
        require_metadata_text(metadata, "quantities.column", &quantity.column)?;
        require_metadata_text(metadata, "quantities.meaning", &quantity.meaning)?;
        require_metadata_text(metadata, "quantities.unit", &quantity.unit)?;
        validate_optional_metadata_text(
            metadata,
            "quantities.source_precision",
            quantity.source_precision.as_deref(),
        )?;
        if !seen_columns.insert(quantity.column.as_str()) {
            return Err(invalid_metadata(
                metadata,
                "quantities.column",
                format!("duplicate column '{}'", quantity.column),
            ));
        }
        if let Some(uncertainty) = &quantity.uncertainty {
            if !uncertainty.value.is_finite() || uncertainty.value < 0.0 {
                return Err(invalid_metadata(
                    metadata,
                    "quantities.uncertainty.value",
                    "uncertainty must be finite and non-negative",
                ));
            }
            require_metadata_text(metadata, "quantities.uncertainty.unit", &uncertainty.unit)?;
        }
    }

    let expected = R::schema();
    if metadata.quantities.len() != expected.len() {
        return Err(FrozenReferenceError::SchemaMismatch {
            dataset_id: metadata.dataset_id.clone(),
            column: "<schema>".to_owned(),
            message: format!(
                "expected {} columns, metadata declares {}",
                expected.len(),
                metadata.quantities.len()
            ),
        });
    }
    for expected_column in expected {
        let Some(actual) = metadata
            .quantities
            .iter()
            .find(|quantity| quantity.column == expected_column.column)
        else {
            return Err(FrozenReferenceError::SchemaMismatch {
                dataset_id: metadata.dataset_id.clone(),
                column: expected_column.column.to_owned(),
                message: "required typed column is missing".to_owned(),
            });
        };
        if actual.unit != expected_column.unit {
            return Err(FrozenReferenceError::SchemaMismatch {
                dataset_id: metadata.dataset_id.clone(),
                column: expected_column.column.to_owned(),
                message: format!(
                    "expected unit '{}', found '{}'",
                    expected_column.unit, actual.unit
                ),
            });
        }
    }
    Ok(())
}

/// Rejects empty required metadata text fields.
fn require_metadata_text(
    metadata: &FrozenReferenceMetadata,
    field: &str,
    value: &str,
) -> Result<(), FrozenReferenceError> {
    if value.trim().is_empty() {
        return Err(invalid_metadata(metadata, field, "must not be empty"));
    }
    Ok(())
}

/// Requires optional metadata text to be either absent or non-empty.
fn validate_optional_metadata_text(
    metadata: &FrozenReferenceMetadata,
    field: &str,
    value: Option<&str>,
) -> Result<(), FrozenReferenceError> {
    if value.is_some_and(|text| text.trim().is_empty()) {
        return Err(invalid_metadata(
            metadata,
            field,
            "must be omitted rather than stored as empty text",
        ));
    }
    Ok(())
}

/// Builds an [`FrozenReferenceError::InvalidMetadata`] referencing the dataset.
fn invalid_metadata(
    metadata: &FrozenReferenceMetadata,
    field: impl Into<String>,
    message: impl Into<String>,
) -> FrozenReferenceError {
    FrozenReferenceError::InvalidMetadata {
        dataset_id: metadata.dataset_id.clone(),
        field: field.into(),
        message: message.into(),
    }
}

/// First narrow row type used to verify the loader architecture.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct TemperaturePressureReference {
    /// Synthetic temperature point, K.
    pub(crate) temperature_k: f64,
    /// Synthetic pressure point, Pa.
    pub(crate) pressure_pa: f64,
}

const TEMPERATURE_PRESSURE_SCHEMA: [ReferenceColumnSchema; 2] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
];

impl FrozenReferenceRow for TemperaturePressureReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &TEMPERATURE_PRESSURE_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.is_empty() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "dataset must contain at least one row".to_owned(),
            });
        }

        for (index, row) in rows.iter().enumerate() {
            validate_positive_finite(dataset_id, index, "temperature_k", row.temperature_k)?;
            validate_positive_finite(dataset_id, index, "pressure_pa", row.pressure_pa)?;
            if let Some(previous) = index.checked_sub(1).map(|i| rows[i].temperature_k) {
                if row.temperature_k == previous {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "temperature_k".to_owned(),
                        message: "duplicate temperature".to_owned(),
                    });
                }
                if row.temperature_k < previous {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "temperature_k".to_owned(),
                        message: "temperatures must be strictly increasing".to_owned(),
                    });
                }
            }
        }
        Ok(())
    }
}

/// Validates that a numeric row field is finite and strictly positive.
fn validate_positive_finite(
    dataset_id: &str,
    row: usize,
    field: &str,
    value: f64,
) -> Result<(), FrozenReferenceError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(FrozenReferenceError::InvalidRows {
            dataset_id: dataset_id.to_owned(),
            row: Some(row),
            field: field.to_owned(),
            message: "must be finite and positive".to_owned(),
        });
    }
    Ok(())
}

/// I5 row for the low-pressure liquid/vapor saturation boundary of ordinary
/// water. It intentionally has its own type rather than reusing the synthetic
/// temperature/pressure record: its monotonic-pressure and triple-point
/// invariants are properties of this physical table, not of every future
/// reference dataset with the same two numeric fields.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct WaterSaturationPressureReference {
    /// Liquid-water saturation temperature, K (above the triple point).
    pub(crate) temperature_k: f64,
    /// Liquid-water saturation pressure, Pa.
    pub(crate) pressure_pa: f64,
}

const WATER_SATURATION_SCHEMA: [ReferenceColumnSchema; 2] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
];

/// The IAPWS liquid/vapor relation is only valid above the triple-point
/// transition to stable liquid water. The first frozen point is deliberately
/// 275 K, safely above this boundary.
const WATER_TRIPLE_POINT_K: f64 = 273.16;

impl FrozenReferenceRow for WaterSaturationPressureReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &WATER_SATURATION_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.is_empty() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "dataset must contain at least one row".to_owned(),
            });
        }

        for (index, row) in rows.iter().enumerate() {
            validate_positive_finite(dataset_id, index, "temperature_k", row.temperature_k)?;
            validate_positive_finite(dataset_id, index, "pressure_pa", row.pressure_pa)?;
            if row.temperature_k <= WATER_TRIPLE_POINT_K {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "temperature_k".to_owned(),
                    message: format!(
                        "liquid-water saturation rows must be above {WATER_TRIPLE_POINT_K} K"
                    ),
                });
            }
            if let Some(previous) = index.checked_sub(1).map(|i| rows[i]) {
                if row.temperature_k == previous.temperature_k {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "temperature_k".to_owned(),
                        message: "duplicate temperature".to_owned(),
                    });
                }
                if row.temperature_k < previous.temperature_k {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "temperature_k".to_owned(),
                        message: "temperatures must be strictly increasing".to_owned(),
                    });
                }
                if row.pressure_pa <= previous.pressure_pa {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "pressure_pa".to_owned(),
                        message: "liquid-water saturation pressures must be strictly increasing"
                            .to_owned(),
                    });
                }
            }
        }
        Ok(())
    }
}

/// I5 row for the low-pressure sublimation boundary of ordinary water ice Ih.
///
/// This is intentionally distinct from [`WaterSaturationPressureReference`].
/// Both tables contain temperature and pressure columns, but their validity
/// ranges, stable condensed state, and monotonicity semantics are different.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct WaterIceSublimationPressureReference {
    /// Ice-Ih sublimation temperature, K (between 50 K and the triple point).
    pub(crate) temperature_k: f64,
    /// Ice-Ih sublimation pressure, Pa.
    pub(crate) pressure_pa: f64,
}

/// One tabulated ideal-gas thermochemical row for natural-abundance CO2.
///
/// The NIST source prints all three functions dimensionlessly. Keeping the
/// frozen values in their source units avoids a hidden gas-constant conversion
/// in evidence data; consumers perform the explicit SI conversion they need.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistCo2IdealGasReference {
    /// Thermodynamic temperature, K.
    pub(crate) temperature_k: f64,
    /// Standard isobaric heat capacity divided by the source gas constant.
    pub(crate) heat_capacity_over_r: f64,
    /// Standard entropy divided by the source gas constant.
    pub(crate) entropy_over_r: f64,
    /// Standard enthalpy divided by `R*T`.
    pub(crate) enthalpy_over_rt: f64,
}

const NIST_CO2_IDEAL_GAS_SCHEMA: [ReferenceColumnSchema; 4] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "heat_capacity_over_r",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "entropy_over_r",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "enthalpy_over_rt",
        unit: "1",
    },
];

impl FrozenReferenceRow for NistCo2IdealGasReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &NIST_CO2_IDEAL_GAS_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() < 2 {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "CO2 ideal-gas evidence needs at least two temperature rows".to_owned(),
            });
        }
        let mut previous_temperature = 0.0_f64;
        for (index, row) in rows.iter().enumerate() {
            if !row.temperature_k.is_finite()
                || row.temperature_k <= previous_temperature
                || !row.heat_capacity_over_r.is_finite()
                || !row.entropy_over_r.is_finite()
                || !row.enthalpy_over_rt.is_finite()
            {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "thermochemical_row".to_owned(),
                    message: "CO2 rows must be finite and strictly temperature ordered".to_owned(),
                });
            }
            previous_temperature = row.temperature_k;
        }
        Ok(())
    }
}

const WATER_ICE_SUBLIMATION_SCHEMA: [ReferenceColumnSchema; 2] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
];

const WATER_ICE_SUBLIMATION_MIN_TEMPERATURE_K: f64 = 50.0;
const WATER_ICE_TRIPLE_POINT_K: f64 = 273.16;

impl FrozenReferenceRow for WaterIceSublimationPressureReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &WATER_ICE_SUBLIMATION_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.is_empty() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "dataset must contain at least one row".to_owned(),
            });
        }

        for (index, row) in rows.iter().enumerate() {
            validate_positive_finite(dataset_id, index, "temperature_k", row.temperature_k)?;
            validate_positive_finite(dataset_id, index, "pressure_pa", row.pressure_pa)?;
            if !(WATER_ICE_SUBLIMATION_MIN_TEMPERATURE_K..=WATER_ICE_TRIPLE_POINT_K)
                .contains(&row.temperature_k)
            {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "temperature_k".to_owned(),
                    message: format!(
                        "ice-Ih sublimation rows must lie in [{WATER_ICE_SUBLIMATION_MIN_TEMPERATURE_K}, {WATER_ICE_TRIPLE_POINT_K}] K"
                    ),
                });
            }
            if let Some(previous) = index.checked_sub(1).map(|i| rows[i]) {
                if row.temperature_k <= previous.temperature_k {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "temperature_k".to_owned(),
                        message: "temperatures must be strictly increasing".to_owned(),
                    });
                }
                if row.pressure_pa <= previous.pressure_pa {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: "pressure_pa".to_owned(),
                        message: "ice-Ih sublimation pressures must be strictly increasing"
                            .to_owned(),
                    });
                }
            }
        }
        Ok(())
    }
}

/// I5 primary NIST-JANAF species data for the Boudouard reaction
/// `2 CO(g) <=> CO2(g) + C(gr)`.
///
/// The frozen file deliberately stores the published CO and CO2 values rather
/// than a pre-combined reaction constant. Carbon is JANAF's reference graphite
/// state (`Delta_f G = 0`, `log Kf = 0`) and is recorded in the metadata.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct JanafBoudouardReference {
    /// Temperature of the JANAF datum, K.
    pub(crate) temperature_k: f64,
    /// Standard molar Gibbs energy of formation of CO, kJ/mol.
    pub(crate) co_delta_f_g_kj_mol: f64,
    /// Base-10 log of the CO formation equilibrium constant.
    pub(crate) co_log10_kf: f64,
    /// Standard molar Gibbs energy of formation of CO2, kJ/mol.
    pub(crate) co2_delta_f_g_kj_mol: f64,
    /// Base-10 log of the CO2 formation equilibrium constant.
    pub(crate) co2_log10_kf: f64,
}

const JANAF_BOUDOUARD_SCHEMA: [ReferenceColumnSchema; 5] = [
    ReferenceColumnSchema {
        column: "temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "co_delta_f_g_kj_mol",
        unit: "kJ/mol",
    },
    ReferenceColumnSchema {
        column: "co_log10_kf",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "co2_delta_f_g_kj_mol",
        unit: "kJ/mol",
    },
    ReferenceColumnSchema {
        column: "co2_log10_kf",
        unit: "1",
    },
];

impl FrozenReferenceRow for JanafBoudouardReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &JANAF_BOUDOUARD_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.is_empty() {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "dataset must contain at least one row".to_owned(),
            });
        }
        for (index, row) in rows.iter().enumerate() {
            validate_positive_finite(dataset_id, index, "temperature_k", row.temperature_k)?;
            for (field, value) in [
                ("co_delta_f_g_kj_mol", row.co_delta_f_g_kj_mol),
                ("co_log10_kf", row.co_log10_kf),
                ("co2_delta_f_g_kj_mol", row.co2_delta_f_g_kj_mol),
                ("co2_log10_kf", row.co2_log10_kf),
            ] {
                if !value.is_finite() {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: field.to_owned(),
                        message: "must be finite".to_owned(),
                    });
                }
            }
            if let Some(previous) = index.checked_sub(1).map(|i| rows[i])
                && row.temperature_k <= previous.temperature_k
            {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "temperature_k".to_owned(),
                    message: "temperatures must be strictly increasing".to_owned(),
                });
            }
        }
        Ok(())
    }
}
