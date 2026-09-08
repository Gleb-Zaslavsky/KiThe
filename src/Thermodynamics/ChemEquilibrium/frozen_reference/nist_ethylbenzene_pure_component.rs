//! Frozen NIST pure-component evidence for the ethylbenzene ternary review.
//!
//! This is deliberately an evidence and oracle module, not a replacement
//! thermochemistry handler.  It supplies reviewed numerical anchors plus
//! independent Antoine and vaporization-enthalpy correlations.  A bounded
//! `G0(T)` closure may be built only after its gas and liquid reference-state
//! conventions have been reviewed as one compatible route.

use std::path::PathBuf;

use serde::Deserialize;

use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, FrozenReferenceError, FrozenReferenceRow, ReferenceColumnSchema,
};

/// NIST liquid heat-capacity measurement retained without averaging it into a
/// synthetic correlation.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistEthylbenzeneHeatCapacityAnchor {
    /// Temperature of the liquid Cp measurement, K.
    pub(crate) temperature_k: f64,
    /// Liquid isobaric heat capacity, J/(mol·K).
    pub(crate) heat_capacity_j_mol_k: f64,
}

/// One NIST vaporization-enthalpy characterization value.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistEthylbenzeneVaporizationPoint {
    /// Temperature of the vaporization measurement, K.
    pub(crate) temperature_k: f64,
    /// Molar enthalpy of vaporization, J/mol.
    pub(crate) vaporization_enthalpy_j_mol: f64,
    /// Optional measurement uncertainty, J/mol.
    #[serde(default)]
    pub(crate) uncertainty_j_mol: Option<f64>,
}

/// Majer-Svoboda vaporization-enthalpy correlation frozen as an independent
/// pure-component characterization route.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistEthylbenzeneMajerSvobodaCorrelation {
    /// Leading coefficient of the vaporization-enthalpy correlation, J/mol.
    pub(crate) coefficient_j_mol: f64,
    /// Dimensionless exponent `beta`.
    pub(crate) beta: f64,
    /// Critical temperature used in the reduced-temperature term, K.
    pub(crate) critical_temperature_k: f64,
    /// Lower bound of the reviewed validity interval, K.
    pub(crate) lower_temperature_k: f64,
    /// Upper bound of the reviewed validity interval, K.
    pub(crate) upper_temperature_k: f64,
}

/// NIST Antoine correlation for `log10(P/bar)`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistEthylbenzeneAntoineCorrelation {
    /// Antoine constant `A`.
    pub(crate) a: f64,
    /// Antoine constant `B`, K.
    pub(crate) b_k: f64,
    /// Antoine constant `C`, K.
    pub(crate) c_k: f64,
    /// Lower bound of the reviewed validity interval, K.
    pub(crate) lower_temperature_k: f64,
    /// Upper bound of the reviewed validity interval, K.
    pub(crate) upper_temperature_k: f64,
}

/// One complete approved numerical-seed row for ethylbenzene.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistEthylbenzenePureComponentReference {
    /// Compound name (must be `ethylbenzene`).
    pub(crate) compound: String,
    /// Molecular formula (must be `C8H10`).
    pub(crate) formula: String,
    /// CAS registry number (must be `100-41-4`).
    pub(crate) cas_registry_number: String,
    /// Standard reference temperature of the formation anchors, K (298.15).
    pub(crate) reference_temperature_k: f64,
    /// Liquid standard formation enthalpy, J/mol.
    pub(crate) liquid_formation_enthalpy_j_mol: f64,
    /// Uncertainty of the liquid formation enthalpy, J/mol.
    pub(crate) liquid_formation_enthalpy_uncertainty_j_mol: f64,
    /// Liquid standard molar entropy, J/(mol·K).
    pub(crate) liquid_standard_entropy_j_mol_k: f64,
    /// Gas standard formation enthalpy, J/mol.
    pub(crate) gas_formation_enthalpy_j_mol: f64,
    /// Uncertainty of the gas formation enthalpy, J/mol.
    pub(crate) gas_formation_enthalpy_uncertainty_j_mol: f64,
    /// Approved liquid heat-capacity anchors (six values).
    pub(crate) liquid_heat_capacity_anchors: Vec<NistEthylbenzeneHeatCapacityAnchor>,
    /// Approved vaporization-enthalpy points (eleven values).
    pub(crate) vaporization_enthalpy_points: Vec<NistEthylbenzeneVaporizationPoint>,
    /// Frozen Majer-Svoboda vaporization-enthalpy correlation.
    pub(crate) majer_svoboda: NistEthylbenzeneMajerSvobodaCorrelation,
    /// Frozen Antoine saturation-pressure correlation.
    pub(crate) antoine: NistEthylbenzeneAntoineCorrelation,
}

const ETHYLBENZENE_SEED_SCHEMA: [ReferenceColumnSchema; 13] = [
    ReferenceColumnSchema {
        column: "compound",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "formula",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "cas_registry_number",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "reference_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "liquid_formation_enthalpy_j_mol",
        unit: "J/mol",
    },
    ReferenceColumnSchema {
        column: "liquid_formation_enthalpy_uncertainty_j_mol",
        unit: "J/mol",
    },
    ReferenceColumnSchema {
        column: "liquid_standard_entropy_j_mol_k",
        unit: "J/mol/K",
    },
    ReferenceColumnSchema {
        column: "gas_formation_enthalpy_j_mol",
        unit: "J/mol",
    },
    ReferenceColumnSchema {
        column: "gas_formation_enthalpy_uncertainty_j_mol",
        unit: "J/mol",
    },
    ReferenceColumnSchema {
        column: "liquid_heat_capacity_anchors",
        unit: "J/mol/K",
    },
    ReferenceColumnSchema {
        column: "vaporization_enthalpy_points",
        unit: "J/mol",
    },
    ReferenceColumnSchema {
        column: "majer_svoboda",
        unit: "J/mol",
    },
    ReferenceColumnSchema {
        column: "antoine",
        unit: "bar",
    },
];

impl FrozenReferenceRow for NistEthylbenzenePureComponentReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &ETHYLBENZENE_SEED_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() != 1 {
            return Err(invalid_row(
                dataset_id,
                None,
                "rows",
                "ethylbenzene numerical seed evidence requires exactly one typed row",
            ));
        }
        let row = &rows[0];
        if row.compound != "ethylbenzene"
            || row.formula != "C8H10"
            || row.cas_registry_number != "100-41-4"
        {
            return Err(invalid_row(
                dataset_id,
                Some(0),
                "identity",
                "expected ethylbenzene, formula C8H10, CAS 100-41-4",
            ));
        }
        for (field, value) in [
            ("reference_temperature_k", row.reference_temperature_k),
            (
                "liquid_formation_enthalpy_j_mol",
                row.liquid_formation_enthalpy_j_mol,
            ),
            (
                "liquid_formation_enthalpy_uncertainty_j_mol",
                row.liquid_formation_enthalpy_uncertainty_j_mol,
            ),
            (
                "liquid_standard_entropy_j_mol_k",
                row.liquid_standard_entropy_j_mol_k,
            ),
            (
                "gas_formation_enthalpy_j_mol",
                row.gas_formation_enthalpy_j_mol,
            ),
            (
                "gas_formation_enthalpy_uncertainty_j_mol",
                row.gas_formation_enthalpy_uncertainty_j_mol,
            ),
        ] {
            if !value.is_finite() {
                return Err(invalid_row(dataset_id, Some(0), field, "must be finite"));
            }
        }
        if row.reference_temperature_k != 298.15
            || row.liquid_formation_enthalpy_uncertainty_j_mol < 0.0
            || row.gas_formation_enthalpy_uncertainty_j_mol < 0.0
        {
            return Err(invalid_row(
                dataset_id,
                Some(0),
                "standard_state_anchors",
                "reference temperature must be 298.15 K and uncertainties must be non-negative",
            ));
        }
        if row.liquid_heat_capacity_anchors.len() != 6 {
            return Err(invalid_row(
                dataset_id,
                Some(0),
                "liquid_heat_capacity_anchors",
                "expected the six approved liquid Cp anchors",
            ));
        }
        for anchor in &row.liquid_heat_capacity_anchors {
            if !anchor.temperature_k.is_finite()
                || !anchor.heat_capacity_j_mol_k.is_finite()
                || anchor.temperature_k <= 0.0
                || anchor.heat_capacity_j_mol_k <= 0.0
            {
                return Err(invalid_row(
                    dataset_id,
                    Some(0),
                    "liquid_heat_capacity_anchors",
                    "temperature and heat capacity must be finite and positive",
                ));
            }
        }
        if row.vaporization_enthalpy_points.len() != 11 {
            return Err(invalid_row(
                dataset_id,
                Some(0),
                "vaporization_enthalpy_points",
                "expected the eleven approved vaporization-enthalpy points",
            ));
        }
        let mut previous_temperature = 0.0_f64;
        for point in &row.vaporization_enthalpy_points {
            if !point.temperature_k.is_finite()
                || !point.vaporization_enthalpy_j_mol.is_finite()
                || point.temperature_k <= previous_temperature
                || point.vaporization_enthalpy_j_mol <= 0.0
                || point
                    .uncertainty_j_mol
                    .is_some_and(|value| !value.is_finite() || value < 0.0)
            {
                return Err(invalid_row(
                    dataset_id,
                    Some(0),
                    "vaporization_enthalpy_points",
                    "points must be finite, positive, and strictly ordered by temperature",
                ));
            }
            previous_temperature = point.temperature_k;
        }
        let correlation = row.majer_svoboda;
        if !correlation.coefficient_j_mol.is_finite()
            || !correlation.beta.is_finite()
            || !correlation.critical_temperature_k.is_finite()
            || !correlation.lower_temperature_k.is_finite()
            || !correlation.upper_temperature_k.is_finite()
            || correlation.coefficient_j_mol <= 0.0
            || correlation.beta <= 0.0
            || correlation.critical_temperature_k <= correlation.upper_temperature_k
            || correlation.lower_temperature_k >= correlation.upper_temperature_k
        {
            return Err(invalid_row(
                dataset_id,
                Some(0),
                "majer_svoboda",
                "correlation must define a positive finite subcritical temperature interval",
            ));
        }
        let antoine = row.antoine;
        if !antoine.a.is_finite()
            || !antoine.b_k.is_finite()
            || !antoine.c_k.is_finite()
            || !antoine.lower_temperature_k.is_finite()
            || !antoine.upper_temperature_k.is_finite()
            || antoine.b_k <= 0.0
            || antoine.lower_temperature_k >= antoine.upper_temperature_k
        {
            return Err(invalid_row(
                dataset_id,
                Some(0),
                "antoine",
                "Antoine correlation must define finite coefficients and an ordered interval",
            ));
        }
        Ok(())
    }
}

/// Builds an [`FrozenReferenceError::InvalidRows`] for this seed schema.
fn invalid_row(
    dataset_id: &str,
    row: Option<usize>,
    field: impl Into<String>,
    message: impl Into<String>,
) -> FrozenReferenceError {
    FrozenReferenceError::InvalidRows {
        dataset_id: dataset_id.to_owned(),
        row,
        field: field.into(),
        message: message.into(),
    }
}

/// Returns the NIST WebBook frozen-data directory path.
fn seed_directory() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src")
        .join("Thermodynamics")
        .join("ChemEquilibrium")
        .join("frozen_reference")
        .join("data")
        .join("nist_webbook")
}

/// Loads the approved NIST numerical seeds without network access.
pub(crate) fn load_nist_ethylbenzene_pure_component_seed()
-> Result<FrozenReferenceDataset<NistEthylbenzenePureComponentReference>, FrozenReferenceError> {
    let directory = seed_directory();
    FrozenReferenceDataset::load(
        directory.join("ethylbenzene_pure_component_seed.metadata.json"),
        directory.join("ethylbenzene_pure_component_seed.rows.json"),
    )
}

/// Evaluates the frozen NIST Antoine correlation and returns saturation pressure
/// in Pa.  This is an independent pure-component oracle, not a VLE fit.
pub(crate) fn nist_ethylbenzene_antoine_pressure_pa(
    correlation: NistEthylbenzeneAntoineCorrelation,
    temperature_k: f64,
) -> Result<f64, String> {
    if !temperature_k.is_finite()
        || temperature_k < correlation.lower_temperature_k
        || temperature_k > correlation.upper_temperature_k
    {
        return Err(format!(
            "ethylbenzene Antoine correlation is valid on [{:.2}, {:.2}] K, got {temperature_k}",
            correlation.lower_temperature_k, correlation.upper_temperature_k
        ));
    }
    let denominator = temperature_k + correlation.c_k;
    if denominator <= 0.0 || !denominator.is_finite() {
        return Err(format!(
            "ethylbenzene Antoine denominator must be finite and positive, got {denominator}"
        ));
    }
    let pressure_pa = 100_000.0 * 10_f64.powf(correlation.a - correlation.b_k / denominator);
    if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(format!(
            "ethylbenzene Antoine correlation produced invalid pressure {pressure_pa} Pa"
        ));
    }
    Ok(pressure_pa)
}

/// Evaluates the frozen Majer-Svoboda vaporization-enthalpy correlation in
/// J/mol inside its reviewed validity interval.
pub(crate) fn nist_ethylbenzene_majer_svoboda_vaporization_enthalpy_j_mol(
    correlation: NistEthylbenzeneMajerSvobodaCorrelation,
    temperature_k: f64,
) -> Result<f64, String> {
    if !temperature_k.is_finite()
        || temperature_k < correlation.lower_temperature_k
        || temperature_k > correlation.upper_temperature_k
    {
        return Err(format!(
            "ethylbenzene Majer-Svoboda correlation is valid on [{:.2}, {:.2}] K, got {temperature_k}",
            correlation.lower_temperature_k, correlation.upper_temperature_k
        ));
    }
    let reduced_temperature = temperature_k / correlation.critical_temperature_k;
    let residual = 1.0 - reduced_temperature;
    let enthalpy = correlation.coefficient_j_mol
        * (-correlation.beta * reduced_temperature).exp()
        * residual.powf(correlation.beta);
    if !enthalpy.is_finite() || enthalpy <= 0.0 {
        return Err(format!(
            "ethylbenzene Majer-Svoboda correlation produced invalid vaporization enthalpy {enthalpy} J/mol"
        ));
    }
    Ok(enthalpy)
}
