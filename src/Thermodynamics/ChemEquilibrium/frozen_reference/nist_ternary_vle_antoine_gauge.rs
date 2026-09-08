//! Test-only `P,T` thermodynamic gauge for the ternary NIST VLE candidate.
//!
//! The selected universe has three independently conserved molecular
//! identities and no interconversion reactions.  For that deliberately narrow
//! case, pure-component saturation pressure fixes each gas/liquid standard
//! Gibbs *difference*.  We choose `G0_gas = 0` and derive
//! `G0_liquid = R T ln(Psat / p0)` at an explicit `p0 = 1 bar`.
//!
//! This is not a production thermochemistry library and is explicitly
//! unsuitable for `P,H`: it owns no absolute enthalpy or heat-capacity model.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::rc::Rc;

use nalgebra::DMatrix;
use serde::Deserialize;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{GibbsFn, R};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, SolveError, compute_reaction_basis,
};
use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, FrozenReferenceError, FrozenReferenceRow, ReferenceColumnSchema,
};

/// The one explicit standard pressure used by every frozen Antoine gauge row.
pub(crate) const FROZEN_ANTOINE_REFERENCE_PRESSURE_PA: f64 = 100_000.0;

/// State for one of the two gauge closures of a molecular identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FrozenAntoineGaugeState {
    /// The ideal-gas standard Gibbs gauge is zero.
    Gas,
    /// The ideal-liquid standard Gibbs value comes from pure saturation.
    Liquid,
}

/// One exact pure-component Antoine record, independent of ternary VLE rows.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FrozenAntoineStandardState {
    /// Compound name (must be unique across the three identities).
    pub(crate) compound: String,
    /// Molecular formula of the compound.
    pub(crate) formula: String,
    /// CAS registry number of the compound.
    pub(crate) cas_registry_number: String,
    /// Element-to-atom-count composition map.
    pub(crate) elemental_composition: BTreeMap<String, f64>,
    /// Antoine constant `A`.
    pub(crate) a: f64,
    /// Antoine constant `B`, K.
    pub(crate) b_k: f64,
    /// Antoine constant `C`, K.
    pub(crate) c_k: f64,
    /// Lower bound of the valid temperature interval, K.
    pub(crate) lower_temperature_k: f64,
    /// Upper bound of the valid temperature interval, K.
    pub(crate) upper_temperature_k: f64,
    /// Pressure unit of the Antoine relation (must be `bar`).
    pub(crate) pressure_unit: String,
    /// Source route from which the coefficients were transcribed.
    pub(crate) source_route: String,
}

/// One published ternary liquid-composition boiling-temperature row.
///
/// The source has no vapour-composition measurement. `source_temperature_k`
/// is characterization evidence only; the independent Raoult route derives
/// its own bubble temperature and vapour composition.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NistTernaryVleInteriorReference {
    /// Liquid-phase toluene mole fraction.
    pub(crate) liquid_toluene_mole_fraction: f64,
    /// Liquid-phase ethylbenzene mole fraction.
    pub(crate) liquid_ethylbenzene_mole_fraction: f64,
    /// Liquid-phase chlorobenzene mole fraction.
    pub(crate) liquid_chlorobenzene_mole_fraction: f64,
    /// Total pressure of the published row, Pa.
    pub(crate) pressure_pa: f64,
    /// Published boiling-temperature characterization, K.
    pub(crate) source_temperature_k: f64,
    /// Uncertainty of the published temperature, K.
    pub(crate) source_temperature_uncertainty_k: f64,
}

/// Result of the independent ideal-Raoult bubble calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TernaryRaoultBubbleSolution {
    /// Bubble temperature found by the independent route, K.
    pub(crate) temperature_k: f64,
    /// Total bubble pressure, Pa.
    pub(crate) pressure_pa: f64,
    /// Derived vapor-phase mole fractions.
    pub(crate) vapor_mole_fractions: [f64; 3],
    /// Number of iterations taken by the bisection solver.
    pub(crate) iterations: usize,
}

/// Result of the independent ideal-Raoult dew calculation.
///
/// At the dew boundary the supplied vapour composition is the bulk gas
/// composition and the returned liquid composition is incipient. This is kept
/// separate from the bubble result because interchanging `x` and `y` is a
/// common and physically consequential VLE test error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TernaryRaoultDewSolution {
    /// Dew temperature found by the independent route, K.
    pub(crate) temperature_k: f64,
    /// Total dew pressure, Pa.
    pub(crate) pressure_pa: f64,
    /// Incipient liquid-phase mole fractions.
    pub(crate) liquid_mole_fractions: [f64; 3],
    /// Number of iterations taken by the bisection solver.
    pub(crate) iterations: usize,
}

/// Physical phase classification returned by the independent flash oracle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TernaryRaoultFlashPhaseClass {
    /// The Rachford-Rice root lies strictly inside `0 < beta < 1`.
    TwoPhase,
    /// The Rachford-Rice liquid endpoint proves a single liquid phase.
    AllLiquid,
    /// The Rachford-Rice vapor endpoint proves a single gas phase.
    AllVapor,
}

/// Independent ideal-Raoult flash result for one ternary bulk inventory.
///
/// `liquid_mole_fractions` and `vapor_mole_fractions` are present only for a
/// genuine two-phase split. A single-phase result retains the corresponding
/// bulk composition in its physical phase and never invents an incipient
/// second-phase composition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TernaryRaoultFlashSolution {
    /// Physical phase classification of the flash result.
    pub(crate) phase_class: TernaryRaoultFlashPhaseClass,
    /// Flash temperature, K.
    pub(crate) temperature_k: f64,
    /// Flash pressure, Pa.
    pub(crate) pressure_pa: f64,
    /// Vapor fraction `beta`, present only for a two-phase split.
    pub(crate) vapor_fraction: Option<f64>,
    /// Liquid mole fractions, present only for a two-phase split.
    pub(crate) liquid_mole_fractions: Option<[f64; 3]>,
    /// Vapor mole fractions, present only for a two-phase split.
    pub(crate) vapor_mole_fractions: Option<[f64; 3]>,
    /// Number of iterations taken by the flash solver.
    pub(crate) iterations: usize,
}

const TERNARY_ANTOINE_SCHEMA: [ReferenceColumnSchema; 11] = [
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
        column: "elemental_composition",
        unit: "atoms",
    },
    ReferenceColumnSchema {
        column: "a",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "b_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "c_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "lower_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "upper_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "pressure_unit",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "source_route",
        unit: "text",
    },
];

const TERNARY_VLE_INTERIOR_SCHEMA: [ReferenceColumnSchema; 6] = [
    ReferenceColumnSchema {
        column: "liquid_toluene_mole_fraction",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "liquid_ethylbenzene_mole_fraction",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "liquid_chlorobenzene_mole_fraction",
        unit: "1",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
    ReferenceColumnSchema {
        column: "source_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "source_temperature_uncertainty_k",
        unit: "K",
    },
];

impl FrozenReferenceRow for FrozenAntoineStandardState {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &TERNARY_ANTOINE_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() != 3 {
            return Err(invalid_row(
                dataset_id,
                None,
                "rows",
                "ternary Antoine gauge requires exactly three molecular identities",
            ));
        }
        let mut identities = BTreeSet::new();
        for (index, row) in rows.iter().enumerate() {
            if row.compound.trim().is_empty()
                || row.formula.trim().is_empty()
                || row.cas_registry_number.trim().is_empty()
                || row.source_route.trim().is_empty()
                || row.pressure_unit != "bar"
            {
                return Err(invalid_row(
                    dataset_id,
                    Some(index),
                    "identity",
                    "identity/source text must be non-empty and Antoine pressure unit must be bar",
                ));
            }
            if !identities.insert(row.compound.as_str()) {
                return Err(invalid_row(
                    dataset_id,
                    Some(index),
                    "compound",
                    "each molecular identity must occur exactly once",
                ));
            }
            for (field, value) in [
                ("a", row.a),
                ("b_k", row.b_k),
                ("c_k", row.c_k),
                ("lower_temperature_k", row.lower_temperature_k),
                ("upper_temperature_k", row.upper_temperature_k),
            ] {
                if !value.is_finite() {
                    return Err(invalid_row(
                        dataset_id,
                        Some(index),
                        field,
                        "must be finite",
                    ));
                }
            }
            if row.b_k <= 0.0 || row.lower_temperature_k >= row.upper_temperature_k {
                return Err(invalid_row(
                    dataset_id,
                    Some(index),
                    "temperature_interval",
                    "Antoine coefficient B must be positive and the interval must be ordered",
                ));
            }
            if row.elemental_composition.is_empty()
                || row.elemental_composition.iter().any(|(element, amount)| {
                    element.trim().is_empty() || !amount.is_finite() || *amount <= 0.0
                })
            {
                return Err(invalid_row(
                    dataset_id,
                    Some(index),
                    "elemental_composition",
                    "element labels and atom counts must be finite and positive",
                ));
            }
        }
        Ok(())
    }
}

impl FrozenReferenceRow for NistTernaryVleInteriorReference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &TERNARY_VLE_INTERIOR_SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() < 3 {
            return Err(invalid_row(
                dataset_id,
                None,
                "rows",
                "ternary VLE characterization needs at least three genuine interior rows",
            ));
        }
        for (index, row) in rows.iter().enumerate() {
            let composition = [
                row.liquid_toluene_mole_fraction,
                row.liquid_ethylbenzene_mole_fraction,
                row.liquid_chlorobenzene_mole_fraction,
            ];
            if composition
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0 || *value >= 1.0)
                || (composition.iter().sum::<f64>() - 1.0).abs() > 1.0e-12
            {
                return Err(invalid_row(
                    dataset_id,
                    Some(index),
                    "liquid_composition",
                    "all selected rows must be strict ternary-simplex points",
                ));
            }
            for (field, value) in [
                ("pressure_pa", row.pressure_pa),
                ("source_temperature_k", row.source_temperature_k),
                (
                    "source_temperature_uncertainty_k",
                    row.source_temperature_uncertainty_k,
                ),
            ] {
                if !value.is_finite() || value <= 0.0 {
                    return Err(invalid_row(
                        dataset_id,
                        Some(index),
                        field,
                        "must be finite and positive",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Builds an [`FrozenReferenceError::InvalidRows`] for the ternary evidence.
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
fn gauge_directory() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src")
        .join("Thermodynamics")
        .join("ChemEquilibrium")
        .join("frozen_reference")
        .join("data")
        .join("nist_webbook")
}

/// Loads the three frozen pure-component Antoine records without a network or
/// a production-library lookup.
pub(crate) fn load_nist_ternary_vle_antoine_gauge()
-> Result<FrozenReferenceDataset<FrozenAntoineStandardState>, FrozenReferenceError> {
    let directory = gauge_directory();
    FrozenReferenceDataset::load(
        directory.join("toluene_ethylbenzene_chlorobenzene_antoine.metadata.json"),
        directory.join("toluene_ethylbenzene_chlorobenzene_antoine.rows.json"),
    )
}

/// Loads selected published ThermoML rows. No derived vapour composition is
/// stored in this evidence file.
pub(crate) fn load_nist_ternary_vle_interior_rows()
-> Result<FrozenReferenceDataset<NistTernaryVleInteriorReference>, FrozenReferenceError> {
    let directory = gauge_directory().join("..").join("nist_thermoml");
    FrozenReferenceDataset::load(
        directory.join("toluene_ethylbenzene_chlorobenzene_selected.metadata.json"),
        directory.join("toluene_ethylbenzene_chlorobenzene_selected.rows.json"),
    )
}

/// Builds the six standard-Gibbs closures used by the frozen ternary gauge.
///
/// The component order is intentionally fixed: all three gas components come
/// first, followed by the matching liquid components.  The production runner
/// sees ordinary Gibbs closures; it has no ternary-specific equilibrium code.
/// The caller is still responsible for evaluating only the common validated
/// Antoine interval.
pub(crate) fn ternary_gauge_gibbs_functions(
    records: &[FrozenAntoineStandardState],
) -> Result<Vec<GibbsFn>, ReactionExtentError> {
    if records.len() != 3 {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "ternary Antoine gauge requires exactly three records, got {}",
            records.len()
        )));
    }

    let mut gibbs = (0..3)
        .map(|_| Rc::new(|_| 0.0) as GibbsFn)
        .collect::<Vec<_>>();
    gibbs.extend(records.iter().cloned().map(|record| {
        Rc::new(move |temperature| {
            gauge_standard_gibbs_j_mol(&record, FrozenAntoineGaugeState::Liquid, temperature)
                .expect("ternary gauge closures are evaluated only inside their validated interval")
        }) as GibbsFn
    }));
    Ok(gibbs)
}

/// Intersects all record domains and rejects an empty overlap.
pub(crate) fn common_temperature_interval(
    records: &[FrozenAntoineStandardState],
) -> Result<(f64, f64), ReactionExtentError> {
    let lower = records
        .iter()
        .map(|record| record.lower_temperature_k)
        .reduce(f64::max)
        .ok_or_else(|| ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_records",
            message: "at least one Antoine record is required".to_owned(),
        })?;
    let upper = records
        .iter()
        .map(|record| record.upper_temperature_k)
        .reduce(f64::min)
        .expect("non-empty records already validated");
    if lower > upper {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_common_interval",
            message: format!(
                "frozen Antoine records have no common interval: [{lower}, {upper}] K"
            ),
        });
    }
    Ok((lower, upper))
}

/// Evaluates one frozen `log10(Psat/bar)` Antoine relation in Pa.
pub(crate) fn saturation_pressure_pa(
    record: &FrozenAntoineStandardState,
    temperature_k: f64,
) -> Result<f64, ReactionExtentError> {
    if !temperature_k.is_finite()
        || temperature_k < record.lower_temperature_k
        || temperature_k > record.upper_temperature_k
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_temperature",
            message: format!(
                "{} Antoine record is valid on [{:.2}, {:.2}] K, got {temperature_k}",
                record.compound, record.lower_temperature_k, record.upper_temperature_k
            ),
        });
    }
    let denominator = temperature_k + record.c_k;
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_denominator",
            message: format!(
                "{} Antoine denominator is invalid at {temperature_k} K",
                record.compound
            ),
        });
    }
    let pressure_pa =
        FROZEN_ANTOINE_REFERENCE_PRESSURE_PA * 10_f64.powf(record.a - record.b_k / denominator);
    if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_pressure",
            message: format!(
                "{} Antoine relation returned {pressure_pa} Pa",
                record.compound
            ),
        });
    }
    Ok(pressure_pa)
}

/// Evaluates the bounded test-only standard Gibbs gauge in J/mol.
pub(crate) fn gauge_standard_gibbs_j_mol(
    record: &FrozenAntoineStandardState,
    state: FrozenAntoineGaugeState,
    temperature_k: f64,
) -> Result<f64, ReactionExtentError> {
    match state {
        FrozenAntoineGaugeState::Gas => {
            // Call the pressure relation for its domain check even though the
            // selected gas gauge itself is exactly zero.
            saturation_pressure_pa(record, temperature_k)?;
            Ok(0.0)
        }
        FrozenAntoineGaugeState::Liquid => {
            let saturation_pressure = saturation_pressure_pa(record, temperature_k)?;
            Ok(R * temperature_k
                * (saturation_pressure / FROZEN_ANTOINE_REFERENCE_PRESSURE_PA).ln())
        }
    }
}

/// Recovers pure saturation pressure from the gauge chemical-potential
/// difference.  This is a construction round-trip, not another data source.
pub(crate) fn saturation_pressure_from_gauge_pa(
    gas_standard_gibbs_j_mol: f64,
    liquid_standard_gibbs_j_mol: f64,
    temperature_k: f64,
) -> Result<f64, ReactionExtentError> {
    if !gas_standard_gibbs_j_mol.is_finite()
        || !liquid_standard_gibbs_j_mol.is_finite()
        || !temperature_k.is_finite()
        || temperature_k <= 0.0
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_gauge",
            message: "Gibbs values must be finite and temperature must be positive".to_owned(),
        });
    }
    let pressure = FROZEN_ANTOINE_REFERENCE_PRESSURE_PA
        * ((liquid_standard_gibbs_j_mol - gas_standard_gibbs_j_mol) / (R * temperature_k)).exp();
    if !pressure.is_finite() || pressure <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_gauge",
            message: format!("gauge round-trip returned invalid pressure {pressure} Pa"),
        });
    }
    Ok(pressure)
}

/// Builds the deterministic molecule-by-element matrix in source row order.
pub(crate) fn molecular_element_matrix(
    records: &[FrozenAntoineStandardState],
) -> Result<(Vec<String>, DMatrix<f64>), ReactionExtentError> {
    if records.is_empty() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "frozen_antoine_records",
            message: "at least one Antoine record is required".to_owned(),
        });
    }
    let labels = records
        .iter()
        .flat_map(|record| record.elemental_composition.keys().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let mut matrix = DMatrix::zeros(records.len(), labels.len());
    for (row, record) in records.iter().enumerate() {
        for (column, label) in labels.iter().enumerate() {
            matrix[(row, column)] = record
                .elemental_composition
                .get(label)
                .copied()
                .unwrap_or(0.0);
        }
    }
    Ok((labels, matrix))
}

/// Verifies that duplicating each molecular identity into gas and liquid gives
/// exactly one transfer reaction per identity and no chemical conversion.
pub(crate) fn phase_transfer_reaction_basis(
    records: &[FrozenAntoineStandardState],
) -> Result<(usize, usize), ReactionExtentError> {
    let (_, molecular) = molecular_element_matrix(records)?;
    let mut duplicated = DMatrix::zeros(molecular.nrows() * 2, molecular.ncols());
    for row in 0..molecular.nrows() {
        duplicated.row_mut(row).copy_from(&molecular.row(row));
        duplicated
            .row_mut(row + molecular.nrows())
            .copy_from(&molecular.row(row));
    }
    let basis = compute_reaction_basis(&duplicated, 1.0e-10)?;
    Ok((basis.rank, basis.num_reactions))
}

/// Computes the ideal-Raoult bubble pressure from a liquid simplex point.
pub(crate) fn raoult_bubble_pressure_pa(
    records: &[FrozenAntoineStandardState],
    temperature_k: f64,
    liquid_mole_fractions: [f64; 3],
) -> Result<f64, ReactionExtentError> {
    if records.len() != 3 {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "ternary Raoult route requires three Antoine records, got {}",
            records.len()
        )));
    }
    if liquid_mole_fractions
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        || (liquid_mole_fractions.iter().sum::<f64>() - 1.0).abs() > 1.0e-12
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_liquid_composition",
            message: "liquid mole fractions must be finite, non-negative, and sum to one"
                .to_owned(),
        });
    }
    records
        .iter()
        .zip(liquid_mole_fractions)
        .try_fold(0.0, |pressure, (record, mole_fraction)| {
            Ok(pressure + mole_fraction * saturation_pressure_pa(record, temperature_k)?)
        })
}

/// Solves the independent ideal-Raoult bubble-temperature equation with
/// bounded bisection over the Antoine common interval.
pub(crate) fn solve_raoult_bubble_temperature(
    records: &[FrozenAntoineStandardState],
    liquid_mole_fractions: [f64; 3],
    pressure_pa: f64,
) -> Result<TernaryRaoultBubbleSolution, ReactionExtentError> {
    if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_pressure",
            message: "bubble pressure must be finite and positive".to_owned(),
        });
    }
    let (mut lower, mut upper) = common_temperature_interval(records)?;
    let mut lower_value =
        raoult_bubble_pressure_pa(records, lower, liquid_mole_fractions)? - pressure_pa;
    let upper_value =
        raoult_bubble_pressure_pa(records, upper, liquid_mole_fractions)? - pressure_pa;
    if lower_value == 0.0 {
        return raoult_solution(records, liquid_mole_fractions, pressure_pa, lower, 0);
    }
    if upper_value == 0.0 {
        return raoult_solution(records, liquid_mole_fractions, pressure_pa, upper, 0);
    }
    if lower_value.signum() == upper_value.signum() {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "ternary_raoult_bubble_temperature",
            message: format!(
                "pressure {pressure_pa} Pa is not bracketed by the frozen common interval [{lower}, {upper}] K"
            ),
        });
    }
    for iteration in 1..=128 {
        let temperature = 0.5 * (lower + upper);
        let value =
            raoult_bubble_pressure_pa(records, temperature, liquid_mole_fractions)? - pressure_pa;
        if value.abs() <= pressure_pa * 1.0e-12 || (upper - lower) <= 1.0e-10 {
            return raoult_solution(
                records,
                liquid_mole_fractions,
                pressure_pa,
                temperature,
                iteration,
            );
        }
        if value.signum() == lower_value.signum() {
            lower = temperature;
            lower_value = value;
        } else {
            upper = temperature;
        }
    }
    Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
}

/// Solves the independent ideal-Raoult dew-temperature equation with bounded
/// bisection over the common Antoine interval.
///
/// The root is `sum_i y_i * P / Psat_i(T) = 1`. No Gibbs closure, TPD result,
/// production solver state, or experimental ternary row participates in this
/// calculation.
pub(crate) fn solve_raoult_dew_temperature(
    records: &[FrozenAntoineStandardState],
    vapor_mole_fractions: [f64; 3],
    pressure_pa: f64,
) -> Result<TernaryRaoultDewSolution, ReactionExtentError> {
    if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_dew_pressure",
            message: "dew pressure must be finite and positive".to_owned(),
        });
    }
    validate_ternary_simplex("ternary_raoult_dew_vapor_composition", vapor_mole_fractions)?;
    let (mut lower, mut upper) = common_temperature_interval(records)?;
    let mut lower_value = raoult_dew_residual(records, lower, vapor_mole_fractions, pressure_pa)?;
    let upper_value = raoult_dew_residual(records, upper, vapor_mole_fractions, pressure_pa)?;
    if lower_value.abs() <= 1.0e-12 {
        return raoult_dew_solution(records, vapor_mole_fractions, pressure_pa, lower, 0);
    }
    if upper_value.abs() <= 1.0e-12 {
        return raoult_dew_solution(records, vapor_mole_fractions, pressure_pa, upper, 0);
    }
    if lower_value.signum() == upper_value.signum() {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "ternary_raoult_dew_temperature",
            message: format!(
                "pressure {pressure_pa} Pa is not bracketed by the frozen common interval [{lower}, {upper}] K"
            ),
        });
    }
    for iteration in 1..=128 {
        let temperature = 0.5 * (lower + upper);
        let value = raoult_dew_residual(records, temperature, vapor_mole_fractions, pressure_pa)?;
        if value.abs() <= 1.0e-12 || (upper - lower) <= 1.0e-10 {
            return raoult_dew_solution(
                records,
                vapor_mole_fractions,
                pressure_pa,
                temperature,
                iteration,
            );
        }
        if value.signum() == lower_value.signum() {
            lower = temperature;
            lower_value = value;
        } else {
            upper = temperature;
        }
    }
    Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
}

/// Solves an independent ideal-Raoult Rachford-Rice flash for a bulk ternary
/// inventory. It shares only pure `Psat(T)` evidence with the gauge; no
/// production residual, Gibbs closure, or phase-control state participates.
pub(crate) fn solve_raoult_flash(
    records: &[FrozenAntoineStandardState],
    temperature_k: f64,
    pressure_pa: f64,
    bulk_mole_fractions: [f64; 3],
) -> Result<TernaryRaoultFlashSolution, ReactionExtentError> {
    if records.len() != 3 {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "ternary Raoult flash requires three Antoine records, got {}",
            records.len()
        )));
    }
    validate_ternary_simplex("ternary_raoult_flash_bulk", bulk_mole_fractions)?;
    if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_flash_pressure",
            message: "flash pressure must be finite and positive".to_owned(),
        });
    }
    let equilibrium_ratios = records
        .iter()
        .map(|record| Ok(saturation_pressure_pa(record, temperature_k)? / pressure_pa))
        .collect::<Result<Vec<_>, ReactionExtentError>>()?;
    if equilibrium_ratios
        .iter()
        .any(|ratio| !ratio.is_finite() || *ratio <= 0.0)
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_flash_equilibrium_ratios",
            message: "all Raoult equilibrium ratios must be finite and positive".to_owned(),
        });
    }
    let endpoint_tolerance = 1.0e-12;
    let lower_value = rachford_rice_value(&bulk_mole_fractions, &equilibrium_ratios, 0.0)?;
    let upper_value = rachford_rice_value(&bulk_mole_fractions, &equilibrium_ratios, 1.0)?;
    if lower_value <= endpoint_tolerance {
        return Ok(TernaryRaoultFlashSolution {
            phase_class: TernaryRaoultFlashPhaseClass::AllLiquid,
            temperature_k,
            pressure_pa,
            vapor_fraction: None,
            liquid_mole_fractions: Some(bulk_mole_fractions),
            vapor_mole_fractions: None,
            iterations: 0,
        });
    }
    if upper_value >= -endpoint_tolerance {
        return Ok(TernaryRaoultFlashSolution {
            phase_class: TernaryRaoultFlashPhaseClass::AllVapor,
            temperature_k,
            pressure_pa,
            vapor_fraction: None,
            liquid_mole_fractions: None,
            vapor_mole_fractions: Some(bulk_mole_fractions),
            iterations: 0,
        });
    }

    let mut lower = 0.0;
    let mut upper = 1.0;
    let mut lower_value = lower_value;
    for iteration in 1..=128 {
        let vapor_fraction = 0.5 * (lower + upper);
        let value = rachford_rice_value(&bulk_mole_fractions, &equilibrium_ratios, vapor_fraction)?;
        if value.abs() <= endpoint_tolerance || (upper - lower) <= 1.0e-12 {
            return flash_two_phase_solution(
                temperature_k,
                pressure_pa,
                bulk_mole_fractions,
                &equilibrium_ratios,
                vapor_fraction,
                iteration,
            );
        }
        if value.signum() == lower_value.signum() {
            lower = vapor_fraction;
            lower_value = value;
        } else {
            upper = vapor_fraction;
        }
    }
    Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
}

/// Validates that three mole fractions are finite, non-negative, and sum to one.
fn validate_ternary_simplex(
    field: &'static str,
    mole_fractions: [f64; 3],
) -> Result<(), ReactionExtentError> {
    if mole_fractions
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        || (mole_fractions.iter().sum::<f64>() - 1.0).abs() > 1.0e-12
    {
        return Err(ReactionExtentError::InvalidProblem {
            field,
            message: "mole fractions must be finite, non-negative, and sum to one".to_owned(),
        });
    }
    Ok(())
}

/// Evaluates the Rachford-Rice function for a given vapor fraction.
fn rachford_rice_value(
    bulk_mole_fractions: &[f64; 3],
    equilibrium_ratios: &[f64],
    vapor_fraction: f64,
) -> Result<f64, ReactionExtentError> {
    if !vapor_fraction.is_finite() || !(0.0..=1.0).contains(&vapor_fraction) {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_flash_vapor_fraction",
            message: "vapor fraction must lie in the closed interval [0, 1]".to_owned(),
        });
    }
    bulk_mole_fractions
        .iter()
        .zip(equilibrium_ratios)
        .try_fold(0.0, |value, (&bulk, &ratio)| {
            let denominator = 1.0 + vapor_fraction * (ratio - 1.0);
            if !denominator.is_finite() || denominator <= 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "ternary_raoult_flash_denominator",
                    message: "Rachford-Rice denominator must be finite and positive".to_owned(),
                });
            }
            Ok(value + bulk * (ratio - 1.0) / denominator)
        })
}

/// Builds a validated two-phase flash result and checks mass reconstruction.
fn flash_two_phase_solution(
    temperature_k: f64,
    pressure_pa: f64,
    bulk_mole_fractions: [f64; 3],
    equilibrium_ratios: &[f64],
    vapor_fraction: f64,
    iterations: usize,
) -> Result<TernaryRaoultFlashSolution, ReactionExtentError> {
    let mut liquid_mole_fractions = [0.0; 3];
    let mut vapor_mole_fractions = [0.0; 3];
    for index in 0..3 {
        liquid_mole_fractions[index] =
            bulk_mole_fractions[index] / (1.0 + vapor_fraction * (equilibrium_ratios[index] - 1.0));
        vapor_mole_fractions[index] = equilibrium_ratios[index] * liquid_mole_fractions[index];
    }
    validate_ternary_simplex(
        "ternary_raoult_flash_liquid_composition",
        liquid_mole_fractions,
    )?;
    validate_ternary_simplex(
        "ternary_raoult_flash_vapor_composition",
        vapor_mole_fractions,
    )?;
    for index in 0..3 {
        let reconstructed = vapor_fraction * vapor_mole_fractions[index]
            + (1.0 - vapor_fraction) * liquid_mole_fractions[index];
        if (reconstructed - bulk_mole_fractions[index]).abs() > 1.0e-10 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "ternary_raoult_flash_conservation",
                message: format!(
                    "component {index} reconstructs {reconstructed}, expected {}",
                    bulk_mole_fractions[index]
                ),
            });
        }
    }
    Ok(TernaryRaoultFlashSolution {
        phase_class: TernaryRaoultFlashPhaseClass::TwoPhase,
        temperature_k,
        pressure_pa,
        vapor_fraction: Some(vapor_fraction),
        liquid_mole_fractions: Some(liquid_mole_fractions),
        vapor_mole_fractions: Some(vapor_mole_fractions),
        iterations,
    })
}

/// Evaluates the ideal-Raoult dew residual `sum_i y_i * P / Psat_i(T) - 1`.
fn raoult_dew_residual(
    records: &[FrozenAntoineStandardState],
    temperature_k: f64,
    vapor_mole_fractions: [f64; 3],
    pressure_pa: f64,
) -> Result<f64, ReactionExtentError> {
    records
        .iter()
        .zip(vapor_mole_fractions)
        .try_fold(0.0, |residual, (record, vapor)| {
            Ok(residual + vapor * pressure_pa / saturation_pressure_pa(record, temperature_k)?)
        })
        .map(|sum| sum - 1.0)
}

/// Assembles a bubble solution with derived vapor composition at a root.
fn raoult_solution(
    records: &[FrozenAntoineStandardState],
    liquid_mole_fractions: [f64; 3],
    pressure_pa: f64,
    temperature_k: f64,
    iterations: usize,
) -> Result<TernaryRaoultBubbleSolution, ReactionExtentError> {
    let mut vapor_mole_fractions = [0.0; 3];
    for (index, record) in records.iter().enumerate() {
        vapor_mole_fractions[index] = liquid_mole_fractions[index]
            * saturation_pressure_pa(record, temperature_k)?
            / pressure_pa;
    }
    if vapor_mole_fractions
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        || (vapor_mole_fractions.iter().sum::<f64>() - 1.0).abs() > 1.0e-10
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ternary_raoult_vapor_composition",
            message: "derived vapor composition must be finite, strictly interior, and sum to one"
                .to_owned(),
        });
    }
    Ok(TernaryRaoultBubbleSolution {
        temperature_k,
        pressure_pa,
        vapor_mole_fractions,
        iterations,
    })
}

/// Assembles a dew solution with incipient liquid composition at a root.
fn raoult_dew_solution(
    records: &[FrozenAntoineStandardState],
    vapor_mole_fractions: [f64; 3],
    pressure_pa: f64,
    temperature_k: f64,
    iterations: usize,
) -> Result<TernaryRaoultDewSolution, ReactionExtentError> {
    let mut liquid_mole_fractions = [0.0; 3];
    for (index, record) in records.iter().enumerate() {
        liquid_mole_fractions[index] = vapor_mole_fractions[index] * pressure_pa
            / saturation_pressure_pa(record, temperature_k)?;
    }
    validate_ternary_simplex(
        "ternary_raoult_dew_liquid_composition",
        liquid_mole_fractions,
    )?;
    Ok(TernaryRaoultDewSolution {
        temperature_k,
        pressure_pa,
        liquid_mole_fractions,
        iterations,
    })
}
