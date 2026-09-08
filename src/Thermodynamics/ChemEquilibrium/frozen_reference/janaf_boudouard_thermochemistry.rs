//! Typed I5 comparison for JANAF Boudouard reaction thermochemistry.
//!
//! This module is deliberately below phase-boundary validation. It compares
//! only standard reaction Gibbs energy and `log10(Kp)` for
//! `2 CO(g) <=> CO2(g) + C(gr)`, using frozen primary JANAF species values and
//! local KiThe standard-state closures. TPD, active sets, pressure boundaries,
//! and `P,H` lifecycle are explicitly outside this first layer.

use std::fmt;

use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, JanafBoudouardReference,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::MOLAR_GAS_CONSTANT;

/// JANAF standard pressure used by C-093, C-095, and C-002.
pub(crate) const JANAF_STANDARD_PRESSURE_PA: f64 = 100_000.0;

/// Reaction quantities derived solely from one frozen JANAF species row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct JanafBoudouardDerived {
    /// Standard reaction Gibbs energy, J/mol.
    pub(crate) delta_r_g_j_mol: f64,
    /// `log10(Kp)` computed from the Gibbs route.
    pub(crate) log10_kp_from_gibbs: f64,
    /// `log10(Kp)` computed from the published log-Kf route.
    pub(crate) log10_kp_from_log_kf: f64,
}

impl JanafBoudouardDerived {
    /// Derives reaction data from primary CO/CO2 columns. Graphite is JANAF's
    /// reference state, so both of its formation quantities are exactly zero.
    pub(crate) fn from_reference(row: JanafBoudouardReference) -> Self {
        let delta_r_g_j_mol = (row.co2_delta_f_g_kj_mol - 2.0 * row.co_delta_f_g_kj_mol) * 1_000.0;
        let log10_kp_from_gibbs =
            -delta_r_g_j_mol / (MOLAR_GAS_CONSTANT * row.temperature_k * std::f64::consts::LN_10);
        let log10_kp_from_log_kf = row.co2_log10_kf - 2.0 * row.co_log10_kf;
        Self {
            delta_r_g_j_mol,
            log10_kp_from_gibbs,
            log10_kp_from_log_kf,
        }
    }

    /// Difference between the two frozen JANAF `log10(Kp)` representations;
    /// used to prove the frozen table is internally self-consistent.
    pub(crate) fn log10_route_delta(self) -> f64 {
        self.log10_kp_from_gibbs - self.log10_kp_from_log_kf
    }
}

/// One comparison of frozen JANAF and local KiThe reaction thermochemistry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct JanafBoudouardThermochemistryRow {
    /// Temperature of the comparison, K.
    pub(crate) temperature_k: f64,
    /// Reaction quantities derived from the frozen JANAF species row.
    pub(crate) janaf: JanafBoudouardDerived,
    /// Local KiThe standard reaction Gibbs energy, J/mol.
    pub(crate) kithe_delta_r_g_j_mol: f64,
    /// Local KiThe `log10(Kp)`.
    pub(crate) kithe_log10_kp: f64,
}

impl JanafBoudouardThermochemistryRow {
    /// Absolute error of the local reaction Gibbs energy versus JANAF, J/mol.
    pub(crate) fn delta_g_error_j_mol(&self) -> f64 {
        self.kithe_delta_r_g_j_mol - self.janaf.delta_r_g_j_mol
    }

    /// Absolute error of the local `log10(Kp)` versus the JANAF Gibbs route.
    pub(crate) fn log10_k_error(&self) -> f64 {
        self.kithe_log10_kp - self.janaf.log10_kp_from_gibbs
    }
}

/// Typed output of the I4-local to I5-JANAF characterization.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct JanafBoudouardThermochemistryReport {
    dataset_id: String,
    source_identifier: Option<String>,
    standard_pressure_pa: f64,
    rows: Vec<JanafBoudouardThermochemistryRow>,
}

impl JanafBoudouardThermochemistryReport {
    /// Keeps local row order aligned with the immutable JANAF source table.
    pub(crate) fn new(
        dataset: &FrozenReferenceDataset<JanafBoudouardReference>,
        rows: Vec<JanafBoudouardThermochemistryRow>,
        standard_pressure_pa: f64,
    ) -> Result<Self, String> {
        if !standard_pressure_pa.is_finite() || standard_pressure_pa <= 0.0 {
            return Err(
                "JANAF comparison standard pressure must be finite and positive".to_owned(),
            );
        }
        if rows.len() != dataset.rows().len() {
            return Err(format!(
                "JANAF Boudouard comparison has {} rows for {} frozen reference rows",
                rows.len(),
                dataset.rows().len()
            ));
        }
        for (index, (comparison, reference)) in rows.iter().zip(dataset.rows()).enumerate() {
            if comparison.temperature_k != reference.temperature_k {
                return Err(format!(
                    "JANAF Boudouard comparison row {index} does not preserve frozen temperature identity"
                ));
            }
            for (field, value) in [
                ("JANAF reaction Gibbs", comparison.janaf.delta_r_g_j_mol),
                (
                    "JANAF Gibbs-derived log10 Kp",
                    comparison.janaf.log10_kp_from_gibbs,
                ),
                (
                    "JANAF log-Kf-derived log10 Kp",
                    comparison.janaf.log10_kp_from_log_kf,
                ),
                ("KiThe reaction Gibbs", comparison.kithe_delta_r_g_j_mol),
                ("KiThe log10 Kp", comparison.kithe_log10_kp),
            ] {
                if !value.is_finite() {
                    return Err(format!(
                        "JANAF Boudouard row {index} has non-finite {field}"
                    ));
                }
            }
        }
        Ok(Self {
            dataset_id: dataset.metadata().dataset_id.clone(),
            source_identifier: dataset.metadata().source.stable_identifier.clone(),
            standard_pressure_pa,
            rows,
        })
    }

    pub(crate) fn rows(&self) -> &[JanafBoudouardThermochemistryRow] {
        &self.rows
    }

    /// Largest internal JANAF Gibbs-versus-log-Kf `log10(Kp)` discrepancy.
    pub(crate) fn max_janaf_log10_route_delta(&self) -> f64 {
        self.rows
            .iter()
            .map(|row| row.janaf.log10_route_delta().abs())
            .fold(0.0_f64, f64::max)
    }

    /// Largest absolute reaction-Gibbs error versus JANAF, J/mol.
    pub(crate) fn max_delta_g_error_j_mol(&self) -> f64 {
        self.rows
            .iter()
            .map(|row| row.delta_g_error_j_mol().abs())
            .fold(0.0_f64, f64::max)
    }

    /// Root-mean-square of the reaction-Gibbs errors versus JANAF, J/mol.
    pub(crate) fn rms_delta_g_error_j_mol(&self) -> f64 {
        root_mean_square(
            self.rows
                .iter()
                .map(JanafBoudouardThermochemistryRow::delta_g_error_j_mol),
        )
    }

    /// Signed mean of the reaction-Gibbs errors versus JANAF (bias indicator).
    pub(crate) fn mean_signed_delta_g_error_j_mol(&self) -> f64 {
        mean(
            self.rows
                .iter()
                .map(JanafBoudouardThermochemistryRow::delta_g_error_j_mol),
        )
    }

    /// Largest absolute `log10(Kp)` error versus JANAF.
    pub(crate) fn max_log10_k_error(&self) -> f64 {
        self.rows
            .iter()
            .map(|row| row.log10_k_error().abs())
            .fold(0.0_f64, f64::max)
    }

    /// Root-mean-square of the `log10(Kp)` errors versus JANAF.
    pub(crate) fn rms_log10_k_error(&self) -> f64 {
        root_mean_square(
            self.rows
                .iter()
                .map(JanafBoudouardThermochemistryRow::log10_k_error),
        )
    }

    /// Signed mean of the `log10(Kp)` errors versus JANAF (bias indicator).
    pub(crate) fn mean_signed_log10_k_error(&self) -> f64 {
        mean(
            self.rows
                .iter()
                .map(JanafBoudouardThermochemistryRow::log10_k_error),
        )
    }

    /// Validates the two frozen JANAF representations only. It is deliberately
    /// separate from the external KiThe comparison, whose tolerance remains
    /// `CharacterizationOnly` until a thermochemistry review is complete.
    pub(crate) fn validate_janaf_internal_routes(
        &self,
        max_log10_delta: f64,
    ) -> Result<(), String> {
        if !max_log10_delta.is_finite() || max_log10_delta <= 0.0 {
            return Err("JANAF log10-K route tolerance must be finite and positive".to_owned());
        }
        for (index, row) in self.rows.iter().enumerate() {
            let delta = row.janaf.log10_route_delta();
            if delta.abs() > max_log10_delta {
                return Err(format!(
                    "JANAF Boudouard row {index} at {} K has Gibbs-vs-logKf delta {delta:e}, exceeding {max_log10_delta:e}",
                    row.temperature_k
                ));
            }
        }
        Ok(())
    }
}

/// Source/model identity for the first Boudouard I5 comparison.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct JanafBoudouardComparisonContract {
    expected_dataset_id: &'static str,
    expected_source_identifier: &'static str,
    standard_pressure_pa: f64,
}

impl JanafBoudouardComparisonContract {
    pub(crate) const fn characterization_only() -> Self {
        Self {
            expected_dataset_id: "janaf.boudouard.reaction_thermodynamics.v1",
            expected_source_identifier: "NIST-JANAF-Boudouard-C093-C095-C002",
            standard_pressure_pa: JANAF_STANDARD_PRESSURE_PA,
        }
    }

    pub(crate) fn standard_pressure_pa(&self) -> f64 {
        self.standard_pressure_pa
    }

    pub(crate) fn validate_dataset(
        &self,
        dataset: &FrozenReferenceDataset<JanafBoudouardReference>,
    ) -> Result<(), String> {
        if dataset.metadata().dataset_id != self.expected_dataset_id {
            return Err(format!(
                "JANAF Boudouard contract expects dataset '{}', got '{}'",
                self.expected_dataset_id,
                dataset.metadata().dataset_id
            ));
        }
        if dataset.metadata().source.stable_identifier.as_deref()
            != Some(self.expected_source_identifier)
        {
            return Err(format!(
                "JANAF Boudouard contract expects source '{}', got {:?}",
                self.expected_source_identifier,
                dataset.metadata().source.stable_identifier
            ));
        }
        Ok(())
    }

    pub(crate) fn validate_report(
        &self,
        report: &JanafBoudouardThermochemistryReport,
    ) -> Result<(), String> {
        if report.dataset_id != self.expected_dataset_id
            || report.source_identifier.as_deref() != Some(self.expected_source_identifier)
            || report.standard_pressure_pa != self.standard_pressure_pa
        {
            return Err(
                "JANAF Boudouard report provenance or pressure convention mismatches its contract"
                    .to_owned(),
            );
        }
        Ok(())
    }
}

impl fmt::Display for JanafBoudouardThermochemistryReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            formatter,
            "{:<7} {:>15} {:>15} {:>13} {:>13} {:>13} {:>13}",
            "T [K]",
            "JANAF dGr [kJ]",
            "KiThe dGr [kJ]",
            "dG err [J]",
            "JANAF log10K",
            "KiThe log10K",
            "dlog10K",
        )?;
        for row in &self.rows {
            writeln!(
                formatter,
                "{:<7.1} {:>15.6} {:>15.6} {:>13.3} {:>13.6} {:>13.6} {:>13.6}",
                row.temperature_k,
                row.janaf.delta_r_g_j_mol / 1_000.0,
                row.kithe_delta_r_g_j_mol / 1_000.0,
                row.delta_g_error_j_mol(),
                row.janaf.log10_kp_from_gibbs,
                row.kithe_log10_kp,
                row.log10_k_error(),
            )?;
        }
        write!(
            formatter,
            "summary: p0={:.0} Pa, max |JANAF route dlog10K|={:.6e}, max |dG|={:.6e} J/mol, rms dG={:.6e} J/mol, mean signed dG={:.6e} J/mol, max |dlog10K|={:.6e}, rms dlog10K={:.6e}, mean signed dlog10K={:.6e}",
            self.standard_pressure_pa,
            self.max_janaf_log10_route_delta(),
            self.max_delta_g_error_j_mol(),
            self.rms_delta_g_error_j_mol(),
            self.mean_signed_delta_g_error_j_mol(),
            self.max_log10_k_error(),
            self.rms_log10_k_error(),
            self.mean_signed_log10_k_error(),
        )
    }
}

/// Root-mean-square of an iterable of values.
fn root_mean_square(values: impl Iterator<Item = f64>) -> f64 {
    let (sum_squares, count) = values.fold((0.0_f64, 0_usize), |(sum, count), value| {
        (sum + value * value, count + 1)
    });
    (sum_squares / count as f64).sqrt()
}

/// Arithmetic mean of an iterable of values.
fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let (sum, count) = values.fold((0.0_f64, 0_usize), |(sum, count), value| {
        (sum + value, count + 1)
    });
    sum / count as f64
}

#[cfg(test)]
mod tests {
    use super::{JANAF_STANDARD_PRESSURE_PA, JanafBoudouardDerived};
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::JanafBoudouardReference;

    #[test]
    fn derived_janaf_reaction_uses_primary_species_columns() {
        let row = JanafBoudouardReference {
            temperature_k: 1_000.0,
            co_delta_f_g_kj_mol: -200.275,
            co_log10_kf: 10.461,
            co2_delta_f_g_kj_mol: -395.886,
            co2_log10_kf: 20.679,
        };
        let derived = JanafBoudouardDerived::from_reference(row);
        assert!((derived.delta_r_g_j_mol - 4_664.0).abs() <= 1e-9);
        assert!((derived.log10_kp_from_log_kf + 0.243).abs() <= 1e-12);
        assert!(derived.log10_route_delta().abs() < 0.005);
        assert_eq!(JANAF_STANDARD_PRESSURE_PA, 100_000.0);
    }
}
