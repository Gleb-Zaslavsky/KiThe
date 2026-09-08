//! Typed evidence and comparison policy for the frozen IAPWS water table.
//!
//! This module is intentionally separate from the generic frozen-data loader.
//! The loader validates files and provenance; this module states what one
//! particular physical comparison means for the current KiThe model.

use std::fmt;

use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, WaterSaturationPressureReference,
};

/// Root evidence emitted by either independent or canonical boundary search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct WaterSaturationRoot {
    /// Total equilibrium pressure at the converged boundary root, Pa.
    pub(crate) total_pressure_pa: f64,
    /// Partial pressure of water vapor, Pa.
    pub(crate) partial_water_pressure_pa: f64,
    /// Residual magnitude at the converged root.
    pub(crate) residual: f64,
    /// Number of iterations taken by the boundary solver.
    pub(crate) iterations: usize,
}

/// One temperature row in the IAPWS-versus-KiThe characterization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct WaterSaturationComparisonRow {
    /// Frozen saturation temperature, K.
    pub(crate) temperature_k: f64,
    /// Published IAPWS saturation pressure at this temperature, Pa.
    pub(crate) iapws_pressure_pa: f64,
    /// Root from an independent I1/I2 saturation-pressure equation.
    pub(crate) independent: WaterSaturationRoot,
    /// Root from the canonical TPD phase-boundary path.
    pub(crate) canonical_tpd: WaterSaturationRoot,
}

impl WaterSaturationComparisonRow {
    /// Relative I1/I2 deviation from the frozen IAPWS pressure.
    pub(crate) fn independent_external_error(&self) -> f64 {
        relative_delta(
            self.independent.partial_water_pressure_pa,
            self.iapws_pressure_pa,
        )
    }

    /// Relative canonical-TPD deviation from the frozen IAPWS pressure.
    pub(crate) fn canonical_external_error(&self) -> f64 {
        relative_delta(
            self.canonical_tpd.partial_water_pressure_pa,
            self.iapws_pressure_pa,
        )
    }

    /// Relative difference between the two internal boundary roots.
    pub(crate) fn internal_root_error(&self) -> f64 {
        relative_delta(
            self.canonical_tpd.partial_water_pressure_pa,
            self.independent.partial_water_pressure_pa,
        )
    }
}

/// Structured output of the frozen-IAPWS boundary diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct WaterSaturationComparisonReport {
    dataset_id: String,
    source_identifier: Option<String>,
    rows: Vec<WaterSaturationComparisonRow>,
}

/// One three-way internal boundary comparison at a fixed temperature.
///
/// This deliberately omits IAPWS values: it proves symmetry between two
/// canonical phase-stability views and the independent local I1/I2 equation.
/// The external comparison remains owned by [`WaterSaturationComparisonReport`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct WaterSaturationSymmetryRow {
    /// Fixed temperature of the three-way boundary comparison, K.
    pub(crate) temperature_k: f64,
    /// Root from an independent I1/I2 saturation-pressure equation.
    pub(crate) independent: WaterSaturationRoot,
    /// Root solving the liquid side against the gas phase.
    pub(crate) liquid_from_gas: WaterSaturationRoot,
    /// Root solving the gas side against the liquid phase.
    pub(crate) gas_from_liquid: WaterSaturationRoot,
}

impl WaterSaturationSymmetryRow {
    /// Relative I1/I2-to-liquid-TPD root deviation.
    pub(crate) fn independent_liquid_error(&self) -> f64 {
        relative_delta(
            self.liquid_from_gas.partial_water_pressure_pa,
            self.independent.partial_water_pressure_pa,
        )
    }

    /// Relative I1/I2-to-gas-TPD root deviation.
    pub(crate) fn independent_gas_error(&self) -> f64 {
        relative_delta(
            self.gas_from_liquid.partial_water_pressure_pa,
            self.independent.partial_water_pressure_pa,
        )
    }

    /// Relative liquid-side-to-gas-side canonical root deviation.
    pub(crate) fn canonical_symmetry_error(&self) -> f64 {
        relative_delta(
            self.gas_from_liquid.partial_water_pressure_pa,
            self.liquid_from_gas.partial_water_pressure_pa,
        )
    }
}

/// Typed internal-evidence report for the two sides of the water boundary.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct WaterSaturationSymmetryReport {
    rows: Vec<WaterSaturationSymmetryRow>,
}

impl WaterSaturationSymmetryReport {
    /// Builds a report after rejecting nonphysical root evidence.
    pub(crate) fn new(rows: Vec<WaterSaturationSymmetryRow>) -> Result<Self, String> {
        if rows.is_empty() {
            return Err("water saturation symmetry report requires at least one row".to_owned());
        }
        let mut previous_temperature = None;
        for (index, row) in rows.iter().enumerate() {
            if !row.temperature_k.is_finite()
                || row.temperature_k <= 0.0
                || previous_temperature.is_some_and(|previous| row.temperature_k <= previous)
            {
                return Err(format!(
                    "water saturation symmetry row {index} has an invalid temperature {} K",
                    row.temperature_k
                ));
            }
            previous_temperature = Some(row.temperature_k);
            for (label, root) in [
                ("independent", row.independent),
                ("liquid-from-gas", row.liquid_from_gas),
                ("gas-from-liquid", row.gas_from_liquid),
            ] {
                for (value_label, value, must_be_positive) in [
                    ("total pressure", root.total_pressure_pa, true),
                    (
                        "water partial pressure",
                        root.partial_water_pressure_pa,
                        true,
                    ),
                    ("residual", root.residual, false),
                ] {
                    if !value.is_finite() || (must_be_positive && value <= 0.0) {
                        return Err(format!(
                            "water saturation symmetry row {index} has invalid {label} {value_label}"
                        ));
                    }
                }
            }
        }
        Ok(Self { rows })
    }

    /// Rows in the source temperature order used by the diagnostic.
    pub(crate) fn rows(&self) -> &[WaterSaturationSymmetryRow] {
        &self.rows
    }

    /// The largest absolute deviation across all pairwise internal roots.
    pub(crate) fn max_internal_relative_error(&self) -> f64 {
        self.rows
            .iter()
            .flat_map(|row| {
                [
                    row.independent_liquid_error(),
                    row.independent_gas_error(),
                    row.canonical_symmetry_error(),
                ]
            })
            .map(f64::abs)
            .fold(0.0_f64, f64::max)
    }

    /// Enforces one strict tolerance for all three descriptions of the same
    /// physical phase boundary.
    pub(crate) fn validate_internal_roots(&self, max_relative_delta: f64) -> Result<(), String> {
        if !max_relative_delta.is_finite() || max_relative_delta <= 0.0 {
            return Err("internal-root tolerance must be finite and positive".to_owned());
        }
        for (index, row) in self.rows.iter().enumerate() {
            for (label, error) in [
                ("I1/I2-vs-liquid-TPD", row.independent_liquid_error()),
                ("I1/I2-vs-gas-TPD", row.independent_gas_error()),
                ("liquid-TPD-vs-gas-TPD", row.canonical_symmetry_error()),
            ] {
                if error.abs() > max_relative_delta {
                    return Err(format!(
                        "water saturation symmetry row {index} at {} K has {label} relative delta {error:e}, exceeding {max_relative_delta:e}",
                        row.temperature_k
                    ));
                }
            }
        }
        Ok(())
    }
}

impl WaterSaturationComparisonReport {
    /// Collects one validated diagnostic row for every frozen reference row.
    pub(crate) fn new(
        dataset: &FrozenReferenceDataset<WaterSaturationPressureReference>,
        rows: Vec<WaterSaturationComparisonRow>,
    ) -> Result<Self, String> {
        if rows.len() != dataset.rows().len() {
            return Err(format!(
                "IAPWS comparison has {} rows for {} frozen reference rows",
                rows.len(),
                dataset.rows().len()
            ));
        }
        for (index, (comparison, reference)) in rows.iter().zip(dataset.rows()).enumerate() {
            if comparison.temperature_k != reference.temperature_k
                || comparison.iapws_pressure_pa != reference.pressure_pa
            {
                return Err(format!(
                    "IAPWS comparison row {index} does not preserve frozen temperature/pressure identity"
                ));
            }
            for (label, value) in [
                (
                    "independent total pressure",
                    comparison.independent.total_pressure_pa,
                ),
                (
                    "independent water partial pressure",
                    comparison.independent.partial_water_pressure_pa,
                ),
                (
                    "canonical total pressure",
                    comparison.canonical_tpd.total_pressure_pa,
                ),
                (
                    "canonical water partial pressure",
                    comparison.canonical_tpd.partial_water_pressure_pa,
                ),
                ("independent residual", comparison.independent.residual),
                ("canonical residual", comparison.canonical_tpd.residual),
            ] {
                if !value.is_finite() {
                    return Err(format!(
                        "IAPWS comparison row {index} has non-finite {label}"
                    ));
                }
            }
        }

        Ok(Self {
            dataset_id: dataset.metadata().dataset_id.clone(),
            source_identifier: dataset.metadata().source.stable_identifier.clone(),
            rows,
        })
    }

    /// Frozen dataset identity this report characterizes.
    pub(crate) fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Row-level diagnostic evidence in frozen source order.
    pub(crate) fn rows(&self) -> &[WaterSaturationComparisonRow] {
        &self.rows
    }

    /// Largest absolute external relative error across all rows and both roots.
    pub(crate) fn max_external_relative_error(&self) -> f64 {
        self.rows
            .iter()
            .flat_map(|row| {
                [
                    row.independent_external_error(),
                    row.canonical_external_error(),
                ]
            })
            .map(f64::abs)
            .fold(0.0_f64, f64::max)
    }

    /// Root-mean-square of the independent-root external relative errors.
    pub(crate) fn independent_external_rms(&self) -> f64 {
        root_mean_square(self.rows.iter().map(|row| row.independent_external_error()))
    }

    /// Root-mean-square of the canonical-TPD external relative errors.
    pub(crate) fn canonical_external_rms(&self) -> f64 {
        root_mean_square(self.rows.iter().map(|row| row.canonical_external_error()))
    }

    /// Largest I1/I2-versus-TPD relative deviation across all rows.
    pub(crate) fn max_internal_root_relative_error(&self) -> f64 {
        self.rows
            .iter()
            .map(|row| row.internal_root_error().abs())
            .fold(0.0_f64, f64::max)
    }

    /// Enforces the numerical independence contract, which is intentionally
    /// stricter than the still-unreviewed external model comparison.
    pub(crate) fn validate_internal_roots(&self, max_relative_delta: f64) -> Result<(), String> {
        if !max_relative_delta.is_finite() || max_relative_delta <= 0.0 {
            return Err("internal-root tolerance must be finite and positive".to_owned());
        }
        for (index, row) in self.rows.iter().enumerate() {
            let error = row.internal_root_error();
            if error.abs() > max_relative_delta {
                return Err(format!(
                    "IAPWS row {index} at {} K has I1/I2-vs-TPD relative delta {error:e}, exceeding {max_relative_delta:e}",
                    row.temperature_k
                ));
            }
        }
        Ok(())
    }
}

/// How this particular external table may be used by the current model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExternalComparisonMode {
    /// Record model discrepancy without treating it as a solver acceptance bound.
    CharacterizationOnly,
}

/// Source- and model-specific policy for the IAPWS water comparison.
///
/// It deliberately does not contain generic nonlinear-solver tolerances. The
/// present NASA/ideal-phase model is not an implementation of IAPWS SR1-86.
/// In this deliberately low-pressure range, the observed external difference
/// is recorded as a model/thermochemistry discrepancy: likely contributors
/// include consistency of the local gas/condensed standard Gibbs functions,
/// reference-state and polynomial approximation effects, and to a lesser
/// degree ideal-gas versus real-fluid behavior. This contract does not assign
/// quantitative responsibility to any one contributor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IapwsWaterSaturationComparisonContract {
    expected_dataset_id: &'static str,
    expected_source_identifier: &'static str,
    model_scope: &'static str,
    mode: ExternalComparisonMode,
}

impl IapwsWaterSaturationComparisonContract {
    pub(crate) const fn characterization_only() -> Self {
        Self {
            expected_dataset_id: "iapws.water_liquid_saturation.low_pressure.v1",
            expected_source_identifier: "IAPWS-SR1-86(1992)",
            model_scope: "local NASA gas + NASA condensed thermochemistry; ideal-gas activity; pure condensed liquid",
            mode: ExternalComparisonMode::CharacterizationOnly,
        }
    }

    pub(crate) fn mode(&self) -> ExternalComparisonMode {
        self.mode
    }

    pub(crate) fn model_scope(&self) -> &str {
        self.model_scope
    }

    /// Refuses to apply an IAPWS-specific interpretation to another dataset.
    pub(crate) fn validate_dataset(
        &self,
        dataset: &FrozenReferenceDataset<WaterSaturationPressureReference>,
    ) -> Result<(), String> {
        if dataset.metadata().dataset_id != self.expected_dataset_id {
            return Err(format!(
                "IAPWS comparison contract expects dataset '{}', got '{}'",
                self.expected_dataset_id,
                dataset.metadata().dataset_id
            ));
        }
        if dataset.metadata().source.stable_identifier.as_deref()
            != Some(self.expected_source_identifier)
        {
            return Err(format!(
                "IAPWS comparison contract expects source '{}', got {:?}",
                self.expected_source_identifier,
                dataset.metadata().source.stable_identifier
            ));
        }
        Ok(())
    }

    /// Validates policy identity without turning a characterization into a
    /// synthetic external-accuracy threshold.
    pub(crate) fn validate_report(
        &self,
        report: &WaterSaturationComparisonReport,
    ) -> Result<(), String> {
        if report.dataset_id != self.expected_dataset_id
            || report.source_identifier.as_deref() != Some(self.expected_source_identifier)
        {
            return Err(
                "IAPWS comparison report provenance does not match its contract".to_owned(),
            );
        }
        Ok(())
    }
}

impl fmt::Display for WaterSaturationComparisonReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            formatter,
            "{:<7} {:>13} {:>13} {:>13} {:>12} {:>12} {:>12}",
            "T [K]",
            "IAPWS [Pa]",
            "I1/I2 pH2O",
            "TPD pH2O",
            "I1 err [%]",
            "TPD err [%]",
            "internal [%]"
        )?;
        for row in &self.rows {
            writeln!(
                formatter,
                "{:<7.1} {:>13.3} {:>13.3} {:>13.3} {:>12.5} {:>12.5} {:>12.5}",
                row.temperature_k,
                row.iapws_pressure_pa,
                row.independent.partial_water_pressure_pa,
                row.canonical_tpd.partial_water_pressure_pa,
                100.0 * row.independent_external_error(),
                100.0 * row.canonical_external_error(),
                100.0 * row.internal_root_error(),
            )?;
        }
        write!(
            formatter,
            "summary: max external relative error={:.6e}, rms I1/I2 external={:.6e}, rms TPD external={:.6e}, max I1/I2-vs-TPD={:.6e}",
            self.max_external_relative_error(),
            self.independent_external_rms(),
            self.canonical_external_rms(),
            self.max_internal_root_relative_error(),
        )
    }
}

impl fmt::Display for WaterSaturationSymmetryReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            formatter,
            "{:<7} {:>13} {:>13} {:>13} {:>13} {:>13} {:>13}",
            "T [K]",
            "I1/I2 pH2O",
            "liq|gas pH2O",
            "gas|liq pH2O",
            "I1-liq [%]",
            "I1-gas [%]",
            "liq-gas [%]",
        )?;
        for row in &self.rows {
            writeln!(
                formatter,
                "{:<7.1} {:>13.3} {:>13.3} {:>13.3} {:>13.5} {:>13.5} {:>13.5}",
                row.temperature_k,
                row.independent.partial_water_pressure_pa,
                row.liquid_from_gas.partial_water_pressure_pa,
                row.gas_from_liquid.partial_water_pressure_pa,
                100.0 * row.independent_liquid_error(),
                100.0 * row.independent_gas_error(),
                100.0 * row.canonical_symmetry_error(),
            )?;
        }
        write!(
            formatter,
            "summary: max three-way internal relative error={:.6e}",
            self.max_internal_relative_error(),
        )
    }
}

/// Relative deviation of an actual value from a reference (`(a-r)/r`).
fn relative_delta(actual: f64, reference: f64) -> f64 {
    (actual - reference) / reference
}

/// Root-mean-square of an iterable of values.
fn root_mean_square(values: impl IntoIterator<Item = f64>) -> f64 {
    let values = values.into_iter().collect::<Vec<_>>();
    debug_assert!(!values.is_empty());
    (values.iter().map(|value| value * value).sum::<f64>() / values.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comparison_row_keeps_external_and_internal_errors_distinct() {
        let row = WaterSaturationComparisonRow {
            temperature_k: 300.0,
            iapws_pressure_pa: 100.0,
            independent: WaterSaturationRoot {
                total_pressure_pa: 200.0,
                partial_water_pressure_pa: 99.0,
                residual: 0.0,
                iterations: 10,
            },
            canonical_tpd: WaterSaturationRoot {
                total_pressure_pa: 200.0,
                partial_water_pressure_pa: 99.000_001,
                residual: 0.0,
                iterations: 11,
            },
        };
        assert!((row.independent_external_error() + 0.01).abs() < 1e-12);
        assert!(row.internal_root_error().abs() < 2e-8);
    }

    #[test]
    fn symmetry_report_keeps_all_three_internal_root_axes() {
        let root = WaterSaturationRoot {
            total_pressure_pa: 100.0,
            partial_water_pressure_pa: 100.0,
            residual: 0.0,
            iterations: 5,
        };
        let report = WaterSaturationSymmetryReport::new(vec![WaterSaturationSymmetryRow {
            temperature_k: 300.0,
            independent: root,
            liquid_from_gas: root,
            gas_from_liquid: root,
        }])
        .unwrap();

        assert_eq!(report.rows().len(), 1);
        assert_eq!(report.max_internal_relative_error(), 0.0);
        report.validate_internal_roots(1e-12).unwrap();
        assert!(report.to_string().contains("liq|gas pH2O"));
    }
}
