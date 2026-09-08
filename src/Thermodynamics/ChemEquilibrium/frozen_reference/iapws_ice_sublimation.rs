//! Typed evidence for the frozen IAPWS ice-Ih sublimation-pressure table.
//!
//! The generic frozen loader owns JSON parsing and provenance validation. This
//! module owns only the physical meaning of the `H2O(s, ice Ih) <=> H2O(g)`
//! comparison: three independent internal roots must agree strictly, while
//! the current local NASA/ideal-phase model is characterized against IAPWS
//! without an invented external-accuracy threshold.

use std::fmt;

use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, WaterIceSublimationPressureReference,
};

/// Root evidence produced by an independent or phase-stability boundary path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct IceSublimationRoot {
    /// Sublimation pressure found by this independent boundary path, Pa.
    pub(crate) pressure_pa: f64,
    /// Residual magnitude at the converged root.
    pub(crate) residual: f64,
    /// Number of iterations taken by the boundary solver.
    pub(crate) iterations: usize,
}

/// All independent descriptions of one ice-Ih sublimation boundary point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct IceSublimationComparisonRow {
    /// Frozen temperature of the ice-Ih sublimation point, K.
    pub(crate) temperature_k: f64,
    /// Published IAPWS sublimation pressure at this temperature, Pa.
    pub(crate) iapws_pressure_pa: f64,
    /// Root from an independent sublimation-pressure equation.
    pub(crate) independent: IceSublimationRoot,
    /// Root from a phase-boundary path solving the ice side against the gas.
    pub(crate) ice_from_gas: IceSublimationRoot,
    /// Root from a phase-boundary path solving the gas side against the ice.
    pub(crate) gas_from_ice: IceSublimationRoot,
}

impl IceSublimationComparisonRow {
    /// Relative deviation of the independent root from the IAPWS value.
    pub(crate) fn independent_external_error(&self) -> f64 {
        relative_delta(self.independent.pressure_pa, self.iapws_pressure_pa)
    }

    /// Relative deviation of the ice-from-gas TPD root from the IAPWS value.
    pub(crate) fn ice_tpd_external_error(&self) -> f64 {
        relative_delta(self.ice_from_gas.pressure_pa, self.iapws_pressure_pa)
    }

    /// Relative deviation of the gas-from-ice TPD root from the IAPWS value.
    pub(crate) fn gas_tpd_external_error(&self) -> f64 {
        relative_delta(self.gas_from_ice.pressure_pa, self.iapws_pressure_pa)
    }

    /// Largest pairwise relative delta among the three internal roots; used to
    /// prove that independent and TPD boundary descriptions agree strictly.
    pub(crate) fn max_internal_error(&self) -> f64 {
        [
            relative_delta(self.independent.pressure_pa, self.ice_from_gas.pressure_pa),
            relative_delta(self.independent.pressure_pa, self.gas_from_ice.pressure_pa),
            relative_delta(self.ice_from_gas.pressure_pa, self.gas_from_ice.pressure_pa),
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0_f64, f64::max)
    }
}

/// Printable, typed IAPWS ice-Ih characterization evidence.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct IceSublimationComparisonReport {
    dataset_id: String,
    source_identifier: Option<String>,
    rows: Vec<IceSublimationComparisonRow>,
}

impl IceSublimationComparisonReport {
    /// Creates a report only when every frozen row has exactly one finite
    /// internal comparison row in source order.
    pub(crate) fn new(
        dataset: &FrozenReferenceDataset<WaterIceSublimationPressureReference>,
        rows: Vec<IceSublimationComparisonRow>,
    ) -> Result<Self, String> {
        if rows.len() != dataset.rows().len() {
            return Err(format!(
                "IAPWS ice comparison has {} rows for {} frozen reference rows",
                rows.len(),
                dataset.rows().len()
            ));
        }
        for (index, (comparison, reference)) in rows.iter().zip(dataset.rows()).enumerate() {
            if comparison.temperature_k != reference.temperature_k
                || comparison.iapws_pressure_pa != reference.pressure_pa
            {
                return Err(format!(
                    "IAPWS ice comparison row {index} does not preserve frozen temperature/pressure identity"
                ));
            }
            for (label, root) in [
                ("independent", comparison.independent),
                ("ice-from-gas", comparison.ice_from_gas),
                ("gas-from-ice", comparison.gas_from_ice),
            ] {
                if !root.pressure_pa.is_finite()
                    || root.pressure_pa <= 0.0
                    || !root.residual.is_finite()
                {
                    return Err(format!(
                        "IAPWS ice comparison row {index} has invalid {label} root evidence"
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

    pub(crate) fn rows(&self) -> &[IceSublimationComparisonRow] {
        &self.rows
    }

    /// Largest absolute external relative error across all rows and all roots.
    pub(crate) fn max_external_relative_error(&self) -> f64 {
        self.rows
            .iter()
            .flat_map(|row| {
                [
                    row.independent_external_error(),
                    row.ice_tpd_external_error(),
                    row.gas_tpd_external_error(),
                ]
            })
            .map(f64::abs)
            .fold(0.0_f64, f64::max)
    }

    /// Root-mean-square of the independent-root external relative errors.
    pub(crate) fn independent_external_rms(&self) -> f64 {
        root_mean_square(self.rows.iter().map(|row| row.independent_external_error()))
    }

    /// Signed mean makes a systematic local thermochemistry bias visible.
    pub(crate) fn independent_external_mean_signed(&self) -> f64 {
        self.rows
            .iter()
            .map(IceSublimationComparisonRow::independent_external_error)
            .sum::<f64>()
            / self.rows.len() as f64
    }

    /// Largest three-way internal relative error across all rows.
    pub(crate) fn max_internal_relative_error(&self) -> f64 {
        self.rows
            .iter()
            .map(IceSublimationComparisonRow::max_internal_error)
            .fold(0.0_f64, f64::max)
    }

    /// Asserts that every row's three-way internal agreement stays within the
    /// given positive tolerance, returning the first offending row otherwise.
    pub(crate) fn validate_internal_roots(&self, max_relative_delta: f64) -> Result<(), String> {
        if !max_relative_delta.is_finite() || max_relative_delta <= 0.0 {
            return Err("internal-root tolerance must be finite and positive".to_owned());
        }
        for (index, row) in self.rows.iter().enumerate() {
            let error = row.max_internal_error();
            if error > max_relative_delta {
                return Err(format!(
                    "IAPWS ice row {index} at {} K has three-way internal relative delta {error:e}, exceeding {max_relative_delta:e}",
                    row.temperature_k
                ));
            }
        }
        Ok(())
    }
}

/// Source/model-specific IAPWS ice-Ih comparison policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IapwsIceSublimationComparisonContract {
    expected_dataset_id: &'static str,
    expected_source_identifier: &'static str,
}

impl IapwsIceSublimationComparisonContract {
    pub(crate) const fn characterization_only() -> Self {
        Self {
            expected_dataset_id: "iapws.water_ice_ih_sublimation.low_pressure.v1",
            expected_source_identifier: "IAPWS-R14-08(2011)",
        }
    }

    pub(crate) fn validate_dataset(
        &self,
        dataset: &FrozenReferenceDataset<WaterIceSublimationPressureReference>,
    ) -> Result<(), String> {
        if dataset.metadata().dataset_id != self.expected_dataset_id {
            return Err(format!(
                "IAPWS ice comparison expects dataset '{}', got '{}'",
                self.expected_dataset_id,
                dataset.metadata().dataset_id
            ));
        }
        if dataset.metadata().source.stable_identifier.as_deref()
            != Some(self.expected_source_identifier)
        {
            return Err(format!(
                "IAPWS ice comparison expects source '{}', got {:?}",
                self.expected_source_identifier,
                dataset.metadata().source.stable_identifier
            ));
        }
        Ok(())
    }

    pub(crate) fn validate_report(
        &self,
        report: &IceSublimationComparisonReport,
    ) -> Result<(), String> {
        if report.dataset_id != self.expected_dataset_id
            || report.source_identifier.as_deref() != Some(self.expected_source_identifier)
        {
            return Err(
                "IAPWS ice comparison report provenance does not match its contract".to_owned(),
            );
        }
        Ok(())
    }
}

impl fmt::Display for IceSublimationComparisonReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            formatter,
            "{:<7} {:>13} {:>13} {:>13} {:>13} {:>12} {:>12}",
            "T [K]",
            "IAPWS [Pa]",
            "I1/I2 [Pa]",
            "ice|gas [Pa]",
            "gas|ice [Pa]",
            "external [%]",
            "internal [%]",
        )?;
        for row in &self.rows {
            writeln!(
                formatter,
                "{:<7.1} {:>13.6} {:>13.6} {:>13.6} {:>13.6} {:>12.5} {:>12.5}",
                row.temperature_k,
                row.iapws_pressure_pa,
                row.independent.pressure_pa,
                row.ice_from_gas.pressure_pa,
                row.gas_from_ice.pressure_pa,
                100.0 * row.independent_external_error(),
                100.0 * row.max_internal_error(),
            )?;
        }
        write!(
            formatter,
            "summary: max external relative error={:.6e}, rms external={:.6e}, mean signed external={:.6e}, max three-way internal relative error={:.6e}",
            self.max_external_relative_error(),
            self.independent_external_rms(),
            self.independent_external_mean_signed(),
            self.max_internal_relative_error(),
        )
    }
}

/// Relative deviation of an actual value from a reference (`(a-r)/r`).
fn relative_delta(actual: f64, reference: f64) -> f64 {
    (actual - reference) / reference
}

/// Root-mean-square of an iterable of values.
fn root_mean_square(values: impl Iterator<Item = f64>) -> f64 {
    let (sum_squares, count) = values.fold((0.0_f64, 0_usize), |(sum, count), value| {
        (sum + value * value, count + 1)
    });
    (sum_squares / count as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::{IceSublimationComparisonReport, IceSublimationComparisonRow, IceSublimationRoot};
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenReferenceDataset, WaterIceSublimationPressureReference,
    };

    #[test]
    fn ice_report_rejects_nonpositive_root_pressure() {
        let metadata = r#"{
            "dataset_format_version": 1,
            "dataset_id": "iapws.water_ice_ih_sublimation.low_pressure.v1",
            "title": "test",
            "evidence_kind": "frozen_external",
            "data_file": "test.json",
            "source": {
                "organization": "IAPWS", "name": "test", "citation": "test",
                "stable_identifier": "IAPWS-R14-08(2011)"
            },
            "transcription": "test",
            "quantities": [
                { "column": "temperature_k", "meaning": "T", "unit": "K" },
                { "column": "pressure_pa", "meaning": "p", "unit": "Pa" }
            ]
        }"#;
        let rows = r#"{
            "dataset_id": "iapws.water_ice_ih_sublimation.low_pressure.v1",
            "rows": [{ "temperature_k": 250.0, "pressure_pa": 76.0 }]
        }"#;
        let dataset =
            FrozenReferenceDataset::<WaterIceSublimationPressureReference>::from_json_strs(
                metadata,
                rows,
                "test.json",
            )
            .unwrap();
        let root = IceSublimationRoot {
            pressure_pa: 0.0,
            residual: 0.0,
            iterations: 0,
        };
        let error = IceSublimationComparisonReport::new(
            &dataset,
            vec![IceSublimationComparisonRow {
                temperature_k: 250.0,
                iapws_pressure_pa: 76.0,
                independent: root,
                ice_from_gas: root,
                gas_from_ice: root,
            }],
        )
        .unwrap_err();
        assert!(
            error.contains("invalid independent root evidence"),
            "{error}"
        );
    }
}
