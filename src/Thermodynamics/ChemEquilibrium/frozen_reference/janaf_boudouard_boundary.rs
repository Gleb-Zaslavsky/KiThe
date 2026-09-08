//! Typed I5 pressure-boundary evidence for the Boudouard reaction.
//!
//! This module extends the frozen JANAF thermochemistry layer without changing
//! it: primary CO/CO2/graphite rows remain the sole external source. Here they
//! are converted into an analytical pressure boundary for a prescribed 50/50
//! ideal-gas reference state and compared with two local KiThe root routes.
//! `P,H`, phase-control policy, and root searching live in the companion test
//! module, not in this read-only evidence model.

use std::fmt;

use crate::Thermodynamics::ChemEquilibrium::frozen_reference::JanafBoudouardReference;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_thermochemistry::{
    JANAF_STANDARD_PRESSURE_PA, JanafBoudouardDerived,
};

/// Fixed controlled gas composition for the first Boudouard boundary story.
pub(crate) const BOUDOUARD_BOUNDARY_Y_CO: f64 = 0.5;
/// Fixed controlled gas composition for the first Boudouard boundary story.
pub(crate) const BOUDOUARD_BOUNDARY_Y_CO2: f64 = 0.5;

/// One external analytical Boudouard boundary derived from a frozen primary
/// JANAF species row, never from a KiThe solver result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct JanafBoudouardBoundaryReference {
    /// Temperature of the derived boundary point, K.
    pub(crate) temperature_k: f64,
    /// Prescribed ideal-gas CO mole fraction (0.5).
    pub(crate) gas_y_co: f64,
    /// Prescribed ideal-gas CO2 mole fraction (0.5).
    pub(crate) gas_y_co2: f64,
    /// Base-10 log of the JANAF-derived equilibrium constant.
    pub(crate) janaf_log10_kp: f64,
    /// JANAF-derived equilibrium constant `10^(log10_kp)`.
    pub(crate) janaf_kp: f64,
    /// Analytical boundary pressure `(y_CO2/y_CO^2) * p0 / Kp`, Pa.
    pub(crate) analytical_boundary_pressure_pa: f64,
}

impl JanafBoudouardBoundaryReference {
    /// Derives `P = (y_CO2 / y_CO^2) * p0 / Kp` for
    /// `2 CO(g) -> CO2(g) + C(gr)`, where pure graphite has unit activity.
    pub(crate) fn from_frozen_reference(
        reference: JanafBoudouardReference,
    ) -> Result<Self, String> {
        let derived = JanafBoudouardDerived::from_reference(reference);
        let janaf_kp = 10_f64.powf(derived.log10_kp_from_gibbs);
        let boundary = (BOUDOUARD_BOUNDARY_Y_CO2 / BOUDOUARD_BOUNDARY_Y_CO.powi(2))
            * JANAF_STANDARD_PRESSURE_PA
            / janaf_kp;
        let result = Self {
            temperature_k: reference.temperature_k,
            gas_y_co: BOUDOUARD_BOUNDARY_Y_CO,
            gas_y_co2: BOUDOUARD_BOUNDARY_Y_CO2,
            janaf_log10_kp: derived.log10_kp_from_gibbs,
            janaf_kp,
            analytical_boundary_pressure_pa: boundary,
        };
        result.validate()?;
        Ok(result)
    }

    /// Validates the prescribed composition and all derived external values.
    pub(crate) fn validate(&self) -> Result<(), String> {
        for (label, value, positive) in [
            ("temperature", self.temperature_k, true),
            ("gas_y_co", self.gas_y_co, true),
            ("gas_y_co2", self.gas_y_co2, true),
            ("janaf_log10_kp", self.janaf_log10_kp, false),
            ("janaf_kp", self.janaf_kp, true),
            (
                "analytical_boundary_pressure_pa",
                self.analytical_boundary_pressure_pa,
                true,
            ),
        ] {
            if !value.is_finite() || (positive && value <= 0.0) {
                return Err(format!(
                    "JANAF Boudouard boundary has invalid {label}={value:e}"
                ));
            }
        }
        if (self.gas_y_co + self.gas_y_co2 - 1.0).abs() > 1e-12 {
            return Err(format!(
                "JANAF Boudouard gas composition must sum to one, got {}",
                self.gas_y_co + self.gas_y_co2
            ));
        }
        Ok(())
    }
}

/// Root evidence emitted by either local pressure-boundary route.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct JanafBoudouardPressureRoot {
    /// Boundary pressure found by this local route, Pa.
    pub(crate) pressure_pa: f64,
    /// Residual magnitude at the converged root.
    pub(crate) residual: f64,
    /// Number of iterations taken by the root solver.
    pub(crate) iterations: usize,
}

/// One three-way Boudouard boundary comparison at a frozen temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct JanafBoudouardBoundaryRow {
    /// Analytical external boundary derived from frozen JANAF data.
    pub(crate) external: JanafBoudouardBoundaryReference,
    /// Root from an independent I1/I2 pressure-boundary equation.
    pub(crate) independent: JanafBoudouardPressureRoot,
    /// Root from the canonical TPD phase-boundary path.
    pub(crate) canonical_tpd: JanafBoudouardPressureRoot,
}

impl JanafBoudouardBoundaryRow {
    /// Relative I1/I2 difference from JANAF's analytical pressure boundary.
    pub(crate) fn independent_external_error(&self) -> f64 {
        relative_delta(
            self.independent.pressure_pa,
            self.external.analytical_boundary_pressure_pa,
        )
    }

    /// Relative canonical-TPD difference from JANAF's analytical boundary.
    pub(crate) fn canonical_external_error(&self) -> f64 {
        relative_delta(
            self.canonical_tpd.pressure_pa,
            self.external.analytical_boundary_pressure_pa,
        )
    }

    /// Relative disagreement of the two internal KiThe routes.
    pub(crate) fn internal_root_error(&self) -> f64 {
        relative_delta(self.canonical_tpd.pressure_pa, self.independent.pressure_pa)
    }
}

/// Characterization-only report for the external JANAF versus local Boudouard
/// pressure-boundary evidence.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct JanafBoudouardBoundaryReport {
    rows: Vec<JanafBoudouardBoundaryRow>,
}

impl JanafBoudouardBoundaryReport {
    /// Builds a report after proving every row is finite, positive, and
    /// strictly ordered by its frozen JANAF temperature identity.
    pub(crate) fn new(rows: Vec<JanafBoudouardBoundaryRow>) -> Result<Self, String> {
        if rows.is_empty() {
            return Err("JANAF Boudouard boundary report requires at least one row".to_owned());
        }
        let mut previous_temperature = None;
        for (index, row) in rows.iter().enumerate() {
            row.external.validate()?;
            if previous_temperature.is_some_and(|previous| row.external.temperature_k <= previous) {
                return Err(format!(
                    "JANAF Boudouard boundary row {index} does not preserve increasing frozen temperatures"
                ));
            }
            previous_temperature = Some(row.external.temperature_k);
            for (route, root) in [
                ("independent", row.independent),
                ("canonical TPD", row.canonical_tpd),
            ] {
                if !root.pressure_pa.is_finite()
                    || root.pressure_pa <= 0.0
                    || !root.residual.is_finite()
                {
                    return Err(format!(
                        "JANAF Boudouard boundary row {index} has invalid {route} root {root:?}"
                    ));
                }
            }
        }
        Ok(Self { rows })
    }

    /// Rows in the requested frozen-temperature order.
    pub(crate) fn rows(&self) -> &[JanafBoudouardBoundaryRow] {
        &self.rows
    }

    /// Largest absolute external relative error across all rows and both routes.
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
    pub(crate) fn rms_independent_external_error(&self) -> f64 {
        root_mean_square(
            self.rows
                .iter()
                .map(JanafBoudouardBoundaryRow::independent_external_error),
        )
    }

    /// Signed mean of the independent-root external errors (bias indicator).
    pub(crate) fn mean_signed_independent_external_error(&self) -> f64 {
        mean(
            self.rows
                .iter()
                .map(JanafBoudouardBoundaryRow::independent_external_error),
        )
    }

    /// Largest I1/I2-versus-TPD root deviation across all rows.
    pub(crate) fn max_internal_root_relative_error(&self) -> f64 {
        self.rows
            .iter()
            .map(|row| row.internal_root_error().abs())
            .fold(0.0_f64, f64::max)
    }

    /// Applies only to two local formulations of the same selected records;
    /// it deliberately does not constrain the external JANAF discrepancy.
    pub(crate) fn validate_internal_roots(&self, max_relative_delta: f64) -> Result<(), String> {
        if !max_relative_delta.is_finite() || max_relative_delta <= 0.0 {
            return Err("Boudouard internal root tolerance must be finite and positive".to_owned());
        }
        for (index, row) in self.rows.iter().enumerate() {
            let delta = row.internal_root_error();
            if delta.abs() > max_relative_delta {
                return Err(format!(
                    "Boudouard row {index} at {} K has I1/I2-vs-TPD root delta {delta:e}, exceeding {max_relative_delta:e}",
                    row.external.temperature_k
                ));
            }
        }
        Ok(())
    }
}

impl fmt::Display for JanafBoudouardBoundaryReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            formatter,
            "{:<7} {:>12} {:>15} {:>15} {:>15} {:>12} {:>12} {:>12}",
            "T [K]",
            "JANAF logK",
            "JANAF P [Pa]",
            "I1/I2 P [Pa]",
            "TPD P [Pa]",
            "I1 err [%]",
            "TPD err [%]",
            "internal [%]",
        )?;
        for row in &self.rows {
            writeln!(
                formatter,
                "{:<7.1} {:>12.6} {:>15.4} {:>15.4} {:>15.4} {:>12.5} {:>12.5} {:>12.5}",
                row.external.temperature_k,
                row.external.janaf_log10_kp,
                row.external.analytical_boundary_pressure_pa,
                row.independent.pressure_pa,
                row.canonical_tpd.pressure_pa,
                100.0 * row.independent_external_error(),
                100.0 * row.canonical_external_error(),
                100.0 * row.internal_root_error(),
            )?;
        }
        write!(
            formatter,
            "summary: max external relative error={:.6e}, rms I1/I2 external={:.6e}, mean signed I1/I2 external={:.6e}, max I1/I2-vs-TPD={:.6e}",
            self.max_external_relative_error(),
            self.rms_independent_external_error(),
            self.mean_signed_independent_external_error(),
            self.max_internal_root_relative_error(),
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

/// Arithmetic mean of an iterable of values.
fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let (sum, count) = values.fold((0.0_f64, 0_usize), |(sum, count), value| {
        (sum + value, count + 1)
    });
    sum / count as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_800_k_boundary_uses_janaf_gibbs_route_and_one_bar_exactly() {
        let boundary =
            JanafBoudouardBoundaryReference::from_frozen_reference(JanafBoudouardReference {
                temperature_k: 800.0,
                co_delta_f_g_kj_mol: -182.497,
                co_log10_kf: 11.916,
                co2_delta_f_g_kj_mol: -395.586,
                co2_log10_kf: 25.829,
            })
            .unwrap();
        assert!((boundary.janaf_log10_kp - 1.997_414).abs() < 1e-6);
        assert!((boundary.analytical_boundary_pressure_pa - 2_012.0).abs() < 20.0);
    }
}
