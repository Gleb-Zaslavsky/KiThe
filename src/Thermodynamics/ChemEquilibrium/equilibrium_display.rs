//! Presentation-only filtering and numeric formatting for equilibrium reports.
//!
//! Solver and validation snapshots always retain every finite physical value.
//! This module operates on copies of those rows solely for terminal, GUI, or
//! export display, so hiding trace species can never change a conservation
//! check, a continuation seed, or an accepted equilibrium result.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_presentation::{
    EquilibriumComponentPresentationRow, EquilibriumPhasePresentationRow,
};

/// Numeric notation used by user-facing equilibrium tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumNumberStyle {
    /// Mantissa/exponent form, useful across many orders of magnitude.
    Scientific { decimals: usize },
    /// Ordinary fixed-point notation for compact values.
    Fixed { decimals: usize },
    /// Engineering notation with an SI prefix attached to the supplied unit.
    EngineeringSi { decimals: usize },
}

/// Unit used when displaying local phase mole fractions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoleFractionDisplay {
    /// Show a dimensionless fraction, for example `2.500e-3`.
    Fraction,
    /// Show a percentage, for example `0.250 %`.
    Percent,
}

/// Display-only threshold and formatting settings.
///
/// `trace_moles` controls row visibility only. The full report remains intact
/// and callers can always choose [`Self::visible_components`] with a different
/// policy later.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumDisplayPolicy {
    trace_moles: f64,
    number_style: EquilibriumNumberStyle,
    mole_fraction_display: MoleFractionDisplay,
}

impl Default for EquilibriumDisplayPolicy {
    fn default() -> Self {
        Self {
            trace_moles: 0.0,
            number_style: EquilibriumNumberStyle::Scientific { decimals: 6 },
            mole_fraction_display: MoleFractionDisplay::Fraction,
        }
    }
}

impl EquilibriumDisplayPolicy {
    /// Builds a validated policy. A zero threshold preserves every row.
    pub fn new(
        trace_moles: f64,
        number_style: EquilibriumNumberStyle,
        mole_fraction_display: MoleFractionDisplay,
    ) -> Result<Self, DisplayPolicyError> {
        if !trace_moles.is_finite() || trace_moles < 0.0 {
            return Err(DisplayPolicyError::InvalidTraceThreshold { trace_moles });
        }
        Ok(Self {
            trace_moles,
            number_style,
            mole_fraction_display,
        })
    }

    /// Trace cutoff in physical mol; this does not alter the source report.
    pub fn trace_moles(&self) -> f64 {
        self.trace_moles
    }

    /// Selected numeric rendering policy.
    pub fn number_style(&self) -> EquilibriumNumberStyle {
        self.number_style
    }

    /// Selected local mole-fraction rendering policy.
    pub fn mole_fraction_display(&self) -> MoleFractionDisplay {
        self.mole_fraction_display
    }

    /// Borrows only rows that are visible under the display threshold.
    ///
    /// A component exactly at the threshold remains visible. `numerical_moles`
    /// is deliberately ignored: trace seeding is a solver diagnostic, whereas
    /// users normally want to filter the published physical state.
    pub fn visible_components<'a>(
        &self,
        components: &'a [EquilibriumComponentPresentationRow],
    ) -> Vec<&'a EquilibriumComponentPresentationRow> {
        components
            .iter()
            .filter(|row| row.physical_moles >= self.trace_moles)
            .collect()
    }

    /// Returns visible amount-column indices for a raw range series.
    ///
    /// The function deliberately returns indices rather than rebuilding a
    /// partial series. A caller can create a compact plot/table while the
    /// original labelled rows stay available for audit, export, and later
    /// threshold changes. This contract is for non-negative physical amount
    /// series such as component moles or phase totals, not residuals.
    pub fn visible_amount_columns(
        &self,
        labels: &[String],
        rows: &[Vec<f64>],
    ) -> Result<Vec<usize>, DisplayPolicyError> {
        if labels.is_empty() || rows.iter().any(|row| row.len() != labels.len()) {
            return Err(DisplayPolicyError::InvalidAmountSeriesShape);
        }
        Ok((0..labels.len())
            .filter(|&column| rows.iter().any(|row| row[column] >= self.trace_moles))
            .collect())
    }

    /// Formats a physical amount without mutating the underlying value.
    pub fn format_moles(&self, moles: f64) -> String {
        self.format_quantity(moles, "mol")
    }

    /// Formats standard molar Gibbs energy without changing its provenance row.
    pub fn format_standard_gibbs(&self, joules_per_mol: f64) -> String {
        self.format_quantity(joules_per_mol, "J/mol")
    }

    /// Formats a local phase mole fraction according to the requested unit.
    pub fn format_mole_fraction(&self, fraction: f64) -> String {
        match self.mole_fraction_display {
            MoleFractionDisplay::Fraction => self.format_dimensionless(fraction, "-"),
            MoleFractionDisplay::Percent => self.format_dimensionless(fraction * 100.0, "%"),
        }
    }

    /// Formats a generic value plus its unit using the selected number style.
    pub fn format_quantity(&self, value: f64, unit: &str) -> String {
        format_quantity(value, unit, self.number_style)
    }

    fn format_dimensionless(&self, value: f64, unit: &str) -> String {
        let style = match self.number_style {
            EquilibriumNumberStyle::EngineeringSi { decimals } => {
                EquilibriumNumberStyle::Fixed { decimals }
            }
            style => style,
        };
        format_quantity(value, unit, style)
    }

    /// Projects phase rows into formatted display rows while keeping the source
    /// report's physical and numerical totals untouched.
    pub fn format_phases(
        &self,
        phases: &[EquilibriumPhasePresentationRow],
    ) -> Vec<EquilibriumFormattedPhaseRow> {
        phases
            .iter()
            .map(|row| EquilibriumFormattedPhaseRow {
                phase: row.phase.clone(),
                status: row.status.clone(),
                physical_total_moles: self.format_moles(row.physical_total_moles),
                numerical_total_moles: self.format_moles(row.numerical_total_moles),
            })
            .collect()
    }

    /// Projects visible component rows into strings suitable for a table cell.
    pub fn format_visible_components(
        &self,
        components: &[EquilibriumComponentPresentationRow],
    ) -> Vec<EquilibriumFormattedComponentRow> {
        self.visible_components(components)
            .into_iter()
            .map(|row| EquilibriumFormattedComponentRow {
                component: row.component.clone(),
                phase: row.phase.clone(),
                substance: row.substance.clone(),
                physical_moles: self.format_moles(row.physical_moles),
                numerical_moles: self.format_moles(row.numerical_moles),
                mole_fraction: self.format_mole_fraction(row.mole_fraction),
                initial_moles: self.format_moles(row.initial_moles),
                standard_gibbs: self.format_standard_gibbs(row.standard_gibbs_j_per_mol),
                library: row.library.clone(),
                record_key: row.record_key.clone(),
                lookup_priority: row.lookup_priority.clone(),
            })
            .collect()
    }
}

/// Formatted phase row for CLI/GUI/export surfaces.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumFormattedPhaseRow {
    pub phase: String,
    pub status: String,
    pub physical_total_moles: String,
    pub numerical_total_moles: String,
}

/// Formatted component row for CLI/GUI/export surfaces.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumFormattedComponentRow {
    pub component: String,
    pub phase: String,
    pub substance: String,
    pub physical_moles: String,
    pub numerical_moles: String,
    pub mole_fraction: String,
    pub initial_moles: String,
    pub standard_gibbs: String,
    pub library: String,
    pub record_key: String,
    pub lookup_priority: String,
}

/// Rejected display-policy input.
#[derive(Debug, Clone, PartialEq)]
pub enum DisplayPolicyError {
    /// Trace filtering only accepts a finite non-negative physical amount.
    InvalidTraceThreshold { trace_moles: f64 },
    /// A presentation-only range filter received rows not aligned to labels.
    InvalidAmountSeriesShape,
}

impl std::fmt::Display for DisplayPolicyError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidTraceThreshold { trace_moles } => write!(
                formatter,
                "trace display threshold must be finite and non-negative, got {trace_moles}"
            ),
            Self::InvalidAmountSeriesShape => write!(
                formatter,
                "amount-series labels and row widths must be non-empty and aligned"
            ),
        }
    }
}

impl std::error::Error for DisplayPolicyError {}

fn format_quantity(value: f64, unit: &str, style: EquilibriumNumberStyle) -> String {
    if value.is_nan() {
        return append_unit("NaN".to_string(), unit);
    }
    if value == f64::INFINITY {
        return append_unit("+inf".to_string(), unit);
    }
    if value == f64::NEG_INFINITY {
        return append_unit("-inf".to_string(), unit);
    }
    match style {
        EquilibriumNumberStyle::Scientific { decimals } => {
            append_unit(format!("{value:.decimals$e}"), unit)
        }
        EquilibriumNumberStyle::Fixed { decimals } => {
            append_unit(format!("{value:.decimals$}"), unit)
        }
        EquilibriumNumberStyle::EngineeringSi { decimals } => {
            let (scaled, prefix) = engineering_scale(value);
            append_unit(format!("{scaled:.decimals$}"), &format!("{prefix}{unit}"))
        }
    }
}

fn append_unit(value: String, unit: &str) -> String {
    if unit.is_empty() || unit == "-" {
        value
    } else {
        format!("{value} {unit}")
    }
}

fn engineering_scale(value: f64) -> (f64, &'static str) {
    if value == 0.0 {
        return (0.0, "");
    }
    const SCALES: &[(i32, &str)] = &[
        (-12, "p"),
        (-9, "n"),
        (-6, "u"),
        (-3, "m"),
        (0, ""),
        (3, "k"),
        (6, "M"),
        (9, "G"),
        (12, "T"),
    ];
    let exponent = value.abs().log10().floor() as i32;
    let requested = exponent.div_euclid(3) * 3;
    let (scale, prefix) = SCALES
        .iter()
        .copied()
        .min_by_key(|(candidate, _)| (candidate - requested).abs())
        .unwrap_or((0, ""));
    (value / 10_f64.powi(scale), prefix)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn component(name: &str, moles: f64, fraction: f64) -> EquilibriumComponentPresentationRow {
        EquilibriumComponentPresentationRow {
            component: format!("gas::{name}"),
            phase: "gas".to_string(),
            substance: name.to_string(),
            physical_moles: moles,
            numerical_moles: moles.max(1e-30),
            mole_fraction: fraction,
            initial_moles: moles,
            standard_gibbs_j_per_mol: -12_500.0,
            library: "NASA_gas".to_string(),
            record_key: name.to_string(),
            lookup_priority: "NASA_gas".to_string(),
        }
    }

    #[test]
    fn trace_filter_never_mutates_or_drops_the_source_rows() {
        let rows = vec![
            component("major", 1.0, 0.999),
            component("trace", 1e-12, 1e-12),
        ];
        let policy = EquilibriumDisplayPolicy::new(
            1e-9,
            EquilibriumNumberStyle::Scientific { decimals: 3 },
            MoleFractionDisplay::Fraction,
        )
        .unwrap();

        let visible = policy.visible_components(&rows);
        assert_eq!(visible.len(), 1);
        assert_eq!(visible[0].substance, "major");
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[1].physical_moles, 1e-12);
    }

    #[test]
    fn unit_and_fraction_formatting_are_presentation_only() {
        let policy = EquilibriumDisplayPolicy::new(
            0.0,
            EquilibriumNumberStyle::EngineeringSi { decimals: 3 },
            MoleFractionDisplay::Percent,
        )
        .unwrap();
        assert_eq!(policy.format_moles(0.00125), "1.250 mmol");
        assert_eq!(policy.format_mole_fraction(0.0025), "0.250 %");
        assert_eq!(policy.format_standard_gibbs(-12_500.0), "-12.500 kJ/mol");
        assert!(
            EquilibriumDisplayPolicy::new(
                -1.0,
                EquilibriumNumberStyle::Fixed { decimals: 2 },
                MoleFractionDisplay::Fraction,
            )
            .is_err()
        );
    }

    #[test]
    fn range_column_filter_returns_indices_without_rewriting_raw_series() {
        let labels = vec!["major".to_string(), "trace".to_string()];
        let rows = vec![vec![1.0, 1e-12], vec![0.8, 2e-12]];
        let policy = EquilibriumDisplayPolicy::new(
            1e-9,
            EquilibriumNumberStyle::Scientific { decimals: 3 },
            MoleFractionDisplay::Fraction,
        )
        .unwrap();
        assert_eq!(
            policy.visible_amount_columns(&labels, &rows).unwrap(),
            vec![0]
        );
        assert_eq!(rows[0][1], 1e-12);
        assert!(
            policy
                .visible_amount_columns(&labels, &[vec![1.0]])
                .is_err()
        );
    }
}
