//! Column provenance for experimental kinetics columns.
//!
//! The goal is to keep a readable chain of how a column was produced:
//! raw import, binding, and transformation steps such as smoothing or resampling.
//! Each step can also carry a quality report when the operation preserves
//! row alignment and a before/after comparison is meaningful.

use crate::Kinetics::experimental_kinetics::ndarray_statistics::FilterQualityReport;

/// Full provenance chain for one column.
#[derive(Debug, Clone)]
pub struct ColumnProvenance {
    /// Root columns that ultimately feed this column.
    pub root_columns: Vec<String>,
    /// Ordered transformation history for the current column.
    pub steps: Vec<ColumnTransformStep>,
}

/// One provenance step together with optional quality metadata.
#[derive(Debug, Clone)]
pub struct ColumnTransformStep {
    pub input_columns: Vec<String>,
    pub output_column: String,
    pub kind: ColumnTransformKind,
    pub quality: Option<TransformQuality>,
    pub operation_id: Option<usize>,
    pub reversible: bool,
}

/// Comparison target used by quality reports.
#[derive(Debug, Clone)]
pub enum QualityReference {
    PreviousStep,
    RawRoot,
    NamedColumn(String),
}

/// Quality status for a transformation.
#[derive(Debug, Clone)]
pub enum QualityStatus {
    Computed,
    NotApplicable(String),
    Failed(String),
}

/// Filter quality together with the comparison target and status.
#[derive(Debug, Clone)]
pub struct TransformQuality {
    pub compared_to: QualityReference,
    pub report: FilterQualityReport,
    pub status: QualityStatus,
}

/// The kind of transformation applied to a column.
///
/// This is intentionally compact and human-readable so it can be shown in GUI
/// tooltips, side panels, or debug logs without extra decoding.
#[derive(Debug, Clone)]
pub enum ColumnTransformKind {
    Raw {
        source: String,
    },
    Import {
        source: String,
    },
    Binding {
        role: String,
        unit: String,
    },
    RollingMean {
        window: usize,
    },
    Hampel {
        window: usize,
        n_sigma: f64,
        strategy: String,
    },
    SavitzkyGolay {
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
        mode: String,
    },
    Lowess {
        frac: f64,
    },
    SplineResample {
        time_col: String,
        new_time_col: String,
        n_points: usize,
        kind: String,
    },
    LsqSplineResample {
        time_col: String,
        new_time_col: String,
        n_points: usize,
        degree: usize,
        n_internal_knots: usize,
        solver: String,
    },
    Fitting {
        model: String,
        x_col: String,
        y_col: String,
        method: String,
        parameters: Vec<(String, f64)>,
        r2: f64,
        tolerance: f64,
        max_iter: usize,
    },
    Manual {
        operation: String,
        details: Option<String>,
    },
}

impl ColumnProvenance {
    /// Create provenance that starts from raw data.
    pub fn raw(column_name: impl Into<String>) -> Self {
        let column_name = column_name.into();
        Self {
            root_columns: vec![column_name.clone()],
            steps: vec![ColumnTransformStep::new(
                vec![column_name.clone()],
                column_name.clone(),
                ColumnTransformKind::Raw {
                    source: column_name,
                },
                None,
                None,
                false,
            )],
        }
    }

    /// Create provenance for a column imported from an external file.
    pub fn imported(column_name: impl Into<String>, source: impl Into<String>) -> Self {
        let column_name = column_name.into();
        let source = source.into();
        Self {
            root_columns: vec![source.clone()],
            steps: vec![ColumnTransformStep::new(
                vec![source.clone()],
                column_name.clone(),
                ColumnTransformKind::Import { source },
                None,
                None,
                false,
            )],
        }
    }

    /// Create provenance for a column that originated from a manual or external action.
    pub fn manual(
        column_name: impl Into<String>,
        operation: impl Into<String>,
        details: Option<String>,
    ) -> Self {
        let column_name = column_name.into();
        Self {
            root_columns: vec![column_name.clone()],
            steps: vec![ColumnTransformStep::new(
                vec![column_name.clone()],
                column_name.clone(),
                ColumnTransformKind::Manual {
                    operation: operation.into(),
                    details,
                },
                None,
                None,
                true,
            )],
        }
    }

    /// Clone an existing provenance chain for a new column name.
    pub fn inherited(column_name: impl Into<String>, source: &ColumnProvenance) -> Self {
        let column_name = column_name.into();
        let mut cloned = source.clone();
        cloned.rename_output(column_name);
        cloned
    }

    /// Append one new step to the chain.
    pub fn append_step(&mut self, step: ColumnTransformStep) {
        for root in &step.input_columns {
            if !self.root_columns.contains(root) {
                self.root_columns.push(root.clone());
            }
        }
        self.steps.push(step);
    }

    /// Append a new step using the common transformation parameters.
    pub fn push_step(
        &mut self,
        input_columns: Vec<String>,
        output_column: impl Into<String>,
        kind: ColumnTransformKind,
        quality: Option<TransformQuality>,
        operation_id: Option<usize>,
        reversible: bool,
    ) {
        self.append_step(ColumnTransformStep::new(
            input_columns,
            output_column,
            kind,
            quality,
            operation_id,
            reversible,
        ));
    }

    /// Update the output column name on all steps.
    pub fn rename_output(&mut self, new_name: impl Into<String>) {
        let new_name = new_name.into();
        for step in &mut self.steps {
            step.output_column = new_name.clone();
        }
    }

    /// Human-readable history for GUI panels and logs.
    pub fn feed_lines(&self) -> Vec<String> {
        let mut lines = vec![format!("roots: {}", self.root_columns.join(", "))];
        lines.extend(
            self.steps
                .iter()
                .enumerate()
                .map(|(idx, step)| format!("{:02}: {}", idx, step.feed_line())),
        );
        lines
    }

    /// Human-readable history as a multiline block.
    pub fn feed_text(&self) -> String {
        if self.steps.is_empty() {
            format!(
                "Column provenance is empty for roots: {}",
                self.root_columns.join(", ")
            )
        } else {
            self.feed_lines().join("\n")
        }
    }
}

impl ColumnTransformStep {
    pub fn new(
        input_columns: Vec<String>,
        output_column: impl Into<String>,
        kind: ColumnTransformKind,
        quality: Option<TransformQuality>,
        operation_id: Option<usize>,
        reversible: bool,
    ) -> Self {
        Self {
            input_columns,
            output_column: output_column.into(),
            kind,
            quality,
            operation_id,
            reversible,
        }
    }

    pub fn feed_line(&self) -> String {
        let quality = self
            .quality
            .as_ref()
            .map(|q| format!("\n{}", q.pretty_lines().join("\n")))
            .unwrap_or_default();

        format!(
            "{} -> {} | input: [{}] | reversible: {} | op_id: {:?}{}",
            self.kind,
            self.output_column,
            self.input_columns.join(", "),
            self.reversible,
            self.operation_id,
            quality
        )
    }
}

impl TransformQuality {
    pub fn computed(report: FilterQualityReport, compared_to: QualityReference) -> Self {
        Self {
            compared_to,
            report,
            status: QualityStatus::Computed,
        }
    }

    pub fn pretty_lines(&self) -> Vec<String> {
        let direction = match &self.compared_to {
            QualityReference::PreviousStep => "previous step".to_string(),
            QualityReference::RawRoot => "raw root".to_string(),
            QualityReference::NamedColumn(name) => format!("column {}", name),
        };
        let status = self.status.to_string();
        vec![
            format!(
                "    quality: comparing {} -> {} ({})",
                self.report.column_raw, self.report.column_filtered, direction
            ),
            format!("      status: {}", status),
            format!("      RMSE: {:.6} (lower is better)", self.report.rmse),
            format!(
                "      NRMSE: {:.6} (lower is better)",
                self.report.normalized_rmse
            ),
            format!(
                "      Corr: {:.6} (higher is better)",
                self.report.correlation
            ),
            format!(
                "      Roughness 1st: {:.6} (lower is better)",
                self.report.roughness_ratio_1st
            ),
            format!(
                "      Roughness 2nd: {:.6} (lower is better)",
                self.report.roughness_ratio_2nd
            ),
        ]
    }
}

impl std::fmt::Display for QualityReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QualityReference::PreviousStep => write!(f, "previous_step"),
            QualityReference::RawRoot => write!(f, "raw_root"),
            QualityReference::NamedColumn(name) => write!(f, "column({})", name),
        }
    }
}

impl std::fmt::Display for QualityStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QualityStatus::Computed => write!(f, "computed"),
            QualityStatus::NotApplicable(reason) => write!(f, "not_applicable({})", reason),
            QualityStatus::Failed(reason) => write!(f, "failed({})", reason),
        }
    }
}

impl std::fmt::Display for TransformQuality {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "comparing {} vs {} [{}] | RMSE={:.6} | NRMSE={:.6} | Corr={:.6} | Roughness 1st={:.6} | Roughness 2nd={:.6}",
            self.report.column_raw,
            self.report.column_filtered,
            self.status,
            self.report.rmse,
            self.report.normalized_rmse,
            self.report.correlation,
            self.report.roughness_ratio_1st,
            self.report.roughness_ratio_2nd
        )
    }
}

impl std::fmt::Display for ColumnTransformKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ColumnTransformKind::Raw { source } => write!(f, "raw({})", source),
            ColumnTransformKind::Import { source } => write!(f, "import({})", source),
            ColumnTransformKind::Binding { role, unit } => {
                write!(f, "binding(role={}, unit={})", role, unit)
            }
            ColumnTransformKind::RollingMean { window } => {
                write!(f, "rolling_mean(window={})", window)
            }
            ColumnTransformKind::Hampel {
                window,
                n_sigma,
                strategy,
            } => write!(
                f,
                "hampel(window={}, n_sigma={}, strategy={})",
                window, n_sigma, strategy
            ),
            ColumnTransformKind::SavitzkyGolay {
                window,
                poly_order,
                deriv,
                delta,
                mode,
            } => write!(
                f,
                "savitzky_golay(window={}, poly_order={}, deriv={}, delta={}, mode={})",
                window, poly_order, deriv, delta, mode
            ),
            ColumnTransformKind::Lowess { frac } => write!(f, "lowess(frac={})", frac),
            ColumnTransformKind::SplineResample {
                time_col,
                new_time_col,
                n_points,
                kind,
            } => write!(
                f,
                "spline_resample(time={}, new_time={}, n_points={}, kind={})",
                time_col, new_time_col, n_points, kind
            ),
            ColumnTransformKind::LsqSplineResample {
                time_col,
                new_time_col,
                n_points,
                degree,
                n_internal_knots,
                solver,
            } => write!(
                f,
                "lsq_spline_resample(time={}, new_time={}, n_points={}, degree={}, knots={}, solver={})",
                time_col, new_time_col, n_points, degree, n_internal_knots, solver
            ),
            ColumnTransformKind::Fitting {
                model,
                x_col,
                y_col,
                method,
                parameters,
                r2,
                tolerance,
                max_iter,
            } => {
                let params = parameters
                    .iter()
                    .map(|(name, value)| format!("{}={:.6}", name, value))
                    .collect::<Vec<_>>()
                    .join(", ");
                write!(
                    f,
                    "fitting(model={}, x={}, y={}, method={}, params=[{}], R2={:.6}, tol={}, max_iter={})",
                    model, x_col, y_col, method, params, r2, tolerance, max_iter
                )
            }
            ColumnTransformKind::Manual { operation, details } => match details {
                Some(details) => write!(f, "manual({}: {})", operation, details),
                None => write!(f, "manual({})", operation),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::ndarray_statistics::FilterQualityReport;

    #[test]
    fn transform_quality_display_includes_all_report_metrics() {
        let report = FilterQualityReport {
            column_raw: "raw".to_string(),
            column_filtered: "filtered".to_string(),
            roughness_ratio_1st: 0.11,
            roughness_ratio_2nd: 0.22,
            rmse: 0.33,
            normalized_rmse: 0.44,
            correlation: 0.55,
        };
        let quality = TransformQuality::computed(report, QualityReference::PreviousStep);
        let rendered = quality.to_string();

        assert!(rendered.contains("comparing raw vs filtered"));
        assert!(rendered.contains("[computed]"));
        assert!(rendered.contains("Roughness 1st=0.110000"));
        assert!(rendered.contains("Roughness 2nd=0.220000"));
        assert!(rendered.contains("RMSE=0.330000"));
        assert!(rendered.contains("NRMSE=0.440000"));
        assert!(rendered.contains("Corr=0.550000"));
    }

    #[test]
    fn transform_quality_pretty_lines_are_human_readable() {
        let report = FilterQualityReport {
            column_raw: "mass".to_string(),
            column_filtered: "sg_mass".to_string(),
            roughness_ratio_1st: 0.171947,
            roughness_ratio_2nd: 0.122832,
            rmse: 4.442255,
            normalized_rmse: 0.114925,
            correlation: 0.993374,
        };
        let quality = TransformQuality::computed(report, QualityReference::PreviousStep);
        let rendered = quality.pretty_lines().join("\n");

        assert!(rendered.contains("quality: comparing mass -> sg_mass"));
        assert!(rendered.contains("status: computed"));
        assert!(rendered.contains("RMSE: 4.442255"));
        assert!(rendered.contains("NRMSE: 0.114925"));
        assert!(rendered.contains("Corr: 0.993374"));
        assert!(rendered.contains("Roughness 1st: 0.171947"));
        assert!(rendered.contains("Roughness 2nd: 0.122832"));
    }
}
