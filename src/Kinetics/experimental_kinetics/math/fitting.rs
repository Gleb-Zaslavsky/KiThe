use crate::Kinetics::experimental_kinetics::column_provenance::{
    ColumnProvenance, ColumnTransformKind, ColumnTransformStep,
};
use crate::Kinetics::experimental_kinetics::experiment_series_main::{TGAExperiment, TGASeries};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    AffectedColumns, ColumnMeta, ColumnNature, ColumnOrigin, TGADataset, TGADomainError, TGASchema,
    Unit,
};
use RustedSciThe::numerical::optimization::fitted_model::{FittedModel, MissingInterval};
use RustedSciThe::numerical::optimization::universal_fitting::{
    Method, UniversalFitting, UniversalFittingResult,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use polars::prelude::*;
use std::collections::HashMap;
use thiserror::Error;

/// A small container for a named kinetic fitting model.
#[derive(Clone, Debug, PartialEq)]
pub struct FittingModel {
    /// Pre-exponential or rate constant.
    pub k: f64,
    /// Additional model parameters for multi-parameter expressions.
    pub params: HashMap<String, f64>,
}

/// Named kinetic models supported by the convenience API.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FittingModelName {
    /// Autocatalytic model.
    Autocatalytic,
    /// First-order kinetic model.
    FirstOrder,
    /// Second-order kinetic model.
    SecondOrder,
    /// Third-order kinetic model.
    ThirdOrder,
    /// Sestak-Berggren model.
    SestakBerggren,
    /// Johnson-Mehl-Avrami model.
    JohnsonMehlAvrami,
    /// Acceleratory model.
    Acceleratory,
    /// Deceleratory model.
    Deceleratory,
    /// Truncated Sestak-Berggren model.
    SestacBergrenTrunc,
    /// Parameterized Sestak-Berggren model.
    SestacBergrenParam,
    /// Single exponential decay.
    DecExp,
    /// Sum of two exponentials.
    TwoExp,
    /// Sum of three exponentials.
    ThreeExp,
    /// Linear model.
    Linear,
    /// Polynomial model of degree `n`.
    Polynom { n: usize },
}

impl FittingModelName {
    /// Stable human-readable model name for logs, provenance and GUI labels.
    pub fn display_name(&self) -> String {
        match self {
            Self::Autocatalytic => "Autocatalytic".to_string(),
            Self::FirstOrder => "FirstOrder".to_string(),
            Self::SecondOrder => "SecondOrder".to_string(),
            Self::ThirdOrder => "ThirdOrder".to_string(),
            Self::SestakBerggren => "SestakBerggren".to_string(),
            Self::JohnsonMehlAvrami => "JohnsonMehlAvrami".to_string(),
            Self::Acceleratory => "Acceleratory".to_string(),
            Self::Deceleratory => "Deceleratory".to_string(),
            Self::SestacBergrenTrunc => "SestacBergrenTrunc".to_string(),
            Self::SestacBergrenParam => "SestacBergrenParam".to_string(),
            Self::DecExp => "DecExp".to_string(),
            Self::TwoExp => "TwoExp".to_string(),
            Self::ThreeExp => "ThreeExp".to_string(),
            Self::Linear => "Linear".to_string(),
            Self::Polynom { n } => format!("Polynom{}", n),
        }
    }

    /// Stable lowercase fragment suitable for generated column names.
    pub fn output_name_fragment(&self) -> String {
        match self {
            Self::Polynom { n } => format!("polynom{}", n),
            _ => self.display_name().to_ascii_lowercase(),
        }
    }

    /// Models exposed by the first GUI fitting pass.
    pub fn all_for_gui() -> Vec<Self> {
        vec![
            Self::Linear,
            Self::DecExp,
            Self::TwoExp,
            Self::ThreeExp,
            Self::FirstOrder,
            Self::SecondOrder,
            Self::ThirdOrder,
            Self::Autocatalytic,
            Self::SestacBergrenTrunc,
            Self::SestacBergrenParam,
            Self::SestakBerggren,
            Self::JohnsonMehlAvrami,
            Self::Acceleratory,
            Self::Deceleratory,
            // Polynomial degree is controlled by a dedicated GUI field.
            Self::Polynom { n: 2 },
        ]
    }

    /// Return the symbolic equation for this kinetic model.
    pub fn equation_str(&self) -> String {
        match self {
            Self::Autocatalytic => "k*x*(1-x)".to_string(),
            Self::FirstOrder => "k*(1-x)".to_string(),
            Self::SecondOrder => "k*(1-x)^2".to_string(),
            Self::ThirdOrder => "k*(1-x)^3".to_string(),
            Self::SestakBerggren => "k*x^n*(1-x)^m*(-ln(1-x))^p".to_string(),
            Self::JohnsonMehlAvrami => "m*(1-x)*(-ln(1-x))^(1-1/m)".to_string(),
            Self::Acceleratory => "m*(1-x)^(1-1/m)".to_string(),
            Self::Deceleratory => "(1-x)^m".to_string(),
            Self::SestacBergrenTrunc => "x^m*(1-x)^n".to_string(),
            Self::SestacBergrenParam => "c*x^m*(1-x)^n".to_string(),
            Self::DecExp => "c*exp(-k*x)".to_string(),
            Self::TwoExp => "a*exp(-p*x) + b*exp(-k*x)".to_string(),
            Self::ThreeExp => "a*exp(-p*x) + b*exp(-k*x) + c*exp(-r*x)".to_string(),
            Self::Linear => "k*x + b".to_string(),
            Self::Polynom { n } => Expr::polyval(*n, "x").0.to_string(),
        }
    }

    pub fn equation(&self) -> Expr {
        Expr::parse_expression(&self.equation_str())
    }

    /// Polynomial degree if this model is a polynomial.
    pub fn polynomial_degree(&self) -> Option<usize> {
        match self {
            Self::Polynom { n } => Some(*n),
            _ => None,
        }
    }

    /// Return the ordered parameter names for this model.
    pub fn vec_of_params(&self) -> Vec<String> {
        match self {
            Self::Autocatalytic => vec!["k".to_string()],
            Self::FirstOrder => vec!["k".to_string()],
            Self::SecondOrder => vec!["k".to_string()],
            Self::ThirdOrder => vec!["k".to_string()],
            Self::SestakBerggren => vec![
                "k".to_string(),
                "n".to_string(),
                "m".to_string(),
                "p".to_string(),
            ],
            Self::JohnsonMehlAvrami => vec!["m".to_string()],
            Self::Acceleratory => vec!["m".to_string()],
            Self::Deceleratory => vec!["m".to_string()],
            Self::SestacBergrenTrunc => vec!["m".to_string(), "n".to_string()],
            Self::SestacBergrenParam => {
                vec!["m".to_string(), "n".to_string(), "c".to_string()]
            }
            Self::DecExp => vec!["c".to_string(), "k".to_string()],
            Self::TwoExp => vec![
                "a".to_string(),
                "p".to_string(),
                "b".to_string(),
                "k".to_string(),
            ],
            Self::ThreeExp => vec![
                "a".to_string(),
                "p".to_string(),
                "b".to_string(),
                "k".to_string(),
                "c".to_string(),
                "r".to_string(),
            ],
            Self::Linear => vec!["k".to_string(), "b".to_string()],
            Self::Polynom { n } => Expr::polyval(*n, "x").1,
        }
    }

    /// Return a simple initial guess vector.
    pub fn vec_of_initial_guess(&self) -> Vec<f64> {
        match self {
            Self::Autocatalytic => vec![0.5],
            Self::FirstOrder => vec![0.5],
            Self::SecondOrder => vec![0.5],
            Self::ThirdOrder => vec![0.5],
            Self::SestakBerggren => vec![1.0, 1.0, 1.0, 1.0],
            Self::JohnsonMehlAvrami => vec![1.0],
            Self::Acceleratory => vec![1.0],
            Self::Deceleratory => vec![1.0],
            Self::SestacBergrenTrunc => vec![1.0, 1.0],
            Self::SestacBergrenParam => vec![1.0, 1.0, 1.0],
            Self::DecExp => vec![1.0, 0.5],
            Self::TwoExp => vec![1.0, 0.5, 1.0, 0.25],
            Self::ThreeExp => vec![1.0, 0.5, 1.0, 0.25, 1.0, 0.1],
            Self::Linear => vec![0.5, 0.5],
            Self::Polynom { n } => vec![1.0; Expr::polyval(*n, "x").1.len()],
        }
    }

    /// Return the nonlinear parameter names used by the VarPro backend.
    ///
    /// Only separable exponential models are handled by VarPro here, because
    /// the amplitudes are linear coefficients and the decay rates are the true
    /// nonlinear unknowns.
    pub fn varpro_unknowns(&self) -> Option<Vec<String>> {
        match self {
            Self::DecExp => Some(vec!["k".to_string()]),
            Self::TwoExp => Some(vec!["p".to_string(), "k".to_string()]),
            Self::ThreeExp => Some(vec!["p".to_string(), "k".to_string(), "r".to_string()]),
            _ => None,
        }
    }

    /// Return the separable basis functions used by the VarPro backend.
    ///
    /// The basis order matches the order of the linear coefficients in the
    /// final solution map.
    pub fn varpro_basis_strs(&self) -> Option<Vec<(&'static str, &'static str)>> {
        match self {
            Self::DecExp => Some(vec![("k", "exp(-k*x)")]),
            Self::TwoExp => Some(vec![("p", "exp(-p*x)"), ("k", "exp(-k*x)")]),
            Self::ThreeExp => Some(vec![
                ("p", "exp(-p*x)"),
                ("k", "exp(-k*x)"),
                ("r", "exp(-r*x)"),
            ]),
            _ => None,
        }
    }

    /// Return the linear coefficient names used by the VarPro backend.
    ///
    /// These are the amplitudes that get solved internally by projection and
    /// later renamed back to the public kinetic-model names.
    pub fn varpro_linear_names(&self) -> Option<Vec<String>> {
        match self {
            Self::DecExp => Some(vec!["c".to_string()]),
            Self::TwoExp => Some(vec!["a".to_string(), "b".to_string()]),
            Self::ThreeExp => Some(vec!["a".to_string(), "b".to_string(), "c".to_string()]),
            _ => None,
        }
    }

    /// Return the default nonlinear initial guess for the VarPro backend.
    pub fn varpro_initial_guess(&self) -> Option<Vec<f64>> {
        match self {
            Self::DecExp => Some(vec![0.5]),
            Self::TwoExp => Some(vec![0.5, 0.25]),
            Self::ThreeExp => Some(vec![0.5, 0.25, 0.1]),
            _ => None,
        }
    }

    /// Return the symbolic independent variable name used by these models.
    pub fn arg_name(&self) -> &'static str {
        "x"
    }

    /// Return the preferred numerical backend for this model.
    ///
    /// The separable exponential models benefit from VarPro, while the rest of
    /// the kinetic library stays on the classic symbolic LM path.
    pub fn preferred_method(&self) -> Method {
        match self {
            Self::DecExp | Self::TwoExp | Self::ThreeExp => Method::VARPRO,
            _ => Method::LM,
        }
    }

    /// Return `true` if the model should be treated as a simple LM fit.
    pub fn uses_classic_lm(&self) -> bool {
        matches!(self.preferred_method(), Method::LM)
    }
}

/// Result of a kinetic fit.
#[derive(Clone, Debug, Default)]
pub struct Fit {
    /// Selected model.
    pub kinetic_model: FittingModelName,
    /// Independent variable samples.
    pub x_data: Vec<f64>,
    /// Observed data samples.
    pub y_data: Vec<f64>,
    /// Optional initial guess.
    pub initial_guess: Option<Vec<f64>>,
    /// Target tolerance.
    pub tolerance: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Runtime evaluator retained so the fitted model can be evaluated later.
    pub fitted_model: Option<FittedModel>,
    /// Fitted parameter map.
    pub map_of_solutions: Option<HashMap<String, f64>>,
    /// Coefficient of determination.
    pub r2: f64,
}

/// Request for fitting one pair of columns inside one experiment.
#[derive(Clone, Debug, PartialEq)]
pub struct FitColumnRequest {
    pub experiment_id: String,
    pub x_col: String,
    pub y_col: String,
    pub output_col: Option<String>,
    pub kinetic_model: FittingModelName,
    pub initial_guess: Option<Vec<f64>>,
    pub tolerance: f64,
    pub max_iter: usize,
    /// Optional inclusive X range used to choose observations for fitting.
    pub training_range: Option<(f64, f64)>,
    /// Optional open interval excluded from training, e.g. a missing-data gap.
    pub excluded_range: Option<(f64, f64)>,
}

impl FitColumnRequest {
    pub fn new(
        experiment_id: impl Into<String>,
        x_col: impl Into<String>,
        y_col: impl Into<String>,
        kinetic_model: FittingModelName,
    ) -> Self {
        Self {
            experiment_id: experiment_id.into(),
            x_col: x_col.into(),
            y_col: y_col.into(),
            output_col: None,
            kinetic_model,
            initial_guess: None,
            tolerance: 1e-6,
            max_iter: 300,
            training_range: None,
            excluded_range: None,
        }
    }

    pub fn with_output_col(mut self, output_col: impl Into<String>) -> Self {
        self.output_col = Some(output_col.into());
        self
    }

    pub fn with_initial_guess(mut self, initial_guess: Vec<f64>) -> Self {
        self.initial_guess = Some(initial_guess);
        self
    }

    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iter = max_iter;
        self
    }

    pub fn with_training_range(mut self, start: f64, end: f64) -> Self {
        self.training_range = Some((start, end));
        self
    }

    pub fn excluding_range(mut self, start: f64, end: f64) -> Self {
        self.excluded_range = Some((start, end));
        self
    }
}

/// Result of fitting a dataframe column and publishing fitted y back to the series.
#[derive(Clone, Debug, PartialEq)]
pub struct FitColumnResult {
    pub experiment_id: String,
    pub x_col: String,
    pub y_col: String,
    pub output_col: String,
    pub kinetic_model: FittingModelName,
    pub method: String,
    pub coefficients: Vec<(String, f64)>,
    pub r2: f64,
    pub tolerance: f64,
    pub max_iter: usize,
}

/// A successful, non-mutating fit session that can be evaluated later.
#[derive(Clone, Debug)]
pub struct FitColumnSession {
    pub result: FitColumnResult,
    pub fitted_model: FittedModel,
    pub training_range: (f64, f64),
    pub excluded_range: Option<(f64, f64)>,
    pub left_rmse: Option<f64>,
    pub right_rmse: Option<f64>,
    pub source_history_len: usize,
    x_unit: Unit,
    x_nature: ColumnNature,
    source_x: Vec<f64>,
    source_y: Vec<f64>,
    pub training_x: Vec<f64>,
    pub training_y: Vec<f64>,
    y_unit: Unit,
    y_nature: ColumnNature,
}

/// How a set of fitted columns should be evaluated on one common X grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FitEvaluationMode {
    /// Evaluate a uniform grid wholly outside the fitted model's training domain.
    /// `points` is the number of coordinates, including both range endpoints.
    Extrapolate { start: f64, end: f64, points: usize },
    /// Evaluate exactly this many coordinates strictly inside an unobserved interval.
    FillGap {
        start: f64,
        end: f64,
        missing_points: usize,
    },
}

/// Fully evaluated, complete reconstruction ready for explicit insertion.
#[derive(Clone, Debug)]
pub struct FittedEvaluationBatch {
    pub source_experiment_id: String,
    pub source_history_len: usize,
    pub x_col: String,
    pub x_unit: Unit,
    pub x_nature: ColumnNature,
    pub mode: FitEvaluationMode,
    /// X values generated by the fitted models, separate from retained source rows.
    pub evaluation_x: Vec<f64>,
    pub source_x: Vec<f64>,
    pub x_values: Vec<f64>,
    pub columns: Vec<EvaluatedColumn>,
}

/// Evaluate all fit sessions on one exact shared grid and combine source samples with predictions.
pub fn evaluate_fit_batch(
    sessions: &[FitColumnSession],
    mode: FitEvaluationMode,
) -> Result<FittedEvaluationBatch, TGADomainError> {
    let first = sessions.first().ok_or_else(|| {
        TGADomainError::InvalidOperation("At least one fitted column is required".to_string())
    })?;
    let source_id = &first.result.experiment_id;
    let x_col = &first.result.x_col;
    let training_range = first.training_range;
    let excluded_range = first.excluded_range;
    let source_x = &first.source_x;

    let mut names = std::collections::HashSet::with_capacity(sessions.len());
    for session in sessions {
        if &session.result.experiment_id != source_id
            || &session.result.x_col != x_col
            || session.training_range != training_range
            || session.excluded_range != excluded_range
            || session.source_history_len != first.source_history_len
            || session.source_x != *source_x
        {
            return Err(TGADomainError::InvalidOperation(
                "Fit sessions must share the experiment, X column, source grid, fit range, and excluded interval".to_string(),
            ));
        }
        if !names.insert(session.result.output_col.as_str()) {
            return Err(TGADomainError::InvalidOperation(format!(
                "Duplicate reconstruction output column '{}'",
                session.result.output_col
            )));
        }
    }

    let (evaluation_x, evaluated_values) = match mode {
        FitEvaluationMode::Extrapolate { start, end, points } => {
            let values = first
                .fitted_model
                .extrapolate_uniform_range(start, end, points)
                .map_err(|err| TGADomainError::InvalidOperation(err.to_string()))?;
            let x = values.x;
            let y = sessions
                .iter()
                .map(|session| {
                    let result = session
                        .fitted_model
                        .extrapolate_at(&x)
                        .map_err(|err| TGADomainError::InvalidOperation(err.to_string()))?;
                    Ok(result.values.column(0).iter().copied().collect::<Vec<_>>())
                })
                .collect::<Result<Vec<_>, TGADomainError>>()?;
            (x, y)
        }
        FitEvaluationMode::FillGap {
            start,
            end,
            missing_points,
        } => {
            let interval = MissingInterval {
                start,
                end,
                missing_points,
            };
            let first_values = first
                .fitted_model
                .fill_missing_intervals(&[interval])
                .map_err(|err| TGADomainError::InvalidOperation(err.to_string()))?
                .into_iter()
                .next()
                .ok_or_else(|| {
                    TGADomainError::InvalidOperation(
                        "No reconstruction points were generated".to_string(),
                    )
                })?;
            let x = first_values.x;
            let mut y_values = vec![first_values.values.column(0).iter().copied().collect()];
            for session in sessions.iter().skip(1) {
                let reconstructed = session
                    .fitted_model
                    .fill_missing_intervals(&[interval])
                    .map_err(|err| TGADomainError::InvalidOperation(err.to_string()))?
                    .into_iter()
                    .next()
                    .ok_or_else(|| {
                        TGADomainError::InvalidOperation(
                            "No reconstruction points were generated".to_string(),
                        )
                    })?;
                if reconstructed.x != x {
                    return Err(TGADomainError::InvalidOperation(
                        "Fitted models produced different missing-point grids".to_string(),
                    ));
                }
                y_values.push(reconstructed.values.column(0).iter().copied().collect());
            }
            (x, y_values)
        }
    };

    if evaluation_x.iter().any(|value| !value.is_finite())
        || evaluation_x.windows(2).any(|pair| pair[1] <= pair[0])
    {
        return Err(TGADomainError::InvalidOperation(
            "Evaluation grid must be finite and strictly increasing".to_string(),
        ));
    }
    if source_x.iter().any(|source| {
        evaluation_x
            .binary_search_by(|value| value.total_cmp(source))
            .is_ok()
    }) {
        return Err(TGADomainError::InvalidOperation(
            "Evaluation grid overlaps an existing X value".to_string(),
        ));
    }

    let mut indexed_x = source_x
        .iter()
        .enumerate()
        .map(|(index, &x)| (x, Some(index), None))
        .chain(
            evaluation_x
                .iter()
                .enumerate()
                .map(|(index, &x)| (x, None, Some(index))),
        )
        .collect::<Vec<_>>();
    indexed_x.sort_by(|a, b| a.0.total_cmp(&b.0));
    if indexed_x.windows(2).any(|pair| pair[0].0 >= pair[1].0) {
        return Err(TGADomainError::InvalidOperation(
            "Combined source and evaluation X values must be unique".to_string(),
        ));
    }
    let x_values = indexed_x.iter().map(|entry| entry.0).collect();

    let columns = sessions
        .iter()
        .zip(evaluated_values)
        .map(|(session, predicted)| {
            let values = indexed_x
                .iter()
                .map(
                    |(_, source_index, evaluation_index)| match (source_index, evaluation_index) {
                        (Some(index), None) => session.source_y[*index],
                        (None, Some(index)) => predicted[*index],
                        _ => f64::NAN,
                    },
                )
                .collect::<Vec<_>>();
            if values.iter().any(|value| !value.is_finite()) {
                return Err(TGADomainError::InvalidOperation(
                    "Reconstructed output contains non-finite values".to_string(),
                ));
            }
            Ok(EvaluatedColumn {
                name: session.result.output_col.clone(),
                unit: session.y_unit,
                nature: session.y_nature,
                values,
                fit: session.result.clone(),
            })
        })
        .collect::<Result<Vec<_>, TGADomainError>>()?;

    Ok(FittedEvaluationBatch {
        source_experiment_id: source_id.clone(),
        source_history_len: first.source_history_len,
        x_col: x_col.clone(),
        x_unit: first.x_unit,
        x_nature: first.x_nature,
        mode,
        evaluation_x,
        source_x: source_x.clone(),
        x_values,
        columns,
    })
}

#[derive(Clone, Debug)]
pub struct EvaluatedColumn {
    pub name: String,
    pub unit: Unit,
    pub nature: ColumnNature,
    pub values: Vec<f64>,
    pub fit: FitColumnResult,
}

/// Errors returned by the kinetic fitting convenience layer.
#[derive(Debug, Error)]
pub enum KineticFittingError {
    /// The x data is empty.
    #[error("x data cannot be empty")]
    XDataEmpty,
    /// The y data is empty.
    #[error("y data cannot be empty")]
    YDataEmpty,
    /// The initial guess is missing.
    #[error("initial guess cannot be empty")]
    InitialGuessMissing,
    /// The solver backend failed.
    #[error("fit failed: {0}")]
    FitFailed(String),
    /// The solver result was missing a map of solutions.
    #[error("fit completed but did not produce a solution map")]
    MissingSolutionMap,
    /// The solver completed but did not provide a coefficient of determination.
    #[error("fit completed but did not provide R-squared")]
    MissingR2,
    /// A model selected the VarPro backend without providing its VarPro definition.
    #[error("model is missing required VarPro configuration")]
    InvalidModelConfiguration,
}

impl Fit {
    /// Create a new kinetic fitting builder for a selected model.
    pub fn new(kinetic_model: FittingModelName) -> Self {
        Self {
            kinetic_model,
            x_data: Vec::new(),
            y_data: Vec::new(),
            initial_guess: None,
            tolerance: 1e-6,
            max_iter: 300,
            fitted_model: None,
            map_of_solutions: None,
            r2: 0.0,
        }
    }

    /// Set the x and y data together.
    pub fn with_data(mut self, x_data: Vec<f64>, y_data: Vec<f64>) -> Self {
        self.x_data = x_data;
        self.y_data = y_data;
        self
    }

    /// Set the initial guess for the parameters.
    pub fn with_initial_guess(mut self, initial_guess: Vec<f64>) -> Self {
        self.initial_guess = Some(initial_guess);
        self
    }

    /// Set the numerical tolerance.
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Set the maximum number of iterations.
    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iter = max_iter;
        self
    }

    /// Run the selected kinetic fit and store the result in this builder.
    pub fn fit(mut self) -> Result<Self, KineticFittingError> {
        self.validate()?;

        let initial_guess = self.effective_initial_guess();
        let method = self.kinetic_model.preferred_method();

        let universal = match method {
            Method::LM => UniversalFitting::new()
                .with_method(Method::LM)
                .with_data(self.x_data.clone(), self.y_data.clone())
                .with_equation(self.kinetic_model.equation())
                .with_arg(self.kinetic_model.arg_name().to_string())
                .with_unknowns(self.kinetic_model.vec_of_params())
                .with_initial_guess(initial_guess)
                .with_tolerance(self.tolerance)
                .with_max_iterations(self.max_iter)
                .build(),
            Method::GAVIN => UniversalFitting::new()
                .with_gavin()
                .with_data(self.x_data.clone(), self.y_data.clone())
                .with_equation(self.kinetic_model.equation())
                .with_arg(self.kinetic_model.arg_name().to_string())
                .with_unknowns(self.kinetic_model.vec_of_params())
                .with_initial_guess(initial_guess)
                .with_tolerance(self.tolerance)
                .with_max_iterations(self.max_iter)
                .build(),
            Method::VARPRO => {
                let nonlinear_unknowns = self
                    .kinetic_model
                    .varpro_unknowns()
                    .ok_or(KineticFittingError::InvalidModelConfiguration)?;
                let basis_functions = self
                    .kinetic_model
                    .varpro_basis_strs()
                    .ok_or(KineticFittingError::InvalidModelConfiguration)?;
                let mut universal = UniversalFitting::new()
                    .with_method(Method::VARPRO)
                    .with_data(self.x_data.clone(), self.y_data.clone())
                    .with_arg(self.kinetic_model.arg_name().to_string())
                    .with_parameters(nonlinear_unknowns)
                    .with_initial_guess(initial_guess)
                    .with_tolerance(self.tolerance)
                    .with_max_iterations(self.max_iter);

                for (parameter_name, basis) in basis_functions {
                    universal = universal.with_basis_str(parameter_name, basis);
                }

                universal.build()
            }
        }
        .map_err(|err| KineticFittingError::FitFailed(err.to_string()))?;

        self.fitted_model = Some(universal.fitted_model().clone());

        match universal {
            UniversalFittingResult::LM(fit) => {
                self.map_of_solutions = fit.solution_map();
                self.r2 = fit.r_squared().ok_or(KineticFittingError::MissingR2)?;
                if self.map_of_solutions.is_none() {
                    return Err(KineticFittingError::MissingSolutionMap);
                }
                Ok(self)
            }
            UniversalFittingResult::VARPRO(fit) => {
                let map = fit
                    .solution_map()
                    .ok_or(KineticFittingError::MissingSolutionMap)?;
                self.map_of_solutions = Some(self.remap_varpro_solution(map));
                self.r2 = fit.r_squared().ok_or(KineticFittingError::MissingR2)?;
                Ok(self)
            }
            UniversalFittingResult::GAVIN(fit) => {
                self.map_of_solutions = Some(
                    fit.parameter_names
                        .clone()
                        .into_iter()
                        .zip(fit.fit.parameters.iter().copied())
                        .collect(),
                );
                self.r2 = fit.fit.r_squared.ok_or(KineticFittingError::MissingR2)?;
                Ok(self)
            }
        }
    }

    /// Store the fitted coefficients as a name-to-value map.
    pub fn with_solution_map(mut self, map: HashMap<String, f64>) -> Self {
        self.map_of_solutions = Some(map);
        self
    }

    /// Return the fitted coefficients as a name-to-value map.
    pub fn solution_map(&self) -> Option<HashMap<String, f64>> {
        self.map_of_solutions.clone()
    }

    /// Short alias for [`solution_map`](Self::solution_map).
    pub fn get_map_of_solutions(&self) -> Option<HashMap<String, f64>> {
        self.solution_map()
    }

    /// Store the coefficient of determination of the last fit.
    pub fn with_r_squared(mut self, r2: f64) -> Self {
        self.r2 = r2;
        self
    }

    /// Return the coefficient of determination of the last fit.
    pub fn r_squared(&self) -> f64 {
        self.r2
    }

    /// Short alias for [`r_squared`](Self::r_squared).
    pub fn r2(&self) -> f64 {
        self.r_squared()
    }

    /// Evaluate the fitted model at the stored x values.
    pub fn evaluate_solution(&self) -> Result<Vec<f64>, KineticFittingError> {
        let model = self
            .fitted_model
            .as_ref()
            .ok_or(KineticFittingError::MissingSolutionMap)?;
        let values = model
            .evaluate_at(&self.x_data)
            .map_err(|err| KineticFittingError::FitFailed(err.to_string()))?;
        Ok(values.values.column(0).iter().copied().collect())
    }

    /// Evaluate the retained fit at arbitrary finite coordinates.
    pub fn evaluate_at(&self, x: &[f64]) -> Result<Vec<f64>, KineticFittingError> {
        let model = self
            .fitted_model
            .as_ref()
            .ok_or(KineticFittingError::MissingSolutionMap)?;
        let values = model
            .evaluate_at(x)
            .map_err(|err| KineticFittingError::FitFailed(err.to_string()))?;
        Ok(values.values.column(0).iter().copied().collect())
    }

    fn validate(&self) -> Result<(), KineticFittingError> {
        if self.x_data.is_empty() {
            return Err(KineticFittingError::XDataEmpty);
        }
        if self.y_data.is_empty() {
            return Err(KineticFittingError::YDataEmpty);
        }
        let expected_guess_len = match self.kinetic_model.preferred_method() {
            Method::LM => self.kinetic_model.vec_of_params().len(),
            Method::GAVIN => self.kinetic_model.vec_of_params().len(),
            Method::VARPRO => self
                .kinetic_model
                .varpro_unknowns()
                .map(|unknowns| unknowns.len())
                .unwrap_or_else(|| self.kinetic_model.vec_of_params().len()),
        };
        let guess_len = self.effective_initial_guess().len();
        if guess_len != expected_guess_len {
            return Err(KineticFittingError::FitFailed(
                "initial guess length does not match parameter count".to_string(),
            ));
        }
        Ok(())
    }

    /// Return the actual initial guess used by the solver.
    ///
    /// If the caller does not provide an explicit guess, the model-specific
    /// default vector is used instead.
    fn effective_initial_guess(&self) -> Vec<f64> {
        match self.kinetic_model.preferred_method() {
            Method::LM => self
                .initial_guess
                .clone()
                .unwrap_or_else(|| self.kinetic_model.vec_of_initial_guess()),
            Method::GAVIN => self
                .initial_guess
                .clone()
                .unwrap_or_else(|| self.kinetic_model.vec_of_initial_guess()),
            Method::VARPRO => self
                .initial_guess
                .clone()
                .map(|guess| self.varpro_initial_guess_from_any_guess(&guess))
                .or_else(|| self.kinetic_model.varpro_initial_guess())
                .unwrap_or_else(|| self.kinetic_model.vec_of_initial_guess()),
        }
    }

    /// Extract the nonlinear seed values from any accepted VarPro guess shape.
    fn varpro_initial_guess_from_any_guess(&self, guess: &[f64]) -> Vec<f64> {
        if let Some(varpro_unknowns) = self.kinetic_model.varpro_unknowns() {
            if guess.len() == varpro_unknowns.len() {
                return guess.to_vec();
            }
        }

        match self.kinetic_model {
            FittingModelName::DecExp if guess.len() >= 2 => vec![guess[1]],
            FittingModelName::TwoExp if guess.len() >= 4 => vec![guess[1], guess[3]],
            FittingModelName::ThreeExp if guess.len() >= 6 => vec![guess[1], guess[3], guess[5]],
            _ => guess.to_vec(),
        }
    }

    /// Rename VarPro coefficient keys so the public API stays model-centric.
    fn remap_varpro_solution(&self, mut map: HashMap<String, f64>) -> HashMap<String, f64> {
        if let Some(linear_names) = self.kinetic_model.varpro_linear_names() {
            for (idx, name) in linear_names.into_iter().enumerate() {
                if let Some(value) = map.remove(&format!("c{idx}")) {
                    map.insert(name, value);
                }
            }
        }
        map
    }
}

impl TGASeries {
    /// Fit one Y column without changing the source experiment.
    ///
    /// If an interval is excluded, the one model is fitted jointly to every
    /// remaining point in the training window, including both sides of the gap.
    pub fn fit_column_for_evaluation(
        &self,
        request: FitColumnRequest,
    ) -> Result<FitColumnSession, TGADomainError> {
        let output_col = request
            .output_col
            .clone()
            .unwrap_or_else(|| default_fit_output_name(request.kinetic_model, &request.y_col));
        if output_col.trim().is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "Fit output name cannot be empty".to_string(),
            ));
        }
        let exp = self.get_experiment_by_id(&request.experiment_id)?;
        let source_history_len = exp.dataset.history_of_operations.vector_of_operations.len();
        let x_meta = exp
            .dataset
            .schema
            .columns
            .get(&request.x_col)
            .ok_or_else(|| TGADomainError::ColumnNotFound(request.x_col.clone()))?;
        let y_meta = exp
            .dataset
            .schema
            .columns
            .get(&request.y_col)
            .ok_or_else(|| TGADomainError::ColumnNotFound(request.y_col.clone()))?;
        if request.x_col == output_col {
            return Err(TGADomainError::InvalidOperation(
                "Reconstruction output cannot overwrite the X column".to_string(),
            ));
        }
        let (source_x, source_y, _) = materialize_fit_columns(exp, &request.x_col, &request.y_col)?;
        if source_x.windows(2).any(|pair| pair[1] <= pair[0]) {
            return Err(TGADomainError::InvalidOperation(format!(
                "X column '{}' must be strictly increasing before fitting",
                request.x_col
            )));
        }
        let source_min = source_x
            .first()
            .copied()
            .ok_or_else(|| TGADomainError::EmptyColumn(request.x_col.clone()))?;
        let source_max = source_x
            .last()
            .copied()
            .ok_or_else(|| TGADomainError::EmptyColumn(request.x_col.clone()))?;
        let training_range = request.training_range.unwrap_or((source_min, source_max));
        validate_finite_increasing_range(training_range, "training")?;
        if training_range.0 < source_min || training_range.1 > source_max {
            return Err(TGADomainError::InvalidOperation(
                "Training range must be inside the source X domain".to_string(),
            ));
        }
        if let Some((gap_start, gap_end)) = request.excluded_range {
            validate_finite_increasing_range((gap_start, gap_end), "excluded")?;
            if gap_start < training_range.0 || gap_end > training_range.1 {
                return Err(TGADomainError::InvalidOperation(
                    "Excluded interval must be inside the training range".to_string(),
                ));
            }
            if source_x.iter().any(|x| *x > gap_start && *x < gap_end) {
                return Err(TGADomainError::InvalidOperation(
                    "The excluded interval contains recorded observations; choose an actually missing interval".to_string(),
                ));
            }
        }

        let mut training_x = Vec::new();
        let mut training_y = Vec::new();
        let mut left_x = Vec::new();
        let mut left_y = Vec::new();
        let mut right_x = Vec::new();
        let mut right_y = Vec::new();
        for (&x, &y) in source_x.iter().zip(&source_y) {
            if x < training_range.0 || x > training_range.1 {
                continue;
            }
            if let Some((gap_start, gap_end)) = request.excluded_range {
                if x > gap_start && x < gap_end {
                    continue;
                }
                if x <= gap_start {
                    left_x.push(x);
                    left_y.push(y);
                }
                if x >= gap_end {
                    right_x.push(x);
                    right_y.push(y);
                }
            }
            training_x.push(x);
            training_y.push(y);
        }
        if training_x.len() < request.kinetic_model.vec_of_params().len().max(2) {
            return Err(TGADomainError::InvalidOperation(
                "The selected fitting range has too few observations".to_string(),
            ));
        }
        if request.excluded_range.is_some() && (left_x.is_empty() || right_x.is_empty()) {
            return Err(TGADomainError::InvalidOperation(
                "Gap reconstruction requires observations on both sides of the missing interval"
                    .to_string(),
            ));
        }

        let mut fit = Fit::new(request.kinetic_model)
            .with_data(training_x.clone(), training_y.clone())
            .with_tolerance(request.tolerance)
            .with_max_iterations(request.max_iter);
        if let Some(guess) = request.initial_guess.clone() {
            fit = fit.with_initial_guess(guess);
        }
        let fit = fit.fit().map_err(fitting_error_to_domain)?;
        let result = FitColumnResult {
            experiment_id: request.experiment_id,
            x_col: request.x_col,
            y_col: request.y_col,
            output_col,
            kinetic_model: request.kinetic_model,
            method: format!("{:?}", request.kinetic_model.preferred_method()),
            coefficients: ordered_coefficients(request.kinetic_model, &fit),
            r2: fit.r_squared(),
            tolerance: request.tolerance,
            max_iter: request.max_iter,
        };
        let model = fit.fitted_model.clone().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "Fitting completed without a reusable evaluation model".to_string(),
            )
        })?;
        Ok(FitColumnSession {
            result,
            fitted_model: model,
            training_range,
            excluded_range: request.excluded_range,
            left_rmse: request
                .excluded_range
                .map(|_| model_rmse(&fit, &left_x, &left_y))
                .transpose()?,
            right_rmse: request
                .excluded_range
                .map(|_| model_rmse(&fit, &right_x, &right_y))
                .transpose()?,
            source_history_len,
            x_unit: x_meta.unit,
            x_nature: x_meta.nature,
            source_x,
            source_y,
            training_x,
            training_y,
            y_unit: y_meta.unit,
            y_nature: y_meta.nature,
        })
    }

    /// Insert a fully evaluated batch as a new experiment, leaving its source untouched.
    pub fn push_fitted_evaluation_batch(
        &mut self,
        batch: &FittedEvaluationBatch,
        new_id: impl Into<String>,
    ) -> Result<(), TGADomainError> {
        let new_id = new_id.into();
        if new_id.trim().is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "Experiment ID cannot be empty".to_string(),
            ));
        }
        if self.get_experiment_by_id(&new_id).is_ok() {
            return Err(TGADomainError::DuplicateExperimentId(new_id));
        }
        let source = self
            .get_experiment_by_id(&batch.source_experiment_id)?
            .clone();
        if source
            .dataset
            .history_of_operations
            .vector_of_operations
            .len()
            != batch.source_history_len
        {
            return Err(TGADomainError::InvalidOperation(
                "Source experiment changed after fitting; refit before committing".to_string(),
            ));
        }
        if batch.x_values.is_empty()
            || batch.x_values.iter().any(|value| !value.is_finite())
            || batch.x_values.windows(2).any(|pair| pair[1] <= pair[0])
            || batch.x_col.trim().is_empty()
            || batch.columns.is_empty()
            || batch.columns.iter().any(|column| {
                column.values.len() != batch.x_values.len()
                    || column.values.iter().any(|value| !value.is_finite())
                    || !column.fit.r2.is_finite()
                    || column
                        .fit
                        .coefficients
                        .iter()
                        .any(|(_, value)| !value.is_finite())
            })
        {
            return Err(TGADomainError::InvalidOperation(
                "Reconstruction batch must contain finite, aligned columns on a strictly increasing grid".to_string(),
            ));
        }
        validate_finite_increasing_range(
            match batch.mode {
                FitEvaluationMode::Extrapolate { start, end, .. }
                | FitEvaluationMode::FillGap { start, end, .. } => (start, end),
            },
            "evaluation",
        )?;
        if batch.evaluation_x.is_empty()
            || batch.evaluation_x.iter().any(|value| !value.is_finite())
            || batch.evaluation_x.windows(2).any(|pair| pair[1] <= pair[0])
        {
            return Err(TGADomainError::InvalidOperation(
                "Evaluation X values must be finite, nonempty, and strictly increasing".to_string(),
            ));
        }
        match batch.mode {
            FitEvaluationMode::Extrapolate { start, end, points } => {
                if points != batch.evaluation_x.len()
                    || batch.evaluation_x.iter().any(|x| *x < start || *x > end)
                {
                    return Err(TGADomainError::InvalidOperation(
                        "Extrapolation grid does not match its requested range and point count"
                            .to_string(),
                    ));
                }
            }
            FitEvaluationMode::FillGap {
                start,
                end,
                missing_points,
            } => {
                if missing_points != batch.evaluation_x.len()
                    || batch.evaluation_x.iter().any(|x| *x <= start || *x >= end)
                {
                    return Err(TGADomainError::InvalidOperation(
                        "Gap grid does not match its requested interval and missing-point count"
                            .to_string(),
                    ));
                }
            }
        }
        let source_y_columns = batch
            .columns
            .iter()
            .map(|column| column.fit.y_col.clone())
            .collect::<Vec<_>>();
        let x_meta = source
            .dataset
            .schema
            .columns
            .get(&batch.x_col)
            .ok_or_else(|| TGADomainError::ColumnNotFound(batch.x_col.clone()))?;
        if x_meta.unit != batch.x_unit || x_meta.nature != batch.x_nature {
            return Err(TGADomainError::InvalidOperation(
                "X column metadata changed after fitting".to_string(),
            ));
        }
        let (current_x, source_values) =
            materialize_fit_columns_many(&source, &batch.x_col, &source_y_columns)?;
        if current_x != batch.source_x
            || batch.source_x.iter().any(|x| {
                batch
                    .evaluation_x
                    .binary_search_by(|value| value.total_cmp(x))
                    .is_ok()
            })
        {
            return Err(TGADomainError::InvalidOperation(
                "Reconstruction source grid is stale or overlaps evaluation points".to_string(),
            ));
        }
        let mut expected_grid = batch
            .source_x
            .iter()
            .copied()
            .chain(batch.evaluation_x.iter().copied())
            .collect::<Vec<_>>();
        expected_grid.sort_by(f64::total_cmp);
        if expected_grid != batch.x_values {
            return Err(TGADomainError::InvalidOperation(
                "Reconstruction X values do not match the source and evaluation grids".to_string(),
            ));
        }
        for column in &batch.columns {
            if column.fit.experiment_id != batch.source_experiment_id
                || column.fit.x_col != batch.x_col
                || column.fit.output_col != column.name
            {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Fitted metadata for output '{}' does not match the reconstruction batch",
                    column.name
                )));
            }
            let y_meta = source
                .dataset
                .schema
                .columns
                .get(&column.fit.y_col)
                .ok_or_else(|| TGADomainError::ColumnNotFound(column.fit.y_col.clone()))?;
            if y_meta.unit != column.unit || y_meta.nature != column.nature {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Source metadata for '{}' changed after fitting",
                    column.fit.y_col
                )));
            }
            let source_values = source_values
                .get(&column.fit.y_col)
                .ok_or_else(|| TGADomainError::ColumnNotFound(column.fit.y_col.clone()))?;
            for (x, source_value) in current_x.iter().zip(source_values) {
                let index = batch
                    .x_values
                    .binary_search_by(|value| value.total_cmp(x))
                    .map_err(|_| {
                        TGADomainError::InvalidOperation(
                            "Source coordinate is absent from reconstructed grid".to_string(),
                        )
                    })?;
                if column.values[index] != *source_value {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Source values for '{}' changed after evaluation",
                        column.fit.y_col
                    )));
                }
            }
        }
        let mut names = std::collections::HashSet::new();
        names.insert(batch.x_col.as_str());
        for column in &batch.columns {
            if column.name.trim().is_empty() || !names.insert(column.name.as_str()) {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Invalid or duplicate reconstruction column name '{}'",
                    column.name
                )));
            }
        }

        let mut df_columns: Vec<Column> =
            vec![Series::new(batch.x_col.clone().into(), batch.x_values.as_slice()).into()];
        for column in &batch.columns {
            df_columns
                .push(Series::new(column.name.clone().into(), column.values.as_slice()).into());
        }
        let df = DataFrame::new(batch.x_values.len(), df_columns)?;
        let mut dataset = TGADataset::default();
        dataset.frame = df.lazy();
        dataset.schema = TGASchema::default();
        let x_provenance = source
            .dataset
            .column_provenance(&batch.x_col)
            .cloned()
            .unwrap_or_else(|| ColumnProvenance::raw(batch.x_col.clone()));
        dataset.schema.columns.insert(
            batch.x_col.clone(),
            ColumnMeta::new(
                batch.x_col.clone(),
                batch.x_unit,
                ColumnOrigin::Imported,
                batch.x_nature,
                x_provenance,
            ),
        );
        for column in &batch.columns {
            let provenance = fitted_column_provenance(&batch, column);
            dataset.schema.columns.insert(
                column.name.clone(),
                ColumnMeta::new(
                    column.name.clone(),
                    column.unit,
                    ColumnOrigin::NumericDerived,
                    column.nature,
                    provenance,
                ),
            );
        }
        let mut bind_nature = |name: &str, nature: ColumnNature| {
            let slot = match nature {
                ColumnNature::Time => &mut dataset.schema.time,
                ColumnNature::Temperature => &mut dataset.schema.temperature,
                ColumnNature::Mass => &mut dataset.schema.mass,
                ColumnNature::Conversion => &mut dataset.schema.alpha,
                ColumnNature::DimensionlessMass => &mut dataset.schema.eta,
                ColumnNature::MassRate => &mut dataset.schema.dm_dt,
                ColumnNature::ConversionRate => &mut dataset.schema.dalpha_dt,
                ColumnNature::DimensionlessMassRate => &mut dataset.schema.deta_dt,
                ColumnNature::TemperatureRate => &mut dataset.schema.dT_dt,
                ColumnNature::ActivationEnergy => &mut dataset.schema.E,
                ColumnNature::R2 => &mut dataset.schema.R2,
                _ => return,
            };
            if slot.is_none() {
                *slot = Some(name.to_string());
            }
        };
        bind_nature(&batch.x_col, batch.x_nature);
        for column in &batch.columns {
            bind_nature(&column.name, column.nature);
        }
        dataset.log_operation(
            "insert_fitted_evaluation_batch",
            AffectedColumns::All,
            None,
            format!(
                "Created from '{}' using {:?} on X='{}'; fitted columns: {}",
                batch.source_experiment_id,
                batch.mode,
                batch.x_col,
                batch
                    .columns
                    .iter()
                    .map(|c| c.name.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            false,
        );
        let mut meta = source.meta;
        meta.id = new_id;
        meta.comment = Some(format!(
            "Reconstructed from {}. {}",
            batch.source_experiment_id,
            meta.comment.unwrap_or_default()
        ));
        self.push(TGAExperiment { dataset, meta })
    }

    /// Fit two numeric columns from one experiment and publish the fitted curve.
    ///
    /// The operation is intentionally scoped to a single experiment: fitting `x`
    /// from one run against `y` from another would silently destroy the physical
    /// meaning of the result. The fitted y values are written back as a new
    /// `FittedCurve` column and carry provenance with model, coefficients and R2.
    pub fn fit_column(
        &mut self,
        request: FitColumnRequest,
    ) -> Result<FitColumnResult, TGADomainError> {
        let output_col = request
            .output_col
            .clone()
            .unwrap_or_else(|| default_fit_output_name(request.kinetic_model, &request.y_col));

        let (x_data, y_data, y_unit) = self.try_apply_by_id(&request.experiment_id, |exp| {
            exp.dataset.ensure_column_name_is_available(&output_col)?;
            materialize_fit_columns(exp, &request.x_col, &request.y_col)
        })?;
        let expected_fit_len = x_data.len();

        let mut fit = Fit::new(request.kinetic_model)
            .with_data(x_data, y_data)
            .with_tolerance(request.tolerance)
            .with_max_iterations(request.max_iter);
        if let Some(initial_guess) = request.initial_guess.clone() {
            fit = fit.with_initial_guess(initial_guess);
        }

        let fit = fit.fit().map_err(fitting_error_to_domain)?;
        let fitted_y = fit.evaluate_solution().map_err(fitting_error_to_domain)?;
        let coefficients = ordered_coefficients(request.kinetic_model, &fit);
        let r2 = fit.r_squared();
        if fitted_y.len() != expected_fit_len
            || fitted_y.iter().any(|value| !value.is_finite())
            || !r2.is_finite()
            || coefficients.iter().any(|(_, value)| !value.is_finite())
        {
            return Err(TGADomainError::InvalidOperation(
                "Fit produced non-finite values or an output length different from its input"
                    .to_string(),
            ));
        }
        let method = format!("{:?}", request.kinetic_model.preferred_method());

        let result = FitColumnResult {
            experiment_id: request.experiment_id.clone(),
            x_col: request.x_col.clone(),
            y_col: request.y_col.clone(),
            output_col: output_col.clone(),
            kinetic_model: request.kinetic_model,
            method,
            coefficients,
            r2,
            tolerance: request.tolerance,
            max_iter: request.max_iter,
        };

        let commit_result = result.clone();
        self.try_transform_by_id(&request.experiment_id, move |exp| {
            let dataset = exp.dataset.add_column_from_vec(
                &commit_result.output_col,
                y_unit,
                ColumnNature::FittedCurve,
                fitted_y,
            )?;
            let mut updated = TGAExperiment {
                dataset,
                meta: exp.meta,
            };
            attach_fit_provenance(&mut updated, &commit_result, y_unit);
            Ok(updated)
        })?;

        Ok(result)
    }
}

fn materialize_fit_columns_many(
    exp: &TGAExperiment,
    x_col: &str,
    y_cols: &[String],
) -> Result<(Vec<f64>, HashMap<String, Vec<f64>>), TGADomainError> {
    let x_meta = exp
        .dataset
        .schema
        .columns
        .get(x_col)
        .ok_or_else(|| TGADomainError::ColumnNotFound(x_col.to_string()))?;
    let mut unique_y = Vec::with_capacity(y_cols.len());
    for name in y_cols {
        if !exp.dataset.schema.columns.contains_key(name) {
            return Err(TGADomainError::ColumnNotFound(name.clone()));
        }
        if name != x_col && !unique_y.contains(name) {
            unique_y.push(name.clone());
        }
    }
    let mut expressions = Vec::with_capacity(unique_y.len() + 1);
    expressions.push(col(x_col));
    expressions.extend(unique_y.iter().map(|name| col(name)));
    let df = exp.dataset.frame.clone().select(expressions).collect()?;
    if df.height() == 0 {
        return Err(TGADomainError::EmptyColumn(x_col.to_string()));
    }
    let mut x_values = Vec::with_capacity(df.height());
    for (index, value) in df.column(&x_meta.name)?.f64()?.iter().enumerate() {
        let value = value.ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Fit input column '{}' contains null at row {}",
                x_col, index
            ))
        })?;
        if !value.is_finite() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Fit input column '{}' contains non-finite value {} at row {}",
                x_col, value, index
            )));
        }
        x_values.push(value);
    }
    let mut y_values = HashMap::with_capacity(y_cols.len());
    for (offset, name) in unique_y.iter().enumerate() {
        let series = df
            .select_at_idx(offset + 1)
            .ok_or_else(|| TGADomainError::ColumnNotFound(name.clone()))?;
        let mut values = Vec::with_capacity(df.height());
        for (index, value) in series.f64()?.iter().enumerate() {
            let value = value.ok_or_else(|| {
                TGADomainError::InvalidOperation(format!(
                    "Fit input column '{}' contains null at row {}",
                    name, index
                ))
            })?;
            if !value.is_finite() {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Fit input column '{}' contains non-finite value {} at row {}",
                    name, value, index
                )));
            }
            values.push(value);
        }
        y_values.insert(name.clone(), values);
    }
    if y_cols.iter().any(|name| name == x_col) {
        y_values.insert(x_col.to_string(), x_values.clone());
    }
    Ok((x_values, y_values))
}

fn materialize_fit_columns(
    exp: &TGAExperiment,
    x_col: &str,
    y_col: &str,
) -> Result<(Vec<f64>, Vec<f64>, Unit), TGADomainError> {
    let x_meta = exp
        .dataset
        .schema
        .columns
        .get(x_col)
        .ok_or_else(|| TGADomainError::ColumnNotFound(x_col.to_string()))?;
    let y_meta = exp
        .dataset
        .schema
        .columns
        .get(y_col)
        .ok_or_else(|| TGADomainError::ColumnNotFound(y_col.to_string()))?;
    let y_unit = y_meta.unit;

    let df = exp
        .dataset
        .frame
        .clone()
        .select([col(x_col), col(y_col)])
        .collect()?;
    if df.height() == 0 {
        return Err(TGADomainError::EmptyColumn(format!("{x_col}/{y_col}")));
    }

    let x = df.column(&x_meta.name)?.f64()?;
    let y = df.column(&y_meta.name)?.f64()?;
    let mut x_data = Vec::with_capacity(df.height());
    let mut y_data = Vec::with_capacity(df.height());

    for (idx, (x_i, y_i)) in x.iter().zip(y.iter()).enumerate() {
        let x_i = x_i.ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Fit input column '{}' contains null at row {}",
                x_col, idx
            ))
        })?;
        let y_i = y_i.ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Fit input column '{}' contains null at row {}",
                y_col, idx
            ))
        })?;
        if !x_i.is_finite() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Fit input column '{}' contains non-finite value {} at row {}",
                x_col, x_i, idx
            )));
        }
        if !y_i.is_finite() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Fit input column '{}' contains non-finite value {} at row {}",
                y_col, y_i, idx
            )));
        }
        x_data.push(x_i);
        y_data.push(y_i);
    }

    if x_data.len() != y_data.len() {
        return Err(TGADomainError::InvalidOperation(format!(
            "Fit input length mismatch: '{}' has {}, '{}' has {}",
            x_col,
            x_data.len(),
            y_col,
            y_data.len()
        )));
    }

    Ok((x_data, y_data, y_unit))
}

fn attach_fit_provenance(exp: &mut TGAExperiment, result: &FitColumnResult, y_unit: Unit) {
    let operation_id = exp
        .dataset
        .history_of_operations
        .vector_of_operations
        .len()
        .checked_sub(1);
    let mut provenance = exp
        .dataset
        .column_provenance(&result.y_col)
        .cloned()
        .unwrap_or_else(|| ColumnProvenance::raw(result.y_col.clone()));
    provenance.rename_output(result.output_col.clone());
    provenance.append_step(ColumnTransformStep::new(
        vec![result.x_col.clone(), result.y_col.clone()],
        result.output_col.clone(),
        ColumnTransformKind::Fitting {
            model: result.kinetic_model.display_name(),
            x_col: result.x_col.clone(),
            y_col: result.y_col.clone(),
            method: result.method.clone(),
            parameters: result.coefficients.clone(),
            r2: result.r2,
            tolerance: result.tolerance,
            max_iter: result.max_iter,
        },
        None,
        operation_id,
        true,
    ));

    if let Some(meta) = exp.dataset.schema.columns.get_mut(&result.output_col) {
        meta.unit = y_unit;
        meta.origin = ColumnOrigin::NumericDerived;
        meta.nature = ColumnNature::FittedCurve;
        meta.provenance = provenance;
    }
}

fn fitted_column_provenance(
    batch: &FittedEvaluationBatch,
    column: &EvaluatedColumn,
) -> ColumnProvenance {
    let mut provenance = ColumnProvenance {
        root_columns: vec![
            batch.source_experiment_id.clone(),
            batch.x_col.clone(),
            column.fit.y_col.clone(),
        ],
        steps: Vec::new(),
    };
    let evaluation = match batch.mode {
        FitEvaluationMode::Extrapolate { start, end, points } => {
            format!("extrapolation [{start}, {end}], {points} points")
        }
        FitEvaluationMode::FillGap {
            start,
            end,
            missing_points,
        } => {
            format!("gap [{start}, {end}], {missing_points} missing points")
        }
    };
    provenance.append_step(ColumnTransformStep::new(
        vec![batch.x_col.clone(), column.fit.y_col.clone()],
        column.name.clone(),
        ColumnTransformKind::Fitting {
            model: column.fit.kinetic_model.display_name(),
            x_col: batch.x_col.clone(),
            y_col: column.fit.y_col.clone(),
            method: format!("{}; {evaluation}", column.fit.method),
            parameters: column.fit.coefficients.clone(),
            r2: column.fit.r2,
            tolerance: column.fit.tolerance,
            max_iter: column.fit.max_iter,
        },
        None,
        None,
        false,
    ));
    provenance
}

fn validate_finite_increasing_range(range: (f64, f64), label: &str) -> Result<(), TGADomainError> {
    if !range.0.is_finite() || !range.1.is_finite() || range.0 >= range.1 {
        return Err(TGADomainError::InvalidOperation(format!(
            "{label} range must have finite, strictly increasing bounds"
        )));
    }
    Ok(())
}

fn model_rmse(fit: &Fit, x: &[f64], y: &[f64]) -> Result<f64, TGADomainError> {
    if x.is_empty() || x.len() != y.len() {
        return Err(TGADomainError::InvalidOperation(
            "Cannot calculate side residuals without aligned observations".to_string(),
        ));
    }
    let predicted = fit.evaluate_at(x).map_err(fitting_error_to_domain)?;
    let mse = predicted
        .iter()
        .zip(y)
        .map(|(pred, observed)| (pred - observed).powi(2))
        .sum::<f64>()
        / y.len() as f64;
    Ok(mse.sqrt())
}

fn ordered_coefficients(model: FittingModelName, fit: &Fit) -> Vec<(String, f64)> {
    let Some(map) = fit.solution_map() else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for name in model.vec_of_params() {
        if let Some(value) = map.get(&name) {
            out.push((name, *value));
        }
    }
    for (name, value) in map {
        if !out.iter().any(|(existing, _)| existing == &name) {
            out.push((name, value));
        }
    }
    out
}

fn default_fit_output_name(model: FittingModelName, y_col: &str) -> String {
    format!(
        "fit_{}_{}",
        model.output_name_fragment(),
        sanitize_column_fragment(y_col)
    )
}

fn sanitize_column_fragment(raw: &str) -> String {
    let sanitized = raw
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect::<String>()
        .trim_matches('_')
        .to_string();
    if sanitized.is_empty() {
        "y".to_string()
    } else {
        sanitized
    }
}

fn fitting_error_to_domain(err: KineticFittingError) -> TGADomainError {
    TGADomainError::InvalidOperation(err.to_string())
}

impl Default for FittingModelName {
    fn default() -> Self {
        Self::Linear
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::experiment_series_main::{
        ExperimentMeta, TGASeries,
    };
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnNature;
    use crate::Kinetics::experimental_kinetics::testing_mod::{
        AdvancedTGAConfig, ExperimentMode, KineticModel, VirtualTGA,
    };
    use approx::assert_relative_eq;

    fn linspace(start: f64, end: f64, len: usize) -> Vec<f64> {
        let step = if len > 1 {
            (end - start) / (len - 1) as f64
        } else {
            0.0
        };
        (0..len).map(|idx| start + idx as f64 * step).collect()
    }

    fn fit_case(
        model: FittingModelName,
        x_data: Vec<f64>,
        y_data: Vec<f64>,
        initial_guess: Vec<f64>,
    ) -> Fit {
        Fit::new(model)
            .with_data(x_data, y_data)
            .with_initial_guess(initial_guess)
            .fit()
            .expect("kinetic model should fit")
    }

    fn assert_map_value(map: &HashMap<String, f64>, key: &str, expected: f64, eps: f64) {
        assert!(
            map.contains_key(key),
            "solution map should contain key {key}"
        );
        assert_relative_eq!(map[key], expected, epsilon = eps);
    }
    #[allow(dead_code)]
    fn assert_expression_values(
        model: FittingModelName,
        x_data: &[f64],
        params: &[(&str, f64)],
        expected: &[f64],
    ) {
        let param_map: HashMap<String, f64> = params
            .iter()
            .map(|(name, value)| ((*name).to_string(), *value))
            .collect();
        let expr = model.equation().set_variable_from_map(&param_map);
        let eval = expr.lambdify1D();
        let actual: Vec<f64> = x_data.iter().map(|&x| eval(x)).collect();

        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(actual, expected, epsilon = 1e-10);
        }
    }

    const R: f64 = 8.314;

    fn arrhenius_rate(k0: f64, e: f64, temperature: f64) -> f64 {
        k0 * (-e / (R * temperature)).exp()
    }

    fn series_from_advanced_config(cfg: &AdvancedTGAConfig) -> TGASeries {
        let mut series = TGASeries::new();
        series
            .create_series_from_synthetic_data(cfg)
            .expect("advanced synthetic experiment should be inserted");
        series
    }

    fn single_component_fit_series() -> (TGASeries, String, f64, f64) {
        let temperature = 600.0;
        let m0 = 12.0;
        let k0 = 2.0e5;
        let e = 80_000.0;
        let cfg = AdvancedTGAConfig {
            n_points: 120,
            dt: 0.5,
            experiment_mode: ExperimentMode::Isothermal {
                temperatures: vec![temperature],
            },
            kinetic_model: KineticModel::ArrheniusSingle { m0, k0, e, r: R },
            mass_noise: None,
            temp_noise: None,
            spikes: None,
            seed: 11,
        };
        let exp_id = format!("T = {:.1} K", temperature);
        let k_eff = arrhenius_rate(k0, e, temperature);
        (series_from_advanced_config(&cfg), exp_id, m0, k_eff)
    }

    fn two_component_fit_series() -> (TGASeries, String, Vec<(f64, f64)>) {
        let temperature = 620.0;
        let m01 = 7.0;
        let k01 = 1.4e4;
        let e1 = 70_000.0;
        let m02 = 5.0;
        let k02 = 1.2e3;
        let e2 = 65_000.0;
        let cfg = AdvancedTGAConfig {
            n_points: 160,
            dt: 0.5,
            experiment_mode: ExperimentMode::Isothermal {
                temperatures: vec![temperature],
            },
            kinetic_model: KineticModel::ArrheniusTwoComponent {
                m01,
                k01,
                e1,
                m02,
                k02,
                e2,
                r: R,
            },
            mass_noise: None,
            temp_noise: None,
            spikes: None,
            seed: 12,
        };
        let exp_id = format!("T = {:.1} K", temperature);
        let mut expected = vec![
            (m01, arrhenius_rate(k01, e1, temperature)),
            (m02, arrhenius_rate(k02, e2, temperature)),
        ];
        expected.sort_by(|left, right| left.1.partial_cmp(&right.1).unwrap());
        (series_from_advanced_config(&cfg), exp_id, expected)
    }

    fn synthetic_series_with_invalid_mass() -> (TGASeries, String) {
        let (mut series, exp_id, _, _) = single_component_fit_series();
        series
            .try_mutate_by_id(&exp_id, |exp| {
                exp.dataset.frame =
                    exp.dataset
                        .frame
                        .clone()
                        .with_columns([when(col("time").eq(lit(2.5)))
                            .then(lit(f64::NAN))
                            .otherwise(col("mass"))
                            .alias("mass")]);
                Ok(())
            })
            .expect("test setup should inject invalid mass value");
        (series, exp_id)
    }

    fn empty_synthetic_series() -> (TGASeries, String) {
        let exp_id = "empty_fit".to_string();
        let mut series = TGASeries::new();
        series
            .create_from_synthetic_data(
                &crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA {
                    time: Vec::new(),
                    temperature: Vec::new(),
                    mass: Vec::new(),
                },
                ExperimentMeta {
                    id: exp_id.clone(),
                    heating_rate: None,
                    isothermal_temperature: Some(600.0),
                    comment: Some("empty fitting negative test".to_string()),
                },
            )
            .expect("empty synthetic dataframe is valid test setup");
        (series, exp_id)
    }

    fn fit_request(exp_id: &str, model: FittingModelName) -> FitColumnRequest {
        FitColumnRequest::new(exp_id, "time", "mass", model)
    }

    fn assert_invalid_operation_contains(err: TGADomainError, needle: &str) {
        match err {
            TGADomainError::InvalidOperation(msg) => {
                assert!(
                    msg.contains(needle),
                    "expected invalid operation to contain '{needle}', got '{msg}'"
                );
            }
            other => panic!("expected InvalidOperation, got {other:?}"),
        }
    }

    fn sorted_two_exp_pairs(coefficients: &[(String, f64)]) -> Vec<(f64, f64)> {
        let mut map = HashMap::new();
        for (name, value) in coefficients {
            map.insert(name.as_str(), *value);
        }
        let mut pairs = vec![
            (*map.get("a").unwrap(), *map.get("p").unwrap()),
            (*map.get("b").unwrap(), *map.get("k").unwrap()),
        ];
        pairs.sort_by(|left, right| left.1.partial_cmp(&right.1).unwrap());
        pairs
    }

    macro_rules! kinetic_test {
        ($name:ident, $model:expr, $x:expr, $y:expr, $guess:expr, $( $key:expr => $expected:expr ),+ $(,)?) => {
            #[test]
            fn $name() {
                let fit = fit_case($model, $x, $y, $guess);
                let map = fit.solution_map().expect("fit should expose a solution map");
                $(assert_map_value(&map, $key, $expected, 1e-6);)+
                assert!(fit.r_squared() > 0.999);
                let y_pred = fit.evaluate_solution().expect("fit should evaluate");
                assert_eq!(y_pred.len(), fit.x_data.len());
            }
        };
    }

    kinetic_test!(
        autocatalytic_fits,
        FittingModelName::Autocatalytic,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| 2.0 * v * (1.0 - v)).collect::<Vec<_>>()
        },
        vec![1.5],
        "k" => 2.0
    );

    kinetic_test!(
        first_order_fits,
        FittingModelName::FirstOrder,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| 1.8 * (1.0 - v)).collect::<Vec<_>>()
        },
        vec![1.0],
        "k" => 1.8
    );

    kinetic_test!(
        second_order_fits,
        FittingModelName::SecondOrder,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| 1.4 * (1.0 - v).powi(2)).collect::<Vec<_>>()
        },
        vec![1.0],
        "k" => 1.4
    );

    kinetic_test!(
        third_order_fits,
        FittingModelName::ThirdOrder,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| 1.2 * (1.0 - v).powi(3)).collect::<Vec<_>>()
        },
        vec![1.0],
        "k" => 1.2
    );

    #[test]
    fn sestak_berggren_smoke_test() {
        let model = FittingModelName::SestakBerggren;
        let expr = model.equation();
        assert!(!expr.to_string().is_empty());
        assert_eq!(model.vec_of_params(), vec!["k", "n", "m", "p"]);
        assert_eq!(model.vec_of_initial_guess().len(), 4);
    }

    #[test]
    fn jma_smoke_test() {
        let model = FittingModelName::JohnsonMehlAvrami;
        let expr = model.equation();
        assert!(!expr.to_string().is_empty());
        assert_eq!(model.vec_of_params(), vec!["m"]);
        assert_eq!(model.vec_of_initial_guess().len(), 1);
    }

    #[test]
    fn exponential_models_prefer_varpro_backend() {
        assert_eq!(FittingModelName::DecExp.preferred_method(), Method::VARPRO);
        assert_eq!(FittingModelName::TwoExp.preferred_method(), Method::VARPRO);
        assert_eq!(
            FittingModelName::ThreeExp.preferred_method(),
            Method::VARPRO
        );
        assert!(FittingModelName::Linear.uses_classic_lm());
        assert!(FittingModelName::FirstOrder.uses_classic_lm());
    }

    #[test]
    fn dec_exp_uses_default_initial_guess_when_not_provided() {
        let x = linspace(0.0, 5.0, 50);
        let y = x
            .iter()
            .map(|&v| 2.0 * (-0.6 * v).exp())
            .collect::<Vec<_>>();

        let fit = Fit::new(FittingModelName::DecExp)
            .with_data(x, y)
            .fit()
            .expect("default initial guess should be enough for the VarPro path");

        let map = fit
            .solution_map()
            .expect("fit should expose a solution map");
        assert_relative_eq!(map["c"], 2.0, epsilon = 1e-6);
        assert_relative_eq!(map["k"], 0.6, epsilon = 1e-6);
        assert!(fit.r_squared() > 0.999_999);
    }

    #[test]
    fn varpro_fit_retains_requested_tolerance_and_iteration_limit() {
        let x = linspace(0.0, 5.0, 50);
        let y = x
            .iter()
            .map(|&value| 2.0 * (-0.6 * value).exp())
            .collect::<Vec<_>>();
        let tolerance = 2.5e-9;
        let max_iter = 137;

        let fit = Fit::new(FittingModelName::DecExp)
            .with_data(x, y)
            .with_tolerance(tolerance)
            .with_max_iterations(max_iter)
            .fit()
            .expect("VARPRO should fit the exact exponential data");

        assert_eq!(fit.tolerance, tolerance);
        assert_eq!(fit.max_iter, max_iter);
        assert!(fit.r_squared() > 0.999_999);
    }

    #[test]
    fn sestak_berggren_expression_parameters_are_ordered() {
        assert_eq!(
            FittingModelName::SestakBerggren.vec_of_params(),
            vec![
                "k".to_string(),
                "n".to_string(),
                "m".to_string(),
                "p".to_string()
            ]
        );
    }

    #[test]
    fn jma_expression_parameters_are_ordered() {
        assert_eq!(
            FittingModelName::JohnsonMehlAvrami.vec_of_params(),
            vec!["m".to_string()]
        );
    }

    kinetic_test!(
        acceleratory_fits,
        FittingModelName::Acceleratory,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| 1.6 * (1.0 - v).powf(1.0 - 1.0 / 1.6)).collect::<Vec<_>>()
        },
        vec![1.0],
        "m" => 1.6
    );

    kinetic_test!(
        deceleratory_fits,
        FittingModelName::Deceleratory,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| (1.0 - v).powf(1.4)).collect::<Vec<_>>()
        },
        vec![1.0],
        "m" => 1.4
    );

    kinetic_test!(
        sestac_bergren_trunc_fits,
        FittingModelName::SestacBergrenTrunc,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter().map(|&v| v.powf(1.1) * (1.0 - v).powf(0.8)).collect::<Vec<_>>()
        },
        vec![1.0, 1.0],
        "m" => 1.1,
        "n" => 0.8
    );

    kinetic_test!(
        sestac_bergren_param_fits,
        FittingModelName::SestacBergrenParam,
        linspace(0.05, 0.85, 50),
        {
            let x = linspace(0.05, 0.85, 50);
            x.iter()
                .map(|&v| 1.5 * v.powf(0.9) * (1.0 - v).powf(0.7))
                .collect::<Vec<_>>()
        },
        vec![1.0, 1.0, 1.0],
        "m" => 0.9,
        "n" => 0.7,
        "c" => 1.5
    );

    kinetic_test!(
        dec_exp_fits,
        FittingModelName::DecExp,
        linspace(0.0, 5.0, 50),
        {
            let x = linspace(0.0, 5.0, 50);
            x.iter().map(|&v| 2.0 * (-0.6 * v).exp()).collect::<Vec<_>>()
        },
        vec![1.5, 0.3],
        "c" => 2.0,
        "k" => 0.6
    );

    kinetic_test!(
        two_exp_fits,
        FittingModelName::TwoExp,
        linspace(0.0, 5.0, 50),
        {
            let x = linspace(0.0, 5.0, 50);
            x.iter()
                .map(|&v| 1.1 * (-0.8 * v).exp() + 0.7 * (-0.25 * v).exp())
                .collect::<Vec<_>>()
        },
        vec![1.0, 0.5, 1.0, 0.2],
        "a" => 1.1,
        "p" => 0.8,
        "b" => 0.7,
        "k" => 0.25
    );

    kinetic_test!(
        three_exp_fits,
        FittingModelName::ThreeExp,
        linspace(0.0, 5.0, 60),
        {
            let x = linspace(0.0, 5.0, 60);
            x.iter()
                .map(|&v| 1.0 * (-0.9 * v).exp() + 0.8 * (-0.3 * v).exp() + 0.5 * (-0.1 * v).exp())
                .collect::<Vec<_>>()
        },
        vec![1.0, 0.5, 1.0, 0.2, 0.4, 0.1],
        "a" => 1.0,
        "p" => 0.9,
        "b" => 0.8,
        "k" => 0.3,
        "c" => 0.5,
        "r" => 0.1
    );

    kinetic_test!(
        linear_fits,
        FittingModelName::Linear,
        linspace(0.0, 10.0, 50),
        {
            let x = linspace(0.0, 10.0, 50);
            x.iter().map(|&v| 2.5 * v + 1.25).collect::<Vec<_>>()
        },
        vec![1.0, 1.0],
        "k" => 2.5,
        "b" => 1.25
    );

    kinetic_test!(
        polynomial_fits,
        FittingModelName::Polynom { n: 3 },
        linspace(0.0, 3.0, 50),
        {
            let x = linspace(0.0, 3.0, 50);
            x.iter()
                .map(|&v| 5.0 * v.powi(3) + 2.0 * v.powi(2) + 3.0 * v + 1.0)
                .collect::<Vec<_>>()
        },
        vec![1.0; 4],
        "c3" => 5.0,
        "c2" => 2.0,
        "c1" => 3.0,
        "c0" => 1.0
    );

    #[test]
    fn series_fit_column_dec_exp_adds_fitted_column_and_provenance() {
        let (mut series, exp_id, expected_c, expected_k) = single_component_fit_series();

        let result = series
            .fit_column(fit_request(&exp_id, FittingModelName::DecExp).with_output_col("fit_mass"))
            .expect("series-level single exponential fit should converge");

        assert_eq!(result.output_col, "fit_mass");
        assert!(result.r2 > 0.999_999);
        let coeffs: HashMap<_, _> = result.coefficients.iter().cloned().collect();
        assert_relative_eq!(coeffs["c"], expected_c, epsilon = 1e-8);
        assert_relative_eq!(coeffs["k"], expected_k, epsilon = 1e-8);

        let exp = series.get_experiment_by_id(&exp_id).unwrap();
        let meta = exp.dataset.schema.columns.get("fit_mass").unwrap();
        assert_eq!(meta.nature, ColumnNature::FittedCurve);
        let df = exp
            .dataset
            .frame
            .clone()
            .select([col("mass"), col("fit_mass")])
            .collect()
            .unwrap();
        assert_eq!(
            df.column("mass").unwrap().len(),
            df.column("fit_mass").unwrap().len()
        );

        let provenance = exp
            .dataset
            .column_provenance_text("fit_mass")
            .expect("fitted column should have provenance");
        assert!(provenance.contains("fitting(model=DecExp"));
        assert!(provenance.contains("R2="));
        assert!(provenance.contains("params=[c="));

        series.undo_last(&exp_id).unwrap();
        assert!(
            !series
                .get_experiment_by_id(&exp_id)
                .unwrap()
                .dataset
                .schema
                .columns
                .contains_key("fit_mass")
        );
    }

    #[test]
    fn series_fit_column_two_exp_handles_component_permutation() {
        let (mut series, exp_id, expected) = two_component_fit_series();

        let result = series
            .fit_column(
                fit_request(&exp_id, FittingModelName::TwoExp)
                    .with_initial_guess(vec![1.0, 0.5, 1.0, 0.2]),
            )
            .expect("series-level two-exponential fit should converge");

        assert_eq!(result.output_col, "fit_twoexp_mass");
        assert!(result.r2 > 0.999_999);

        let actual = sorted_two_exp_pairs(&result.coefficients);
        for ((actual_a, actual_k), (expected_a, expected_k)) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(actual_a, expected_a, epsilon = 1e-8);
            assert_relative_eq!(actual_k, expected_k, epsilon = 1e-8);
        }

        let exp = series.get_experiment_by_id(&exp_id).unwrap();
        assert!(exp.dataset.schema.columns.contains_key("fit_twoexp_mass"));
    }

    #[test]
    fn series_fit_column_rejects_non_finite_input() {
        let (mut series, exp_id) = synthetic_series_with_invalid_mass();

        let err = series
            .fit_column(fit_request(&exp_id, FittingModelName::DecExp))
            .unwrap_err();

        assert!(format!("{:?}", err).contains("non-finite"));
    }

    #[test]
    fn series_fit_column_negative_cases_are_reported() {
        let (mut series, exp_id, _, _) = single_component_fit_series();

        assert!(matches!(
            series
                .fit_column(fit_request("missing", FittingModelName::DecExp))
                .unwrap_err(),
            TGADomainError::ExperimentNotFound(_)
        ));
        assert!(matches!(
            series
                .fit_column(FitColumnRequest::new(
                    &exp_id,
                    "missing_time",
                    "mass",
                    FittingModelName::DecExp
                ))
                .unwrap_err(),
            TGADomainError::ColumnNotFound(_)
        ));

        series
            .fit_column(fit_request(&exp_id, FittingModelName::DecExp).with_output_col("fit_mass"))
            .expect("first fit should create output column");
        assert!(matches!(
            series
                .fit_column(
                    fit_request(&exp_id, FittingModelName::DecExp).with_output_col("fit_mass")
                )
                .unwrap_err(),
            TGADomainError::ColumnAlreadyExists(_)
        ));

        let before = series.get_experiment_by_id(&exp_id).unwrap();
        let before_frame = before.dataset.frame.clone().collect().unwrap();
        let before_provenance = before
            .dataset
            .schema
            .columns
            .iter()
            .map(|(name, meta)| (name.clone(), format!("{:?}", meta.provenance)))
            .collect::<std::collections::HashMap<_, _>>();
        let before_columns = before.dataset.schema.columns.len();
        let before_history_lines = before.dataset.history_of_operations.feed_lines();
        let before_history = before
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();
        let before_undo_depth = before.dataset.undo_depth();
        let err = series
            .fit_column(
                fit_request(&exp_id, FittingModelName::DecExp)
                    .with_output_col("fit_bad_guess")
                    .with_initial_guess(Vec::new()),
            )
            .unwrap_err();
        assert_invalid_operation_contains(err, "initial guess length");
        let after = series.get_experiment_by_id(&exp_id).unwrap();
        assert_eq!(after.dataset.frame.clone().collect().unwrap(), before_frame);
        assert_eq!(after.dataset.schema.columns.len(), before_columns);
        let after_provenance = after
            .dataset
            .schema
            .columns
            .iter()
            .map(|(name, meta)| (name.clone(), format!("{:?}", meta.provenance)))
            .collect::<std::collections::HashMap<_, _>>();
        assert_eq!(after_provenance, before_provenance);
        assert_eq!(
            after.dataset.history_of_operations.feed_lines(),
            before_history_lines
        );
        assert_eq!(
            after
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            before_history
        );
        assert_eq!(after.dataset.undo_depth(), before_undo_depth);
        assert!(!after.dataset.schema.columns.contains_key("fit_bad_guess"));
    }

    #[test]
    fn series_fit_column_rejects_empty_input() {
        let (mut series, exp_id) = empty_synthetic_series();

        assert!(matches!(
            series
                .fit_column(fit_request(&exp_id, FittingModelName::DecExp))
                .unwrap_err(),
            TGADomainError::EmptyColumn(_)
        ));
    }

    #[test]
    fn staged_gap_reconstruction_fits_both_sides_and_commits_complete_new_experiment() {
        let observed_x = vec![0.0, 1.0, 2.0, 3.0, 7.0, 8.0, 9.0, 10.0];
        let virtual_tga = VirtualTGA {
            time: observed_x.clone(),
            temperature: observed_x.iter().map(|x| 300.0 + 4.0 * x).collect(),
            mass: observed_x.iter().map(|x| 12.0 - 0.5 * x).collect(),
        };
        let mut series = TGASeries::new();
        series
            .create_from_synthetic_data(
                &virtual_tga,
                ExperimentMeta {
                    id: "gapped_source".into(),
                    ..ExperimentMeta::default()
                },
            )
            .unwrap();
        let source_history = series
            .get_experiment_by_id("gapped_source")
            .unwrap()
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();

        let temperature = series
            .fit_column_for_evaluation(
                FitColumnRequest::new(
                    "gapped_source",
                    "time",
                    "temperature",
                    FittingModelName::Linear,
                )
                .with_output_col("temperature_rebuilt")
                .excluding_range(3.0, 7.0),
            )
            .unwrap();
        let mass = series
            .fit_column_for_evaluation(
                FitColumnRequest::new("gapped_source", "time", "mass", FittingModelName::Linear)
                    .with_output_col("mass_rebuilt")
                    .excluding_range(3.0, 7.0),
            )
            .unwrap();
        let mut duplicate_output = mass.clone();
        duplicate_output.result.output_col = temperature.result.output_col.clone();
        assert!(
            evaluate_fit_batch(
                &[temperature.clone(), duplicate_output],
                FitEvaluationMode::FillGap {
                    start: 3.0,
                    end: 7.0,
                    missing_points: 3,
                },
            )
            .is_err()
        );
        assert!(temperature.result.r2 > 0.999_999);
        assert!(mass.result.r2 > 0.999_999);
        assert!(temperature.left_rmse.unwrap() < 1e-8);
        assert!(temperature.right_rmse.unwrap() < 1e-8);

        let batch = evaluate_fit_batch(
            &[temperature, mass],
            FitEvaluationMode::FillGap {
                start: 3.0,
                end: 7.0,
                missing_points: 3,
            },
        )
        .unwrap();
        assert_eq!(batch.columns.len(), 2);
        assert!(batch.x_values.windows(2).all(|pair| pair[1] > pair[0]));
        assert_eq!(
            batch
                .x_values
                .iter()
                .filter(|x| **x > 3.0 && **x < 7.0)
                .count(),
            3
        );
        assert_eq!(batch.columns[0].values.len(), batch.columns[1].values.len());
        assert!(
            batch
                .columns
                .iter()
                .all(|column| column.values.iter().all(|value| value.is_finite()))
        );

        // Preparing and evaluating a fit is non-mutating; the explicit commit creates a sibling run.
        let source = series.get_experiment_by_id("gapped_source").unwrap();
        assert_eq!(
            source
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            source_history
        );
        assert_eq!(
            source.dataset.frame.clone().collect().unwrap().height(),
            observed_x.len()
        );
        let replay = batch.clone();
        let mut invalid_batch = batch.clone();
        invalid_batch.x_values[0] = f64::NAN;
        let count_before_invalid_commit = series.len();
        assert!(
            series
                .push_fitted_evaluation_batch(&invalid_batch, "invalid_reconstruction")
                .is_err()
        );
        assert_eq!(series.len(), count_before_invalid_commit);
        series
            .push_fitted_evaluation_batch(&batch, "gapped_reconstruction")
            .unwrap();
        let series_len = series.len();
        assert!(matches!(
            series.push_fitted_evaluation_batch(&replay, "gapped_reconstruction"),
            Err(TGADomainError::DuplicateExperimentId(_))
        ));
        assert_eq!(series.len(), series_len);
        let reconstructed = series
            .get_experiment_by_id("gapped_reconstruction")
            .unwrap();
        let frame = reconstructed.dataset.frame.clone().collect().unwrap();
        assert_eq!(frame.height(), observed_x.len() + 3);
        assert_eq!(frame.width(), 3);
        assert!(
            frame
                .columns()
                .iter()
                .all(|column| column.null_count() == 0)
        );
        assert_eq!(
            reconstructed
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            1
        );
        assert!(
            reconstructed
                .dataset
                .history_of_operations
                .vector_of_operations[0]
                .description
                .contains("FillGap { start: 3.0, end: 7.0, missing_points: 3 }")
        );
        assert!(
            reconstructed
                .dataset
                .column_provenance_text("mass_rebuilt")
                .unwrap()
                .contains("gap [3")
        );
        assert_eq!(
            series
                .get_experiment_by_id("gapped_source")
                .unwrap()
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap()
                .height(),
            observed_x.len()
        );
    }

    #[test]
    fn evaluated_batch_is_rejected_if_source_history_changes_before_commit() {
        let (mut series, exp_id, _, _) = single_component_fit_series();
        let source = series.get_experiment_by_id(&exp_id).unwrap();
        let x = source
            .dataset
            .frame
            .clone()
            .select([col("time")])
            .collect()
            .unwrap();
        let max_x = x.column("time").unwrap().f64().unwrap().max().unwrap();
        let session = series
            .fit_column_for_evaluation(
                fit_request(&exp_id, FittingModelName::DecExp).with_output_col("mass_rebuilt"),
            )
            .unwrap();
        let batch = evaluate_fit_batch(
            &[session],
            FitEvaluationMode::Extrapolate {
                start: max_x + 1.0,
                end: max_x + 2.0,
                points: 3,
            },
        )
        .unwrap();
        series
            .try_mutate_by_id(&exp_id, |experiment| {
                experiment.dataset.log_operation(
                    "test_source_revision",
                    AffectedColumns::All,
                    None,
                    "Simulate a source edit after fitting".to_string(),
                    false,
                );
                Ok(())
            })
            .unwrap();
        let count_before = series.len();
        assert!(matches!(
            series.push_fitted_evaluation_batch(&batch, "stale_reconstruction"),
            Err(TGADomainError::InvalidOperation(message)) if message.contains("changed after fitting")
        ));
        assert_eq!(series.len(), count_before);
        assert!(series.get_experiment_by_id("stale_reconstruction").is_err());
    }
}
