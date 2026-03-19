//! # Integral isoconversional kinetic methods: OFW, KAS, Starink
//!
//! ## What this module does
//! Implements three classical **integral isoconversional** methods for
//! model-free kinetic analysis of non-isothermal TGA data:
//!
//! | Method | Full name | Year |
//! |--------|-----------|------|
//! | OFW    | Ozawa–Flynn–Wall | 1965/1966 |
//! | KAS    | Kissinger–Akahira–Sunose | 1957/1971 |
//! | Starink | Starink | 1996 |
//!
//! All three share the same linear regression backbone: at each fixed
//! conversion level α a straight line is fitted to `log(f(β, Tα))` vs `1/Tα`
//! across experiments run at different heating rates β.  The slope of that
//! line is proportional to `-Eα/R`.
//!
//! Also hosts the [`VyazovkinMethod`] adapter that delegates to
//! `Vyazovkin::VyazovkinSolver` through the same `KineticMethod` trait.
//!
//! ## Pipeline
//! ```text
//! KineticDataView  (one ExperimentData per heating rate)
//!      │
//!      ▼
//! ConversionGridBuilder  →  ConversionGrid  (ndarray matrices, shape [n_exp × n_α])
//!      │
//!      ▼
//! IntegralIsoconversionalSolver::solve
//!      │
//!      ▼
//! IsoconversionalResult  (Vec<IsoLayerResult>, one per α point)
//! ```
//!
//! ## Main data structures
//! - [`IntIsoconversionalMethod`] — enum selecting KAS / OFW / Starink.
//! - [`IntegralIsoConfig`] — three numbers that fully parameterise a method:
//!   `exponent` (power of T in the denominator), `log_kind` (ln or log₁₀),
//!   `slope_factor` (converts regression slope to Eα in J/mol).
//! - [`IntegralIsoconversionalSolver`] — the core solver; built from a config
//!   via factory methods `kas()`, `ofw()`, `starink()`.
//! - [`IsoLayerResult`] — result for one α layer: `eta`, `ea`, optional `k`,
//!   and the full [`LinearRegressionResult`].
//! - [`IsoconversionalResult`] — collects all layers plus the method name;
//!   can be exported to a Polars `DataFrame` or pretty-printed.
//! - [`OFW`], [`KAS`], [`VyazovkinMethod`] — zero-size structs implementing
//!   `KineticMethod` as entry points for the unified pipeline.
//!
//! ## Math
//! All three methods approximate the Arrhenius temperature integral
//! `∫ exp(-E/RT) dT` and linearise the result:
//!
//! | Method | Linearised form | Slope |
//! |--------|-----------------|-------|
//! | OFW    | `ln(β)` vs `1/T` | `-1.052 E/R` |
//! | KAS    | `ln(β/T²)` vs `1/T` | `-E/R` |
//! | Starink| `ln(β/T^1.92)` vs `1/T` | `-E/(1.0008 R)` |
//!
//! The `exponent` and `slope_factor` fields of [`IntegralIsoConfig`] encode
//! exactly these differences, so a single `solve` loop handles all three.
//!
//! ## Non-trivial technique
//! The `1/T` column is pre-cached in `ConversionGrid::inv_temperature` so the
//! inner loop never performs a division.  The y-values (`ln(β/T^b)`) are
//! computed on-the-fly per layer rather than stored, keeping memory O(n_exp)
//! instead of O(n_exp × n_α).

use super::super::kinetic_methods::{
    ConversionGrid, ConversionGridBuilder, GridInterpolation, KineticDataView, KineticMethod,
    KineticRequirements, TGADomainError,
};
use super::Vyazovkin::VyazovkinSolver;
use super::kinetic_regression::{LinearRegressionResult, linear_regression};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    AffectedColumns, ColumnMeta, ColumnNature, ColumnOrigin, History, TGADataset, TGASchema, Unit,
};
use log::info;
use ndarray::Array1;
use polars::prelude::*;
use std::collections::HashMap;
use std::time::Instant;
pub const R: f64 = 8.314;
//====================================================================================================================
//          METHODS
//======================================================================================================================
/// Selects which integral isoconversional method the unified solver will run.
#[derive(Clone, Copy, Debug)]
pub enum IntIsoconversionalMethod {
    /// Kissinger–Akahira–Sunose: regresses `ln(β/T²)` vs `1/T`.
    KAS,
    /// Ozawa–Flynn–Wall: regresses `ln(β)` vs `1/T`.
    OFW,
    /// Starink: regresses `ln(β/T^1.92)` vs `1/T` with a refined slope factor.
    Starink,
}

//========================================================================================================================
//HANDLING RESULT
//=====================================================================================================================
// "результат слоя α".
/// Result for a single isoconversional layer (one fixed α value).
///
/// Produced by every solver in this crate; the `regression` field is a
/// placeholder with default values for methods that do not use linear
/// regression (e.g. Vyazovkin).
#[derive(Clone, Debug)]
pub struct IsoLayerResult {
    /// Conversion value α this layer corresponds to.
    pub eta: f64,
    /// Activation energy Eα (J/mol) estimated at this conversion.
    pub ea: f64,
    /// Pre-exponential factor (optional; not computed by all methods).
    pub k: Option<f64>,
    /// Linear regression diagnostics (slope, intercept, R²).
    pub regression: LinearRegressionResult,
}

/// Aggregated output of any isoconversional method.
///
/// Contains one [`IsoLayerResult`] per α grid point and the human-readable
/// method name.  Can be exported to a Polars `DataFrame` via `to_dataframe`
/// or inspected with `pretty_print_and_assert`.
// Русское название: "итоговый результат метода".
#[derive(Clone, Debug)]
pub struct IsoconversionalResult {
    /// Short name of the method that produced this result (e.g. `"KAS"`).
    pub method: &'static str,
    /// One result entry per α layer, ordered by increasing α.
    pub layers: Vec<IsoLayerResult>,
}

impl IsoconversionalResult {
    /// Converts the result to a Polars eager `DataFrame` with columns
    /// `alpha`, `Ea` (J/mol), and `R2`.
    pub fn to_dataframe(&self) -> DataFrame {
        let eta: Vec<f64> = self.layers.iter().map(|l| l.eta).collect();

        let ea: Vec<f64> = self.layers.iter().map(|l| l.ea).collect();

        let r2: Vec<f64> = self.layers.iter().map(|l| l.regression.r2).collect();

        DataFrame::new_infer_height(vec![
            Series::new("eta".into(), eta).into(),
            Series::new("Ea".into(), ea).into(),
            Series::new("R2".into(), r2).into(),
        ])
        .unwrap()
    }

    /// Converts the result to a TGADataset with proper schema metadata.
    /// Columns:
    /// - alpha (dimensionless conversion)
    /// - Ea (activation energy, J/mol; unit left as Unknown)
    /// - R2 (dimensionless)
    ///
    /// Column origin is set to NumericDerived and ColumnNature is filled.
    /// An operation record is added to the dataset history.
    pub fn to_tga_dataset(&self) -> Result<TGADataset, TGADomainError> {
        let eta: Vec<f64> = self.layers.iter().map(|l| l.eta).collect();
        let ea: Vec<f64> = self.layers.iter().map(|l| l.ea).collect();
        let r2: Vec<f64> = self.layers.iter().map(|l| l.regression.r2).collect();
        info!("found activasion energy vector of lengh {}", ea.len());
        for (i, eta) in eta.iter().enumerate() {
            println!("eta {:.2},Ea {:2}, r2 {:2}", eta, ea[i], r2[i]);
        }

        let df = DataFrame::new_infer_height(vec![
            Series::new("eta".into(), eta).into(),
            Series::new("Ea".into(), ea).into(),
            Series::new("R2".into(), r2).into(),
        ])?;

        let mut columns: HashMap<String, ColumnMeta> = HashMap::new();
        columns.insert(
            "eta".to_string(),
            ColumnMeta {
                name: "eta".to_string(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::NumericDerived,
                nature: ColumnNature::Conversion,
            },
        );
        columns.insert(
            "Ea".to_string(),
            ColumnMeta {
                name: "Ea".to_string(),
                unit: Unit::Unknown,
                origin: ColumnOrigin::NumericDerived,
                nature: ColumnNature::ActivationEnergy,
            },
        );
        columns.insert(
            "R2".to_string(),
            ColumnMeta {
                name: "R2".to_string(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::NumericDerived,
                nature: ColumnNature::R2,
            },
        );

        let schema = TGASchema {
            columns,
            time: None,
            temperature: None,
            mass: None,
            alpha: Some("alpha".to_string()),
            dm_dt: None,
            eta: None,
            deta_dt: None,
            dalpha_dt: None,
            dT_dt: None,
        };

        let mut dataset = TGADataset {
            frame: df.lazy(),
            schema,
            oneframeplot: None,
            history_of_operations: History::new(),
        };

        dataset.log_operation(
            "create_from_isoconversional_result",
            AffectedColumns::Specific(vec![
                "alpha".to_string(),
                "Ea".to_string(),
                "R2".to_string(),
            ]),
            None,
            format!(
                "Created dataset from isoconversional result ({})",
                self.method
            ),
            false,
        );

        Ok(dataset)
    }

    /// Pretty‑print a table of eta and Ea values within the given conversion range.
    /// Only layers with eta between `eta_min` and `eta_max` are considered.
    /// The output is limited to approximately `n` rows (sampled evenly).
    /// Pretty‑print with assertion that R² exceeds a threshold.
    pub fn pretty_print_and_assert(
        &self,
        eta_min: f64,
        eta_max: f64,
        n: usize,
        threshold: Option<f64>,
    ) {
        use tabled::{Table, Tabled, settings::Style};

        #[derive(Tabled)]
        struct EaRow {
            #[tabled(rename = "η")]
            eta: String,
            #[tabled(rename = "Ea (kJ/mol)")]
            ea: String,
            #[tabled(rename = "R²")]
            r2: String,
        }

        let filtered: Vec<&IsoLayerResult> = self
            .layers
            .iter()
            .filter(|layer| layer.eta >= eta_min && layer.eta <= eta_max)
            .collect();

        if filtered.is_empty() {
            println!("No layers in the range [{}, {}]", eta_min, eta_max);
            return;
        }

        let step = if filtered.len() <= n {
            1
        } else {
            filtered.len() / n
        };

        let sampled: Vec<&IsoLayerResult> =
            filtered.iter().step_by(step).take(n).copied().collect();

        let rows: Vec<EaRow> = sampled
            .clone()
            .into_iter()
            .map(|layer| EaRow {
                eta: format!("{:.4}", layer.eta),
                ea: format!("{:.2}", layer.ea / 1000.0),
                r2: format!("{:.4}", layer.regression.r2),
            })
            .collect();

        let table = Table::new(rows).with(Style::rounded()).to_string();
        println!("{}", table);
        if let Some(threshold) = threshold {
            // Assert R² > threshold for each sampled layer
            for layer in &sampled {
                assert!(
                    layer.regression.r2 > threshold,
                    "R² = {} is not greater than {} for η = {}",
                    layer.regression.r2,
                    threshold,
                    layer.eta
                );
            }
        }
    }
}

//========================================================================================
//               ADJUSTING PROBLEM
/// Selects whether the y-axis of the regression uses natural log or log base-10.
#[derive(Clone, Debug)]
pub enum LogKind {
    /// Natural logarithm (used by KAS and Starink).
    Ln,
    /// Base-10 logarithm (used by the original Ozawa formulation).
    Log10,
}

/// Configuration that fully parameterises one integral isoconversional method.
///
/// The three fields encode the method-specific constants so that a single
/// `solve` loop handles KAS, OFW, and Starink without branching.
#[derive(Clone, Debug)]
pub struct IntegralIsoConfig {
    /// Power of T in the denominator: `β / T^exponent`.
    /// KAS = 2.0, Starink = 1.92, OFW = 0.0.
    pub exponent: f64,
    /// Whether to take `ln` or `log10` of `β / T^exponent`.
    pub log_kind: LogKind,
    /// Multiplier that converts the regression slope to Eα (J/mol):
    /// `Ea = slope * slope_factor`.
    pub slope_factor: f64,
}

impl Default for IntegralIsoConfig {
    fn default() -> Self {
        Self {
            exponent: 0.0,
            log_kind: LogKind::Ln,
            slope_factor: 0.0,
        }
    }
}

/// Core solver shared by KAS, OFW, and Starink.
///
/// Parameterised entirely by [`IntegralIsoConfig`]; use the factory methods
/// [`Self::kas`], [`Self::ofw`], [`Self::starink`] to get a correctly
/// configured instance.
pub struct IntegralIsoconversionalSolver {
    pub config: IntegralIsoConfig,
    pub method: IntIsoconversionalMethod,
}

impl IntegralIsoconversionalSolver {
    /// Constructs a solver for the given method variant using the appropriate
    /// pre-set [`IntegralIsoConfig`].
    pub fn new(method: IntIsoconversionalMethod) -> Self {
        match method {
            IntIsoconversionalMethod::KAS => Self::kas(),

            IntIsoconversionalMethod::Starink => Self::starink(),

            IntIsoconversionalMethod::OFW => Self::ofw(),
        }
    }

    fn method_name(&self) -> &'static str {
        match self.method {
            IntIsoconversionalMethod::KAS => "KAS",
            IntIsoconversionalMethod::OFW => "Ozawa–Flynn–Wall",
            IntIsoconversionalMethod::Starink => "Starink",
        }
    }

    /// Runs the isoconversional regression over all α layers in `grid`.
    ///
    /// For each layer j the method:
    /// 1. Reads `1/T` (pre-cached in `grid.inv_temperature`) and `T` for each experiment.
    /// 2. Computes `y = log(β / T^exponent)` using the config.
    /// 3. Fits `y = slope * (1/T) + intercept` via OLS.
    /// 4. Converts the slope to Eα via `slope_factor`.
    pub fn solve(&self, grid: &ConversionGrid) -> Result<IsoconversionalResult, TGADomainError> {
        let n_eta = grid.eta.len();

        let mut layers = Vec::with_capacity(n_eta);

        for j in 0..n_eta {
            let eta = grid.eta[j];

            let inv_t = grid.inv_temperature.column(j);

            let temp = grid.temperature.column(j);

            // collect heating rates with validation instead of unwrap to avoid panics
            let mut beta: Vec<f64> = Vec::with_capacity(grid.meta.len());
            for m in &grid.meta {
                match m.heating_rate {
                    Some(b) if b.is_finite() && b > 0.0 => beta.push(b),
                    _ => {
                        return Err(TGADomainError::InvalidOperation(
                            "Heating rate missing or invalid".into(),
                        ));
                    }
                }
            }

            let x = inv_t.to_owned();

            let mut y = Array1::<f64>::zeros(beta.len());

            for (i, b) in beta.iter().enumerate() {
                let t = temp[i];
                if !t.is_finite() || t <= 0.0 {
                    return Err(TGADomainError::InvalidOperation(
                        "Non-physical temperature encountered in grid".into(),
                    ));
                }

                let val = b / t.powf(self.config.exponent);

                let logv = match self.config.log_kind {
                    LogKind::Ln => val.ln(),

                    LogKind::Log10 => val.log10(),
                };

                if !logv.is_finite() {
                    return Err(TGADomainError::InvalidOperation(
                        "Log argument non-finite in isoconversional solver".into(),
                    ));
                }

                y[i] = logv;
            }

            // y = ln(beta/T^b)
            // x = 1/T
            // y = k*x + b
            // k = -C*E/R
            // E = k*(-R/C)
            let reg = linear_regression(&x, &y);

            let ea = reg.slope * self.config.slope_factor;

            layers.push(IsoLayerResult {
                eta,

                ea,
                k: None,

                regression: reg,
            });
        }

        Ok(IsoconversionalResult {
            method: self.method_name(),

            layers,
        })
    }

    //Готовые фабрики методов
    /// Factory: returns a solver configured for the KAS method.
    /// Regression: `ln(β/T²)` vs `1/T`; slope factor = `-R`.
    pub fn kas() -> Self {
        Self {
            config: IntegralIsoConfig {
                exponent: 2.0,

                log_kind: LogKind::Ln,

                slope_factor: -R,
            },
            method: IntIsoconversionalMethod::KAS,
        }
    }

    /// Factory: returns a solver configured for the Starink method.
    /// Regression: `ln(β/T^1.92)` vs `1/T`; slope factor = `-R/1.0008`.
    pub fn starink() -> Self {
        Self {
            config: IntegralIsoConfig {
                exponent: 1.92,

                log_kind: LogKind::Ln,

                slope_factor: -R / 1.0008,
            },
            method: IntIsoconversionalMethod::Starink,
        }
    }

    /// Factory: returns a solver configured for the OFW method.
    /// Regression: `ln(β)` vs `1/T`; slope factor = `-R/1.052`.
    pub fn ofw() -> Self {
        Self {
            config: IntegralIsoConfig {
                exponent: 0.0,

                log_kind: LogKind::Ln,

                slope_factor: -R / 1.052,
            },
            method: IntIsoconversionalMethod::OFW,
        }
    }
}

/// `KineticMethod` entry point for the Ozawa–Flynn–Wall method.
///
/// Builds a `ConversionGrid` from the view and delegates to
/// `IntegralIsoconversionalSolver::ofw()`.
pub struct OFW {}

impl KineticMethod for OFW {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        "Ozawa–Flynn–Wall"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        // self.check_input(data)?;
        let now = Instant::now();
        let grid = ConversionGridBuilder::new()
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(data)?;
        println!("grid built in {} ms", now.elapsed().as_millis());
        let now = Instant::now();
        let solver = IntegralIsoconversionalSolver::new(IntIsoconversionalMethod::OFW);

        let out = solver.solve(&grid);
        println!("solver completed in {} ms", now.elapsed().as_millis());
        out
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,

            needs_conversion: true,

            needs_conversion_rate: false,

            needs_temperature: true,

            needs_heating_rate: true,
        }
    }
    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        vec![
            ColumnNature::Conversion,
            ColumnNature::Temperature,
            ColumnNature::Time,
            ColumnNature::ConversionRate,
        ]
    }
}

/// `KineticMethod` entry point for the Kissinger–Akahira–Sunose method.
///
/// Builds a `ConversionGrid` from the view and delegates to
/// `IntegralIsoconversionalSolver::kas()`.
pub struct KAS {}

impl KineticMethod for KAS {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        "KAS"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        // self.check_input(data)?;

        let grid = ConversionGridBuilder::new()
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(data)?;

        let solver = IntegralIsoconversionalSolver::new(IntIsoconversionalMethod::KAS);

        solver.solve(&grid)
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,

            needs_conversion: true,

            needs_conversion_rate: false,

            needs_temperature: true,

            needs_heating_rate: true,
        }
    }
    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        vec![
            ColumnNature::Conversion,
            ColumnNature::Temperature,
            ColumnNature::Time,
            ColumnNature::ConversionRate,
        ]
    }
}

pub struct Starink {}

impl KineticMethod for Starink {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        "Ozawa–Flynn–Wall"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        // self.check_input(data)?;
        let now = Instant::now();
        let grid = ConversionGridBuilder::new()
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(data)?;
        println!("grid built in {} ms", now.elapsed().as_millis());
        let now = Instant::now();
        let solver = IntegralIsoconversionalSolver::new(IntIsoconversionalMethod::Starink);

        let out = solver.solve(&grid);
        println!("solver completed in {} ms", now.elapsed().as_millis());
        out
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,

            needs_conversion: true,

            needs_conversion_rate: false,

            needs_temperature: true,

            needs_heating_rate: true,
        }
    }
    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        vec![
            ColumnNature::Conversion,
            ColumnNature::Temperature,
            ColumnNature::Time,
            ColumnNature::ConversionRate,
        ]
    }
}
/// `KineticMethod` entry point for the Vyazovkin non-linear method.
///
/// Builds a `ConversionGrid` **with the dt matrix** and delegates to
/// `VyazovkinSolver::default()`.
pub struct VyazovkinMethod {}

impl KineticMethod for VyazovkinMethod {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        "Vyazovkin"
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,

            needs_conversion: true,

            needs_conversion_rate: false,

            needs_temperature: true,

            needs_heating_rate: true,
        }
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        let grid = ConversionGridBuilder::default()
            .with_dt_matrix()
            .build_nonisothermal(data)?;

        VyazovkinSolver::default().solve(&grid)
    }

    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        vec![
            ColumnNature::Conversion,
            ColumnNature::Temperature,
            ColumnNature::Time,
            ColumnNature::ConversionRate,
        ]
    }
}
