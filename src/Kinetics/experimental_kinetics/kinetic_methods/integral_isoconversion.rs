//! integral isoconversional kinetic methods:
//! Ozawa-Flinn-Wall
//! KAS
//! Starink
//!
//! PIPELINE:
//! UnitedDataset
//!      │
//!      ▼
//!ConversionGridBuilder
//!      │
//!      ▼
//!ConversionGrid (ndarray)
//!      │
//!      ▼
//!IntegralIsoconversionalSolver
//!      │
//!      ▼
//!IsoconversionalResult
//!
use super::super::kinetic_methods::{
    ConversionGrid, ConversionGridBuilder, GridInterpolation, KineticDataView, KineticMethod,
    KineticRequirements, TGADomainError,
};
use super::Vyazovkin::VyazovkinSolver;
use super::kinetic_regression::{LinearRegressionResult, linear_regression};
use ndarray::Array1;
use polars::prelude::*;
use std::time::Instant;
pub const R: f64 = 8.314;
//====================================================================================================================
//          METHODS
//======================================================================================================================
#[derive(Clone, Copy, Debug)]
pub enum IntIsoconversionalMethod {
    KAS,
    OFW,
    Starink,
}

//========================================================================================================================
//HANDLING RESULT
//=====================================================================================================================
// "результат слоя α".
#[derive(Clone, Debug)]
pub struct IsoLayerResult {
    pub eta: f64,

    pub ea: f64,
    pub k: Option<f64>,
    pub regression: LinearRegressionResult,
}

// Итоговый результат метода
#[derive(Clone, Debug)]
pub struct IsoconversionalResult {
    pub method: &'static str,

    pub layers: Vec<IsoLayerResult>,
}

impl IsoconversionalResult {
    /// to polars eager DataFrame
    pub fn to_dataframe(&self) -> DataFrame {
        let eta: Vec<f64> = self.layers.iter().map(|l| l.eta).collect();

        let ea: Vec<f64> = self.layers.iter().map(|l| l.ea).collect();

        let r2: Vec<f64> = self.layers.iter().map(|l| l.regression.r2).collect();

        DataFrame::new_infer_height(vec![
            Series::new("alpha".into(), eta).into(),
            Series::new("Ea".into(), ea).into(),
            Series::new("R2".into(), r2).into(),
        ])
        .unwrap()
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
#[derive(Clone, Debug)]
pub enum LogKind {
    Ln,
    Log10,
}
#[derive(Clone, Debug)]
pub struct IntegralIsoConfig {
    pub exponent: f64,

    pub log_kind: LogKind,

    /// slope -> Ea conversion factor
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

pub struct IntegralIsoconversionalSolver {
    pub config: IntegralIsoConfig,
    pub method: IntIsoconversionalMethod,
}

impl IntegralIsoconversionalSolver {
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
    //KAS
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

    //Starink
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

    //OFW
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
            .build(data)?;
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
}

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
            .build(data)?;

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
}

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
            .build(data)?;

        VyazovkinSolver::default().solve(&grid)
    }
}
