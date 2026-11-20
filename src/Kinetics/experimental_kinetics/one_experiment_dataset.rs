//! # TGA data preprocessing pipeline
//!
//! Recommended order:
//! 1. Read raw data (CSV / TXT)
//! 2. Bind semantic columns (time, mass, temperature)
//! 3. Calibrate raw signals if needed
//! 4. Dimensionless transforms
//! 5. Smoothing / filtering
//! 6. Derivatives
//! 7. Trim invalid edges
//!
//! ## Русский ориентир
//! 1. Читаем сырые данные
//! 2. Привязываем смысловые колонки: время, масса, температура
//! 3. При необходимости калибруем сигналы
//! 4. Строим безразмерные величины
//! 5. Сглаживаем и фильтруем
//! 6. Считаем производные
//! 7. Обрезаем некорректные края
/// # Recommended TGA data processing pipeline
///
/// KiThe follows an explicit, step-by-step processing model.
/// KiThe использует явную пошаговую модель обработки данных.
/// No operation implicitly removes rows or modifies unrelated columns.
/// Ни одна операция не удаляет строки и не меняет посторонние колонки неявно.
///
/// ## General principles
///
/// 1. **All column operations preserve table length**
///    - rolling mean, differentiation, filtering produce `null`
///    - row removal is always explicit
///
/// 2. **`null` values are valid intermediate state**
///    - smoothing and differentiation naturally introduce `null`
///    - users must explicitly remove invalid edges
///
/// 3. **Row removal is explicit**
///    - use `trim_edges(left, right)` or `filter_rows(...)`
///    - no operation removes rows automatically
///
/// 4. **Raw experimental columns should be preserved**
///    - smoothing/filtering is usually applied to derived columns
///
/// ## Typical workflow
///
/// ```text
/// raw data (time, temperature, mass)
///     ↓
/// optional calibration / scaling
///     ↓
/// optional smoothing (rolling mean, Hampel, etc.)
///     ↓
/// trim_edges (to remove null introduced by smoothing)
///     ↓
/// relative mass α / conversion η
///     ↓
/// differentiation (rates)
///     ↓
/// trim_edges (to remove null from finite differences)
/// ```
///
/// ## Notes
///
/// - Hampel filter assumes the input column has no null values.
/// - It is recommended to call `trim_edges` before Hampel filtering.
/// - Smoothing and filtering do not modify column metadata (units, origin).
///
use polars::error::{PolarsError, PolarsResult};

use polars::prelude::LazyFrame;
use polars::prelude::*;
use polars::series::Series;

use crate::Kinetics::experimental_kinetics::column_provenance::{
    ColumnProvenance, ColumnTransformKind, QualityReference, TransformQuality,
};
use log::info;
use std::collections::{HashMap, HashSet, VecDeque};
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;
//===============================================================
// One datafame plot data
#[derive(Clone, Debug)]
pub struct OneFramePlot {
    pub x: Option<String>,
    pub y: Option<String>,
}

//===========================================

pub struct UnaryOp {
    pub func: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    pub output_unit: Unit,
    pub domain_check: Option<fn(&TGADataset, &str) -> Result<(), TGADomainError>>,
}
//=================================================================
// ERROR HANDLING
//=================================================================
#[derive(Debug, Clone)]
pub enum TGADomainError {
    PolarsError(PolarsError),
    InvalidFilterParameters {
        filter: String,
        details: String,
    },
    NonNumericColumn {
        column: String,
        actual_type: String,
    },
    InvalidColumnValues {
        column: String,
        details: String,
        indices: Vec<usize>,
    },
    IncompatibleUnits {
        operation: String,
        expected: String,
        actual: Unit,
    },
    NoConvergence {
        method: String,
        details: String,
    },
    InsufficientData {
        available: usize,
        required: usize,
    },
    DegenerateRegression {
        details: String,
    },
    FittingInputLengthMismatch {
        x_len: usize,
        y_len: usize,
    },
    InvalidFittingConfiguration {
        details: String,
    },
    ConversionRangeNotCovered {
        experiment_id: String,
        requested_min: f64,
        requested_max: f64,
        available_min: f64,
        available_max: f64,
    },
    OperationContext {
        experiment_id: Option<String>,
        column: Option<String>,
        operation: String,
        source: Box<TGADomainError>,
    },
    EmptyColumn(String),
    NonPositiveValues {
        column: String,
    },
    ColumnNotFound(String),
    AmbiguousColumnNature {
        nature: ColumnNature,
        columns: Vec<String>,
    },
    ExperimentDataLengthMismatch {
        experiment_id: String,
        column: String,
        expected: usize,
        actual: usize,
    },
    DatasetColumnLengthMismatch {
        column: String,
        expected: usize,
        actual: usize,
    },
    MassNotBound,
    TimeNotBound,
    NoDimensionless,
    AlphaNotFound,
    TemperatureNotBound,
    DtDtNotFound,
    DmDtNotFound,
    DalphaDtNotFound,

    EtaNotFound,
    DetaDtNotFound,

    UnsupportedRateUnit,
    InvalidReferenceMass,
    InvalidUnitForConversion,
    InvalidUnitForBinding {
        field: &'static str,
        unit: Unit,
    },
    ColumnAlreadyExists(String),
    DuplicateExperimentId(String),
    ExperimentNotFound(String),
    NotImplemented,
    InvalidOperation(String),
    OutOfRange,
    TimeNonMonotonic(String),
    InvalidConversionRange,
}

impl TGADomainError {
    pub fn with_context(
        self,
        operation: impl Into<String>,
        experiment_id: Option<&str>,
        column: Option<&str>,
    ) -> Self {
        Self::OperationContext {
            experiment_id: experiment_id.map(str::to_owned),
            column: column.map(str::to_owned),
            operation: operation.into(),
            source: Box::new(self),
        }
    }
}

impl std::fmt::Display for TGADomainError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OperationContext {
                experiment_id,
                column,
                operation,
                source,
            } => {
                write!(f, "{operation} failed")?;
                if let Some(experiment_id) = experiment_id {
                    write!(f, " for experiment '{experiment_id}'")?;
                }
                if let Some(column) = column {
                    write!(f, " on column '{column}'")?;
                }
                write!(f, ": {source}")
            }
            Self::InvalidFilterParameters { filter, details } => {
                write!(f, "Invalid parameters for {filter}: {details}")
            }
            Self::InvalidColumnValues {
                column,
                details,
                indices,
            } => write!(
                f,
                "Invalid values in column '{column}': {details} (rows: {indices:?})"
            ),
            Self::NoConvergence { method, details } => {
                write!(f, "{method} did not produce a valid result: {details}")
            }
            Self::InsufficientData {
                available,
                required,
            } => write!(
                f,
                "Insufficient observations: {available} available, {required} required"
            ),
            Self::DegenerateRegression { details } => {
                write!(f, "Regression is degenerate: {details}")
            }
            Self::FittingInputLengthMismatch { x_len, y_len } => write!(
                f,
                "Fitting X/Y columns have different lengths ({x_len} vs {y_len})"
            ),
            Self::InvalidFittingConfiguration { details } => {
                write!(f, "Invalid fitting configuration: {details}")
            }
            Self::ConversionRangeNotCovered {
                experiment_id,
                requested_min,
                requested_max,
                available_min,
                available_max,
            } => write!(
                f,
                "Experiment '{experiment_id}' covers conversion [{available_min}, {available_max}], but requested [{requested_min}, {requested_max}]"
            ),
            Self::DatasetColumnLengthMismatch {
                column,
                expected,
                actual,
            } => write!(
                f,
                "Column '{column}' has {actual} values, but the dataset has {expected} rows"
            ),
            Self::InvalidOperation(message) => write!(f, "{message}"),
            _ => write!(f, "{self:?}"),
        }
    }
}

impl std::error::Error for TGADomainError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::PolarsError(source) => Some(source),
            Self::OperationContext { source, .. } => Some(source.as_ref()),
            _ => None,
        }
    }
}

/// Detailed diagnostics for the active time column.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeOrderViolation {
    pub index: usize,
    pub previous: f64,
    pub current: f64,
}

/// Structured report describing why the active time column is or is not monotonic.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeMonotonicityReport {
    pub time_column: String,
    pub len: usize,
    pub null_indices: Vec<usize>,
    pub nan_indices: Vec<usize>,
    pub infinite_indices: Vec<usize>,
    pub duplicate_pairs: Vec<TimeOrderViolation>,
    pub decreasing_pairs: Vec<TimeOrderViolation>,
    pub is_strictly_increasing: bool,
    pub sortable_without_data_loss: bool,
}

impl TimeMonotonicityReport {
    pub fn has_invalid_values(&self) -> bool {
        !self.null_indices.is_empty()
            || !self.nan_indices.is_empty()
            || !self.infinite_indices.is_empty()
    }
}

impl From<PolarsError> for TGADomainError {
    fn from(err: PolarsError) -> Self {
        TGADomainError::PolarsError(err)
    }
}
//=================================================================
//UNITS
//=================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unit {
    Second,
    Minute,
    Hour,

    Kelvin,
    Celsius,

    MilliVolt,
    Milligram,
    MilligramPerSecond,
    MilligramPerMinute,
    MilligramPerHour,
    KelvinPerSecond,
    KelvinPerMinute,
    KelvinPerHour,
    CelsiusPerSecond,
    CelsiusPerMinute,
    CelsiusPerHour,
    PerSecond,
    PerMinute,
    PerHour,
    Gram,

    Dimensionless,
    Unknown,
}

impl Unit {
    fn is_valid_time_binding_unit(self) -> bool {
        matches!(self, Unit::Second | Unit::Minute | Unit::Hour)
    }

    fn is_valid_temperature_binding_unit(self) -> bool {
        matches!(self, Unit::Kelvin | Unit::Celsius)
    }

    fn is_valid_mass_binding_unit(self) -> bool {
        matches!(self, Unit::MilliVolt | Unit::Milligram | Unit::Gram)
    }

    /// Normalize a rate unit to the time base used by the time column.
    ///
    /// KiThe currently supports time in seconds or hours. The numeric derivative
    /// is computed against the stored time values, so this helper only adjusts the
    /// semantic unit label to match the time axis.
    pub fn for_time_base(self, time_unit: Unit) -> Result<Self, TGADomainError> {
        match time_unit {
            Unit::Second => match self {
                Unit::Milligram
                | Unit::MilligramPerSecond
                | Unit::MilligramPerMinute
                | Unit::MilligramPerHour => Ok(Unit::MilligramPerSecond),
                Unit::Kelvin
                | Unit::KelvinPerSecond
                | Unit::KelvinPerMinute
                | Unit::KelvinPerHour => Ok(Unit::KelvinPerSecond),
                Unit::Celsius
                | Unit::CelsiusPerSecond
                | Unit::CelsiusPerMinute
                | Unit::CelsiusPerHour => Ok(Unit::CelsiusPerSecond),
                Unit::Dimensionless | Unit::PerSecond | Unit::PerMinute | Unit::PerHour => {
                    Ok(Unit::PerSecond)
                }
                _ => Err(TGADomainError::UnsupportedRateUnit),
            },
            Unit::Minute => match self {
                Unit::Milligram
                | Unit::MilligramPerSecond
                | Unit::MilligramPerMinute
                | Unit::MilligramPerHour => Ok(Unit::MilligramPerMinute),
                Unit::Kelvin
                | Unit::KelvinPerSecond
                | Unit::KelvinPerMinute
                | Unit::KelvinPerHour => Ok(Unit::KelvinPerMinute),
                Unit::Celsius
                | Unit::CelsiusPerSecond
                | Unit::CelsiusPerMinute
                | Unit::CelsiusPerHour => Ok(Unit::CelsiusPerMinute),
                Unit::Dimensionless | Unit::PerSecond | Unit::PerMinute | Unit::PerHour => {
                    Ok(Unit::PerMinute)
                }
                _ => Err(TGADomainError::UnsupportedRateUnit),
            },
            Unit::Hour => match self {
                Unit::Milligram
                | Unit::MilligramPerSecond
                | Unit::MilligramPerMinute
                | Unit::MilligramPerHour => Ok(Unit::MilligramPerHour),
                Unit::Kelvin
                | Unit::KelvinPerSecond
                | Unit::KelvinPerMinute
                | Unit::KelvinPerHour => Ok(Unit::KelvinPerHour),
                Unit::Celsius
                | Unit::CelsiusPerSecond
                | Unit::CelsiusPerMinute
                | Unit::CelsiusPerHour => Ok(Unit::CelsiusPerHour),
                Unit::Dimensionless | Unit::PerSecond | Unit::PerMinute | Unit::PerHour => {
                    Ok(Unit::PerHour)
                }
                _ => Err(TGADomainError::UnsupportedRateUnit),
            },
            _ => Err(TGADomainError::UnsupportedRateUnit),
        }
    }

    pub fn parse(s: &str) -> Result<Self, PolarsError> {
        match s {
            "s" => Ok(Unit::Second),
            "min" => Ok(Unit::Minute),
            "h" => Ok(Unit::Hour),
            "C" => Ok(Unit::Celsius),
            "K" => Ok(Unit::Kelvin),
            "mV" => Ok(Unit::MilliVolt),
            "mg" => Ok(Unit::Milligram),
            "mg/s" => Ok(Unit::MilligramPerSecond),
            "mg/min" => Ok(Unit::MilligramPerMinute),
            "mg/h" => Ok(Unit::MilligramPerHour),
            "K/s" => Ok(Unit::KelvinPerSecond),
            "K/min" => Ok(Unit::KelvinPerMinute),
            "K/h" => Ok(Unit::KelvinPerHour),
            "C/s" => Ok(Unit::CelsiusPerSecond),
            "C/min" => Ok(Unit::CelsiusPerMinute),
            "C/h" => Ok(Unit::CelsiusPerHour),
            "/s" => Ok(Unit::PerSecond),
            "/min" => Ok(Unit::PerMinute),
            "1/s" => Ok(Unit::PerSecond),
            "1/min" => Ok(Unit::PerMinute),
            "/h" => Ok(Unit::PerHour),
            "1/h" => Ok(Unit::PerHour),
            "g" => Ok(Unit::Gram),
            "Dimensionless" => Ok(Unit::Dimensionless),
            "Unknown" => Ok(Unit::Unknown),
            "-" => Ok(Unit::Dimensionless),

            _ => Err(PolarsError::ComputeError(
                format!("Unknown unit: {}", s).into(),
            )),
        }
    }
    pub fn convert_to_string(&self) -> String {
        match self {
            Unit::Second => "s".to_string(),
            Unit::Minute => "min".to_string(),
            Unit::Hour => "h".to_string(),
            Unit::Kelvin => "K".to_string(),
            Unit::Celsius => "C".to_string(),
            Unit::MilliVolt => "mV".to_string(),
            Unit::Milligram => "mg".to_string(),
            Unit::MilligramPerSecond => "mg/s".to_string(),
            Unit::MilligramPerMinute => "mg/min".to_string(),
            Unit::MilligramPerHour => "mg/h".to_string(),
            Unit::KelvinPerSecond => "K/s".to_string(),
            Unit::KelvinPerMinute => "K/min".to_string(),
            Unit::KelvinPerHour => "K/h".to_string(),
            Unit::CelsiusPerSecond => "C/s".to_string(),
            Unit::CelsiusPerMinute => "C/min".to_string(),
            Unit::CelsiusPerHour => "C/h".to_string(),
            Unit::PerSecond => "/s".to_string(),
            Unit::PerMinute => "/min".to_string(),
            Unit::PerHour => "/h".to_string(),
            Unit::Gram => "g".to_string(),
            Unit::Dimensionless => "Dimensionless".to_string(),
            Unit::Unknown => "Unknown".to_string(),
        }
    }
}

impl std::fmt::Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unit::Second => write!(f, "s"),
            Unit::Minute => write!(f, "min"),
            Unit::Hour => write!(f, "h"),
            Unit::Kelvin => write!(f, "K"),
            Unit::Celsius => write!(f, "C"),
            Unit::MilliVolt => write!(f, "mV"),
            Unit::Milligram => write!(f, "mg"),
            Unit::MilligramPerSecond => write!(f, "mg/s"),
            Unit::MilligramPerMinute => write!(f, "mg/min"),
            Unit::MilligramPerHour => write!(f, "mg/h"),
            Unit::KelvinPerSecond => write!(f, "K/s"),
            Unit::KelvinPerMinute => write!(f, "K/min"),
            Unit::KelvinPerHour => write!(f, "K/h"),
            Unit::CelsiusPerSecond => write!(f, "C/s"),
            Unit::CelsiusPerMinute => write!(f, "C/min"),
            Unit::CelsiusPerHour => write!(f, "C/h"),
            Unit::PerSecond => write!(f, "/s"),
            Unit::PerMinute => write!(f, "/min"),
            Unit::PerHour => write!(f, "/h"),
            Unit::Gram => write!(f, "g"),
            Unit::Dimensionless => write!(f, "Dimensionless"),
            Unit::Unknown => write!(f, "Unknown"),
        }
    }
}

//=====================================================================
//              COLUMN METADATA
//====================================================================
#[derive(Debug, Clone)]
pub enum ColumnRole {
    Time,
    Temperature,
    Mass,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum ColumnNature {
    Time,
    Temperature,
    Mass,
    // alpha = relative mass m/m0
    Conversion,
    // eta = 1 - alpha
    DimensionlessMass,
    MassRate,
    ConversionRate,
    TemperatureRate,
    DimensionlessMassRate,
    Unknown,
    ActivationEnergy,
    PredexFactor,
    R2,
    FittedCurve,
}

/// Controls when a numerical derivative is withheld across a large time gap.
///
/// Derivatives use the bound time column's stored values, so `MaxInterval` is
/// expressed in that column's unit (for example, seconds or hours).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TimeGapPolicy {
    /// Mark a derivative as null when either adjacent interval exceeds this
    /// multiple of the median positive sampling interval.
    MedianStepMultiple(f64),
    /// Mark a derivative as null when either adjacent interval exceeds this
    /// absolute duration in the bound time unit.
    MaxInterval(f64),
    /// Do not reject large intervals; still reject non-positive or non-finite
    /// intervals.
    Allow,
}

impl Default for TimeGapPolicy {
    fn default() -> Self {
        Self::MedianStepMultiple(5.0)
    }
}
/// Columns Origin
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum ColumnOrigin {
    Raw,
    PolarsDerived,
    NumericDerived,
    Imported,
}
/// Column  metadata
#[derive(Debug, Clone)]
pub struct ColumnMeta {
    pub name: String,
    pub unit: Unit,
    pub origin: ColumnOrigin,
    pub nature: ColumnNature,
    pub provenance: ColumnProvenance,
}

impl ColumnMeta {
    pub fn new(
        name: impl Into<String>,
        unit: Unit,
        origin: ColumnOrigin,
        nature: ColumnNature,
        provenance: ColumnProvenance,
    ) -> Self {
        Self {
            name: name.into(),
            unit,
            origin,
            nature,
            provenance,
        }
    }

    pub fn raw(name: impl Into<String>, unit: Unit, nature: ColumnNature) -> Self {
        let name = name.into();
        Self::new(
            name.clone(),
            unit,
            ColumnOrigin::Raw,
            nature,
            ColumnProvenance::raw(name),
        )
    }

    pub fn set_nature(&mut self, nature: ColumnNature) {
        self.nature = nature;
    }

    pub fn get_nature(&self) -> ColumnNature {
        self.nature
    }
}

/// Scheme of TGA dataset
#[derive(Debug, Clone)]
pub struct TGASchema {
    pub columns: HashMap<String, ColumnMeta>,
    pub time: Option<String>,
    pub temperature: Option<String>,
    pub mass: Option<String>,
    pub alpha: Option<String>,
    pub dm_dt: Option<String>,
    pub eta: Option<String>,
    pub deta_dt: Option<String>,
    pub dalpha_dt: Option<String>,
    pub dT_dt: Option<String>,
    pub E: Option<String>,
    pub R2: Option<String>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
#[serde(default)]
struct CsvSchemaBindings {
    time: Option<String>,
    temperature: Option<String>,
    mass: Option<String>,
    alpha: Option<String>,
    dm_dt: Option<String>,
    eta: Option<String>,
    deta_dt: Option<String>,
    dalpha_dt: Option<String>,
    dT_dt: Option<String>,
    E: Option<String>,
    R2: Option<String>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
struct CsvColumnMetadata {
    name: String,
    nature: ColumnNature,
    origin: ColumnOrigin,
}

impl CsvColumnMetadata {
    fn from_schema(schema: &TGASchema, names: &[String]) -> Vec<Self> {
        names
            .iter()
            .map(|name| {
                let meta = schema.columns.get(name);
                Self {
                    name: name.clone(),
                    nature: meta
                        .map(|meta| meta.nature)
                        .unwrap_or(ColumnNature::Unknown),
                    origin: meta
                        .map(|meta| meta.origin)
                        .unwrap_or(ColumnOrigin::Imported),
                }
            })
            .collect()
    }

    fn apply_all(metadata: Vec<Self>, schema: &mut TGASchema) -> PolarsResult<()> {
        let mut seen = HashSet::with_capacity(metadata.len());
        for item in metadata {
            if !seen.insert(item.name.clone()) {
                return Err(PolarsError::ComputeError(
                    format!("CSV metadata repeats column '{}'", item.name).into(),
                ));
            }
            let meta = schema.columns.get_mut(&item.name).ok_or_else(|| {
                PolarsError::ComputeError(
                    format!("CSV metadata refers to missing column '{}'", item.name).into(),
                )
            })?;
            meta.nature = item.nature;
            meta.origin = item.origin;
        }
        Ok(())
    }
}

impl CsvSchemaBindings {
    fn from_schema(schema: &TGASchema, selected: Option<&[&str]>) -> Self {
        let keep = |name: &Option<String>| {
            name.as_ref()
                .filter(|name| selected.is_none_or(|selected| selected.contains(&name.as_str())))
                .cloned()
        };
        Self {
            time: keep(&schema.time),
            temperature: keep(&schema.temperature),
            mass: keep(&schema.mass),
            alpha: keep(&schema.alpha),
            dm_dt: keep(&schema.dm_dt),
            eta: keep(&schema.eta),
            deta_dt: keep(&schema.deta_dt),
            dalpha_dt: keep(&schema.dalpha_dt),
            dT_dt: keep(&schema.dT_dt),
            E: keep(&schema.E),
            R2: keep(&schema.R2),
        }
    }

    fn apply(self, schema: &mut TGASchema) -> PolarsResult<()> {
        for (field, name) in [
            ("time", &self.time),
            ("temperature", &self.temperature),
            ("mass", &self.mass),
            ("alpha", &self.alpha),
            ("dm_dt", &self.dm_dt),
            ("eta", &self.eta),
            ("deta_dt", &self.deta_dt),
            ("dalpha_dt", &self.dalpha_dt),
            ("dT_dt", &self.dT_dt),
            ("E", &self.E),
            ("R2", &self.R2),
        ] {
            if let Some(name) = name {
                let meta = schema.columns.get_mut(name).ok_or_else(|| {
                    PolarsError::ComputeError(
                        format!("CSV binding '{field}' refers to missing column '{name}'").into(),
                    )
                })?;
                meta.nature = match field {
                    "time" => ColumnNature::Time,
                    "temperature" => ColumnNature::Temperature,
                    "mass" => ColumnNature::Mass,
                    "alpha" => ColumnNature::DimensionlessMass,
                    "dm_dt" => ColumnNature::MassRate,
                    "eta" => ColumnNature::Conversion,
                    "deta_dt" => ColumnNature::ConversionRate,
                    "dalpha_dt" => ColumnNature::DimensionlessMassRate,
                    "dT_dt" => ColumnNature::TemperatureRate,
                    "E" => ColumnNature::ActivationEnergy,
                    "R2" => ColumnNature::R2,
                    _ => unreachable!(),
                };
            }
        }
        schema.time = self.time;
        schema.temperature = self.temperature;
        schema.mass = self.mass;
        schema.alpha = self.alpha;
        schema.dm_dt = self.dm_dt;
        schema.eta = self.eta;
        schema.deta_dt = self.deta_dt;
        schema.dalpha_dt = self.dalpha_dt;
        schema.dT_dt = self.dT_dt;
        schema.E = self.E;
        schema.R2 = self.R2;
        schema.refresh_active_column_natures()?;
        Ok(())
    }
}

impl Default for TGASchema {
    fn default() -> Self {
        Self {
            columns: HashMap::new(),
            time: None,
            temperature: None,
            mass: None,
            alpha: None,
            dm_dt: None,
            eta: None,
            deta_dt: None,
            dalpha_dt: None,
            dT_dt: None,
            E: None,
            R2: None,
        }
    }
}

impl TGASchema {
    fn get(&self, name: &str) -> Option<&ColumnMeta> {
        self.columns.get(name)
    }

    fn bind_primary_role(&mut self, role: &ColumnRole, name: &str) {
        for binding in [
            &mut self.time,
            &mut self.temperature,
            &mut self.mass,
            &mut self.alpha,
            &mut self.dm_dt,
            &mut self.eta,
            &mut self.deta_dt,
            &mut self.dalpha_dt,
            &mut self.dT_dt,
            &mut self.E,
            &mut self.R2,
        ] {
            if binding.as_deref() == Some(name) {
                *binding = None;
            }
        }
        match role {
            ColumnRole::Time => self.time = Some(name.to_string()),
            ColumnRole::Temperature => self.temperature = Some(name.to_string()),
            ColumnRole::Mass => self.mass = Some(name.to_string()),
        }
    }

    pub(crate) fn refresh_active_column_natures(&mut self) -> PolarsResult<()> {
        let bindings = [
            (&self.time, ColumnNature::Time, "time"),
            (&self.temperature, ColumnNature::Temperature, "temperature"),
            (&self.mass, ColumnNature::Mass, "mass"),
            (&self.alpha, ColumnNature::DimensionlessMass, "alpha"),
            (&self.dm_dt, ColumnNature::MassRate, "dm_dt"),
            (&self.eta, ColumnNature::Conversion, "eta"),
            (&self.deta_dt, ColumnNature::ConversionRate, "deta_dt"),
            (
                &self.dalpha_dt,
                ColumnNature::DimensionlessMassRate,
                "dalpha_dt",
            ),
            (&self.dT_dt, ColumnNature::TemperatureRate, "dT_dt"),
            (&self.E, ColumnNature::ActivationEnergy, "E"),
            (&self.R2, ColumnNature::R2, "R2"),
        ];
        let mut assigned = HashMap::<&str, (&str, ColumnNature)>::new();
        for (binding, nature, field) in bindings {
            let Some(name) = binding.as_deref() else {
                continue;
            };
            if !self.columns.contains_key(name) {
                return Err(PolarsError::ComputeError(
                    format!("Schema binding '{field}' refers to missing column '{name}'").into(),
                ));
            }
            if let Some((previous_field, previous_nature)) = assigned.insert(name, (field, nature))
                && previous_nature != nature
            {
                return Err(PolarsError::ComputeError(
                    format!(
                        "Column '{name}' cannot be active for both '{previous_field}' and '{field}'"
                    )
                    .into(),
                ));
            }
        }
        for (name, (_, nature)) in assigned {
            let metadata = self.columns.get_mut(name).ok_or_else(|| {
                PolarsError::ComputeError(
                    format!("Validated schema column '{name}' disappeared").into(),
                )
            })?;
            metadata.nature = nature;
        }
        Ok(())
    }

    fn nature_from_header(name: &str) -> ColumnNature {
        let label = name.split(['[', '(']).next().unwrap_or(name).trim();

        match label {
            "t" => return ColumnNature::Time,
            "T" => return ColumnNature::Temperature,
            "m" => return ColumnNature::Mass,
            _ => {}
        }

        match label.to_ascii_lowercase().as_str() {
            "time" => ColumnNature::Time,
            "temp" | "temperature" => ColumnNature::Temperature,
            "mass" | "weight" => ColumnNature::Mass,
            _ => ColumnNature::Unknown,
        }
    }

    // Infer a role only for unique header matches; ambiguous candidates stay
    // Unknown so later nature-based lookups cannot silently select one of them.
    fn infer_bindings_from_headers(&mut self, names: &[String]) {
        let mut time = Vec::new();
        let mut temperature = Vec::new();
        let mut mass = Vec::new();

        for name in names {
            let nature = Self::nature_from_header(name);
            match nature {
                ColumnNature::Time => time.push(name.clone()),
                ColumnNature::Temperature => temperature.push(name.clone()),
                ColumnNature::Mass => mass.push(name.clone()),
                _ => {}
            }
        }

        self.time = self.bind_unique_header_match(time, ColumnNature::Time);
        self.temperature = self.bind_unique_header_match(temperature, ColumnNature::Temperature);
        self.mass = self.bind_unique_header_match(mass, ColumnNature::Mass);
    }

    fn bind_unique_header_match(
        &mut self,
        mut matches: Vec<String>,
        nature: ColumnNature,
    ) -> Option<String> {
        if matches.len() == 1 {
            let name = matches.pop()?;
            if let Some(meta) = self.columns.get_mut(&name) {
                let has_compatible_unit = match nature {
                    ColumnNature::Time => meta.unit.is_valid_time_binding_unit(),
                    ColumnNature::Temperature => meta.unit.is_valid_temperature_binding_unit(),
                    ColumnNature::Mass => meta.unit.is_valid_mass_binding_unit(),
                    _ => false,
                };
                if !has_compatible_unit {
                    return None;
                }
                meta.nature = nature;
            } else {
                return None;
            }
            Some(name)
        } else {
            None
        }
    }

    fn from_columns(names: Vec<String>, units: Vec<Unit>) -> PolarsResult<Self> {
        if names.len() != units.len() {
            return Err(PolarsError::ComputeError(
                format!(
                    "Column metadata has {} names but {} units",
                    names.len(),
                    units.len()
                )
                .into(),
            ));
        }

        let mut columns = HashMap::new();
        for (name, unit) in names.iter().zip(units.iter()) {
            if name.trim().is_empty() {
                return Err(PolarsError::ComputeError(
                    "Column metadata contains an empty column name".into(),
                ));
            }
            if columns.contains_key(name) {
                return Err(PolarsError::ComputeError(
                    format!("Column metadata contains duplicate name '{name}'").into(),
                ));
            }
            columns.insert(
                name.clone(),
                ColumnMeta::new(
                    name.clone(),
                    *unit,
                    ColumnOrigin::Imported,
                    ColumnNature::Unknown,
                    ColumnProvenance::imported(name.clone(), name.clone()),
                ),
            );
        }
        let mut schema = Self {
            columns,
            ..Self::default()
        };
        schema.infer_bindings_from_headers(&names);
        Ok(schema)
    }

    /// Update schema when a column is renamed
    /// 1. Insert new column metadata with same unit/origin as old_name
    /// 2. Keep old column metadata unchanged
    /// 3. Update semantic fields (time, mass, etc.) if they reference old_name
    pub fn update_schema(&mut self, old_name: &str, new_name: &str) {
        // Clone the metadata of the old column if it exists
        if let Some(old_meta) = self.columns.get(old_name) {
            let mut provenance = old_meta.provenance.clone();
            provenance.rename_output(new_name);
            let new_meta = ColumnMeta::new(
                new_name.to_string(),
                old_meta.unit,
                old_meta.origin,
                old_meta.nature,
                provenance,
            );
            // Insert the new column metadata
            self.columns.insert(new_name.to_string(), new_meta);
        }

        for binding in [
            &mut self.time,
            &mut self.temperature,
            &mut self.mass,
            &mut self.alpha,
            &mut self.dm_dt,
            &mut self.eta,
            &mut self.deta_dt,
            &mut self.dalpha_dt,
            &mut self.dT_dt,
            &mut self.E,
            &mut self.R2,
        ] {
            if binding.as_deref() == Some(old_name) {
                *binding = Some(new_name.to_string());
            }
        }
    }

    /// Update the active binding after a numerical transform and mark its output
    /// as derived while retaining the source unit and semantic nature.
    pub fn update_schema_for_numeric_transform(&mut self, input: &str, output: &str) {
        self.update_schema(input, output);
        if let Some(metadata) = self.columns.get_mut(output) {
            metadata.origin = ColumnOrigin::NumericDerived;
        }
    }

    pub fn drop_schema_record(&mut self, name: &str) {
        if let Some(_) = self.columns.get(name) {
            self.columns.remove(name);
        }

        for binding in [
            &mut self.time,
            &mut self.temperature,
            &mut self.mass,
            &mut self.alpha,
            &mut self.dm_dt,
            &mut self.eta,
            &mut self.deta_dt,
            &mut self.dalpha_dt,
            &mut self.dT_dt,
            &mut self.E,
            &mut self.R2,
        ] {
            if binding.as_deref() == Some(name) {
                *binding = None;
            }
        }
    }
}
//=============================================================================
//             HISTORY OF OPERATIONS
//===============================================================================
#[derive(Clone, Debug)]
pub struct History {
    pub vector_of_operations: Vec<OperationRecord>,
    pub columns_has_changed: Vec<Option<String>>,
}

impl Default for History {
    fn default() -> Self {
        Self::new()
    }
}
impl History {
    pub fn new() -> Self {
        Self {
            vector_of_operations: Vec::new(),
            columns_has_changed: Vec::new(),
        }
    }
    pub fn columns_has_changed(&self) -> Vec<Option<String>> {
        self.columns_has_changed.clone()
    }
    pub fn take_column(&mut self, column_name: &str) -> Option<String> {
        // Find the index of the column name in the columns_has_changed vector
        let column_number = self
            .columns_has_changed
            .iter()
            .position(|col| col.as_ref().map_or(false, |name| name == column_name));

        // If found, take the column name at that position
        if let Some(index) = column_number {
            self.columns_has_changed
                .get_mut(index)
                .and_then(|col| col.take())
        } else {
            None
        }
    }

    pub fn feed_lines(&self) -> Vec<String> {
        self.vector_of_operations
            .iter()
            .map(OperationRecord::feed_line)
            .collect()
    }

    pub fn feed_text(&self) -> String {
        let lines = self.feed_lines();
        if lines.is_empty() {
            "History is empty".to_string()
        } else {
            lines.join("\n")
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ColumnTypes {
    Time,
    Temperature,
    Mass,
    Alpha,
    Dmdt,
    Eta,
    Detadt,
    Dalphadt,
    DTdt,
    Undefined,
    All,
}

#[derive(Clone, Debug)]
pub enum AffectedColumns {
    Specific(Vec<String>),
    All,
    Semantic(Vec<ColumnTypes>),
}

impl std::fmt::Display for AffectedColumns {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AffectedColumns::Specific(cols) => write!(f, "specific({})", cols.join(", ")),
            AffectedColumns::All => write!(f, "all"),
            AffectedColumns::Semantic(types) => write!(f, "semantic({:?})", types),
        }
    }
}

#[derive(Clone, Debug)]
pub struct OperationRecord {
    pub timestamp: usize,
    pub operation_name: String,
    pub affected_columns: AffectedColumns,
    pub expr: Option<Expr>,
    pub description: String,
    pub reversible: bool,
}

impl OperationRecord {
    pub fn feed_line(&self) -> String {
        format!(
            "[{}] {} | {} | affected: {} | reversible: {}",
            self.timestamp,
            self.operation_name,
            self.description,
            self.affected_columns,
            self.reversible
        )
    }
}

impl std::fmt::Display for OperationRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.feed_line())
    }
}

#[derive(Clone, Debug)]
pub struct ColumnHistory {
    pub column_name: String,
    pub operations: Vec<OperationRecord>,
}

impl ColumnHistory {
    pub fn new(column_name: String) -> Self {
        Self {
            column_name,
            operations: Vec::new(),
        }
    }

    pub fn reversible_count(&self) -> usize {
        self.operations.iter().filter(|op| op.reversible).count()
    }

    pub fn irreversible_count(&self) -> usize {
        self.operations.iter().filter(|op| !op.reversible).count()
    }

    pub fn has_irreversible(&self) -> bool {
        self.operations.iter().any(|op| !op.reversible)
    }

    pub fn feed_text(&self) -> String {
        if self.operations.is_empty() {
            format!("Column '{}' has no history", self.column_name)
        } else {
            self.operations
                .iter()
                .map(OperationRecord::feed_line)
                .collect::<Vec<_>>()
                .join("\n")
        }
    }
}

#[derive(Clone)]
pub(crate) struct TGADatasetSnapshot {
    frame: LazyFrame,
    schema: TGASchema,
    oneframeplot: Option<OneFramePlot>,
    columns_has_changed: Vec<Option<String>>,
    operation_index: usize,
}

/// A bounded copy-on-write undo stack. Cloning a dataset shares old snapshots;
/// changing the stack copies only the small deque of `Arc` handles.
#[derive(Clone, Default)]
pub(crate) struct UndoStack {
    entries: Arc<VecDeque<Arc<TGADatasetSnapshot>>>,
}

impl UndoStack {
    const MAX_DEPTH: usize = 20;

    fn push(&mut self, snapshot: TGADatasetSnapshot) {
        let entries = Arc::make_mut(&mut self.entries);
        if entries.len() == Self::MAX_DEPTH {
            entries.pop_front();
        }
        entries.push_back(Arc::new(snapshot));
    }

    fn pop(&mut self) -> Option<TGADatasetSnapshot> {
        let snapshot = Arc::make_mut(&mut self.entries).pop_back()?;
        Some(Arc::try_unwrap(snapshot).unwrap_or_else(|shared| (*shared).clone()))
    }

    fn len(&self) -> usize {
        self.entries.len()
    }

    fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub(crate) fn clear(&mut self) {
        self.entries = Arc::new(VecDeque::new());
    }

    #[cfg(test)]
    pub(crate) fn shares_storage_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.entries, &other.entries)
    }
}

//================================================================================
//              TGADATASET - MAIN DATA STRUCTURE
//==================================================================================
/// Main entity Plolars-layer
#[derive(Clone)]
pub struct TGADataset {
    pub(crate) frame: LazyFrame,
    pub(crate) schema: TGASchema,
    pub oneframeplot: Option<OneFramePlot>,
    pub history_of_operations: History,
    pub(crate) undo_stack: UndoStack,
    pub(crate) undo_snapshot_latch: bool,
}

impl Default for TGADataset {
    fn default() -> Self {
        Self {
            frame: LazyFrame::default(),
            schema: TGASchema::default(),
            oneframeplot: None,
            history_of_operations: History::default(),
            undo_stack: UndoStack::default(),
            undo_snapshot_latch: false,
        }
    }
}

impl std::fmt::Debug for TGADataset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TGADataset")
            .field("schema", &DebugSchema(&self.schema))
            .field("frame", &"<LazyFrame>")
            .finish()
    }
}

pub fn expr_into_str(expr: &Expr) -> Expr {
    expr.clone().reverse()
}
/// Helper struct to nicely format the TGASchema for debugging
struct DebugSchema<'a>(&'a TGASchema);

impl<'a> std::fmt::Debug for DebugSchema<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let schema = self.0;
        let mut builder = f.debug_struct("TGASchema");

        // Show all columns with their metadata
        builder.field("columns_count", &schema.columns.len());

        let mut columns_display: Vec<_> = schema
            .columns
            .values()
            .map(|m| format!("{}: {} ({})", m.name, m.unit, format!("{:?}", m.origin)))
            .collect();
        columns_display.sort();
        builder.field("columns", &columns_display);

        // Show bound semantic columns
        builder.field("time", &schema.time);
        builder.field("temperature", &schema.temperature);
        builder.field("mass", &schema.mass);
        builder.field("alpha", &schema.alpha);
        builder.field("dm_dt", &schema.dm_dt);
        builder.field("eta", &schema.eta);
        builder.field("deta_dt", &schema.deta_dt);
        builder.field("dalpha_dt", &schema.dalpha_dt);
        builder.field("dT_dt", &schema.dT_dt);

        builder.finish()
    }
}
impl TGADataset {
    /// Read-only access to the lazy data frame.
    ///
    /// Mutations must go through dataset operations so history and undo remain aligned.
    pub fn frame(&self) -> &LazyFrame {
        &self.frame
    }

    /// Read-only access to column metadata and active bindings.
    pub fn schema(&self) -> &TGASchema {
        &self.schema
    }

    fn ensure_column_is_bindable(&self, name: &str) -> Result<(), TGADomainError> {
        let frame_schema = self.frame.clone().collect_schema()?;
        let dtype = frame_schema
            .get(name)
            .ok_or_else(|| TGADomainError::ColumnNotFound(name.to_string()))?;
        if !self.schema.columns.contains_key(name) {
            return Err(TGADomainError::ColumnNotFound(name.to_string()));
        }
        if !matches!(
            dtype,
            DataType::Int8
                | DataType::Int16
                | DataType::Int32
                | DataType::Int64
                | DataType::UInt8
                | DataType::UInt16
                | DataType::UInt32
                | DataType::UInt64
                | DataType::Float32
                | DataType::Float64
        ) {
            return Err(TGADomainError::NonNumericColumn {
                column: name.to_string(),
                actual_type: dtype.to_string(),
            });
        }
        Ok(())
    }
    //=========================================================================
    // HISTORY HELPERS

    pub fn log_operation(
        &mut self,
        name: &str,
        affected: AffectedColumns,
        expr: Option<Expr>,
        description: String,
        reversible: bool,
    ) {
        let timestamp = self.history_of_operations.vector_of_operations.len();
        self.populate_columns_has_changed(&affected);
        self.history_of_operations
            .vector_of_operations
            .push(OperationRecord {
                timestamp,
                operation_name: name.to_string(),
                affected_columns: affected,
                expr,
                description,
                reversible,
            });
    }

    /// Store a pre-change snapshot so the last mutation can be reverted later.
    /// We keep the snapshot small enough for day-to-day use by limiting the stack depth.
    pub fn push_undo_snapshot(&mut self) {
        if self.undo_snapshot_latch {
            self.undo_snapshot_latch = false;
            return;
        }
        self.undo_stack.push(TGADatasetSnapshot {
            frame: self.frame.clone(),
            schema: self.schema.clone(),
            oneframeplot: self.oneframeplot.clone(),
            columns_has_changed: self.history_of_operations.columns_has_changed.clone(),
            operation_index: self.history_of_operations.vector_of_operations.len(),
        });
    }

    /// Returns `true` when a column name is already occupied in the schema.
    pub(crate) fn has_column_name(&self, name: &str) -> bool {
        self.schema.columns.contains_key(name)
    }

    /// Reject duplicate names before creating a new column.
    pub(crate) fn ensure_column_name_is_available(&self, name: &str) -> Result<(), TGADomainError> {
        if self.has_column_name(name) || self.frame.clone().collect_schema()?.get(name).is_some() {
            Err(TGADomainError::ColumnAlreadyExists(name.to_string()))
        } else {
            Ok(())
        }
    }

    /// Mark the next mutation as part of a composite operation.
    pub(crate) fn suppress_next_undo_snapshot(&mut self) {
        self.undo_snapshot_latch = true;
    }

    /// Seed provenance for columns that are already present in the dataset.
    ///
    /// This keeps older construction paths working while the new metadata model
    /// is being threaded through the codebase.
    pub(crate) fn initialize_column_provenance(&mut self) {
        for (name, meta) in self.schema.columns.iter_mut() {
            if meta.provenance.steps.is_empty() {
                meta.provenance = match meta.origin {
                    ColumnOrigin::Raw => ColumnProvenance::raw(name.clone()),
                    ColumnOrigin::Imported => {
                        ColumnProvenance::imported(name.clone(), name.clone())
                    }
                    ColumnOrigin::PolarsDerived | ColumnOrigin::NumericDerived => {
                        ColumnProvenance::manual(
                            name.clone(),
                            "initialize_column_provenance",
                            Some("backfilled provenance".to_string()),
                        )
                    }
                };
            }
        }
    }

    /// Return provenance for a single column if we have recorded it.
    pub fn column_provenance(&self, column_name: &str) -> Option<&ColumnProvenance> {
        self.schema
            .columns
            .get(column_name)
            .map(|meta| &meta.provenance)
    }

    /// Human-readable provenance for GUI panels and debugging.
    pub fn column_provenance_text(&self, column_name: &str) -> Option<String> {
        self.column_provenance(column_name)
            .map(ColumnProvenance::feed_text)
    }

    /// Record one transformation step for a column.
    ///
    /// If `source_column` and `target_column` are different, provenance is cloned
    /// from the source and extended with the new step. If the column is modified
    /// in place, the step is appended to the existing chain.
    pub(crate) fn record_column_transform(
        &mut self,
        source_column: &str,
        target_column: &str,
        kind: ColumnTransformKind,
        reversible: bool,
        quality: Option<
            crate::Kinetics::experimental_kinetics::math::ndarray_statistics::FilterQualityReport,
        >,
    ) {
        let operation_id = Some(self.history_of_operations.vector_of_operations.len());
        let input_columns = vec![source_column.to_string()];
        let quality = Some(match quality {
            Some(report) => TransformQuality::computed(report, QualityReference::PreviousStep),
            None => TransformQuality::unavailable(
                QualityReference::PreviousStep,
                "quality metrics could not be computed for this transformation",
            ),
        });

        if source_column == target_column {
            let entry = self
                .schema
                .columns
                .entry(target_column.to_string())
                .or_insert_with(|| {
                    ColumnMeta::raw(
                        target_column.to_string(),
                        Unit::Unknown,
                        ColumnNature::Unknown,
                    )
                });
            entry.provenance.push_step(
                input_columns,
                target_column.to_string(),
                kind,
                quality,
                operation_id,
                reversible,
            );
            return;
        }

        let mut provenance = self
            .schema
            .columns
            .get(source_column)
            .map(|meta| meta.provenance.clone())
            .unwrap_or_else(|| ColumnProvenance::raw(source_column));
        provenance.rename_output(target_column);
        provenance.push_step(
            input_columns,
            target_column.to_string(),
            kind,
            quality,
            operation_id,
            reversible,
        );

        if let Some(target_meta) = self.schema.columns.get_mut(target_column) {
            target_meta.provenance = provenance;
        } else {
            self.schema.columns.insert(
                target_column.to_string(),
                ColumnMeta::new(
                    target_column.to_string(),
                    Unit::Unknown,
                    ColumnOrigin::NumericDerived,
                    ColumnNature::Unknown,
                    provenance,
                ),
            );
        }
    }

    /// Returns `true` if there is at least one saved state to undo to.
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// Number of saved snapshots available for undo.
    pub fn undo_depth(&self) -> usize {
        self.undo_stack.len()
    }

    /// Restore the last saved state and add an `undo` record to the audit trail.
    pub fn undo_last(&mut self) -> Result<(), TGADomainError> {
        let Some(snapshot) = self.undo_stack.pop() else {
            return Err(TGADomainError::InvalidOperation(
                "No undo snapshot available".to_string(),
            ));
        };

        let undone = self
            .history_of_operations
            .vector_of_operations
            .get(snapshot.operation_index)
            .cloned();

        self.frame = snapshot.frame;
        self.schema = snapshot.schema;
        self.oneframeplot = snapshot.oneframeplot;
        self.history_of_operations.columns_has_changed = snapshot.columns_has_changed;

        let (affected, description) = if let Some(op) = undone {
            (
                op.affected_columns,
                format!("Reverted operation '{}'", op.operation_name),
            )
        } else {
            (AffectedColumns::All, "Reverted last operation".to_string())
        };

        self.log_operation("undo_last", affected, None, description, true);
        Ok(())
    }

    fn populate_columns_has_changed(&mut self, affected: &AffectedColumns) {
        let cols = match affected {
            AffectedColumns::Specific(cols) => cols.clone(),
            AffectedColumns::All => self.list_of_columns(),
            AffectedColumns::Semantic(types) => {
                let mut result = Vec::new();
                for col in self.list_of_columns() {
                    if types.iter().any(|ct| self.column_matches_type(&col, ct)) {
                        result.push(col);
                    }
                }
                result
            }
        };

        for col in cols {
            if !self
                .history_of_operations
                .columns_has_changed
                .contains(&Some(col.clone()))
            {
                self.history_of_operations
                    .columns_has_changed
                    .push(Some(col));
            }
        }
    }
    ///
    pub fn list_of_columns_to_recalc(&mut self) -> Vec<String> {
        let list_of_columns = self.list_of_columns();
        let mut list_of_columns_to_recalc = Vec::new();
        for col in list_of_columns {
            if let Some(col_to_recalc) = self.take_column(&col) {
                list_of_columns_to_recalc.push(col_to_recalc);
            }
        }
        list_of_columns_to_recalc
    }

    fn column_matches_type(&self, col: &str, ct: &ColumnTypes) -> bool {
        match ct {
            ColumnTypes::Time => self.schema.time.as_deref() == Some(col),
            ColumnTypes::Temperature => self.schema.temperature.as_deref() == Some(col),
            ColumnTypes::Mass => self.schema.mass.as_deref() == Some(col),
            ColumnTypes::Alpha => self.schema.alpha.as_deref() == Some(col),
            ColumnTypes::Dmdt => self.schema.dm_dt.as_deref() == Some(col),
            ColumnTypes::Eta => self.schema.eta.as_deref() == Some(col),
            ColumnTypes::Detadt => self.schema.deta_dt.as_deref() == Some(col),
            ColumnTypes::Dalphadt => self.schema.dalpha_dt.as_deref() == Some(col),
            ColumnTypes::DTdt => self.schema.dT_dt.as_deref() == Some(col),
            ColumnTypes::All => true,
            ColumnTypes::Undefined => false,
        }
    }

    pub fn operations_on_column(&self, col: &str) -> Vec<OperationRecord> {
        let h: Vec<&OperationRecord> = self
            .history_of_operations
            .vector_of_operations
            .iter()
            .filter(|op| match &op.affected_columns {
                AffectedColumns::Specific(cols) => cols.contains(&col.to_string()),
                AffectedColumns::All => true,
                AffectedColumns::Semantic(types) => {
                    types.iter().any(|ct| self.column_matches_type(col, ct))
                }
            })
            .collect();
        h.iter().map(|&op| op.clone()).collect()
    }

    pub fn get_column_history(&self, col: &str) -> ColumnHistory {
        let operations = self.operations_on_column(col).into_iter().collect();
        ColumnHistory {
            column_name: col.to_string(),
            operations,
        }
    }

    pub fn list_of_columns(&self) -> Vec<String> {
        self.schema.columns.keys().cloned().collect()
    }

    pub fn print_operation_history(&self) {
        println!("\n=== Operation History ===");
        for op in &self.history_of_operations.vector_of_operations {
            println!(
                "[{}] {} - {} (reversible: {})",
                op.timestamp, op.operation_name, op.description, op.reversible
            );
            match &op.affected_columns {
                AffectedColumns::Specific(cols) => println!("  Columns: {:?}", cols),
                AffectedColumns::All => println!("  Columns: ALL"),
                AffectedColumns::Semantic(types) => println!("  Semantic: {:?}", types),
            }
        }
    }

    pub fn print_column_history(&self, col: &str) {
        let history = self.get_column_history(col);
        println!("\n=== History for column '{}' ===", col);
        println!("Total operations: {}", history.operations.len());
        println!("Reversible: {}", history.reversible_count());
        println!("Irreversible: {}", history.irreversible_count());
        for op in &history.operations {
            println!(
                "  [{}] {} - {} (reversible: {})",
                op.timestamp, op.operation_name, op.description, op.reversible
            );
        }
    }

    pub fn columns_has_changed(&self) -> Vec<Option<String>> {
        self.history_of_operations.columns_has_changed()
    }
    pub fn take_column(&mut self, column_name: &str) -> Option<String> {
        self.history_of_operations.take_column(column_name)
    }

    //=========================================================================
    // INPUT/OUTPUT
    /// Normalize txt file to csv
    /// ?????????????? txt ????? ? csv ??????
    pub fn normalize_txt_to_csv(input: &Path, output: &Path) -> std::io::Result<()> {
        let input = std::fs::read_to_string(input)?;
        let mut out = String::new();

        for line in input.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let cols: Vec<&str> = line.split_whitespace().collect();
            out.push_str(&cols.join(","));
            out.push('\n');
        }

        std::fs::write(output, out)?;
        Ok(())
    }
    /// creation of TGA entity
    /// ???????? TGA ???????? ?? CSV ?????
    pub fn from_csv(
        path: &str,
        time_col: &str,
        temp_col: &str,
        mass_col: &str,
    ) -> PolarsResult<Self> {
        let requested_roles = [time_col, temp_col, mass_col]
            .into_iter()
            .collect::<HashSet<_>>();
        if requested_roles.len() != 3 {
            return Err(PolarsError::ComputeError(
                "Time, temperature, and mass must use distinct columns".into(),
            ));
        }
        println!("\n column: {} ", time_col);
        println!("\n column: {} ", temp_col);
        println!("\n column: {} ", mass_col);

        let plpath = PlRefPath::new(path);
        let frame = LazyCsvReader::new(plpath)
            .with_has_header(true)
            .finish()?
            .with_columns(vec![
                col(time_col).cast(DataType::Float64).alias(time_col),
                col(temp_col).cast(DataType::Float64).alias(temp_col),
                col(mass_col).cast(DataType::Float64).alias(mass_col),
            ]);

        let mut columns = HashMap::new();

        columns.insert(
            time_col.to_string(),
            ColumnMeta::raw(time_col, Unit::Second, ColumnNature::Time),
        );

        columns.insert(
            temp_col.into(),
            ColumnMeta::raw(temp_col, Unit::Celsius, ColumnNature::Temperature),
        );

        columns.insert(
            mass_col.to_string(),
            ColumnMeta::raw(mass_col, Unit::Milligram, ColumnNature::Mass),
        );

        let mut dataset = Self {
            frame,
            schema: TGASchema {
                columns,
                time: Some(time_col.into()),
                temperature: Some(temp_col.into()),
                mass: Some(mass_col.into()),
                ..TGASchema::default()
            },
            oneframeplot: None,
            history_of_operations: History::new(),
            undo_stack: UndoStack::default(),
            undo_snapshot_latch: false,
        };
        dataset.initialize_column_provenance();

        dataset.log_operation(
            "from_csv",
            AffectedColumns::Specific(vec![
                time_col.to_string(),
                temp_col.to_string(),
                mass_col.to_string(),
            ]),
            None,
            format!(
                "Imported CSV '{}' with bound columns time='{}', temperature='{}', mass='{}'",
                path, time_col, temp_col, mass_col
            ),
            false,
        );

        Ok(dataset)
    }

    /// Add a numeric column from a user-provided Vec<f64>.
    ///
    /// This is useful for fitted curves or any external numeric results.
    /// The function checks that the vector length matches the dataset height,
    /// then appends the column as `NumericDerived` and logs the operation.
    pub fn add_column_from_vec(
        mut self,
        name: &str,
        unit: Unit,
        nature: ColumnNature,
        data: Vec<f64>,
    ) -> Result<Self, TGADomainError> {
        self.ensure_column_name_is_available(name)?;
        let height = self.frame.clone().collect()?.height();
        if data.len() != height {
            return Err(TGADomainError::InvalidOperation(format!(
                "Column '{}' length {} does not match dataset height {}",
                name,
                data.len(),
                height
            )));
        }

        let series = Series::new(name.into(), data);
        self.push_undo_snapshot();
        self.frame = self.frame.with_column(lit(series));

        self.schema.columns.insert(
            name.to_string(),
            ColumnMeta::new(
                name.to_string(),
                unit,
                ColumnOrigin::NumericDerived,
                nature,
                ColumnProvenance::manual(
                    name,
                    "add_column_from_vec",
                    Some("external numeric data".to_string()),
                ),
            ),
        );

        self.log_operation(
            "add_column_from_vec",
            AffectedColumns::Specific(vec![name.to_string()]),
            None,
            format!("Added numeric column '{}' from Vec<f64>", name),
            true,
        );

        Ok(self)
    }

    /// Создание TGA сущности из необработанного CSV файла
    pub fn from_csv_raw(path: &Path) -> PolarsResult<Self> {
        let plpath = PlRefPath::try_from_path(path)?;

        let mut frame = LazyCsvReader::new(plpath).with_has_header(true).finish()?;

        let schema = frame.collect_schema()?;

        let mut columns = HashMap::new();
        for name in schema.iter_names() {
            columns.insert(
                name.to_string(),
                ColumnMeta::raw(name.to_string(), Unit::Unknown, ColumnNature::Unknown),
            );
        }

        let mut dataset = Self {
            frame,
            schema: TGASchema {
                columns,
                ..TGASchema::default()
            },
            oneframeplot: None,
            history_of_operations: History::new(),
            undo_stack: UndoStack::default(),
            undo_snapshot_latch: false,
        };
        dataset.initialize_column_provenance();

        dataset.log_operation(
            "from_csv_raw",
            AffectedColumns::All,
            None,
            format!(
                "Imported raw CSV '{}' with {} columns",
                path.display(),
                dataset.schema.columns.len()
            ),
            false,
        );

        Ok(dataset)
    }

    /// Привязка временной колонки с указанием единиц измерения
    /// Build an appended dataset while leaving both source datasets unchanged.
    /// The input schemas and time bindings must match. A time axis that already
    /// continues is preserved; a reset/overlapping axis is shifted by the median
    /// of the last ten target sample intervals.
    pub fn append_continuation(&self, next: &Self) -> Result<Self, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_deref()
            .ok_or(TGADomainError::TimeNotBound)?;
        let next_time_col = next
            .schema
            .time
            .as_deref()
            .ok_or(TGADomainError::TimeNotBound)?;
        if time_col != next_time_col {
            return Err(TGADomainError::InvalidOperation(format!(
                "Continuation time column '{}' does not match target time column '{}'",
                next_time_col, time_col
            )));
        }
        if self.schema.columns.get(time_col).map(|meta| meta.unit)
            != next.schema.columns.get(next_time_col).map(|meta| meta.unit)
        {
            return Err(TGADomainError::InvalidOperation(
                "Continuation time units must match the target experiment".into(),
            ));
        }

        let target_schema = self.frame.clone().collect_schema()?;
        let next_schema = next.frame.clone().collect_schema()?;
        let target_names = target_schema
            .iter_names()
            .map(|name| name.to_string())
            .collect::<Vec<_>>();
        let next_names = next_schema
            .iter_names()
            .map(|name| name.to_string())
            .collect::<Vec<_>>();
        if target_names != next_names {
            return Err(TGADomainError::InvalidOperation(
                "Continuation columns and their order must match the target experiment".into(),
            ));
        }
        for name in &target_names {
            if target_schema.get(name) != next_schema.get(name) {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Continuation schema differs for column '{name}'"
                )));
            }
        }

        let read_time = |dataset: &Self| -> Result<Vec<f64>, TGADomainError> {
            let frame = dataset
                .frame
                .clone()
                .select([col(time_col).cast(DataType::Float64)])
                .collect()?;
            let time = frame.column(time_col)?.f64()?;
            let values = (0..time.len())
                .map(|index| {
                    time.get(index).filter(|value| value.is_finite()).ok_or_else(|| {
                    TGADomainError::InvalidOperation(format!(
                        "Time column '{time_col}' has null or non-finite value at row {index}"
                    ))
                })
                })
                .collect::<Result<Vec<_>, _>>()?;
            if values.len() < 2 || values.windows(2).any(|pair| pair[1] <= pair[0]) {
                return Err(TGADomainError::TimeNonMonotonic(time_col.to_string()));
            }
            Ok(values)
        };
        let target_time = read_time(self)?;
        let next_time = read_time(next)?;
        let mut steps = target_time
            .windows(2)
            .rev()
            .take(10)
            .map(|pair| pair[1] - pair[0])
            .collect::<Vec<_>>();
        steps.sort_by(f64::total_cmp);
        let step = steps[steps.len() / 2];
        let time_shift = if next_time[0] > target_time[target_time.len() - 1] {
            0.0
        } else {
            target_time[target_time.len() - 1] + step - next_time[0]
        };
        let next_frame = if time_shift == 0.0 {
            next.frame.clone()
        } else {
            next.frame.clone().with_column(
                (col(time_col).cast(DataType::Float64) + lit(time_shift)).alias(time_col),
            )
        };
        let mut result = self.clone();
        result.frame = concat([self.frame.clone(), next_frame], UnionArgs::default())?;
        result.log_operation(
            "append_continuation",
            AffectedColumns::All,
            None,
            format!("Appended continuation; time shift for '{time_col}' is {time_shift}"),
            false,
        );
        Ok(result)
    }

    pub fn bind_time_inplace(&mut self, col: &str, unit: Unit) -> Result<(), TGADomainError> {
        if !unit.is_valid_time_binding_unit() {
            return Err(TGADomainError::InvalidUnitForBinding {
                field: "time",
                unit,
            });
        }
        self.ensure_column_is_bindable(col)?;
        self.push_undo_snapshot();
        let meta = self
            .schema
            .columns
            .get_mut(col)
            .ok_or(TGADomainError::ColumnNotFound(col.into()))?;

        meta.unit = unit;
        meta.nature = ColumnNature::Time;
        self.schema.bind_primary_role(&ColumnRole::Time, col);
        self.record_column_transform(
            col,
            col,
            ColumnTransformKind::Binding {
                role: "time".to_string(),
                unit: format!("{:?}", unit),
            },
            true,
            None,
        );

        self.log_operation(
            "bind_time",
            AffectedColumns::Semantic(vec![ColumnTypes::Time]),
            None,
            format!("Bound time column '{}' with unit {:?}", col, unit),
            true,
        );

        Ok(())
    }

    pub fn bind_time(mut self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        self.bind_time_inplace(col, unit)?;
        Ok(self)
    }

    /// Привязка температурной колонки с указанием единиц измерения
    pub fn bind_temperature_inplace(
        &mut self,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGADomainError> {
        if !unit.is_valid_temperature_binding_unit() {
            return Err(TGADomainError::InvalidUnitForBinding {
                field: "temperature",
                unit,
            });
        }
        self.ensure_column_is_bindable(col)?;
        self.push_undo_snapshot();
        let meta = self
            .schema
            .columns
            .get_mut(col)
            .ok_or(TGADomainError::ColumnNotFound(col.into()))?;

        meta.unit = unit;
        meta.nature = ColumnNature::Temperature;
        self.schema.bind_primary_role(&ColumnRole::Temperature, col);
        self.record_column_transform(
            col,
            col,
            ColumnTransformKind::Binding {
                role: "temperature".to_string(),
                unit: format!("{:?}", unit),
            },
            true,
            None,
        );

        self.log_operation(
            "bind_temperature",
            AffectedColumns::Semantic(vec![ColumnTypes::Temperature]),
            None,
            format!("Bound temperature column '{}' with unit {:?}", col, unit),
            true,
        );

        Ok(())
    }

    pub fn bind_temperature(mut self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        self.bind_temperature_inplace(col, unit)?;
        Ok(self)
    }

    /// Привязка массовой колонки с указанием единиц измерения
    pub fn bind_mass_inplace(&mut self, col: &str, unit: Unit) -> Result<(), TGADomainError> {
        if !unit.is_valid_mass_binding_unit() {
            return Err(TGADomainError::InvalidUnitForBinding {
                field: "mass",
                unit,
            });
        }
        self.ensure_column_is_bindable(col)?;
        self.push_undo_snapshot();
        let meta = self
            .schema
            .columns
            .get_mut(col)
            .ok_or(TGADomainError::ColumnNotFound(col.into()))?;

        meta.unit = unit;
        meta.nature = ColumnNature::Mass;
        self.schema.bind_primary_role(&ColumnRole::Mass, col);
        self.record_column_transform(
            col,
            col,
            ColumnTransformKind::Binding {
                role: "mass".to_string(),
                unit: format!("{:?}", unit),
            },
            true,
            None,
        );

        self.log_operation(
            "bind_mass",
            AffectedColumns::Semantic(vec![ColumnTypes::Mass]),
            None,
            format!("Bound mass column '{}' with unit {:?}", col, unit),
            true,
        );

        Ok(())
    }

    pub fn bind_mass(mut self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        self.bind_mass_inplace(col, unit)?;
        Ok(self)
    }
    /// Привязка колонки с указанием роли, имени и единиц измерения
    pub fn bind_column_inplace(
        &mut self,
        role: ColumnRole,
        name: &str,
        unit: Unit,
    ) -> Result<(), TGADomainError> {
        let role_label = format!("{:?}", &role);
        let valid = match role {
            ColumnRole::Time => unit.is_valid_time_binding_unit(),
            ColumnRole::Temperature => unit.is_valid_temperature_binding_unit(),
            ColumnRole::Mass => unit.is_valid_mass_binding_unit(),
        };
        if !valid {
            let field = match role {
                ColumnRole::Time => "time",
                ColumnRole::Temperature => "temperature",
                ColumnRole::Mass => "mass",
            };
            return Err(TGADomainError::InvalidUnitForBinding { field, unit });
        }
        let nature = match role {
            ColumnRole::Time => ColumnNature::Time,
            ColumnRole::Temperature => ColumnNature::Temperature,
            ColumnRole::Mass => ColumnNature::Mass,
        };
        self.ensure_column_is_bindable(name)?;
        self.push_undo_snapshot();
        let meta = self
            .schema
            .columns
            .get_mut(name)
            .ok_or_else(|| TGADomainError::ColumnNotFound(name.to_string()))?;
        meta.unit = unit;
        meta.nature = nature;

        self.schema.bind_primary_role(&role, name);
        self.record_column_transform(
            name,
            name,
            ColumnTransformKind::Binding {
                role: role_label,
                unit: format!("{:?}", unit),
            },
            true,
            None,
        );

        let affected = match role {
            ColumnRole::Time => AffectedColumns::Semantic(vec![ColumnTypes::Time]),
            ColumnRole::Temperature => AffectedColumns::Semantic(vec![ColumnTypes::Temperature]),
            ColumnRole::Mass => AffectedColumns::Semantic(vec![ColumnTypes::Mass]),
        };

        self.log_operation(
            "bind_column",
            affected,
            None,
            format!("Bound {:?} column '{}' with unit {:?}", role, name, unit),
            true,
        );

        Ok(())
    }

    pub fn bind_column(
        mut self,
        role: ColumnRole,
        name: &str,
        unit: Unit,
    ) -> Result<Self, TGADomainError> {
        self.bind_column_inplace(role, name, unit)?;
        Ok(self)
    }

    /// Экспорт данных в CSV файл с указанием единиц измерения в заголовке
    pub fn to_csv_with_units(&self, path: &Path) -> PolarsResult<()> {
        let df = &mut self.frame.clone().collect()?;
        let mut file = std::fs::File::create(path)?;

        // Meta header
        let meta: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|name| {
                let unit = self
                    .schema
                    .get(name)
                    .map(|m| m.unit.convert_to_string())
                    .unwrap_or("?".into());
                format!("{} [{}]", name, unit)
            })
            .collect();

        writeln!(file, "# KiThe TGA Dataset")?;
        writeln!(file, "# {}", meta.join(", "))?;
        let bindings = CsvSchemaBindings::from_schema(&self.schema, None);
        writeln!(
            file,
            "# Bindings: {}",
            serde_json::to_string(&bindings).map_err(|err| {
                PolarsError::ComputeError(format!("Failed to serialize CSV bindings: {err}").into())
            })?
        )?;
        let names = df
            .get_column_names()
            .iter()
            .map(|name| name.to_string())
            .collect::<Vec<_>>();
        let column_metadata = CsvColumnMetadata::from_schema(&self.schema, &names);
        writeln!(
            file,
            "# Column metadata: {}",
            serde_json::to_string(&column_metadata).map_err(|err| {
                PolarsError::ComputeError(
                    format!("Failed to serialize CSV column metadata: {err}").into(),
                )
            })?
        )?;

        CsvWriter::new(&mut file).include_header(true).finish(df)?;

        Ok(())
    }

    /// Export selected columns to CSV with units encoded in a metadata header, similar to `to_csv_with_units`.
    /// Only the specified columns are exported, and units are taken from the schema.
    /// ??????? ????????? ??????? ? CSV ???? ? ????????? ?????? ?????????.
    /// ??????? ????????? ??????? ? CSV ???? ? ????????? ?????? ?????????.
    /// Экспорт выбранных колонок в CSV файл с указанием единиц измерения
    pub fn to_csv_with_units_selected(&self, path: &Path, columns: &[&str]) -> PolarsResult<()> {
        if columns.is_empty() {
            // Nothing to export
            return Ok(());
        }

        // Validate that requested columns exist in the frame schema
        let schema = self.frame.clone().collect_schema()?;
        for &name in columns {
            if !schema.iter_names().any(|n| n == name) {
                return Err(PolarsError::ColumnNotFound(
                    format!("Column '{}' not found", name).into(),
                ));
            }
        }

        // Build selection in requested order
        let mut select_exprs: Vec<Expr> = Vec::with_capacity(columns.len());
        for &name in columns {
            select_exprs.push(col(name));
            println!("\n column: {} ", name);
        }

        let mut df = self.frame.clone().select(select_exprs).collect()?;

        // Construct metadata header line
        let header_meta: Vec<String> = columns
            .iter()
            .map(|&name| {
                let unit = self
                    .schema
                    .get(name)
                    .map(|m| m.unit.convert_to_string())
                    .unwrap_or_else(|| "?".to_string());
                format!("{} [{}]", name, unit)
            })
            .collect();

        let mut file = std::fs::File::create(path)?;
        writeln!(file, "# KiThe TGA Dataset (selected columns)")?;
        writeln!(file, "# {}", header_meta.join(", "))?;
        let bindings = CsvSchemaBindings::from_schema(&self.schema, Some(columns));
        writeln!(
            file,
            "# Bindings: {}",
            serde_json::to_string(&bindings).map_err(|err| {
                PolarsError::ComputeError(format!("Failed to serialize CSV bindings: {err}").into())
            })?
        )?;
        let column_names = columns
            .iter()
            .map(|name| (*name).to_string())
            .collect::<Vec<_>>();
        let column_metadata = CsvColumnMetadata::from_schema(&self.schema, &column_names);
        writeln!(
            file,
            "# Column metadata: {}",
            serde_json::to_string(&column_metadata).map_err(|err| {
                PolarsError::ComputeError(
                    format!("Failed to serialize CSV column metadata: {err}").into(),
                )
            })?
        )?;

        CsvWriter::new(&mut file)
            .include_header(true)
            .finish(&mut df)?;
        Ok(())
    }

    /// Создание TGA сущности из CSV файла с единицами измерения
    pub fn from_csv_with_units(path: &Path) -> PolarsResult<Self> {
        let text = std::fs::read_to_string(path)?;
        let mut lines = text.lines().peekable();

        lines.next(); // "# KiThe TGA Dataset"
        let meta_line = lines.next().ok_or_else(|| {
            PolarsError::ComputeError("Missing column metadata line in TGA CSV".into())
        })?;
        let meta_line = meta_line.strip_prefix("# ").ok_or_else(|| {
            PolarsError::ComputeError("Malformed column metadata line in TGA CSV".into())
        })?;

        let mut names = Vec::new();
        let mut units = Vec::new();

        for entry in meta_line.split(',') {
            let (name, unit_with_bracket) = entry.trim().split_once('[').ok_or_else(|| {
                PolarsError::ComputeError(
                    format!("Malformed column metadata entry '{entry}'").into(),
                )
            })?;
            let unit = unit_with_bracket.strip_suffix(']').ok_or_else(|| {
                PolarsError::ComputeError(
                    format!("Malformed unit in metadata entry '{entry}'").into(),
                )
            })?;
            names.push(name.trim().to_string());
            units.push(Unit::parse(unit)?);
        }

        let mut bindings = None;
        let mut column_metadata = None;
        while let Some(line) = lines.peek().copied() {
            if let Some(json) = line.strip_prefix("# Bindings: ") {
                lines.next();
                bindings = Some(serde_json::from_str::<CsvSchemaBindings>(json).map_err(
                    |err| {
                        PolarsError::ComputeError(format!("Malformed CSV bindings: {err}").into())
                    },
                )?);
            } else if let Some(json) = line.strip_prefix("# Column metadata: ") {
                lines.next();
                column_metadata = Some(
                    serde_json::from_str::<Vec<CsvColumnMetadata>>(json).map_err(|err| {
                        PolarsError::ComputeError(
                            format!("Malformed CSV column metadata: {err}").into(),
                        )
                    })?,
                );
            } else {
                break;
            }
        }

        let data = lines.collect::<Vec<_>>().join("\n");
        let tmp = tempfile::NamedTempFile::new()?;
        std::fs::write(tmp.path(), data)?;

        let plpath = PlRefPath::try_from_path(tmp.path())?;
        // Materialize while the temporary CSV still exists, then release the
        // file normally; the returned lazy frame owns the in-memory data.
        let frame = LazyCsvReader::new(plpath)
            .with_has_header(true)
            .finish()?
            .collect()?
            .lazy();

        let mut schema = TGASchema::from_columns(names, units)?;
        if let Some(column_metadata) = column_metadata {
            CsvColumnMetadata::apply_all(column_metadata, &mut schema)?;
        }
        if let Some(bindings) = bindings {
            bindings.apply(&mut schema)?;
        }

        let mut dataset = Self {
            frame,
            schema,
            oneframeplot: None,
            history_of_operations: History::new(),
            undo_stack: UndoStack::default(),
            undo_snapshot_latch: false,
        };
        dataset.initialize_column_provenance();

        dataset.log_operation(
            "from_csv_with_units",
            AffectedColumns::All,
            None,
            format!(
                "Imported CSV '{}' with explicit units in header",
                path.display()
            ),
            false,
        );

        Ok(dataset)
    }
    /// There are 2 functions to read and parse from_csv_raw a function to read raw files with no metadata and headers and
    /// from_csv_with_units. This is a wrapper function from_csv_universal which recognizes if there in the file headers and call
    /// from_csv_with_units otherwise from_csv_raw
    /// Универсальная функция чтения CSV файлов (с или без единиц измерения)
    pub fn from_csv_universal(path: &Path) -> PolarsResult<Self> {
        let extension = path.extension().and_then(|ext| ext.to_str()).unwrap_or("");

        // If it's a txt file, normalize it to CSV first
        if extension.eq_ignore_ascii_case("txt") {
            let temp_file = tempfile::NamedTempFile::new()?;
            let temp_path = temp_file.path().to_path_buf();

            // Normalize TXT to CSV
            Self::normalize_txt_to_csv(path, &temp_path)?;

            // Process the normalized CSV file
            let result = Self::from_csv_universal_inner(&temp_path)?;

            // Keep temp file alive until processing is done
            std::mem::forget(temp_file);
            return Ok(result);
        }

        // For non-TXT files, proceed with original logic
        Self::from_csv_universal_inner(path)
    }

    /// Inner function that contains the original logic from from_csv_universal
    fn from_csv_universal_inner(path: &Path) -> PolarsResult<Self> {
        let text = std::fs::read_to_string(path)?;
        let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());

        let first = lines.next().unwrap_or_default();
        let second = lines.next().unwrap_or_default();

        let has_units_header = first.starts_with("# KiThe TGA Dataset")
            && second.starts_with("# ")
            && second.contains('[')
            && second.contains(']');

        if has_units_header {
            Self::from_csv_with_units(path)
        } else {
            Self::from_csv_raw(path)
        }
    }

    //==============================================================
    // RATES
    /// Вычисление производной с заданными параметрами
    pub fn derive_rate0(
        self,
        source_col: &str,
        new_col: &str,
        out_meta: ColumnMeta,
    ) -> Result<Self, TGADomainError> {
        self.derive_rate0_with_gap_policy(source_col, new_col, out_meta, TimeGapPolicy::default())
    }

    /// Derive a column with an explicit policy for large gaps in the bound time axis.
    pub fn derive_rate0_with_gap_policy(
        mut self,
        source_col: &str,
        new_col: &str,
        out_meta: ColumnMeta,
        gap_policy: TimeGapPolicy,
    ) -> Result<Self, TGADomainError> {
        self.ensure_column_name_is_available(new_col)?;
        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();
        let time_unit = self
            .schema
            .columns
            .get(&time)
            .map(|meta| meta.unit)
            .ok_or(TGADomainError::ColumnNotFound(time.clone()))?;

        self.schema
            .columns
            .get(source_col)
            .ok_or(TGADomainError::ColumnNotFound(source_col.into()))?;
        let rate_expr = self.derivative_expr(source_col, &time, new_col, gap_policy)?;
        self.push_undo_snapshot();
        self.frame = self.frame.with_column(rate_expr.clone());

        let out_unit = out_meta.unit.for_time_base(time_unit)?;
        let mut out_meta = out_meta;
        out_meta.unit = out_unit;
        self.schema.columns.insert(new_col.into(), out_meta);

        self.log_operation(
            "derive_rate0",
            AffectedColumns::Specific(vec![new_col.to_string()]),
            Some(rate_expr),
            format!(
                "Computed rate {} from {} with unit {:?}; time-gap policy {:?}",
                new_col, source_col, out_unit, gap_policy
            ),
            true,
        );

        Ok(self)
    }
    /// Вычисление производной с автоматическим определением единиц измерения
    pub fn derive_rate(
        self,
        source_col: &str,
        new_col: &str,
        new_col_nature: ColumnNature,
    ) -> Result<Self, TGADomainError> {
        self.derive_rate_with_gap_policy(
            source_col,
            new_col,
            new_col_nature,
            TimeGapPolicy::default(),
        )
    }

    /// Derive a rate with an explicit policy for large gaps in the bound time axis.
    pub fn derive_rate_with_gap_policy(
        mut self,
        source_col: &str,
        new_col: &str,
        new_col_nature: ColumnNature,
        gap_policy: TimeGapPolicy,
    ) -> Result<Self, TGADomainError> {
        self.ensure_column_name_is_available(new_col)?;
        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();
        let time_unit = self
            .schema
            .columns
            .get(&time)
            .map(|meta| meta.unit)
            .ok_or(TGADomainError::ColumnNotFound(time.clone()))?;

        let src_meta = self
            .schema
            .columns
            .get(source_col)
            .ok_or(TGADomainError::ColumnNotFound(source_col.into()))?;

        let out_unit = src_meta.unit.for_time_base(time_unit)?;
        let rate_expr = self.derivative_expr(source_col, &time, new_col, gap_policy)?;
        self.push_undo_snapshot();
        self.frame = self.frame.with_column(rate_expr.clone());

        self.schema.columns.insert(
            new_col.into(),
            ColumnMeta::new(
                new_col.to_string(),
                out_unit,
                ColumnOrigin::PolarsDerived,
                new_col_nature,
                ColumnProvenance::manual(
                    new_col,
                    "derive_rate",
                    Some(format!("derived from {}", source_col)),
                ),
            ),
        );

        match new_col_nature {
            ColumnNature::MassRate => self.schema.dm_dt = Some(new_col.to_string()),
            ColumnNature::TemperatureRate => self.schema.dT_dt = Some(new_col.to_string()),
            ColumnNature::ConversionRate => self.schema.deta_dt = Some(new_col.to_string()),
            ColumnNature::DimensionlessMassRate => {
                self.schema.dalpha_dt = Some(new_col.to_string())
            }
            _ => {}
        }

        self.log_operation(
            "derive_rate",
            AffectedColumns::Specific(vec![new_col.to_string()]),
            Some(rate_expr),
            format!(
                "Computed rate {} from {} with unit {:?}; time-gap policy {:?}",
                new_col, source_col, out_unit, gap_policy
            ),
            true,
        );

        Ok(self)
    }

    /// Build a nonuniform three-point derivative expression without dropping rows.
    fn derivative_expr(
        &self,
        source_col: &str,
        time_col: &str,
        new_col: &str,
        gap_policy: TimeGapPolicy,
    ) -> Result<Expr, TGADomainError> {
        let max_interval = match gap_policy {
            TimeGapPolicy::MedianStepMultiple(factor) => {
                if !factor.is_finite() || factor <= 0.0 {
                    return Err(TGADomainError::InvalidFilterParameters {
                        filter: "time derivative".into(),
                        details: "median-step gap factor must be finite and positive".into(),
                    });
                }
                let limit = self.median_time_step(time_col)? * factor;
                if !limit.is_finite() || limit <= 0.0 {
                    return Err(TGADomainError::InvalidFilterParameters {
                        filter: "time derivative".into(),
                        details: "computed maximum time interval is not finite and positive".into(),
                    });
                }
                Some(limit)
            }
            TimeGapPolicy::MaxInterval(interval) => {
                if !interval.is_finite() || interval <= 0.0 {
                    return Err(TGADomainError::InvalidFilterParameters {
                        filter: "time derivative".into(),
                        details: "maximum time interval must be finite and positive".into(),
                    });
                }
                Some(interval)
            }
            TimeGapPolicy::Allow => None,
        };

        let previous_step = col(time_col) - col(time_col).shift(lit(1));
        let next_step = col(time_col).shift(lit(-1)) - col(time_col);
        let previous_value = col(source_col).shift(lit(1));
        let current_value = col(source_col);
        let next_value = col(source_col).shift(lit(-1));
        let sum_steps = previous_step.clone() + next_step.clone();

        let mut valid = previous_step
            .clone()
            .gt(lit(0.0))
            .and(previous_step.clone().lt(lit(f64::INFINITY)))
            .and(next_step.clone().gt(lit(0.0)))
            .and(next_step.clone().lt(lit(f64::INFINITY)));
        if let Some(max_interval) = max_interval {
            valid = valid
                .and(previous_step.clone().lt_eq(lit(max_interval)))
                .and(next_step.clone().lt_eq(lit(max_interval)));
        }
        valid = valid
            .and(sum_steps.clone().gt(lit(0.0)))
            .and(sum_steps.clone().lt(lit(f64::INFINITY)));

        let first_weight = -next_step.clone() / (previous_step.clone() * sum_steps.clone());
        let center_weight = (next_step.clone() - previous_step.clone())
            / (previous_step.clone() * next_step.clone());
        let last_weight = previous_step.clone() / (next_step.clone() * sum_steps);
        let derivative = first_weight * previous_value
            + center_weight * current_value
            + last_weight * next_value;
        valid = valid
            .and(derivative.clone().gt(lit(-f64::MAX)))
            .and(derivative.clone().lt(lit(f64::MAX)));

        Ok(when(valid)
            .then(derivative)
            .otherwise(lit(NULL))
            .alias(new_col))
    }

    /// Median positive interval between adjacent finite time values.
    fn median_time_step(&self, time_col: &str) -> Result<f64, TGADomainError> {
        let frame = self
            .frame
            .clone()
            .select([col(time_col).cast(DataType::Float64)])
            .collect()?;
        let values = frame.column(time_col)?.f64()?;
        let mut intervals = Vec::with_capacity(values.len().saturating_sub(1));
        let mut previous: Option<f64> = None;
        for index in 0..values.len() {
            match values.get(index).filter(|value| value.is_finite()) {
                Some(current) => {
                    if let Some(previous) = previous {
                        let interval = current - previous;
                        if interval.is_finite() && interval > 0.0 {
                            intervals.push(interval);
                        }
                    }
                    previous = Some(current);
                }
                None => previous = None,
            }
        }

        if intervals.is_empty() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Cannot infer typical sampling interval from time column '{time_col}'"
            )));
        }
        let middle = intervals.len() / 2;
        intervals.select_nth_unstable_by(middle, f64::total_cmp);
        Ok(intervals[middle])
    }

    /// Вычисление скорости изменения массы
    pub fn derive_mass_rate_from_column(
        self,
        source_col: &str,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let mut updated = self.derive_rate(source_col, new_col, ColumnNature::MassRate)?;
        updated.schema.dm_dt = Some(new_col.to_string());
        Ok(updated)
    }

    /// Вычисление скорости изменения массы
    pub fn derive_mass_rate(self, new_col: &str) -> Result<Self, TGADomainError> {
        let mass = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();
        self.derive_mass_rate_from_column(&mass, new_col)
    }

    /// Вычисление скорости изменения температуры
    pub fn derive_temperature_rate_from_column(
        self,
        source_col: &str,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let start = std::time::Instant::now();
        let mut updated = self.derive_rate(source_col, new_col, ColumnNature::TemperatureRate)?;
        updated.schema.dT_dt = Some(new_col.to_string());
        info!(
            "Derived temperature rate in {:?} ms",
            start.elapsed().as_millis()
        );
        Ok(updated)
    }

    /// Вычисление скорости изменения температуры
    pub fn derive_temperature_rate(self, new_col: &str) -> Result<Self, TGADomainError> {
        let temp = self
            .schema
            .temperature
            .as_ref()
            .ok_or(TGADomainError::TemperatureNotBound)?
            .clone();
        self.derive_temperature_rate_from_column(&temp, new_col)
    }
    /// Вычисление безразмерной скорости изменения
    pub fn derive_dimensionless_rate(
        self,
        col_name: &str,
        out_name: &str,
    ) -> Result<Self, TGADomainError> {
        let out_meta = ColumnMeta::new(
            out_name.to_string(),
            Unit::PerSecond,
            ColumnOrigin::PolarsDerived,
            ColumnNature::DimensionlessMassRate,
            ColumnProvenance::manual(
                out_name,
                "derive_dimensionless_rate",
                Some(format!("derived from {}", col_name)),
            ),
        );
        let (mut updated, is_conversion_rate) = if self.schema.eta.is_some() {
            (self.derive_rate0(col_name, out_name, out_meta)?, true)
        } else if self.schema.alpha.is_some() {
            (self.derive_rate0(col_name, out_name, out_meta)?, false)
        } else {
            return Err(TGADomainError::NoDimensionless);
        };
        if is_conversion_rate {
            updated.schema.deta_dt = Some(out_name.to_string());
        } else {
            updated.schema.dalpha_dt = Some(out_name.to_string());
        }
        Ok(updated)
    }

    /// `deta/dt` for conversion `eta = 1 - alpha`.
    pub fn derive_deta_dt(self, out_name: &str) -> Result<Self, TGADomainError> {
        let start = Instant::now();
        let out_meta = ColumnMeta::new(
            out_name.to_string(),
            Unit::PerSecond,
            ColumnOrigin::PolarsDerived,
            ColumnNature::ConversionRate,
            ColumnProvenance::manual(
                out_name,
                "derive_deta_dt",
                Some("derived from eta".to_string()),
            ),
        );
        let eta_col = self
            .schema
            .eta
            .as_ref()
            .ok_or(TGADomainError::EtaNotFound)?
            .clone();
        let mut updated = self.derive_rate0(&eta_col, out_name, out_meta)?;
        updated.schema.deta_dt = Some(out_name.to_string());
        info!("Derived deta_dt in {:?} ms", start.elapsed().as_millis());
        Ok(updated)
    }

    /// `dalpha/dt` for relative mass `alpha = m / m0`.
    pub fn derive_dalpha_dt(self, out_name: &str) -> Result<Self, TGADomainError> {
        let start = Instant::now();
        let out_meta = ColumnMeta::new(
            out_name.to_string(),
            Unit::PerSecond,
            ColumnOrigin::PolarsDerived,
            ColumnNature::DimensionlessMassRate,
            ColumnProvenance::manual(
                out_name,
                "derive_dalpha_dt",
                Some("derived from alpha".to_string()),
            ),
        );
        let alpha_col = self
            .schema
            .alpha
            .as_ref()
            .ok_or(TGADomainError::AlphaNotFound)?
            .clone();
        let mut updated = self.derive_rate0(&alpha_col, out_name, out_meta)?;
        updated.schema.dalpha_dt = Some(out_name.to_string());
        info!("Derived dalpha_dt in {:?} ms", start.elapsed().as_millis());
        Ok(updated)
    }

    pub fn drop_nulls(&mut self) -> Result<(), TGADomainError> {
        self.push_undo_snapshot();
        self.frame = self.frame.clone().drop_nulls(None);
        self.log_operation(
            "drop_nulls",
            AffectedColumns::All,
            None,
            "Dropped null rows from the dataset".to_string(),
            true,
        );
        Ok(())
    }

    pub fn get_column_by_nature(
        &self,
        nature: ColumnNature,
    ) -> Result<Option<String>, TGADomainError> {
        let schema_binding = match nature {
            ColumnNature::Time => Some(self.schema.time.as_deref()),
            ColumnNature::Temperature => Some(self.schema.temperature.as_deref()),
            ColumnNature::Mass => Some(self.schema.mass.as_deref()),
            ColumnNature::DimensionlessMass => Some(self.schema.alpha.as_deref()),
            ColumnNature::Conversion => Some(self.schema.eta.as_deref()),
            ColumnNature::MassRate => Some(self.schema.dm_dt.as_deref()),
            ColumnNature::ConversionRate => Some(self.schema.deta_dt.as_deref()),
            ColumnNature::TemperatureRate => Some(self.schema.dT_dt.as_deref()),
            ColumnNature::DimensionlessMassRate => Some(self.schema.dalpha_dt.as_deref()),
            ColumnNature::ActivationEnergy => Some(self.schema.E.as_deref()),
            ColumnNature::R2 => Some(self.schema.R2.as_deref()),
            ColumnNature::PredexFactor | ColumnNature::FittedCurve | ColumnNature::Unknown => None,
        };

        if let Some(binding) = schema_binding {
            return Ok(binding.map(str::to_owned));
        }

        let matches = self
            .schema
            .columns
            .values()
            .filter(|meta| meta.nature == nature)
            .map(|meta| meta.name.clone())
            .collect::<Vec<_>>();
        match matches.as_slice() {
            [] => Ok(None),
            [name] => Ok(Some(name.clone())),
            _ => Err(TGADomainError::AmbiguousColumnNature {
                nature,
                columns: matches,
            }),
        }
    }
    pub fn get_columns_by_nature(
        &self,
        nature: Vec<ColumnNature>,
    ) -> Result<Vec<Option<String>>, TGADomainError> {
        nature
            .into_iter()
            .map(|nature| self.get_column_by_nature(nature))
            .collect()
    }
}

#[cfg(test)]
mod schema_constructor_tests {
    use super::{ColumnNature, ColumnOrigin, TGASchema, Unit};

    #[test]
    fn from_columns_rejects_mismatched_name_and_unit_counts() {
        let err = TGASchema::from_columns(vec!["time".into()], Vec::<Unit>::new())
            .err()
            .unwrap();
        assert!(err.to_string().contains("1 names but 0 units"));
    }

    #[test]
    fn from_columns_infers_unique_roles_from_headers_not_positions() {
        let schema = TGASchema::from_columns(
            vec!["mass".into(), "aux".into(), "T".into(), "time".into()],
            vec![Unit::Milligram, Unit::Unknown, Unit::Celsius, Unit::Second],
        )
        .unwrap();

        assert_eq!(schema.time.as_deref(), Some("time"));
        assert_eq!(schema.temperature.as_deref(), Some("T"));
        assert_eq!(schema.mass.as_deref(), Some("mass"));
        assert_eq!(schema.columns["mass"].nature, ColumnNature::Mass);
        assert_eq!(schema.columns["T"].nature, ColumnNature::Temperature);
        assert_eq!(schema.columns["aux"].nature, ColumnNature::Unknown);
        assert!(
            schema
                .columns
                .values()
                .all(|meta| meta.origin == ColumnOrigin::Imported)
        );
    }

    #[test]
    fn from_columns_leaves_ambiguous_roles_unbound() {
        let schema = TGASchema::from_columns(
            vec![
                "time".into(),
                "t".into(),
                "temperature".into(),
                "mass".into(),
            ],
            vec![Unit::Second, Unit::Second, Unit::Celsius, Unit::Milligram],
        )
        .unwrap();

        assert_eq!(schema.time, None);
        assert_eq!(schema.columns["time"].nature, ColumnNature::Unknown);
        assert_eq!(schema.columns["t"].nature, ColumnNature::Unknown);
        assert_eq!(schema.temperature.as_deref(), Some("temperature"));
        assert_eq!(schema.mass.as_deref(), Some("mass"));
    }
}
//==============================================================
//numerical layer
#[derive(Debug)]
pub struct NumericBlock {
    pub grid: Arc<Vec<f64>>, // time or temperature
    pub columns: HashMap<String, Arc<Vec<f64>>>,
}

/// A row-aligned, nullable preview for UI tables. Unlike `materialize`, this
/// preserves source row indices and does not discard incomplete rows.
#[derive(Debug, Clone, PartialEq)]
pub struct PreviewRows {
    pub row_indices: Vec<usize>,
    pub columns: HashMap<String, Vec<Option<f64>>>,
}

impl TGADataset {
    /// Samples the same source rows from every requested column for display.
    pub fn preview_rows(&self, columns: &[&str], max_rows: usize) -> PolarsResult<PreviewRows> {
        let expressions: Vec<_> = columns.iter().map(|name| col(*name)).collect();
        let df = self.frame.clone().select(expressions).collect()?;
        let row_count = df.height();
        let step = if max_rows == 0 {
            row_count.max(1)
        } else {
            (row_count as f64 / max_rows as f64).ceil().max(1.0) as usize
        };
        let row_indices: Vec<_> = (0..row_count).step_by(step).collect();
        let mut sampled = HashMap::with_capacity(columns.len());
        for name in columns {
            let series = df.column(name)?.f64()?;
            let values = row_indices.iter().map(|&row| series.get(row)).collect();
            sampled.insert((*name).to_string(), values);
        }
        Ok(PreviewRows {
            row_indices,
            columns: sampled,
        })
    }

    /// Преобразование данных в численный формат для дальнейшей обработки
    pub fn materialize(&self, grid_col: &str, cols: &[&str]) -> PolarsResult<NumericBlock> {
        let mut requested = HashSet::with_capacity(cols.len() + 1);
        requested.insert(grid_col);
        for name in cols {
            if !requested.insert(*name) {
                return Err(PolarsError::ComputeError(
                    format!("Column '{name}' was requested more than once for materialization")
                        .into(),
                ));
            }
        }

        let mut select_exprs = Vec::with_capacity(cols.len() + 1);
        select_exprs.push(col(grid_col));
        select_exprs.extend(cols.iter().map(|name| col(*name)));

        let df = self.frame.clone().select(select_exprs).collect()?;
        let grid_values = df.column(grid_col)?.f64()?;
        let value_series = cols
            .iter()
            .map(|name| Ok((*name, df.column(name)?.f64()?)))
            .collect::<PolarsResult<Vec<_>>>()?;
        let mut grid = Vec::with_capacity(df.height());
        let mut columns = cols
            .iter()
            .map(|name| (name.to_string(), Vec::with_capacity(df.height())))
            .collect::<HashMap<_, _>>();
        let mut row_values = Vec::with_capacity(cols.len());
        let mut previous_time = None;

        // Apply one complete-case mask to the grid and every requested column.
        for row in 0..df.height() {
            let Some(x) = grid_values.get(row) else {
                continue;
            };
            row_values.clear();
            let mut complete = true;
            for (_, series) in &value_series {
                match series.get(row) {
                    Some(value) => row_values.push(value),
                    None => {
                        complete = false;
                        break;
                    }
                }
            }
            if !complete {
                continue;
            }

            if !x.is_finite() {
                return Err(PolarsError::ComputeError(
                    format!("Grid column '{grid_col}' has a non-finite value at row {row}").into(),
                ));
            }
            for ((name, _), value) in value_series.iter().zip(&row_values) {
                if !value.is_finite() {
                    return Err(PolarsError::ComputeError(
                        format!("Column '{name}' has a non-finite value at row {row}").into(),
                    ));
                }
            }

            if self.schema.time.as_deref() == Some(grid_col) {
                if let Some((previous_row, previous_x)) = previous_time {
                    if x <= previous_x {
                        return Err(PolarsError::ComputeError(
                            format!(
                                "Bound time column '{grid_col}' must be strictly increasing after null rows are omitted; rows {previous_row} and {row} contain {previous_x} then {x}"
                            )
                            .into(),
                        ));
                    }
                }
                previous_time = Some((row, x));
            }

            grid.push(x);
            for ((name, _), value) in value_series.iter().zip(&row_values) {
                columns
                    .get_mut(*name)
                    .ok_or_else(|| {
                        PolarsError::ComputeError(
                            format!("Missing output buffer for column '{name}'").into(),
                        )
                    })?
                    .push(*value);
            }
        }

        Ok(NumericBlock {
            grid: Arc::new(grid),
            columns: columns
                .into_iter()
                .map(|(name, data)| (name, Arc::new(data)))
                .collect(),
        })
    }
}
//======================================================================
/// numeical column representation
pub struct NumericColumn {
    pub data: Arc<Vec<f64>>,
}

impl NumericColumn {
    /*
    /// Преобразование в ndarray представление
    pub fn as_ndarray(&self) -> ArrayView1<f64> {
        ArrayView1::from(self.data.as_slice())
    }

    /// Преобразование в DVector представление
    pub fn as_dvector(&self) -> DVectorView<f64> {
        DVectorView::from_slice(self.data.as_slice(), self.data.len())
    }
    */
}
/// return data back to Polars
impl TGADataset {
    /// Добавление числовой колонки в набор данных
    pub fn add_numeric_column(
        mut self,
        name: &str,
        unit: Unit,
        data: Arc<Vec<f64>>,
        nature: ColumnNature,
    ) -> PolarsResult<Self> {
        self.ensure_column_name_is_available(name)
            .map_err(|err| PolarsError::ComputeError(format!("{:?}", err).into()))?;
        let series = Series::new(name.into(), data.as_ref());
        self.push_undo_snapshot();
        self.frame = self.frame.with_column(lit(series));

        self.schema.columns.insert(
            name.into(),
            ColumnMeta::new(
                name.to_string(),
                unit,
                ColumnOrigin::NumericDerived,
                nature,
                ColumnProvenance::manual(
                    name,
                    "add_numeric_column",
                    Some("external numeric data".to_string()),
                ),
            ),
        );

        self.log_operation(
            "add_numeric_column",
            AffectedColumns::Specific(vec![name.to_string()]),
            None,
            format!("Added numeric column {} with unit {:?}", name, unit),
            true,
        );

        Ok(self)
    }
}
