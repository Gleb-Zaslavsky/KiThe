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
/// # Recommended TGA data processing pipeline
///
/// KiThe follows an explicit, step-by-step processing model.
/// No operation implicitly removes rows or modifies unrelated columns.
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
/// dimensionless mass / conversion
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
use nalgebra::DVectorView;
use ndarray::ArrayView1;
use polars::error::{PolarsError, PolarsResult};

use polars::prelude::LazyFrame;
use polars::prelude::*;
use polars::series::Series;

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

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
    EmptyColumn(String),
    NonPositiveValues { column: String },
    ColumnNotFound(String),
    MassNotBound,
    TimeNotBound,
    NoDimensionless,
    AlphaNotFound,
    InvalidColumnType(String),
    TemperatureNotBound,
    DtDtNotFound,
    DmDtNotFound,
    DalphaDtNotFound,

    EtaNotFound,
    DetaDtNotFound,

    UnsupportedRateUnit,
    InvalidReferenceMass,
    InvalidUnitForConversion,
    NotImplemented,
    InvalidOperation(String),
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
    Hour,

    Kelvin,
    Celsius,

    MilliVolt,
    Milligram,
    MilligramPerSecond,
    KelvinPerSecond,
    CelsiusPerSecond,
    PerSecond,
    Gram,

    Dimensionless,
    Unknown,
}

impl Unit {
    pub fn parse(s: &str) -> Result<Self, PolarsError> {
        match s {
            "s" => Ok(Unit::Second),
            "h" => Ok(Unit::Hour),
            "C" => Ok(Unit::Celsius),
            "K" => Ok(Unit::Kelvin),
            "mV" => Ok(Unit::MilliVolt),
            "mg" => Ok(Unit::Milligram),
            "mg/s" => Ok(Unit::MilligramPerSecond),
            "K/s" => Ok(Unit::KelvinPerSecond),
            "C/s" => Ok(Unit::CelsiusPerSecond),
            "/s" => Ok(Unit::PerSecond),
            "1/s" => Ok(Unit::PerSecond),
            "g" => Ok(Unit::Gram),
            "-" => Ok(Unit::Dimensionless),

            _ => Err(PolarsError::ComputeError(
                format!("Unknown unit: {}", s).into(),
            )),
        }
    }
    pub fn convert_to_string(&self) -> String {
        match self {
            Unit::Second => "s".to_string(),
            Unit::Hour => "h".to_string(),
            Unit::Kelvin => "K".to_string(),
            Unit::Celsius => "C".to_string(),
            Unit::MilliVolt => "mV".to_string(),
            Unit::Milligram => "mg".to_string(),
            Unit::MilligramPerSecond => "mg/s".to_string(),
            Unit::KelvinPerSecond => "K/s".to_string(),
            Unit::CelsiusPerSecond => "C/s".to_string(),
            Unit::PerSecond => "/s".to_string(),
            Unit::Gram => "g".to_string(),
            Unit::Dimensionless => "-".to_string(),
            Unit::Unknown => "-".to_string(),
        }
    }
}

impl std::fmt::Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unit::Second => write!(f, "s"),
            Unit::Hour => write!(f, "h"),
            Unit::Kelvin => write!(f, "K"),
            Unit::Celsius => write!(f, "C"),
            Unit::MilliVolt => write!(f, "mV"),
            Unit::Milligram => write!(f, "mg"),
            Unit::MilligramPerSecond => write!(f, "mg/s"),
            Unit::KelvinPerSecond => write!(f, "K/s"),
            Unit::CelsiusPerSecond => write!(f, "C/s"),
            Unit::PerSecond => write!(f, "/s"),
            Unit::Gram => write!(f, "g"),
            Unit::Dimensionless => write!(f, "-"),
            Unit::Unknown => write!(f, "-"),
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
/// Columns Origin
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
}

impl TGASchema {
    fn get(&self, name: &str) -> Option<&ColumnMeta> {
        self.columns.get(name)
    }

    fn from_columns(names: Vec<String>, units: Vec<Unit>) -> Self {
        let mut columns = HashMap::new();
        for (i, name) in names.iter().enumerate() {
            let origin = if i < 3 {
                ColumnOrigin::Raw
            } else {
                ColumnOrigin::Imported
            };
            columns.insert(
                name.clone(),
                ColumnMeta {
                    name: name.clone(),
                    unit: units[i],
                    origin,
                },
            );
        }
        Self {
            columns,
            time: Some(names[0].clone()),
            temperature: Some(names[1].clone()),
            mass: Some(names[2].clone()),
            alpha: None,
            dm_dt: None,
            eta: None,
            deta_dt: None,
            dalpha_dt: None,
            dT_dt: None,
        }
    }
}
//================================================================================
//              TGADATASET - MAIN DATA STRUCTURE
//==================================================================================
/// Main entity Plolars-layer
#[derive(Clone)]
pub struct TGADataset {
    pub frame: LazyFrame,
    pub schema: TGASchema,
    pub oneframeplot: Option<OneFramePlot>,
}

impl std::fmt::Debug for TGADataset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TGADataset")
            .field("schema", &DebugSchema(&self.schema))
            .field("frame", &"<LazyFrame>")
            .finish()
    }
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
    //=========================================================================
    // INPUT/OUTPUT
    /// Normalize txt file to csv
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
    pub fn from_csv(
        path: &str,
        time_col: &str,
        temp_col: &str,
        mass_col: &str,
    ) -> PolarsResult<Self> {
        let path_buf = PathBuf::from(path);
        let path_arc = Arc::<Path>::from(path_buf.as_path());
        let plpath = PlPath::Local(path_arc);
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
            time_col.into(),
            ColumnMeta {
                name: time_col.into(),
                unit: Unit::Second,
                origin: ColumnOrigin::Raw,
            },
        );

        columns.insert(
            temp_col.into(),
            ColumnMeta {
                name: temp_col.into(),
                unit: Unit::Celsius,
                origin: ColumnOrigin::Raw,
            },
        );

        columns.insert(
            mass_col.into(),
            ColumnMeta {
                name: mass_col.into(),
                unit: Unit::Milligram,
                origin: ColumnOrigin::Raw,
            },
        );

        Ok(Self {
            frame,
            schema: TGASchema {
                columns,
                time: Some(time_col.into()),
                temperature: Some(temp_col.into()),
                mass: Some(mass_col.into()),
                alpha: None,
                dm_dt: None,
                eta: None,
                deta_dt: None,
                dalpha_dt: None,
                dT_dt: None,
            },
            oneframeplot: None,
        })
    }

    pub fn from_csv_raw(path: &Path) -> PolarsResult<Self> {
        let plpath = PlPath::Local(Arc::from(path));

        let mut frame = LazyCsvReader::new(plpath).with_has_header(true).finish()?;

        let schema = frame.collect_schema()?;

        let mut columns = HashMap::new();
        for name in schema.iter_names() {
            columns.insert(
                name.to_string(),
                ColumnMeta {
                    name: name.to_string(),
                    unit: Unit::Unknown,
                    origin: ColumnOrigin::Raw,
                },
            );
        }

        Ok(Self {
            frame,
            schema: TGASchema {
                columns,
                time: None,
                temperature: None,
                mass: None,
                alpha: None,
                dm_dt: None,
                eta: None,
                deta_dt: None,
                dalpha_dt: None,
                dT_dt: None,
            },
            oneframeplot: None,
        })
    }

    pub fn bind_time(mut self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        let meta = self
            .schema
            .columns
            .get_mut(col)
            .ok_or(TGADomainError::ColumnNotFound(col.into()))?;

        meta.unit = unit;
        self.schema.time = Some(col.into());
        Ok(self)
    }

    pub fn bind_temperature(mut self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        let meta = self
            .schema
            .columns
            .get_mut(col)
            .ok_or(TGADomainError::ColumnNotFound(col.into()))?;

        meta.unit = unit;
        self.schema.temperature = Some(col.into());
        Ok(self)
    }

    pub fn bind_mass(mut self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        let meta = self
            .schema
            .columns
            .get_mut(col)
            .ok_or(TGADomainError::ColumnNotFound(col.into()))?;

        meta.unit = unit;
        self.schema.mass = Some(col.into());
        Ok(self)
    }
    pub fn bind_column(
        mut self,
        role: ColumnRole,
        name: &str,
        unit: Unit,
    ) -> Result<Self, TGADomainError> {
        let meta = ColumnMeta {
            name: name.into(),
            unit,
            origin: ColumnOrigin::Raw,
        };

        self.schema.columns.insert(name.into(), meta);

        match role {
            ColumnRole::Time => self.schema.time = Some(name.into()),
            ColumnRole::Temperature => self.schema.temperature = Some(name.into()),
            ColumnRole::Mass => self.schema.mass = Some(name.into()),
        }

        Ok(self)
    }

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

        CsvWriter::new(&mut file).include_header(true).finish(df)?;

        Ok(())
    }

    /// Export selected columns to CSV with units encoded in a metadata header, similar to `to_csv_with_units`.
    /// Only the specified columns are exported, and units are taken from the schema.
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

        CsvWriter::new(&mut file)
            .include_header(true)
            .finish(&mut df)?;
        Ok(())
    }

    pub fn from_csv_with_units(path: &Path) -> PolarsResult<Self> {
        let text = std::fs::read_to_string(path)?;
        let mut lines = text.lines();

        lines.next(); // "# KiThe TGA Dataset"

        let meta = lines
            .next()
            .unwrap()
            .trim_start_matches("# ")
            .split(',')
            .map(|s| s.trim())
            .collect::<Vec<_>>();

        let mut names = Vec::new();
        let mut units = Vec::new();

        for entry in meta {
            let (name, unit) = entry.split_once('[').unwrap();
            names.push(name.trim().to_string());
            units.push(Unit::parse(unit.trim_end_matches(']'))?);
        }

        let data = lines.collect::<Vec<_>>().join("\n");
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), data)?;
        let tmp_path_arc = Arc::<Path>::from(tmp.path());
        let plpath = PlPath::Local(tmp_path_arc);
        let frame = LazyCsvReader::new(plpath).with_has_header(true).finish()?;
        std::mem::forget(tmp); // Keep temp file alive

        let schema = TGASchema::from_columns(names, units);

        Ok(Self {
            frame,
            schema,
            oneframeplot: None,
        })
    }
    /// There are 2 functions to read and parse from_csv_raw a function to read raw files with no metadata and headers and
    /// from_csv_with_units. This is a wrapper function from_csv_universal which recognizes if there in the file headers and call
    /// from_csv_with_units otherwise from_csv_raw
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
    pub fn derive_rate0(
        mut self,
        source_col: &str,
        new_col: &str,
        out_meta: ColumnMeta,
    ) -> Result<Self, TGADomainError> {
        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        let dv = col(source_col).shift(lit(-1)) - col(source_col).shift(lit(1));
        let dt = col(&time).shift(lit(-1)) - col(time).shift(lit(1));
        let rate = (dv / dt).alias(new_col);
        self.frame = self.frame.with_column(rate);

        self.schema.columns.insert(new_col.into(), out_meta);

        Ok(self)
    }
    pub fn derive_rate(mut self, source_col: &str, new_col: &str) -> Result<Self, TGADomainError> {
        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        let src_meta = self
            .schema
            .columns
            .get(source_col)
            // .or_else(|| self.schema.mass.as_ref().filter(|m| m.name == source_col))
            .ok_or(TGADomainError::ColumnNotFound(source_col.into()))?;

        let out_unit = match src_meta.unit {
            Unit::Milligram => Unit::MilligramPerSecond,
            Unit::Dimensionless => Unit::PerSecond,
            Unit::Kelvin => Unit::KelvinPerSecond,
            Unit::Celsius => Unit::CelsiusPerSecond,
            _ => return Err(TGADomainError::UnsupportedRateUnit),
        };
        let dv = col(source_col).shift(lit(-1)) - col(source_col).shift(lit(1));
        let dt = col(&time).shift(lit(-1)) - col(time).shift(lit(1));
        let rate = (dv / dt).alias(new_col);
        self.frame = self.frame.with_column(rate);

        self.schema.columns.insert(
            new_col.into(),
            ColumnMeta {
                name: new_col.into(),
                unit: out_unit,
                origin: ColumnOrigin::PolarsDerived,
            },
        );

        Ok(self)
    }

    pub fn derive_mass_rate(mut self, new_col: &str) -> Result<Self, TGADomainError> {
        let mass = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();
        self.schema.dm_dt = Some(new_col.to_string());
        self.derive_rate(&mass, new_col)
    }

    pub fn derive_temperature_rate(mut self, new_col: &str) -> Result<Self, TGADomainError> {
        let temp = self
            .schema
            .temperature
            .as_ref()
            .ok_or(TGADomainError::TemperatureNotBound)?
            .clone();
        self.schema.dT_dt = Some(new_col.to_string());
        self.derive_rate(&temp, new_col)
    }
    pub fn derive_dimensionless_rate(
        mut self,
        col_name: &str,
        out_name: &str,
    ) -> Result<Self, TGADomainError> {
        let out_meta = ColumnMeta {
            name: out_name.into(),
            unit: Unit::PerSecond,
            origin: ColumnOrigin::PolarsDerived,
        };
        if let Some(_) = self.schema.eta {
            self.schema.deta_dt = Some(out_name.to_string());
        } else if let Some(_) = self.schema.alpha {
            self.schema.dalpha_dt = Some(out_name.to_string());
        } else {
            return Err(TGADomainError::NoDimensionless);
        };

        self.derive_rate0(col_name, out_name, out_meta)
    }

    pub fn derive_deta_dt(mut self, out_name: &str) -> Result<Self, TGADomainError> {
        let out_meta = ColumnMeta {
            name: out_name.into(),
            unit: Unit::PerSecond,
            origin: ColumnOrigin::PolarsDerived,
        };
        let eta_col = self
            .schema
            .eta
            .as_ref()
            .ok_or(TGADomainError::EtaNotFound)?
            .clone();
        self.schema.deta_dt = Some(out_name.to_string());
        self.derive_rate0(&eta_col, out_name, out_meta)
    }

    pub fn derive_dalpha_dt(mut self, out_name: &str) -> Result<Self, TGADomainError> {
        let out_meta = ColumnMeta {
            name: out_name.into(),
            unit: Unit::PerSecond,
            origin: ColumnOrigin::PolarsDerived,
        };
        let alpha_col = self
            .schema
            .alpha
            .as_ref()
            .ok_or(TGADomainError::AlphaNotFound)?
            .clone();
        self.schema.dalpha_dt = Some(out_name.to_string());
        self.derive_rate0(&alpha_col, out_name, out_meta)
    }
}
//==============================================================
//numerical layer
pub struct NumericBlock {
    pub grid: Arc<Vec<f64>>, // time or temperature
    pub columns: HashMap<String, Arc<Vec<f64>>>,
}

impl TGADataset {
    pub fn materialize(&self, grid_col: &str, cols: &[&str]) -> PolarsResult<NumericBlock> {
        let mut select_exprs = Vec::with_capacity(cols.len() + 1);
        select_exprs.push(col(grid_col));

        for &c in cols {
            select_exprs.push(col(c));
        }

        let df = self.frame.clone().select(select_exprs).collect()?;

        let grid = Arc::new(df.column(grid_col)?.f64()?.into_no_null_iter().collect());

        let mut columns = HashMap::new();
        for &c in cols {
            let data = Arc::new(df.column(c)?.f64()?.into_no_null_iter().collect());
            columns.insert(c.to_string(), data);
        }

        Ok(NumericBlock { grid, columns })
    }
}
//======================================================================
/// numeical column representation
pub struct NumericColumn {
    pub data: Arc<Vec<f64>>,
}

impl NumericColumn {
    pub fn as_ndarray(&self) -> ArrayView1<f64> {
        ArrayView1::from(self.data.as_slice())
    }

    pub fn as_dvector(&self) -> DVectorView<f64> {
        DVectorView::from_slice(self.data.as_slice(), self.data.len())
    }
}
/// return data back to Polars
impl TGADataset {
    pub fn add_numeric_column(
        mut self,
        name: &str,
        unit: Unit,
        data: Arc<Vec<f64>>,
    ) -> PolarsResult<Self> {
        let series = Series::new(name.into(), data.as_ref());
        self.frame = self.frame.with_column(lit(series));

        self.schema.columns.insert(
            name.into(),
            ColumnMeta {
                name: name.into(),
                unit,
                origin: ColumnOrigin::NumericDerived,
            },
        );

        Ok(self)
    }
}
