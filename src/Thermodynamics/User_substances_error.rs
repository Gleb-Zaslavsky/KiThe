//! Error handling for SubsData operations
//!
//! This module provides comprehensive error handling for substance data operations,
//! including proper error propagation from thermo_api and transport_api modules.

use crate::Thermodynamics::DBhandlers::Diffusion::PairDiffusionError;
use crate::Thermodynamics::DBhandlers::thermo_api::ThermoError;
use crate::Thermodynamics::DBhandlers::transport_api::TransportError;
use std::error::Error;
use std::fmt;

/// Comprehensive error type for SubsData operations
#[derive(Debug)]
pub enum SubsDataError {
    /// Errors from thermodynamic calculations
    ThermoError(ThermoError),
    /// Errors from transport property calculations
    TransportError(TransportError),
    /// Substance not found in search results
    SubstanceNotFound(String),
    /// Calculator not available for substance
    CalculatorNotAvailable {
        substance: String,
        calc_type: String,
    },
    /// Missing required data
    MissingData { field: String, substance: String },
    /// Invalid temperature
    InvalidTemperature(f64),
    /// Pressure not set when required
    PressureNotSet,
    ///
    TemperatureNotSet,
    /// Molar mass not found for substance
    MolarMassNotFound(String),
    /// Heat capacity not available
    HeatCapacityNotAvailable(String),
    /// Library not found
    LibraryNotFound(String),
    /// Calculation failed
    CalculationFailed {
        substance: String,
        operation: String,
        reason: String,
    },
    /// Coefficient extraction failed
    CoefficientExtractionFailed {
        substance: String,
        temperature: Option<f64>,
    },
    /// Function creation failed
    FunctionCreationFailed {
        substance: String,
        function_type: String,
    },
    /// Symbolic expression creation failed
    SymbolicCreationFailed {
        substance: String,
        expression_type: String,
    },
    /// Matrix operation failed
    MatrixOperationFailed(String),
    /// NIST data retrieval failed
    NistRetrievalFailed { substance: String, reason: String },
}

impl fmt::Display for SubsDataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SubsDataError::ThermoError(e) => write!(f, "Thermodynamic calculation error: {}", e),
            SubsDataError::TransportError(e) => {
                write!(f, "Transport property calculation error: {}", e)
            }
            SubsDataError::SubstanceNotFound(substance) => {
                write!(f, "Substance '{}' not found in search results", substance)
            }
            SubsDataError::CalculatorNotAvailable {
                substance,
                calc_type,
            } => {
                write!(
                    f,
                    "Calculator '{}' not available for substance '{}'",
                    calc_type, substance
                )
            }
            SubsDataError::MissingData { field, substance } => {
                write!(
                    f,
                    "Missing required data '{}' for substance '{}'",
                    field, substance
                )
            }
            SubsDataError::InvalidTemperature(temp) => {
                write!(f, "Invalid temperature: {} K", temp)
            }

            SubsDataError::PressureNotSet => {
                write!(f, "Pressure must be set before performing this operation")
            }
            SubsDataError::TemperatureNotSet => {
                write!(
                    f,
                    "Temperature must be set before performing this operation"
                )
            }
            SubsDataError::MolarMassNotFound(substance) => {
                write!(f, "Molar mass not found for substance '{}'", substance)
            }
            SubsDataError::HeatCapacityNotAvailable(substance) => {
                write!(
                    f,
                    "Heat capacity not available for substance '{}'",
                    substance
                )
            }
            SubsDataError::LibraryNotFound(library) => {
                write!(f, "Library '{}' not found", library)
            }
            SubsDataError::CalculationFailed {
                substance,
                operation,
                reason,
            } => {
                write!(
                    f,
                    "Calculation '{}' failed for substance '{}': {}",
                    operation, substance, reason
                )
            }
            SubsDataError::CoefficientExtractionFailed {
                substance,
                temperature,
            } => {
                write!(
                    f,
                    "Failed to extract coefficients for substance '{}' at {:?} K",
                    substance, temperature
                )
            }
            SubsDataError::FunctionCreationFailed {
                substance,
                function_type,
            } => {
                write!(
                    f,
                    "Failed to create {} function for substance '{}'",
                    function_type, substance
                )
            }
            SubsDataError::SymbolicCreationFailed {
                substance,
                expression_type,
            } => {
                write!(
                    f,
                    "Failed to create {} symbolic expression for substance '{}'",
                    expression_type, substance
                )
            }
            SubsDataError::MatrixOperationFailed(reason) => {
                write!(f, "Matrix operation failed: {}", reason)
            }
            SubsDataError::NistRetrievalFailed { substance, reason } => {
                write!(
                    f,
                    "Failed to retrieve NIST data for substance '{}': {}",
                    substance, reason
                )
            }
        }
    }
}

impl Error for SubsDataError {}

impl From<ThermoError> for SubsDataError {
    fn from(err: ThermoError) -> Self {
        SubsDataError::ThermoError(err)
    }
}

impl From<TransportError> for SubsDataError {
    fn from(err: TransportError) -> Self {
        SubsDataError::TransportError(err)
    }
}
impl From<PairDiffusionError> for TransportError {
    fn from(err: PairDiffusionError) -> Self {
        match err {
            PairDiffusionError::InvalidTemperatureRange => TransportError::InvalidTemperatureRange,
            PairDiffusionError::SymbolicError(msg) => {
                TransportError::CalculationError(format!("Symbolic error: {}", msg))
            }
            PairDiffusionError::InvalidUnit(msg) => TransportError::InvalidUnit(msg),
            PairDiffusionError::SubstanceNotFound(msg) => TransportError::MissingData(msg),
        }
    }
}
impl From<PairDiffusionError> for SubsDataError {
    fn from(err: PairDiffusionError) -> Self {
        SubsDataError::TransportError(err.into())
    }
}
/// Result type for SubsData operations
pub type SubsDataResult<T> = Result<T, SubsDataError>;
////////////////////////////ERROR LOGGING/////////////////////////////////////////////////////////////////////////////
#[derive(Debug, Clone)]
pub enum LogErrorType {
    SubstanceNotFound(String),
    CalculatorNotAvailable(String),

    SerdeError(String),
    MissingData(String),
    InvalidUnit(String),
    InvalidTemperature(String),
    InvalidTemperatureRange(String),
    InvalidDensity(String),
    InvalidPressure(String),
    PressureNotSet(String),
    ParseError(String),

    TemperatureNotSet(String),
    MolarMassNotFound(String),
    HeatCapacityNotAvailable(String),
    LibraryNotFound(String),
    CalculationFailed(String),
    CoefficientExtractionFailed(String),
    FunctionCreationFailed(String),
    SymbolicCreationFailed(String),
    MatrixOperationFailed(String),
    NistRetrievalFailed(String),
    FittingError(String),
}
impl From<ThermoError> for LogErrorType {
    fn from(error: ThermoError) -> Self {
        match error {
            ThermoError::NoCoefficientsFound { temperature, range } => {
                LogErrorType::InvalidTemperatureRange(format!(
                    "No coefficients found for temperature {} K. Valid range: {}",
                    temperature, range
                ))
            }
            ThermoError::InvalidTemperatureRange => LogErrorType::InvalidTemperatureRange(
                "Invalid temperature range in coefficient data".to_string(),
            ),
            ThermoError::SerdeError(msg) => {
                LogErrorType::SerdeError(format!("Failed to deserialize NASA data: {}", msg))
            }
            ThermoError::UnsupportedUnit(unit) => LogErrorType::InvalidUnit(format!(
                "Unsupported unit: {}. Only 'J' and 'cal' are supported",
                unit
            )),
            ThermoError::SymbolicError(msg) => {
                LogErrorType::SerdeError(format!("Symbolic error: {}", msg))
            }
            ThermoError::CalculationError(msg) => {
                LogErrorType::CalculationFailed(format!("Calculation error: {}", msg))
            }
            ThermoError::FittingError(msg) => {
                LogErrorType::FittingError(format!("Fitting error: {}", msg))
            }
        }
    }
}
impl From<TransportError> for LogErrorType {
    fn from(error: TransportError) -> Self {
        match error {
            TransportError::MissingData(msg) => {
                LogErrorType::MissingData(format!("Missing required transport data: {}", msg))
            }
            TransportError::InvalidUnit(msg) => {
                LogErrorType::InvalidUnit(format!("Invalid  transport unit error: {}", msg))
            }
            TransportError::InvalidTemperature(msg) => LogErrorType::InvalidTemperature(format!(
                "Invalid temperature error ( transport): {}",
                msg
            )),
            TransportError::InvalidTemperatureRange => LogErrorType::InvalidTemperatureRange(
                "Invalid temperature range in transport data".to_string(),
            ),
            TransportError::CalculationError(msg) => {
                LogErrorType::CalculationFailed(format!("Transport calculation error: {}", msg))
            }
            TransportError::SerdeError(msg) => {
                LogErrorType::SerdeError(format!("Failed to deserialize transport data: {}", msg))
            }
            TransportError::InvalidDensity(msg) => {
                LogErrorType::InvalidDensity(format!("Invalid density error: {}", msg))
            }
            TransportError::InvalidPressure(msg) => LogErrorType::InvalidPressure(format!(
                "Invalid pressure error ( transport): {}",
                msg
            )),
            TransportError::InvalidMolarMass(msg) => {
                LogErrorType::MolarMassNotFound(format!("Invalid molar mass error: {}", msg))
            }
            TransportError::ParseError(msg) => {
                LogErrorType::ParseError(format!("Parse error: {}", msg))
            }
            TransportError::MissingCoefficients(msg) => LogErrorType::CoefficientExtractionFailed(
                format!("Missing coefficients error: {}", msg),
            ),
            TransportError::SymbolicError(msg) => {
                LogErrorType::SymbolicCreationFailed(format!("Symbolic error: {}", msg))
            }
            TransportError::FittingError(msg) => {
                LogErrorType::FittingError(format!("Fitting error: {}", msg))
            }
        }
    }
}

impl From<SubsDataError> for LogErrorType {
    fn from(error: SubsDataError) -> Self {
        match error {
            SubsDataError::ThermoError(ThermoError) => LogErrorType::from(ThermoError),
            SubsDataError::TransportError(TransportError) => LogErrorType::from(TransportError),
            SubsDataError::SubstanceNotFound(substance) => LogErrorType::SubstanceNotFound(
                format!("Substance '{}' not found in any library", substance),
            ),
            SubsDataError::CalculatorNotAvailable {
                substance,
                calc_type,
            } => LogErrorType::CalculatorNotAvailable(format!(
                "Calculator '{}' not available for substance '{}'",
                calc_type, substance
            )),
            SubsDataError::MissingData { field, substance } => LogErrorType::MissingData(format!(
                "Missing required data '{}' for substance '{}'",
                field, substance
            )),
            SubsDataError::InvalidTemperature(temp) => {
                LogErrorType::InvalidTemperature(format!("Invalid temperature: {} K", temp))
            }
            SubsDataError::PressureNotSet => LogErrorType::PressureNotSet(
                "Pressure must be set before performing this operation".to_string(),
            ),
            SubsDataError::TemperatureNotSet => LogErrorType::TemperatureNotSet(
                "Temperature must be set before performing this operation".to_string(),
            ),
            SubsDataError::MolarMassNotFound(substance) => LogErrorType::MolarMassNotFound(
                format!("Molar mass not found for substance '{}'", substance),
            ),
            SubsDataError::HeatCapacityNotAvailable(substance) => {
                LogErrorType::HeatCapacityNotAvailable(format!(
                    "Heat capacity not available for substance '{}'",
                    substance
                ))
            }
            SubsDataError::LibraryNotFound(library) => {
                LogErrorType::LibraryNotFound(format!("Library '{}' not found", library))
            }
            SubsDataError::CalculationFailed {
                substance,
                operation,
                reason,
            } => LogErrorType::CalculationFailed(format!(
                "Calculation '{}' failed for substance '{}': {}",
                operation, substance, reason
            )),
            SubsDataError::CoefficientExtractionFailed {
                substance,
                temperature,
            } => LogErrorType::CoefficientExtractionFailed(format!(
                "No coefficients found for substance '{}' for temperature {:?} K.",
                substance, temperature
            )),
            SubsDataError::FunctionCreationFailed {
                substance,
                function_type,
            } => LogErrorType::FunctionCreationFailed(format!(
                "Failed to create {} function for substance '{}'",
                function_type, substance
            )),
            SubsDataError::SymbolicCreationFailed {
                substance,
                expression_type,
            } => LogErrorType::SymbolicCreationFailed(format!(
                "Failed to create {} symbolic expression for substance '{}'",
                expression_type, substance
            )),
            SubsDataError::MatrixOperationFailed(reason) => {
                LogErrorType::MatrixOperationFailed(format!("Matrix operation failed: {}", reason))
            }
            SubsDataError::NistRetrievalFailed { substance, reason } => {
                LogErrorType::NistRetrievalFailed(format!(
                    "Failed to retrieve NIST data for substance '{}': {}",
                    substance, reason
                ))
            }
        }
    }
}

impl From<&ThermoError> for LogErrorType {
    fn from(error: &ThermoError) -> Self {
        match error {
            ThermoError::NoCoefficientsFound { temperature, range } => {
                LogErrorType::InvalidTemperatureRange(format!(
                    "No coefficients found for temperature {} K. Valid range: {}",
                    temperature, range
                ))
            }
            ThermoError::InvalidTemperatureRange => LogErrorType::InvalidTemperatureRange(
                "Invalid temperature range in coefficient data".to_string(),
            ),
            ThermoError::SerdeError(msg) => {
                LogErrorType::SerdeError(format!("Failed to deserialize NASA data: {}", msg))
            }
            ThermoError::UnsupportedUnit(unit) => LogErrorType::InvalidUnit(format!(
                "Unsupported unit: {}. Only 'J' and 'cal' are supported",
                unit
            )),
            ThermoError::SymbolicError(msg) => {
                LogErrorType::SerdeError(format!("Symbolic error: {}", msg))
            }
            ThermoError::CalculationError(msg) => {
                LogErrorType::CalculationFailed(format!("Calculation error: {}", msg))
            }
            ThermoError::FittingError(msg) => {
                LogErrorType::FittingError(format!("Fitting error: {}", msg))
            }
        }
    }
}
impl From<&TransportError> for LogErrorType {
    fn from(error: &TransportError) -> Self {
        match error {
            TransportError::MissingData(msg) => {
                LogErrorType::MissingData(format!("Missing required transport data: {}", msg))
            }
            TransportError::InvalidUnit(msg) => {
                LogErrorType::InvalidUnit(format!("Invalid  transport unit error: {}", msg))
            }
            TransportError::InvalidTemperature(msg) => LogErrorType::InvalidTemperature(format!(
                "Invalid temperature error ( transport): {}",
                msg
            )),
            TransportError::CalculationError(msg) => {
                LogErrorType::CalculationFailed(format!("Transport calculation error: {}", msg))
            }
            TransportError::SerdeError(msg) => {
                LogErrorType::SerdeError(format!("Failed to deserialize transport data: {}", msg))
            }
            TransportError::InvalidDensity(msg) => {
                LogErrorType::InvalidDensity(format!("Invalid density error: {}", msg))
            }
            TransportError::InvalidPressure(msg) => LogErrorType::InvalidPressure(format!(
                "Invalid pressure error ( transport): {}",
                msg
            )),
            TransportError::InvalidMolarMass(msg) => {
                LogErrorType::MolarMassNotFound(format!("Invalid molar mass error: {}", msg))
            }
            TransportError::ParseError(msg) => {
                LogErrorType::ParseError(format!("Parse error: {}", msg))
            }
            TransportError::MissingCoefficients(msg) => LogErrorType::CoefficientExtractionFailed(
                format!("Missing coefficients error: {}", msg),
            ),
            TransportError::InvalidTemperatureRange => LogErrorType::InvalidTemperatureRange(
                "Invalid temperature range in transport data".to_string(),
            ),
            TransportError::SymbolicError(msg) => {
                LogErrorType::SymbolicCreationFailed(format!("Symbolic error: {}", msg))
            }
            TransportError::FittingError(msg) => {
                LogErrorType::FittingError(format!("Fitting error: {}", msg))
            }
        }
    }
}
impl From<&SubsDataError> for LogErrorType {
    fn from(error: &SubsDataError) -> Self {
        match error {
            SubsDataError::ThermoError(thermo_error) => LogErrorType::from(thermo_error),
            SubsDataError::TransportError(transport_error) => LogErrorType::from(transport_error),
            SubsDataError::SubstanceNotFound(substance) => LogErrorType::SubstanceNotFound(
                format!("Substance '{}' not found in any library", substance),
            ),
            SubsDataError::CalculatorNotAvailable {
                substance,
                calc_type,
            } => LogErrorType::CalculatorNotAvailable(format!(
                "Calculator '{}' not available for substance '{}'",
                calc_type, substance
            )),
            SubsDataError::MissingData { field, substance } => LogErrorType::MissingData(format!(
                "Missing required data '{}' for substance '{}'",
                field, substance
            )),
            SubsDataError::InvalidTemperature(temp) => {
                LogErrorType::InvalidTemperature(format!("Invalid temperature: {} K", temp))
            }
            SubsDataError::PressureNotSet => LogErrorType::PressureNotSet(
                "Pressure must be set before performing this operation".to_string(),
            ),
            SubsDataError::TemperatureNotSet => LogErrorType::TemperatureNotSet(
                "Temperature must be set before performing this operation".to_string(),
            ),
            SubsDataError::MolarMassNotFound(substance) => LogErrorType::MolarMassNotFound(
                format!("Molar mass not found for substance '{}'", substance),
            ),
            SubsDataError::HeatCapacityNotAvailable(substance) => {
                LogErrorType::HeatCapacityNotAvailable(format!(
                    "Heat capacity not available for substance '{}'",
                    substance
                ))
            }
            SubsDataError::LibraryNotFound(library) => {
                LogErrorType::LibraryNotFound(format!("Library '{}' not found", library))
            }
            SubsDataError::CalculationFailed {
                substance,
                operation,
                reason,
            } => LogErrorType::CalculationFailed(format!(
                "Calculation '{}' failed for substance '{}': {}",
                operation, substance, reason
            )),
            SubsDataError::CoefficientExtractionFailed {
                substance,
                temperature,
            } => LogErrorType::CoefficientExtractionFailed(format!(
                "No coefficients found for substance '{}' for temperature {:?} K.",
                substance, temperature
            )),
            SubsDataError::FunctionCreationFailed {
                substance,
                function_type,
            } => LogErrorType::FunctionCreationFailed(format!(
                "Failed to create {} function for substance '{}'",
                function_type, substance
            )),
            SubsDataError::SymbolicCreationFailed {
                substance,
                expression_type,
            } => LogErrorType::SymbolicCreationFailed(format!(
                "Failed to create {} symbolic expression for substance '{}'",
                expression_type, substance
            )),
            SubsDataError::MatrixOperationFailed(reason) => {
                LogErrorType::MatrixOperationFailed(format!("Matrix operation failed: {}", reason))
            }
            SubsDataError::NistRetrievalFailed { substance, reason } => {
                LogErrorType::NistRetrievalFailed(format!(
                    "Failed to retrieve NIST data for substance '{}': {}",
                    substance, reason
                ))
            }
        }
    }
}

impl From<&PairDiffusionError> for LogErrorType {
    fn from(error: &PairDiffusionError) -> Self {
        match error {
            PairDiffusionError::InvalidTemperatureRange => LogErrorType::InvalidTemperatureRange(
                "Invalid temperature range in diffusion data".to_string(),
            ),

            PairDiffusionError::SymbolicError(msg) => {
                LogErrorType::SymbolicCreationFailed(format!("Symbolic error: {}", msg))
            }
            PairDiffusionError::InvalidUnit(msg) => {
                LogErrorType::InvalidUnit(format!("Invalid unit error: {}", msg))
            }
            PairDiffusionError::SubstanceNotFound(msg) => {
                LogErrorType::MissingData("Substance not found in diffusion library".to_string())
            }
        }
    }
}
#[derive(Debug, Clone)]
pub struct ExceptionLogEntry {
    pub substance: String,
    pub error_type: LogErrorType,
    pub function: String,
    pub message: String,
    pub timestamp: std::time::SystemTime,
}

pub trait ExceptionLogger {
    fn log_exception(&mut self, entry: ExceptionLogEntry);
    fn get_logs(&self) -> &Vec<ExceptionLogEntry>;
    fn clear_logs(&mut self);
    fn print_logs(&self); // Pretty print using table crate
}

#[derive(Debug, Clone)]
pub struct SimpleExceptionLogger {
    logs: Vec<ExceptionLogEntry>,
}

impl SimpleExceptionLogger {
    pub fn new() -> Self {
        Self { logs: Vec::new() }
    }

    pub fn print_pretty_table(&self) {
        use prettytable::{Cell, Row, Table};

        if self.logs.is_empty() {
            println!("No errors logged.");
            return;
        }

        let mut table = Table::new();
        table.add_row(Row::new(vec![
            Cell::new("Substance"),
            Cell::new("Error Type"),
            Cell::new("Function"),
            Cell::new("Message"),
        ]));

        for entry in &self.logs {
            let error_debug = format!("{:?}", entry.error_type);
            let error_type_str = error_debug.split('(').next().unwrap_or("Unknown");
            table.add_row(Row::new(vec![
                Cell::new(&entry.substance),
                Cell::new(error_type_str),
                Cell::new(&entry.function),
                Cell::new(&entry.message),
            ]));
        }

        table.printstd();
    }
}

impl ExceptionLogger for SimpleExceptionLogger {
    fn log_exception(&mut self, entry: ExceptionLogEntry) {
        self.logs.push(entry);
    }

    fn get_logs(&self) -> &Vec<ExceptionLogEntry> {
        &self.logs
    }

    fn clear_logs(&mut self) {
        self.logs.clear();
    }

    fn print_logs(&self) {
        for entry in &self.logs {
            println!(
                "`{:?}` for substance `{}` in `{}`: {}",
                entry.error_type, entry.substance, entry.function, entry.message
            );
        }
    }
}

#[macro_export]
macro_rules! log_error {
    ($logger:expr, $error:expr, $substance:expr, $function:expr) => {{
        use $crate::Thermodynamics::User_substances_error::{ExceptionLogEntry, ExceptionLogger};
        $logger.log_exception(ExceptionLogEntry {
            substance: $substance.to_string(),
            error_type: $error.into(),
            function: $function.to_string(),
            message: $error.to_string(),
            timestamp: std::time::SystemTime::now(),
        });
    }};
}
