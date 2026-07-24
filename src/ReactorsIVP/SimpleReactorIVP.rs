//! # Simple Reactor IVP Module
//!
//! This module provides reactor-domain helpers for the condensed-phase IVP workflow.
//! The production target is constant density, one state per species, a thermal
//! pair, and LSODE2-friendly setup helpers.
//!
//! ## Main Structures
//!
//! - **`SimpleReactorTask`**: Main reactor modeling structure that aggregates
//!   kinetics, thermodynamics, and condensed-phase reactor inputs
//! - **`IVPSolver`**: Wrapper for the assembled IVP equations and results
//! - **`FastElemReact`**: Simple structure for elementary reactions with
//!   Arrhenius parameters
//!
//! ## Key Features
//!
//! - **Constant-density condensed-phase modeling**: Density is supplied
//!   directly by the user
//! - **Typed initial state**: The canonical state is assembled from a single
//!   validated source
//! - **LSODE2-oriented setup**: Solver handoff stays typed and narrow
//! - **Reaction preprocessing**: Stoichiometry, rate expressions, and heat
//!   effects are prepared once
//!
//! ## Mathematical Model
//!
//! The module solves the condensed-phase reactor IVP in the physical
//! coordinate x with one state per species plus the thermal pair.
use crate::Kinetics::User_reactions::KinData;
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::ReactorsBVP::SimpleReactorBVP::ReactorError;
use crate::ReactorsBVP::reactor_BVP_utils::{
    ScalingConfig, create_bounds_map, create_tolerance_map,
};
use crate::ReactorsIVP::solver_backend::{
    ReactorIvpBackendError, ReactorIvpMethod, ReactorIvpSolverConfig,
};

use RustedSciThe::numerical::LSODE2::solver::Lsode2EvaluationTelemetry;
use RustedSciThe::numerical::LSODE2::{
    Lsode2AlgorithmSnapshot, Lsode2Error, Lsode2NativeStatistics, Lsode2ProblemConfig,
    Lsode2SolveSummary, Lsode2Solver,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_ivp_generated::IvpBackendStatistics;
use log::{info, warn};
use nalgebra::{DMatrix, DVector};
use prettytable::{Table as PrettyTable, row};
use tabled::{Table, Tabled};

use std::collections::HashMap;
use std::fmt;

/// Universal gas constant in J/(mol·K)

/// Canonical independent variable for the condensed-phase reactor IVP.
///
/// The condensed reactor model is a spatial problem in the physical coordinate
/// `x`, not a time-integration problem and not a legacy dimensionless `z`
/// variable.
pub const REACTOR_IVP_INDEPENDENT_VARIABLE: &str = "x";

/// Local IVP error type for condensed-phase helpers.
#[derive(Debug, Clone)]
pub enum IvpError {
    MissingData(String),
    InvalidConfiguration(String),
    InvalidNumericValue(String),
    CalculationError(String),
    IndexOutOfBounds(String),
}

impl fmt::Display for IvpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IvpError::MissingData(msg) => write!(f, "Missing data: {}", msg),
            IvpError::InvalidConfiguration(msg) => write!(f, "Invalid configuration: {}", msg),
            IvpError::InvalidNumericValue(msg) => write!(f, "Invalid numeric value: {}", msg),
            IvpError::CalculationError(msg) => write!(f, "Calculation error: {}", msg),
            IvpError::IndexOutOfBounds(msg) => write!(f, "Index out of bounds: {}", msg),
        }
    }
}

impl std::error::Error for IvpError {}

impl From<ReactorError> for IvpError {
    fn from(value: ReactorError) -> Self {
        match value {
            ReactorError::MissingData(msg) => Self::MissingData(msg),
            ReactorError::InvalidConfiguration(msg) => Self::InvalidConfiguration(msg),
            ReactorError::InvalidNumericValue(msg) => Self::InvalidNumericValue(msg),
            ReactorError::CalculationError(msg) => Self::CalculationError(msg),
            ReactorError::ParseError(msg) => Self::CalculationError(msg),
            ReactorError::IndexOutOfBounds(msg) => Self::IndexOutOfBounds(msg),
        }
    }
}

impl From<IvpError> for ReactorError {
    fn from(value: IvpError) -> Self {
        match value {
            IvpError::MissingData(msg) => Self::MissingData(msg),
            IvpError::InvalidConfiguration(msg) => Self::InvalidConfiguration(msg),
            IvpError::InvalidNumericValue(msg) => Self::InvalidNumericValue(msg),
            IvpError::CalculationError(msg) => Self::CalculationError(msg),
            IvpError::IndexOutOfBounds(msg) => Self::IndexOutOfBounds(msg),
        }
    }
}

impl From<ReactorIvpBackendError> for IvpError {
    fn from(value: ReactorIvpBackendError) -> Self {
        match value {
            ReactorIvpBackendError::InvalidConfiguration(msg) => {
                IvpError::InvalidConfiguration(msg)
            }
        }
    }
}

impl From<Lsode2Error> for IvpError {
    fn from(value: Lsode2Error) -> Self {
        match value {
            Lsode2Error::UnsupportedBackend(message) => {
                IvpError::InvalidConfiguration(message.to_string())
            }
            Lsode2Error::InvalidConfig(message) => IvpError::InvalidConfiguration(message),
            Lsode2Error::GeneratedBackend(err) => IvpError::CalculationError(err.to_string()),
            Lsode2Error::NativeStep(message) => IvpError::CalculationError(message),
        }
    }
}

/// Canonical state order for condensed-phase IVP snapshots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IvpCanonicalStateOrder {
    ThermalThenSpecies,
}

impl Default for IvpCanonicalStateOrder {
    fn default() -> Self {
        Self::ThermalThenSpecies
    }
}

/// Typed physical inputs for a condensed-phase reactor task.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorIvpPhysicalConfig {
    pub ro: f64,
    pub Cp: f64,
    pub Lambda: f64,
    pub m: f64,
    pub L: f64,
}

impl ReactorIvpPhysicalConfig {
    /// Validate the finite-positive physical inputs required by the condensed model.
    pub fn validate(&self) -> Result<(), IvpError> {
        let validate = |name: &str, value: f64| -> Result<(), IvpError> {
            if !value.is_finite() {
                return Err(IvpError::InvalidNumericValue(format!(
                    "{} must be finite, got {}",
                    name, value
                )));
            }
            if value <= 0.0 {
                return Err(IvpError::InvalidNumericValue(format!(
                    "{} must be positive, got {}",
                    name, value
                )));
            }
            Ok(())
        };

        validate("ro", self.ro)?;
        validate("Cp", self.Cp)?;
        validate("Lambda", self.Lambda)?;
        validate("m", self.m)?;
        validate("L", self.L)?;
        Ok(())
    }
}

/// Typed initial state for the condensed-phase reactor task.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorIvpInitialStateConfig {
    pub temperature: f64,
    pub heat_flux: f64,
    pub species: HashMap<String, f64>,
}

impl ReactorIvpInitialStateConfig {
    /// Validate the typed initial state against the required canonical entries.
    pub fn validate(&self, species_order: &[String]) -> Result<(), IvpError> {
        if !self.temperature.is_finite() {
            return Err(IvpError::InvalidNumericValue(format!(
                "Initial temperature must be finite, got {}",
                self.temperature
            )));
        }
        if !self.heat_flux.is_finite() {
            return Err(IvpError::InvalidNumericValue(format!(
                "Initial heat flux must be finite, got {}",
                self.heat_flux
            )));
        }
        for substance in species_order {
            let value = self.species.get(substance).ok_or_else(|| {
                IvpError::MissingData(format!(
                    "Missing initial condition for substance `{}`",
                    substance
                ))
            })?;
            if !value.is_finite() {
                return Err(IvpError::InvalidNumericValue(format!(
                    "Initial condition for `{}` must be finite, got {}",
                    substance, value
                )));
            }
        }
        Ok(())
    }
}

/// Typed input bundle for the condensed-phase IVP contract.
#[derive(Debug, Clone)]
pub struct ReactorIvpTaskInputs {
    pub physical: ReactorIvpPhysicalConfig,
    pub scaling: ScalingConfig,
    pub initial_state: ReactorIvpInitialStateConfig,
    pub thermal_effects: Vec<f64>,
}

impl ReactorIvpTaskInputs {
    /// Validate the coupled physical, scaling, and initial-state inputs.
    pub fn validate(&self, species_order: &[String]) -> Result<(), IvpError> {
        self.physical.validate()?;
        self.scaling.validate().map_err(IvpError::from)?;
        self.initial_state.validate(species_order)?;
        if self.thermal_effects.is_empty() {
            return Err(IvpError::MissingData(
                "No thermal effects provided for condensed IVP".to_string(),
            ));
        }
        Ok(())
    }
}

/// Immutable condensed-IVP solve snapshot published after a successful solve.
#[derive(Debug, Clone)]
pub struct IvpSolveSnapshot {
    pub axis: DVector<f64>,
    pub values: DMatrix<f64>,
    pub variable_names: Vec<String>,
    pub summary: Lsode2SolveSummary,
    pub backend_plan: ReactorIvpSolverConfig,
    pub algorithm: Lsode2AlgorithmSnapshot,
    pub statistics: IvpBackendStatistics,
    pub native_statistics: Lsode2NativeStatistics,
    pub evaluation_telemetry: Lsode2EvaluationTelemetry,
}

/// Lightweight solve report derived from a validated snapshot.
///
/// The report intentionally copies only diagnostics and high-level metadata so
/// callers can inspect the solve without cloning the full solution matrix.
#[derive(Debug, Clone)]
pub struct IvpSolveReport {
    pub status: String,
    pub method_label: String,
    pub time_points: usize,
    pub variable_count: usize,
    pub backend_plan: ReactorIvpSolverConfig,
    pub algorithm: Lsode2AlgorithmSnapshot,
    pub statistics: IvpBackendStatistics,
    pub native_statistics: Lsode2NativeStatistics,
    pub evaluation_telemetry: Lsode2EvaluationTelemetry,
}

/// Borrowed result view over a validated condensed-IVP solve snapshot.
///
/// The view keeps result inspection idempotent and avoids cloning the matrix
/// or the LSODE2 summary when the caller only needs read-only access.
#[derive(Debug, Clone, Copy)]
pub struct IvpResultView<'a> {
    snapshot: &'a IvpSolveSnapshot,
}

fn build_solution_preview_rows(
    variable_names: &[String],
    values: &DMatrix<f64>,
) -> Result<Vec<IvpSolutionPreviewRow>, IvpError> {
    if variable_names.len() != values.ncols() {
        return Err(IvpError::CalculationError(format!(
            "variable count {} does not match solution columns {}",
            variable_names.len(),
            values.ncols()
        )));
    }
    if values.nrows() == 0 || values.ncols() == 0 {
        return Err(IvpError::MissingData(
            "Solution matrix is empty and cannot be previewed".to_string(),
        ));
    }

    let mut rows = Vec::with_capacity(variable_names.len());
    for (index, variable) in variable_names.iter().enumerate() {
        let column = values.column(index);
        let last = column.len().saturating_sub(1);
        let middle = column.len() / 2;
        rows.push(IvpSolutionPreviewRow {
            index,
            variable: variable.clone(),
            first: column[0],
            middle: column[middle],
            last: column[last],
        });
    }
    Ok(rows)
}

/// A single key-value row for a diagnostics table.
///
/// Keeping diagnostics in row form makes it easy to feed the same data into
/// console logs, GUI previews, and regression tests.
#[derive(Debug, Clone, PartialEq, Tabled)]
pub struct IvpReportRow {
    pub field: String,
    pub value: String,
}

/// Snapshot of the condensed-phase conservation checks.
///
/// The report mirrors the BVP balance diagnostics closely enough that the IVP
/// path can expose the same kind of physical sanity check after a successful
/// solve.
#[derive(Debug, Clone, PartialEq)]
pub struct IvpConservationReport {
    pub energy_balance_error_abs: f64,
    pub energy_balance_error_rel: f64,
    pub sum_of_mass_fractions: Vec<(usize, f64)>,
    pub atomic_mass_balance_error: Vec<(usize, f64)>,
}

/// One row of a compact solution preview table.
///
/// The preview is intentionally data-only so callers can render it however they
/// want without mutating the stored solve snapshot.
#[derive(Debug, Clone, PartialEq, Tabled)]
pub struct IvpSolutionPreviewRow {
    pub index: usize,
    pub variable: String,
    pub first: f64,
    pub middle: f64,
    pub last: f64,
}

impl IvpSolveSnapshot {
    /// Create a borrowed read-only result view over this snapshot.
    pub fn result_view(&self) -> IvpResultView<'_> {
        IvpResultView { snapshot: self }
    }

    /// Return the number of mesh points stored in the snapshot.
    pub fn time_points(&self) -> usize {
        self.axis.len()
    }

    /// Return the number of state variables stored in the snapshot.
    pub fn variable_count(&self) -> usize {
        self.values.ncols()
    }

    /// Return the solver status string.
    pub fn status(&self) -> &str {
        &self.summary.status
    }

    /// Return the typed LSODE2 backend plan that produced this snapshot.
    pub fn backend_config(&self) -> ReactorIvpSolverConfig {
        self.backend_plan
    }

    /// Return the typed backend plan that produced this snapshot.
    pub fn resolved_backend_plan(&self) -> ReactorIvpSolverConfig {
        self.backend_plan
    }

    /// Return the resolved LSODE2 method family from the published plan.
    pub fn method(&self) -> ReactorIvpMethod {
        self.backend_plan.method
    }

    /// Return the resolved method family from the published summary.
    pub fn method_label(&self) -> &str {
        self.summary.method
    }

    /// Return the published algorithm snapshot.
    pub fn algorithm_snapshot(&self) -> &Lsode2AlgorithmSnapshot {
        &self.algorithm
    }

    /// Return the published backend statistics.
    pub fn backend_statistics(&self) -> &IvpBackendStatistics {
        &self.statistics
    }

    /// Return the published native LSODE2 statistics.
    pub fn native_statistics(&self) -> &Lsode2NativeStatistics {
        &self.native_statistics
    }

    /// Return the published evaluation telemetry.
    pub fn evaluation_telemetry(&self) -> &Lsode2EvaluationTelemetry {
        &self.evaluation_telemetry
    }

    /// Build compact preview rows for the stored solution matrix.
    ///
    /// The helper stays on the validated snapshot so every read path uses the
    /// same canonical metadata and solution data.
    pub fn solution_preview_rows(&self) -> Result<Vec<IvpSolutionPreviewRow>, IvpError> {
        build_solution_preview_rows(&self.variable_names, &self.values)
    }

    /// Build a small diagnostics report without cloning the solution matrix.
    pub fn diagnostics_report(&self) -> IvpSolveReport {
        IvpSolveReport {
            status: self.summary.status.clone(),
            method_label: self.summary.method.to_string(),
            time_points: self.axis.len(),
            variable_count: self.values.ncols(),
            backend_plan: self.backend_plan,
            algorithm: self.algorithm.clone(),
            statistics: self.statistics.clone(),
            native_statistics: self.native_statistics.clone(),
            evaluation_telemetry: self.evaluation_telemetry.clone(),
        }
    }

    /// Return a diagnostics table snapshot as plain rows.
    ///
    /// The rows are intentionally separate from the render step so callers can
    /// inspect or format them without coupling to `println!`.
    pub fn diagnostics_report_rows(&self) -> Vec<IvpReportRow> {
        let report = self.diagnostics_report();
        let algorithm = &report.algorithm;
        let stats = &report.statistics;
        let native = &report.native_statistics;
        let telemetry = &report.evaluation_telemetry;

        vec![
            IvpReportRow {
                field: "status".to_string(),
                value: report.status,
            },
            IvpReportRow {
                field: "method_label".to_string(),
                value: report.method_label,
            },
            IvpReportRow {
                field: "time_points".to_string(),
                value: report.time_points.to_string(),
            },
            IvpReportRow {
                field: "variable_count".to_string(),
                value: report.variable_count.to_string(),
            },
            IvpReportRow {
                field: "controller_mode".to_string(),
                value: algorithm.controller_mode.to_string(),
            },
            IvpReportRow {
                field: "active_family".to_string(),
                value: algorithm.active_family.to_string(),
            },
            IvpReportRow {
                field: "resolved_source".to_string(),
                value: format!("{:?}", report.backend_plan.execution_backend),
            },
            IvpReportRow {
                field: "resolved_structure".to_string(),
                value: format!("{:?}", report.backend_plan.matrix_backend),
            },
            IvpReportRow {
                field: "solve_calls".to_string(),
                value: stats.solve_calls.to_string(),
            },
            IvpReportRow {
                field: "residual_calls".to_string(),
                value: stats.residual_calls.to_string(),
            },
            IvpReportRow {
                field: "jacobian_calls".to_string(),
                value: stats.jacobian_calls.to_string(),
            },
            IvpReportRow {
                field: "native_step_attempts".to_string(),
                value: native.native_step_attempts.to_string(),
            },
            IvpReportRow {
                field: "accepted_steps".to_string(),
                value: telemetry.accepted_steps.to_string(),
            },
        ]
    }
}

impl IvpConservationReport {
    /// Return the report in a table-friendly row format.
    pub fn rows(&self) -> Vec<IvpReportRow> {
        vec![
            IvpReportRow {
                field: "absolute_energy_balance_error".to_string(),
                value: self.energy_balance_error_abs.to_string(),
            },
            IvpReportRow {
                field: "relative_energy_balance_error_percent".to_string(),
                value: self.energy_balance_error_rel.to_string(),
            },
            IvpReportRow {
                field: "mass_fraction_rows_off_1".to_string(),
                value: self.sum_of_mass_fractions.len().to_string(),
            },
            IvpReportRow {
                field: "atomic_balance_rows_off_0".to_string(),
                value: self.atomic_mass_balance_error.len().to_string(),
            },
        ]
    }
}

/// Owned render snapshot for the postprocessed condensed-IVP solution.
///
/// The solver keeps the validated solve snapshot for diagnostics, while this
/// view exposes the physical solution values that should be shown in plots,
/// tables, and exports after postprocessing.
#[derive(Debug, Clone)]
pub struct IvpSolutionRenderData {
    pub x_mesh: DVector<f64>,
    pub solution: DMatrix<f64>,
    pub unknowns: Vec<String>,
    pub arg_name: String,
}

impl IvpSolutionRenderData {
    /// Build compact preview rows for the rendered solution matrix.
    pub fn solution_preview_rows(&self) -> Result<Vec<IvpSolutionPreviewRow>, IvpError> {
        build_solution_preview_rows(&self.unknowns, &self.solution)
    }
}

/// Physical postprocessing output for the condensed-IVP solution.
///
/// This mirrors the BVP pattern: the raw solver state is validated first and
/// then converted to dimensional coordinates for user-facing output.
#[derive(Debug, Clone)]
pub struct IvpPostprocessingReport {
    pub x_mesh: DVector<f64>,
    pub solution: DMatrix<f64>,
}

impl<'a> IvpResultView<'a> {
    /// Return the underlying validated solve snapshot.
    pub fn snapshot(&self) -> &IvpSolveSnapshot {
        self.snapshot
    }

    /// Return diagnostics as a borrowed report snapshot.
    pub fn report(&self) -> IvpSolveReport {
        self.snapshot.diagnostics_report()
    }

    /// Return diagnostics rows for table rendering.
    pub fn report_rows(&self) -> Vec<IvpReportRow> {
        self.snapshot.diagnostics_report_rows()
    }

    /// Return compact preview rows for the solution matrix.
    pub fn solution_preview_rows(&self) -> Result<Vec<IvpSolutionPreviewRow>, IvpError> {
        self.snapshot.solution_preview_rows()
    }
}

/*
    Condensed-phase IVP model notes.

    The IVP path is spatial and uses the physical coordinate `x`.
    The canonical state is ordered as:
        [Teta, q, C0, C1, ...]

    Thermal pair:
        dTeta/dx = q / Lambda
        dq/dx = Pe_q * q - L^2 * Q / dT

    Units follow the condensed-phase task document:
        ro     - density [kg/m^3]
        Cp     - heat capacity [J/(kg*K)]
        Lambda - thermal conductivity [W/(m*K)]
        L      - characteristic length [m]
        m      - condensed-phase flow parameter / mass-flux-like input

    Sign convention:
        species source terms follow the stoichiometric analyzer contract.
        In the current canonical `A -> B` fixture, the assembled RHS evaluates
        with a positive `A` entry and a negative `B` entry.

    Species pair is no longer part of the IVP contract.
    There is no diffusion state, no `J_i`, and no `Pe_D` branch here.
*/
/// Simple structure for elementary chemical reactions with Arrhenius kinetics
///
/// Rate = A * T^n * exp(-E/(R*T)) * ∏[Ci]^νi
#[derive(Debug, Clone)]
pub struct FastElemReact {
    /// Chemical equation (e.g., "A + B => C + D")
    pub eq: String,
    /// Pre-exponential factor (units depend on reaction order)
    pub A: f64,
    /// Temperature exponent (dimensionless)
    pub n: f64,
    /// Activation energy (J/mol)
    pub E: f64,
    /// Heat of reaction (J/kg)
    pub Q: f64,
}

/// Main reactor modeling structure for condensed-phase IVP calculations.`r`n///`r`n/// The task keeps physics, solver backend configuration, and canonical state`r`n/// assembly in one place.
#[derive(Debug, Clone)]
pub struct SimpleReactorTask {
    /// Optional problem identifier
    pub problem_name: Option<String>,
    /// Optional problem description
    pub problem_description: Option<String>,
    /// Kinetic data (reactions, substances, rate constants)
    pub kindata: KinData,
    /// Heat effects for each reaction (J/kg)
    pub thermal_effects: Vec<f64>,
    /// Typed initial conditions for the condensed-phase IVP.
    pub initial_conditions: HashMap<String, f64>,
    /// Heat capacity (J/kg·K)
    pub Cp: f64,
    /// Thermal conductivity (W/m·K)
    pub Lambda: f64,
    pub ro: f64,
    /// Mass flow rate (kg/s)
    pub m: f64,
    /// Scaling parameters for dimensionless transformation
    pub scaling: ScalingConfig,
    /// Characteristic length (m)
    pub L: f64,
    /// Temperature scaling expression: T = dT*(Teta + 1)
    pub T_scaling: Expr,

    /// Thermal Peclet number: Pe_q = L*m*Cp/λ
    pub Pe_q: f64,

    /// Reaction rate expressions for each reaction
    pub map_eq_rate: HashMap<String, Expr>,
    /// System of differential equations (substance -> (variable, equation))
    pub map_of_equations: HashMap<String, (String, Expr)>,
    /// heat release function
    pub heat_release: Expr,
    /// Canonical state order used by the condensed-phase IVP helpers.
    pub canonical_state_order: IvpCanonicalStateOrder,
    /// Typed solver facade configuration for the condensed-phase IVP path.
    pub solver_backend_config: ReactorIvpSolverConfig,
    /// Last successful condensed-phase LSODE2 solve snapshot.
    pub last_solve_snapshot: Option<IvpSolveSnapshot>,
    /// Last computed condensed-phase conservation report.
    pub last_conservation_report: Option<IvpConservationReport>,
    /// IVP solver instance
    pub solver: IVPSolver,
}
/// Boundary Value Problem solver wrapper
///
/// Supports multiple solver backends:

#[derive(Debug, Clone)]
pub struct IVPSolver {
    /// Independent variable name (typically "x" or "z")
    pub arg_name: String,
    /// Domain range (start, end) - typically (0.0, 1.0) for dimensionless
    pub x_range: (f64, f64),
    /// Names of unknown variables ["Teta", "q", "C0", "J0", "C1", "J1", ...]
    pub unknowns: Vec<String>,
    /// System of differential equations dy/dx = f(x,y)
    pub eq_system: Vec<Expr>,
    /// Solution matrix (variables × mesh_points)
    pub solution: Option<DMatrix<f64>>,
    /// Spatial mesh points
    pub x_mesh: Option<DVector<f64>>,
}

impl Default for IVPSolver {
    fn default() -> Self {
        Self {
            arg_name: REACTOR_IVP_INDEPENDENT_VARIABLE.to_string(),
            x_range: (0.0, 1.0),
            unknowns: Vec::new(),
            eq_system: Vec::new(),
            solution: None,
            x_mesh: None,
        }
    }
}
impl IVPSolver {
    /// Get reference to solution matrix
    ///
    /// Returns None if solve hasn't been called yet
    pub fn get_solution(&self) -> Option<&DMatrix<f64>> {
        self.solution.as_ref()
    }

    /// Build a compact preview of the current solution without printing.
    ///
    /// The result is meant for tables, logs, or GUI inspection panels.
    pub fn solution_preview_rows(&self) -> Result<Vec<IvpSolutionPreviewRow>, IvpError> {
        let solution = self.solution.as_ref().ok_or_else(|| {
            IvpError::MissingData("No solution matrix is available for preview".to_string())
        })?;
        if solution.nrows() == 0 || solution.ncols() == 0 {
            return Err(IvpError::MissingData(
                "Solution matrix is empty and cannot be previewed".to_string(),
            ));
        }
        if self.unknowns.len() != solution.ncols() {
            return Err(IvpError::CalculationError(format!(
                "unknowns length {} does not match solution columns {}",
                self.unknowns.len(),
                solution.ncols()
            )));
        }

        let mut rows = Vec::with_capacity(self.unknowns.len());
        for (index, variable) in self.unknowns.iter().enumerate() {
            let column = solution.column(index);
            let last = column.len().saturating_sub(1);
            let middle = column.len() / 2;
            rows.push(IvpSolutionPreviewRow {
                index,
                variable: match variable.as_str() {
                    "Teta" => "T".to_string(),
                    other => other.to_string(),
                },
                first: column[0],
                middle: column[middle],
                last: column[last],
            });
        }
        Ok(rows)
    }

    /// Debug print solution summary
    ///
    /// Prints solution matrix dimensions and sample values for each variable
    pub fn debug_solution(&self) {
        if let Some(solution) = &self.solution {
            println!("\n=== SOLUTION DEBUG ===");
            println!(
                "Solution matrix shape: {} x {}",
                solution.nrows(),
                solution.ncols()
            );
            println!("Unknowns: {:?}", self.unknowns);
            match self.solution_preview_rows() {
                Ok(rows) => {
                    // Render preview rows without mutating the stored solution.
                    println!("{}", Table::new(rows));
                }
                Err(err) => {
                    println!("Preview unavailable: {}", err);
                }
            }
            println!("=== END DEBUG ===\n");
        }
    }
}

impl SimpleReactorTask {
    /// Create new reactor task with default values
    pub fn new() -> Self {
        Self {
            problem_name: None,
            problem_description: None,
            kindata: KinData::new(),
            thermal_effects: Vec::new(),
            initial_conditions: HashMap::new(),
            Cp: 0.0,
            Lambda: 0.0,
            ro: 0.0,

            m: 0.0,
            scaling: ScalingConfig::default(),
            L: 1.0,
            T_scaling: Expr::Const(0.0),

            Pe_q: 0.0,
            map_eq_rate: HashMap::new(),
            map_of_equations: HashMap::new(),
            heat_release: Expr::Const(0.0),
            canonical_state_order: IvpCanonicalStateOrder::default(),
            solver_backend_config: ReactorIvpSolverConfig::default(),
            last_solve_snapshot: None,
            last_conservation_report: None,
            solver: IVPSolver::default(),
        }
    }

    /// Returns the solved matrix if the reactor has already been solved.
    ///
    /// Keeping the availability check in one place prevents the conservation
    /// and preview helpers from repeating the same `Option` handling.
    fn solution_ref(&self) -> Result<&DMatrix<f64>, IvpError> {
        self.solver.solution.as_ref().ok_or_else(|| {
            IvpError::MissingData(
                "Solver solution is not available; run the condensed IVP solve first".to_string(),
            )
        })
    }

    /// Returns the x-mesh if the solver has already produced one.
    fn x_mesh_ref(&self) -> Result<&DVector<f64>, IvpError> {
        self.solver.x_mesh.as_ref().ok_or_else(|| {
            IvpError::MissingData(
                "Solver x mesh is not available; run the condensed IVP solve first".to_string(),
            )
        })
    }

    /// Returns the molecular weights used by the conservation conversion helpers.
    fn molar_masses_ref(&self) -> Result<&[f64], IvpError> {
        self.kindata
            .stecheodata
            .vec_of_molmasses
            .as_deref()
            .ok_or_else(|| {
                IvpError::MissingData("Molar masses are not available in stecheodata".to_string())
            })
    }

    /// Return the index of a canonical solver variable.
    fn unknown_index(&self, name: &str) -> Result<usize, IvpError> {
        self.solver
            .unknowns
            .iter()
            .position(|candidate| candidate == name)
            .ok_or_else(|| {
                IvpError::MissingData(format!(
                    "Variable `{}` is not available in the solved state",
                    name
                ))
            })
    }

    /// Extract the concentration columns from the full solver matrix.
    ///
    /// The condensed reactor keeps one state per species and stores them with
    /// canonical `C*` names.
    fn get_only_concentrations(&self) -> Result<DMatrix<f64>, IvpError> {
        let solution = self.solution_ref()?;
        let mut concentration_indices = Vec::new();
        for (index, unknown) in self.solver.unknowns.iter().enumerate() {
            if unknown.starts_with('C') {
                concentration_indices.push(index);
            }
        }
        if concentration_indices.is_empty() {
            return Err(IvpError::MissingData(
                "No concentration variables were found in the solved state".to_string(),
            ));
        }

        let mut only_concentrations = DMatrix::zeros(solution.nrows(), concentration_indices.len());
        for (column_index, &solution_index) in concentration_indices.iter().enumerate() {
            only_concentrations.set_column(column_index, &solution.column(solution_index));
        }
        Ok(only_concentrations)
    }

    /// Convert mass-fraction-like state columns to molar fractions.
    ///
    /// The helper mirrors the BVP postprocessing logic so the same conservation
    /// check can be applied to the condensed IVP path.
    pub fn from_mass_to_molar_fractions(
        &self,
        matrix_of_mass_fractions: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, IvpError> {
        let mi = self.molar_masses_ref()?;
        if mi.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(IvpError::InvalidNumericValue(
                "Molar masses must be finite and positive for fraction conversion".to_string(),
            ));
        }

        let mi: Vec<f64> = mi.iter().map(|m| *m / 1000.0).collect();
        let mi = DVector::from_vec(mi).transpose();
        let mut matrix_of_molar_fractions = DMatrix::zeros(
            matrix_of_mass_fractions.nrows(),
            matrix_of_mass_fractions.ncols(),
        );

        for (row_index, row) in matrix_of_mass_fractions.row_iter().enumerate() {
            let minv = mi.map(|x| 1.0 / x);
            let denom = row.dot(&minv);
            let row_new = row
                .iter()
                .enumerate()
                .map(|(species_index, value)| (value / mi[species_index]) / denom)
                .collect::<Vec<f64>>();
            matrix_of_molar_fractions.set_row(row_index, &DVector::from_vec(row_new).transpose());
        }

        Ok(matrix_of_molar_fractions)
    }

    /// Convert mass fractions to molar concentrations.
    pub fn from_mass_fractions_to_molar_concentration(
        &self,
        matrix_of_mass_fractions: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, IvpError> {
        let mi = self.molar_masses_ref()?;
        if mi.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(IvpError::InvalidNumericValue(
                "Molar masses must be finite and positive for concentration conversion".to_string(),
            ));
        }

        let mi: Vec<f64> = mi.iter().map(|m| *m / 1000.0).collect();
        let mi = DVector::from_vec(mi).transpose();
        let mut matrix_of_molar_concentrations = DMatrix::zeros(
            matrix_of_mass_fractions.nrows(),
            matrix_of_mass_fractions.ncols(),
        );

        for (row_index, row) in matrix_of_mass_fractions.row_iter().enumerate() {
            let row_new = row
                .iter()
                .enumerate()
                .map(|(species_index, value)| value / mi[species_index])
                .collect::<Vec<f64>>();
            matrix_of_molar_concentrations
                .set_row(row_index, &DVector::from_vec(row_new).transpose());
        }

        Ok(matrix_of_molar_concentrations)
    }

    /// Backward-compatible wrapper for the legacy misspelled API name.
    pub fn from_mass_fractions_to_molar_conentration(
        &self,
        matrix_of_mass_fractions: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, IvpError> {
        self.from_mass_fractions_to_molar_concentration(matrix_of_mass_fractions)
    }

    /// Build the condensed-phase conservation report from the current solve.
    pub fn conservation_report(&self) -> Result<IvpConservationReport, IvpError> {
        let matrix_of_elements = self
            .kindata
            .stecheodata
            .matrix_of_elements
            .as_ref()
            .ok_or_else(|| {
                IvpError::MissingData(
                    "Element matrix is not available for conservation checks".to_string(),
                )
            })?
            .transpose();

        let concentrations = self.get_only_concentrations()?;
        if concentrations.ncols() != self.kindata.substances.len() {
            return Err(IvpError::CalculationError(format!(
                "Concentration matrix has {} columns but {} substances are tracked",
                concentrations.ncols(),
                self.kindata.substances.len()
            )));
        }

        let mut sum_of_mass_fractions = Vec::new();
        for (row_index, row) in concentrations.row_iter().enumerate() {
            let sum: f64 = row.iter().sum();
            if (sum - 1.0).abs() > 0.01 {
                warn!(
                    "ATTENTION! Sum of concentrations in row {} is {:.3} not 1.0",
                    row_index, sum
                );
                sum_of_mass_fractions.push((row_index, sum));
            }
        }

        let molar_concentrations = self.from_mass_to_molar_fractions(concentrations)?;
        let initial_concentrations: DVector<f64> = molar_concentrations.row(0).transpose().into();
        let initial_vector_of_elements = &matrix_of_elements * initial_concentrations;

        let mut atomic_mass_balance_error = Vec::new();
        for row_index in 0..molar_concentrations.nrows() {
            let concentrations_at_step: DVector<f64> =
                molar_concentrations.row(row_index).transpose().into();
            let vector_of_elements = &matrix_of_elements * concentrations_at_step;
            let error = (&initial_vector_of_elements - &vector_of_elements).norm();
            if error > 0.01 {
                warn!(
                    "ATTENTION! Mass balance error in step {} is {:.3}",
                    row_index, error
                );
                atomic_mass_balance_error.push((row_index, error));
            }
        }

        let solution = self.solution_ref()?;
        let q_index = self.unknown_index("q")?;
        let t_index = self.unknown_index("Teta")?;
        let q_profile: Vec<f64> = solution.column(q_index).iter().copied().collect();
        let t_profile: Vec<f64> = solution.column(t_index).iter().copied().collect();
        let x_mesh: Vec<f64> = self.x_mesh_ref()?.iter().map(|x| x * self.L).collect();
        let heat_release_values = self.heat_release_profile_values()?;

        let (total_heat_release, _) =
            estimate_error_simpsons_richardson(&x_mesh, &heat_release_values)
                .map_err(|err| IvpError::CalculationError(err.to_string()))?;
        let (integrated_q, _) = estimate_error_simpsons_richardson(&x_mesh, &q_profile)
            .map_err(|err| IvpError::CalculationError(err.to_string()))?;

        let q_f = q_profile
            .last()
            .copied()
            .ok_or_else(|| IvpError::MissingData("Heat flux column is empty".to_string()))?
            * self.scaling.T_scale
            / self.L;
        let q_0 = q_profile
            .first()
            .copied()
            .ok_or_else(|| IvpError::MissingData("Heat flux column is empty".to_string()))?
            * self.scaling.T_scale
            / self.L;
        let dq = q_f - q_0;

        let t_f = t_profile
            .last()
            .copied()
            .ok_or_else(|| IvpError::MissingData("Temperature column is empty".to_string()))?
            * self.scaling.T_scale
            + self.scaling.dT;
        let t_0 = t_profile
            .first()
            .copied()
            .ok_or_else(|| IvpError::MissingData("Temperature column is empty".to_string()))?
            * self.scaling.T_scale
            + self.scaling.dT;
        let sensible_heat = self.m * self.Cp * (t_f - t_0);
        let absolute_energy_error = -dq - total_heat_release + sensible_heat;
        let energy_reference = total_heat_release.abs().max(1e-12);
        let energy_balance_error_rel = 100.0 * absolute_energy_error / energy_reference;

        if !absolute_energy_error.is_finite() || !energy_balance_error_rel.is_finite() {
            return Err(IvpError::CalculationError(
                "Energy conservation produced a non-finite result".to_string(),
            ));
        }

        let q_integral_reference = integrated_q.abs().max(1e-12);
        let q_integral_rel_error = 100.0 * (integrated_q.abs() / q_integral_reference);
        info!(
            "Condensed IVP conservation summary: integrated heat release {:.6}, integrated q {:.6}, relative q integration {:.3} %",
            total_heat_release, integrated_q, q_integral_rel_error
        );

        Ok(IvpConservationReport {
            energy_balance_error_abs: absolute_energy_error,
            energy_balance_error_rel,
            sum_of_mass_fractions,
            atomic_mass_balance_error,
        })
    }

    /// Return the latest computed conservation report, if one is available.
    pub fn latest_conservation_report(&self) -> Result<&IvpConservationReport, IvpError> {
        self.last_conservation_report.as_ref().ok_or_else(|| {
            IvpError::MissingData(
                "No conservation report is available; run the condensed IVP solve first"
                    .to_string(),
            )
        })
    }

    /// Refresh and store the current conservation report.
    fn refresh_conservation_report(&mut self) -> Result<(), IvpError> {
        let report = self.conservation_report()?;
        self.last_conservation_report = Some(report);
        Ok(())
    }

    /// Print the latest conservation summary as a compact table.
    fn pretty_print_conservation_report(&self, report: &IvpConservationReport) {
        let mut table = PrettyTable::new();
        table.add_row(row!["Parameter", "Value"]);
        for row_data in report.rows() {
            table.add_row(row![row_data.field, row_data.value]);
        }
        info!("{}", table);
    }

    /// Return the heat-release profile evaluated on the solved state.
    fn heat_release_profile_values(&self) -> Result<Vec<f64>, IvpError> {
        let solution = self.solution_ref()?;
        let unknowns = self
            .solver
            .unknowns
            .iter()
            .map(|unknown| unknown.as_str())
            .collect::<Vec<_>>();
        let mut values = Vec::with_capacity(solution.nrows());
        for solution_for_timestep in solution.row_iter() {
            let solution_for_timestep = solution_for_timestep.iter().copied().collect::<Vec<f64>>();
            let value = self
                .heat_release
                .eval_expression(unknowns.as_slice(), &solution_for_timestep);
            values.push(value);
        }
        Ok(values)
    }

    /// Create tolerance map from simplified config for this reactor's substances
    pub fn create_tolerance_map_for_system(
        &self,
        tolerance_config: HashMap<String, f64>,
    ) -> HashMap<String, f64> {
        create_tolerance_map(tolerance_config, &self.kindata.substances)
    }

    /// Create bounds map from simplified config for this reactor's substances
    pub fn create_bounds_map_for_system(
        &self,
        bounds_config: HashMap<String, (f64, f64)>,
    ) -> HashMap<String, (f64, f64)> {
        create_bounds_map(bounds_config, &self.kindata.substances)
    }

    /// Set scaling parameters using ScalingConfig
    pub fn set_scaling(&mut self, scaling: ScalingConfig) -> Result<(), ReactorError> {
        scaling.validate()?;
        self.scaling = scaling;
        Ok(())
    }

    /// Set scaling parameters from individual values
    pub fn set_scaling_values(
        &mut self,
        dT: f64,
        L: f64,
        T_scale: f64,
    ) -> Result<(), ReactorError> {
        let scaling = ScalingConfig::new(dT, L, T_scale);
        self.set_scaling(scaling)
    }
    /////////////////////////////////SETTERS////////////////////////////////////////////////////////////////////////////////
    /// Set problem name for identification
    pub fn set_problem_name(&mut self, name: &str) {
        self.problem_name = Some(name.to_string());
    }

    /// Set problem description
    pub fn set_problem_description(&mut self, description: &str) {
        self.problem_description = Some(description.to_string());
    }
    /// Set the typed initial conditions for the condensed-phase IVP.
    ///
    /// The map should contain `T`, `q`, and one entry per substance.
    pub fn set_initial_conditions(&mut self, conditions: HashMap<String, f64>) {
        self.initial_conditions = conditions;
    }

    /// Builder-style setter for condensed-phase initial conditions.
    pub fn with_initial_conditions(mut self, conditions: HashMap<String, f64>) -> Self {
        self.set_initial_conditions(conditions);
        self
    }

    /// Set the constant condensed-phase density.
    pub fn set_density(&mut self, ro: f64) -> Result<(), IvpError> {
        if !ro.is_finite() || ro <= 0.0 {
            return Err(IvpError::InvalidNumericValue(format!(
                "Density must be finite and positive, got {}",
                ro
            )));
        }
        self.ro = ro;
        Ok(())
    }

    /// Builder-style setter for the constant condensed-phase density.
    pub fn with_density(mut self, ro: f64) -> Result<Self, IvpError> {
        self.set_density(ro)?;
        Ok(self)
    }

    /// Set the typed solver-backend contract for the condensed-phase IVP.
    pub fn set_solver_backend_config(&mut self, config: ReactorIvpSolverConfig) {
        self.solver_backend_config = config;
    }

    /// Builder-style setter for the typed solver-backend contract.
    pub fn with_solver_backend_config(mut self, config: ReactorIvpSolverConfig) -> Self {
        self.set_solver_backend_config(config);
        self
    }

    /// Collect the typed condensed-phase input bundle from the current task state.
    ///
    /// This is the canonical reactor-facing configuration object used by the
    /// condensed IVP path. Legacy scalar fields remain on the task for
    /// compatibility, but new code should prefer the typed bundle.
    pub fn condensed_inputs(&self) -> Result<ReactorIvpTaskInputs, IvpError> {
        let mut species = HashMap::with_capacity(self.kindata.substances.len());
        for substance in &self.kindata.substances {
            let value = *self.initial_conditions.get(substance).ok_or_else(|| {
                IvpError::MissingData(format!(
                    "Missing initial condition for substance `{}`",
                    substance
                ))
            })?;
            species.insert(substance.clone(), value);
        }

        let initial_state = ReactorIvpInitialStateConfig {
            temperature: *self.initial_conditions.get("T").ok_or_else(|| {
                IvpError::MissingData("Missing `T` in initial conditions".to_string())
            })?,
            heat_flux: *self.initial_conditions.get("q").ok_or_else(|| {
                IvpError::MissingData("Missing `q` in initial conditions".to_string())
            })?,
            species,
        };

        Ok(ReactorIvpTaskInputs {
            physical: ReactorIvpPhysicalConfig {
                ro: self.ro,
                Cp: self.Cp,
                Lambda: self.Lambda,
                m: self.m,
                L: self.scaling.L,
            },
            scaling: self.scaling.clone(),
            initial_state,
            thermal_effects: self.thermal_effects.clone(),
        })
    }

    /// Return `true` when the task has enough data for the condensed-phase path.
    pub fn prefers_condensed_setup(&self) -> bool {
        !self.initial_conditions.is_empty()
    }

    /// Return the canonical state vector names for the condensed-phase IVP.
    pub fn canonical_state_names(&self) -> Vec<String> {
        let mut names = vec!["Teta".to_string(), "q".to_string()];
        names.extend(
            self.kindata
                .substances
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("C{}", idx)),
        );
        names
    }

    /// Build the canonical initial-state vector.
    ///
    /// The stored map keeps dimensional input values. The thermal entry is
    /// converted to the canonical `Teta` coordinate here so that
    /// equation assembly and initial state stay in one canonical order.
    pub fn build_initial_state_vector(&self) -> Result<DVector<f64>, IvpError> {
        let typed_inputs = self.condensed_inputs()?;
        typed_inputs.validate(&self.kindata.substances)?;

        let t0 = typed_inputs.initial_state.temperature;
        let q0 = typed_inputs.initial_state.heat_flux;
        let mut values = Vec::with_capacity(2 + self.kindata.substances.len());
        values.push((t0 - self.scaling.dT) / self.scaling.T_scale);
        values.push(q0);

        for substance in &self.kindata.substances {
            let value = *typed_inputs
                .initial_state
                .species
                .get(substance)
                .ok_or_else(|| {
                    IvpError::MissingData(format!(
                        "Missing initial condition for substance `{}`",
                        substance
                    ))
                })?;
            values.push(value);
        }

        Ok(DVector::from_vec(values))
    }

    /// Validate the condensed-phase physical contract before assembling equations.
    pub fn validate_condensed_phase_contract(&self) -> Result<(), IvpError> {
        let typed_inputs = self.condensed_inputs()?;
        typed_inputs.validate(&self.kindata.substances)?;

        if self.kindata.substances.is_empty() {
            return Err(IvpError::MissingData(
                "No substances found in kinetic data".to_string(),
            ));
        }
        if self.kindata.vec_of_equations.is_empty() {
            return Err(IvpError::MissingData(
                "No reactions found in kinetic data".to_string(),
            ));
        }
        if self.thermal_effects.len() != self.kindata.vec_of_equations.len() {
            return Err(IvpError::InvalidConfiguration(
                "Thermal effects length must match number of reactions".to_string(),
            ));
        }
        Ok(())
    }

    /// Condensed-phase setup workflow.
    pub fn setup_condensed_ivp(&mut self) -> Result<(), IvpError> {
        self.validate_condensed_phase_contract()?;
        self.scaling_processing().map_err(IvpError::from)?;
        self.kinetic_processing().map_err(IvpError::from)?;
        self.peclet_numbers().map_err(IvpError::from)?;
        self.solver.arg_name = REACTOR_IVP_INDEPENDENT_VARIABLE.to_string();
        self.solver.x_range = (
            self.solver_backend_config.x0,
            self.solver_backend_config.x_bound,
        );
        self.create_condensed_ivp_equations()
            .map_err(IvpError::from)?;
        // let _handoff = self.build_lsode2_problem_config()?;
        self.check_condensed_before_solution()?;
        Ok(())
    }

    /// Validate the canonical condensed-phase snapshot before solve handoff.
    pub fn check_condensed_before_solution(&self) -> Result<(), IvpError> {
        let expected_len = self.kindata.substances.len() + 2;
        if self.solver.eq_system.len() != expected_len {
            return Err(IvpError::InvalidConfiguration(format!(
                "eq_system length {} != expected {}",
                self.solver.eq_system.len(),
                expected_len
            )));
        }
        if self.solver.unknowns.len() != expected_len {
            return Err(IvpError::InvalidConfiguration(format!(
                "unknowns length {} != expected {}",
                self.solver.unknowns.len(),
                expected_len
            )));
        }
        if self.map_of_equations.len() != expected_len {
            return Err(IvpError::InvalidConfiguration(format!(
                "map_of_equations length {} != expected {}",
                self.map_of_equations.len(),
                expected_len
            )));
        }
        Ok(())
    }

    /// Build the typed LSODE2 problem configuration for the condensed-phase IVP.
    ///
    /// This is the handoff boundary between reactor physics and the solver
    /// backend. The returned config should be enough to instantiate
    /// `RustedSciThe::numerical::LSODE2::Lsode2Solver` without re-reading any
    /// additional reactor state.
    pub fn build_lsode2_problem_config(&self) -> Result<Lsode2ProblemConfig, IvpError> {
        self.check_condensed_before_solution()?;
        let initial_state = self.build_initial_state_vector()?;
        self.solver_backend_config
            .to_rusted_problem_config(
                self.solver.eq_system.clone(),
                self.solver.unknowns.clone(),
                REACTOR_IVP_INDEPENDENT_VARIABLE.to_string(),
                initial_state,
            )
            .map_err(IvpError::from)
    }

    /// Build an LSODE2 solver from the currently assembled condensed-phase state.
    pub fn build_condensed_lsode2_solver(&self) -> Result<Lsode2Solver, IvpError> {
        let config = self.build_lsode2_problem_config()?;
        Lsode2Solver::new(config).map_err(|err| {
            IvpError::CalculationError(format!("failed to build LSODE2 solver: {err}"))
        })
    }

    /// Assemble the condensed-phase symbolic RHS on the canonical state order.
    pub fn create_condensed_ivp_equations(&mut self) -> Result<(), IvpError> {
        let substances = &self.kindata.substances;
        let n = substances.len();
        if n == 0 {
            return Err(IvpError::MissingData(
                "No substances found in kinetic data".to_string(),
            ));
        }

        let equations = self.kindata.vec_of_equations.clone();
        let k = equations.len();
        if k == 0 {
            return Err(IvpError::MissingData(
                "No reactions found in kinetic data".to_string(),
            ));
        }

        let k_sym_vec = self.kindata.K_sym_vec.as_ref().ok_or_else(|| {
            IvpError::MissingData("Symbolic rate constants not calculated".to_string())
        })?;
        if k_sym_vec.len() != k {
            return Err(IvpError::InvalidConfiguration(
                "Mismatch between number of reactions and rate constants".to_string(),
            ));
        }

        let stoich_matrix = &self.kindata.stecheodata.stecheo_matrx;
        let reactant_powers = &self.kindata.stecheodata.stecheo_reags;
        if stoich_matrix.len() != k || reactant_powers.len() != k {
            return Err(IvpError::InvalidConfiguration(
                "Stoichiometric matrix size mismatch".to_string(),
            ));
        }

        let molar_masses = self
            .kindata
            .stecheodata
            .vec_of_molmasses
            .as_ref()
            .ok_or_else(|| IvpError::MissingData("Molar masses not calculated".to_string()))?;
        if molar_masses.len() != n {
            return Err(IvpError::InvalidConfiguration(format!(
                "Molar mass count {} must match substance count {}",
                molar_masses.len(),
                n
            )));
        }

        let heat_effects = if self.thermal_effects.len() == k {
            self.thermal_effects.clone()
        } else {
            return Err(IvpError::InvalidConfiguration(
                "Thermal effects length must match number of reactions".to_string(),
            ));
        };

        let (ci_expr, ci_names) = Expr::IndexedVars(n, "C");
        let ci_names: Vec<String> = ci_names.iter().map(|name| name.replace('_', "")).collect();

        let mut map_eq_rate: HashMap<String, Expr> = HashMap::with_capacity(k);
        let mut rates: Vec<Expr> = Vec::with_capacity(k);
        let ro_m = Expr::Const(self.ro);

        for reaction_index in 0..k {
            let mut rate_expr = k_sym_vec[reaction_index].clone();
            for species_index in 0..n {
                let power = reactant_powers
                    .get(reaction_index)
                    .and_then(|row| row.get(species_index))
                    .ok_or_else(|| {
                        IvpError::IndexOutOfBounds(format!(
                            "reactant_powers[{}][{}] out of bounds",
                            reaction_index, species_index
                        ))
                    })?;
                let ci = ci_expr.get(species_index).ok_or_else(|| {
                    IvpError::IndexOutOfBounds(format!("C[{}] out of bounds", species_index))
                })?;
                let molar_mass = Expr::Const(molar_masses[species_index] / 1000.0);
                rate_expr =
                    rate_expr * (ro_m.clone() * ci.clone() / molar_mass).pow(Expr::Const(*power));
            }
            let rate_expr = rate_expr.simplify_();
            map_eq_rate.insert(equations[reaction_index].clone(), rate_expr.clone());
            rates.push(rate_expr);
        }

        let mut unknowns = Vec::with_capacity(n + 2);
        let mut rhs = Vec::with_capacity(n + 2);
        let mut map_of_equations: HashMap<String, (String, Expr)> = HashMap::with_capacity(n + 2);

        unknowns.push("Teta".to_string());
        let q = Expr::Var("q".to_string());
        let rhs_teta = q.clone() / Expr::Const(self.Lambda);
        rhs.push(rhs_teta.clone());
        map_of_equations.insert("Teta".to_string(), ("Teta".to_string(), rhs_teta));

        unknowns.push("q".to_string());
        let mut heat_release = Expr::Const(0.0);
        for (reaction_index, heat_effect) in heat_effects.iter().enumerate() {
            heat_release = (heat_release
                + Expr::Const(*heat_effect) * rates[reaction_index].clone())
            .simplify_();
        }
        self.heat_release = heat_release.clone();
        let rhs_q = q * Expr::Const(self.Pe_q)
            - heat_release * Expr::Const(self.L.powi(2) / self.scaling.T_scale);
        rhs.push(rhs_q.clone());
        map_of_equations.insert("q".to_string(), ("q".to_string(), rhs_q));

        for species_index in 0..n {
            let substance = &substances[species_index];
            let c_name = ci_names
                .get(species_index)
                .ok_or_else(|| {
                    IvpError::IndexOutOfBounds(format!("C[{}] out of bounds", species_index))
                })?
                .clone();
            unknowns.push(c_name.clone());

            let mut mass_source = Expr::Const(0.0);
            for reaction_index in 0..k {
                let stoich_coeff = stoich_matrix
                    .get(reaction_index)
                    .and_then(|row| row.get(species_index))
                    .ok_or_else(|| {
                        IvpError::IndexOutOfBounds(format!(
                            "stoich_matrix[{}][{}] out of bounds",
                            reaction_index, species_index
                        ))
                    })?;
                let molar_mass = Expr::Const(molar_masses[species_index] / 1000.0);
                mass_source = mass_source
                    + rates[reaction_index].clone() * Expr::Const(*stoich_coeff) * molar_mass;
            }

            let rhs_species = (-mass_source * Expr::Const(self.L.powi(2))).simplify_();
            rhs.push(rhs_species.clone());
            map_of_equations.insert(substance.clone(), (c_name, rhs_species));
        }

        self.map_eq_rate = map_eq_rate;
        self.solver.eq_system = rhs;
        self.solver.unknowns = unknowns;
        self.map_of_equations = map_of_equations;

        Ok(())
    }

    /// Validate a condensed-phase LSODE2 result snapshot before publishing it.
    pub(crate) fn validate_condensed_solve_snapshot(
        &self,
        axis: &DVector<f64>,
        values: &DMatrix<f64>,
        summary: &Lsode2SolveSummary,
    ) -> Result<(), IvpError> {
        if axis.is_empty() {
            return Err(IvpError::CalculationError(
                "LSODE2 returned an empty axis".to_string(),
            ));
        }
        if values.nrows() != axis.len() {
            return Err(IvpError::CalculationError(format!(
                "solution row count {} must match axis length {}",
                values.nrows(),
                axis.len()
            )));
        }
        if values.ncols() != self.solver.unknowns.len() {
            return Err(IvpError::CalculationError(format!(
                "solution column count {} must match variable count {}",
                values.ncols(),
                self.solver.unknowns.len()
            )));
        }
        if summary.time_points != axis.len() {
            return Err(IvpError::CalculationError(format!(
                "summary time_points {} must match axis length {}",
                summary.time_points,
                axis.len()
            )));
        }
        if summary.variable_count != values.ncols() {
            return Err(IvpError::CalculationError(format!(
                "summary variable_count {} must match solution column count {}",
                summary.variable_count,
                values.ncols()
            )));
        }
        if summary.status.trim().is_empty() {
            return Err(IvpError::CalculationError(
                "LSODE2 returned an empty status".to_string(),
            ));
        }
        if summary.method.trim().is_empty() {
            return Err(IvpError::CalculationError(
                "LSODE2 returned an empty method label".to_string(),
            ));
        }
        if !Self::is_accepted_solve_status(&summary.status) {
            return Err(IvpError::CalculationError(format!(
                "LSODE2 returned an unexpected termination status: {}",
                summary.status
            )));
        }
        if summary.algorithm.controller_mode.trim().is_empty()
            || summary.algorithm.active_family.trim().is_empty()
            || summary.algorithm.mused_family.trim().is_empty()
            || summary.algorithm.mcur_family.trim().is_empty()
            || summary.algorithm.preferred_family.trim().is_empty()
            || summary.algorithm.switch_reason.trim().is_empty()
            || summary.algorithm.note.trim().is_empty()
        {
            return Err(IvpError::CalculationError(
                "LSODE2 returned incomplete algorithm diagnostics".to_string(),
            ));
        }
        if summary.jacobian_backend.trim().is_empty()
            || summary.linear_solver_backend.trim().is_empty()
            || summary.linear_solver_reason.trim().is_empty()
            || summary.resolved_source.trim().is_empty()
            || summary.resolved_structure.trim().is_empty()
        {
            return Err(IvpError::CalculationError(
                "LSODE2 returned incomplete backend diagnostics".to_string(),
            ));
        }
        if !summary.max_abs_solution.is_finite() || summary.max_abs_solution < 0.0 {
            return Err(IvpError::CalculationError(
                "LSODE2 returned an invalid solution magnitude".to_string(),
            ));
        }
        for (idx, value) in axis.iter().enumerate() {
            if !value.is_finite() {
                return Err(IvpError::CalculationError(format!(
                    "axis value at index {} is not finite: {}",
                    idx, value
                )));
            }
        }
        for (idx, value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(IvpError::CalculationError(format!(
                    "solution value at linear index {} is not finite: {}",
                    idx, value
                )));
            }
        }
        for pair in axis.iter().zip(axis.iter().skip(1)) {
            if pair.1 < pair.0 {
                return Err(IvpError::CalculationError(
                    "LSODE2 axis must be monotonic increasing".to_string(),
                ));
            }
        }
        if let Some(final_t) = summary.final_t {
            let last_axis = axis[axis.len() - 1];
            let tolerance = 1e-10_f64.max(1e-10 * last_axis.abs());
            if !final_t.is_finite() || (final_t - last_axis).abs() > tolerance {
                return Err(IvpError::CalculationError(format!(
                    "summary final_t {} does not match axis tail {}",
                    final_t, last_axis
                )));
            }
        }
        if let Some(final_y) = &summary.final_y {
            if final_y.len() != values.ncols() {
                return Err(IvpError::CalculationError(format!(
                    "summary final_y length {} must match solution column count {}",
                    final_y.len(),
                    values.ncols()
                )));
            }
            for (idx, value) in final_y.iter().enumerate() {
                if !value.is_finite() {
                    return Err(IvpError::CalculationError(format!(
                        "summary final_y value at index {} is not finite: {}",
                        idx, value
                    )));
                }
            }
        }
        Ok(())
    }

    /// Accept the statuses that correspond to successful or intentionally stopped solves.
    fn is_accepted_solve_status(status: &str) -> bool {
        status.starts_with("finished") || status == "stopped_by_condition"
    }

    /// Solve the condensed-phase IVP and publish a validated typed snapshot.
    ///
    /// The task keeps the last successful result in memory so repeated failed
    /// reruns do not destroy the previous usable output.
    pub fn solve_condensed_ivp(&mut self) -> Result<IvpSolveSnapshot, IvpError> {
        self.setup_condensed_ivp()?;
        let mut solver = self.build_condensed_lsode2_solver()?;
        let summary = solver.solve_with_summary().map_err(IvpError::from)?;
        let (axis, values) = solver.get_result();
        let algorithm = summary.algorithm.clone();
        let statistics = summary.statistics.clone();
        let native_statistics = summary.native_statistics.clone();
        let evaluation_telemetry = summary.evaluation_telemetry.clone();

        self.validate_condensed_solve_snapshot(&axis, &values, &summary)?;

        self.solver.x_mesh = Some(axis.clone());
        self.solver.solution = Some(values.clone());
        let snapshot = IvpSolveSnapshot {
            axis,
            values,
            variable_names: self.solver.unknowns.clone(),
            summary,
            backend_plan: self.solver_backend_config,
            algorithm,
            statistics,
            native_statistics,
            evaluation_telemetry,
        };
        self.last_solve_snapshot = Some(snapshot.clone());
        match self.refresh_conservation_report() {
            Ok(()) => {
                if let Some(report) = self.last_conservation_report.as_ref() {
                    self.pretty_print_conservation_report(report);
                }
            }
            Err(err) => {
                warn!(
                    "Condensed IVP conservation check could not be completed: {}",
                    err
                );
            }
        }
        if let Err(err) = self.postprocessing() {
            return Err(err);
        }
        Ok(snapshot)
    }

    /// Solve the condensed-phase IVP using the currently configured backend.
    ///
    /// This is the narrow public entry point intended for callers that already
    /// prepared the task state through the builder or setup methods.
    pub fn solve(&mut self) -> Result<IvpSolveSnapshot, IvpError> {
        self.solve_condensed_ivp()
    }

    /// Solve the condensed-phase IVP with a temporary backend override.
    ///
    /// The caller provides a typed LSODE2 backend plan while the task keeps the
    /// physical reactor state intact. On failure the previous backend plan is
    /// restored so a rejected override does not poison the task.
    pub fn solve_with_config(
        &mut self,
        config: ReactorIvpSolverConfig,
    ) -> Result<IvpSolveSnapshot, IvpError> {
        let previous = self.solver_backend_config;
        self.solver_backend_config = config;
        match self.solve() {
            Ok(snapshot) => Ok(snapshot),
            Err(err) => {
                self.solver_backend_config = previous;
                Err(err)
            }
        }
    }

    /// Return a diagnostics-only view of the most recent validated solve.
    pub fn latest_solve_report(&self) -> Result<IvpSolveReport, IvpError> {
        let snapshot = self.last_solve_snapshot.as_ref().ok_or_else(|| {
            IvpError::MissingData("No successful IVP solve snapshot is available".to_string())
        })?;
        Ok(snapshot.diagnostics_report())
    }

    /// Return the latest diagnostics as plain rows for table rendering.
    pub fn latest_solve_report_rows(&self) -> Result<Vec<IvpReportRow>, IvpError> {
        let snapshot = self.last_solve_snapshot.as_ref().ok_or_else(|| {
            IvpError::MissingData("No successful IVP solve snapshot is available".to_string())
        })?;
        Ok(snapshot.diagnostics_report_rows())
    }

    /// Return a borrowed read-only view of the most recent validated solve.
    pub fn latest_result_view(&self) -> Result<IvpResultView<'_>, IvpError> {
        let snapshot = self.last_solve_snapshot.as_ref().ok_or_else(|| {
            IvpError::MissingData("No successful IVP solve snapshot is available".to_string())
        })?;
        Ok(snapshot.result_view())
    }

    /// Build an owned render snapshot from the current solver state.
    ///
    /// The render snapshot is expected to be called after postprocessing so the
    /// user-facing view reflects physical units instead of canonical scaled
    /// variables.
    pub fn solution_render_data(&self) -> Result<IvpSolutionRenderData, IvpError> {
        let unknowns = self
            .solver
            .unknowns
            .iter()
            .map(|name| match name.as_str() {
                "Teta" => "T".to_string(),
                other => other.to_string(),
            })
            .collect::<Vec<_>>();
        Ok(IvpSolutionRenderData {
            x_mesh: self.x_mesh_ref()?.clone(),
            solution: self.solution_ref()?.clone(),
            unknowns,
            arg_name: self.solver.arg_name.clone(),
        })
    }

    /// Build a dimensional postprocessing snapshot without mutating the solver.
    ///
    /// The raw LSODE2 solve works in the canonical IVP coordinates, but the
    /// GUI and export paths should show dimensional values after the solve.
    pub fn postprocessing_report(&self) -> Result<IvpPostprocessingReport, IvpError> {
        let mut x_mesh = self.x_mesh_ref()?.clone();
        let mut solution = self.solution_ref()?.clone();

        if self.solver.unknowns.len() != solution.ncols() {
            return Err(IvpError::CalculationError(format!(
                "Solution has {} columns but {} unknowns are tracked",
                solution.ncols(),
                self.solver.unknowns.len()
            )));
        }

        let dT = self.scaling.dT;
        let T_scale = self.scaling.T_scale;
        let L = self.L;

        x_mesh.iter_mut().for_each(|x| *x *= L);
        for (index, mut column) in solution.column_iter_mut().enumerate() {
            let variable = self.solver.unknowns.get(index).ok_or_else(|| {
                IvpError::IndexOutOfBounds(format!(
                    "Solution column {} is missing a matching variable name",
                    index
                ))
            })?;
            if variable == "Teta" {
                column
                    .iter_mut()
                    .for_each(|value| *value = *value * T_scale + dT);
            } else if variable == "q" {
                column
                    .iter_mut()
                    .for_each(|value| *value = *value * T_scale / L);
            }
        }

        Ok(IvpPostprocessingReport { x_mesh, solution })
    }

    /// Convert the solved IVP state to dimensional output in place.
    ///
    /// The raw canonical snapshot is still preserved in `last_solve_snapshot`,
    /// but the live solver state is rewritten to the physical coordinate system
    /// so plots and previews follow the same contract as the BVP path.
    pub fn postprocessing(&mut self) -> Result<(), IvpError> {
        let report = self.postprocessing_report()?;
        self.solver.x_mesh = Some(report.x_mesh);
        self.solver.solution = Some(report.solution);
        Ok(())
    }

    /// Set transport properties
    pub fn set_transport_properties(&mut self, lambda: f64, cp: f64) {
        self.Lambda = lambda;
        self.Cp = cp;
    }

    /// Set thermal effects (heat of reaction) for each reaction
    pub fn set_thermal_effects(&mut self, thermal_effects: Vec<f64>) {
        self.thermal_effects = thermal_effects;
    }

    /// Set elementary reactions from FastElemReact structures
    ///
    /// Convenient method to quickly set up elementary reactions with Arrhenius kinetics
    pub fn fast_react_set(&mut self, vec_of_maps: Vec<FastElemReact>) -> Result<(), ReactorError> {
        let mut eq_vec: Vec<String> = Vec::new();
        let mut elementary_reaction_vec = Vec::new();
        let mut Q_vec = Vec::new();
        for (idx, map_of_reactiondata) in vec_of_maps.iter().enumerate() {
            // Check equation string
            let eq = if !map_of_reactiondata.eq.is_empty() {
                map_of_reactiondata.eq.clone()
            } else {
                return Err(ReactorError::MissingData(format!(
                    "No equation in input hashmap at index {}",
                    idx
                )));
            };

            // Check Arrhenius parameters
            let A = map_of_reactiondata.A;
            let n = map_of_reactiondata.n;
            let E = map_of_reactiondata.E;
            let Q = map_of_reactiondata.Q;
            // Check for NaN (in case of uninitialized f64)
            if A.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'A' in input hashmap at index {}",
                    idx
                )));
            }
            if n.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'n' in input hashmap at index {}",
                    idx
                )));
            }
            if E.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'E' in input hashmap at index {}",
                    idx
                )));
            }
            if Q.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'Q' in input hashmap at index {}",
                    idx
                )));
            }
            let arrenius = vec![A, n, E];
            let reactdata = ReactionData::new_elementary(eq.clone(), arrenius, None);
            eq_vec.push(eq);
            Q_vec.push(Q);
            elementary_reaction_vec.push(reactdata);
        }

        let mut kindata = KinData::new();
        kindata.vec_of_equations = eq_vec;
        kindata.vec_of_reaction_data = Some(elementary_reaction_vec);
        self.kindata = kindata;
        self.thermal_effects = Q_vec;
        Ok(())
    }
    ///////////////////////////////////////////VALIDATION////////////////////////////////////////////////
    ///////////////////////////////////////////KINETICS AND THERMAL PREPROCESSING////////////////////////////////////////////////

    ///
    pub fn kinetic_processing(&mut self) -> Result<(), ReactorError> {
        let kd = &mut self.kindata;
        // stoichiometry and element matrix
        kd.analyze_reactions()
            .map_err(|err| ReactorError::CalculationError(err.to_string()))?;
        // in elementary reactions there are only Arrhenius parameters - no concentration or pressure dependencies
        kd.calc_sym_constants(None, None, Some(self.T_scaling.clone()))
            .map_err(|err| ReactorError::CalculationError(err.to_string()))?;
        Ok(())
    }

    ///
    ///
    /// Process scaling parameters for dimensionless equations
    /// Teta = (T - dT)/T_scaling
    /// Creates temperature scaling: T = (Teta*T_scaling + dT)
    /// Sets characteristic length L for spatial scaling: z = x/L
    pub fn scaling_processing(&mut self) -> Result<(), ReactorError> {
        // Validate scaling parameters
        self.scaling.validate()?;

        // Create temperature scaling expression: T = dT*(Teta + 1)
        let dT = self.scaling.dT;
        let T_scale = self.scaling.T_scale;
        let Teta = Expr::Var("Teta".to_owned());
        self.T_scaling = Teta.clone() * Expr::Const(T_scale) + Expr::Const(dT);

        // Set characteristic length.
        self.L = self.scaling.L;
        Ok(())
    }
    pub fn peclet_numbers(&mut self) -> Result<(), ReactorError> {
        if self.Lambda <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Lambda must be positive".to_string(),
            ));
        }
        if self.m <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Mass flow rate must be positive".to_string(),
            ));
        }
        if self.L.is_nan() {
            return Err(ReactorError::InvalidConfiguration("missing L".to_string()));
        }
        if self.m.is_nan() {
            return Err(ReactorError::InvalidConfiguration("missing m".to_string()));
        }
        if self.Cp.is_nan() {
            return Err(ReactorError::InvalidConfiguration("missing Cp".to_string()));
        }
        if self.Lambda.is_nan() {
            return Err(ReactorError::InvalidConfiguration(
                "missing Lambda".to_string(),
            ));
        }

        let Pe_q = (self.L * self.m * self.Cp) / self.Lambda;

        self.Pe_q = Pe_q;

        Ok(())
    }
}

/// Integrate sampled values with the composite trapezoidal rule.
fn trapezoidal(y: &Vec<f64>, x: &Vec<f64>) -> Result<f64, String> {
    if x.len() != y.len() {
        return Err("x and y must have the same length".to_string());
    }
    if x.len() < 2 {
        return Err("At least 2 points required for integration".to_string());
    }

    let mut integral = 0.0;
    for (x_pair, y_pair) in x.windows(2).zip(y.windows(2)) {
        let h = x_pair[1] - x_pair[0];
        integral += h * (y_pair[0] + y_pair[1]) / 2.0;
    }
    Ok(integral)
}

/// Integrate sampled values with Simpson's rule when the spacing is uniform.
///
/// If the local spacing is not uniform, the helper falls back to a trapezoid
/// step for the current segment.
fn simpsons(y: &Vec<f64>, x: &Vec<f64>) -> Result<f64, String> {
    if x.len() != y.len() {
        return Err("x and y must have the same length".to_string());
    }
    if x.len() < 2 {
        return Err("At least 2 points required for integration".to_string());
    }
    if x.len() == 2 {
        return trapezoidal(y, x);
    }

    let mut integral = 0.0;
    let mut i = 0;
    while i + 2 < x.len() {
        let h1 = x[i + 1] - x[i];
        let h2 = x[i + 2] - x[i + 1];
        if (h1 - h2).abs() < 1e-10 {
            let h = h1 + h2;
            integral += h / 6.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
        } else {
            integral += 0.5 * (y[i] + y[i + 1]) * h1;
            if i + 2 == x.len() - 1 {
                integral += 0.5 * (y[i + 1] + y[i + 2]) * h2;
            }
        }
        i += 2;
    }

    if i + 1 < x.len() {
        integral += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
    }

    Ok(integral)
}

/// Estimate the trapezoidal integration error with Richardson extrapolation.
#[allow(dead_code)]
fn estimate_error_richardson(x: &Vec<f64>, y: &Vec<f64>) -> Result<(f64, f64), String> {
    let i_fine = trapezoidal(y, x)?;
    let x_coarse: Vec<f64> = x.iter().step_by(2).copied().collect();
    let y_coarse: Vec<f64> = y.iter().step_by(2).copied().collect();

    if x_coarse.len() < 2 {
        return Ok((i_fine, 0.0));
    }

    let i_coarse = trapezoidal(&y_coarse, &x_coarse)?;
    Ok((i_fine, (i_coarse - i_fine) / 3.0))
}

/// Estimate the Simpson integration error with Richardson extrapolation.
fn estimate_error_simpsons_richardson(x: &Vec<f64>, y: &Vec<f64>) -> Result<(f64, f64), String> {
    let i_fine = simpsons(y, x)?;
    let x_coarse: Vec<f64> = x.iter().step_by(2).copied().collect();
    let y_coarse: Vec<f64> = y.iter().step_by(2).copied().collect();

    if x_coarse.len() < 3 {
        return Ok((i_fine, 0.0));
    }

    let i_coarse = simpsons(&y_coarse, &x_coarse)?;
    Ok((i_fine, (i_coarse - i_fine) / 15.0))
}

/////////////////////////////////////////TESTS/////////////////////////////////////////////////////////
