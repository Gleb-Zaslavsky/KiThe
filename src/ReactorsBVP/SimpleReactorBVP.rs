//! # Simple Reactor BVP Module
//!
//! This module provides a comprehensive framework for modeling chemical reactors using boundary value problems (BVP).
//! It implements dimensionless reactor equations for mass and heat transfer with chemical reactions.
//!
//! ## Main Structures
//!
//! - **`SimpleReactorTask`**: Main reactor modeling structure that aggregates kinetics, thermodynamics, and transport properties
//! - **`BVPSolver`**: Wrapper for different BVP solvers (NRBVP damped Newton-Raphson, BVPsci adaptive)
//! - **`ToleranceConfig`**: Helper for setting solver tolerances across all variables (C, J, Teta, q)
//! - **`BoundsConfig`**: Helper for setting variable bounds across all variables
//! - **`FastElemReact`**: Simple structure for elementary reactions with Arrhenius parameters
//!
//! ## Key Features
//!
//! - **Automatic dimensionless scaling**: Converts dimensional equations to dimensionless form for numerical stability
//! - **Flexible tolerance/bounds setup**: Automatically expands simple configs to full variable maps (C0,C1,... J0,J1,...)
//! - **Multiple solver backends**: Supports both damped Newton-Raphson and adaptive mesh refinement solvers
//! - **Transport property calculation**: Automatic calculation of Peclet numbers and Lewis numbers
//! - **Post-processing**: Converts dimensionless results back to dimensional form
//!
//! ## Mathematical Model
//!
//! The module solves dimensionless reactor equations:
//! - **Mass balance**: dCi/dz = Ji/(D*ρ), dJi/dz = Pe_D*Ji - l²*Gi
//! - **Energy balance**: dTeta/dz = q/λ, dq/dz = Pe_q*q - l²*Q/dT
//!
//! Where z = x/L (dimensionless coordinate), Teta = (T-T₀)/dT (dimensionless temperature)
//! PAY ATTENTION TO THE DIMENSION OF INPUT PARAMETERS
use crate::Kinetics::User_reactions::KinData;
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::ReactorsBVP::reactor_BVP_utils::{
    BoundsConfig, ScalingConfig, ToleranceConfig, create_bounds_map, create_tolerance_map,
};
use crate::ReactorsBVP::solver_backend::ReactorBvpSolverConfig;
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{
    DampedSolverOptions, NRBVP, SolverParams,
};
use RustedSciThe::numerical::BVP_sci::BVP_sci_symb::BVPwrap as BVPsci;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use log::{info, warn};

use nalgebra::{DMatrix, DVector};

use std::collections::HashMap;
use std::fmt;

/// Universal gas constant in J/(mol·K)
pub const R_G: f64 = 8.314;

#[derive(Debug, Clone)]
pub struct SolutionQuality {
    pub energy_balane_error_abs: f64,
    pub energy_balane_error_rel: f64,

    /// steps where sum of mass fractions is larger than the threshold
    pub sum_of_mass_fractions: Vec<(usize, f64)>,
    pub atomic_mass_balance_error: Vec<(usize, f64)>,
}

impl Default for SolutionQuality {
    fn default() -> Self {
        Self {
            energy_balane_error_abs: 0.0,
            energy_balane_error_rel: 0.0,

            sum_of_mass_fractions: Vec::new(),
            atomic_mass_balance_error: Vec::new(),
        }
    }
}
/// Error types for reactor modeling operations
#[derive(Debug)]
pub enum ReactorError {
    /// Missing required data (e.g., molar masses, kinetic parameters)
    MissingData(String),
    /// Invalid configuration (e.g., negative transport coefficients)
    InvalidConfiguration(String),
    /// Numeric input is non-finite or otherwise unusable for reactor physics
    InvalidNumericValue(String),
    /// Numerical calculation errors
    CalculationError(String),
    /// Parsing errors for chemical equations or parameters
    ParseError(String),
    /// Array/vector index out of bounds
    IndexOutOfBounds(String),
}

impl fmt::Display for ReactorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ReactorError::MissingData(msg) => write!(f, "Missing data: {}", msg),
            ReactorError::InvalidConfiguration(msg) => write!(f, "Invalid configuration: {}", msg),
            ReactorError::InvalidNumericValue(msg) => write!(f, "Invalid numeric value: {}", msg),
            ReactorError::CalculationError(msg) => write!(f, "Calculation error: {}", msg),
            ReactorError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            ReactorError::IndexOutOfBounds(msg) => write!(f, "Index out of bounds: {}", msg),
        }
    }
}

impl std::error::Error for ReactorError {}

/// Policy for assembling symbolic RHS expressions in `create_bvp_equations`.
///
/// `Auto` keeps small systems sequential and switches to parallel assembly for
/// larger ones where Rayon overhead is likely to pay off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymbolicRhsAssemblyPolicy {
    Auto,
    Sequential,
    Parallel,
}
/*
       (1/l)*  d(Lambda*( (1/l)* dT/d(x/l)  )  )/d(x/l) - c*m*( (1/l)* dT/d(x/l) + Q =0
        define: z = x/l
        d(Lambda* dT/dz   )/dz - c*m*l* dT/dz + l^2*Q =0
        define  Teta = (T-dT)/dT
         d(Lambda* dTeta/dz   )/dz - c*m*l* dTeta/dz + l^2*Q/dT =0
        define: q = Lambda* dTeta/dz
         d(Lambda* dTeta/dz   )/dz - c*m*l* dT/dz + l^2*Q/dT =0
        got equations:
        dT/dz = q/Lambda
        dq/dz - (c*m*l)/Lambda *q + (l^2)*Q/dT =0
        and finally:
        dT/dz = q/Lambda
        dq/dz - Pe_q*q + (l^2)*Q/dT =0

        (1/l)*  d(D*ro*( (1/l)* dCi/d(x/l)  )  )/d(x/l) - m*( (1/l)* dCi/d(x/l) + Gi =0
        define: z = x/l
        d(D*ro* dCi/dz   )/dz - m*l* dT/dz + l^2*Gi =0
        define: Ji = D*ro dC/dz
         d(Ji)/dz - m*l* dС/dz + l^2*Gi =0
        got equations:
        dCi/dz = Ji/D*ro;
        dJi/dz - (m*l)/ro*D *Ji + (l^2)*Gi =0

*/
/// Canonical in-memory description of an elementary reaction used by the BVP
/// pipeline.
///
/// `eq` stores the symbolic reaction equation, `A`, `n`, and `E` are the
/// Arrhenius parameters, and `Q` is the heat release contribution used by the
/// reactor energy balance.
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

/// Main reactor modeling structure that aggregates all reactor properties and methods
///
/// This structure handles the complete workflow from kinetic data to BVP solution:
/// 1. Kinetic preprocessing (stoichiometry, rate expressions)
/// 2. Transport property calculations (Peclet numbers, diffusion)
/// 3. Dimensionless scaling and equation setup
/// 4. BVP solving with multiple solver options
/// 5. Post-processing and visualization
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
    /// Pressure (Pa)
    pub P: f64,
    /// Mean temperature (K)
    pub Tm: f64,
    /// Heat capacity (J/kg·K)
    pub Cp: f64,
    /// Boundary conditions for substances and temperature
    pub boundary_condition: HashMap<String, f64>,
    /// Thermal conductivity (W/m·K)
    pub Lambda: f64,
    /// Diffusion coefficients for each substance (m²/s)
    pub Diffusion: HashMap<String, f64>,
    /// Mass flow rate (kg/s)
    pub m: f64,
    /// Scaling parameters for dimensionless transformation
    pub scaling: ScalingConfig,
    /// Characteristic length (m)
    pub L: f64,
    /// Temperature scaling expression: T = dT*(Teta + 1)
    pub T_scaling: Expr,
    /// Mean molar mass (kg/mol)
    pub M: f64,
    /// Transport coefficients D*ρ for each substance
    pub D_ro_map: HashMap<String, f64>,
    /// Thermal Peclet number: Pe_q = L*m*Cp/λ
    pub Pe_q: f64,
    /// Mass Peclet numbers for each substance: Pe_D = m*L/(D*ρ)
    pub Pe_D: Vec<f64>,
    /// Reaction rate expressions for each reaction
    pub map_eq_rate: HashMap<String, Expr>,
    /// System of differential equations (substance -> (variable, equation))
    pub map_of_equations: HashMap<String, (String, Expr)>,
    /// heat release function
    pub heat_release: Expr,
    /// Policy for symbolic RHS assembly in `create_bvp_equations`
    pub symbolic_rhs_assembly_policy: SymbolicRhsAssemblyPolicy,
    /// BVP solver instance
    pub solver: BVPSolver,
}
/// Boundary Value Problem solver wrapper
///
/// Supports multiple solver backends:
/// - NRBVP: Damped Newton-Raphson with finite differences
/// - BVPsci: Adaptive mesh refinement with collocation
#[derive(Debug, Clone)]
pub struct BVPSolver {
    /// Independent variable name (typically "x" or "z")
    pub arg_name: String,
    /// Domain range (start, end) - typically (0.0, 1.0) for dimensionless
    pub x_range: (f64, f64),
    /// Names of unknown variables ["Teta", "q", "C0", "J0", "C1", "J1", ...]
    pub unknowns: Vec<String>,
    /// System of differential equations dy/dx = f(x,y)
    pub eq_system: Vec<Expr>,
    /// Boundary conditions: variable -> (boundary_index, value)
    pub BorderConditions: HashMap<String, (usize, f64)>,
    /// Solution matrix (variables × mesh_points)
    pub solution: Option<DMatrix<f64>>,
    /// Spatial mesh points
    pub x_mesh: Option<DVector<f64>>,
    /// struct that stores balance errors - filled at the end of solution
    pub quality: SolutionQuality,
}

impl Default for BVPSolver {
    fn default() -> Self {
        Self {
            arg_name: "x".to_string(),
            x_range: (0.0, 1.0),
            unknowns: Vec::new(),
            eq_system: Vec::new(),
            BorderConditions: HashMap::new(),
            solution: None,
            x_mesh: None,
            quality: SolutionQuality::default(),
        }
    }
}

/// Validate the matrix returned by the BVP backend.
///
/// The backend may complete without a structured failure but still provide an
/// unusable solution buffer. This helper keeps the public contract strict by
/// rejecting empty or non-finite matrices before they enter reactor state.
pub(crate) fn validate_bvp_solution_matrix(solution: &DMatrix<f64>) -> Result<(), ReactorError> {
    if solution.nrows() == 0 || solution.ncols() == 0 {
        return Err(ReactorError::CalculationError(
            "BVP solver returned an empty solution matrix".to_string(),
        ));
    }
    if let Some((idx, value)) = solution
        .iter()
        .copied()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(ReactorError::InvalidNumericValue(format!(
            "BVP solver returned non-finite value {} at flat index {}",
            value, idx
        )));
    }
    Ok(())
}

/// Owned handoff payload for building an NRBVP backend snapshot.
///
/// The solver owns the canonical reactor state, while this struct packages the
/// backend-specific runtime settings so both direct and parser-driven solve
/// paths can assemble the same NRBVP object through one helper.
pub(crate) struct NrbvpHandoffConfig {
    pub initial_guess: DMatrix<f64>,
    pub t0: f64,
    pub t_end: f64,
    pub n_steps: usize,
    pub solver_backend_config: ReactorBvpSolverConfig,
    /// Optional solver-ready damped options built natively by RustedSciThe.
    ///
    /// When present, this takes precedence over the local facade conversion so
    /// parser-driven paths can hand the backend a ready-made solver contract.
    pub solver_options: Option<DampedSolverOptions>,
    pub scheme: String,
    pub strategy: String,
    pub strategy_params: Option<SolverParams>,
    pub linear_sys_method: Option<String>,
    pub method: String,
    pub abs_tolerance: f64,
    pub rel_tolerance: Option<HashMap<String, f64>>,
    pub max_iterations: usize,
    pub bounds: Option<HashMap<String, (f64, f64)>>,
    pub loglevel: Option<String>,
    pub dont_save_log: bool,
}

impl NrbvpHandoffConfig {
    /// Build a handoff config with the fields required by the backend.
    pub(crate) fn new(
        initial_guess: DMatrix<f64>,
        t0: f64,
        t_end: f64,
        n_steps: usize,
        scheme: String,
        strategy: String,
        strategy_params: Option<SolverParams>,
        linear_sys_method: Option<String>,
        method: String,
        abs_tolerance: f64,
        rel_tolerance: Option<HashMap<String, f64>>,
        max_iterations: usize,
        bounds: Option<HashMap<String, (f64, f64)>>,
        loglevel: Option<String>,
        dont_save_log: bool,
    ) -> Self {
        Self {
            initial_guess,
            t0,
            t_end,
            n_steps,
            solver_backend_config: ReactorBvpSolverConfig::default_lambdify(),
            solver_options: None,
            scheme,
            strategy,
            strategy_params,
            linear_sys_method,
            method,
            abs_tolerance,
            rel_tolerance,
            max_iterations,
            bounds,
            loglevel,
            dont_save_log,
        }
    }

    /// Override the generated callback backend while preserving legacy solver settings.
    pub(crate) fn with_solver_backend_config(
        mut self,
        solver_backend_config: ReactorBvpSolverConfig,
    ) -> Self {
        self.solver_backend_config = solver_backend_config;
        self
    }

    /// Attach a solver-ready damped options object and bypass the local facade conversion.
    pub(crate) fn with_solver_options(mut self, solver_options: DampedSolverOptions) -> Self {
        self.solver_options = Some(solver_options);
        self
    }
}

/// Apply legacy Newton/runtime settings on top of the generated-backend facade.
///
/// The matrix/execution/symbolic route comes from `solver_backend_config`, while
/// older public methods still pass scheme, strategy, tolerances, and bounds as
/// positional values. Keeping the merge here prevents parser/direct paths from
/// drifting while the public API is migrated gradually.
fn build_damped_options(config: &NrbvpHandoffConfig) -> Result<DampedSolverOptions, ReactorError> {
    let mut options = if let Some(solver_options) = config.solver_options.clone() {
        solver_options
    } else {
        config
            .solver_backend_config
            .to_rusted_options()
            .map_err(ReactorError::InvalidConfiguration)?
    };
    options = options
        .with_scheme_name(config.scheme.clone())
        .with_strategy_params(config.strategy_params.clone())
        .with_abs_tolerance(config.abs_tolerance)
        .with_max_iterations(config.max_iterations)
        .with_loglevel(config.loglevel.clone());
    options.strategy = config.strategy.clone();
    options.linear_sys_method = config.linear_sys_method.clone();
    options.method = config.method.clone();
    if let Some(rel_tolerance) = config.rel_tolerance.clone() {
        options = options.with_rel_tolerance(rel_tolerance);
    }
    if let Some(bounds) = config.bounds.clone() {
        options = options.with_bounds(bounds);
    }
    Ok(options)
}

impl BVPSolver {
    /// Convert the solver boundary-condition map into the backend shape.
    ///
    /// The backend expects each physical boundary value wrapped in a single-item
    /// vector, so we build that representation in one place and reuse it across
    /// solver entry points.
    pub(crate) fn backend_boundary_conditions(&self) -> HashMap<String, Vec<(usize, f64)>> {
        let mut backend_bc = HashMap::with_capacity(self.BorderConditions.len());
        for (name, &(boundary_index, value)) in &self.BorderConditions {
            backend_bc.insert(name.clone(), vec![(boundary_index, value)]);
        }
        backend_bc
    }

    /// Build the owned NRBVP backend object from the current solver snapshot.
    ///
    /// This keeps the ownership boundary explicit: the solver still owns the
    /// canonical equation system, but the backend gets its own owned copy at one
    /// well-defined handoff point.
    pub(crate) fn build_nrbvp_backend(
        &self,
        config: NrbvpHandoffConfig,
    ) -> Result<NRBVP, ReactorError> {
        let options = build_damped_options(&config)?;
        let mut backend = NRBVP::new_with_options(
            self.eq_system.clone(),
            config.initial_guess,
            self.unknowns.clone(),
            self.arg_name.clone(),
            self.backend_boundary_conditions(),
            config.t0,
            config.t_end,
            config.n_steps,
            options,
        );
        backend.dont_save_log(config.dont_save_log);
        Ok(backend)
    }

    /// Solve BVP using damped Newton-Raphson method (NRBVP)
    ///
    /// This is the main solver method with full parameter control
    pub fn solve_NRBVP(
        &mut self,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        scheme: String,
        strategy: String,
        strategy_params: Option<SolverParams>,
        linear_sys_method: Option<String>,
        method: String,
        abs_tolerance: f64,
        rel_tolerance: Option<HashMap<String, f64>>,
        max_iterations: usize,
        Bounds: Option<HashMap<String, (f64, f64)>>,
        loglevel: Option<String>,
    ) -> Result<(), ReactorError> {
        self.solve_NRBVP_impl(
            initial_guess,
            n_steps,
            scheme,
            strategy,
            strategy_params,
            linear_sys_method,
            method,
            abs_tolerance,
            rel_tolerance,
            max_iterations,
            Bounds,
            loglevel,
        )
    }

    /// Solve BVP with simplified tolerance configuration
    ///
    /// Uses ToleranceConfig to automatically generate tolerances for all variables
    pub fn solve_NRBVP_with_tolerance_config(
        &mut self,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        scheme: String,
        strategy: String,
        strategy_params: Option<SolverParams>,
        linear_sys_method: Option<String>,
        method: String,
        abs_tolerance: f64,
        tolerance_config: ToleranceConfig,
        substances: &[String],
        max_iterations: usize,
        Bounds: Option<HashMap<String, (f64, f64)>>,
        loglevel: Option<String>,
    ) -> Result<(), ReactorError> {
        let rel_tolerance = Some(tolerance_config.to_full_tolerance_map(substances));
        self.solve_NRBVP_impl(
            initial_guess,
            n_steps,
            scheme,
            strategy,
            strategy_params,
            linear_sys_method,
            method,
            abs_tolerance,
            rel_tolerance,
            max_iterations,
            Bounds,
            loglevel,
        )
    }

    /// Solve BVP with both tolerance and bounds configurations
    ///
    /// Most convenient method - automatically generates both tolerances and bounds
    /// for all variables from simple configs
    pub fn solve_NRBVP_with_configs(
        &mut self,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        scheme: String,
        strategy: String,
        strategy_params: Option<SolverParams>,
        linear_sys_method: Option<String>,
        method: String,
        abs_tolerance: f64,
        tolerance_config: ToleranceConfig,
        bounds_config: BoundsConfig,
        substances: &[String],
        max_iterations: usize,
        loglevel: Option<String>,
    ) -> Result<(), ReactorError> {
        let rel_tolerance = Some(tolerance_config.to_full_tolerance_map(substances));
        let bounds = Some(bounds_config.to_full_bounds_map(substances));
        self.solve_NRBVP_impl(
            initial_guess,
            n_steps,
            scheme,
            strategy,
            strategy_params,
            linear_sys_method,
            method,
            abs_tolerance,
            rel_tolerance,
            max_iterations,
            bounds,
            loglevel,
        )
    }

    /// Internal implementation for NRBVP solver
    ///
    /// All public solve methods delegate to this implementation
    fn solve_NRBVP_impl(
        &mut self,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        scheme: String,
        strategy: String,
        strategy_params: Option<SolverParams>,
        linear_sys_method: Option<String>,
        method: String,
        abs_tolerance: f64,
        rel_tolerance: Option<HashMap<String, f64>>,
        max_iterations: usize,
        Bounds: Option<HashMap<String, (f64, f64)>>,
        loglevel: Option<String>,
    ) -> Result<(), ReactorError> {
        info!("starting solver!");
        let mut bvp = self.build_nrbvp_backend(NrbvpHandoffConfig::new(
            initial_guess,
            self.x_range.0,
            self.x_range.1,
            n_steps,
            scheme,
            strategy,
            strategy_params,
            linear_sys_method,
            method,
            abs_tolerance,
            rel_tolerance,
            max_iterations,
            Bounds,
            loglevel,
            false,
        ))?;
        let solve_result = bvp
            .try_solve()
            .map_err(|err| ReactorError::CalculationError(format!("{err:?}")))?;
        let converged = solve_result.is_some();
        if !converged {
            warn!("BVP solver did not converge to a solution");
        }

        let solution = bvp.get_result().ok_or_else(|| {
            ReactorError::CalculationError(
                "BVP solver finished without producing a solution matrix".to_string(),
            )
        })?;
        validate_bvp_solution_matrix(&solution)?;
        if bvp.x_mesh.is_empty() {
            return Err(ReactorError::CalculationError(
                "BVP solver finished without producing an x mesh".to_string(),
            ));
        }
        if bvp.x_mesh.iter().any(|value| !value.is_finite()) {
            return Err(ReactorError::InvalidNumericValue(
                "BVP solver x mesh contains non-finite values".to_string(),
            ));
        }

        // Store solution for later access.
        self.solution = Some(solution);
        self.x_mesh = Some(bvp.x_mesh);

        if converged {
            Ok(())
        } else {
            Err(ReactorError::CalculationError(
                "BVP solver did not converge to a solution".to_string(),
            ))
        }
    }
    /// Solve BVP using adaptive mesh refinement (BVPsci)
    ///
    /// Alternative solver with automatic mesh adaptation
    pub fn solve_BVPsci(
        &mut self,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        max_nodes: usize,
        tol: f64,
    ) {
        let BC = self.backend_boundary_conditions();

        let mut bvp = BVPsci::new(
            None,
            Some(0.0 as f64),
            Some(1.0 as f64),
            Some(n_steps),
            self.eq_system.clone(),
            self.unknowns.clone(),
            vec![],
            None,
            BC,
            "x".to_string(),
            tol,
            max_nodes,
            initial_guess,
        );
        bvp.solve();
        if let Some(solution) = bvp.get_result() {
            self.solution = Some(solution);
            self.x_mesh = Some(bvp.mesh());
        }
    }
    /// Get reference to solution matrix
    ///
    /// Returns None if solve hasn't been called yet
    pub fn get_solution(&self) -> Option<&DMatrix<f64>> {
        self.solution.as_ref()
    }

    /// Debug print solution summary
    ///
    /// Prints solution matrix dimensions and sample values for each variable
    pub fn debug_solution(&self) {
        if let Some(solution) = &self.solution {
            info!("\n=== SOLUTION DEBUG ===");
            info!(
                "Solution matrix shape: {} x {}",
                solution.nrows(),
                solution.ncols()
            );
            info!("Unknowns: {:?}", self.unknowns);

            // Print the first and last few values for each variable.
            for (i, var_name) in self.unknowns.iter().enumerate() {
                if i < solution.ncols() {
                    let col = solution.column(i);
                    info!(
                        "{}: [first...{:.6}, {:.6}, {:.6}, ... last: {:.6}, {:.6}, {:.6}]",
                        var_name,
                        col[0],
                        col[1.min(col.len() - 1)],
                        col[2.min(col.len() - 1)],
                        col[col.len() - 3],
                        col[col.len() - 2],
                        col[col.len() - 1]
                    );
                }
            }
            info!("=== END DEBUG ===\n");
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
            P: 0.0,
            Tm: 0.0,
            Cp: 0.0,
            boundary_condition: HashMap::new(),
            Lambda: 0.0,
            Diffusion: HashMap::new(),

            m: 0.0,
            scaling: ScalingConfig::default(),
            L: 1.0,
            T_scaling: Expr::Const(0.0),
            M: 0.0,
            D_ro_map: HashMap::new(),
            Pe_q: 0.0,
            Pe_D: Vec::new(),
            map_eq_rate: HashMap::new(),
            map_of_equations: HashMap::new(),
            heat_release: Expr::Const(0.0),
            symbolic_rhs_assembly_policy: SymbolicRhsAssemblyPolicy::Auto,
            solver: BVPSolver::default(),
        }
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

    /// Set the symbolic RHS assembly policy in place.
    pub fn set_symbolic_rhs_assembly_policy(&mut self, policy: SymbolicRhsAssemblyPolicy) {
        self.symbolic_rhs_assembly_policy = policy;
    }

    /// Builder-style helper for configuring symbolic RHS assembly.
    ///
    /// This is convenient when the caller wants to configure the reactor in a
    /// fluent style before calling `setup_bvp`.
    pub fn with_symbolic_rhs_assembly_policy(mut self, policy: SymbolicRhsAssemblyPolicy) -> Self {
        self.symbolic_rhs_assembly_policy = policy;
        self
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
    /// Set boundary conditions for substances and temperature
    ///
    /// Keys should include substance names and "T" for temperature
    pub fn set_boundary_conditions(&mut self, conditions: HashMap<String, f64>) {
        self.boundary_condition = conditions;
    }

    /// Read-only access to the diffusion coefficients snapshot.
    ///
    /// New code should prefer this getter instead of reaching into the public
    /// `Diffusion` field directly.
    pub fn diffusion_coefficients(&self) -> &HashMap<String, f64> {
        &self.Diffusion
    }

    /// Read-only access to the raw boundary-condition snapshot.
    ///
    /// This keeps the legacy field available while giving external code a
    /// cleaner place to read from.
    pub fn boundary_conditions(&self) -> &HashMap<String, f64> {
        &self.boundary_condition
    }

    /// Read-only access to the cached transport coefficients `D * rho`.
    ///
    /// The cache is populated by `transport_coefficients()` and then reused by
    /// Lewis and Peclet number calculations.
    pub fn transport_cache(&self) -> &HashMap<String, f64> {
        &self.D_ro_map
    }

    /// Read-only access to the mass Peclet numbers.
    ///
    /// `Pe_D` is a domain-standard symbol, so the method keeps the scientific
    /// naming while still offering a stable accessor for new code.
    pub fn mass_peclet_numbers(&self) -> &[f64] {
        &self.Pe_D
    }

    /// Read-only access to the thermal Peclet number.
    pub fn thermal_peclet_number(&self) -> f64 {
        self.Pe_q
    }

    /// Read-only access to the temperature scaling shift `dT`.
    pub fn temperature_shift(&self) -> f64 {
        self.scaling.dT
    }

    /// Read-only access to the characteristic reactor length.
    pub fn characteristic_length(&self) -> f64 {
        self.scaling.L
    }

    /// Read-only access to the temperature scaling factor.
    pub fn temperature_scale(&self) -> f64 {
        self.scaling.T_scale
    }

    /// Read-only access to the full scaling configuration snapshot.
    ///
    /// New callers should prefer this over stitching `dT`, `L`, and `T_scale`
    /// together manually from the public fields.
    pub fn scaling_config(&self) -> &ScalingConfig {
        &self.scaling
    }

    /// Set all reactor parameters at once.
    ///
    /// `thermal_effects` are reaction heats in J/mol, `P` is pressure in Pa,
    /// `Tm` is temperature in K, `Cp` is heat capacity in J/(kg·K), `Lambda`
    /// is thermal conductivity in W/(m·K), `Diffusion` stores species
    /// diffusivities in m²/s, `m` is mass flow in kg/s, and `scaling`
    /// controls the dimensionless transformation.
    pub fn set_parameters(
        &mut self,
        thermal_effects: Vec<f64>,
        P: f64,
        Tm: f64,
        Cp: f64,
        boundary_condition: HashMap<String, f64>,
        Lambda: f64,
        Diffusion: HashMap<String, f64>,
        m: f64,

        scaling: ScalingConfig,
    ) {
        self.thermal_effects = thermal_effects;
        self.P = P;
        self.Tm = Tm;
        self.Cp = Cp;
        self.boundary_condition = boundary_condition;
        self.Lambda = Lambda;
        self.Diffusion = Diffusion;
        self.m = m;

        self.scaling = scaling;
    }
    /// Complete BVP setup workflow
    ///
    /// Orchestrates the entire setup process:
    /// 1. Kinetic preprocessing
    /// 2. Scaling and transport calculations  
    /// 3. Equation system assembly
    /// 4. Boundary condition setup
    pub fn setup_bvp(&mut self) -> Result<(), ReactorError> {
        //
        self.check_task()?;
        info!("task checked!");
        // Process
        self.scaling_processing()?;
        info!("scaling processed!");
        //
        // Process kinetics
        self.kinetic_processing()?;
        info!("kinetics processed!");

        // Calculate mean molar mass
        self.mean_molar_mass()?;
        info!("mean molar mass calculated");

        // dbg!("here");
        // Calculate Peclet numbers
        self.peclet_numbers()?;
        info!("Peclet numbers calculated");
        // Create BVP equations
        self.create_bvp_equations()?;
        info!("BVP equations created");
        self.set_solver_BC()?;
        info!("Boundary conditions created!");
        self.check_before_solution()?;
        info!("BVP setup completed!");

        Ok(())
    }

    /// Set transport properties.
    ///
    /// `lambda` is thermal conductivity in W/(m·K), `cp` is heat capacity in
    /// J/(kg·K), and `diffusion` maps species names to diffusivities in m²/s.
    pub fn set_transport_properties(
        &mut self,
        lambda: f64,
        cp: f64,
        diffusion: HashMap<String, f64>,
    ) {
        self.Lambda = lambda;
        self.Cp = cp;
        self.Diffusion = diffusion;
    }

    /// Set operating conditions.
    ///
    /// `pressure` is in Pa, `temperature` is in K, and `mass_flow` is in kg/s.
    pub fn set_operating_conditions(&mut self, pressure: f64, temperature: f64, mass_flow: f64) {
        self.P = pressure;
        self.Tm = temperature;
        self.m = mass_flow;
    }

    /// Set thermal effects (heat of reaction) for each reaction
    pub fn set_thermal_effects(&mut self, thermal_effects: Vec<f64>) {
        self.thermal_effects = thermal_effects;
    }

    /// Set elementary reactions from FastElemReact structures
    ///
    /// Convenient method to quickly set up elementary reactions with Arrhenius kinetics
    pub fn fast_react_set(&mut self, vec_of_maps: Vec<FastElemReact>) -> Result<(), ReactorError> {
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
            // Check for non-finite values (in case of uninitialized or corrupted f64)
            if !A.is_finite() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'A' in input hashmap at index {}",
                    idx
                )));
            }
            if !n.is_finite() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'n' in input hashmap at index {}",
                    idx
                )));
            }
            if !E.is_finite() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'E' in input hashmap at index {}",
                    idx
                )));
            }
            if !Q.is_finite() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'Q' in input hashmap at index {}",
                    idx
                )));
            }
            let arrenius = vec![A, n, E];
            let reactdata = ReactionData::new_elementary(eq.clone(), arrenius, None);
            Q_vec.push(Q);
            elementary_reaction_vec.push(reactdata);
        }

        let mut kindata = KinData::new();
        kindata
            .set_reaction_data_directly(elementary_reaction_vec, None)
            .map_err(|err| ReactorError::CalculationError(err.to_string()))?;
        self.kindata = kindata;
        self.thermal_effects = Q_vec;
        Ok(())
    }
    ///////////////////////////////////////////VALIDATION////////////////////////////////////////////////
    /// Validate reactor task configuration
    ///
    /// Checks:
    /// - All physical properties are positive
    /// - Diffusion coefficients exist for all substances
    /// - Thermal effects match number of reactions
    /// - Scaling parameters are valid
    /// - Boundary conditions are complete
    pub fn check_task(&self) -> Result<(), ReactorError> {
        let ensure_finite_positive = |name: &str, value: f64| -> Result<(), ReactorError> {
            if !value.is_finite() {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "{} must be finite, got {}",
                    name, value
                )));
            }
            if value <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "{} must be positive, got {}",
                    name, value
                )));
            }
            Ok(())
        };

        // Check basic properties
        ensure_finite_positive("P", self.P)?;
        ensure_finite_positive("Tm", self.Tm)?;
        ensure_finite_positive("Cp", self.Cp)?;
        ensure_finite_positive("Lambda", self.Lambda)?;
        ensure_finite_positive("m", self.m)?;

        ensure_finite_positive("M", self.M)?;

        // Check diffusion entries match substances
        if self.Diffusion.len() != self.kindata.substances.len() {
            return Err(ReactorError::InvalidConfiguration(
                "Diffusion entries must match number of substances".to_string(),
            ));
        }
        for substance in &self.kindata.substances {
            if !self.Diffusion.contains_key(substance) {
                return Err(ReactorError::MissingData(format!(
                    "Missing diffusion coefficient for {}",
                    substance
                )));
            }
            let diffusion = *self.Diffusion.get(substance).ok_or_else(|| {
                ReactorError::MissingData(format!(
                    "Missing diffusion coefficient for {}",
                    substance
                ))
            })?;
            if !diffusion.is_finite() || diffusion <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Diffusion coefficient for {} must be finite and positive, got {}",
                    substance, diffusion
                )));
            }
        }

        // Check thermal effects length
        if self.thermal_effects.len() != self.kindata.vec_of_equations.len() {
            return Err(ReactorError::InvalidConfiguration(
                "Thermal effects length must match number of reactions".to_string(),
            ));
        }
        for (idx, effect) in self.thermal_effects.iter().enumerate() {
            if !effect.is_finite() {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Thermal effect at index {} must be finite, got {}",
                    idx, effect
                )));
            }
        }

        // Validate scaling parameters
        self.scaling.validate()?;

        // Check boundary conditions
        let temperature = self.boundary_condition.get("T").ok_or_else(|| {
            ReactorError::MissingData("Missing T in boundary conditions".to_string())
        })?;
        if !temperature.is_finite() || *temperature <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Boundary temperature must be finite and positive, got {}",
                temperature
            )));
        }

        let mut boundary_fraction_sum = 0.0;
        for substance in &self.kindata.substances {
            let fraction = *self.boundary_condition.get(substance).ok_or_else(|| {
                ReactorError::MissingData(format!("Missing boundary condition for {}", substance))
            })?;
            if !fraction.is_finite() || fraction < 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Boundary fraction for {} must be finite and non-negative, got {}",
                    substance, fraction
                )));
            }
            boundary_fraction_sum += fraction;
        }
        if !boundary_fraction_sum.is_finite() || boundary_fraction_sum <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Boundary fractions must sum to a positive finite value, got {}",
                boundary_fraction_sum
            )));
        }
        if (boundary_fraction_sum - 1.0).abs() > 1e-6 {
            return Err(ReactorError::InvalidConfiguration(format!(
                "Boundary fractions must be normalized to 1.0, got {}",
                boundary_fraction_sum
            )));
        }

        Ok(())
    }
    /// Validate system before solving
    ///
    /// Checks that all arrays have consistent dimensions:
    /// - Equation system length = 2*n_substances + 2
    /// - Unknown variables, equations, and boundary conditions match
    /// - Peclet numbers are calculated and positive
    pub fn check_before_solution(&self) -> Result<(), ReactorError> {
        let n_substances = self.kindata.substances.len();
        let expected_len = 2 * n_substances + 2;

        if self.solver.eq_system.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "eq_system length {} != expected {}",
                self.solver.eq_system.len(),
                expected_len
            )));
        }
        if self.solver.unknowns.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "unknowns length {} != expected {}",
                self.solver.unknowns.len(),
                expected_len
            )));
        }
        if self.map_of_equations.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "map_of_equations length {} != expected {}",
                self.map_of_equations.len(),
                expected_len
            )));
        }
        if self.solver.BorderConditions.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "BorderConditions length {} != expected {}",
                self.solver.BorderConditions.len(),
                expected_len
            )));
        }
        if self.Pe_D.len() != n_substances {
            return Err(ReactorError::InvalidConfiguration(format!(
                "Pe_D length {} != substances {}",
                self.Pe_D.len(),
                n_substances
            )));
        }
        if !self.M.is_finite() || self.M <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "M must be finite and positive, got {}",
                self.M
            )));
        }
        if !self.Pe_q.is_finite() || self.Pe_q <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Pe_q must be finite and positive, got {}",
                self.Pe_q
            )));
        }
        for (idx, pe_d) in self.Pe_D.iter().enumerate() {
            if !pe_d.is_finite() || *pe_d <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Pe_D[{}] must be finite and positive, got {}",
                    idx, pe_d
                )));
            }
        }
        Ok(())
    }
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
    pub fn mean_molar_mass(&mut self) -> Result<(), ReactorError> {
        let mut mean_mass_inv = 0.0;
        let molar_masses = self
            .kindata
            .stecheodata
            .vec_of_molmasses
            .as_ref()
            .ok_or_else(|| ReactorError::MissingData("Molar masses not calculated".to_string()))?;
        // M_mean = sum_i(xi*Mi)
        //   x(i) = (ω(i) / M(i)) / ∑(ω(j) / M(j)) where ω(j) - is mass fruction
        // so M_men = ∑(xi*Mi) = ∑ω(i)  / ∑(ω(j) / M(j))= 1/∑(ω(j) / M(j))
        if !self.M.is_finite() || self.M <= 0.0 {
            for (i, substance) in self.kindata.substances.iter().enumerate() {
                if let Some(conc) = self.boundary_condition.get(substance) {
                    let mol_mass = molar_masses.get(i).ok_or_else(|| {
                        ReactorError::IndexOutOfBounds(format!(
                            "Molar mass index {} out of bounds",
                            i
                        ))
                    })?;
                    if !conc.is_finite() || *conc < 0.0 {
                        return Err(ReactorError::InvalidNumericValue(format!(
                            "Boundary fraction for {} must be finite and non-negative, got {}",
                            substance, conc
                        )));
                    }
                    if !mol_mass.is_finite() || *mol_mass <= 0.0 {
                        return Err(ReactorError::InvalidNumericValue(format!(
                            "Molar mass for {} must be finite and positive, got {}",
                            substance, mol_mass
                        )));
                    }
                    mean_mass_inv += conc / mol_mass;
                }
            }
            if !mean_mass_inv.is_finite() || mean_mass_inv <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Mean molar mass denominator must be finite and positive, got {}",
                    mean_mass_inv
                )));
            }
            let mean_mass = 1.0 / mean_mass_inv;
            if !mean_mass.is_finite() || mean_mass <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Mean molar mass must be finite and positive, got {}",
                    mean_mass
                )));
            }
            self.M = mean_mass / 1000.0; // from g/mol to kg/mol
        }
        Ok(())
    }
    ///
    /// Recompute the cached transport coefficients in place.
    ///
    /// The method updates `self.D_ro_map` directly so later hot-path reads can
    /// borrow the cached data instead of rebuilding or cloning the map.
    pub fn transport_coefficients(&mut self) -> Result<(), ReactorError> {
        /*
           D = A*T^1.5/P;
           D0  = A*T0^1.5/P0;
           D/D0 = (P0/P)*(T/T0)^1.5;
           D = D0*(P0/P)*(T/T0)^1.5;
           PV = (m/M)RT
           ro = P/(R*T)
           D*ro = (P/(R*T))*D0*(P0/P)*(T/T0)^1.5=(P0/(R*T0))*D0*(T/T0)^0.5=
           D0*ro0*(T/T0)^0.5 ;
           finally
           D*ro = D0*ro0*(T/T0)^0.5
        */
        if !self.P.is_finite() || self.P <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Pressure must be finite and positive, got {}",
                self.P
            )));
        }
        if !self.Tm.is_finite() || self.Tm <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Temperature must be finite and positive, got {}",
                self.Tm
            )));
        }
        if !self.M.is_finite() || self.M <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Mean molar mass must be finite and positive, got {}",
                self.M
            )));
        }
        // ro at standard conditions
        // PV = (m/M)*RT => ro = m/V = P*M/RT;
        let ro0 = self.M * self.P / (R_G * 298.15);
        if !ro0.is_finite() || ro0 <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Reference density must be finite and positive, got {}",
                ro0
            )));
        }
        let mut D_ro_map = HashMap::new();
        for subs in self.kindata.substances.iter() {
            if let Some(D_i) = self.Diffusion.get(subs) {
                if !D_i.is_finite() || *D_i <= 0.0 {
                    return Err(ReactorError::InvalidNumericValue(format!(
                        "Diffusion coefficient for {} must be finite and positive, got {}",
                        subs, D_i
                    )));
                }
                let D_ro = D_i * ro0 * (self.Tm / 298.15).powf(0.5);
                if !D_ro.is_finite() || D_ro <= 0.0 {
                    return Err(ReactorError::InvalidNumericValue(format!(
                        "D*ro for {} must be finite and positive, got {}",
                        subs, D_ro
                    )));
                }
                D_ro_map.insert(subs.clone(), D_ro);
            }
        }
        self.D_ro_map = D_ro_map;
        Ok(())
    }
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

        // Set characteristic length
        self.L = self.scaling.L;
        Ok(())
    }
    pub fn peclet_numbers(&mut self) -> Result<(), ReactorError> {
        if !self.Lambda.is_finite() || self.Lambda <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Lambda must be positive".to_string(),
            ));
        }
        if !self.m.is_finite() || self.m <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Mass flow rate must be positive".to_string(),
            ));
        }
        if !self.L.is_finite() || !self.Cp.is_finite() {
            return Err(ReactorError::InvalidNumericValue(
                "Peclet inputs must be finite".to_string(),
            ));
        }

        let Pe_q = (self.L * self.m * self.Cp) / self.Lambda;
        if !Pe_q.is_finite() || Pe_q <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Pe_q must be finite and positive, got {}",
                Pe_q
            )));
        }
        let mut Pe_D = Vec::new();
        self.transport_coefficients()?;

        for subs in self.kindata.substances.iter() {
            let ro_D_i = self.D_ro_map.get(subs).ok_or_else(|| {
                ReactorError::MissingData(format!(
                    "Transport coefficient for substance '{}' not found",
                    subs
                ))
            })?;
            if *ro_D_i <= 0.0 {
                return Err(ReactorError::InvalidConfiguration(format!(
                    "Transport coefficient for '{}' must be positive",
                    subs
                )));
            }
            let Pe_D_i = (self.m * self.L) / ro_D_i;
            if !Pe_D_i.is_finite() || Pe_D_i <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Pe_D for '{}' must be finite and positive, got {}",
                    subs, Pe_D_i
                )));
            }
            Pe_D.push(Pe_D_i);
        }

        self.Pe_q = Pe_q;
        self.Pe_D = Pe_D;
        Ok(())
    }
    /// Return the ideal-gas density for the current reactor state.
    ///
    /// The value is only meaningful when the thermodynamic inputs are finite
    /// and strictly positive, so we fail fast instead of silently returning
    /// `Inf` or `NaN`.
    pub fn ideal_gas_density(&self) -> Result<f64, ReactorError> {
        if !self.M.is_finite() || self.M <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Mean molar mass must be finite and positive, got {}",
                self.M
            )));
        }
        if !self.P.is_finite() || self.P <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Pressure must be finite and positive, got {}",
                self.P
            )));
        }
        if !self.Tm.is_finite() || self.Tm <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Temperature must be finite and positive, got {}",
                self.Tm
            )));
        }
        let density = self.M * self.P / (R_G * self.Tm);
        if !density.is_finite() || density <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Ideal gas density must be finite and positive, got {}",
                density
            )));
        }
        Ok(density)
    }
    /// Compute the average Lewis number from the current transport snapshot.
    ///
    /// The method borrows the diffusion map directly so it stays cheap and
    /// does not clone the full transport state on every call.
    pub fn Le_number(&self) -> Result<f64, ReactorError> {
        if self.D_ro_map.is_empty() {
            return Err(ReactorError::MissingData(
                "Cannot compute Lewis number without diffusion coefficients".to_string(),
            ));
        }
        if !self.Lambda.is_finite() || self.Lambda <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Thermal conductivity must be finite and positive, got {}",
                self.Lambda
            )));
        }
        if !self.Cp.is_finite() || self.Cp <= 0.0 {
            return Err(ReactorError::InvalidNumericValue(format!(
                "Heat capacity must be finite and positive, got {}",
                self.Cp
            )));
        }

        let mut le_sum = 0.0;
        let mut le_count = 0usize;

        for diffusion in self.D_ro_map.values() {
            if !diffusion.is_finite() || *diffusion <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Diffusion coefficient must be finite and positive, got {}",
                    diffusion
                )));
            }
            le_sum += self.Lambda / (self.Cp * diffusion);
            le_count += 1;
        }

        Ok(le_sum / le_count as f64)
    }
}

/////////////////////////////////////////TESTS/////////////////////////////////////////////////////////
