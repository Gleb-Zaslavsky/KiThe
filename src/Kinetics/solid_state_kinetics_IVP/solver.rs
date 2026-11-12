use super::models::*;
use super::options::SolidStateIvpSolverOptions;
use super::solution::IVPSolution;
use RustedSciThe::Utils::plots::{plots_ref, plots_terminal};
use RustedSciThe::numerical::BDF::BDF_api::BdfSolverOptions;
use RustedSciThe::numerical::BE::BeSolverOptions;
use RustedSciThe::numerical::LSODE2::{Lsode2ProblemConfig, Lsode2StopCondition};
use RustedSciThe::numerical::Radau::api::{RadauConfig, RadauOutputPolicy, RadauProblem};
use RustedSciThe::numerical::StiffAPI::{
    PreparedStiffSolver, StiffPreset, StiffSolveOutcome, StiffSolverKind, StiffSolverRequest,
    StiffTelemetrySnapshot,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};

/// Initial Value Problem solver for solid-state kinetic models with temperature dependence.
///
/// This struct provides a complete framework for solving kinetic problems involving
/// solid-state reactions with Arrhenius temperature dependence. It combines kinetic
/// models with thermal profiles and numerical ODE solving.
///
/// # Usage Pattern
/// 1. `new()` - Create instance with solver type
/// 2. `set_problem()` - Set thermal and kinetic parameters
/// 3. `set_model()` - Choose and configure kinetic model
/// 4. `check_task()` - Validate configuration
/// 5. `solve()` - Solve the ODE system
/// 6. `gnuplot()` or `save_result()` - Output results
///
/// # Fields
/// - Thermal parameters: `T0` (initial temperature), `beta` (heating rate), `E` (activation energy), `A` (pre-exponential factor)
/// - Time parameter: `t_final` (final integration time)
/// - Kinetic model: Combined Arrhenius and solid-state kinetic expressions
///
/// # Examples
/// ```rust, ignore
/// use your_crate::{KineticModelIVP, KineticModelNames};
/// use RustedSciThe::numerical::StiffAPI::StiffSolverKind;
///
/// let mut ivp = KineticModelIVP::new(StiffSolverKind::Bdf);
/// ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6)?;
/// ivp.set_model(KineticModelNames::JMA, vec![2.0])?;
/// ivp.solve()?;
/// ```
pub struct KineticModelIVP {
    /// Prepared StiffAPI solver (initialized after solve() is called).
    solver: Option<PreparedStiffSolver>,
    /// Materialized time-by-state trajectory, independent of native solver type.
    result: Option<(DVector<f64>, DMatrix<f64>)>,
    telemetry: Option<StiffTelemetrySnapshot>,
    /// Numerical solver parameters (step size, tolerance, etc.)
    solver_options: SolidStateIvpSolverOptions,
    /// Final integration time [s]
    t_final: f64,
    /// Heating rate [K/s]
    beta: f64,
    /// Initial temperature [K]
    T0: f64,
    /// Activation energy [J/mol]
    E: f64,
    /// Pre-exponential factor [1/s]
    A: f64,
    /// Selected kinetic model with parameters
    kinmodel: Option<KineticModel>,
    /// Combined kinetic expression (Arrhenius × kinetic model)
    /// Type of numerical solver to use
    solvertype: StiffSolverKind,
    stop_conversion_at: Option<f64>,
}

impl KineticModelIVP {
    /// Creates a new KineticModelIVP instance with default parameters.
    ///
    /// # Arguments
    /// * `solvertype` - The numerical solver type to use for ODE integration
    ///
    /// # Returns
    /// A new instance with default solver parameters and uninitialized problem data
    pub fn new(solvertype: StiffSolverKind) -> Self {
        Self {
            solver: None,
            result: None,
            telemetry: None,
            solver_options: SolidStateIvpSolverOptions::default(),
            t_final: 0.0,
            beta: 0.0,
            T0: 0.0,
            E: 0.0,
            A: 0.0,
            kinmodel: None,
            solvertype,
            stop_conversion_at: None,
        }
    }

    /// Convenience constructor for the default LSODE2 workflow.
    pub fn new_lsode2() -> Self {
        Self::new(StiffSolverKind::Lsode2)
    }

    pub fn set_stop_condition(&mut self, stop_conversion_at: Option<f64>) {
        self.stop_conversion_at = stop_conversion_at;
        self.invalidate_solution();
    }
    /// Sets the thermal and kinetic parameters for the problem.
    ///
    /// # Arguments
    /// * `t_final` - Final integration time [s] (must be > 0)
    /// * `beta` - Heating rate [K/s] (must be > 0)
    /// * `T0` - Initial temperature [K] (must be > 0)
    /// * `E` - Activation energy [J/mol] (must be > 0)
    /// * `A` - Pre-exponential factor [1/s] (must be > 0)
    ///
    /// # Returns
    /// * `Ok(())` - Parameters set successfully
    /// * `Err(String)` - Error message if any parameter is invalid (≤ 0)
    pub fn set_problem(
        &mut self,
        t_final: f64,
        beta: f64,
        T0: f64,
        E: f64,
        A: f64,
    ) -> Result<(), String> {
        if t_final <= 0.0 {
            return Err("t_final must be positive".to_string());
        }
        if beta <= 0.0 {
            return Err("beta must be positive".to_string());
        }
        if T0 <= 0.0 {
            return Err("T0 must be positive".to_string());
        }
        if E <= 0.0 {
            return Err("E must be positive".to_string());
        }
        if A <= 0.0 {
            return Err("A must be positive".to_string());
        }

        self.t_final = t_final;
        self.beta = beta;
        self.T0 = T0;
        self.E = E;
        self.A = A;
        self.invalidate_solution();
        Ok(())
    }
    pub fn set_solver_options(&mut self, options: SolidStateIvpSolverOptions) {
        self.solver_options = options;
        self.invalidate_solution();
    }

    pub fn set_arrhenius(&mut self, e: f64, a: f64) -> Result<(), String> {
        if e <= 0.0 {
            return Err("E must be positive".to_string());
        }
        if a <= 0.0 {
            return Err("A must be positive".to_string());
        }
        self.E = e;
        self.A = a;
        self.invalidate_solution();
        Ok(())
    }

    pub fn set_heating_program(&mut self, T0: f64, beta: f64) -> Result<(), String> {
        if beta <= 0.0 {
            return Err("beta must be positive".to_string());
        }
        if T0 <= 0.0 {
            return Err("T0 must be positive".to_string());
        }

        self.beta = beta;
        self.T0 = T0;
        self.invalidate_solution();

        Ok(())
    }

    pub fn set_time(&mut self, t_final: f64) -> Result<(), String> {
        if t_final <= 0.0 {
            return Err("t_final must be positive".to_string());
        }

        self.t_final = t_final;
        self.invalidate_solution();

        Ok(())
    }
    /// Sets the kinetic model and its parameters.
    ///
    /// This method validates and stores the kinetic model; current thermal and
    /// Arrhenius inputs are combined with it when a solve is prepared.
    ///
    /// # Arguments
    /// * `name` - The kinetic model type to use
    /// * `params` - Model parameters (empty vector for parameter-free models)
    ///
    /// # Returns
    /// * `Ok(())` - Model set successfully
    /// * `Err(String)` - Error message if parameters are invalid
    pub fn set_model(&mut self, name: KineticModelNames, params: Vec<f64>) -> Result<(), String> {
        self.kinmodel = Some(KineticModel::from_name_and_params(name, params)?);
        self.invalidate_solution();
        Ok(())
    }

    /// Validates that all required parameters and models are properly set.
    ///
    /// This method should be called before `solve()` to ensure the problem
    /// is completely configured.
    ///
    /// # Returns
    /// * `Ok(())` - All parameters are valid and set
    /// * `Err(String)` - Detailed error message indicating what is missing or invalid
    pub fn check_task(&self) -> Result<(), String> {
        if self.t_final <= 0.0 {
            return Err("t_final not set or invalid".to_string());
        }
        if self.beta <= 0.0 {
            return Err("beta not set or invalid".to_string());
        }
        if self.E <= 0.0 {
            return Err("E not set or invalid".to_string());
        }
        if self.A <= 0.0 {
            return Err("A not set or invalid".to_string());
        }
        if self.T0 <= 0.0 {
            return Err("T0 not set or invalid".to_string());
        }
        if self.kinmodel.is_none() {
            return Err("kinetic model not set".to_string());
        }
        Ok(())
    }
    /// Solves the kinetic ODE system.
    ///
    /// This method automatically calls `check_task()` first, then initializes
    /// and runs the ODE solver with the configured kinetic expression.
    ///
    /// # Returns
    /// * `Ok(())` - Solution computed successfully
    /// * `Err(String)` - Error message if validation fails or solver encounters issues
    pub fn solve(&mut self) -> Result<(), String> {
        self.check_task()?;
        self.invalidate_solution();
        self.solver_options.validate()?;
        let stop_target = self.stop_conversion_at.unwrap_or(0.999);
        if !stop_target.is_finite() || !(0.0..=1.0).contains(&stop_target) {
            return Err("stop conversion must be finite and within [0, 1]".into());
        }
        let equations = vec![self.kinetic_expression()?];
        let variables = vec!["a".to_owned()];
        let initial_state = DVector::from_vec(vec![1e-5]);
        let options = &self.solver_options;
        if matches!(&self.solvertype, StiffSolverKind::Lsode2) && options.max_steps.is_some() {
            return Err("LSODE2 does not expose a max_steps control through StiffAPI".into());
        }
        if matches!(&self.solvertype, StiffSolverKind::BackwardEuler)
            && options.first_step.is_some()
        {
            return Err("Backward Euler does not expose a first_step control".into());
        }

        let request = match self.solvertype {
            StiffSolverKind::Bdf => {
                let mut options = BdfSolverOptions::for_bdf(
                    equations,
                    variables,
                    "t".to_owned(),
                    0.0,
                    initial_state,
                    self.t_final,
                    options.max_step,
                    options.rtol,
                    options.atol,
                    None,
                    false,
                    options.first_step,
                );
                if let Some(max_steps) = self.solver_options.max_steps {
                    options.max_steps = max_steps;
                }
                StiffSolverRequest::Bdf { options }
            }
            StiffSolverKind::Radau => {
                let problem = RadauProblem::new(equations, variables, "t");
                let mut config = RadauConfig {
                    t0: 0.0,
                    t_bound: self.t_final,
                    rtol: options.rtol,
                    atol: options.atol,
                    first_step: options.first_step,
                    max_step: options.max_step,
                    output: RadauOutputPolicy::Dense,
                    ..RadauConfig::default()
                };
                if let Some(max_steps) = options.max_steps {
                    config.max_steps = max_steps;
                }
                StiffSolverRequest::Radau {
                    problem,
                    config,
                    initial_state: vec![1e-5],
                    parameters: Vec::new(),
                }
            }
            StiffSolverKind::BackwardEuler => {
                let mut be_options = BeSolverOptions::new(
                    equations,
                    variables,
                    "t".to_owned(),
                    options.newton_tolerance,
                    options.max_newton_iterations,
                    Some(options.max_step),
                    0.0,
                    self.t_final,
                    initial_state,
                );
                if let Some(max_steps) = options.max_steps {
                    be_options = be_options.with_max_steps(max_steps);
                }
                StiffSolverRequest::BackwardEuler {
                    options: be_options,
                    equation_parameters: Vec::new(),
                    parameter_values: Vec::new(),
                }
            }
            StiffSolverKind::Lsode2 => {
                let mut config = Lsode2ProblemConfig::new(
                    equations,
                    variables,
                    "t".to_owned(),
                    0.0,
                    initial_state,
                    self.t_final,
                    options.max_step,
                    options.rtol,
                    options.atol,
                );
                config.stop_conditions = vec![Lsode2StopCondition::greater_equal("a", stop_target)];
                StiffSolverRequest::Lsode2 {
                    config,
                    logging_level: None,
                }
            }
        };

        let preset = StiffPreset::solve_1d().with_telemetry(options.telemetry);
        let mut prepared = request
            .prepare(Some(&preset))
            .map_err(|error| format!("StiffAPI preparation failed: {error}"))?;
        let outcome = prepared
            .try_solve()
            .map_err(|error| format!("StiffAPI solve failed: {error}"))?;
        let (mut trajectory, telemetry) = match outcome {
            StiffSolveOutcome::Radau(solution) => (
                solution
                    .output_trajectory()
                    .map_err(|error| format!("Radau result extraction failed: {error}"))?,
                StiffTelemetrySnapshot::from_radau_solution(&solution),
            ),
            _ => {
                let view = prepared
                    .trajectory_view()
                    .map_err(|error| format!("StiffAPI result extraction failed: {error}"))?;
                let telemetry = prepared
                    .telemetry_snapshot()
                    .map_err(|error| format!("StiffAPI telemetry snapshot failed: {error}"))?;
                ((view.axis.clone(), view.values.clone()), telemetry)
            }
        };
        if let Some(first) =
            (0..trajectory.1.nrows()).find(|&row| trajectory.1[(row, 0)] >= stop_target)
        {
            let retained = first + 1;
            trajectory.0 = trajectory.0.rows(0, retained).into_owned();
            trajectory.1 = trajectory.1.rows(0, retained).into_owned();
        }
        self.result = Some(trajectory);
        self.telemetry = (options.telemetry != RustedSciThe::numerical::SolverTelemetryMode::Off)
            .then_some(telemetry);
        self.solver = Some(prepared);
        Ok(())
    }
    /// Plots the solution using gnuplot.
    ///
    /// # Returns
    /// * `Ok(())` - Plot generated successfully
    /// * `Err(String)` - Error if solver not initialized (call `solve()` first)
    pub fn plot(&self) -> Result<(), String> {
        let (time, values) = self.get_result()?;
        plots_ref("t", &["a".to_owned()], &time, &values);
        Ok(())
    }

    pub fn plot_in_terminal(&self) {
        // Skip terminal plotting when the solver has not produced a result yet.
        if let Ok((t, a)) = self.get_result() {
            plots_terminal("t".to_string(), vec!["a".to_string()], t, a)
        }
    }

    pub fn get_result(&self) -> Result<(DVector<f64>, DMatrix<f64>), String> {
        self.result
            .clone()
            .ok_or_else(|| "Solver not initialized. Call solve() first.".to_owned())
    }

    /// Native StiffAPI telemetry for the last successful solve, when enabled.
    pub fn telemetry_snapshot(&self) -> Option<&StiffTelemetrySnapshot> {
        self.telemetry.as_ref()
    }
    pub fn get_solution(&self) -> Result<IVPSolution, String> {
        let (t, a_matrix) = self.get_result()?;
        let model = self.kinmodel.as_ref().ok_or("kinetic model not set")?;
        IVPSolution::from_trajectory(t, a_matrix, self.T0, self.beta, self.E, self.A, model)
    }
    /// Saves the solution results to file.
    ///
    /// # Returns
    /// * `Ok(())` - Results saved successfully
    /// * `Err(String)` - Error if solver not initialized or file I/O fails
    pub fn save_result(&self) -> Result<(), String> {
        let solver = self
            .solver
            .as_ref()
            .ok_or("Solver not initialized. Call solve() first.")?;
        if let Some(bdf) = solver.as_bdf() {
            bdf.save_result()
                .map_err(|error| format!("Failed to save result: {error}"))?;
        }
        Ok(())
    }

    fn invalidate_solution(&mut self) {
        self.solver = None;
        self.result = None;
        self.telemetry = None;
    }

    fn arrhenius(&self) -> Expr {
        let E = Expr::Const(self.E);
        let A = Expr::Const(self.A);
        let beta = Expr::Const(self.beta);
        let R_sym = Expr::Const(8.31446261815324);
        let T0 = Expr::Const(self.T0);
        let t = Expr::Var("t".to_string());
        let k = A * Expr::Exp(Box::new(-E / (R_sym * (T0 + t * beta))));
        k
    }

    fn kinetic_expression(&self) -> Result<Expr, String> {
        let model = self
            .kinmodel
            .as_ref()
            .ok_or_else(|| "kinetic model not set".to_owned())?;
        Ok(self.arrhenius() * model.get_expr())
    }
}

impl Default for KineticModelIVP {
    /// Uses LSODE2 as the ergonomic default IVP solver.
    fn default() -> Self {
        Self::new_lsode2()
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
